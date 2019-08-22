/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
Modified: 2001 AlansFixes
**********/

/*
 * NIiter(ckt,maxIter)
 *
 *  This subroutine performs the actual numerical iteration.
 *  It uses the sparse matrix stored in the circuit struct
 *  along with the matrix loading program, the load data, the
 *  convergence test function, and the convergence parameters
 */

#include "ngspice/ngspice.h"
#include "ngspice/trandefs.h"
#include "ngspice/cktdefs.h"
#include "ngspice/smpdefs.h"
#include "ngspice/sperror.h"


/* NIiter() - return value is non-zero for convergence failure */

int
NIiter(CKTcircuit *ckt, int maxIter)
{
    double startTime, *OldCKTstate0 = NULL;
    int error, i, j;

    int iterno = 0;
    int ipass = 0;

#ifdef STEPDEBUG
    SPICE_debug(("entering...ckt->CTKmode=0x%08x\n", (uint32_t)(ckt->CKTmode)));
#endif

    /* some convergence issues that get resolved by increasing max iter */
    if (maxIter < 100)
        maxIter = 100;

    if ((ckt->CKTmode & MODETRANOP) && (ckt->CKTmode & MODEUIC)) {
        SWAP(double *, ckt->CKTrhs, ckt->CKTrhsOld);
        error = CKTload(ckt);
        if (error)
            return(error);
        return(OK);
    }

#ifdef WANT_SENSE2
    if (ckt->CKTsenInfo) {
        error = NIsenReinit(ckt);
        if (error)
            return(error);
    }
#endif

    if (ckt->CKTniState & NIUNINITIALIZED) {
        error = NIreinit(ckt);
        if (error) {
#ifdef STEPDEBUG
            SPICE_debug(("re-init returned error \n"));
#endif
            return(error);
        }
    }

    /* OldCKTstate0 = TMALLOC(double, ckt->CKTnumStates + 1); */

    for (;;) {

        ckt->CKTnoncon = 0;

#ifdef NEWPRED
        if (!(ckt->CKTmode & MODEINITPRED))
#endif
        {

#ifdef STEPDEBUG
            SPICE_debug(("#################### iterno=%d: noncon=%d, CKTmode=0x%08x\n", iterno, ckt->CKTnoncon, (uint32_t)(ckt->CKTmode)));
#endif
            error = CKTload(ckt);
            iterno++;
            if (error) {
                ckt->CKTstat->STATnumIter += iterno;
#ifdef STEPDEBUG
                SPICE_debug(("load returned error \n"));
#endif
                FREE(OldCKTstate0);
                return (error);
            }

#ifdef STEPDEBUG
            //SPICE_debug(("after loading, before solving\n"));
            //CKTdump(ckt, 0.0, NULL);
#endif

            if (!(ckt->CKTniState & NIDIDPREORDER)) {
                error = SMPpreOrder(ckt->CKTmatrix);
                if (error) {
                    ckt->CKTstat->STATnumIter += iterno;
#ifdef STEPDEBUG
                    SPICE_debug(("pre-order returned error \n"));
#endif
                    FREE(OldCKTstate0);
                    return(error); /* badly formed matrix */
                }
                ckt->CKTniState |= NIDIDPREORDER;
            }

            if ((ckt->CKTmode & MODEINITJCT) ||
                ((ckt->CKTmode & MODEINITTRAN) && (iterno == 1)))
            {
                ckt->CKTniState |= NISHOULDREORDER;
            }

            if (ckt->CKTniState & NISHOULDREORDER) {
                startTime = SPfrontEnd->IFseconds();
                error = SMPreorder(ckt->CKTmatrix, ckt->CKTpivotAbsTol,
                                   ckt->CKTpivotRelTol, ckt->CKTdiagGmin);
                ckt->CKTstat->STATreorderTime +=
                    SPfrontEnd->IFseconds() - startTime;
                if (error) {
                    /* new feature - we can now find out something about what is
                     * wrong - so we ask for the troublesome entry
                     */
                    SMPgetError(ckt->CKTmatrix, &i, &j);
                    SPfrontEnd->IFerrorf (ERR_WARNING, "singular matrix:  check nodes %s and %s\n", NODENAME(ckt, i), NODENAME(ckt, j));
                    ckt->CKTstat->STATnumIter += iterno;
#ifdef STEPDEBUG
                    SPICE_debug(("reorder returned error \n"));
#endif
                    FREE(OldCKTstate0);
                    return(error); /* can't handle these errors - pass up! */
                }
                ckt->CKTniState &= ~NISHOULDREORDER;
            } else {
                startTime = SPfrontEnd->IFseconds();
                error = SMPluFac(ckt->CKTmatrix, ckt->CKTpivotAbsTol,
                                 ckt->CKTdiagGmin);
                ckt->CKTstat->STATdecompTime +=
                    SPfrontEnd->IFseconds() - startTime;
                if (error) {
                    if (error == E_SINGULAR) {
                        ckt->CKTniState |= NISHOULDREORDER;
                        DEBUGMSG(" forced reordering....\n");
                        continue;
                    }
                    /* CKTload(ckt); */
                    /* SMPprint(ckt->CKTmatrix, stdout); */
                    /* seems to be singular - pass the bad news up */
                    ckt->CKTstat->STATnumIter += iterno;
#ifdef STEPDEBUG
                    SPICE_debug(("lufac returned error \n"));
#endif
                    FREE(OldCKTstate0);
                    return(error);
                }
            }

            /* moved it to here as if xspice is included then CKTload changes
               CKTnumStates the first time it is run */
            if (!OldCKTstate0)
                OldCKTstate0 = TMALLOC(double, ckt->CKTnumStates + 1);
            memcpy(OldCKTstate0, ckt->CKTstate0,
                   (size_t) ckt->CKTnumStates * sizeof(double));

            startTime = SPfrontEnd->IFseconds();
            SMPsolve(ckt->CKTmatrix, ckt->CKTrhs, ckt->CKTrhsSpare);
            ckt->CKTstat->STATsolveTime +=
                SPfrontEnd->IFseconds() - startTime;
#ifdef STEPDEBUG
            /*XXXX*/
            if (ckt->CKTrhs[0] != 0.0)
                SPICE_debug(("NIiter: CKTrhs[0] = %g\n", ckt->CKTrhs[0]));
            if (ckt->CKTrhsSpare[0] != 0.0)
                SPICE_debug(("NIiter: CKTrhsSpare[0] = %g\n", ckt->CKTrhsSpare[0]));
            if (ckt->CKTrhsOld[0] != 0.0)
                SPICE_debug(("NIiter: CKTrhsOld[0] = %g\n", ckt->CKTrhsOld[0]));
            /*XXXX*/
#endif
            ckt->CKTrhs[0] = 0;
            ckt->CKTrhsSpare[0] = 0;
            ckt->CKTrhsOld[0] = 0;

            if (iterno > maxIter) {
                ckt->CKTstat->STATnumIter += iterno;
                /* we don't use this info during transient analysis */
                if (ckt->CKTcurrentAnalysis != DOING_TRAN) {
                    FREE(errMsg);
                    errMsg = copy("Too many iterations without convergence");
#ifdef STEPDEBUG
                    fprintf(stderr, "too many iterations without convergence: %d iter's (max iter == %d)\n",
                    iterno, maxIter);
#endif
                }
                FREE(OldCKTstate0);
                return(E_ITERLIM);
            }

            if ((ckt->CKTnoncon == 0) && (iterno != 1))
                ckt->CKTnoncon = NIconvTest(ckt);
            else
                ckt->CKTnoncon = 1;

#ifdef STEPDEBUG
            SPICE_debug(("noncon is %d\n", ckt->CKTnoncon));
#endif
        } /* if (!(ckt->CKTmode & MODEINITPRED)) */

        if ((ckt->CKTnodeDamping != 0) && (ckt->CKTnoncon != 0) &&
            ((ckt->CKTmode & MODETRANOP) || (ckt->CKTmode & MODEDCOP)) &&
            (iterno > 1)) {

#ifdef STEPDEBUG
            SPICE_debug(("enter some if() block.\n"));
#endif

            CKTnode *node;
            double diff, maxdiff = 0;
            for (node = ckt->CKTnodes->next; node; node = node->next)
                if (node->type == SP_VOLTAGE) {
                    diff = ckt->CKTrhs[node->number] - ckt->CKTrhsOld[node->number];
                    if (maxdiff < diff)
                        maxdiff = diff;
                }

            if (maxdiff > 10) {
                double damp_factor = 10 / maxdiff;
                if (damp_factor < 0.1)
                    damp_factor = 0.1;
                for (node = ckt->CKTnodes->next; node; node = node->next) {
                    diff = ckt->CKTrhs[node->number] - ckt->CKTrhsOld[node->number];
                    ckt->CKTrhs[node->number] =
                        ckt->CKTrhsOld[node->number] + (damp_factor * diff);
                }
                for (i = 0; i < ckt->CKTnumStates; i++) {
                    diff = ckt->CKTstate0[i] - OldCKTstate0[i];
                    ckt->CKTstate0[i] = OldCKTstate0[i] + (damp_factor * diff);
                }
            }
        }

        if (ckt->CKTmode & MODEINITFLOAT) {
#ifdef STEPDEBUG
            SPICE_debug(("enter if (ckt->CKTmode & MODEINITFLOAT): noncon=%d\n", ckt->CKTnoncon));
#endif
            if ((ckt->CKTmode & MODEDC) && ckt->CKThadNodeset) {
#ifdef STEPDEBUG
                SPICE_debug(("enter if if ((ckt->CKTmode & MODEDC) && ckt->CKThadNodeset)\n"));
#endif
                
                if (ipass) {
#ifdef STEPDEBUG
                    SPICE_debug(("enter if if (ipass): ipass=%d\n", ipass));
#endif
                    ckt->CKTnoncon = ipass;
                }
                ipass = 0;
            }
            if (ckt->CKTnoncon == 0) {
                ckt->CKTstat->STATnumIter += iterno;
                FREE(OldCKTstate0);
#ifdef STEPDEBUG
                SPICE_debug(("return OK!!!!!!!!!!!!!!!!!!!\n"));
#endif
                return(OK);
            }
            
        } else if (ckt->CKTmode & MODEINITJCT) {
            ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFIX;
            ckt->CKTniState |= NISHOULDREORDER;
        } else if (ckt->CKTmode & MODEINITFIX) {
            if (ckt->CKTnoncon == 0)
                ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFLOAT;
#ifdef STEPDEBUG
            SPICE_debug(("ipass=1\n"));
#endif
            ipass = 1;
        } else if (ckt->CKTmode & MODEINITSMSIG) {
            ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFLOAT;
        } else if (ckt->CKTmode & MODEINITTRAN) {
            if (iterno <= 1)
                ckt->CKTniState |= NISHOULDREORDER;
            ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFLOAT;
        } else if (ckt->CKTmode & MODEINITPRED) {
            ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFLOAT;
        } else {
            ckt->CKTstat->STATnumIter += iterno;
#ifdef STEPDEBUG
            SPICE_debug(("bad initf state \n"));
#endif
            FREE(OldCKTstate0);
            return(E_INTERN);
            /* impossible - no such INITF flag! */
        }

        /* build up the lvnim1 array from the lvn array */
#ifdef STEPDEBUG
        SPICE_debug(("swapping CKTrhs and CKTrhsOld...\n"));
#endif
        SWAP(double *, ckt->CKTrhs, ckt->CKTrhsOld);
#ifdef STEPDEBUG
        //SPICE_debug(("after loading, after solving\n"));
        //CKTdump(ckt, 0.0, NULL);
#endif
    } /* end for(;;) */
    /*NOTREACHED*/
}
