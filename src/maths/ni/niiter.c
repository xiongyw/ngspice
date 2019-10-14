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

extern bool g_near;
extern int g_num_workers;

int fake_near(CKTcircuit *ckt, int maxIter, int num_workers);
int near(CKTcircuit *ckt, int maxIter, int num_workers);


/* NIiter() - return value is non-zero for convergence failure */

int NIiter(CKTcircuit *ckt, int maxIter)
{
    double startTime, *OldCKTstate0 = NULL;
    int error, i, j;

    int iterno = 0;
    int ipass = 0;

    SPICE_debug(("entering...ckt->CTKmode=0x%08x\n", (uint32_t)(ckt->CKTmode)));

    /* some convergence issues that get resolved by increasing max iter */
    if (maxIter < 100)
        maxIter = 100;


    

    /*
     * test EA for DCOP 
     */
    if (g_near && (ckt->CKTmode & MODEDCOP)) {
        return near(ckt, maxIter, g_num_workers);
        //return fake_near(ckt, maxIter, g_num_workers);
    }


    

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

            SPICE_debug(("#################### iterno=%d: CKTnoncon=%d, CKTmode=0x%08x\n", iterno, ckt->CKTnoncon, (uint32_t)(ckt->CKTmode)));
            error = CKTload(ckt);
            iterno++;
            if (error) {
                ckt->CKTstat->STATnumIter += iterno;

                SPICE_debug(("load returned error \n"));

                FREE(OldCKTstate0);
                return (error);
            }


            SPICE_debug(("%d: after CKTload(), CKTmode=0x%08x:\n", iterno, ckt->CKTmode));
            for (i = 0; i < SMPmatSize(ckt->CKTmatrix) + 1; ++ i) {
                SPICE_debug(("rhs, rhsOld, spare[%d]=%f, %f, %f\n", i, ckt->CKTrhs[i], ckt->CKTrhsOld[i], ckt->CKTrhsSpare[i]));
            }



            if (!(ckt->CKTniState & NIDIDPREORDER)) {
                error = SMPpreOrder(ckt->CKTmatrix);
                if (error) {
                    ckt->CKTstat->STATnumIter += iterno;
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
                    FREE(OldCKTstate0);
                    return(error);
                }
            }

            /* moved it to here as if xspice is included then CKTload changes
               CKTnumStates the first time it is run */
            if (!OldCKTstate0)
                OldCKTstate0 = TMALLOC(double, ckt->CKTnumStates + 1);
            memcpy(OldCKTstate0, ckt->CKTstate0, (size_t) ckt->CKTnumStates * sizeof(double));

            startTime = SPfrontEnd->IFseconds();
            SMPsolve(ckt->CKTmatrix, ckt->CKTrhs, ckt->CKTrhsSpare);
            ckt->CKTstat->STATsolveTime += SPfrontEnd->IFseconds() - startTime;

            SPICE_debug(("%d: after SMPsolve(), CKTmode=0x%08x, CKTnoncon=%d:\n", iterno, ckt->CKTmode, ckt->CKTnoncon));
            for (i = 0; i < SMPmatSize(ckt->CKTmatrix) + 1; ++ i) {
                SPICE_debug(("rhs, rhsOld, spare[%d]=%f, %f, %f\n", i, ckt->CKTrhs[i], ckt->CKTrhsOld[i], ckt->CKTrhsSpare[i]));
            }

            /*XXXX*/
            if (ckt->CKTrhs[0] != 0.0)
                SPICE_debug(("NIiter: CKTrhs[0] = %g\n", ckt->CKTrhs[0]));
            if (ckt->CKTrhsSpare[0] != 0.0)
                SPICE_debug(("NIiter: CKTrhsSpare[0] = %g\n", ckt->CKTrhsSpare[0]));
            if (ckt->CKTrhsOld[0] != 0.0)
                SPICE_debug(("NIiter: CKTrhsOld[0] = %g\n", ckt->CKTrhsOld[0]));
            /*XXXX*/

            ckt->CKTrhs[0] = 0;
            ckt->CKTrhsSpare[0] = 0;
            ckt->CKTrhsOld[0] = 0;

            if (iterno > maxIter) {
                ckt->CKTstat->STATnumIter += iterno;
                /* we don't use this info during transient analysis */
                if (ckt->CKTcurrentAnalysis != DOING_TRAN) {
                    FREE(errMsg);
                    errMsg = copy("Too many iterations without convergence");

                    fprintf(stderr, "too many iterations without convergence: %d iter's (max iter == %d)\n",
                    iterno, maxIter);

                }
                FREE(OldCKTstate0);
                return(E_ITERLIM);
            }


            if ((ckt->CKTnoncon == 0) && (iterno != 1)) {
                ckt->CKTnoncon = NIconvTest(ckt);
                SPICE_debug(("NIconvTest()=%d\n", ckt->CKTnoncon));
            } else {
                ckt->CKTnoncon = 1;
            }


        } /* if (!(ckt->CKTmode & MODEINITPRED)) */

        if ((ckt->CKTnodeDamping != 0) && (ckt->CKTnoncon != 0) &&
            ((ckt->CKTmode & MODETRANOP) || (ckt->CKTmode & MODEDCOP)) &&
            (iterno > 1)) {

            SPICE_debug(("enter some if() block.\n"));

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
            SPICE_debug(("enter if (ckt->CKTmode & MODEINITFLOAT): noncon=%d\n", ckt->CKTnoncon));
            if ((ckt->CKTmode & MODEDC) && ckt->CKThadNodeset) {
                SPICE_debug(("enter if if ((ckt->CKTmode & MODEDC) && ckt->CKThadNodeset)\n"));

                if (ipass) {
                    SPICE_debug(("enter if if (ipass): ipass=%d\n", ipass));
                    ckt->CKTnoncon = ipass;
                }
                ipass = 0;
            }
            if (ckt->CKTnoncon == 0) {
                ckt->CKTstat->STATnumIter += iterno;
                FREE(OldCKTstate0);

                printf("converged! iterno=%d\n", iterno);

                return(OK);
            }

        } else if (ckt->CKTmode & MODEINITJCT) {
            ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFIX;
            ckt->CKTniState |= NISHOULDREORDER;
        } else if (ckt->CKTmode & MODEINITFIX) {
            if (ckt->CKTnoncon == 0)
                ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFLOAT;

            SPICE_debug(("ipass=1\n"));
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

            SPICE_debug(("bad initf state \n"));

            FREE(OldCKTstate0);
            return(E_INTERN);
            /* impossible - no such INITF flag! */
        }

        /* build up the lvnim1 array from the lvn array */

        SPICE_debug(("swapping CKTrhs and CKTrhsOld...\n"));

        SWAP(double *, ckt->CKTrhs, ckt->CKTrhsOld);

        //SPICE_debug(("after loading, after solving\n"));
        //CKTdump(ckt, 0.0, NULL);

    } /* end for(;;) */
    /*NOTREACHED*/
}








#include <assert.h>


/* rhsOld is used as the current OP */
void set_rhsOld(CKTcircuit* ckt, double* x0)
{
    int size = SMPmatSize(ckt->CKTmatrix);
    memcpy(ckt->CKTrhsOld, x0, sizeof(double) * (size + 1));
}

typedef struct {
    double* x;
    double* dx;
    double dist;
    bool is_conv;
    double* state0;
}worker_io; // workers' input/output

void deallocate_worker_io(worker_io* w, int N)
{
    for (int j = 0; j < N; j ++) {
        free(w[j].x);
        free(w[j].dx);
        free(w[j].state0);
    }
    free(w);
}


worker_io* allocate_worker_io(int N, int matrix_size, int CKTnumStates)
{
    int i;
    worker_io* w = (worker_io*)malloc(sizeof(worker_io) * N);
    if (!w) return NULL;
    for (i = 0; i < N; i ++) {
        w[i].x = (double*)malloc(sizeof(double) * (matrix_size + 1));
        if (!w[i].x) {
            deallocate_worker_io(w, i);
            return NULL;
        }
        
        w[i].dx = (double*)malloc(sizeof(double) * (matrix_size + 1));
        if (!w[i].dx) {
            free(w[i].x);
            deallocate_worker_io(w, i);
            return NULL;
        }
        
        w[i].dist = 0;
        w[i].is_conv = FALSE;

        
        w[i].state0 = (double*)malloc(sizeof(double) * CKTnumStates);
        
    }
    
    return w;

}



/*
 * for a given ckt and OP x, get x' (store in ckt->CKTrhs), and calculate dist(x, x') and noncon.
 */
int worker_whip(CKTcircuit* ckt, worker_io* io, bool is_1st_iter)
{
    // set OP
    int size = SMPmatSize(ckt->CKTmatrix);
    memcpy(ckt->CKTrhs, io->x, sizeof(double) * (size + 1));
    SWAP(double *, ckt->CKTrhs, ckt->CKTrhsOld);

    // set state0
    memcpy(ckt->CKTstate0, io->state0, (size_t) ckt->CKTnumStates * sizeof(double));

    SPICE_debug(("before CKTload()\n"));
    for (int i = 0; i < size + 1; ++ i) {
        SPICE_debug(("rhs, rhsOld, spare[%d]=%f, %f, %f\n", i, ckt->CKTrhs[i], ckt->CKTrhsOld[i], ckt->CKTrhsSpare[i]));
    }

    // load
    ckt->CKTnoncon = 0;
    CKTload(ckt);
    
    SPICE_debug(("after CKTload()\n"));
    for (int i = 0; i < size + 1; ++ i) {
        SPICE_debug(("rhs, rhsOld, spare[%d]=%f, %f, %f\n", i, ckt->CKTrhs[i], ckt->CKTrhsOld[i], ckt->CKTrhsSpare[i]));
    }

    // solve
    if (is_1st_iter) {
        SMPpreOrder(ckt->CKTmatrix);
        SMPreorder(ckt->CKTmatrix, ckt->CKTpivotAbsTol, ckt->CKTpivotRelTol, ckt->CKTdiagGmin);
        SMPsolve(ckt->CKTmatrix, ckt->CKTrhs, ckt->CKTrhsSpare);
    } else {
        SMPluFac(ckt->CKTmatrix, ckt->CKTpivotAbsTol, ckt->CKTdiagGmin);
        SMPsolve(ckt->CKTmatrix, ckt->CKTrhs, ckt->CKTrhsSpare);
    }

    SPICE_debug(("after solve:\n"));
    for (int i = 0; i < size + 1; ++ i) {
        SPICE_debug(("rhs, rhsOld, spare[%d]=%f, %f, %f\n", i, ckt->CKTrhs[i], ckt->CKTrhsOld[i], ckt->CKTrhsSpare[i]));
    }


    // check converge
    if (ckt->CKTnoncon == 0) {
        ckt->CKTnoncon = NIconvTest(ckt);
    } else {
        ckt->CKTnoncon = 1;
    }

    /*
     * output
     */
    double tmp = 0.;
    for (int i = 0; i < size + 1; i ++) {
        io->dx[i] = ckt->CKTrhs[i] - io->x[i];
        tmp +=  io->dx[i] * io->dx[i];
    }
    io->dist = tmp;
    io->is_conv = (ckt->CKTnoncon == 0)? TRUE : FALSE;

    // store state0
    memcpy(io->state0, ckt->CKTstate0, (size_t) ckt->CKTnumStates * sizeof(double));

    return 0;
}

int set_worker_io(worker_io* w, int N, const double* x, const double* dx, double ratio, int matrix_size,
                  const double* state0, int CKTnumStates)
{
    double *center = NULL;
    double *step = NULL;
    
    center = (double*)malloc(sizeof(double) * (matrix_size + 1));
    assert(center);
    step = (double*)malloc(sizeof(double) * (matrix_size + 1));
    assert(step);

    for (int i = 0; i < matrix_size + 1; i ++) {
        center[i] = x[i] + dx[i];
        step[i] = dx[i] * 2.0 * ratio / N;
    }

    for (int i = 0; i < N; i ++) {

        // candidates
        for (int j = 0; j < matrix_size + 1; j ++) {
            if (N == 1) { // for testing
                w[i].x[j] = center[j];
            } else {
                w[i].x[j] = center[j] + step[j] * (i - N / 2);
            }
            //SPICE_debug(("%d: x[%d]=%f\n", i, j, w[i].x[j]));
        }

        // state0
        memcpy(w[i].state0, state0, sizeof(double) * CKTnumStates);
    }

    free(center);
    free(step);
    
    return 0;
}

#if (0)
simple nonlinear circuit
.global gnd
.model 1n4148 d (is=0.1pa, rs=16 cjo=2pf tt=12n bv=100 ibv=1na)
r 2 0 10k
d 1 2 1n4148
v 1 0 dc 10
.control
debug
op
.endc
    
       1        d       2
       +-------|>|------+
       |     1N4148     |
       |                |
    +  -                x
    v(10)               x r(10k)
    -  -                x
       |                |
       |                |
       +----------------+
       0
#endif
/* near() only works for OP: Evolutionary Algorithm for Newton-Raphson iterations. */
int near(CKTcircuit *ckt, int maxIter, int N /* number of workers*/)
{
    int error, i, j, k;

    int iterno = 0;
    double ratio = 1.;
    int nconv = 0;  // number of converged workers

    SPICE_debug(("entering...ckt->CTKmode=0x%08x\n", (uint32_t)(ckt->CKTmode)));

    // allocate worker_io
    int size = SMPmatSize(ckt->CKTmatrix);
    worker_io* w = allocate_worker_io(N + 1, size, ckt->CKTnumStates); // slot 0 stores the best candidate
    assert(w);

    /*
     * it seems that the CKTmode & SMPxxx() calls change in the following way (at least for OP):
     * 
     * =========================================================================
     *  iterations  |    CKTmode            | SMPxxx() calls
     * =========================================================================
     *  1st         | 0x210 (INITJCT)       | preorder + reorder + solve
     * -------------+-----------------------+-----------------------------------
     *  rest        | 0x410 (INITFIX)       |
     * -------------+-----------------------| luFac + solve
     *  last        | 0x110 (INITFLOAT)     | 
     * =========================================================================
     *
     * NB: 
     * - CKTmode affects how CKTload() works (fixme: how?)
     * - SMPxxx() calls solve the matrix
     */


    /*
     * the 1st iteration, use w[0] as io
     */
    for (i = 0; i < size + 1; ++ i) {
        w[0].x[i] = 0.0;
    }
    memcpy(w[0].state0, ckt->CKTstate0, sizeof(double)*ckt->CKTnumStates);
    ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITJCT;
    worker_whip(ckt, w, TRUE);

    SPICE_debug(("w[0].dist=%.16f, w[0].is_conv=%d\n", w[0].dist, w[0].is_conv));
    for (i = 0; i < size + 1; ++ i) {
        SPICE_debug(("w[0].dx[%d]=%.16f\n", i, w[0].dx[i]));
    }

    ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFIX;

    /*
     * the rest iteration
     */
    for (k = 1; k < maxIter; k ++) {
        
        SPICE_debug(("######################>> iterno=%d: noncon=%d, CKTmode=0x%08x\n", k, ckt->CKTnoncon, (uint32_t)(ckt->CKTmode)));

        set_worker_io(w + 1, N, w[0].x, w[0].dx, ratio, size, w[0].state0, ckt->CKTnumStates);

        for (i = 1; i < N + 1; i ++) {
            
            worker_whip(ckt, w + i, FALSE);
        
            SPICE_debug(("w[%d].dist=%.16f, w[%d].is_conv=%d\n", i, w[i].dist, i, w[i].is_conv));
            for (j = 0; j < size + 1; ++ j) {
                SPICE_debug(("w[%d].dx[%d]=%10.9f\n", i, w[i].dx[j]));
            }
            
            if (w[i].is_conv) {
                printf("ckt->CKTmode=0x%08x: w[%3d] converged, dist=%.16f! N=%d, iterno=%d\n", ckt->CKTmode, i, w[i].dist, N, k+1);
                nconv ++;
            }
        }

        // find the minimum distance
        int min_idx = 1;
        double min = w[1].dist;
        for (i = 2; i < N + 1; ++ i) {
            if (w[i].dist < min) {
                min = w[i].dist;
                min_idx = i;
            }
        }
        SPICE_debug(("k=%d: min_idx=%d, min dist=%.16f\n", k, min_idx, w[min_idx].dist));

        // copy the best candidate to the first slot
        memcpy(w[0].x, w[min_idx].x, sizeof(double) * (size + 1));
        memcpy(w[0].dx, w[min_idx].dx, sizeof(double) * (size + 1));
        w[0].dist = w[min_idx].dist;
        w[0].is_conv = w[min_idx].is_conv;
        memcpy(w[0].state0, w[min_idx].state0, sizeof(double) * ckt->CKTnumStates);


        if (nconv != 0) {
            if (ckt->CKTmode & MODEINITFIX) {
                ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFLOAT;
            } else if (ckt->CKTmode & MODEINITFLOAT) {
                printf("MODEINITFLOAT: converged!!!!!!!!!!!!!! N=%d, iterno=%d\n", N, k+1);
                break;
            }
        }
        
    }


    if((ckt->CKTmode & MODEINITFLOAT) && nconv != 0) {
        SPICE_debug(("final result: w[0].dist=%.16f, w[0].is_conv=%d\n", w[0].dist, w[0].is_conv));
        for (i = 0; i < size + 1; ++ i) {
            SPICE_debug(("w[0].dx[%d]=%.16f\n", i, w[0].dx[i]));
        }

        /* copy final result from w[0] to CKT */
        for (i = 0; i < size + 1; i ++) {
            ckt->CKTrhs[i] = w[0].x[i] + w[0].dx[1];
            ckt->CKTrhsOld[i] = w[0].x[i];
        }
        memcpy(ckt->CKTstate0, w[0].state0, sizeof(double) * ckt->CKTnumStates);

        ckt->CKTstat->STATnumIter += k;

        deallocate_worker_io(w, N + 1);

        return OK;
    } else {

        deallocate_worker_io(w, N + 1);
        return(E_ITERLIM);
    }



    
    deallocate_worker_io(w, N + 1);
    return OK;

}



int fake_near(CKTcircuit *ckt, int maxIter, int N /* number of workers*/)
{
    double *OldCKTstate0 = NULL;
    int error, i, j, k;

    int iterno = 0;

    SPICE_debug(("entering...ckt->CTKmode=0x%08x\n", (uint32_t)(ckt->CKTmode)));

    /* make a copy of CTKrhsOld[] */
    int size = SMPmatSize(ckt->CKTmatrix);
    double *x0 = (double*)malloc(sizeof(double) * (size + 1));
    assert(x0);
    memcpy(x0, ckt->CKTrhsOld, sizeof(double) * (size + 1));

    for (i = 0; i < size + 1; ++ i) {
        SPICE_debug(("rhs, rhsOld, spare[%d]=%f, %f, %f\n", i, ckt->CKTrhs[i], ckt->CKTrhsOld[i], ckt->CKTrhsSpare[i]));
    }

    /*
     * it seems that the CKTmode & SMPxxx() calls change in the following way (at least for OP):
     * 
     * =========================================================================
     *  iterations  |    CKTmode            | SMPxxx() calls
     * =========================================================================
     *  1st         | 0x210 (INITJCT)       | preorder + reorder + solve
     * -------------+-----------------------+-----------------------------------
     *  rest        | 0x410 (INITFIX)       |
     * -------------+-----------------------| luFac + solve
     *  last        | 0x110 (INITFLOAT)     | 
     * =========================================================================
     *
     * NB: 
     * - CKTmode affects how CKTload() works (fixme: how?)
     * - SMPxxx() calls solve the matrix
     */


    /*
     * the 1st iteration
     */

    for (i = 0; i < size + 1; ++ i) {
        ckt->CKTrhs[i] = 0.0;
        ckt->CKTrhsOld[i] = 0.0;
        ckt->CKTrhsSpare[i] = 0.0;
    }

    SPICE_debug(("before CKTload(), CKTmode=0x%08x:\n", ckt->CKTmode));
    for (i = 0; i < size + 1; ++ i) {
        SPICE_debug(("rhs, rhsOld, spare[%d]=%f, %f, %f\n", i, ckt->CKTrhs[i], ckt->CKTrhsOld[i], ckt->CKTrhsSpare[i]));
    }

    // load
    error = CKTload(ckt);
    if (error) {
        ckt->CKTstat->STATnumIter += iterno;
    
        SPICE_debug(("load returned error \n"));
    
        FREE(OldCKTstate0);
        return (error);
    }

    SPICE_debug(("after CKTload(), CKTmode=0x%08x:\n", ckt->CKTmode));
    for (i = 0; i < size + 1; ++ i) {
        SPICE_debug(("rhs, rhsOld, spare[%d]=%f, %f, %f\n", i, ckt->CKTrhs[i], ckt->CKTrhsOld[i], ckt->CKTrhsSpare[i]));
    }

    // preorder: call only once
    error = SMPpreOrder(ckt->CKTmatrix);
    if (error) {
        ckt->CKTstat->STATnumIter += iterno;
    
        SPICE_debug(("pre-order returned error \n"));
    
        FREE(OldCKTstate0);
        return(error); /* badly formed matrix */
    }

    // reorder: spOrderAndFactor()
    error = SMPreorder(ckt->CKTmatrix, ckt->CKTpivotAbsTol, ckt->CKTpivotRelTol, ckt->CKTdiagGmin);
    if (error) {
        /* new feature - we can now find out something about what is
         * wrong - so we ask for the troublesome entry
         */
        SMPgetError(ckt->CKTmatrix, &i, &j);
        SPfrontEnd->IFerrorf (ERR_WARNING, "singular matrix:  check nodes %s and %s\n", NODENAME(ckt, i), NODENAME(ckt, j));
        ckt->CKTstat->STATnumIter += iterno;
    
        SPICE_debug(("reorder returned error \n"));
    
        FREE(OldCKTstate0);
        return(error); /* can't handle these errors - pass up! */
    }

    // solve
    SMPsolve(ckt->CKTmatrix, ckt->CKTrhs, ckt->CKTrhsSpare);
    SPICE_debug(("after solve:\n"));
    for (i = 0; i < size + 1; ++ i) {
        SPICE_debug(("rhs, rhsOld, spare[%d]=%f, %f, %f; noncon=%d\n", i, ckt->CKTrhs[i], ckt->CKTrhsOld[i], ckt->CKTrhsSpare[i], ckt->CKTnoncon));
    }

    if (ckt->CKTmode & MODEINITJCT) {
        ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFIX;
        ckt->CKTniState |= NISHOULDREORDER;
    }

    SWAP(double *, ckt->CKTrhs, ckt->CKTrhsOld);

    for (k = 1; k < maxIter; k ++) {


        SPICE_debug(("######################>> iterno=%d: noncon=%d, CKTmode=0x%08x\n", k, ckt->CKTnoncon, (uint32_t)(ckt->CKTmode)));

        ckt->CKTnoncon = 0;

        SPICE_debug(("%d: before CKTload(), CKTmode=0x%08x:\n", k, ckt->CKTmode));
        for (i = 0; i < size + 1; ++ i) {
            SPICE_debug(("rhs, rhsOld, spare[%d]=%f, %f, %f\n", i, ckt->CKTrhs[i], ckt->CKTrhsOld[i], ckt->CKTrhsSpare[i]));
        }
        
        CKTload(ckt);
    
        SPICE_debug(("%d: after CKTload(), CKTmode=0x%08x:\n", k, ckt->CKTmode));
        for (i = 0; i < size + 1; ++ i) {
            SPICE_debug(("rhs, rhsOld, spare[%d]=%f, %f, %f\n", i, ckt->CKTrhs[i], ckt->CKTrhsOld[i], ckt->CKTrhsSpare[i]));
        }
    
        SMPluFac(ckt->CKTmatrix, ckt->CKTpivotAbsTol, ckt->CKTdiagGmin);
        SMPsolve(ckt->CKTmatrix, ckt->CKTrhs, ckt->CKTrhsSpare);
        
        SPICE_debug(("%d: after solve, noncon=%d\n", k, ckt->CKTnoncon));
        for (i = 0; i < size + 1; ++ i) {
            SPICE_debug(("rhs, rhsOld, spare[%d]=%f, %f, %f\n", i, ckt->CKTrhs[i], ckt->CKTrhsOld[i], ckt->CKTrhsSpare[i]));
        }


        // check convergence
        if (ckt->CKTnoncon == 0) {
            ckt->CKTnoncon = NIconvTest(ckt);
            SPICE_debug(("NIconvTest=%d\n", ckt->CKTnoncon));
        }

        if (ckt->CKTnoncon == 0) {
            if (ckt->CKTmode & MODEINITFIX) {
                ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFLOAT;
                SPICE_debug(("ready for the last iteration. CKTmode=0x%08x\n", ckt->CKTmode));
            } else if(ckt->CKTmode & MODEINITFLOAT) {
                printf("EA: CKTmode=0x%08x, Converged!, iterno=%d\n", ckt->CKTmode, k + 1);
                ckt->CKTstat->STATnumIter += k;
                return (OK);
            }
        }
        
        SWAP(double *, ckt->CKTrhs, ckt->CKTrhsOld);
        
    }

    return(E_ITERLIM);

}

