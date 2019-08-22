/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
Modified: Apr 2000 - Paolo Nenzi
**********/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "resdefs.h"
#include "ngspice/sperror.h"


/* actually load the current resistance value into the sparse matrix
 * previously provided */
int
RESload(GENmodel *inModel, CKTcircuit *ckt)
{

#ifdef STEPDEBUG
    SPICE_debug(("entering...\n"));
#endif

    RESmodel *model = (RESmodel *)inModel;

    /*  loop through all the resistor models */
    for( ; model != NULL; model = RESnextModel(model)) {
#ifdef STEPDEBUG
        SPICE_debug(("RESmodel: %s\n", model->gen.GENmodName));
#endif
        RESinstance *here;

        /* loop through all the instances of the model */
        for (here = RESinstances(model); here != NULL ;
                here = RESnextInstance(here)) {

#ifdef STEPDEBUG
            SPICE_debug(("  RESinstances: %s %d %d %.2f\n", here->gen.GENname, here->RESposNode, here->RESnegNode, here->RESresist));
#endif
            // CKTrhsOld are the right hand side vector from the previous iteration which 
            // now contains the solution (nodal voltage) for that iteration:
            // i = v * g
            here->REScurrent = (*(ckt->CKTrhsOld+here->RESposNode) -
                                *(ckt->CKTrhsOld+here->RESnegNode)) * here->RESconduct;

            // element stamp
            *(here->RESposPosPtr) += here->RESconduct;
            *(here->RESnegNegPtr) += here->RESconduct;
            *(here->RESposNegPtr) -= here->RESconduct;
            *(here->RESnegPosPtr) -= here->RESconduct;
        }
    }
    return(OK);
}


/* actually load the current resistance value into the sparse matrix
 * previously provided */
int
RESacload(GENmodel *inModel, CKTcircuit *ckt)
{
    RESmodel *model = (RESmodel *)inModel;
    double g;

    NG_IGNORE(ckt);

    /*  loop through all the resistor models */
    for( ; model != NULL; model = RESnextModel(model)) {
        RESinstance *here;

        /* loop through all the instances of the model */
        for (here = RESinstances(model); here != NULL ;
             here = RESnextInstance(here)) {

            if (here->RESacresGiven)
                g = here->RESacConduct;
            else
                g = here->RESconduct;

            *(here->RESposPosPtr) += g;
            *(here->RESnegNegPtr) += g;
            *(here->RESposNegPtr) -= g;
            *(here->RESnegPosPtr) -= g;
        }
    }
    return(OK);
}
