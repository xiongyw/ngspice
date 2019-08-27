/*
 * Copyright (c) 2014, NVIDIA Corporation. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, 
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, 
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, 
 *    this list of conditions and the following disclaimer in the documentation and/or 
 *    other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to 
 *    endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, 
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, 
 * OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "ngspice/ngspice.h"
#include "ngspice/config.h"
#include "resdefs.h"

extern "C"
__global__ void cuRESload_kernel (RESparamGPUstruct, double *, int, int *, double *) ;

extern "C"
int
cuRESload
(
GENmodel *inModel, CKTcircuit *ckt
)
{
    RESmodel *model = (RESmodel *)inModel ;
    int thread_x, thread_y, block_x ;

    cudaError_t status ;

#ifdef STEPDEBUG
    SPICE_debug(("entering...\n"));
#endif

    /*  loop through all the resistor models */
    for ( ; model != NULL ; model = model->RESnextModel)
    {
#ifdef STEPDEBUG
        SPICE_debug(("  RESmodName=%s\n", model->RESmodName));
#endif
        /* Determining how many blocks should exist in the kernel */
        thread_x = 1 ;
        thread_y = 256 ;
        if (model->n_instances % thread_y != 0)
            block_x = (int)((model->n_instances + thread_y - 1) / thread_y) ;
        else
            block_x = model->n_instances / thread_y ;

        dim3 thread (thread_x, thread_y) ;

        /* Kernel launch */
        status = cudaGetLastError () ; // clear error status
        
#ifdef STEPDEBUG
        SPICE_debug(("  calling cuRESload_kernel() for model %s: block=%d, thread=(%d,%d)\n", model->RESmodName, block_x, thread_x, thread_y));
#endif
        cuRESload_kernel <<< block_x, thread >>> (model->RESparamGPU, ckt->d_CKTrhsOld, model->n_instances,
                                                  model->d_PositionVector, ckt->d_CKTloadOutput) ;

        cudaDeviceSynchronize () ;

        status = cudaGetLastError () ; // check for launch error
        if (status != cudaSuccess)
        {
            fprintf (stderr, "Kernel launch failure in the Resistor Model\n\n") ;
            return (E_NOMEM) ;
        }
    }

    return (OK) ;
}

extern "C"
__global__
void
cuRESload_kernel
(
RESparamGPUstruct RESentry, double *CKTrhsOld, int n_instances, int *d_PositionVector, double * d_CKTloadOutput
)
{
    double m, difference, factor ;

    int instance_ID ;

    instance_ID = threadIdx.y + blockDim.y * blockIdx.x ;
    if (instance_ID < n_instances)
    {
        if (threadIdx.x == 0)
        {
            if (!(RESentry.d_REStc1GivenArray [instance_ID]))
                RESentry.d_REStc1Array [instance_ID] = 0.0 ;
            
            if (!(RESentry.d_REStc2GivenArray [instance_ID]))
                RESentry.d_REStc2Array [instance_ID] = 0.0 ;
            
            if (!(RESentry.d_RESmGivenArray [instance_ID]))
                RESentry.d_RESmArray [instance_ID] = 1.0 ;

            RESentry.d_REScurrentArray [instance_ID] = (CKTrhsOld [RESentry.d_RESposNodeArray [instance_ID]] -
                                                    CKTrhsOld [RESentry.d_RESnegNodeArray [instance_ID]]) *
                                                    RESentry.d_RESconductArray [instance_ID] ;
            
            difference = (RESentry.d_REStempArray [instance_ID] + RESentry.d_RESdtempArray [instance_ID]) - 300.15 ;
            factor = 1.0 + (RESentry.d_REStc1Array [instance_ID]) * difference +
                     (RESentry.d_REStc2Array [instance_ID]) * difference * difference ;
            
            m = (RESentry.d_RESmArray [instance_ID]) / factor ;
            
            d_CKTloadOutput [d_PositionVector [instance_ID]] = m * RESentry.d_RESconductArray [instance_ID] ;
        }
    }

    return ;
}
