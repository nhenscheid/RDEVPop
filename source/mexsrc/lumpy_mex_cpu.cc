// ---------------------------------------------------------------------------
// **********************OAIQSIM SIMULATION TOOLBOX***************************
// file: lumpy_mex_cpu.cc
// purpose: This is a mex CPU method to evaluate a lumpy background with specified lump centers 
//
// Author:  Nick Henscheid
// Date:    3-2019
// Contact: nhenscheid@math.arizona.edu
// This software is in the public domain, furnished "as is", without technical
// support, and with no warranty, express or implied, as to its usefulness for
// any purpose.
// ---------------------------------------------------------------------------
#include "mex.h"
#include "lumpy_bgnd_cpu.h"

using namespace CPU;
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{  
    int dim          = (int)mxGetScalar(prhs[0]);
	int K            = (int)mxGetScalar(prhs[1]); 
	int nEval        = (int)mxGetScalar(prhs[2]);  
    bool unif        = (bool)mxGetScalar(prhs[3]);   // Uniform cov or nonuniform cov
    numtype* centers = (numtype*)mxGetData(prhs[4]);
	numtype* evalpts = (numtype*)mxGetData(prhs[5]); 
	numtype B0       = (numtype)mxGetScalar(prhs[6]);
	numtype* b0      = (numtype*)mxGetData(prhs[7]);
	numtype* cov     = (numtype*)mxGetData(prhs[8]);
    if(nlhs>0)
    	plhs[0]     = mxCreateNumericMatrix(nEval,1,mxDOUBLE_CLASS,mxREAL);
	numtype* result  = (numtype*)mxGetData(plhs[0]);

    //Create LumpyBgnd object and evaluate the result
    LumpyBgnd<numtype> L(unif,dim,K,centers,B0,b0,cov);
	L.eval(nEval,evalpts,result);
}

