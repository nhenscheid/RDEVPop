// ---------------------------------------------------------------------------
// **********************OAIQSIM SIMULATION TOOLBOX***************************
// file: lumpy_gpu_mex.c
// purpose: This is the mex function interface to the GPU lumpy background computation
//
// Author:  Nick Henscheid
// Date:    9-2016, 3-2019
// Contact: nhenscheid@math.arizona.edu
// This software is in the public domain, furnished "as is", without technical
// support, and with no warranty, express or implied, as to its usefulness for
// any purpose.
// ---------------------------------------------------------------------------

#include "mex.h"
#include "lumpy_kernel_gpu.cuh"
#include "util_gpu.cuh"  // For CHECK() 

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int dim        = (int)mxGetScalar(prhs[0]);  //Dimension (e.g. 2,3)
	int K          = (int)mxGetScalar(prhs[1]);  //Number of lumps
	int nEval      = (int)mxGetScalar(prhs[2]);   //Number of eval points'
    int nParam     = (int)mxGetScalar(prhs[3]);
	float* centers = (float*)mxGetData(prhs[4]);  //Lump centers
	float* evalpts = (float*)mxGetData(prhs[5]);  //Field evaluation points
	float B0       = (float)mxGetScalar(prhs[6]); //Constant/"DC offset"
	float* theta   = (float*)mxGetData(prhs[7]);  //Vector of lump parameters
	
	plhs[0]         = mxCreateNumericMatrix(nEval,1,mxSINGLE_CLASS,mxREAL);
	float* result  = (float*)mxGetData(plhs[0]);
    
    float *centers_d,*evalpts_d,*result_d,*theta_d;
    
	CHECK(cudaMalloc((void**)&(centers_d),dim*K*sizeof(float)));
    CHECK(cudaMalloc((void**)&(evalpts_d),dim*nEval*sizeof(float)));	
    CHECK(cudaMalloc((void**)&(result_d),nEval*sizeof(float)));
    CHECK(cudaMalloc((void**)&(theta_d),K*nParam*sizeof(float)));
    
	CHECK(cudaMemcpy(centers_d,centers,dim*K*sizeof(float),cudaMemcpyHostToDevice));
	CHECK(cudaMemcpy(evalpts_d,evalpts,dim*nEval*sizeof(float),cudaMemcpyHostToDevice));
	CHECK(cudaMemcpy(result_d,result,nEval*sizeof(float),cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(theta_d,theta,K*nParam*sizeof(float),cudaMemcpyHostToDevice));

	// Cuda parameters (1d grid of 1d blocks)
	dim3 block(512,1);
	dim3 grid((nEval+block.x-1)/block.x,1);
	//Lauch kernel
    
    if(nParam == 2){
        // Isotropic 
        IsoLumpyKernel<float><<<grid,block>>>(dim,K,nEval,centers_d,evalpts_d,result_d,B0,theta_d);
    } else if(dim == 2){
        // 2D anisotropic
        AnisoLumpyKernel2D<float><<<grid,block>>>(K,nEval,centers_d,evalpts_d,result_d,B0,theta_d);
    } else if(dim == 3){
        // 3D anisotropic
        AnisoLumpyKernel3D<float><<<grid,block>>>(K,nEval,centers_d,evalpts_d,result_d,B0,theta_d);
    } else {
        mexErrMsgTxt("Error: incorrect settings! (check nParam and/or dim)\n");   
    }

	CHECK(cudaMemcpy(result,result_d,nEval*sizeof(float),cudaMemcpyDeviceToHost));
    CHECK(cudaDeviceReset());
}
