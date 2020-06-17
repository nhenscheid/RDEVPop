// ------------------------------------------------------------------------
// *************************** OAIQSIM TOOLBOX ****************************
// File:    lumpy_kernel_gpu.cu
// Purpose: LumpyKernel computes generalized lumpy background functions 
//          u(r) = B0 + sum_{j=1}^K phi(r-r_j;theta_j)
//          The standard example is the "Gaussian" lumpy background where 
//          phi(r) = (b0/pi*rb^2)*exp(-|x|^2/rb^2)
//          The "lump function" phi is defined as a device function 
// Inputs:  dim is the ambient dimension (any dim>0 possible, in theory)
//          K is the number of lumps
//          nEval is the number of evaluation points 
//          nParam is the number of parameters (number of columns in theta)
//          centers is an nEval-by-dim (row-major) array of centers (r_j)
//          evalpts is an nEval-by-dim (row-major) array of eval points (r)
//          result is an nEval-by-1 array to store the result 
//          theta is an K-by-nParam (row-major) array of parameters 
//          (theta_j)
//          The template parameter T is typically either float or double
// References: "Factors Influencing Lesion Detection in Medical Imaging" 
//              J.P. Rolland, 1990 and 
//             "Foundations of Image Science", H.H. Barrett and K.J. Myers 
// Author:  Nick Henscheid
// Date:    9-2016
// Contact: nph@email.arizona.edu
// This software is in the public domain, furnished "as is", without 
// technical support, and with no warranty, express or implied, as to its 
// usefulness for any purpose.
// ------------------------------------------------------------------------
#ifndef __OAIQSIM_OBJECTS_LUMPY_KERNEL_GPU_H__
#define __OAIQSIM_OBJECTS_LUMPY_KERNEL_GPU_H__
#include <math.h>

template<typename T>
__device__ T IsoGaussianLump(int dim, T* R, T* R0, T* theta)
{
    // Evaluates a single (Isotropic) Gaussian Lump with center R0 at 
    // position R with param theta = [b0,rb^2] i.e. theta[0] = b0, 
    // theta[1] = rb^2. b0 and rb^2 are the lump scale, and width 
    // (rb^2 = 2*var) (see documentation) 
    T b0  = theta[0];
    T rb2 = theta[1];
    T exponent = 0;
    
    if((b0!=0)&&(rb2>0)){
        //T A = b0/pow(2.0*M_PI*rb2,((T)dim)/2.0); 
        for(int idim = 0;idim<dim;++idim){
            exponent += (R[idim]-R0[idim])*(R[idim]-R0[idim]);
        }
        return b0*exp(-0.5*exponent/rb2); 
    } else{
        return 0;
    }
}

template<typename T>
__device__ T AnisoGaussianLump2D(T* R, T* R0, T* theta)
{
    // Evaluates a single anisotropic Gaussian Lump with center R0 at 
    // position R with param theta = [b0,P] i.e. theta[0] = b0, 
    // theta[1] = rb^2. b0 and P are the lump scale and precision matrix, respectively 
    // i.e. l(R) \propto b0exp(-0.5(R-R0)'*P*(R-R0))
    T b0  = theta[0];
    T* P  = theta+1; // Column-major symmetric matrix i.e. P = [p11,p21,p22] in 2D
    T exponent = 0;
    //T det = P[0]*P[2]-P[1]*P[1];
    T xs  = R[0]-R0[0];
    T ys  = R[1]-R0[1];
    
    if(b0!=0){
        //T A = b0*sqrt(det)/(2*M_PI);  
        exponent = P[0]*xs*xs + 2.0*P[1]*xs*ys + P[2]*ys*ys;
        return b0*exp(-0.5*exponent); 
    } else{
        return 0;
    }
}

template<typename T>
__device__ T AnisoGaussianLump3D(T* R, T* R0, T* theta)
{
    // Evaluates a single anisotropic Gaussian Lump with center R0 at 
    // position R with param theta = [b0,P] i.e. theta[0] = b0, 
    // theta[1] = rb^2. b0 and P are the lump scale and precision matrix, respectively 
    // i.e. l(R) \propto b0exp(-0.5(R-R0)'*P*(R-R0))
    T b0  = theta[0];
    T* P  = theta+1; // Column-major symmetric matrix i.e. P = [p11,p21,p31,p22,p32,p33] in 3D
    T exponent = 0;
    //T det = P[0]*(P[3]*P[5]-P[4]*P[4]) - P[1]*(P[1]*P[5]-P[4]*P[2]) + P[2]*(P[1]*P[4]-P[3]*P[2]);
    T xs  = R[0]-R0[0];
    T ys  = R[1]-R0[1];
    T zs  = R[2]-R0[2];
    
    if(b0!=0){
        //T A = b0*sqrt(det)/powf(2*M_PI,(T)1.5);  
        exponent = P[0]*xs*xs + P[3]*ys*ys + P[5]*zs*zs 
                 + 2.0*P[1]*xs*ys + 2.0*P[2]*xs*zs + 2.0*P[4]*ys*zs;
        return b0*exp(-0.5*exponent); 
    } else{
        return 0;
    }
}
// *************************************** //
// ******* Isotropic Lumpy Kernel ******** //
// *************************************** //
template<typename T>
__global__ void IsoLumpyKernel(int dim,int K,int nEval,
                            T* centers,T* evalpts,T* result,T B0,T* theta)
{
    unsigned int tId = blockIdx.x*blockDim.x + threadIdx.x; //1d thread structure
    unsigned int idx = tId*dim;
    int nParam = 2;
    
	if(tId<nEval){
		result[tId] = B0;  // Constant DC shift
		for(int iLump=0;iLump<K;++iLump) //THIS SHOULD BE PARALLELIZED (or at least localized)
        {
            result[tId] += IsoGaussianLump(dim,evalpts+idx,centers+iLump,theta+nParam*iLump);
        }
	}
}

// ******************************************** //
// ******* 2D Anisotropic Lumpy Kernel ******** //
// ******************************************** //

template<typename T>
__global__ void AnisoLumpyKernel2D(int K,int nEval,
                                   T* centers,T* evalpts,T* result,T B0,T* theta)
{
    unsigned int tId = blockIdx.x*blockDim.x + threadIdx.x; //1d thread structure
    unsigned int idx = 2*tId;
    int nParam = 4;  // theta = [b0,p11,p21,p22]
    
	if(tId<nEval){
		result[tId] = B0;  // Constant DC shift
		for(int iLump=0;iLump<K;++iLump) //THIS SHOULD BE PARALLELIZED
        {
            result[tId] += AnisoGaussianLump2D(evalpts+idx,centers+2*iLump,theta+nParam*iLump);
        }
	}
}

// ******************************************** //
// ******* 3D Anisotropic Lumpy Kernel ******** //
// ******************************************** //

template<typename T>
__global__ void AnisoLumpyKernel3D(int K,int nEval,
                                   T* centers,T* evalpts,T* result,T B0,T* theta)
{
    unsigned int tId = blockIdx.x*blockDim.x + threadIdx.x; //1d thread structure
    unsigned int idx = 3*tId;
    int nParam = 7; // theta = [b0,p11,p21,p31,p22,p32,p33];
    
	if(tId<nEval){
		result[tId] = B0;  // Constant DC shift
		for(int iLump=0;iLump<K;++iLump) //THIS SHOULD BE PARALLELIZED
        {
            result[tId] += AnisoGaussianLump3D(evalpts+idx,centers+3*iLump,theta+nParam*iLump);
        }
	}
}


#endif