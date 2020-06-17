#ifndef __UTIL_GPU_H__
#define __UTIL_GPU_H__

// Standard includes
#include<iostream>
#include<fstream>
#include<vector>
#include<stdio.h>
#include<cuda.h>
#include<string>
#include<math.h>
#include<float.h>
#include<stdlib.h>
#include<sys/time.h>
// Thrust
#include <thrust/complex.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

//Include util_cpu to get all the non-gpu-specific utilities
#include<util_cpu.h>  

// Namespace is to avoid typedef conflicts with CPU code.
namespace GPU{
    typedef double numtype;
    typedef thrust::complex<numtype> complexnumtype;
    const   thrust::complex<numtype> I(0.0,1.0);   // imaginary unit.  Allows for "1.0+1.0*I" notation.
// **********************************//
    //PI and CVAC are defined in util_cpu!
//#define PI   (numtype)M_PI
//#define CVAC (numtype)299792458   //Speed of light in vacuum

} // End namespace    
    
//Shorthand notation for device/host functions
#define __dh__ __device__ __host__
#define __d__ __device__


//Macro for "verbose mode"
#define VERB(call) if(VMODE){call;}

//Macro for cublas style array indexing (column major)
#define IDX2C(i,j,ld)((j)*(ld)+(i))


// *****************************************
// ************* CUDA Macros ***************
// *****************************************
//Wrapper for all CUDA API function calls
#define CHECK(call)															\
{																			\
	const cudaError_t error = call;											\
	if(error != cudaSuccess)												\
	{																		\
		printf("Error: %s:%d, ",__FILE__,__LINE__);							\
		printf("code:%d, reason: %s\n", error, cudaGetErrorString(error));	\
		exit(1);															\
	}																		\
}	

#define CUBLASCHECK(call)													\
{																			\
	const cublasStatus_t status = call;										\
	if (status != CUBLAS_STATUS_SUCCESS)									\
    {																		\
        printf("cuBLAS Error: %s:%d, \n",__FILE__,__LINE__);	            \
		printf("code:%d \n", status);                                       \
        exit(1);															\
    }																		\
}																			\

#define MGPUSYNC(NGPU)                                                      \
{                                                                           \
	for(int igpu=0;igpu<NGPU;++igpu)                                        \
	{                                                                       \
		CHECK(cudaSetDevice(igpu));                                         \
		CHECK(cudaDeviceSynchronize());                                     \
	}                                                                       \
}                                                                           \

// Get clock time
double cpuSecond(){
    struct timeval tp;
    gettimeofday(&tp,NULL);
    return ((double)tp.tv_sec + (double)tp.tv_usec*1e-6);
}
        
inline void MemStats(void){
	int numGPU;
	CHECK(cudaGetDeviceCount(&numGPU));
	std::cout<<"You have "<<numGPU<<" GPU devices available."<<std::endl;
	size_t freeMem;
	size_t totalMem;
	for(int i=0;i<numGPU;++i){
		CHECK(cudaSetDevice(i));
		CHECK(cudaMemGetInfo(&freeMem,&totalMem));
		std::cout<<"For device "<< i << ", total mem: "<<totalMem<<", free mem:"<<freeMem<<std::endl;
	}
}

class Managed 
{
public:
  void *operator new(size_t len) {
    void *ptr;
    CHECK(cudaMallocManaged(&ptr, len));
    CHECK(cudaDeviceSynchronize());
    return ptr;
  }

  void operator delete(void *ptr) {
    CHECK(cudaDeviceSynchronize());
    CHECK(cudaFree(ptr));
  }
};

// Misc GPU utility functions


// Initialize an array of T to a scalar value.
namespace GPU{
    template<typename T>
    __global__ void initArray(int n,T* x,T a)
    {
        int idx = blockIdx.x*blockDim.x + threadIdx.x;
        if(idx<n)
            x[idx] = a;
    }
}


// Misc Math Functions




// NEED TO CONVERT TO GPU CODE
// //Compute the p norm of an array. T should be a floating point type, e.g. float or double.
// //This should get replaced with something better.
// template<typename T>
// inline T PNorm(T* x,int N,T p){
// 	T result = 0;
// 	for(int i=0;i<N;++i){
// 		result += pow(abs(x[i]),p);
// 	}
// 	result = pow(result,1/p);
// 	return result;
// }
// //Compute p distance of two arrays 
// template<typename T>
// inline T PDist(T* x,T* y,int N,T p){
// 	T result = 0;
// 	for(int i=0;i<N;++i){
// 		result += pow(abs(x[i]-y[i]),p);
// 	}
// 	result = pow(result,1/p);
// 	return result;
// }



#endif /* !UTIL_CUH */

