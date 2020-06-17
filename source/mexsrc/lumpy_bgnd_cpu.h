// ------------------------------------------------------------------------
// *************************** OAIQSIM TOOLBOX ****************************
// File:    lumpy_bgnd_cpu.h
// Purpose: This header contains the CPU class definition and compute 
//          kernel for the standard lumpy background random field.
// Author:  Nick Henscheid
// Date:    9-2016, 2-2019, 4-2019
// Contact: nph@email.arizona.edu
// References: "Factors Influencing Lesion Detection in Medical Imaging" 
//              J.P. Rolland, 1990 and 
//             "Foundations of Image Science", H.H. Barrett and K.J. Myers 
// This software is in the public domain, furnished "as is", without 
// technical support, and with no warranty, express or implied, as to its 
// usefulness for any purpose.
// ------------------------------------------------------------------------
#ifndef __OAIQSIM_OBJECTS_LUMPY_KERNEL_CPU_H__
#define __OAIQSIM_OBJECTS_LUMPY_KERNEL_CPU_H__
#include <math.h>
#include <util_cpu.h>

namespace CPU{

// Forward declarations
template<typename T> 
T LumpyKernelUniform(T*,T*,int,int,T,T*,T*);

template<typename T> 
T LumpyKernelNonuniform(T*,T*,int,int,T,T*,T*);

template<typename T>
T GaussianLump2D(T*,T*,T,T*);

template<typename T>
T GaussianLump3D(T*,T*,T,T*);

template<typename T>
bool CheckInput(T*,int,int);

// LumpyBgnd object 
template<typename T>
class LumpyBgnd{
    public:
        // Default constructor
        LumpyBgnd(){
            unif = 1; //Uniform P
            dim = 2; K = 1;
            B0 = 0; 
            b0 = (T*)malloc(sizeof(T)); b0[0] = 1.0;
            P = (T*)malloc(3*sizeof(T)); P[0] = 1.0;P[1] = 0.0;P[2] = 1.0;
            centers = (T*)malloc(2*sizeof(T)); centers[0] = 0.5;centers[1] = 0.5;
        };
        // Construct default object of dimension dim
        LumpyBgnd(int dim_i){
            if(dim_i>0){
                unif = 1;
                dim = dim_i; K = 1;
                B0 = 0; 
                b0 = (T*)malloc(sizeof(T)); b0[0] = 1.0;
                P = (T*)malloc(dim*sizeof(T)); 
                centers = (T*)malloc(dim*sizeof(T));
                for(int i=0;i<dim;++i){
                    P[i] = 1.0;
                    centers[i] = 0.5;
                }                
            }
        };
        LumpyBgnd(bool type,int dim_i,int K_i,T* centers_i,T B0_i,T* b0_i,T* P_i){
            if((dim_i>0)&&(K_i>0)){
                unif = type;
                dim = dim_i; K = K_i;
                centers = centers_i;
                B0 = B0_i; b0 = b0_i; P = P_i;
            }
        }
        
        void set_dim(int d){if(dim==2||dim==3){dim = d;}}
        int  dimension(void){return dim;}
        void set_K(int k){K=k;}
        int  num_lumps(void){return K;}
        void set_B0(numtype B){B0 = B;}
        numtype DC_val(void){return B0;}
        void set_centers(numtype* p){if(p){centers = p;}}// if(p){} checks for null ptr
        void set_b0(numtype* p){if(p){b0 = p;}}
        void set_precision(numtype* p){if(p){P = p;}}
        
        void eval(int nEval,T* evalpts,T* result)
        {
            if(unif==1){
                for(int ix=0;ix<nEval;++ix){
                    result[ix] = LumpyKernelUniform<T>(centers,evalpts+dim*ix,dim,K,B0,b0,P);
                }
            }else{
                for(int ix=0;ix<nEval;++ix){
                    result[ix] = LumpyKernelNonuniform<T>(centers,evalpts+dim*ix,dim,K,B0,b0,P);
                }
            }
        }
        
        T operator()(T* evalpt){
            // Evaluate the texture for a single point
            return LumpyKernel(centers,evalpt,dim,K,B0,b0,P);
        }
        
    private:
        int dim;
        int K;
        T B0;
        T* centers;
        T* b0;
        T* P;
        bool unif; // Uniform vs. non-uniform  
};// LumpyBgnd

template<typename T> 
T LumpyKernelUniform(T* centers,T* eval,int dim,int K,T B0,T* b0,T* P)
{
    // Should check if eval satisfies a domain condition?
    T result   = B0;
    if(dim==2){// Could probably use templates to do this dim switch
        for(int iLump=0;iLump<K;++iLump){ 
            // Should figure out how to limit the loop to only those lumps that
            // are "close" to the evaluation point.  Pseudogrid?
            result += GaussianLump2D(eval,centers+dim*iLump,b0[iLump],P);
        }
    }else if(dim==3){
        for(int iLump=0;iLump<K;++iLump){
            // Should figure out how to limit the loop to only those lumps that
            // are "close" to the evaluation point.  Pseudogrid?
            result += GaussianLump3D(eval,centers+dim*iLump,b0[iLump],P);
        }
    }
    return result;
}

template<typename T> 
T LumpyKernelNonuniform(T* centers,T* eval,int dim,int K,T B0,T* b0,T* P)
{
    // Should check if eval satisfies a domain condition?
    T result   = B0;
    if(dim==2){ // Could probably use templates to do this dim switch
        for(int iLump=0;iLump<K;++iLump){
            // Should figure out how to limit the loop to only those lumps that
            // are "close" to the evaluation point.  Pseudogrid?
            result += GaussianLump2D(eval,centers+dim*iLump,b0[iLump],P + iLump*dim*(dim+1)/2);
        }
    }else if(dim==3){
        for(int iLump=0;iLump<K;++iLump){
            // Should figure out how to limit the loop to only those lumps that
            // are "close" to the evaluation point.  Pseudogrid?
            result += GaussianLump3D(eval,centers+dim*iLump,b0[iLump],P + iLump*dim*(dim+1)/2);
        }
    }
 
    return result;
}

template<typename T>
T GaussianLump2D(T* R, T* R0, T b0, T* P)
{
    // Evaluates a single *2D* anisotropic Gaussian Lump with center R0 at 
    // position R with amplitude b0 and precision matrix P
    // P[0],...,theta[d*(d+1)/2-1] gives the precision matrix in symmetric
    // matrix column-major form.  Note that P is the precision matrix, i.e.
    // l(R) = b0*exp(-0.5(R-R0)'*P*(R-R0))
    T xs  = R[0]-R0[0];
    T ys  = R[1]-R0[1];
    T exponent = P[0]*xs*xs + 2.0*P[1]*xs*ys + P[2]*ys*ys;
    return b0*exp(-0.5*exponent);
}

template<typename T>
T GaussianLump3D(T* R, T* R0, T b0, T* P)
{
    // Evaluates a single *3D* anisotropic Gaussian Lump with center R0 at 
    // position R with amplitude b0 and precision matrix P
    // P[0],...,theta[d*(d+1)/2-1] gives the precision matrix in symmetric
    // matrix column-major form.  Note that P is the precision matrix, i.e.
    // l(R) = b0*exp(-0.5(R-R0)'*P*(R-R0))
    T xs  = R[0]-R0[0];
    T ys  = R[1]-R0[1];
    T zs  = R[2]-R0[2];
    
    T exponent = P[0]*xs*xs + P[3]*ys*ys + P[5]*zs*zs 
               + 2.0*P[1]*xs*ys + 2.0*P[2]*xs*zs + 2.0*P[4]*ys*zs;
    return b0*exp(-0.5*exponent);
}


}//namespace CPU
#endif