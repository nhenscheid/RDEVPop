// ************************ OAIQSIM TOOLBOX ******************************
// File:    util_cpu.h
// Purpose: Contains utility functions for CPU code
// Author:  Nick Henscheid
// Date:    4-2017, 11-2018
// Contact: nhenscheid@math.arizona.edu
// References: 
// This software is in the public domain, furnished "as is", without 
// technical support, and with no warranty, express or implied, as to its 
// usefulness for any purpose.
// ------------------------------------------------------------------------
#ifndef __UTIL_H__
#define __UTIL_H__

#include<iostream>
#include<fstream>
#include<stdio.h>
#include<string>    // For std::string
#include<string.h>  // For strcmp, etc.
#include<math.h>
#include<float.h>
#include<stdlib.h>
#include<complex>   // std::complex (superseded by thrust/complex on GPU)
#include<vector>    // std::vector<T>
#include<algorithm> // for std::max etc

// Namespace is to avoid typedef conflicts with GPU code.
namespace CPU{
    typedef double numtype;  // Should probably put this in a more central location
    typedef std::complex<numtype> complexnumtype;
    const   complexnumtype I(0.0,1.0);   // imaginary unit.  Allows for "1.0+1.0*I" notation.
    
    #define PI   (numtype)3.14159265358979323846
    #define CVAC (numtype)299792458   //Speed of light in vacuum
}

//Macro for "verbose mode"  //WHERE IS THIS USED...?
#define VERB(call) if(VMODE){call;}

//Macro for cublas style array indexing (column major)
#define IDX2C(i,j,ld)((j)*(ld)+(i))

                                   
// ***************************************
// *************  File I/O ***************
// ***************************************
// Get number of lines in a txt file
inline int GetNumLines(const char* fname){
	int N = 0;
	std::string line;
	std::ifstream file(fname);
	if(!file){ 
		std::cout<<"Cannot open file.\n";
		return 0;
	}
	while(std::getline(file,line)){
		++N;
	}
	return N;
}

// Load file into array of T's
template<typename T>
inline void LoadFile(T* x, int N, const char* fname) {
	int i;
	std::ifstream in(fname);

	if (!in) {
    	std::cout << "Cannot open file.\n";
    	return;
	}

	for(i=0;i<N;++i){
		in >> x[i];
	}
	
  	in.close();
}

// Save array of T's to a file
template<typename T>
inline void SaveFile(T* x, int N, const char* fname) {
	int i;
    int precision = 16; // Double precision saving
	std::ofstream out;
	out.open(fname);
	out.precision(precision);
	if (!out) {
    	std::cout << "Cannot open file.\n";
    	return;
	}

	for(i=0;i<N;++i){
		out << x[i];
		out << "\n";
	}
	
  	out.close();
}

// Save N arrays of T's to a csv
// The array is assumed to be unrolled i.e. x = [x1;x2;...;xM]
// Column headings are contained in the std::vector<string> cnames
// ASSUMING THAT EACH COLUMN HAS SAME NUMBER OF ROWS

template<typename T>
inline void SaveCSV(T* x,int M, int N, const char* fname, std::vector<std::string> cnames) {
	int i,j;
    int precision = 16; // Double precision saving
	std::ofstream out;
	out.open(fname);
	out.precision(precision);
	if (!out) {
    	std::cout << "Cannot open file.\n";
    	return;
	}

    for(int icol = 0;icol<N-1;++icol){
        out << cnames[icol] << ",";
    }
    out<< cnames[N-1] <<"\n";
	for(i=0;i<M;++i){ //rows
        for(j=0;j<N-1;++j){ //cols
    		out << x[j*M + i]<<",";
        }
        out << x[(N-1)*M + i] << "\n";
	}
	
  	out.close();
}

// ***************************************
// ***********  Misc. Math Ops ***********
// ***************************************
//Compute the p norm of an array. T should be a floating point type, e.g. float or double.
//This should get replaced with something better.
template<typename T>
inline T PNorm(T* x,int N,T p){
	T result = 0;
	for(int i=0;i<N;++i){
		result += pow(abs(x[i]),p);
	}
	result = pow(result,1/p);
	return result;
}
//Compute p distance of two arrays 
template<typename T>
inline T PDist(T* x,T* y,int N,T p){
	T result = 0;
	for(int i=0;i<N;++i){
		result += pow(abs(x[i]-y[i]),p);
	}
	result = pow(result,1/p);
	return result;
}


#endif /* !UTIL_CUH */

