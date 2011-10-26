//===========================================================================
//                                                                           
// File: perftest.C                                                          
//                                                                           
// Created: Thu May 23 10:09:35 2002                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: perftest.C,v 1.3 2006-02-15 12:38:38 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/utils/timeutils.h"
#include "GoTools/utils/errormacros.h"
#include <iostream>
#include <algorithm>


using namespace Go;
using namespace std;


template <int NN>
struct Count { enum { N = NN }; };

template <>
struct Count<0> { enum { N = 0}; };

template <int N, typename FI, typename OI>
inline void copy_n_ct(FI from, Count<N>, OI to)
{
    *to = *from;
    copy_n_ct(++from, Count<N-1>(), ++to);
}

template <typename FI, typename OI>
inline void copy_n_ct(FI, Count<0>, OI)
{}


template <int N>
inline void copyint_n_ct(int* from, int* to)
{
    *to++ = *from++;
    copyint_n_ct<N-1>(from, to);
}

template <>
inline void copyint_n_ct<0> (int*, int*)
{}


int main(int argc, char** argv)
{
    if (argc < 2) {
	cout << "Usage: " << argv[0] << " number_of_iterations\n\n" << flush;
	return 1;
    }
    const int testsize = 8;

    int array1[testsize];
    int array2[testsize];

    int numiter = atoi(argv[1]);
    double t0 = getCurrentTime();
    for (int i = 0; i < numiter; ++i) {
	copy_n_ct(array1, Count<testsize>(), array2);
    }
    double t1 = getCurrentTime();
    for (int i = 0; i < numiter; ++i) {
  	copyint_n_ct<testsize>(array1, array2);
    }
    double t2 = getCurrentTime();
    for (int i = 0; i < numiter; ++i) {
	copy(array1, array1 + testsize, array2);
    }
    double t3 = getCurrentTime();
    for (int i = 0; i < numiter; ++i) {
	for (int j = 0; j < testsize; ++j) {
	    array2[j] = array1[j];
	}
    }
    double t4 = getCurrentTime();
    for (int i = 0; i < numiter; ++i) {
	array2[0] = array1[0];
	array2[1] = array1[1];
	array2[2] = array1[2];
	array2[3] = array1[3];
	array2[4] = array1[4];
	array2[5] = array1[5];
	array2[6] = array1[6];
	array2[7] = array1[7];
    }
    double* a1 = reinterpret_cast<double*>(array1);
    double* a2 = reinterpret_cast<double*>(array2);
    const int ts2 = testsize/2;
    double t5 = getCurrentTime();
    for (int i = 0; i < numiter; ++i) {
	for (int j = 0; j < ts2; ++j) {
	    a2[j] = a1[j];
	}
    }
    double t6 = getCurrentTime();

    double d1 = (t1-t0)/double(numiter);
    double d2 = (t2-t1)/double(numiter);
    double d3 = (t3-t2)/double(numiter);
    double d4 = (t4-t3)/double(numiter);
    double d5 = (t5-t4)/double(numiter);
    double d6 = (t6-t5)/double(numiter);

    cout << d1 << "   " << d2 << "   " << d3 << "   " << d4 << "   " << d5 << "   " << d6 <<endl;
}
