//===========================================================================
//                                                                           
// File: test_matrix.C                                                       
//                                                                           
// Created: Wed Mar 10 10:41:26 2004                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: test_matrix.C,v 1.2 2006-04-06 14:10:41 afr Exp $
//                                                                           
//===========================================================================

#include "GoTools/utils/MatrixXD.h"

using namespace std;
using namespace Go;

template <int N>
void printAllSubs(const MatrixXD<double, N>& m)
{
    cout << m << "det = " << m.det() << "\n\n";

    for (int i = 0; i < N; ++i) {
	for (int j = 0; j < N; ++j) {
	    MatrixXD<double, N-1> sm = m.submatrix(i,j);
	    printAllSubs(sm);
	}
    }
}

template <>
void printAllSubs(const MatrixXD<double, 2>& m)
{
    cout << m << "det = " << m.det() << "\n\n";
}

int main()
{
    const int N = 5;
    int n = 1;
    MatrixXD<double, N> m;
    for (int i = 0; i < N; ++i) {
	for (int j = 0; j < N; ++j) {
	    m(i,j) = n++;
	}
    }

    printAllSubs(m);
}
