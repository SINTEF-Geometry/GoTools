//===========================================================================
//                                                                           
// File: compmeancurv.C                                                      
//                                                                           
// Created: Thu Feb 19 11:10:48 2004                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: compmeancurv.C,v 1.1 2004-03-26 08:34:11 afr Exp $
//                                                                           
//===========================================================================

#include <iostream>
#include <vector>

using namespace std;


// Taking second order differences
void comp_finite_diff(int M, int N, int dir,
		      const vector<double>& pts,
		      vector<double>& diffs)
{
    double h;
    int off;
    if (dir == 0) {
	// U-direction
	h = 1.0/(M-1.0);
	off = 1;
    } else {
	// V-direction
	h = 1.0/(N-1.0);
	off = M;
    }

    diffs.clear();
    diffs.resize(M*N, 0.0);
    for (int i = 1; i < M-1; ++i) {
	for (int j = 1; j < N-1; ++j) {
	    diffs[i*M + j] = pts[i*M + j + off] - pts[i*M + j - off];
	    diffs[i*M + j] /= 2*h;
	}
    }
    
}


int main()
{
    // Read the grid
    int M, N;
    cin >> M >> N;
    vector<double> f(M*N);
    for (int i = 0; i < M*N; ++i) {
	cin >> f[i];
    }

    vector<double> fu, fv, fuu, fuv, fvv, sum;

    comp_finite_diff(M, N, 0, f, fu);
    comp_finite_diff(M, N, 1, f, fv);
    comp_finite_diff(M, N, 0, fu, fuu);
    comp_finite_diff(M, N, 1, fu, fuv);
    comp_finite_diff(M, N, 1, fv, fvv);

    sum.resize(M*N);
    for (int i = 0; i < M*N; ++i) {
	sum[i] = (1+fv[i]*fv[i])*fuu[i]
	    + 2*fu[i]*fv[i]*fuv[i]
	    + (1+fu[i]*fu[i])*fvv[i];
	int c = i%M;
	int r = i/M;
	if ((r > 1) && (r < M-2) && (c > 1) && (c < N-2)) {
	    cout << sum[i] << " ";
	}
	if ((i+1)%M == 0) cout << endl;
    }
}
