//===========================================================================
//                                                                           
// File: surfeval2.C                                                          
//                                                                           
// Created: Nov. 10, 2008
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
//===========================================================================


#include <fstream>
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/errormacros.h"

using namespace Go;
using namespace std;

int main(int argc, char* argv[] )
{
    ALWAYS_ERROR_IF(argc < 3, "Usage: " << argv[0]
		    << " inputsurf inputpoints" << endl);

    // Open input surface file
    ifstream is(argv[1]);
    ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

    ObjectHeader head;
    is >> head;

    // Read surface from file
    SplineSurface sf;
    is >> sf;

    // Get points
    ifstream pts(argv[2]);
    ALWAYS_ERROR_IF(pts.bad(), "Bad or no input filename");
    int n;
    pts >> n;
    vector<double> pt(n*2);
    for (int i = 0; i < n; ++i) {
	pts >> pt[2*i] >> pt[2*i+1];
    }

    vector<double>::const_iterator coefs = sf.coefs_begin();
    int kk1 = sf.order_u();
    int kk2 = sf.order_v();
    int nn1 = sf.numCoefs_u();
    int dim = sf.dimension();
    std::vector<Point> p(3, Point(dim));
    std::vector<Point> p2(3, Point(dim));
    int kp, kq, kr, kd;
    for (int j = 0; j < n; ++j) {
	sf.point(p, pt[2*j], pt[2*j+1], 1);
	cout << p[0][0] << " " << p[0][1] << " " << p[0][2] << " " ;
	cout << p[1][0] << " " << p[1][1] << " " << p[1][2] << " " ;
	cout << p[2][0] << " " << p[2][1] << " " << p[2][2] << endl;

	vector<double> b0, bu, bv;
	sf.computeBasis(&pt[2*j], b0, bu, bv);
	int left1 = sf.basis_u().knotInterval(pt[2*j]) - kk1 + 1;
	int left2 = sf.basis_v().knotInterval(pt[2*j+1]) - kk2 + 1;
	
	p2[0].setValue(0.0);
	p2[1].setValue(0.0);
	p2[2].setValue(0.0);
	for (kq=left2, kr=0; kq<left2+kk2; kq++)
	    for (kp=left1; kp<left1+kk1; kp++, kr++)
		for (kd=0; kd<dim; kd++)
		{
		    double val = coefs[(kq*nn1+kp)*dim+kd];
		    p2[0].begin()[kd] += val*b0[kr];
		    p2[1].begin()[kd] += val*bu[kr];
		    p2[2].begin()[kd] += val*bv[kr];
		}
	cout << p2[0][0] << " " << p2[0][1] << " " << p2[0][2] << " " ;
	cout << p2[1][0] << " " << p2[1][1] << " " << p2[1][2] << " " ;
	cout << p2[2][0] << " " << p2[2][1] << " " << p2[2][2] << endl;
	cout << endl;
    }
}





