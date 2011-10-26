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
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/errormacros.h"

using namespace Go;
using namespace std;

int main(int argc, char* argv[] )
{
    ALWAYS_ERROR_IF(argc < 4, "Usage: " << argv[0]
		    << " inputvol inputpoints show_derivs?" << endl);

    // Open input surface file
    ifstream is(argv[1]);
    ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

    ObjectHeader head;
    is >> head;

    // Read surface from file
    SplineVolume vol;
    is >> vol;

    // Get points
    ifstream pts(argv[2]);
    ALWAYS_ERROR_IF(pts.bad(), "Bad or no input filename");
    int n;
    pts >> n;
    vector<double> pt(n*3);
    for (int i = 0; i < n; ++i) {
	pts >> pt[3*i] >> pt[3*i+1] >> pt[3*i+2];
    }

    bool showDerivs = (atoi(argv[3]) != 0);

    vector<double>::const_iterator coefs = vol.coefs_begin();
    int kk1 = vol.order(0);
    int kk2 = vol.order(1);
    int kk3 = vol.order(2);
    int nn1 = vol.numCoefs(0);
    int nn2 = vol.numCoefs(1);
    int dim = vol.dimension();
    std::vector<Point> p(4, Point(dim));
    //Point p(dim);
    std::vector<Point> p2(4, Point(dim));
    int kp, kq, kh, kr, kd;
    for (int j = 0; j < n; ++j) {
 	vol.point(p, pt[3*j], pt[3*j+1], pt[3*j+2], 1);
 	cout << p[0][0] << " " << p[0][1] << " " << p[0][2] << " " ;
 	cout << p[1][0] << " " << p[1][1] << " " << p[1][2] << " " ;
	cout << p[2][0] << " " << p[2][1] << " " << p[2][2] << " " ;
 	cout << p[3][0] << " " << p[3][1] << " " << p[3][2] << endl;

	vector<double> b0, bu, bv, bw;
	vol.computeBasis(&pt[3*j], b0, bu, bv, bw);
	int left1 = vol.basis(0).knotInterval(pt[3*j]) - kk1 + 1;
	int left2 = vol.basis(1).knotInterval(pt[3*j+1]) - kk2 + 1;
	int left3 = vol.basis(2).knotInterval(pt[3*j+2]) - kk3 + 1;
	
	BasisPts res;
	vol.computeBasis(pt[3*j], pt[3*j+1], pt[3*j+2], res);

	p2[0].setValue(0.0);
	p2[1].setValue(0.0);
	p2[2].setValue(0.0);
	p2[3].setValue(0.0);
	for (kh=left3, kr=0; kh<left3+kk3; kh++)
	    for (kq=left2; kq<left2+kk2; kq++)
		for (kp=left1; kp<left1+kk1; kp++, kr++)
		    for (kd=0; kd<dim; kd++)
		    {
			double val = coefs[((kq+kh*nn2)*nn1+kp)*dim+kd];
			p2[0].begin()[kd] += val*b0[kr];
			p2[1].begin()[kd] += val*bu[kr];
			p2[2].begin()[kd] += val*bv[kr];
			p2[3].begin()[kd] += val*bw[kr];
		    }
	cout << p2[0][0] << " " << p2[0][1] << " " << p2[0][2] << " " ;
	if (showDerivs)
	  {
	    cout << p2[1][0] << " " << p2[1][1] << " " << p2[1][2] << " " ;
	    cout << p2[2][0] << " " << p2[2][1] << " " << p2[2][2] << " " ;
	    cout << p2[3][0] << " " << p2[3][1] << " " << p2[3][2];
	  }
	cout << endl;
	cout << endl;
    }

    ofstream os("tmp_vol.g2");
    vol.writeStandardHeader(os);
    vol.write(os);
}





