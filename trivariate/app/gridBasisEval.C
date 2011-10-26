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

typedef std::vector<double>  Dvector; // for convenience to shorten
typedef std::vector<Dvector> Dmatrix; // function argument lists...

int main(int argc, char* argv[] )
{
    ALWAYS_ERROR_IF(argc != 5, "Usage: " << " inputfile nmb_u nmb_v nmb_w" << endl);

    // Open input surface file
    ifstream is(argv[1]);
    ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

    ObjectHeader head;
    is >> head;

    // Read volume from file
    SplineVolume vol;
    is >> vol;

    int nmb_u = atoi(argv[2]);
    int nmb_v = atoi(argv[3]);
    int nmb_w = atoi(argv[4]);

    // Evaluate points
    int dim = vol.dimension();
    int size = nmb_u*nmb_v*nmb_w*dim;
    std::vector< double> points(size), der_u(size), der_v(size), der_w(size);
    std::vector<double> par_u(nmb_u), par_v(nmb_v), par_w(nmb_w);
    vol.gridEvaluator(nmb_u, nmb_v, nmb_w, points, der_u, der_v, der_w,
		       par_u, par_v, par_w);

    // Evaluate basis functions
    std::vector<BasisDerivs> result;
    vol.computeBasisGrid(par_u, par_v, par_w, result);
    
    // Altenative evaluation of basis functions
    Dmatrix basisValues, basisDerivs_u, basisDerivs_v, basisDerivs_w;
    vol.computeBasisGrid(par_u, par_v, par_w, basisValues, basisDerivs_u, basisDerivs_v,
			 basisDerivs_w);
    // Write output
    vector<double>::const_iterator coefs = vol.coefs_begin();
    int kk1 = vol.order(0);
    int kk2 = vol.order(1);
    int kk3 = vol.order(2);
    int nn1 = vol.numCoefs(0);
    int nn2 = vol.numCoefs(1);
    int nn3 = vol.numCoefs(2);
    int nn = nn1*nn2*nn3;

    int i1, i2, i3, idx;
    for (i3=0, idx=0; i3<nmb_w; ++i3)
	for (i2=0; i2<nmb_v; ++i2)
	    for (i1=0; i1<nmb_u; ++i1, ++idx)
	    {
		cout << result[idx].param[0] << " "  << result[idx].param[1] << " ";
		cout << result[idx].param[2] << endl;  
		cout << points[idx*dim] << " " << points[idx*dim+1] << " " << points[idx*dim+2] << " ";
		cout << der_u[idx*dim] << " " << der_u[idx*dim+1] << " " << der_u[idx*dim+2] << " ";
		cout << der_v[idx*dim] << " " << der_v[idx*dim+1] << " " << der_v[idx*dim+2] << " ";
		cout << der_w[idx*dim] << " " << der_w[idx*dim+1] << " " << der_w[idx*dim+2] << endl;

		vector<double> b0, bu, bv, bw;
		vol.computeBasis(result[idx].param, b0, bu, bv, bw);

		int left1 = result[idx].left_idx[0] - kk1 + 1;
		int left2 = result[idx].left_idx[1] - kk2 + 1;
		int left3 = result[idx].left_idx[2] - kk3 + 1;
		int kp, kq, kh, kr, kd;
		std::vector<Point> p1(4, Point(dim));
		p1[0].setValue(0.0);
		p1[1].setValue(0.0);
		p1[2].setValue(0.0);
		p1[3].setValue(0.0);
		for (kh=left3, kr=0; kh<left3+kk3; kh++)
		    for (kq=left2; kq<left2+kk2; kq++)
			for (kp=left1; kp<left1+kk1; kp++, kr++)
			    for (kd=0; kd<dim; kd++)
			    {
				double val = coefs[((kq+kh*nn2)*nn1+kp)*dim+kd];
				p1[0].begin()[kd] += val*b0[kr];
				p1[1].begin()[kd] += val*bu[kr];
				p1[2].begin()[kd] += val*bv[kr];
				p1[3].begin()[kd] += val*bw[kr];
		    }
		cout << p1[0][0] << " " << p1[0][1] << " " << p1[0][2] << " " ;
		cout << p1[1][0] << " " << p1[1][1] << " " << p1[1][2] << " " ;
		cout << p1[2][0] << " " << p1[2][1] << " " << p1[2][2] << " " ;
		cout << p1[3][0] << " " << p1[3][1] << " " << p1[3][2] << endl;

		std::vector<Point> p2(4, Point(dim));
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
				p2[0].begin()[kd] += val*result[idx].basisValues[kr];
				p2[1].begin()[kd] += val*result[idx].basisDerivs_u[kr];
				p2[2].begin()[kd] += val*result[idx].basisDerivs_v[kr];
				p2[3].begin()[kd] += val*result[idx].basisDerivs_w[kr];
			    }
		cout << p2[0][0] << " " << p2[0][1] << " " << p2[0][2] << " " ;
		cout << p2[1][0] << " " << p2[1][1] << " " << p2[1][2] << " " ;
		cout << p2[2][0] << " " << p2[2][1] << " " << p2[2][2] << " " ;
		cout << p2[3][0] << " " << p2[3][1] << " " << p2[3][2] << endl;

		std::vector<Point> p3(4, Point(dim));
		p3[0].setValue(0.0);
		p3[1].setValue(0.0);
		p3[2].setValue(0.0);
		p3[3].setValue(0.0);
		for (kr=0; kr<nn; kr++)
		    for (kd=0; kd<dim; kd++)
		    {
			double val = coefs[kr*dim+kd];
			p3[0].begin()[kd] += val*basisValues[idx][kr];
			p3[1].begin()[kd] += val*basisDerivs_u[idx][kr];
			p3[2].begin()[kd] += val*basisDerivs_v[idx][kr];
			p3[3].begin()[kd] += val*basisDerivs_w[idx][kr];
		    }
		cout << p3[0][0] << " " << p3[0][1] << " " << p3[0][2] << " " ;
		cout << p3[1][0] << " " << p3[1][1] << " " << p3[1][2] << " " ;
		cout << p3[2][0] << " " << p3[2][1] << " " << p3[2][2] << " " ;
		cout << p3[3][0] << " " << p3[3][1] << " " << p3[3][2] << endl;
		cout << endl;
	    }


    ofstream os("tmp_vol.g2");
    vol.writeStandardHeader(os);
    vol.write(os);
}





