//===========================================================================
//                                                                           
// File: interpcurve.C                                                       
//                                                                           
// Created: Mon Oct 16 18:46:34 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: interpcurve.C,v 1.7 2006-04-06 14:29:32 afr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/geometry/SplineInterpolator.h"
#include "GoTools/geometry/SplineApproximator.h"
#include "GoTools/geometry/SplineCurve.h"
#include <vector>

using namespace std;

int main()
{
    // Read data from standard input
    int dim, numpt;
    cin >> dim >> numpt;
    std::vector<double> param(numpt);
    std::vector<double> data(numpt * dim);
    for (int i = 0; i < numpt; ++i) {
	param[i] = static_cast<double>(i);
	for (int dd = 0; dd < dim; ++dd) {
	    cin >> data[i*dim + dd];
	}
    }

    // Interpolate and print out
     Go::SplineInterpolator interpol;
     interpol.setFreeConditions();
//      interpol.setNaturalConditions();
//      interpol.setHermiteConditions(Go::Point(0.0, 0.0, 1.0),
//  				  Go::Point(0.0, 0.0, -1.0));

//     Go::SplineApproximator interpol;
//     interpol.setNumCoefs(static_cast<int>(floor(max(4.0, numpt*0.2))));

    Go::SplineCurve cv;
    cv.interpolate(interpol, numpt, dim, &param[0], &data[0]);
//      std::vector<double> coefs;
//      interpol.interpolate(numpt, dim, param.begin(), data.begin(), coefs);

//      Go::SplineCurve cv(interpol.basis().numCoefs(),
//  			 interpol.basis().order(),
//  			 interpol.basis().begin(),
//  			 coefs.begin(),
//  			 dim);

    cout.precision(15);
    cv.writeStandardHeader(cout);
    cout << cv;

    vector<Go::Point> p(3, Go::Point(1));
    for (int i = 0; i < numpt; ++i) {
	cv.point(p, param[i], 2);
	cout << "Pts: " << p[0] << '\n' << p[1] << '\n' << p[2] << endl;
    }
}
