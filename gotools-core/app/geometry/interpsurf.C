//===========================================================================
//                                                                           
// File: interpsurf.C                                                        
//                                                                           
// Created: Tue Oct 17 16:44:37 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: interpsurf.C,v 1.5 2005-06-17 05:51:58 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================



#include "GoTools/geometry/SplineInterpolator.h"
#include "GoTools/geometry/SplineApproximator.h"
#include "GoTools/geometry/SplineSurface.h"
#include <vector>

using namespace std;

int main()
{
    // Read data from standard input
    int dim, numpt1, numpt2;
    cin >> dim >> numpt1 >> numpt2;
    std::vector<double> param1(numpt1);
    std::vector<double> param2(numpt2);
    std::vector<double> data(numpt1*numpt2 * dim);
    for (int i = 0; i < numpt1; ++i)
	cin >> param1[i];
     for (int i = 0; i < numpt2; ++i)
	cin >> param2[i];
    for (int i = 0; i < numpt1*numpt2; ++i) {
	for (int dd = 0; dd < dim; ++dd) {
	    cin >> data[i*dim + dd];
	}
    }

    // Interpolate and print out
    Go::SplineInterpolator interpol1;
    interpol1.setNaturalConditions();
    Go::SplineInterpolator interpol2;
    interpol2.setNaturalConditions();
//      Go::SplineApproximator interpol1;
//      interpol1.setNumCoefs(4);
//      Go::SplineApproximator interpol2;
//      interpol2.setNumCoefs(4);

    Go::SplineSurface sf;
    sf.interpolate(interpol1, interpol2,
		   numpt1, numpt2, dim,
		   &param1[0], &param2[0], &data[0]);

    cout.precision(15);
    sf.writeStandardHeader(cout);
    cout << sf;
}
