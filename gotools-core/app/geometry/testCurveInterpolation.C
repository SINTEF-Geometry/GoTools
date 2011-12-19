#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/CurveInterpolator.h"
#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <iostream>
#include <fstream>


using namespace Go;
using std::ifstream;
using std::ofstream;
using std::vector;


int main(int argc, char* argv[])
{
    if (argc != 3) {
	std::cout << "Usage:  input file, output file"
		  << std::endl;
	return -1;
    }

    SplineCurve cv_in;
    ObjectHeader header;

    ifstream infile(argv[1]);
    infile >> header >> cv_in;
    ofstream outfile(argv[2]);

    // Fetch Greville parameters
    BsplineBasis basis = cv_in.basis();
    int nmb = basis.numCoefs();
    vector<double> par(nmb);
    int ki;
    for (ki=0; ki<nmb; ++ki)
      par[ki] = basis.grevilleParameter(ki);
    
    // Evaluate the curve in a regular grid
    vector<double> points;
    cv_in.gridEvaluator(points, par);
   
    // Fetch weights
    vector<double> weights;
    if (cv_in.rational())
      cv_in.getWeights(weights);

    // Make surface
    shared_ptr<SplineCurve> cv_out = 
      shared_ptr<SplineCurve>(CurveInterpolator::regularInterpolation(basis,
								      par,
								      points,
								      cv_in.dimension(),
								      cv_in.rational(),
								      weights));

    cv_out->writeStandardHeader(outfile);
    cv_out->write(outfile);
}
