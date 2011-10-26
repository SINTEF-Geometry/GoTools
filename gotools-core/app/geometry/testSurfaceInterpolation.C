#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SurfaceInterpolator.h"
#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <iostream>
#include <fstream>


using namespace Go;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::shared_ptr;


int main(int argc, char* argv[])
{
    if (argc != 3) {
	std::cout << "Usage:  input file, output file"
		  << std::endl;
	return -1;
    }

    SplineSurface sf_in;
    ObjectHeader header;

    ifstream infile(argv[1]);
    infile >> header >> sf_in;
    ofstream outfile(argv[2]);

    // Fetch Greville parameters
    BsplineBasis basis_u = sf_in.basis_u();
    int nmb_u = basis_u.numCoefs();
    BsplineBasis basis_v = sf_in.basis_v();
    int nmb_v = basis_v.numCoefs();
    vector<double> par_u(nmb_u);
    vector<double> par_v(nmb_v);
    int ki;
    for (ki=0; ki<nmb_u; ++ki)
      par_u[ki] = basis_u.grevilleParameter(ki);
    for (ki=0; ki<nmb_v; ++ki)
      par_v[ki] = basis_v.grevilleParameter(ki);
    
    // Evaluate the surface in a regular grid
    vector<double> points;
    sf_in.gridEvaluator(points, par_u, par_v);
   
    // Fetch weights
    vector<double> weights;
    if (sf_in.rational())
      sf_in.getWeights(weights);

    // Make surface
    shared_ptr<SplineSurface> sf_out = 
      shared_ptr<SplineSurface>(SurfaceInterpolator::regularInterpolation(basis_u,
									 basis_v,
									 par_u,
									 par_v,
									 points,
									 sf_in.dimension(),
									 sf_in.rational(),
									 weights));

    sf_out->writeStandardHeader(outfile);
    sf_out->write(outfile);
}
