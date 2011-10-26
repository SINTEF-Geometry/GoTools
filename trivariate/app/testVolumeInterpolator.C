#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/VolumeInterpolator.h"
#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <iostream>
#include <fstream>


using namespace Go;
using std::ifstream;
using std::ofstream;
using std::shared_ptr;
using std::vector;


int main(int argc, char* argv[])
{
    if (argc != 3) {
	std::cout << "Usage:  input file, output file"
		  << std::endl;
	return -1;
    }

    SplineVolume vol_in;
    ObjectHeader header;

    ifstream infile(argv[1]);
    infile >> header >> vol_in;
    ofstream outfile(argv[2]);

    // Fetch Greville parameters
    BsplineBasis basis_u = vol_in.basis(0);
    int nmb_u = basis_u.numCoefs();
    BsplineBasis basis_v = vol_in.basis(1);
    int nmb_v = basis_v.numCoefs();
    BsplineBasis basis_w = vol_in.basis(2);
    int nmb_w = basis_w.numCoefs();
    vector<double> par_u(nmb_u);
    vector<double> par_v(nmb_v);
    vector<double> par_w(nmb_w);
    int ki;
    for (ki=0; ki<nmb_u; ++ki)
      par_u[ki] = basis_u.grevilleParameter(ki);
    for (ki=0; ki<nmb_v; ++ki)
      par_v[ki] = basis_v.grevilleParameter(ki);
    for (ki=0; ki<nmb_w; ++ki)
      par_w[ki] = basis_w.grevilleParameter(ki);
    
    // Evaluate the volume in a regular grid
    vector<double> points;
    vol_in.gridEvaluator(par_u, par_v, par_w, points);
   
    // Fetch weights
    vector<double> weights;
    if (vol_in.rational())
      vol_in.getWeights(weights);

    // Make volume
    shared_ptr<SplineVolume> vol_out = 
      shared_ptr<SplineVolume>(VolumeInterpolator::regularInterpolation(basis_u,
									basis_v,
									basis_w,
									par_u,
									par_v,
									par_w,
									points,
									vol_in.dimension(),
									vol_in.rational(),
									weights));

    vol_out->writeStandardHeader(outfile);
    vol_out->write(outfile);
}
