#include <fstream>
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/SmoothVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/errormacros.h"
#include "GoTools/trivariate/CoonsPatchVolumeGen.h"

using namespace Go;
using namespace std;

//===========================================================================
//                                                                           
// File: createCoonsVolume
//                                                                           
//                                                                           
// Description:
//
// This program demonstrates the use of the static function 'createCoonsPatch'
// in namespace 'CoonsPatchVolumeGen'.
// The function can create a new 'SplineVolume' representing the coons patch of
// six SplineSurfaces, the six faces of the volume.
//
// This program reads the 6 boundary surfaces needed to create the volume. 
// The surfaces must be non-rational, lie in the same space, and have
// conciding edge curves where needed. The input file must satisfy these
// conditions. This is not tested.
// The B-spline bases of the surfaces might be reversed, swapped, have order
// raised, knots inserted and parameter interval rescaled to [0,1].
// The input filename is data/volume_boundaries.g2. This file is created
// by the example program createVolumeBoundaries in compositemodel.
//
// The input surfaces are not all particularily "nice". Thus, a Coons 
// approach is not expected to give a very good parameterization. To 
// improve the representation, the volume is smoothed, i.e. the inner
// coefficients of the volume is modified to minimize a smoothing
// funcitonal. The volume is written to the file data/volume1.g2 prior to
// smoothing and data/volume1.g2 after smoothing
// 
//===========================================================================

int main(int argc, char* argv[] )
{
   // Open input volume file.
  ifstream is("data/volume_boundaries.g2");
  ALWAYS_ERROR_IF(is.bad(), "Bad or no volume input filename");

  // Read surfaces
  vector<shared_ptr<SplineSurface> > bdsf(6);
  for (int ki=0; ki<6; ++ki)
    {
        ObjectHeader head;
	head.read(is);
	bdsf[ki] = shared_ptr<SplineSurface>(new SplineSurface());
	bdsf[ki]->read(is);
    }

  std::cout << "Create Coons volume" << std::endl;

  // Create Coons volume
  shared_ptr<SplineVolume> vol(CoonsPatchVolumeGen::createCoonsPatch(bdsf[0],
								     bdsf[1],
								     bdsf[2],
								     bdsf[3],
								     bdsf[4],
								     bdsf[5]));

  // Write to file
  ofstream of1("data/volume1.g2");
  vol->writeStandardHeader(of1);  
  vol->write(of1);

  std::cout << "Volume smoothing" << std::endl;

  // Perform smoothing
  SmoothVolume smooth(true);  // Create engine. Do not modify input volume

  // Set all volume coefficients to be free, except those at the boundaries
  int numcfs[3];
  for (int ki=0; ki<3; ++ki)
    numcfs[ki] = vol->numCoefs(ki);
  vector<CoefStatus> coefstat(numcfs[0]*numcfs[1]*numcfs[2], CoefFree);

  // Fix coefficients at the boundaries
  for (int kr=0; kr<numcfs[2]; kr+=(numcfs[2]-1))
    for (int kj=0; kj<numcfs[1]; ++kj)
      for (int ki=0; ki<numcfs[0]; ++ki)
	coefstat[(kj+kr*numcfs[2])*numcfs[1]+ki] = CoefKnown;

  for (int kr=0; kr<numcfs[2]; ++kr)
    for (int kj=0; kj<numcfs[1]; kj+=(numcfs[1]-1))
      for (int ki=0; ki<numcfs[0]; ++ki)
	coefstat[(kj+kr*numcfs[2])*numcfs[1]+ki] = CoefKnown;

  for (int kr=0; kr<numcfs[2]; ++kr)
    for (int kj=0; kj<numcfs[1]; ++kj)
      for (int ki=0; ki<numcfs[0]; ki+=(numcfs[0]-1))
	coefstat[(kj+kr*numcfs[2])*numcfs[1]+ki] = CoefKnown;

  // Attach data to smoothing engine
  smooth.attach(vol, coefstat);

  // Define weights. Only smoothing is performed, thus only smoothing
  // weights are set
  double wgt1 = 1.0e-5; // Minimize expression in 1. derivative. Use with care
  double wgt2 = 0.5;    // Minimize expression in 2. derivative. 
  double wgt3 = 1.0 - wgt2 - wgt1; // Minimize expression in 3. derivative. 
  smooth.setOptimize();

  // Perform smoothing and fetch result
  shared_ptr<SplineVolume> vol2;
  int status = smooth.equationSolve(vol2);
  std::cout << "Volume smoothing status: " << status << std::endl;

  // Write to file
  ofstream of2("data/volume2.g2");
  vol->writeStandardHeader(of2);  
  vol->write(of2);
}

