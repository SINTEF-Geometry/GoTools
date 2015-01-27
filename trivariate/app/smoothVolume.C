#include <fstream>
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/SmoothVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/errormacros.h"
#include "GoTools/trivariate/CoonsPatchVolumeGen.h"

using namespace Go;
using namespace std;

int main(int argc, char* argv[] )
{
  if (argc != 5) {
    std::cout << "Usage: input volumes, output volumes, output surfaces 1, output surfaces 2" << std::endl;
    return -1;
  }

  std::ifstream is(argv[1]);
  std::ofstream os(argv[2]);
  std::ofstream os2(argv[3]);
  std::ofstream os3(argv[4]);

  while (!is.eof())
    {
      ObjectHeader head;
      is >> head;

      shared_ptr<SplineVolume> vol(new SplineVolume());
      vol->read(is);

      // Write constant parameter surfaces to file. Visualization purpose
      int kj;
      double start = vol->startparam(0);
      double end = vol->endparam(0);
      double del = (end - start)/(double)4;
      for (kj=0; kj<5; ++kj, start+=del)
	{
	  shared_ptr<SplineSurface> tmp(vol->constParamSurface(start, 0));
	  tmp->writeStandardHeader(os2);
	  tmp->write(os2);
	}

      start = vol->startparam(1);
      end = vol->endparam(1);
      del = (end - start)/(double)4;
      for (kj=0; kj<5; ++kj, start+=del)
	{
	  shared_ptr<SplineSurface> tmp(vol->constParamSurface(start, 1));
	  tmp->writeStandardHeader(os2);
	  tmp->write(os2);
	}
      start = vol->startparam(2);
      end = vol->endparam(2);
      del = (end - start)/(double)4;
      for (kj=0; kj<5; ++kj, start+=del)
	{
	  shared_ptr<SplineSurface> tmp(vol->constParamSurface(start, 2));
	  tmp->writeStandardHeader(os2);
	  tmp->write(os2);
	}

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
	    coefstat[(kj+kr*numcfs[1])*numcfs[0]+ki] = CoefKnown;

      for (int kr=0; kr<numcfs[2]; ++kr)
	for (int kj=0; kj<numcfs[1]; kj+=(numcfs[1]-1))
	  for (int ki=0; ki<numcfs[0]; ++ki)
	    coefstat[(kj+kr*numcfs[1])*numcfs[0]+ki] = CoefKnown;

      for (int kr=0; kr<numcfs[2]; ++kr)
	for (int kj=0; kj<numcfs[1]; ++kj)
	  for (int ki=0; ki<numcfs[0]; ki+=(numcfs[0]-1))
	    coefstat[(kj+kr*numcfs[1])*numcfs[0]+ki] = CoefKnown;

      // Attach data to smoothing engine
      smooth.attach(vol, coefstat);

      // Define weights. Only smoothing is performed, thus only smoothing
      // weights are set
      double wgt1 = 1.0e-5; // Minimize expression in 1. derivative. Use with care
      double wgt2 = 0.5;    // Minimize expression in 2. derivative. 
      double wgt3 = 1.0 - wgt2 - wgt1; // Minimize expression in 3. derivative. 
      smooth.setOptimize(wgt1, wgt2, wgt3);

      // Perform smoothing and fetch result
      shared_ptr<SplineVolume> vol2;
      int status = smooth.equationSolve(vol2);
      std::cout << "Volume smoothing status: " << status << std::endl;

      // Write to file
      vol2->writeStandardHeader(os);  
      vol2->write(os);

      // Write constant parameter surfaces to file. Visualization purpose
      start = vol2->startparam(0);
      end = vol2->endparam(0);
      del = (end - start)/(double)4;
      for (kj=0; kj<5; ++kj, start+=del)
	{
	  shared_ptr<SplineSurface> tmp(vol2->constParamSurface(start, 0));
	  tmp->writeStandardHeader(os3);
	  tmp->write(os3);
	}

      start = vol2->startparam(1);
      end = vol2->endparam(1);
      del = (end - start)/(double)4;
      for (kj=0; kj<5; ++kj, start+=del)
	{
	  shared_ptr<SplineSurface> tmp(vol2->constParamSurface(start, 1));
	  tmp->writeStandardHeader(os3);
	  tmp->write(os3);
	}
      start = vol2->startparam(2);
      end = vol2->endparam(2);
      del = (end - start)/(double)4;
      for (kj=0; kj<5; ++kj, start+=del)
	{
	  shared_ptr<SplineSurface> tmp(vol2->constParamSurface(start, 2));
	  tmp->writeStandardHeader(os3);
	  tmp->write(os3);
	}
    }
}

