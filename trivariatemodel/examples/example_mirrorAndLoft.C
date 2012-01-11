//===========================================================================
//                                                                           
// File: example_mirrorAndLoft
//                                                                           
//===========================================================================

#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/trivariate/LoftVolumeCreator.h"
#include <fstream>

using namespace Go;


// Description:
//  
// This program demonstrates how to create a volume model from a set
// of spline surfaces. The surfaces must be connected. The surface set
// is mirrored around a given plane, and lofting between corresponding 
// surfaces is performed.
//
// Input/Output
// The file containing the input surface set and the plane specifiction 
// is hardcoded. Tolerances are also hardcoded
// The resulting volume model is written to a specifified file


int main( int argc, char* argv[] )
{
  // The input file contains a connected surface set where each surface
  // is expected to be a spline surface, rational or not
  std::string input("data/nut_top.g2");
  
  // The output file will contain a set of spline volumes. To visualize
  // using goview, use the application program trivariate/app/makeShields 
  // to feth the boundary surfaces of all volumes
  std::string output("data/nut.g2");

  // Prepare input and output files
  std::ifstream infile(input.c_str());
  std::ofstream outfile(output.c_str());

  // Specify the plane around which to mirror the surface set
  Point pnt(0.0, 0.4, 0.0);
  Point norm(0.0, 1.0, 0.0);

  // Tolerances used in topology build
  // The neighbour tolerance is used in topology build if more than
  // one surface are given. Two surfaces lying closer than this tolerance
  // are assumed to be neighbours 
  double neighbour = 0.01;
  // Two surface lying closer than the neighbour tolerance, but more distant
  // than the gap tolerance are assumed to be neighbours, but the surface
  // set is not found to be C0 continuous within the given tolerance (gap)
  // The gap tolerance must be smaller than the neighbouring tolerance
  double gap = 0.0001;
  // Two neighbouring surfaces where the angle between some corresponding
  // surface normals are larger than the bend tolerances are found to
  // create an intential corner. Angular tolerances are given in radians.
  double bend = 0.5;
  // If the angle between corresponding surface normals are larger than the
  // kink tolerance but smaller than the bend tolerance the surface set is 
  // found to be intentially, but not truely G1. If all angles are less than
  // the kink tolerance, the surface set is seen as G1
  // The kink tolerance must be smaller than the bend tolerance
  double kink = 0.01;
  // Approximation tolerance. Not used in this context
  double approxtol = 0.01;

  // Create a factory class to read/create composite models, most often
  // one or more surfaces where the adjacency relationship between the
  // surfaces are known
  std::cout << "Reading input data" << std::endl;
  CompositeModelFactory factory(approxtol, gap, neighbour, kink, bend);

  // Read data from file. At this stage, it is not known whether the
  // model consists of curves or surfaces
  shared_ptr<CompositeModel> model = shared_ptr<CompositeModel>(factory.createFromG2(infile));
  // A surface model inherits composite model, check if we have
  // a surface model
  shared_ptr<SurfaceModel> sfmodel = dynamic_pointer_cast<SurfaceModel,CompositeModel>(model);
  if (!sfmodel.get())
    {
      std::cout << "No surface set is specified. Stop execution" << std::endl;
      exit(-1);
    }

  // Check the input surface set
  std::cout << "Check surface set"  << std::endl;
  bool isOK = sfmodel->allSplines();
  if (!isOK)
    {
      std::cout << "Not all surfaces are splines. Stopping computation" << std::endl;
      exit(-1);
    }

  // Check the number of connected parts in the surface set
  // The tolerances specified previously are stored in the surface model
  // and used for testing
  std::vector<shared_ptr<SurfaceModel> > models = sfmodel->getConnectedModels();
  if (models.size() != 1)
    {
      std::cout << "Not a connected model. Stop execution" << std::endl;
      exit(-1);
    }

  // Check for gaps in the model
  ftCurve gaps = sfmodel->getGaps();
  if (gaps.numSegments() > 0)
    {
      // There are at least one gap in the model. 
      std::cout << "Gap(s) in the model. Stop execution" << std::endl;
      exit(-1);
    }

  // Check if the surface set has a corner-to-corner configuration.
  // The mirror and loft operations are not dependent on this configuration,
  // but the topology build for volumes is currently limited with respect to 
  // configurations. 
  bool corner_config = sfmodel->isCornerToCorner();
  if (!corner_config)
    {
      std::cout << "Not a corner-to-corner configuration. Stop execution" << std::endl;
      exit(-1);
    }
      
  // Ensure common spline spaces and corresponding coefficients
  // This is not necessary for the operations in this application, but
  // ensures a volume model with corresponding coefficients at the common
  // boundary for adjacent volumes
  std::cout << "Ensure common spline space"  << std::endl;
  sfmodel->makeCommonSplineSpaces();

   // Fetch the number of surfaces in the surface set
  int nmb_sfs = sfmodel->nmbEntities();
  vector<shared_ptr<ParamSurface> > mirrored;  // Storage for mirrored surfaces
  vector<shared_ptr<ftVolume> > blocks;  // Storage for volumes including
                                         // topology information
  std::cout << "Perform mirror and loft"  << std::endl;
  for (int ki=0; ki<nmb_sfs; ++ki)
    {
      // Fetch the current surface
      shared_ptr<ParamSurface> curr = sfmodel->getSurface(ki);

      // Mirror around the given plane
      shared_ptr<ParamSurface> mod = 
	shared_ptr<ParamSurface>(curr->mirrorSurface(pnt, norm));
      mirrored.push_back(mod);

      // Represent the original and the mirrored surface as spline surfaces
      vector<shared_ptr<SplineSurface> > sfs(2);
      sfs[0] = dynamic_pointer_cast<SplineSurface,ParamSurface>(curr);
      sfs[1] = dynamic_pointer_cast<SplineSurface,ParamSurface>(mod);

      if (sfs[0].get() && sfs[1].get())
	{
	  // Perform lofting
	  shared_ptr<ParamVolume> vol = 
	    shared_ptr<ParamVolume>(LoftVolumeCreator::loftVolume(sfs.begin(), 2));
	  // Add data structures for volume topology
	  shared_ptr<ftVolume> ftvol = 
	    shared_ptr<ftVolume>(new ftVolume(vol, gap, kink));
	  blocks.push_back(ftvol);
	}
    }

  // Create a volume model. The topology build performs coincidence testing
  // between possible adjacent volumes. This may be time consuming
  std::cout << "Volume topology build"  << std::endl;
  shared_ptr<VolumeModel> volmodel = 
    shared_ptr<VolumeModel>(new VolumeModel(blocks, gap, neighbour, 
					    kink, 10.0*kink));


  // Fetch the number of volumes
  int nmb_vol = volmodel->nmbEntities();

  // Write the volumes to the output file
  // With the current file format, the topology information is lost
  for (int ki=0; ki<nmb_vol; ++ki)
    {
      shared_ptr<ParamVolume> vol = volmodel->getVolume(ki);
       vol->writeStandardHeader(outfile);
       vol->write(outfile);
     }
}
  
