//===========================================================================
//                                                                           
// File: multiPatchSweep
//                                                                           
//===========================================================================

#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/trivariate/SweepVolumeCreator.h"
#include <fstream>

using namespace Go;
using std::cout;


//===========================================================================
//                                                                           
// Description:
//  
// The idea of this program is to sweep a set surfaces to create a multi patch 
// volume model.
// The surface set will be created by the example programs in compositemodel:
// createSplitDisc and createBlockStructuredDisc
//
// Input/Output
// Input is the file data/split_disc4.g2
// The output volumes are written to data/swept_vol.g2
//                                                                           
//===========================================================================


int main( int argc, char* argv[] )
{
  // Prepare for input and output files
  std::string infile("data/split_disc4.g2");
  std::string outfile("data/swept_vol.g2");
  std::ifstream input("infile.c_str()");

  // Read input data
  // Define tolerances
  // The neighbour tolerance is used in topology build if more than
  // one surface are given. Two surfaces lying closer than this tolerance
  // are assumed to be neighbours 
  double neighbour = 0.001;
  // Two surface lying closer than the neighbour tolerance, but more distant
  // than the gap tolerance are assumed to be neighbours, but the surface
  // set is not found to be C0 continuous within the given tolerance (gap)
  double gap = 0.0001;
  // Two neighbouring surfaces where the angle between some corresponding
  // surface normals are larger than the bend tolerances are found to
  // create an intential corner. Angular tolerances are given in radians.
  double bend = 0.01;
  // If the angle between corresponding surface normals are larger than the
  // kink tolerance but smaller than the bend tolerance the surface set is 
  // found to be intentially, but not truely G1. If all angles are less than
  // the kink tolerance, the surface set is seen as G1
  double kink = 0.001;
  // Tolerance intended for approximations
  double approxtol = 0.0001;

  // Create a factory class to read/create composite models, most often
  // one or more surfaces where the adjacency relationship between the
  // surfaces are known
  std::cout << "Reading input data" << std::endl;
  CompositeModelFactory factory(approxtol, gap, neighbour, kink, bend);

  // Read data from file. At this stage, it is not known whether the
  // model consists of curves or surfaces
  shared_ptr<CompositeModel> model(factory.createFromG2(input));

  // A surface model inherits composite model, check if we have
  // a surface model
  shared_ptr<SurfaceModel> sfmodel = 
    dynamic_pointer_cast<SurfaceModel,CompositeModel>(model);

if (sfmodel.get())
  {
  
    // Make sure that neighbouring surfaces have the same spline space
    // and coincidence of corresponding coefficients
    // First check if all surfaces are non-trimmed spline surfaces, i.e.
    // the previous opeation succeeded
    bool isOK = sfmodel->allSplines();
    if (!isOK)
      {
	std::cout << "Not all surfaces are splines. Stopping computation" << std::endl;
	exit(-1);
      }
    
    // Create linear swept volumes
    std::cout << "Create swept volumes" << std::endl;
    vector<shared_ptr<ftVolume> > blocks;  // Storage for volumes including
                                           // topology information
    double sweep_length_fac = 5.0;
    for (ki=0; ki<nmb; ++ki)
      {
	// Fetch current spline surface
	shared_ptr<ParamSurface> sf = model4->getSurface(ki);
	shared_ptr<SplineSurface> surf = 
	  dynamic_pointer_cast<SplineSurface,ParamSurface>(sf);

	if (surf.get())
	  {
	    // Create swept volume
	    shared_ptr<SplineCurve> cv(new SplineCurve(centre, 
						       centre+sweep_length_fac*z_axis));
	    shared_ptr<ParamVolume> vol = 
	      shared_ptr<ParamVolume>(SweepVolumeCreator::linearSweptVolume(*surf, 
									    *cv, 
									    centre));
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
  std::ofstream of(outfile.c_str());
  for (ki=0; ki<nmb_vol; ++ki)
    {
      shared_ptr<ParamVolume> vol = volmodel->getVolume(ki);
       vol->writeStandardHeader(of);
       vol->write(of);
     }
  
}

  
		      
