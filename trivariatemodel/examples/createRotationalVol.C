//===========================================================================
//                                                                           
// File: createRotationalVol
//                                                                           
//===========================================================================

#include <fstream>
#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/trivariatemodel/VolumeModelCreator.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"

//===========================================================================
///                                                                           
/// Description:
///  
/// Create a block structured volume model from a face set
/// describing a boundary represented solid representing a
/// rotational object.
///
/// Input/Output:
///
/// Input to the geometry construction is the file data/rotational_model.g2
/// The output volumes are stored in the file data/rotational_model_out.g2
///                                                                           
//===========================================================================

using namespace Go;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::vector;

int main(int argc, char* argv[] )
{
  // Prepare for input and output files
  std::string infile("data/rotational_model.g2");
  std::string outfile("data/rotational_model_out.g2");

  ifstream input(infile.c_str());

  ofstream of(outfile.c_str());

  // Define tolerances
  // The neighbour tolerance is used in topology build if more than
  // one surface are given. Two surfaces lying closer than this tolerance
  // are assumed to be neighbours 
  double neighbour = 0.01;
  // Two surface lying closer than the neighbour tolerance, but more distant
  // than the gap tolerance are assumed to be neighbours, but the surface
  // set is not found to be C0 continuous within the given tolerance (gap)
  double gap = 0.001;
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
  double approxtol = 0.001;

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
    // Check if the model is rotational. In that case intersect the model
    // with a planar surface with one boundary at the rotational axis that
    // is larger that one sector of the rotational model. The resulting
    // planar trimmed surface is represented as a block structured model 
    // where each block is a NURBS surface. Finally, perform rotational 
    // sweep to create volume blocks.
    shared_ptr<VolumeModel> volmod;
    bool rotational = 
      VolumeModelCreator::createRotationalModel(sfmodel, volmod);
    std::cout << "Rotational: " << rotational << std::endl;

    if (rotational)
      {
	// The model can be created by rotational sweep and a result is
	// obtained. Write the result to file. Note that the file contains
	// volumes. To visualize them use trivariate/app/makeShields to
	// extract the boundary surfaces of each volume
        int nmb_vols = volmod->nmbEntities();
	for (int kr=0; kr<nmb_vols; ++kr)
	  {
	    shared_ptr<ParamVolume> curr_vol2 = volmod->getVolume(kr);
	    curr_vol2->writeStandardHeader(of);
	    curr_vol2->write(of);
	  }
    }
  }

}

