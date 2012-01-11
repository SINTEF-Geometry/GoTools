//===========================================================================
//                                                                           
// File: face2splineSet
//                                                                           
// Description:
//  
// This program demonstrates how to create a set of spline surfaces,
// meeting in a corner-to-corner configuration and with corresponding
// coefficients at common boundaries, from one possibly trimmed face
//
// The program reads a bounded surface from a file, splits this surface
// into several bounded surfaces where each surface has (at most) 4 boundary
// curves. Finally, each bounded surface is approximated by a spline surface
// within a given tolerance and C0 continuities at common boundaries is
// ensured.
//
// Input/Output
// The file containing the input bounded surface is hardcoded. Tolerances
// are also hardcoded
// The result surfaces at different stages are written to specified files
//
//   
//===========================================================================

#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/RegularizeFace.h"
#include <fstream>

using std::vector;
using namespace Go;

int main( int argc, char* argv[] )
{
  // The input file contains one surface with a hole. It can be replaced
  // with other surfaces with inner and outer trimming.
  // Note that the functionality shown in this example program is still
  // under development, and may fail for some input surfaces
  std::string input_face("data/plane_sf_diamond_hole.g2");

  // The first output file contains a set of trimmed surfaces where each
  // surface has 4 boundaries and only outer trimming. The surfaces meet
  // in a corner-to-corner configuration, i.e. no T-joints
  std::string output_sfs1("data/trimmed_sfs.g2");

  // The second output file contains spline surfaces approximating the
  // 4-sided trimmed surfaces
  std::string output_sfs2("data/spline_sfs1.g2");

  // The last output file contains spline surfaces where the spline spaces
  // of the surfaces are identical along commong boundaries and where
  // corresponding coefficients belonging to neighbouring surfaces are
  // identical
  std::string output_sfs3("data/spline_sfs2.g2");

  // Prepare input and output files
  std::ifstream infile(input_face.c_str());
  std::ofstream of1(output_sfs1.c_str());
  std::ofstream of2(output_sfs2.c_str());
  std::ofstream of3(output_sfs3.c_str());

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
  CompositeModel *model = factory.createFromG2(infile);

  // A surface model inherits composite model, check if we have
  // a surface model
  SurfaceModel *sfmodel = dynamic_cast<SurfaceModel*>(model);

   if (sfmodel)
  {
    // The input file contained one or more surfaces. Fetch the 
    // face corresponding the first surface in this file
    shared_ptr<ftSurface> face = sfmodel->getFace(0);

    // Create class for splitting of one given face into one or more
    // faces with 4 boundaries
    std::cout << "Replace input surface with 4-sided surfaces"  << std::endl;
    RegularizeFace reg(face, gap, kink, neighbour);

    // Performe the split and fetch output faces
    // The gap tolerance is used when creating new trimming curves
    vector<shared_ptr<ftSurface> > sub_faces = reg.getRegularFaces();

    // Write corresponding surfaces to file
    for (size_t ki=0; ki<sub_faces.size(); ++ki)
      {
	// Fetch given surface from face set
	shared_ptr<ParamSurface> surf = sub_faces[ki]->surface();
	surf->writeStandardHeader(of1);
	surf->write(of1);
      }

    // Replace by spline surfaces
    // Boundary curves are approximated by spline curves, or fetched from
    // an already created adjacent spline surface. Spline surfaces are 
    // initially created as Coons patches interpolating the boundary curves
    // then updated with respect a point set fetched from the initial 
    // trimmed surface
    // The gap tolerance is used in the curve and surface approximations
    // First create a surface model representing the current face set
    std::cout << "Replace trimmed surfaces by spline surfaces"  << std::endl;
    shared_ptr<SurfaceModel> model2 =
      shared_ptr<SurfaceModel>(new SurfaceModel(approxtol, gap, neighbour,
						kink, bend, sub_faces,
						true));

    // Replace each trimmed surface by a spline surface
    model2->replaceRegularSurfaces();

    // Write current surface set to a file
    int nmb = model2->nmbEntities();
    for (int kr=0; kr<nmb; ++kr)
      {
	shared_ptr<ParamSurface> sf = model2->getSurface(kr);
	sf->writeStandardHeader(of2);
	sf->write(of2);
      }

    // Make sure that neighbouring surfaces have the same spline space
    // and coincidence of corresponding coefficients
    // First check if all surfaces are non-trimmed spline surfaces, i.e.
    // the previous opeation succeeded
    bool isOK = model2->allSplines();
    if (!isOK)
      {
	std::cout << "Not all surfaces are splines. Stopping computation" << std::endl;
	exit(-1);
      }
    
    // Ensure common spline spaces and corresponding coefficients
    std::cout << "Ensure common spline space"  << std::endl;
     model2->makeCommonSplineSpaces();

    // Write result to file
    for (int kr=0; kr<nmb; ++kr)
      {
	shared_ptr<ParamSurface> sf = model2->getSurface(kr);
	sf->writeStandardHeader(of3);
	sf->write(of3);
      }
  }
}


