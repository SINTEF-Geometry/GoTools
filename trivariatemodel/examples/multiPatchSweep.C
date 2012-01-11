//===========================================================================
//                                                                           
// File: multiPatchSweep
//                                                                           
// Description:
//  
// The idea of this program is to create a set of spline surfaces and
// sweep each of those surfaces to create a multi patch volume model.
// The surface set will be created first by a number of split operations
// and regularization of a face set to approximate a number of bounded
// surfaces by spline surfaces in a corner-to-corner configuaration. Then
// the surfaces are swept to create volumes
// The construction uses planar, rectangular surfaces and a truncated sylinder,
// but the operations performed using these surfaces, do not depend on that
// level of regularity.
//
// Input/Output
// Input to the geometry construction is hardcoded
// The current surfaces and volumes are written to g2-files as we go along
//
//   
//===========================================================================

#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/RegularizeFaceSet.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SweepSurfaceCreator.h"
#include "GoTools/trivariate/SweepVolumeCreator.h"
#include <fstream>

using namespace Go;
using std::cout;

int main( int argc, char* argv[] )
{
  // Prepare for output files
  std::string outfile1("data/rect_sf.g2");
  std::string outfile2("data/cylinder_sf.g2");
  std::string outfile3("data/disc.g2");
  std::string outfile4("data/rect_tube.g2");
  std::string outfile5("data/split_disc1.g2");
  std::string outfile6("data/split_disc2.g2");
  std::string outfile7("data/split_disc3.g2");
  std::string outfile8("data/split_disc4.g2");
  std::string outfile9("data/swept_vol.g2");
  
  // First a rectangular spline surface is created as the base surface
  // for a disc represented as a bounded surface
  // Define data
  std::cout << "Creating underlying surface" << std::endl;
  int dim = 3;
  int numcoefs = 2;
  int order = 2;
  double knots[4] = {0.0, 0.0, 1.0, 1.0};
  double coefs[12] = {-2.0, -2.0, 0.0,
		      2.0, -2.0, 0.0,
		      -2.0, 2.0, 0.0,
		      2.0, 2.0, 0.0};

  // Create spline rectangle
  shared_ptr<SplineSurface> rect(new SplineSurface(numcoefs, numcoefs,
						   order, order,
						   knots, knots, coefs,
						   dim));

  std::ofstream of1(outfile1.c_str());
  rect->writeStandardHeader(of1);
  rect->write(of1);

  // Create cylinder
  std::cout << "Creating cylinder " << std::endl;
  double radius = 1.5;
  Point centre(0.0, 0.0, 0.0);
  Point z_axis(0.0, 0.0, 1.0);
  Point x_axis(1.0, 0.0, 0.0);

  shared_ptr<Cylinder> cylinder(new Cylinder(radius, centre, z_axis, x_axis));

  // Restrict the cylinder height. The 1. parameter directions is
  // running from 0 to 2*pi, while the 2. is infinite
  cylinder->setParameterBounds(0.0, -1, 2*M_PI, 1);

  // Represent the cylinder as a NURBS surface
  shared_ptr<SplineSurface> cyl(cylinder->createSplineSurface());

  std::ofstream of2(outfile2.c_str());
  cyl->writeStandardHeader(of2);
  cyl->write(of2);

  // Represent each surface as one surface set to perform a Boolean operation
  // between these two surface set
  // First define topology tolerances
  double neighbour = 0.01;
  double gap = 0.0001;
  double bend = 0.5;
  double kink = 0.01;
  // Approximation tolerance. Not used in this context
  double approxtol = 0.01;

  vector<shared_ptr<ParamSurface> > set1(1);
  set1[0] = rect;
  shared_ptr<SurfaceModel> model1(new SurfaceModel(approxtol, gap, neighbour,
						   kink, bend, set1));
  vector<shared_ptr<ParamSurface> > set2(1);
  set2[0] = cyl;
  shared_ptr<SurfaceModel> model2(new SurfaceModel(approxtol, gap, neighbour,
						   kink, bend, set2));
  
  // Split both surface models with respect to the intersection curves
  // between them, and fetch the circular disc
 std:cout << "Intersecting base rectangle with cylinder" << std::endl;
  vector<shared_ptr<SurfaceModel> > sub_models = 
    model1->splitSurfaceModels(model2);
  if (sub_models.size() != 4)
    {
      std::cout << "Unexpected number of pieces" << std::endl;
      exit(-1);
    }

  // By construction, the number of models after split should be 4 and
  // the first one is the part of the first model lying inside the
  // second surface set
  shared_ptr<SurfaceModel> disc = sub_models[0];

  std::ofstream of3(outfile3.c_str());
  disc->getSurface(0)->writeStandardHeader(of2);
  disc->getSurface(0)->write(of2);

  // Create a quadratic tube intersecting the disc in the inner in order
  // to split the disc into two parts
  // Start by creating the lower boundary curve in each of the surface
  // defining the tube
  std::cout << "Creating rectangular cube as swept surfaces" << std::endl;
  vector<Point> pnts(4);
  pnts[0] = Point(-0.5, -0.5, -1);
  pnts[1] = Point(0.5, -0.5, -1);
  pnts[2] = Point(0.5, 0.5, -1);
  pnts[3] = Point(-0.5, 0.5, -1);
  
  // Each curve is a linear spline curve on the parameter interval
  // [0,1] which interpolates two tube corners
  vector<shared_ptr<SplineCurve> > cvs(4);
  int ki, kj;
  for (ki=0; ki<4; ++ki)
    {
      kj = (ki+1)%4;
      cvs[ki] = 
	shared_ptr<SplineCurve>(new SplineCurve(pnts[ki], pnts[kj]));
    }

  // Perform a linear sweep of each curve in the direction of the
  // cylinder axis to create surfaces
  // Create swept surfaces
  SweepSurfaceCreator sweep;
  vector<shared_ptr<ParamSurface> > swept_sfs(4);
  for (ki=0; ki<4; ++ki)
    {
      shared_ptr<SplineCurve> axis(new SplineCurve(pnts[ki],
						   pnts[ki]+2.0*z_axis));
      swept_sfs[ki] = 
	shared_ptr<SplineSurface>(sweep.linearSweptSurface(*cvs[ki],
							   *axis, 
							   pnts[ki]));
    }

  std::ofstream of4(outfile4.c_str());
  for (ki=0; ki<4; ++ki)
    {
      swept_sfs[ki]->writeStandardHeader(of2);
      swept_sfs[ki]->write(of2);
    }

  // Create surface model containing the swept surfaces
  shared_ptr<SurfaceModel> model3(new SurfaceModel(approxtol, gap, neighbour,
						   kink, bend, swept_sfs));

  // Split the dist with respect to the tube model and the other way around
  std::cout << "Split disc by intersecting with rectangular tube" << std::endl;
  vector<shared_ptr<SurfaceModel> > sub_models2 = 
    disc->splitSurfaceModels(model3);
  if (sub_models2.size() != 4)
    {
      std::cout << "Unexpected number of pieces" << std::endl;
      exit(-1);
    }
 
  // Now we want to keept both pieces of the disc surface. They are stored 
  // the first two surface models. Join both parts into one model
  sub_models2[0]->append(sub_models2[1]);
  
  std::ofstream of5(outfile5.c_str());
  int nmb = sub_models2[0]->nmbEntities();
  for (ki=0; ki<nmb; ++ki)
    {
      shared_ptr<ParamSurface> sf = sub_models2[0]->getSurface(ki);
      sf->writeStandardHeader(of5);
      sf->write(of5);
    }
  
  // Create class for splitting a given face set into a face set where
  // all associated surfaces has 4 boundaries
  std::cout << "Replace current surfaces with 4-sided surfaces"  << std::endl;
  RegularizeFaceSet reg(sub_models2[0]);

  // Perform splitting
  // Note that this functionality is under developement and unexpected results
  // may occur. 
  shared_ptr<SurfaceModel> model4 = reg.getRegularModel();
  nmb = model4->nmbEntities();
  std::ofstream of6(outfile6.c_str());
    for (ki=0; ki<nmb; ++ki)
      {
	shared_ptr<ParamSurface> sf = model4->getSurface(ki);
	sf->writeStandardHeader(of6);
	sf->write(of6);
      }

    // Replace by spline surfaces
    // Boundary curves are approximated by spline curves, or fetched from
    // an already created adjacent spline surface. Spline surfaces are 
    // initially created as Coons patches interpolating the boundary curves
    // then updated with respect a point set fetched from the initial 
    // trimmed surface
    // The gap tolerance is used in the curve and surface approximations
    std::cout << "Replace trimmed surfaces by spline surfaces"  << std::endl;
    model4->replaceRegularSurfaces();

    // Write current surface set to a file
    std::ofstream of7(outfile7.c_str());
     nmb = model4->nmbEntities();
    for (ki=0; ki<nmb; ++ki)
      {
	shared_ptr<ParamSurface> sf = model4->getSurface(ki);
	sf->writeStandardHeader(of7);
	sf->write(of7);
      }

    // Make sure that neighbouring surfaces have the same spline space
    // and coincidence of corresponding coefficients
    // First check if all surfaces are non-trimmed spline surfaces, i.e.
    // the previous opeation succeeded
    bool isOK = model4->allSplines();
    if (!isOK)
      {
	std::cout << "Not all surfaces are splines. Stopping computation" << std::endl;
	exit(-1);
      }
    
    // Ensure common spline spaces and corresponding coefficients
    std::cout << "Ensure common spline space"  << std::endl;
     model4->makeCommonSplineSpaces();

    // Write result to file
    std::ofstream of8(outfile8.c_str());
    for (ki=0; ki<nmb; ++ki)
      {
	shared_ptr<ParamSurface> sf = model4->getSurface(ki);
	sf->writeStandardHeader(of8);
	sf->write(of8);
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
  std::ofstream of9(outfile9.c_str());
  for (ki=0; ki<nmb_vol; ++ki)
    {
      shared_ptr<ParamVolume> vol = volmodel->getVolume(ki);
       vol->writeStandardHeader(of9);
       vol->write(of9);
     }
  
}

  
		      
