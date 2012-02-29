//===========================================================================
//                                                                           
// File: createVolumeBoundaries.C
//                                                                           
//===========================================================================

#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/geometry/Sphere.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SweepSurfaceCreator.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/Sphere.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/geometry/LoopUtils.h"
#include "GoTools/geometry/SurfaceInterpolator.h"
#include "GoTools/creators/ApproxCurve.h"
#include "GoTools/creators/ApproxSurf.h"
#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/creators/LoftSurfaceCreator.h"
#include "GoTools/compositemodel/AdaptSurface.h"
#include "GoTools/geometry/SplineInterpolator.h"
#include <fstream>

using namespace Go;
using std::cout;



//===========================================================================
///                                                                           
/// Description:
///  
/// This programs creates a set of B-spline surfaces intended as the boundary
/// surfaces for a spline volume. A number of different methods are used in the
/// surface construction. Thus, this expample program illustrates some of the
/// possibilities for surface construction.
///
/// Input/Output:
///
/// Input to the geometry construction is partly hardcoded, partly fetched
/// from given data files.
/// Current surfaces written to g2-files as we go along. The final surface set
/// is written to the file data/volume_boundaries.g2
///                                                                           
//===========================================================================
int main( int argc, char* argv[] )
{
  // Input files
  std::string infile1("data/point_sequence.g2");
  std::string infile2("data/point_set.g2");
  
  // Prepare for output files
  std::string outfile1("data/bdsf1.g2");
  std::string outfile2("data/bdsf2.g2");
  std::string outfile3("data/bdsf3.g2");
  std::string outfile4("data/bdsf4.g2");
  std::string outfile5("data/bdsf5.g2");
  std::string outfile6("data/bdsf6.g2");
  std::string outfile7("data/volume_boundaries.g2");
 
  // Read file containing a point sequence
  std::ifstream if1(infile1.c_str());
  ObjectHeader header;
  header.read(if1);
  PointCloud<3> points;   // Points in 3-dimensional space
  points.read(if1);

  // Fetch raw data
  int nmb_pnts = points.numPoints();
  int dim = points.dimension();
  vector<double> pts(points.rawData(), points.rawData()+dim*nmb_pnts);

  std::cout << "Approximating point sequence" << std::endl;

  // Parameterize. The point sequence has some changes in curvature.
  // Thus, we prefer centripetal parameterization
  vector<double> param(nmb_pnts);
  param[0] = 0.0;
  for (int ki=1; ki<nmb_pnts; ++ki)
    {
      double dist = sqrt(Utils::distance_squared(pts.begin()+(ki-1)*dim,
					  pts.begin()+ki*dim, 
					  pts.begin()+ki*dim));
      param[ki] = param[ki-1]+sqrt(dist);
    }

  // Create curve approximation
  int max_iter = 10;   // Maximum number of iterations performing parameter 
  // iterations and refinement of the initial spline space
  int order = 4;       // Polynomial order of resulting curve
  double maxdist, avdist;  // Maximum and average approximation error
  double epsge = 1.0e-3;

  // Prepare for approximation
  ApproxCurve approx(pts, param, dim, epsge, order, order);
  approx.setSmooth(1.0e-6);  // Apply medium size weights to smoothing

  // Perform approximation and fetch result
  shared_ptr<SplineCurve> crv = approx.getApproxCurve(maxdist, avdist, max_iter);
  std::cout << "Curve approximation, maximum distance: " << maxdist;
  std::cout << ", average distace: " << avdist << std::endl;
  
  std::cout << "1. boundary surface is a swept surface" << std::endl;

  // Make first boundary surface by sweeping this curve in the y-direction
  Point pos = crv->ParamCurve::point(crv->startparam());  // Point on curve
  Point y_axis(0.0, 1.0, 0.0);   // Sweep direction
  double len = 3.0;            // Length of surface in sweep direction
  shared_ptr<SplineCurve> crv2(new SplineCurve(pos, pos+len*y_axis));  // Linear spline curve
  shared_ptr<SplineSurface> bdsf1 =
    shared_ptr<SplineSurface>(SweepSurfaceCreator::linearSweptSurface(*crv,
								      *crv2,
								      pos));
  
  // Write current result to file
    std::ofstream of1(outfile1.c_str());
    crv->writeStandardHeader(of1);
    crv->write(of1);
    bdsf1->writeStandardHeader(of1);
    bdsf1->write(of1);
 
    std::cout << "Second boundary surface" << std::endl;

    // Create the second (opposite) boundary surface as a part of a sphere
    // First create sphere
    Point centre(2.0, 1.5, 2.0);
    double radius = 3;
    Point x_axis(1.0, 0.0, 0.0);
    Point up_axis(0.0, 1.0, 0.0);
    shared_ptr<Sphere> sphere(new Sphere(radius, centre, up_axis, x_axis));

    // Limit the sphere
    sphere->setParameterBounds(1.25*M_PI, -0.2*M_PI, 1.75*M_PI, 0.25*M_PI);

    // Create NURBS representation of the sphere
    shared_ptr<SplineSurface> sphere_sf(sphere->createSplineSurface());

    // The volume will later be created as a Coons volume. Rational boundary
    // surfaces are not allowed in this construction. Thus, we need to create
    // a non-rational approximation. 
    // Start by creating a surface representing the initial spline space
    // The initial approximating surface is set to have the same spline space
    // as the rational surface. It also have the same coefficients, but the
    // weights are excluded. The coefficients are not used in the construction
    shared_ptr<SplineSurface> sphere_sf2(new SplineSurface(sphere_sf->basis_u(),
							   sphere_sf->basis_v(),
							   sphere_sf->coefs_begin(),
							   sphere_sf->dimension(),
							   false));

    // The NURBS representation of the spherical surface is quadratic. We want
    // a cubic surface.
    sphere_sf2->raiseOrder(1, 1);  // Raise the order with one in both parameter directions

    // Modify surface to approximate the sphere using the same tolerance as
    // for the point approxamation
    shared_ptr<SplineSurface> bdsf2 = AdaptSurface::approxInSplineSpace(sphere_sf, 
									sphere_sf2, 
									epsge);

    // Write to file
    std::ofstream of2(outfile2.c_str());
    bdsf2->writeStandardHeader(of2);
    bdsf2->write(of2);

    std::cout << "The third boundary surface" << std::endl;

    // Make the third boundary surface by lofting
    // First fetch the corner points in the maximum y-direction for both
    // existing surfaces.
    Point pt1sf1 = bdsf1->ParamSurface::point(bdsf1->startparam_u(),
					      bdsf1->endparam_v());
    Point pt2sf1 = bdsf1->ParamSurface::point(bdsf1->endparam_u(),
					      bdsf1->endparam_v());
    Point pt1sf2 = bdsf2->ParamSurface::point(bdsf2->startparam_u(),
					      bdsf2->endparam_v());
    Point pt2sf2 = bdsf2->ParamSurface::point(bdsf2->endparam_u(),
					      bdsf2->endparam_v());
    
    // Create intermediate curve by Hermite interpolation. The curve
    // interpolates the midpoint between corresponding corner points,
    // and is assigned a tangent pointing along the y-axis
    SplineInterpolator interpol;  // Create an empty SplineInterpolator.
    // Set tangent conditions in endpoints
    double len2 = 1.2;
    interpol.setHermiteConditions(len2*y_axis, -len2*y_axis);
    vector<double> midcvpos;
    vector<double> par(2);
    Point mid1 = 0.5*(pt1sf1+pt1sf2) + 0.5*y_axis;
    Point mid2 = 0.5*(pt2sf1+pt2sf2) + 0.3*y_axis;
    midcvpos.insert(midcvpos.end(), mid1.begin(), mid1.end());
    midcvpos.insert(midcvpos.end(), mid2.begin(), mid2.end());
    par[0] = 0.0;
    par[1] = mid1.dist(mid2);
    // Create an empty spline curve and create content by interpolation
    shared_ptr<SplineCurve> mid_cv(new SplineCurve());
    mid_cv->interpolate(interpol, 2, 3, &par[0], &midcvpos[0]);

    // Collect curves
    vector<shared_ptr<SplineCurve> > loft_cvs1(3);
    loft_cvs1[0] = shared_ptr<SplineCurve>(bdsf1->constParamCurve(bdsf1->endparam_v(),
								  true));
    loft_cvs1[1] = mid_cv;
    loft_cvs1[2] = shared_ptr<SplineCurve>(bdsf2->constParamCurve(bdsf2->endparam_v(),
								  true));
    // Loft surface.
    shared_ptr<SplineSurface> bdsf3(LoftSurfaceCreator::loftSurface(loft_cvs1, 3));
    // Write to file
     std::ofstream of3(outfile3.c_str());
     loft_cvs1[0]->writeStandardHeader(of3);
     loft_cvs1[0]->write(of3);
     loft_cvs1[1]->writeStandardHeader(of3);
     loft_cvs1[1]->write(of3);
     loft_cvs1[2]->writeStandardHeader(of3);
     loft_cvs1[2]->write(of3);
     bdsf3->writeStandardHeader(of3);
     bdsf3->write(of3);
   
     std::cout << "The fourth boundary surface" << std::endl;

    // Fourth surface, again by loft
    // Collect curves
    vector<shared_ptr<SplineCurve> > loft_cvs2(2);
    loft_cvs2[0] = shared_ptr<SplineCurve>(bdsf1->constParamCurve(bdsf1->startparam_v(),
								  true));
    loft_cvs2[1] = shared_ptr<SplineCurve>(bdsf2->constParamCurve(bdsf2->startparam_v(),
								  true));
    // Loft surface.
    shared_ptr<SplineSurface> bdsf4(LoftSurfaceCreator::loftSurface(loft_cvs2, 2));

     std::ofstream of4(outfile4.c_str());
     loft_cvs2[0]->writeStandardHeader(of4);
     loft_cvs2[0]->write(of4);
     loft_cvs2[1]->writeStandardHeader(of4);
     loft_cvs2[1]->write(of4);
     bdsf4->writeStandardHeader(of4);
     bdsf4->write(of4);
      
     std::cout << "The fifth boundary surface" << std::endl;

     // Start by making a Coons patch surface interpolating the boundary curves
     // of all adjacent surfaces. First fetch boundary curves
     vector<shared_ptr<ParamCurve> > coons_bd1(4);
     coons_bd1[0] = shared_ptr<ParamCurve>(bdsf1->constParamCurve(bdsf1->startparam_u(),
							      false));
     coons_bd1[1] = shared_ptr<ParamCurve>(bdsf3->constParamCurve(bdsf3->startparam_u(),
								   false));
     coons_bd1[2] = shared_ptr<ParamCurve>(bdsf2->constParamCurve(bdsf2->startparam_u(),
								   false));
     coons_bd1[3] = shared_ptr<ParamCurve>(bdsf4->constParamCurve(bdsf4->startparam_u(),
								   false));

     // The Coons patch constructor assumes ccw orientation of the curves
     bool reorganized = LoopUtils::makeLoopCCW(coons_bd1, epsge);
     // epsge = tolerance for distance between corresponding curve ends
     if (!reorganized)
       std::cout << "Curve loop not continuous" << std::endl;

     // Create surface
     CurveLoop loop1(coons_bd1, epsge);  // Organize curves in loop
     shared_ptr<SplineSurface> bdsf5(CoonsPatchGen::createCoonsPatch(loop1));

     // Modify surface by creating a fake field and add it to the surface.
     // First make a simple function representing the length distribution of 
     // the vector field 
     double knots1[2];
     double knots2[2];
     knots1[0] = bdsf5->startparam_u();
     knots1[1] = bdsf5->endparam_u();
     knots2[0] = bdsf5->startparam_v();
     knots2[1] = bdsf5->endparam_v();
     double coef = 0.0;
     // Make initial distribution surface (constant in both parameter directions)
     shared_ptr<SplineSurface> tmp_len(new SplineSurface(1, 1, 1, 1, &knots1[0],
							 &knots2[0], &coef, 1));

     // Increase the degree to quadratic
     tmp_len->raiseOrder(2, 2);

     // Modify the coefficient in the inner of the surface
     vector<double>::iterator coefs_it = tmp_len->coefs_begin();
     coefs_it[4] = 2.0;

     // Create field by regular interpolation of the surface normals multiplied
     // with the length distribution
     // Evaluate in the greville parameters
    // Fetch Greville parameters
    BsplineBasis basis_u = bdsf5->basis_u();
    int nmb_u = basis_u.numCoefs();
    BsplineBasis basis_v = bdsf5->basis_v();
    int nmb_v = basis_v.numCoefs();
    vector<double> reg_points;
    vector<double> par_u(nmb_u);
    vector<double> par_v(nmb_v);
    for (int kr=0; kr<nmb_u; ++kr)
      par_u[kr] = basis_u.grevilleParameter(kr);
    for (int kr=0; kr<nmb_v; ++kr)
      par_v[kr] = basis_v.grevilleParameter(kr);
    for (int kr=0; kr<nmb_v; ++kr)
      for (int kh=0; kh<nmb_u; ++kh)
	{
	  Point norm;
	  bdsf5->normal(norm, par_u[kh], par_v[kr]);
	  Point len = tmp_len->ParamSurface::point(par_u[kh], par_v[kr]);
	  norm.normalize();
	  norm *= len[0];
	  reg_points.insert(reg_points.end(), norm.begin(), norm.end());
	}


    // Regular interpolation
    vector<double> dummy_wgts;
    shared_ptr<SplineSurface> field(SurfaceInterpolator::regularInterpolation(basis_u,
									      basis_v,
									      par_u,
									      par_v,
									      reg_points,
									      bdsf5->dimension(),
									      false,
									      dummy_wgts));

    // Add field to boundary surface
    bdsf5->add(field.get());

     // Write to file
     std::ofstream of5(outfile5.c_str());
     bdsf5->writeStandardHeader(of5);
     bdsf5->write(of5);

     std::cout << "Last boundary surface" << std::endl;
     
     // The last boundary surface is created as a Coons surface and
     // modified by letting this base surface approximate a given point set
          // Start by making a Coons patch surface interpolating the boundary curves
     // of all adjacent surfaces. First fetch boundary curves
     vector<shared_ptr<ParamCurve> > coons_bd2(4);
     coons_bd2[0] = shared_ptr<ParamCurve>(bdsf1->constParamCurve(bdsf1->endparam_u(),
							      false));
     coons_bd2[1] = shared_ptr<ParamCurve>(bdsf3->constParamCurve(bdsf3->endparam_u(),
								   false));
     coons_bd2[2] = shared_ptr<ParamCurve>(bdsf2->constParamCurve(bdsf2->endparam_u(),
								   false));
     coons_bd2[3] = shared_ptr<ParamCurve>(bdsf4->constParamCurve(bdsf4->endparam_u(),
								   false));

     // The Coons patch constructor assumes ccw orientation of the curves
     reorganized = LoopUtils::makeLoopCCW(coons_bd2, epsge);
     // epsge = tolerance for distance between corresponding curve ends
     if (!reorganized)
       std::cout << "Curve loop not continuous" << std::endl;

     // Create surface
     CurveLoop loop2(coons_bd2, epsge);  // Organize curves in loop
     shared_ptr<SplineSurface> bdsf6(CoonsPatchGen::createCoonsPatch(loop2));

     // Write to file
     std::ofstream of6(outfile6.c_str());
     bdsf6->writeStandardHeader(of6);
     bdsf6->write(of6);

     // Read point set from file
     std::ifstream if2(infile2.c_str());
     header.read(if2);
     PointCloud<3> points2;   // Points in 3-dimensional space
     points2.read(if2);

     // Parameterize the points by projecting them onto the base surface
     // bdsf6
     vector<double> point_set(points2.rawData(), 
			      points2.rawData()+points2.numPoints()*points2.dimension());
     vector<double> parvals;
     shared_ptr<ParamSurface> tmpbd = bdsf6;
     SurfaceTools::parameterizeByBaseSurf(*tmpbd, point_set, parvals); // The function is found
     // in gotools-core/geometry/SurfaceTools

     // Modify base surface
     // First create engine
     double eps2 = 0.1;   // The points do not origine from a smooth surface.
     // Thus it is impossible to approximate them within a tight tolerance.
     ApproxSurf approxsf(bdsf6, point_set, parvals, points2.dimension(),
			 eps2);
     approxsf.setFixBoundary(1);   // Do not modify the surface boundaries
     approxsf.setSmoothingWeight(0.01);
		      
     // Perform modification
     max_iter = 3;  // Do not expect a tight approximation, but allow for some
     // refinement of the spline space
     int nmb_out;
     shared_ptr<SplineSurface> bdsf7 = approxsf.getApproxSurf(maxdist, avdist, 
							      nmb_out, max_iter);

      bdsf7->writeStandardHeader(of6);
      bdsf7->write(of6);
    
      // Collect all boundary surfaces. The Coons volume creator expects
      // surfaces in the order umin, umax, vmin, vmax, wmin, wmax, but does
      // not expect the surfaces to be oriented in a constistent way.
     std::ofstream of7(outfile7.c_str());
     bdsf1->writeStandardHeader(of7);
     bdsf1->write(of7);
     bdsf2->writeStandardHeader(of7);
     bdsf2->write(of7);
     bdsf3->writeStandardHeader(of7);
     bdsf3->write(of7);
     bdsf4->writeStandardHeader(of7);
     bdsf4->write(of7);
     bdsf5->writeStandardHeader(of7);
     bdsf5->write(of7);
     bdsf7->writeStandardHeader(of7);
     bdsf7->write(of7);
      
}

