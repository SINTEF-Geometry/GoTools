#include "GoTools/geometry/GeometryTools.h"

#include "GoTools/geometry/Line.h"

//***************************************************************************
//
// Implementation file of the free function projectCurve defined in
// GeometryTools.h/
//
//***************************************************************************

namespace Go
{

shared_ptr<ParamCurve> projectCurve(shared_ptr<ParamCurve> incurve,
				    const Point& normal,
				    bool planar)
{
  int dim = incurve->dimension();
  ALWAYS_ERROR_IF(dim != 2 && dim != 3, "Dimension must be 2 or 3");
  ALWAYS_ERROR_IF(dim!=2 && dim!=normal.dimension(),
		  "Error in dimension.");

  if (incurve->instanceType() == Class_Line)
    {
      shared_ptr<Line> line = dynamic_pointer_cast<Line, ParamCurve>(incurve);
      MESSAGE_IF(!line->isBounded(), "Line not bounded.");

      shared_ptr<SplineCurve> spline_cv(line->geometryCurve());
      shared_ptr<SplineCurve> proj_cv =
	projectCurve(*spline_cv, normal, planar);

      Point from_pt = proj_cv->ParamCurve::point(line->startparam());
      Point to_pt = proj_cv->ParamCurve::point(line->endparam());
      Point dir_vec = to_pt - from_pt;
      shared_ptr<ParamCurve> proj_line =
	shared_ptr<Line>(new Line(from_pt, dir_vec));
      proj_line->setParameterInterval(line->startparam(), line->endparam());
      return proj_line;
    }
  else if (incurve->instanceType() == Class_SplineCurve)
    {
      shared_ptr<SplineCurve> spline_cv =
	dynamic_pointer_cast<SplineCurve, ParamCurve>(incurve);
      return projectCurve(*spline_cv, normal, planar);
    }
  else {
      shared_ptr<SplineCurve> spline_cv(incurve->geometryCurve());
      if (spline_cv.get() == NULL) {
	  MESSAGE("Unexpected curve type!");
      } else {
	  return projectCurve(*spline_cv, normal, planar);
      }
  }

  return shared_ptr<ParamCurve>();
}


shared_ptr<SplineCurve> projectCurve(const SplineCurve& incurve,
				     const Point& normal,
				     bool planar)
  //*************************************************************************
  //
  // Project a SplineCurve into a given plane. The curve can
  // be returned either as a 3D curve lying in the plane or as
  // a 2D curve. In the latter case, the coordinate system is
  // rotated such that the given plane normal coincides with the
  // z-axis.
  //
  //*************************************************************************
{
  int dim = incurve.dimension();
  ALWAYS_ERROR_IF(dim != 2 && dim != 3, "Dimension must be 2 or 3");
  ALWAYS_ERROR_IF(dim!=2 && dim!=normal.dimension(),
		  "Error in dimension.");

  // Make sure that the given normal is normalized.
  Point norm = normal;
  norm.normalize();

  // Fetch the coefficients of the input curve
  int nmbcoef = incurve.numCoefs();
  bool rational = incurve.rational();
  std::vector<double>::const_iterator coefs;
  if (rational)
    coefs = incurve.rcoefs_begin();
  else
    coefs = incurve.coefs_begin();
  int dim1 = dim + (rational);

  // Find the position of the plane which approximates
  // the input curve best. First compute the mean distance
  // from the curve coefficients to the plane with the given 
  // normal passing through origo.
  double meandist = 0.0, dist;
  int ki, kj;
  std::vector<double>::const_iterator co;
  std::vector<double>::iterator co2, co3;
  for (ki=0, co=coefs; ki<nmbcoef; ki++, co+=dim1)
    {
      for (kj=0, dist=0.0; kj<dim; kj++)
	dist += co[kj]*norm[kj];
      meandist += dist;
    }
  meandist /= (double)nmbcoef;

  // The point in the plane
  Point point(meandist*norm);

  // Project the coefficients into the found plane
  std::vector<double> coef2(dim1*nmbcoef);
  for (ki=0, co=coefs, co2=coef2.begin(); ki<nmbcoef; ki++, co+=dim1, co2+=dim1)
    {
      Point pp(co[0],co[1],co[2]);
      double td = (point - pp)*norm;
      for (kj=0; kj<dim; kj++)
	co2[kj] = co[kj] - td*norm[kj];
      if (rational)
	co2[dim] = co[dim];   // Keep the weight
    }

  SplineCurve *planecrv;
  if (planar)
    {
      // Rotate coordinate system such that the z-axis coincides
      // with the given plane normal, and express the coefficients in 
      // this coordinate system. Remove z in each coefficient.
      std::vector<double> coef3((dim1-1)*nmbcoef);
      double tdum;
      for (co3=coef3.begin(), co2=coef2.begin(), ki=0; ki<nmbcoef; 
	   ki++, co3+=(dim1-1), co2+=dim1)
	{
	  co3[0] = norm[2]*co2[0] + norm[0]*co2[1] + norm[1]*co2[2];
	  co3[1] = norm[1]*co2[0] + norm[2]*co2[1] + norm[0]*co2[2];
	  tdum = norm[0]*co2[0]+norm[1]*co2[1]+norm[2]*co2[2];
	  if (rational)
	    co3[2] = co2[dim];
	}
      
      planecrv = new SplineCurve(nmbcoef, incurve.order(),
				   incurve.basis().begin(),
				   coef3.begin(), 2, rational);
    }
  else
    planecrv = new SplineCurve(nmbcoef, incurve.order(),
			       incurve.basis().begin(),
			       coef2.begin(), 3, rational);

  return shared_ptr<SplineCurve>(planecrv);
}


} // namespace Go
