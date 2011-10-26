#include <vector>
using std::vector;
using std::pair;
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/geometry/SISL_code.h"

#ifdef __BORLANDC__
using std::free;
#endif

//***************************************************************************
//
// Implementation file of the free function intersect2Dcurves defined in
// GoIntersections.h/
//
//***************************************************************************

using std::pair;
using std::vector;

namespace Go
{

void intersectCurvePoint(const ParamCurve* crv, Point pnt, double epsge,
			 vector<double>& intersections, 
			 vector<pair<double, double> >& int_crvs)
  //************************************************************************
  // 
  // Intersect a curve with a point
  //
  //***********************************************************************
{
  // First make sure that the curve is a spline curve
  ALWAYS_ERROR_IF(crv->instanceType() != Class_SplineCurve,
		  " Intersection involving general parametric curves is not implemented.");

  // Make sisl curve and call sisl.
  SISLCurve *pc = Curve2SISL(*(dynamic_cast<const SplineCurve*>(crv)), false);

  int knpt=0, kncrv=0;
  double *par=0;
  SISLIntcurve **vcrv = 0;
  int stat = 0;
  s1871(pc, pnt.begin(), pnt.size(), epsge, &knpt, &par, &kncrv, &vcrv, &stat);
  ALWAYS_ERROR_IF(stat<0,"Error in intersection, code: " << stat);


  // Remember intersections points. 
  if (knpt > 0)
    intersections.insert(intersections.end(), par, par+knpt);

  // Remember intersection curves
  for (int ki=0; ki<kncrv; ++ki)
    int_crvs.push_back(std::make_pair(vcrv[ki]->epar1[0], vcrv[ki]->epar1[vcrv[ki]->ipoint-1]));

  if (kncrv > 0)
    freeIntcrvlist(vcrv, kncrv);

  if (par != 0) free(par);
  if (pc) freeCurve(pc);
}

void intersect2Dcurves(const ParamCurve* cv1, const ParamCurve* cv2, double epsge,
		       vector<pair<double,double> >& intersections,
		       vector<int>& pretopology,
		       vector<pair<pair<double,double>, pair<double,double> > >& int_crvs)
  //************************************************************************
  // 
  // Intersect two 2D spline curves. Collect intersection parameters
  // and pretopology information.
  //
  //***********************************************************************
{

  // First make sure that the curves are spline curves.
  ALWAYS_ERROR_IF(cv1->instanceType() != Class_SplineCurve ||
  	      cv2->instanceType() != Class_SplineCurve,
		  " Intersection between general parametric curves is not implemented.");

  MESSAGE_IF(cv1->dimension() != 2,
		"Dimension different from 2, pretopology not reliable.");

  // Make sisl curves and call sisl.
  SISLCurve *pc1 = Curve2SISL(*(dynamic_cast<const SplineCurve*>(cv1)), false);
  SISLCurve *pc2 = Curve2SISL(*(dynamic_cast<const SplineCurve*>(cv2)), false);

  int kntrack = 0;
  int trackflag = 0;  // Do not make tracks.
  SISLTrack **track =0;
  int knpt=0, kncrv=0;
  double *par1=0, *par2=0;
  int *pretop = 0;
  SISLIntcurve **vcrv = 0;
  int stat = 0;
  sh1857(pc1, pc2, 0.0, epsge, trackflag, &kntrack, &track,
	 &knpt, &par1, &par2, &pretop, &kncrv, &vcrv, &stat);

  ALWAYS_ERROR_IF(stat<0,"Error in intersection, code: " << stat);


  // Remember intersections points. 
  int ki;
  for (ki=0; ki<knpt; ki++)
    {
      intersections.push_back(std::make_pair(par1[ki],par2[ki]));
      pretopology.insert(pretopology.end(), pretop+4*ki, pretop+4*(ki+1));
    }

  // Remember intersection curves
  for (ki=0; ki<kncrv; ++ki)
    int_crvs.push_back(std::make_pair(std::make_pair(vcrv[ki]->epar1[0], vcrv[ki]->epar2[0]), 
				      std::make_pair(vcrv[ki]->epar1[vcrv[ki]->ipoint-1],
						     vcrv[ki]->epar2[vcrv[ki]->ipoint-1])));

  if (kncrv > 0)
    freeIntcrvlist(vcrv, kncrv);

  if (par1 != 0) free(par1);
  if (par2 != 0) free(par2);
  if (pc1 != 0) freeCurve(pc1);
  if (pc2 != 0) freeCurve(pc2);
  if (pretop != 0) free(pretop);
}



} // namespace Go  
