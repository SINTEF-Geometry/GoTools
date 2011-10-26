#include <vector>
using std::vector;
using std::pair;
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/geometry/SISL_code.h"
 
//***************************************************************************
//
// Implementation file of the free function intersectcurves defined in
// GoIntersections.h/
//
//***************************************************************************

namespace Go
{

void intersectcurves(SplineCurve* cv1, SplineCurve* cv2, double epsge,
		     vector<std::pair<double,double> >& intersections)
  //************************************************************************
  // 
  // Intersect two spline curves. Collect intersection parameters.
  //
  //***********************************************************************
{

  // Make sisl curves and call sisl.
  SISLCurve *pc1 = Curve2SISL(*cv1, false);
  SISLCurve *pc2 = Curve2SISL(*cv2, false);

  int kntrack = 0;
  int trackflag = 0;  // Do not make tracks.
  SISLTrack **track =0;
  int knpt=0, kncrv=0;
  double *par1=0, *par2=0;
  int *pretop = 0;
  SISLIntcurve **intcrvs = 0;
  int stat = 0;
  sh1857(pc1, pc2, 0.0, epsge, trackflag, &kntrack, &track,
	 &knpt, &par1, &par2, &pretop, &kncrv, &intcrvs, &stat);

  ALWAYS_ERROR_IF(stat<0,"Error in intersection, code: " << stat);


  // Remember intersections points. The intersection curves are
  // skipped.
  int ki;
  for (ki=0; ki<knpt; ki++)
    {
      intersections.push_back(std::make_pair(par1[ki],par2[ki]));
    }

  if (kncrv > 0)
    freeIntcrvlist(intcrvs, kncrv);

  if (par1 != 0) free(par1);
  if (par2 != 0) free(par2);
  if (pc1 != 0) freeCurve(pc1);
  if (pc2 != 0) freeCurve(pc2);
}



} // namespace Go  
