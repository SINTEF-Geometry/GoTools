#include <vector>
using std::vector;
using std::pair;
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/utils/Point.h"
#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/geometry/SISL_code.h"

#ifdef __BORLANDC__
using std::free;
#endif

//***************************************************************************
//
// Implementation file of the free function intersectCurveSurf defined in
// GoIntersections.h/
//
//***************************************************************************

using std::pair;
using std::vector;

namespace Go
{

 void intersectCurveSurf(const SplineCurve *cv, const SplineSurface *sf,
			 double epsge, 
			 vector<pair<double, Point> >& int_pts,
			 vector<int>& pretopology,
			 vector<pair<pair<double,Point>, 
			 pair<double,Point> > >& int_crvs)
 {
   SISLSurf* sislsf = GoSurf2SISL(*sf, false);
   SISLCurve* sislcv = Curve2SISL(*cv, false);
   int kntrack = 0;
   int trackflag = 0;  // Do not make tracks.
   SISLTrack **track =0;
   int knpt=0, kncrv=0;
   double *par1=0, *par2=0;
   int *pretop = 0;
   SISLIntcurve **vcrv = 0;
   int stat = 0;

   sh1858(sislsf, sislcv, 0.0, epsge, trackflag, &kntrack, &track,
	  &knpt, &par1, &par2, &pretop, &kncrv, &vcrv, &stat);
   ALWAYS_ERROR_IF(stat<0,"Error in intersection, code: " << stat);

   // Remember intersections points. 
   int ki;
   for (ki=0; ki<knpt; ki++)
     {
       int_pts.push_back(std::make_pair(par2[ki],
					Point(par1[2*ki],par1[2*ki+1])));
       pretopology.insert(pretopology.end(), pretop+4*ki+2, pretop+4*(ki+1));
       pretopology.insert(pretopology.end(), pretop+4*ki, pretop+4*ki+2);
     }

   // Remember intersection curves
   for (ki=0; ki<kncrv; ++ki)
     {
       int nmb_pt = vcrv[ki]->ipoint;
       Point par1 = Point(vcrv[ki]->epar1[0],vcrv[ki]->epar1[1]);
       Point par2 = Point(vcrv[ki]->epar1[2*(nmb_pt-1)],
			  vcrv[ki]->epar1[2*nmb_pt-1]);
       int_crvs.push_back(std::make_pair(std::make_pair(vcrv[ki]->epar2[0], 
							par1),
					 std::make_pair(vcrv[ki]->epar2[nmb_pt-1],
							par2)));
     }

   if (kncrv > 0)
     freeIntcrvlist(vcrv, kncrv);

   if (par1 != NULL) free(par1);
   if (par2 != NULL) free(par2);
   if (sislsf != NULL) freeSurf(sislsf);
   if (sislcv != NULL) freeCurve(sislcv);
   if (pretop != NULL) free(pretop);
 }

} // namespace Go  
