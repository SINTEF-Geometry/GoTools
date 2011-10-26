//===========================================================================
//                                                                           
// File: VolumeSpaceCurve
//                                                                           
// Created: Nov 10, 2010
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/trivariate/VolumeSpaceCurve.h"
#include "GoTools/trivariate/ParamVolume.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/utils/Array.h"
#include "GoTools/utils/Volumes.h"

using namespace Go;
using std::vector;
using std::shared_ptr;

//===========================================================================
VolumeSpaceCurve::VolumeSpaceCurve(shared_ptr<ParamVolume> vol, 
				   shared_ptr<ParamCurve> crv)
  : vol_(vol), crv_(crv)
//===========================================================================
{
}

//===========================================================================
VolumeSpaceCurve::~VolumeSpaceCurve()
//===========================================================================
{
}
 
//===========================================================================
Point VolumeSpaceCurve::eval(double t) const
//===========================================================================
{
  vector<Point> result(1);
  eval(t, 0, &result[0]);
  return result[0];
}

//===========================================================================
void VolumeSpaceCurve::eval(double t, int n, Point der[]) const
//===========================================================================
{
  if (n > 1)
    n = 1;
  if (n < 0)
    n = 0;

  // Evaluate curve
  vector<Point> cv_der(n+1);
  crv_->point(cv_der, t, n);

  // Evaluate volume
  vector<Point> vol_der(1+n*3);
  vol_->point(vol_der, cv_der[0][0], cv_der[0][1], cv_der[0][2], n);

  der[0] = vol_der[0];
  if (n == 1)
    {
      der[1] = cv_der[1][0]*vol_der[1] + cv_der[1][1]*vol_der[2] +
	cv_der[1][2]*vol_der[3];
    }
}

//===========================================================================
double VolumeSpaceCurve::start() const
//===========================================================================
{
  return crv_->startparam();
}

//===========================================================================
double VolumeSpaceCurve::end() const
//===========================================================================
{
  return crv_->endparam();
}

//===========================================================================
int VolumeSpaceCurve::dim() const
//===========================================================================
{
// Trivariate parameter area
  return 3;
}

//===========================================================================
bool VolumeSpaceCurve::approximationOK(double par, Point approxpos,
					   double tol1, double tol2) const
//===========================================================================
{
  Point pos = eval(par);
  double dist = pos.dist(approxpos);
  return (dist <= tol1);
}

