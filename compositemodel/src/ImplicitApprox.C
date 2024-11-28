/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

#include "GoTools/compositemodel/ImplicitApprox.h"
#include "GoTools/compositemodel/RevEngPoint.h"
#include "GoTools/implicitization/ImplicitizePointCloudAlgo.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/LineCloud.h"
#include "GoTools/geometry/SweepSurfaceCreator.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/ftPlane.h"
#include "GoTools/compositemodel/ftCurve.h"
#include "GoTools/implicitization/BernsteinPoly.h"
#include "newmat.h"
#include "newmatap.h"
#include "sisl.h"
#include <fstream>

using namespace Go;
using std::vector;

//===========================================================================
ImplicitApprox::ImplicitApprox()
  : eps_(1.0e-12)
//===========================================================================
{
}

//===========================================================================
ImplicitApprox::~ImplicitApprox()
//===========================================================================
{
}

//===========================================================================
void ImplicitApprox::approx(vector<RevEngPoint*> points, int degree)
//===========================================================================
{
  // Extract xyz values
  vector<Vector3D> xyz(points.size());
  for (size_t ki=0; ki<points.size(); ++ki)
    xyz[ki] = points[ki]->getPoint();

  PointCloud3D pointset(xyz);

  // Implicitize
  degree_ = degree;
  ImplicitizePointCloudAlgo implicitize(pointset, degree);
  implicitize.perform();
  
  // Get result
  implicitize.getResultData(implicit_, bc_, sigma_min_);
  
  // Differentiate
  Vector4D bdir1(1.0, 0.0, 0.0, 0.0);
  Vector4D bdir2(0.0, 1.0, 0.0, 0.0);
  Vector4D bdir3(0.0, 0.0, 1.0, 0.0);
  Vector4D bdir4(0.0, 0.0, 0.0, 1.0);
  implicit_.deriv(1, bdir1, deriv1_);
  implicit_.deriv(1, bdir2, deriv2_);
  implicit_.deriv(1, bdir3, deriv3_);
  implicit_.deriv(1, bdir4, deriv4_);

}

//===========================================================================
void ImplicitApprox::approxPoints(vector<Point> points, int degree)
//===========================================================================
{
  // Extract xyz values
  vector<Vector3D> xyz(points.size());
  for (size_t ki=0; ki<points.size(); ++ki)
    xyz[ki] = Vector3D(points[ki].begin());

  PointCloud3D pointset(xyz);

  // Implicitize
  degree_ = degree;
  ImplicitizePointCloudAlgo implicitize(pointset, degree);
  implicitize.perform();
  
  // Get result
  implicitize.getResultData(implicit_, bc_, sigma_min_);
  
  // Differentiate
  Vector4D bdir1(1.0, 0.0, 0.0, 0.0);
  Vector4D bdir2(0.0, 1.0, 0.0, 0.0);
  Vector4D bdir3(0.0, 0.0, 1.0, 0.0);
  Vector4D bdir4(0.0, 0.0, 0.0, 1.0);
  implicit_.deriv(1, bdir1, deriv1_);
  implicit_.deriv(1, bdir2, deriv2_);
  implicit_.deriv(1, bdir3, deriv3_);
  implicit_.deriv(1, bdir4, deriv4_);

}

//===========================================================================
void ImplicitApprox::approx(vector<pair<vector<RevEngPoint*>::iterator,
			    vector<RevEngPoint*>::iterator> >& points,
			    int degree)
//===========================================================================
{
  // Extract xyz values
  vector<Vector3D> xyz;
  for (size_t kj=0; kj<points.size(); ++kj)
    {
      for (auto it=points[kj].first; it!=points[kj].second; ++it)
	{
	  Vector3D curr = (*it)->getPoint();
	  xyz.push_back(curr);
	}
    }
  PointCloud3D pointset(xyz);

  // Implicitize
  degree_ = degree;
  ImplicitizePointCloudAlgo implicitize(pointset, degree);
  implicitize.perform();
  
  // Get result
  implicitize.getResultData(implicit_, bc_, sigma_min_);
  
  // Differentiate
  Vector4D bdir1(1.0, 0.0, 0.0, 0.0);
  Vector4D bdir2(0.0, 1.0, 0.0, 0.0);
  Vector4D bdir3(0.0, 0.0, 1.0, 0.0);
  Vector4D bdir4(0.0, 0.0, 0.0, 1.0);
  implicit_.deriv(1, bdir1, deriv1_);
  implicit_.deriv(1, bdir2, deriv2_);
  implicit_.deriv(1, bdir3, deriv3_);
  implicit_.deriv(1, bdir4, deriv4_);

}

//===========================================================================
void ImplicitApprox::evaluate(Point& pt, double& val, Point& grad)
//===========================================================================
{
  Vector3D xyz(pt[0], pt[1], pt[2]);
  Vector4D bary = bc_.cartToBary(xyz);
  val = implicit_(bary);
  double d1 = deriv1_(bary);
  double d2 = deriv2_(bary);
  double d3 = deriv3_(bary);
  double d4 = deriv4_(bary);
  Vector4D dv(d1,d2,d3,d4);
  Vector4D bary2 = bary+dv;
  Vector3D pt2 = bc_.baryToCart(bary2);
  Vector3D grad2 = pt2 - xyz;
  grad = Point(grad2[0], grad2[1], grad2[2]);
}

//===========================================================================
double ImplicitApprox::estimateDist(RevEngPoint* pt)
//===========================================================================
{
  Vector3D xyz = pt->getPoint();
  Vector4D bary = bc_.cartToBary(xyz);
  double dist0 = implicit_(bary);
  double d1 = deriv1_(bary);
  double d2 = deriv2_(bary);
  double d3 = deriv3_(bary);
  double d4 = deriv4_(bary);
  Vector4D dv(d1,d2,d3,d4);
  Vector4D bary2 = bary+dv;
  Vector3D pt2 = bc_.baryToCart(bary2);
  Vector3D grad = pt2 - xyz;
  double len = grad.length();
  double dist = (len > eps_) ? dist0/len : dist0;

  Point norm = pt->getLocFuncNormal();
  Vector3D norm2(norm[0], norm[1], norm[2]);
  norm2 *= 100;
  Vector3D xyz2 = xyz + norm2;
  Vector3D xyz3 = xyz - norm2;
  Vector4D bary3 = bc_.cartToBary(xyz2);
  Vector4D bary4 = bc_.cartToBary(xyz3);
  BernsteinPoly line = implicit_.pickLine(bary3, bary4);

  // Compute zeroes of bernstein polynomial
  // First make sisl curve
  int ik = degree_ + 1;
  vector<double> et(2*ik, 0.0);  // Knot vector of line curve
  for (int ki=0; ki<ik; ++ki)
    et[ik+ki] = 1.0;
  vector<double> ecoef(line.coefsBegin(), line.coefsEnd());
  SISLCurve *qc = newCurve(ik, ik, &et[0], &ecoef[0], 1, 1, 1);
  double zero = 0.0;
  
  // Intersect
  double eps = 1.0e-6;
  int kstat = 0;
  int kcrv=0, kpt=0;
  double *epar = 0;
  SISLIntcurve **intcv = 0;
  if (qc)
    s1871(qc, &zero, 1, eps, &kpt, &epar, &kcrv, &intcv, &kstat);
  if (qc)
    freeCurve(qc);
  
  // Compute cartesian points and curves associated with intersections
  double dd = std::numeric_limits<double>::max();
  for (int kr=0; kr<kpt; ++kr)
    {
      Vector4D barypt = (1.0 - epar[kr])*bary3 + epar[kr]*bary4;
      Vector3D pos = bc_.baryToCart(barypt);
      double dd2 = xyz.dist(pos);
      if (dd2 < dd)
	dd = dd2;
    }
  if (epar) free(epar);
  if (intcv) freeIntcrvlist(intcv, kcrv);
  return dd; //dist;
}

//===========================================================================
bool ImplicitApprox::projectPoint(Point point, Point dir,
				  Point& projpos, Point& normal)
//===========================================================================
{
  double len = 100.0;
  dir.normalize();
  
  Point xdir(1.0, 0.0, 0.0);
  Point ydir(0.0, 1.0, 0.0);
  Point zdir(0.0, 0.0, 1.0);
  double a1 = xdir.angle(dir);
  double a2 = ydir.angle(dir);
  double a3 = zdir.angle(dir);
  Point dir2;
  if (a1 > std::min(a2, a3))
    dir2 = xdir;
  else if (a2 > a3)
    dir2 = ydir;
  else
    dir2 = zdir;
  Point dir3 = dir%dir2;
  dir2 = dir%dir3;
  dir2.normalize();
  dir3.normalize();
  Point points[3];
  points[0] = point;
  points[1] = point + dir2;
  points[2] = point + dir3;

  Vector3D proj[3];
  int ka;
  for (ka=0; ka<3; ++ka)
    {
      Vector3D xyz(points[ka].begin());
      Point p1 = points[ka] - len*dir;
      Point p2 = points[ka] + len*dir;

      Vector3D cart1(p1.begin());
      Vector3D cart2(p2.begin());
      Vector4D bary1 = bc_.cartToBary(Vector3D(cart1[0], cart1[1], cart1[2]));
      Vector4D bary2 = bc_.cartToBary(Vector3D(cart2[0], cart2[1], cart2[2]));

      // Pick line
      BernsteinPoly line = implicit_.pickLine(bary1, bary2);

      // Compute zeroes of bernstein polynomial
      // First make sisl curve
      int ik = degree_ + 1;
      vector<double> et(2*ik, 0.0);  // Knot vector of line curve
      for (int ki=0; ki<ik; ++ki)
	et[ik+ki] = 1.0;
      vector<double> ecoef(line.coefsBegin(), line.coefsEnd());
      SISLCurve *qc = newCurve(ik, ik, &et[0], &ecoef[0], 1, 1, 1);
      double zero = 0.0;

      // Intersect
      double eps = 1.0e-6;
      int kstat = 0;
      int kcrv=0, kpt=0;
      double *epar = 0;
      SISLIntcurve **intcv = 0;
      if (qc)
	s1871(qc, &zero, 1, eps, &kpt, &epar, &kcrv, &intcv, &kstat);
      if (qc)
	freeCurve(qc);
      if (kpt == 0)
	return false;

      // Compute cartesian points and curves associated with intersections
      double mindist = std::numeric_limits<double>::max();
      for (int kr=0; kr<kpt; ++kr)
	{
	  Vector4D barypt = (1.0 - epar[kr])*bary1 + epar[kr]*bary2;
		
	  Vector3D pos = bc_.baryToCart(barypt);
	  double dist = pos.dist(xyz);
	  if (dist < mindist)
	    {
	      mindist = dist;
	      proj[ka] = pos;
	    }
	}
      if (epar) free(epar);
      if (intcv) freeIntcrvlist(intcv, kcrv);
      
    }

  projpos = Point(proj[0][0], proj[0][1], proj[0][2]);
  Point pt2(proj[1][0], proj[1][1], proj[1][2]);
  Point pt3(proj[2][0], proj[2][1], proj[2][2]);
  
  Point vec1 = pt2 - projpos;
  Point vec2 = pt3 - projpos;
  normal = vec1.cross(vec2);
  normal.normalize_checked();
  return true;
}

//===========================================================================
void ImplicitApprox::visualize(vector<RevEngPoint*> points, std::ostream& os)
//===========================================================================
{
  // View direction
  Point dir = points[0]->getLocFuncNormal();
  dir.normalize();

  BoundingBox bb(3);
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      Vector3D xyz = points[ki]->getPoint();
      bb.addUnionWith(Point(xyz[0], xyz[1], xyz[2]));
    }
  Point low = bb.low();
  Point high = bb.high();
  Point bmid = 0.5*(low + high);
  Point diag = high - low;
  double diaglen = diag.length();

  double gap = 1.0e-6;
  Point xdir(1.0, 0.0, 0.0);
  Point ydir(0.0, 1.0, 0.0);
  Point zdir(0.0, 0.0, 1.0);
  CompositeModelFactory factory(gap, gap, 10.0*gap, 0.01, 0.05);
  shared_ptr<SurfaceModel> boxmod(factory.createFromBox(low-2.0*diag, xdir, ydir, 
  							5*diag[0], 5*diag[1], 5*diag[2]));
    
    // Find the coordinate direction with the largest angle with the view direction
    double a1 = xdir.angle(dir);
    double a2 = ydir.angle(dir);
    double a3 = zdir.angle(dir);
    Point dir2;
    if (a1 > std::min(a2, a3))
      dir2 = xdir;
    else if (a2 > a3)
      dir2 = ydir;
    else
      dir2 = zdir;
    Point dir3 = dir%dir2;
    dir2 = dir%dir3;
    if (dir2*(high-low) < 0.0)
      dir2 *= -1.0;
    if (dir3*(high-low) < 0.0)
      dir3 *= -1.0;
    dir2.normalize();
    dir3.normalize();
    double len = low.dist(high);
    int nmb_sample = 100;
    shared_ptr<SplineCurve> cv1(new SplineCurve(bmid-len*dir2, 0.0, bmid+len*dir2, 1.0));
    shared_ptr<SplineCurve> cv2(new SplineCurve(bmid-len*dir3, 0.0, bmid+len*dir3, 1.0));
    SweepSurfaceCreator sweep;
    shared_ptr<SplineSurface> ssf(sweep.linearSweptSurface(*cv1, *cv2, bmid));
    double del = 1.0/(double)(nmb_sample-1);
    double p1, p2;
    int ki, kj, kr;
    int ik = degree_ + 1;
    vector<double> et(2*ik, 0.0);  // Knot vector of line curve
    for (int ki=0; ki<ik; ++ki)
      et[ik+ki] = 1.0;

    vector<double> sfpoints;
    vector<double> vecs;
    vector<double> linesegs;
    vector<double> der;
    vector<double> der2;
    vector<double> lineder;
    // Evaluate line
    vector<double> tmpline;
    for (kj=0, p2=0.0; kj<nmb_sample; ++kj, p2+=del)
      {
  	for (ki=0, p1=0.0; ki<nmb_sample; ++ki, p1+=del)
  	  {
  	    // Compute barysentric coordinates of end points of line
  	    // First cartesian
  	    Point sfpos = ssf->ParamSurface::point(p1,p2);
  	    Point cart1 = sfpos + len*dir;
  	    Point cart2 = sfpos - len*dir;
  	    tmpline.insert(tmpline.end(), cart1.begin(), cart1.end());
  	    tmpline.insert(tmpline.end(), cart2.begin(), cart2.end());

  	    Vector4D bary1 = bc_.cartToBary(Vector3D(cart1[0], cart1[1], cart1[2]));
  	    Vector4D bary2 = bc_.cartToBary(Vector3D(cart2[0], cart2[1], cart2[2]));

  	    Vector3D tp1 = bc_.baryToCart(bary1);
  	    Vector3D tp2 = bc_.baryToCart(bary2);
	    
  	    // Pick line
  	    BernsteinPoly line = implicit_.pickLine(bary1, bary2);

  	    // Compute zeroes of bernstein polynomial
  	    // First make sisl curve
  	    vector<double> ecoef(line.coefsBegin(), line.coefsEnd());
  	    SISLCurve *qc = newCurve(ik, ik, &et[0], &ecoef[0], 1, 1, 1);
  	    double zero = 0.0;

  	    // Intersect
  	    double eps = 1.0e-6;
  	    int kstat = 0;
  	    int kcrv=0, kpt=0;
  	    double *epar = 0;
  	    SISLIntcurve **intcv = 0;
  	    if (qc)
  	      s1871(qc, &zero, 1, eps, &kpt, &epar, &kcrv, &intcv, &kstat);
  	    if (qc)
  	      freeCurve(qc);

  	    // Compute cartesian points and curves associated with intersections
  	    for (kr=0; kr<kpt; ++kr)
  	      {
  		Vector4D barypt = (1.0 - epar[kr])*bary1 + epar[kr]*bary2;
  		int kb;
  		for (kb=0; kb<4; ++kb)
  		  if (barypt[kb] < -0.001 || barypt[kb] > 1.001)
  		    break;
  		if (kb < 4)
  		  continue;
		
  		Vector3D pos = bc_.baryToCart(barypt);
  		sfpoints.insert(sfpoints.end(), pos.begin(), pos.end());

  	      }
	    if (epar) free(epar);
	    if (intcv) freeIntcrvlist(intcv, kcrv);
  	  }
      }
    
    // Output
    if (sfpoints.size() > 0)
      {
  	PointCloud3D ptcloud(&sfpoints[0], sfpoints.size()/3);
  	os << "400 1 0 4 255 0 0 255" << std::endl;
  	ptcloud.write(os);
      }
}

//===========================================================================
void ImplicitApprox::visualize(vector<Point> points, Point& dir, std::ostream& os)
//===========================================================================
{
  BoundingBox bb(3);
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      bb.addUnionWith(points[ki]);
    }
  Point low = bb.low();
  Point high = bb.high();
  Point bmid = 0.5*(low + high);
  Point diag = high - low;
  double diaglen = diag.length();

  double gap = 1.0e-6;
  Point xdir(1.0, 0.0, 0.0);
  Point ydir(0.0, 1.0, 0.0);
  Point zdir(0.0, 0.0, 1.0);
  CompositeModelFactory factory(gap, gap, 10.0*gap, 0.01, 0.05);
  shared_ptr<SurfaceModel> boxmod(factory.createFromBox(low-2.0*diag, xdir, ydir, 
							5*diag[0], 5*diag[1], 5*diag[2]));
    
    // Find the coordinate direction with the largest angle with the view direction
    double a1 = xdir.angle(dir);
    double a2 = ydir.angle(dir);
    double a3 = zdir.angle(dir);
    Point dir2;
    if (a1 > std::min(a2, a3))
      dir2 = xdir;
    else if (a2 > a3)
      dir2 = ydir;
    else
      dir2 = zdir;
    Point dir3 = dir%dir2;
    dir2 = dir%dir3;
    if (dir2*(high-low) < 0.0)
      dir2 *= -1.0;
    if (dir3*(high-low) < 0.0)
      dir3 *= -1.0;
    dir2.normalize();
    dir3.normalize();
    double len = low.dist(high);
    int nmb_sample = 100;
    shared_ptr<SplineCurve> cv1(new SplineCurve(bmid-len*dir2, 0.0, bmid+len*dir2, 1.0));
    shared_ptr<SplineCurve> cv2(new SplineCurve(bmid-len*dir3, 0.0, bmid+len*dir3, 1.0));
    SweepSurfaceCreator sweep;
    shared_ptr<SplineSurface> ssf(sweep.linearSweptSurface(*cv1, *cv2, bmid));
    double del = 1.0/(double)(nmb_sample-1);
    double p1, p2;
    int ki, kj, kr;
    int ik = degree_ + 1;
    vector<double> et(2*ik, 0.0);  // Knot vector of line curve
    for (int ki=0; ki<ik; ++ki)
      et[ik+ki] = 1.0;

    vector<double> sfpoints;
    vector<double> vecs;
    vector<double> linesegs;
    vector<double> der;
    vector<double> der2;
    vector<double> lineder;
    // Evaluate line
    vector<double> tmpline;
    for (kj=0, p2=0.0; kj<nmb_sample; ++kj, p2+=del)
      {
	for (ki=0, p1=0.0; ki<nmb_sample; ++ki, p1+=del)
	  {
	    // Compute barysentric coordinates of end points of line
	    // First cartesian
	    Point sfpos = ssf->ParamSurface::point(p1,p2);
	    Point cart1 = sfpos + len*dir;
	    Point cart2 = sfpos - len*dir;
	    tmpline.insert(tmpline.end(), cart1.begin(), cart1.end());
	    tmpline.insert(tmpline.end(), cart2.begin(), cart2.end());

	    Vector4D bary1 = bc_.cartToBary(Vector3D(cart1[0], cart1[1], cart1[2]));
	    Vector4D bary2 = bc_.cartToBary(Vector3D(cart2[0], cart2[1], cart2[2]));

	    Vector3D tp1 = bc_.baryToCart(bary1);
	    Vector3D tp2 = bc_.baryToCart(bary2);
	    
	    // Pick line
	    BernsteinPoly line = implicit_.pickLine(bary1, bary2);

	    // Compute zeroes of bernstein polynomial
	    // First make sisl curve
	    vector<double> ecoef(line.coefsBegin(), line.coefsEnd());
	    SISLCurve *qc = newCurve(ik, ik, &et[0], &ecoef[0], 1, 1, 1);
	    double zero = 0.0;

	    // Intersect
	    double eps = 1.0e-6;
	    int kstat = 0;
	    int kcrv=0, kpt=0;
	    double *epar = 0;
	    SISLIntcurve **intcv = 0;
	    if (qc)
	      s1871(qc, &zero, 1, eps, &kpt, &epar, &kcrv, &intcv, &kstat);
	    if (qc)
	      freeCurve(qc);

	    // Compute cartesian points and curves associated with intersections
	    for (kr=0; kr<kpt; ++kr)
	      {
		Vector4D barypt = (1.0 - epar[kr])*bary1 + epar[kr]*bary2;
		int kb;
		for (kb=0; kb<4; ++kb)
		  if (barypt[kb] < -0.001 || barypt[kb] > 1.001)
		    break;
		if (kb < 4)
		  continue;
		
		Vector3D pos = bc_.baryToCart(barypt);
		sfpoints.insert(sfpoints.end(), pos.begin(), pos.end());

	      }
	    if (epar) free(epar);
	    if (intcv) freeIntcrvlist(intcv, kcrv);
	  }
      }
    
    // Output
    if (sfpoints.size() > 0)
      {
	PointCloud3D ptcloud(&sfpoints[0], sfpoints.size()/3);
	os << "400 1 0 4 255 0 0 255" << std::endl;
	ptcloud.write(os);
      }
}


double fc(int deg, int power[], double coef[], double x, double y, double z)
{
  int nn = (deg+1)*(deg+2)*(deg+3)/6;
  double res = 0.0;
  for (int ki=0; ki<nn; ++ki)
    {
      double tmp1 = 1;
      for (int kj=0; kj<power[4*ki]; ++kj)
	tmp1 *= x;
      double tmp2 = 1;
      for (int kj=0; kj<power[4*ki+1]; ++kj)
	tmp2 *= y;
      double tmp3 = 1;
      for (int kj=0; kj<power[4*ki+2]; ++kj)
	tmp3 *= z;
      res += (coef[ki]*tmp1*tmp2*tmp3);
    }
  return res;
}

double f(int deg, int power[], int ki, double x, double y, double z)
{
  double tmp1 = 1;
  for (int kj=0; kj<power[4*ki]; ++kj)
    tmp1 *= x;
  double tmp2 = 1;
  for (int kj=0; kj<power[4*ki+1]; ++kj)
    tmp2 *= y;
  double tmp3 = 1;
  for (int kj=0; kj<power[4*ki+2]; ++kj)
    tmp3 *= z;
  double res = tmp1*tmp2*tmp3;
  return res;
}

double fx(int deg, int power[], int ki, double x, double y, double z)
{
  double res = 0.0;
  if (power[4*ki] > 0)
    {
      double tmp = (double)power[4*ki];
      for (int kj=0; kj<power[4*ki]-1; ++kj)
	tmp *= x;
      for (int kj=0; kj<power[4*ki+1]; ++kj)
	tmp *= y;
      for (int kj=0; kj<power[4*ki+2]; ++kj)
	tmp *= z;
      res = tmp;
    }

  return res;
}

double fy(int deg, int power[], int ki, double x, double y, double z)
{
  double res = 0.0;
  if (power[4*ki+1] > 0)
    {
      double tmp = (double)power[4*ki+1];
      for (int kj=0; kj<power[4*ki+1]-1; ++kj)
	tmp *= y;
      for (int kj=0; kj<power[4*ki]; ++kj)
	tmp *= x;
      for (int kj=0; kj<power[4*ki+2]; ++kj)
	tmp *= z;
      res = tmp;
    }
  return res;
}

double fz(int deg, int power[], int ki, double x, double y, double z)
{
  double res = 0.0;
  if (power[4*ki+2] > 0)
    {
      double tmp = (double)power[4*ki+2];
      for (int kj=0; kj<power[4*ki+2]-1; ++kj)
	tmp *= z;
      for (int kj=0; kj<power[4*ki]; ++kj)
	tmp *= x;
      for (int kj=0; kj<power[4*ki+1]; ++kj)
	tmp *= y;
      res = tmp;
    }
  return res;
}


//===========================================================================
void ImplicitApprox::polynomialSurf(vector<Point>& pos_and_der, int degree,
				    vector<double>& coefs)
//===========================================================================
{
  // Assemble matrix
  int nmbvar = (degree+1)*(degree+2)*(degree+3)/6;
  size_t nmb_pts = pos_and_der.size()/3;
  vector<vector<double> > M(nmb_pts);
  vector<vector<double> > N1(nmb_pts);
  vector<vector<double> > N2(nmb_pts);
  for (size_t kr=0; kr<nmb_pts; ++kr)
    {
      M[kr].resize(nmbvar, 0.0);
      N1[kr].resize(nmbvar, 0.0);
      N2[kr].resize(nmbvar, 0.0);
    }

  vector<int> power(4*nmbvar);
  for (int ki=0, kr=0; ki<=degree; ++ki)
    for (int kj=0; kj<=degree-ki; ++kj)
      for (int kh=0; kh<=degree-ki-kj; ++kh, kr+=4)
      {
	power[kr] = ki;
	power[kr+1] = kj;
	power[kr+2] = kh;
	power[kr+3] = degree-ki-kj-kh;
      }
  
  for (size_t kr=0, kh=0; kr<pos_and_der.size(); kr+=3, ++kh)
    {
      for (int ki=0; ki<nmbvar; ++ki)
	{
	  double m1 = f(degree, &power[0], ki, pos_and_der[kr][0],
			pos_and_der[kr][1], pos_and_der[kr][2]);
	  Point tmp1(3);
	  tmp1[0] = fx(degree, &power[0], ki, pos_and_der[kr][0],
		       pos_and_der[kr][1], pos_and_der[kr][2]);
	  tmp1[1] = fy(degree, &power[0], ki, pos_and_der[kr][0],
		       pos_and_der[kr][1], pos_and_der[kr][2]);
	  tmp1[2] = fz(degree, &power[0], ki, pos_and_der[kr][0],
		       pos_and_der[kr][1], pos_and_der[kr][2]);
	  double n1 = tmp1*pos_and_der[kr+1];
	  double n2 = tmp1*pos_and_der[kr+2];
	  M[kh][ki] += m1;
	  N1[kh][ki] += n1;
	  N2[kh][ki] += n2;
	}
    }


  double lambda = 0.2;
  NEWMAT::Matrix mat;
  mat.ReSize(nmb_pts,nmbvar);
  for (int ki=0; ki<nmb_pts; ++ki)
    for (int kj=0; kj<nmbvar; ++kj)
      mat.element(ki,kj) = M[ki][kj] + lambda*(N1[ki][kj] + N2[ki][kj]);

  static NEWMAT::DiagonalMatrix diag;
  static NEWMAT::Matrix V;
  try {
    NEWMAT::SVD(mat, diag, mat, V);
  } catch(...) {
    std::cout << "Exception in SVD" << std::endl;
    return;
  }

  int write_info = 1;
  if (write_info)
    {
      std::cout << "Singular values:" << std::endl;
      for (int ki = 0; ki < nmbvar; ++ki)
	std::cout << ki << "\t" << diag.element(ki, ki) << std::endl;
      
      // Write out info about singular values
      double s_min = diag.element(nmbvar-1,nmbvar-1);
      double s_max = diag.element(0, 0);
      std::cout << "Implicitization:" << std::endl
		<< "s_min = " << s_min << std::endl
		<< "s_max = " << s_max << std::endl
		<< "Ratio of s_min/s_max = " << s_min/s_max << std::endl;
      std::cout << "Ratio s_min/s_next_min = " << diag.element(nmbvar-1,nmbvar-1)/diag.element(nmbvar-2,nmbvar-2) << std::endl;
    }

   coefs.resize(nmbvar);
  for (int ki=0; ki<nmbvar; ++ki)
    coefs[ki] = V.element(ki, nmbvar-1);

}

//===========================================================================
void ImplicitApprox::polynomialSurfAccuracy(vector<Point>& pos_and_der, 
					    int degree, vector<double>& coefs,
					    double& maxfield, double& avfield,
					    double& maxdist, double& avdist,
					    int& ndiv, double& maxang,
					    double& avang)
//===========================================================================
{
  int nmbvar = (degree+1)*(degree+2)*(degree+3)/6;
  size_t nmb_pts = pos_and_der.size()/3;
  double eps = 1.0e-10;
  double min_grad = 1.0e-9;

  vector<int> power(4*nmbvar);
  for (int ki=0, kr=0; ki<=degree; ++ki)
    for (int kj=0; kj<=degree-ki; ++kj)
      for (int kh=0; kh<=degree-ki-kj; ++kh, kr+=4)
      {
	power[kr] = ki;
	power[kr+1] = kj;
	power[kr+2] = kh;
	power[kr+3] = degree-ki-kj-kh;
      }
  
  // Test accuracy
  maxfield = 0.0;
  avfield = 0.0;
  double avgradlen = 0.0;
  double mingradlen = 1.0e8;
  double maxgradlen = 0.0;
  double minang = 1.0e8;
  maxang = 0.0;
  avang = 0.0;
  double mindist = 1.0e8;
  maxdist = 0.0;
  avdist = 0.0;
  vector<Point> out;
  ndiv = 0;
  for (size_t kr=0; kr<nmb_pts; kr+=3)
    {
      double field = 0.0;
      double dx=0.0, dy=0.0, dz=0.0;
      for (int ki=0; ki<nmbvar; ++ki)
	{
	  field += coefs[ki]*f(degree, &power[0], ki, pos_and_der[kr][0],
			       pos_and_der[kr][1], pos_and_der[kr][2]);
	  dx += coefs[ki]*fx(degree, &power[0], ki, pos_and_der[kr][0],
			     pos_and_der[kr][1], pos_and_der[kr][2]);
	  dy += coefs[ki]*fy(degree, &power[0], ki, pos_and_der[kr][0],
			     pos_and_der[kr][1], pos_and_der[kr][2]);
	  dz += coefs[ki]*fz(degree, &power[0], ki, pos_and_der[kr][0],
			     pos_and_der[kr][1], pos_and_der[kr][2]);
	}
      Point grad(dx, dy, dz);
      Point normc = pos_and_der[kr+1].cross(pos_and_der[kr+2]);
      double ang = normc.angle(grad);
      minang = std::min(minang, ang);
      maxang = std::max(maxang, ang);
      avang += ang;
      double gradlen = grad.length();
      double edist = (gradlen < 1.0e-17) ? 0.0 : fabs(field)/gradlen;
      maxfield = std::max(maxfield, fabs(field));
      avfield += fabs(field);
      mingradlen = std::min(mingradlen, gradlen);
      maxgradlen = std::max(maxgradlen, gradlen);
      avgradlen += gradlen;
      if (gradlen > 1.0e-10)
	grad.normalize();

      double delta = 1.0e-9;
      Point norm = grad; //tan1[kr].cross(tan2[kr]);
      norm.normalize();
      double t0 = 0.0;
      Point pos0 = pos_and_der[kr];
      double dt = dx*norm[0] + dy*norm[1] + dz*norm[2];
      double tdel = -field/dt;
      double field0;
      for (int ka=0; ka<10; ++ka)
	{
	  if (fabs(tdel) < delta)
	    break;
	  t0 += tdel;
	  pos0 = pos_and_der[kr] + t0*norm;
	  double dx0 = 0.0, dy0 = 0.0, dz0 = 0.0;
	  field0 = 0.0;
	  for (int ki=0; ki<nmbvar; ++ki)
	    {
	      field0 += coefs[ki]*f(degree, &power[0], ki, pos0[0],pos0[1],pos0[2]);
	      dx0 += coefs[ki]*fx(degree, &power[0], ki, pos0[0],pos0[1],pos0[2]);
	      dy0 += coefs[ki]*fy(degree, &power[0], ki, pos0[0],pos0[1],pos0[2]);
	      dz0 += coefs[ki]*fz(degree, &power[0], ki, pos0[0],pos0[1],pos0[2]);
	    }
	  dt = dx0*norm[0] + dy0*norm[1] + dz0*norm[2];
	  tdel = -field0/dt;
	  int stop_break = 1;
	}

      if (fabs(field0) > fabs(field) || fabs(field0) > eps || fabs(t0) > 1.0 ||
	  gradlen < min_grad)
	{
	  out.push_back(pos_and_der[kr]);
	  ndiv++;
	}
      else
	{
	  maxdist = std::max(maxdist, fabs(t0));
	  mindist = std::min(mindist, fabs(t0));
	  avdist += fabs(t0);
	}
   }
  avfield /= (double)nmb_pts;
  avang /= (double)nmb_pts;
  avgradlen /= (double)nmb_pts;
  avdist /= (double)(nmb_pts-ndiv);
  
  std::cout << "Maximum field: " << maxfield << std::endl;
  std::cout << "Average field: " << avfield << std::endl;
  std::cout << "Minimum distance: " << mindist << std::endl;
  std::cout << "Maximum  distance: " << maxdist << std::endl;
  std::cout << "Average  distance: " << avdist << std::endl;
  std::cout << "Num divergent: " << ndiv << std::endl;
  std::cout << "Minimum angle difference: " << minang << std::endl;
  std::cout << "Maximum angle difference: " << maxang << std::endl;
  std::cout << "Avarage angle difference: " << avang << std::endl;
  std::cout << "Minimum gradient length: " << mingradlen << std::endl;
  std::cout << "Maximum gradient length: " << maxgradlen << std::endl;
  std::cout << "Average gradient length: " << avgradlen << std::endl;

  std::ofstream ofo("out.g2");
  ofo << "400 1 0 4 50 50 155 255" << std::endl;
  ofo << out.size() << std::endl;
  for (int ka=0; ka<(int)out.size(); ++ka)
    ofo << out[ka] << std::endl;

  int stop_break = 1;
  
}
