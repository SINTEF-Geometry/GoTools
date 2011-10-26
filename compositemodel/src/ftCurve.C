//===========================================================================
//                                                                           
// File: ftCurve.C                                                           
//                                                                           
// Created: Mon Apr  3 15:21:58 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: ftCurve.C,v 1.3 2009-06-12 08:55:15 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include <string>
#include "GoTools/compositemodel/ftCurve.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/compositemodel/ftPlane.h"
#include "GoTools/compositemodel/ftFaceBase.h"
#include "GoTools/geometry/orientCurves.h"
#include "GoTools/tesselator/CurveTesselator.h"

using std::vector;
using std::endl;
using std::string;
using std::shared_ptr;

namespace Go
{


//===========================================================================
ftCurveSegment::ftCurveSegment(ftCurveType type,
			       tpJointType joint,
			       ftFaceBase* f0,
			       ftFaceBase* f1,
			       shared_ptr<ParamCurve> paramcv0,
			       shared_ptr<ParamCurve> paramcv1,
			       shared_ptr<ParamCurve> spacecv,
			       double eps_geo)
//===========================================================================
    : segment_type_(type), joint_(joint), space_curve_(spacecv)
{

//     ALWAYS_ERROR_IF(space_curve_.get() == 0,
// 		"You must specify a space curve for the segment.",
// 		InputError());

    underlying_face_[0] = f0;
    underlying_face_[1] = f1;
    parameter_curve_[0] = paramcv0;
    parameter_curve_[1] = paramcv1;

    // There must always be a spacecurve
    if (space_curve_.get() == 0) {
	redefineSpaceCurve(eps_geo);
    }

}


//===========================================================================
double ftCurveSegment::startOfSegment() const
//===========================================================================
{
    return space_curve_->startparam();
}

//===========================================================================
double ftCurveSegment::endOfSegment() const
//===========================================================================
{
    return space_curve_->endparam();
}

//===========================================================================
void ftCurveSegment::point(double t, Point& pt) const
//===========================================================================
{
    space_curve_->point(pt, t);
}

//===========================================================================
void ftCurveSegment::tangent(double t, Point& tan) const
//===========================================================================
{
    vector<Point> pts(2, Point(0.0, 0.0, 0.0));
    space_curve_->point(pts, t, 1);
    tan = pts[1];
}

//===========================================================================
void ftCurveSegment::paramcurvePoint(int number, double t,
				     Point& pt) const
//===========================================================================
{
  if (parameter_curve_[number].get() != 0)
    parameter_curve_[number]->point(pt, t);
  else
    {
      // Find the point on the curve
      Point pts(space_curve_->dimension());
      space_curve_->point(pts, t);
      Point clo_pt = pts;
      double clo_u, clo_v, clo_dist;
      // Find the closest point on the surface
      shared_ptr<ParamSurface> srf = underlying_face_[number]->surface();
      // cout << srf.get() << endl;
      srf->closestPoint(pts, clo_u, clo_v, clo_pt, clo_dist, 1e-10);
      pt.setValue(clo_u, clo_v);
    }
}

//===========================================================================
void ftCurveSegment::paramcurveTangent(int number, double t,
				       Point& tan) const
//===========================================================================
{
  if (parameter_curve_[number].get() != 0)
    {
      vector<Point> pts(2, Point(0.0, 0.0));
      parameter_curve_[number]->point(pts, t, 1);
      tan = pts[1];
    }
  else
    {
      // Find the point on the curve
      vector<Point> pts(2, Point(0.0, 0.0, 0.0));
      space_curve_->point(pts, t, 1);
      Point clo_pt = pts[0];
      double clo_u, clo_v, clo_dist;
      // Find the closest point on the surface
      shared_ptr<ParamSurface> srf = underlying_face_[number]->surface();
      // cout << srf.get() << endl;
      srf->closestPoint(pts[0], clo_u, clo_v, clo_pt, clo_dist, 1e-10);

      // Evaluate surface
      vector<Point> pts2(3, Point(0.0, 0.0, 0.0));
      srf->point(pts2, clo_u, clo_v, 1);

      double fufu = pts2[1]*pts2[1];
      double fufv = pts2[1]*pts2[2];
      double fvfv = pts2[2]*pts2[2];
      double sdfu = pts[1]*pts2[1];
      double sdfv = pts[1]*pts2[2];

      // Solve min(s-Fu*x-Fv*y)^2 to find x and y.
      // s(t) = F(u(t),v(t)), s(t) is the curve in geometry space,
      // F(u,v) is the surface. (u(t),v(t)) is the curve in the
      // parameter space of F, x = u'(t) and y = v'(t).
      double det = fufu*fvfv - fufv*fufv;
      if (fabs(det) < 10.0e-11)
	tan = Point(0.0, 0.0);
      else
	{
	  double x = (fvfv*sdfu - fufv*sdfv)/det;
	  double y = (fufu*sdfv - fufv*sdfu)/det;
	  tan.setValue(x,y);
	}
    }
}

//===========================================================================
void ftCurveSegment::reverse()
//===========================================================================
{
    // @ After reversal, the joint types will be wrong.
    // This must be handled by the calling function, for
    // instance ftCurve::reverse().
    //    joint_ = JOINT_NONE;
    if (parameter_curve_[0].get())
	parameter_curve_[0]->reverseParameterDirection();
    if (parameter_curve_[1].get())
	parameter_curve_[1]->reverseParameterDirection();
    if (space_curve_.get())
	space_curve_->reverseParameterDirection();
  
}

//===========================================================================
Point ftCurveSegment::startPoint() const
//===========================================================================
{
    Point pt;
    point(startOfSegment(), pt);
    return pt;
}

//===========================================================================
Point ftCurveSegment::endPoint() const
//===========================================================================
{
    Point pt;
    point(endOfSegment(), pt);
    return pt;
}

// @@ Value eps currently not being used. Using hardcoded 1e-10.
//===========================================================================
void ftCurveSegment::normal(double t, int side, Point& normal,
			    double eps) const
//===========================================================================
{
  if (parameter_curve_[side].get() != 0)
    {
      Point par_pt;
      parameter_curve_[side]->point(par_pt, t);
      normal = underlying_face_[side]->normal(par_pt[0], par_pt[1]);
    }
  else
    {
      // Find the point on the curve
      Point pts = space_curve_->point(t);
      Point clo_pt;
      double clo_u, clo_v, clo_dist;
      // Find the closest point on the surface
      shared_ptr<ParamSurface> srf = underlying_face_[side]->surface();
      // cout << srf.get() << endl;
      srf->closestPoint(pts, clo_u, clo_v, clo_pt, clo_dist, 1e-10);

      // Evaluate surface
      normal = underlying_face_[side]->normal(clo_u, clo_v);
    }
}


//===========================================================================
double ftCurveSegment::arcLength(double t1, double t2) const
//===========================================================================
{
    MESSAGE("Arc length function is a slow hack");
    int num_samp = 10;
    double l = 0;
    Point old_pt = space_curve_->point(t1);
    Point new_pt; 
    double inc = (t2 - t1)/(double)(num_samp - 1);
    for (int i = 1; i < num_samp; ++i) {
	new_pt = space_curve_->point(t1 + i*inc);
	l += old_pt.dist(new_pt);
	old_pt = new_pt;
    }
    return l;
}

//===========================================================================
shared_ptr<LineStrip> ftCurveSegment::tesselate(int resolution) const
//===========================================================================
{
  CurveTesselator tesselator(*space_curve_.get());
  tesselator.changeRes(resolution);
  shared_ptr<LineStrip> strip =  tesselator.getMesh();
  return strip;
}

//===========================================================================
double ftCurve::arcLength(double t1, int seg1, double t2, int seg2) const
//===========================================================================
{
    double length = 0.0;
    double a, b;
    for (int i = seg1; i <= seg2; ++i) {
	if (i == seg1)
	    a = t1;
	else
	    a = startOfSegment(i);
	if (i == seg2)
	    b = t2;
	else
	    b = endOfSegment(i);
	length += segments_[i].arcLength(a, b);
    }
    return length;
}

//===========================================================================
void ftCurve::reparametrize(double eps_go)
//===========================================================================
{
    // We will reparametrize every segment from 0 to its length
    for (size_t ki = 0; ki < segments_.size(); ++ki)
	segments_[ki].reparametrize(eps_go);
}

//===========================================================================
void ftCurve::joinSegments(double gap_tol, double neighbour_tol,
			   double kink_tol, double bend_tol)
//===========================================================================
{
    int n = (int)segments_.size();
    for (int i = 0; i < n; ++i) {
	int j = (i+1) % n; // Next segment to test against
	Point p1, p2, t1, t2;
	segments_[i].point(segments_[i].endOfSegment(), p1);
	segments_[j].point(segments_[j].startOfSegment(), p2);
	double dist = p1.dist(p2);
	if (dist > neighbour_tol) {
	    segments_[i].setJointAfter(JOINT_DISC);
	} else if (dist > gap_tol) {
	    segments_[i].setJointAfter(JOINT_GAP);
	} else {
	    // The points coincide
	    segments_[i].tangent(segments_[i].endOfSegment(), t1);
	    segments_[j].tangent(segments_[j].startOfSegment(), t2);
	    double ang_dist = t1.angle(t2);
	    if (ang_dist > bend_tol) {
		segments_[i].setJointAfter(JOINT_G0);
	    } else if (ang_dist > kink_tol) {
		segments_[i].setJointAfter(JOINT_KINK);
	    } else {
		// The tangents coincide, too!
		segments_[i].setJointAfter(JOINT_G1);
	    }
	}
    }
    if (n > 0)
	if (segments_[n-1].jointAfter() == JOINT_DISC)
	    segments_[n-1].setJointAfter(JOINT_NONE);
}

//===========================================================================
void ftCurve::orientSegments(double neighbour_tol, bool assume_manifold)
//===========================================================================
{
    // @ Note that all segments must have a space curve defined
    int num_seg = numSegments();
    vector<ParamCurve*> curves(num_seg);
    for (int i = 0; i < num_seg; ++i) {
	curves[i] = segments_[i].spaceCurve().get();
	ASSERT(curves[i] != 0);
    }
    vector<int> perm(num_seg);
    vector<bool> flip(num_seg);
    orientCurves(curves, perm, flip, neighbour_tol, assume_manifold);
    // Making the new segment vector
    vector<ftCurveSegment> new_segments;
    new_segments.reserve(num_seg);
    for (int i = 0; i < num_seg; ++i) {
	new_segments.push_back(segments_[perm[i]]);
	if (flip[i]) {
	    new_segments[i].reverse();
	}
    }
    segments_.swap(new_segments);
}

//===========================================================================
void ftCurve::tesselate(vector<shared_ptr<LineStrip> >& meshes) const
//===========================================================================
{
  int res = 100;
  tesselate(res, meshes);
}

//===========================================================================
  void ftCurve::tesselate(int resolution,
			  vector<shared_ptr<LineStrip> >& meshes) const
//===========================================================================
{
  meshes.clear();
  for (size_t ki=0; ki<segments_.size(); ++ki)
    {
      shared_ptr<LineStrip> strip = segments_[ki].tesselate(resolution);
      meshes.push_back(strip);
    }
  
}

//===========================================================================
  void ftCurve::tesselate(double density,
			  vector<shared_ptr<LineStrip> >& meshes) const
//===========================================================================
{
  int min_nmb = 5;
  int max_nmb = 10000;
  meshes.clear();
  for (size_t ki=0; ki<segments_.size(); ++ki)
    {
      double len = segments_[ki].spaceCurve()->estimatedCurveLength();
      int resolution = (int)(len/density);
      resolution = std::max(min_nmb, std::min(resolution, max_nmb));
      shared_ptr<LineStrip> strip = segments_[ki].tesselate(resolution);
      meshes.push_back(strip);
    }
  
}

//===========================================================================
void ftCurve::write(std::ostream& os) const
//===========================================================================
{
    string ctypes[] = {"CURVE_INTERSECTION", "CURVE_EDGE",
			"CURVE_KINK", "CURVE_CORNER", "CURVE_GAP" };
    string jtypes[] = {"JOINT_G1", "JOINT_KINK", "JOINT_G0",
			"JOINT_GAP", "JOINT_DISC", "JOINT_NONE"};
    int n = (int)segments_.size();
    os << "Curve type is " << (type_==CURVE_NOTYPE ? 
			       string("CURVE_NOTYPE") : ctypes[type_])
	 << " with " << n << " segments." << endl;
    Point p1, p2;
    for (int i = 0; i < n; ++i) {
	os << "Segment " << i << endl;
	ftCurveType ty = segments_[i].segmentType();
	os << "Segment type is " << (ty==CURVE_NOTYPE ? 
				     string("CURVE_NOTYPE") : ctypes[ty])
	   << endl;
	segments_[i].point(segments_[i].startOfSegment(), p1);
	segments_[i].point(segments_[i].endOfSegment(), p2);
	os << p1 << p2;
	os << "Joint to next segment is "
	     << jtypes[segments_[i].jointAfter()] << endl;
    }
}


//===========================================================================
void ftCurve::writeSpaceCurve(std::ostream& os) const
//===========================================================================
{
    int n = (int)segments_.size();
    for (int i = 0; i < n; ++i) {
	segments_[i].spaceCurve()->writeStandardHeader(os);
	segments_[i].spaceCurve()->write(os);
	os << endl << endl;
    }
}



} // namespace Go




