//===========================================================================
//                                                                           
// File: ftCurve.h                                                           
//                                                                           
// Created: Thu Mar 23 11:34:00 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: ftCurve.h,v 1.6 2009-06-12 08:55:14 vsk Exp $
//                                                                           
// Description: ftCurve and helper ftCurveSegment classes
//                                                                           
//===========================================================================

#ifndef _FTCURVE_H
#define _FTCURVE_H


#include <algorithm>
#include <vector>
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/topology/tpJointType.h"
#include "GoTools/compositemodel/ftFaceBase.h"

namespace Go
{

  class LineStrip;
class ftFaceBase;

/** Indicates type of curve described by a ftCurve. */
enum ftCurveType
{
    CURVE_NOTYPE = -1,   ///< Untyped curve
    CURVE_INTERSECTION,  ///< Curve mostly in interior of surfaces
    CURVE_EDGE,          ///< Curve lying on the edge of surfaces
    CURVE_KINK,          ///< Between surfaces joined by kink
    CURVE_CORNER,        ///< Between surfaces meeting at a corner line
    CURVE_GAP,           ///< Between surfaces joined by gap
    CURVE_SINGULAR       ///< Curve mostly in interior of one surface
};




//===========================================================================
/** ftCurveSegment -  One segment describing a curve lying on one or two surfaces.
 * The curve can be of type intersection, edge curve, kink curve, corner curve or
 * gap curve. The curve is an output of a computation corresponding to the curve
 * type.
 */
//===========================================================================
class GO_API ftCurveSegment
{

public:
    /** Constructor. The segment owns the nurbscurves, but not the faces
     */
    ftCurveSegment(ftCurveType type,
		   tpJointType joint,
		   ftFaceBase* f0,
		   ftFaceBase* f1,
		   shared_ptr<ParamCurve> paramcv0,
		   shared_ptr<ParamCurve> paramcv1,
		   shared_ptr<ParamCurve> spacecv,
		   double eps_geo = 1.0e-6);

    /// Fetch info about the type of this curve segment
    /// -1 = no specified type, 0 = intersection curve, 1 = edge curve,
    /// 2 = kink curve, 3 = corner curve, 4 = gap curve, 5 = singularity curve
    inline ftCurveType segmentType() const;
    
    /// Set joint type between this curve segment and the next
    /// 0 = g1 continuity, 1 = not quite g1, 2 = g0 continuity, 3 = gap,
    /// 4 = discontinuous, 5 = last segment
    inline void setJointAfter(tpJointType jt);

    /// Fetch joint type between this curve segment and the next
    inline tpJointType jointAfter() const;

    /// Fetch face number number associated with this curve segment. 
    inline ftFaceBase* face(int number) const;

    /// Return start parameter of segment
    double startOfSegment() const;

    /// Return end parameter of segment
    double endOfSegment() const;

    /** Evaluate a point on the curve. */
    void point(double t, Point& pt) const;
    /** Get the curve tangent in a given parameter. */
    void tangent(double t, Point& tan) const;
    /// Comput point on curve in the parameter domain of face number in the
    /// parameter t. If the parameter curve does not exist, the closest point
    /// in the specified face to the space curve point corresponding to t is
    /// returned.
    void paramcurvePoint(int number, double t, Point& pt) const;
    /// Tangent to the curve in the parameter domain of face number
    void paramcurveTangent(int number, double t, Point& pt) const;

    /// Reverse the direction of the curve segment. NB! The joint types will not be
    /// correct
    void reverse();

    /// Internal use
    std::vector<ftCurveSegment> chopOff(const BoundingBox& box, bool& eraseme);

    /// Return the space curve corresponding to this segment
    shared_ptr<ParamCurve> spaceCurve() const;

    /// Return the curve in the parameter domain of face number
    shared_ptr<ParamCurve> parameterCurve(int number) const;

    /// Start point of segment
    Point startPoint() const;

    /// End point of segment
    Point endPoint() const;

    /** The normal of the surface corresponding to this curve. */
    void normal(double t, int side, Point& normal, double eps) const;

    /// The arc length of this curve segment limited by the interval [t1, t2]
    double arcLength(double t1, double t2) const;
    
    /** The eps_go argument is used as a tolerance if ever this function has to 
     * regenerate the spacecurve from its underlying surface and parametric curve */
    void reparametrize(double eps_go = 1.0e-6);

    /// Tesselate curve segment to return a linear approximation
    shared_ptr<LineStrip> tesselate(int resolution) const;

    /// Debug purpose
    void deleteSpaceCurveRepresentation() { // @@ Debug purpose - remove later!
	if ((parameter_curve_[0].get() != 0 && underlying_face_[0] != 0)||
	    (parameter_curve_[1].get() != 0 && underlying_face_[1] != 0)) {
	    std::cout << "Removing curve " << std::endl;
	    space_curve_.reset();
	} else {
	    std::cerr << "Parameter curve not present - could not ";
	    std::cerr << " delete space curve. " << std::endl;
	}
    }

protected:
    ftCurveType segment_type_;
    tpJointType joint_; // Use 1 if edge/intersect, 2 if kink/gap
    /// Surface(s)  on which curve is lying
    ftFaceBase* underlying_face_[2];    
    /// Curves in surface domain
    shared_ptr<ParamCurve> parameter_curve_[2]; 
    /// Curve in space
    shared_ptr<ParamCurve> space_curve_;     

    // make sure the space_curve_ is a SplineCurve
    void redefineSpaceCurve(double eps_go);
};


//===========================================================================
inline ftCurveType ftCurveSegment::segmentType() const
//===========================================================================
{
    return segment_type_;
}

//===========================================================================
inline void ftCurveSegment::setJointAfter(tpJointType jt)
//===========================================================================
{
    joint_ = jt;
}

//===========================================================================
inline tpJointType ftCurveSegment::jointAfter() const
//===========================================================================
{
    return joint_;
}

//===========================================================================
inline ftFaceBase* ftCurveSegment::face(int number) const
//===========================================================================
{
    return underlying_face_[number];
}

//===========================================================================
inline shared_ptr<ParamCurve> ftCurveSegment::spaceCurve() const
//===========================================================================
{
    return space_curve_;
}

//===========================================================================
inline shared_ptr<ParamCurve> ftCurveSegment::parameterCurve(int number) const
//===========================================================================
{
    return parameter_curve_[number];
}

//===========================================================================
/** ftCurve - contains a number of (possibly disjoint) curve segments.
 * The ftCurve class has interfaces for most operations (such as evaluation)
 * that one might want to perform on them
 *
 * \author Atgeirr F Rasmussen <atgeirr@sintef.no>
 */
//===========================================================================
class GO_API ftCurve
{
public:
    /// Constructor
  ftCurve() {}

  /// Constructor of curve given curve type
  /// -1 = no specified type, 0 = intersection curve, 1 = edge curve,
  /// 2 = kink curve, 3 = corner curve, 4 = gap curve, 5 = singularity curve
    ftCurve(ftCurveType t) : type_(t){}

  /// Append curve
  ftCurve& operator += (const ftCurve& cv);

  /// Append the ftCurve cv after this curve
  void appendCurve(const ftCurve& cv);

  /// Place the ftCurve cv before this curve
  void prependCurve(const ftCurve& cv);

  /// A new segment after this curve
  void appendSegment(const ftCurveSegment& seg);

  /// A new segment prior to this curve
  void prependSegment(const ftCurveSegment& seg);

    /** Flag the connections between segments according to their continuity.
     * This function assumes that the segments have already been sorted and oriented.
     * Sorting and orienting of segments can be done with orientSegments() */
    void joinSegments(double gap_tol, double neighbour_tol,
		      double kink_tol, double bend_tol);

    /** Assures that the segments are consistently oriented and sorted such that 
     * the end point of segment 'k' is coincident (within the specified tolerance) with
     * the start point of segment 'k+1' if those segments are connected.  If we cannot
     * a priori be sure that our curve is describing an 1-manifold (no 'branches' on the
     * curve), the last argument should be set to 'false'.  This may effect performance,
     * though not severely unless the number of segments is very high. */
    void orientSegments(double neighbour_tol, bool assume_manifold = true);
    
    /// Reverse the curve and maintain correct joint types between segments
    void reverse();

    /// Remove the parts of this curve lying outside the specified bounding box
    void chopOff(const BoundingBox& box);

    // Curve info
    /// Number of segments in ftCurve
    int numSegments() const;

    /// Start parameter of segment number segment
    double startOfSegment(int segment) const;
    /// End parameter of segment number segment
    double endOfSegment(int segment) const;
    /// Fetch segment number segment
    const ftCurveSegment& segment(int segment) const;
    /// Type of curve
    /// -1 = no specified type, 0 = intersection curve, 1 = edge curve,
    /// 2 = kink curve, 3 = corner curve, 4 = gap curve, 5 = singularity curve    
    ftCurveType curveType() const;
    /// Type of joint between segment number segment and the next segment
    /// 0 = g1 continuity, 1 = not quite g1, 2 = g0 continuity, 3 = gap,
    /// 4 = discontinuous, 5 = last segment
    tpJointType jointAfter(int segment) const;
    /// Least continous joint between segments
    tpJointType worstJointType() const;
    /// Number of disjoint subcurves
    int numDisjointSubcurves() const;

    // Evaluation
    /// Compute position given segment number and parameter within segment
    void point(double t, int segment, Point& pt) const;
    /// Compute tangent given segment number and parameter within segment
    void tangent(double t, int segment, Point& tan) const;
    /// Compute normal given segment number, parameter within segment and the
    /// number of the associated face
    void normal(double t, int segment, int side, Point& normal,
		double eps) const;

    /// Curve length
    double arcLength(double t1, int seg1, double t2, int seg2) const;

    /// Reparametrization
    /** The eps_go argument is used as a tolerance if ever this function has to 
     * regenerate the spacecurve from its underlying surface and parametric curve */
    void reparametrize(double eps_go = 1.0e-6);

    /// Tesselation with default number of line segments
    void tesselate(std::vector<shared_ptr<LineStrip> >& meshes) const;
    /// Tesselation, number of line segments for each ftCurveSegment given
    void tesselate(int resolution, 
		   std::vector<shared_ptr<LineStrip> >& meshes) const;
    /// Tesselation, number of line segments for each ftCurveSegment depends on
    /// the length of the segment and the given density
    void tesselate(double density, 
		   std::vector<shared_ptr<LineStrip> >& meshes) const;

    // For debugging
    /// Debug
    void write(std::ostream& os) const;
    /// Debug
    void writeSpaceCurve(std::ostream& os) const;
    /// Debug
    void deleteSpaceCurveRepresentation() { // @@ Debug purpose - remove later!
	for (size_t i = 0; i < segments_.size(); ++i) {
	    segments_[i].deleteSpaceCurveRepresentation();
	}
    }

protected:
    ftCurveType type_;
    std::vector<ftCurveSegment> segments_;
};


//===========================================================================
inline ftCurve& ftCurve::operator += (const ftCurve& cv)
//===========================================================================
{
    appendCurve(cv);
    return *this;
}

//===========================================================================
inline void ftCurve::appendCurve(const ftCurve& cv)
//===========================================================================
{
    segments_.insert(segments_.end(),
		     cv.segments_.begin(),
		     cv.segments_.end());
}

//===========================================================================
inline void ftCurve::prependCurve(const ftCurve& cv)
//===========================================================================
{
    segments_.insert(segments_.begin(),
		     cv.segments_.begin(),
		     cv.segments_.end());
}

//===========================================================================
inline void ftCurve::appendSegment(const ftCurveSegment& seg)
//===========================================================================
{
    segments_.push_back(seg);
}

//===========================================================================
inline void ftCurve::prependSegment(const ftCurveSegment& seg)
//===========================================================================
{
    segments_.insert(segments_.begin(), seg);
}

//===========================================================================
inline void ftCurve::reverse()
//===========================================================================

{
    std::reverse(segments_.begin(), segments_.end());
    size_t i;
    for (i = 0; i< segments_.size(); ++i)
        segments_[i].reverse();
    tpJointType lastJoint = segments_[0].jointAfter();
    for (i = 1; i< segments_.size(); ++i)
        segments_[i-1].setJointAfter(segments_[i].jointAfter());
    segments_[segments_.size()-1].setJointAfter(lastJoint);
}


//===========================================================================
inline void ftCurve::chopOff(const BoundingBox& box)
//===========================================================================
{
    std::vector<ftCurveSegment>::iterator it = segments_.begin();
    bool erase;
    while (it != segments_.end()) {
	std::vector<ftCurveSegment> new_segs = it->chopOff(box, erase);
	if (erase) {
	    it = segments_.erase(it);
	} else {
	    segments_.insert(it + 1, new_segs.begin(), new_segs.end());
	    it = it + new_segs.size() + 1;
	}
    }
}

//===========================================================================
inline int ftCurve::numSegments() const
//===========================================================================
{
    return (int)segments_.size();
}

//===========================================================================
inline double ftCurve::startOfSegment(int segment) const
//===========================================================================
{
    return segments_[segment].startOfSegment();
}

//===========================================================================
inline double ftCurve::endOfSegment(int segment) const
//===========================================================================
{
    return segments_[segment].endOfSegment();
}

//===========================================================================
inline const ftCurveSegment& ftCurve::segment(int segment) const
//===========================================================================
{
    return segments_[segment];
}

//===========================================================================
inline ftCurveType ftCurve::curveType() const
//===========================================================================
{
    return type_;
}

//===========================================================================
inline tpJointType ftCurve::jointAfter(int segment) const
//===========================================================================
{
    return segments_[segment].jointAfter();
}

//===========================================================================
inline tpJointType ftCurve::worstJointType() const
//===========================================================================
{
    tpJointType worst = JOINT_G1;   // is equal to 0
    for (size_t i = 0; i < segments_.size(); ++i)
	if ((int)(segments_[i].jointAfter()) > (int)worst)
	    worst = segments_[i].jointAfter();
    return worst;
}

//===========================================================================
inline int ftCurve::numDisjointSubcurves() const
//===========================================================================
{
    int n = 0;
    for (size_t i = 0; i < segments_.size(); ++i)
	if (segments_[i].jointAfter() >= JOINT_DISC)
	    ++n;
    return n;
}

//===========================================================================
inline void ftCurve::point(double t, int segment, Point& pt) const
//===========================================================================
{
    segments_[segment].point(t, pt);
}

//===========================================================================
inline void ftCurve::tangent(double t, int segment, Point& tan) const
//===========================================================================
{
    segments_[segment].tangent(t, tan);
}

//===========================================================================
inline void ftCurve::normal(double t, int segment, int side,
			    Point& normal, double eps) const
//===========================================================================
{
    segments_[segment].normal(t, side, normal, eps);
}



} // namespace Go


#endif // _FTCURVE_H

