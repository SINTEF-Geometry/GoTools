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
		   std::shared_ptr<ParamCurve> paramcv0,
		   std::shared_ptr<ParamCurve> paramcv1,
		   std::shared_ptr<ParamCurve> spacecv,
		   double eps_geo = 1.0e-6);

    inline ftCurveType segmentType() const;
    inline void setJointAfter(tpJointType jt);
    inline tpJointType jointAfter() const;
    inline ftFaceBase* face(int number) const;
    double startOfSegment() const;
    double endOfSegment() const;
    /** Evaluate a point on the curve. */
    void point(double t, Point& pt) const;
    /** Get the curve tangent in a given parameter. */
    void tangent(double t, Point& tan) const;
    /// Throws if the parameter plane curve asked for does not exist
    void paramcurvePoint(int number, double t, Point& pt) const;
    /// Throws if the parameter plane curve asked for does not exist
    void paramcurveTangent(int number, double t, Point& pt) const;
    void reverse();
    std::vector<ftCurveSegment> chopOff(const BoundingBox& box, bool& eraseme);

    std::shared_ptr<ParamCurve> spaceCurve() const;

    std::shared_ptr<ParamCurve> parameterCurve(int number) const;

    Point startPoint() const;
    Point endPoint() const;

    /** The normal of the surface corresponding to this curve. */
    void normal(double t, int side, Point& normal, double eps) const;
    double arcLength(double t1, double t2) const;
    
    /** The eps_go argument is used as a tolerance if ever this function has to 
     * regenerate the spacecurve from its underlying surface and parametric curve */
    void reparametrize(double eps_go = 1.0e-6);

    std::shared_ptr<LineStrip> tesselate(int resolution) const;

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
    tpJointType joint_;
    // Use 1 if edge/intersect, 2 if kink/gap
    ftFaceBase* underlying_face_[2];    // Surface(s)  on which curve is lying
    std::shared_ptr<ParamCurve> parameter_curve_[2]; // Curves in surface domain
    std::shared_ptr<ParamCurve> space_curve_;        // Curve in space

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
inline std::shared_ptr<ParamCurve> ftCurveSegment::spaceCurve() const
//===========================================================================
{
    return space_curve_;
}

//===========================================================================
inline std::shared_ptr<ParamCurve> ftCurveSegment::parameterCurve(int number) const
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
    ftCurve(ftCurveType t) : type_(t){}
    ftCurve& operator += (const ftCurve& cv);
    void appendCurve(const ftCurve& cv);
    void prependCurve(const ftCurve& cv);
    void appendSegment(const ftCurveSegment& seg);
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
    
    void reverse();
    void chopOff(const BoundingBox& box);

    // Curve info
    int numSegments() const;
    double startOfSegment(int segment) const;
    double endOfSegment(int segment) const;
    const ftCurveSegment& segment(int segment) const;
    ftCurveType curveType() const;
    tpJointType jointAfter(int segment) const;
    tpJointType worstJointType() const;
    int numDisjointSubcurves() const;

    // Evaluation
    void point(double t, int segment, Point& pt) const;
    void tangent(double t, int segment, Point& tan) const;
    void normal(double t, int segment, int side, Point& normal,
		double eps) const;

    /// Curve length
    double arcLength(double t1, int seg1, double t2, int seg2) const;

    /// Reparametrization
    /** The eps_go argument is used as a tolerance if ever this function has to 
     * regenerate the spacecurve from its underlying surface and parametric curve */
    void reparametrize(double eps_go = 1.0e-6);

    // Tesselation
    void tesselate(std::vector<std::shared_ptr<LineStrip> >& meshes) const;
    void tesselate(int resolution, 
		   std::vector<std::shared_ptr<LineStrip> >& meshes) const;
    void tesselate(double density, 
		   std::vector<std::shared_ptr<LineStrip> >& meshes) const;

    // For debugging
    void write(std::ostream& os) const;
    void writeSpaceCurve(std::ostream& os) const;
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

