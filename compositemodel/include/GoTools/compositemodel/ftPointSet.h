//===========================================================================
//                                                                           
// File: ftPointSet.h                                                       
//                                                                           
// Created: Wed Mar 14 2001                                         
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description: ftPointSet and helper ftSamplePoint classes
//                                                                           
//===========================================================================

#ifndef _FTPOINTSET_H
#define _FTPOINTSET_H


//===========================================================================
//===========================================================================

#include "GoTools/utils/Array.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/parametrization/PrOrganizedPoints.h"
#include <list>
#include <vector>             // Standard library STL vector

namespace Go
{

class ftSamplePoint;
 class ftFaceBase;
 class ftEdgeBase;
 class ftSurfaceSetPoint;

// We use a shared_ptr, thus allowing the use of sub classes of ftSamplePoint.
typedef std::list<shared_ptr<ftSamplePoint> > PointList;
typedef ftSamplePoint* PointIter;

/** ftSamplePoint -  One point in a set of sample points
 * 
 */
//===========================================================================
class ftSamplePoint
//===========================================================================
{
public:
    /// Constructor
    ftSamplePoint(Vector3D xyz, int bnd);
    /// Destructor
    virtual ~ftSamplePoint() {}

    /// Set point
    void setPoint(Vector3D xyz)
    { xyz_ = xyz; }
    /// Set parameter value
    void setPar(Vector2D uv)
    { uv_ = uv; }
    /// Set dist
    void setDist(double d)
    { dist_ = d; }
    /// Set index
    void setIndex(int i)
    { index_ = i; }

    /// Return pointer to sub class entity if the point is of that type
    virtual ftSurfaceSetPoint* asSurfaceSetPoint()
	{
	    return 0;
	}

    /// Check if a face is associated this point
    virtual bool containsFace(ftFaceBase* face) const
	{
	    return false;
	}

    /// Add a new neighbouring point
    void addNeighbour(PointIter next);
    /// If neighbour exists it is removed from neighbour vector.
    void removeNeighbour(PointIter neighbour);

    /// Am I on the boundary?
    bool isOnBoundary() const
    { return at_boundary_ == 1; }
    /// Am I on sub surface boundary?
    bool isOnSubSurfaceBoundary() const
    { return ((at_boundary_ == 2) || (at_boundary_ == 1)); }
    /// Get number of neighbours.
    int getNmbNeighbour() const
    { return (int)next_.size();}
    /// Fetch all neighbouring points
    const std::vector<PointIter>& getNeighbours() const
    { return next_; }

    /// Check if the point pnt is a neighbour to this point
    bool isConnected(PointIter pnt)
	{
	    for (size_t ki=0; ki<next_.size(); ++ki)
		if (next_[ki] == pnt)
		    return true;
	    return false;
	}

    /// Return 3D point value
    Vector3D getPoint() const
    { return xyz_; }
    /// Return parameter value
    Vector2D getPar() const
    { return uv_; }
    /// Return distance
    double getDist() const
    { return dist_; }
    /// Return index
    int getIndex() const
    { return index_; }
    /// Order the neighbours for one given point
    void orderNeighbours(ftSamplePoint* nextpoint, bool forward);
    /// Get first neighbour
    PointIter getFirstNeighbour()
      { return next_[0]; }

    /// Distance between sample points
    double pntDist(ftSamplePoint* other) const;

    /// Fetch all triangles containing this point
    void getAttachedTriangles(std::vector<std::vector<int> >& triangles) const;

    /// Debug
    virtual
      void write2Dval(std::ostream& os) const;
    
 protected:
    Vector3D xyz_;
    Vector2D uv_;
    double dist_;
    int index_;
    int at_boundary_; // 0: inner point, 1: boundary point (on merged surface),
                      // 2: boundary point on subface, inner point on merged surface.
    std::vector<PointIter> next_;


};  // End of ftSamplePoint
  

//===========================================================================
/** ftPointSet -  A set of sample points used as data in surface approximation
 * 
 */
//===========================================================================
class ftPointSet : public PrOrganizedPoints
//===========================================================================
//===========================================================================
{
public:
    /// Constructor
    ftPointSet();
    /// Destructor
    virtual ~ftPointSet();

    /// Get the number of points
    int size() const
    { return (int)index_to_iter_.size(); }

    /// Operator overload []
    ftSamplePoint* operator[](int idx)
    {
    	nb_ordered_ = false;
    	return (index_to_iter_[idx]);
    }
  
    /// Operator overload []
    const ftSamplePoint* operator[](int idx) const
    {
    	return (index_to_iter_[idx]);
    }
  
    /// Add a new point to the point set
    PointIter addEntry(shared_ptr<ftSamplePoint> point)
    {
	points_.push_back(point);
	points_.back()->setIndex((int)index_to_iter_.size());
	index_to_iter_.push_back((--points_.end())->get());
	nb_ordered_ = false;
	return (--points_.end())->get();
    }

    /// Remove a point from the point set
    void removePoint(PointIter point);

    /// Mark the first point around the boundary of this point set
    void setFirst(PointIter point)
    {
      first_ = point;
    }

    /// Mark the second point around the boundary of this point set
    void setSecond(PointIter point)
    {
      second_ = point;
    }

    /// Mark the last point added to this point set
    PointIter lastAdded()
    {
	return (--points_.end())->get();
    }

    /// Return the maximum distance in the pointset
    double getMaxDist() const;
    /// Return the medium distance in the pointset
    double getMeanDist() const;

    /// Compute the distances from the points in the point set
    /// to the given surface, at the same parameter value.
    void computeParametricDist(shared_ptr<ParamSurface> surf);
    /// Compute the distances from the points in the point set
    /// to the given surface.
    void computeDist(shared_ptr<ParamSurface> surf);
    /// Compute the distances from the points in the point set
    /// to the given surface, and reparametrize.
    void computeDistAndRepar(shared_ptr<ParamSurface> surf);

    /// Reparameterize points at the boundary of the point set
    double reparBdy(shared_ptr<ParamSurface> surf, bool use_seed = true);
    /// Reparameterize points in the inner of the point set
    double reparInnerPoints(shared_ptr<ParamSurface> surf, bool use_seed = true);

    /// Reorganize the sequence of neighbours
    void orderNeighbours();

    /// Extend the point set with points from another set. Avoid points with no
    /// connectivity
    void append(shared_ptr<ftPointSet> triang);

    /// Remove identical boundary nodes
    void cleanNodeIdentity(double tol);

    /// Given to faces, rearrange points on the common boundary between these
    /// two faces
    void mergeBoundary(shared_ptr<ftFaceBase> face1, int range1_idx1, 
		       int range1_idx2, shared_ptr<ftFaceBase> face2,
		       int range2_idx1, int range2_idx2, double eps);

    /// Fetch all triangles in the connectivity graph
    void getTriangles(std::vector<std::vector<int> >& triangles) const;

    /// Get the position of all points
    void getPoints(std::vector<Vector3D>& positions) const;

    /// Write point set to stream
    void write(std::ostream& os) const;

    /// Write parameter points to stream
     void write2D(std::ostream& os) const;

     /// Debug
    void printPoints(std::ostream& os) const;


    // From PrOrganizedPoints:
    /// Number of points in point set
    virtual int getNumNodes() const;
    /// Get content of point number i
    virtual Vector3D get3dNode(int i) const;
    /// Change content on point number i
    virtual void set3dNode(int i, const Vector3D& p);
    /// Fetch all neighbours to point number i
    virtual void getNeighbours(int i, std::vector<int>& neighbours) const;
    /// Check if point number i lies at the boundary of the point set
    virtual bool isBoundary(int i) const;

    /// Fetch 1. parameter of point number i
    virtual double getU(int i) const;
    /// Fetch 2. parameter of point number i
    virtual double getV(int i) const;
    /// Set 1. parameter of point number i
    virtual void setU(int i, double u);
    /// Set 2. parameter of point number i
    virtual void setV(int i, double v);
protected:
    PointList points_;
    std::vector<PointIter> index_to_iter_;
    mutable bool nb_ordered_;
    PointIter  first_;  // first_ and second_ are consecutive (CCW) points on
    PointIter  second_; // boundary. Defines direction.

 private:
    void addConnectivityInfo(PointIter pnt, PointIter pnt2, ftFaceBase* other_face);

    void mergeBoundaryEdges(std::vector<shared_ptr<ftEdgeBase> >& edges,
			    std::vector<shared_ptr<ParamCurve> >& crvs,
			    double tol) const;

};  // End of ftPointSet


} // namespace Go


#endif // _FTPOINTSET_H
