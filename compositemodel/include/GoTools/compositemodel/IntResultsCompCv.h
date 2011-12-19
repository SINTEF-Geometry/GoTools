//===========================================================================
//                                                                           
// File: IntResultsCompCv
//                                                                           
// Created: July 2009
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _INTRESULTSCOMPCV_H
#define _INTRESULTSCOMPCV_H

#include "GoTools/compositemodel/IntResultsModel.h"
#include "GoTools/geometry/PointOnCurve.h"


namespace Go
{

  class CompositeCurve;

//===========================================================================
/** Storage of results from intersection operations related to a CompositeCurve
*/
//
//===========================================================================

class IntResultsCompCv : public IntResultsModel
{
 public:

  /// Constructor
  // Note that only pointers are set to the initiating CompositeCurve.
  // If this curve is deleted, this pointer will be obsolete
  IntResultsCompCv(CompositeCurve *compcv, const ftLine& line);

  /// Constructor
  IntResultsCompCv(CompositeCurve *compcv, const ftPlane& plane);

  /// Destructor
  ~IntResultsCompCv();

  /// Add an intersection point
  void addIntPt(shared_ptr<ParamCurve> cv, double* parval);

  /// Add an intersection curve
  void addIntCv(shared_ptr<ParamCurve> cv, double* startpar, double* endpar);

  /// Check if any intersection curves are found
  virtual bool hasIntCurves() const
  {
    return (int_seg_1cv_.size() > 0);
  }

  /// Return the number of intersection curve segments
  virtual int nmbCurveSegments() const
  {
      return (int)int_seg_1cv_.size();
  }

  /// Check if any intersection point are found
  virtual bool hasIntPoints() const
  {
    return (int_pts_1cv_.size() > 0);
  }
  
  /// Return the number of intersection points
  virtual int nmbIntPoints() const 
  {
    return (int)int_pts_1cv_.size();
  }

  /// Fetch all intersection points
  void getIntersectionPoints(std::vector<PointOnCurve>& int_points) const;

  /// Fetch all intersection curves
  void 
    getIntersectionCurves(std::vector<std::pair<PointOnCurve, PointOnCurve> >& int_crvs) const;

  /// Check if any intersections are found
  virtual bool hasIntersections() const
  {
    return (int_seg_1cv_.size() > 0 || int_pts_1cv_.size() > 0);
  }

  /// Tesselation
  /// Tesselate with respect to a default value
  virtual void tesselate(std::vector<shared_ptr<LineStrip> >& meshes,
			 PointCloud3D& points) const;

  /// Tesselate with respect to a given resolution, the same for each
  /// intersection curve
  virtual void tesselate(int resolution, 
			 std::vector<shared_ptr<LineStrip> >& meshes,
			 PointCloud3D& points) const;

  /// Tesselate with respect to a given density
  virtual void tesselate(double density, 
			 std::vector<shared_ptr<LineStrip> >& meshes,
			 PointCloud3D& points) const;
 private:
  CompositeCurve *compcv1_;
  CompositeCurve *compcv2_;

  std::vector<PointOnCurve> int_pts_1cv_;
  std::vector<std::pair<PointOnCurve, PointOnCurve> > int_seg_1cv_;
};

} // namespace Go



#endif // _INTRESULTSCOMPCV_H

