//===========================================================================
//                                                                           
// File: IntResultsSfModel
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

#ifndef _INTRESULTSSFMODEL_H
#define _INTRESULTSSfMODEL_H

#include "GoTools/compositemodel/IntResultsModel.h"
#include "GoTools/geometry/PointCloud.h"


namespace Go
{

  class SurfaceModel;
  class LineStrip;

//===========================================================================
/** Storage of results from intersection operations related to a SurfaceModel
*/
//
//===========================================================================

class IntResultsSfModel : public IntResultsModel
{
 public:

  /// Constructor
  // Note that only pointers are set to the initiating SurfaceModel
  // If this model is deleted, this pointer will be obsolete
  IntResultsSfModel(SurfaceModel *sfmodel, const ftLine& line);

  /// Constructor
  IntResultsSfModel(SurfaceModel *sfmodel, const ftPlane& plane);

  /// Destructor
  ~IntResultsSfModel();

  /// Add an intersection point
  void addIntPts(std::vector<ftPoint>& intpts);

  /// Add an intersection curve
  void addIntCvs(ftCurve& cvs);

  /// Check if any intersection curves are found
  virtual bool hasIntCurves() const
  {
    return (intcvs_.numSegments() > 0);
  }
  
  /// Return the number of intersection curve segments
  virtual int nmbCurveSegments() const
  {
    return intcvs_.numSegments();
  }

  /// Fetch all intersection curve segments joined in one curve
  const ftCurve& getIntersectionCurves()
  {
    return intcvs_;
  }

  /// Check if any intersection points are found
  virtual bool hasIntPoints() const
  {
    return (int_pts_.size() > 0);
  }
  
  /// Return the number of intersection points
  int nmbIntPoints() const 
  {
    return (int)int_pts_.size();
  }

  /// Fetch all intersection points
  std::vector<ftPoint>& getIntersectionPoints()
    {
      return int_pts_;
    }

  /// Check if any intersections are found
  virtual bool hasIntersections() const
  {
    return (intcvs_.numSegments() > 0 || int_pts_.size() > 0);
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
  SurfaceModel *sfmodel1_;
  SurfaceModel *sfmodel2_;

  std::vector<ftPoint> int_pts_;
  ftCurve intcvs_;
};

} // namespace Go



#endif // _INTRESULTSCOMPCV_H

