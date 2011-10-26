//===========================================================================
//                                                                           
// File: IntResultsModel
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

#ifndef _INTRESULTSMODEL_H
#define _INTRESULTSMODEL_H

#include "GoTools/compositemodel/ftPoint.h"
#include "GoTools/compositemodel/ftCurve.h"
#include "GoTools/compositemodel/ftPlane.h"
#include "GoTools/compositemodel/ftLine.h"
#include "GoTools/geometry/PointCloud.h"
#include <vector>


namespace Go
{

enum IntersectionType
{
  No_Type = 0,
  SurfaceModel_SurfaceModel,
  SurfaceModel_Plane,
  SurfaceModel_Line,
  CompositeCurve_CompositeCurve,
  CompositeCurve_Plane,
  CompositeCurve_Line
};

class LineStrip;

//===========================================================================
/** Storage of results from intersection operations related to a CompositeModel
*/
//
//===========================================================================
class IntResultsModel
{
 public:

  /// Constructor
  IntResultsModel(IntersectionType type);

  /// Destructor
  ~IntResultsModel()
    {}

  /// Plane involved in the intersection
  void addPlaneInfo(const ftPlane& plane)
  {
    plane_ = plane;
  }

  /// Line involved in the intersection
  void addLineInfo(const ftLine& line)
  {
    line_ = line;
  }

  /// Check if any intersection curves are found
  virtual bool hasIntCurves() const = 0;

  /// Return the number of intersection curve segments
  virtual int nmbCurveSegments() const = 0;

  /// Check if any intersection point are found
  virtual bool hasIntPoints() const = 0;
  
  /// Return the number of intersection points
  virtual int nmbIntPoints() const = 0;

  /// Check if any intersections are found
  virtual bool hasIntersections() const = 0;
  
  /// Tesselation
  /// Tesselate with respect to a default value
  virtual void tesselate(std::vector<std::shared_ptr<LineStrip> >& meshes,
			 PointCloud3D& points) const = 0;

  /// Tesselate with respect to a given resolution, the same for each
  /// intersection curve
  virtual void tesselate(int resolution, 
			 std::vector<std::shared_ptr<LineStrip> >& meshes,
			 PointCloud3D& points) const = 0;

  /// Tesselate with respect to a given density
  virtual void tesselate(double density, 
			 std::vector<std::shared_ptr<LineStrip> >& meshes,
			 PointCloud3D& points) const = 0;

 protected:
  IntersectionType type_;
  ftLine line_;   // Is set if a line is involved in the current intersection
  ftPlane plane_; // Is set if a plane is involved in the current intersection
  int numpar_;    // Number of parameter values associated with an intersection result
};

} // namespace Go



#endif // _INTRESULTSMODEL_H
