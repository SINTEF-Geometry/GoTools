//===========================================================================
//                                                                           
// File: CompositeModel.h                                                   
//                                                                           
// Created: February 2007
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _COMPOSITEMODEL_H
#define _COMPOSITEMODEL_H

#include "GoTools/utils/Point.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/topology/tpTolerances.h"
#include "GoTools/tesselator/GeneralMesh.h"
#include "GoTools/compositemodel/ftLine.h"
#include "GoTools/compositemodel/IntResultsModel.h"
#include <vector>
#include <set>


namespace Go
{

enum closestPointLevel
  {
    LOCAL_SEARCH = 1,   // Iterate for a closest point in each sub object
    SEMI_LOCAL_SEARCH,  // Use iteration with some pre processing
    GLOBAL_SEARCH       // Global closest point for each object. Slow, but secure
  };

 class SurfaceModel;
 class IntResultsModel;
 class LineCloud;
 
//===========================================================================
/** Abstract base class for a surface model or a curve model.
 */
//===========================================================================

class GO_API CompositeModel  
{
 public:
  /// Constructor with topology tolerances
  /// \param gap If the distance between two points are less than 'gap' they
  ///            are viewed as identical.
  /// \param neighbour Maximum distance between curves or surfaces viewed as adjacent.
  /// \param kink If two adjacent curves or surfaces meet with an angle less than
  ///             'kink', they are seen as G1 continous. (Angles in radians)
  /// \param bend If two surfaces meet along a common boundary and corresponding
  ///             surface normals form an angle which is larger than 'bend', there 
  ///             is an intentional sharp edge between the surfaces.
  ///             Similarily if two curves meet with an angle larger than 'bend', there 
  ///             is an intentional corner. (Angles in radians)
  CompositeModel(double gap, double neighbour, double kink, double bend);

  /// Destructor
  virtual ~CompositeModel();

  /// Set or reset topology tolerances
  /// \param gap If the distance between two points are less than 'gap' they
  ///            are viewed as identical.
  /// \param neighbour Maximum distance between curves or surfaces viewed as adjacent.
  /// \param kink If two adjacent curves or surfaces meet with an angle less than
  ///             'kink', they are seen as G1 continous. (Angles in radians)
  /// \param bend If two surfaces meet along a common boundary and corresponding
  ///             surface normals form an angle which is larger than 'bend', there 
  ///             is an intentional sharp edge between the surfaces.
  ///             Similarily if two curves meet with an angle larger than 'bend', there 
  ///             is an intentional corner. (Angles in radians)
  void setTolerances(double gap, double neighbour, double kink, double bend);

  /// Return topology tolerances
  tpTolerances getTolerances()
      {
	  return toptol_;
      }

  /// Return surface model pointer
  /// \return Pointer to surface model
  virtual SurfaceModel* asSurfaceModel()
      {
	  return 0;
      }

  /// Make a copy of the current model
  /// \return Pointer to new CompositeModel
  virtual CompositeModel* clone() const { return 0; }

  /// Number of simple entities
  /// \return Number of simple entities
  virtual int nmbEntities() const = 0;

  /// Evaluate position
  /// \param idx Index of curve
  /// \param par[] Parameter value
  /// \retval pnt Result
  virtual void evaluate(int idx,      // Index of surface
			double par[], // Parameter value
			Point& pnt) const = 0;  // Result

  /// Evaluate position and a number of derivatives
  /// \param idx Index
  /// \param par[] Parameter value
  /// \param nder Number of derivatives to compute, 0=only position
  /// \retval der Result
  virtual void evaluate(int idx,      // Index
			double par[], // Parameter value
			int nder,   // Number of derivatives to compute, 0=only position
			std::vector<Point>& der) const = 0;  // Result

  /// Compute one closest point
  /// \param pnt Input point
  /// \retval clo_pnt Found closest point
  /// \retval idx Index of curve or surface where the closest point is found
  /// \retval clo_par[] Parameter value corresponding to the closest point
  /// \retval dist Distance between input point and found closest point
  virtual void
    closestPoint(Point& pnt,     // Input point
		 Point& clo_pnt, // Found closest point
		 int& idx,         // Index of surface where the closest point is found
		 double clo_par[],  // Parameter value corresponding to the closest point
		 double& dist) = 0;  // Distance between input point and found closest point

  /// Intersection with a line. Expected output is points, probably one point. Curves 
  /// can occur in special configurations.
  /// \param line The line.
  /// \return Pointer to an IntResultsModel.  
     virtual shared_ptr<IntResultsModel> intersect(const ftLine& line) = 0;

  /// Intersection with a plane.
  /// \param plane The plane.
  /// \return Pointer to an IntResultsModel.  
     virtual shared_ptr<IntResultsModel> intersect_plane(const ftPlane& plane) = 0;

  /// Extremal point(s) in a given direction
  /// \param dir Direction 
  /// \param ext_pnt Found extremal point
  /// \param idx Index of curve or surface where the extremal point is found
  /// \param ext_par[] Parameter value of extremal point
  virtual void
    extremalPoint(Point& dir,     // Direction
		  Point& ext_pnt, // Found extremal point
		  int& idx,  // Index of curve or surface where the closest point is found
		  double ext_par[]) = 0;  // Parameter value of extremal point

  /// Bounding box of the entire model
  /// \return Bounding box
  virtual BoundingBox boundingBox() = 0;

  /// Bounding box corresponding to one entity
  /// \param idx Index of entity
  /// \return Bounding box
  virtual BoundingBox boundingBox(int idx) const = 0;  // Index of entity

  /// Whether one particular entity is degenerate
  /// \param idx Index of entity
  /// \return Whether the curve is degenerated
  virtual bool isDegenerate(int idx) const = 0;

  /// Curvature of an entity (Curve or Surface),
  /// Only Curve is implemented yet.
  /// \param idx Index of entity
  /// \param par Parameter value at which to compute curvature
  /// \return The curvature.
  virtual double curvature(int idx, // Index of entity
			   double *par) const = 0;  // Parameter value at which to compute curvature

  /// Turn parameter directions of one entity
  /// \param idx Index of entity
  virtual void turn(int idx) = 0;  // Turn parameter directions of one entity

  /// Turn parameter directions of all entities
  virtual void turn() = 0;   // Turn parameter directions of all entities

/*   // Draw. This function does not draw itself, but produces information which openGl can use */
/*   virtual void draw(/\* Some appropriate parameter list *\/) const = 0; */

  /// Tesselate model with respect to a default parameter
  /// \retval meshes Tesselated model
  virtual void tesselate(std::vector<shared_ptr<GeneralMesh> >& meshes) const = 0;

  /// Tesselate model with respect to a given resolution
  /// \param resolution[] Given resolution
  /// \retval meshes Tesselated model
  virtual
  void tesselate(int resolution[],
		 std::vector<shared_ptr<GeneralMesh> >& meshes) const = 0;

  /// Tesselate model with respect to a given tesselation density
  /// \param density Tesselation density
  /// \retval meshes Tesselated model
  virtual
  void tesselate(double density,
		 std::vector<shared_ptr<GeneralMesh> >& meshes) const = 0;

  /// Return a tesselation of the control polygon of this entity
  /// \retval ctr_pol Tesselated control polygon of this entity
  virtual 
    void tesselatedCtrPolygon(std::vector<shared_ptr<LineCloud> >& ctr_pol) const = 0;

 protected:
  // Tolerances used in topology analysis
  // These tolerances needs to be stored with the class as a topology
  // structure may become obsolete if the tolerances change
  tpTolerances toptol_;
  mutable int closest_idx_;

  double boxVecDist(const BoundingBox& b, const Point& v) const;

  bool boxExtreme(const BoundingBox& box, const Point& dir, 
		  const Point& curr_pnt) const;
};

} // namespace Go



#endif // _COMPOSITEMODEL_H

