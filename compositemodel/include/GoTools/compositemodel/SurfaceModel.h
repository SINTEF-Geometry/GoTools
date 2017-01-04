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

#ifndef _SURFACEMODEL_H
#define _SURFACEMODEL_H

#include "GoTools/compositemodel/CompositeModel.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/compositemodel/ftFaceBase.h"
#include "GoTools/compositemodel/CellDivision.h"
//#include "GoTools/topology/tpTopologyTable.h"
#include "GoTools/compositemodel/ftCurve.h"
#include "GoTools/compositemodel/ftPoint.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/compositemodel/ftEdgeBase.h"
//#include "GoTools/compositemodel/Loop.h"
#include "GoTools/compositemodel/ftEdge.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/compositemodel/ftPlane.h"
#include "GoTools/compositemodel/ftLine.h"
#include "GoTools/compositemodel/FaceUtilities.h"
#include <vector>

namespace Go
{




 class ftPointSet;
 class IntResultsSfModel;
 class Loop;
 class Body;
 struct SamplePointData;

//===========================================================================
/** A surface set or shell including topological information
 */
// Note that the functions below may throw exceptions. More information will be
// added regarding the functions that may throw during implementation.
//
//===========================================================================

class GO_API SurfaceModel : public CompositeModel  
{
 public:
  /// Constructor taking a vector of faces
  /// \param approxtol Approximation error tolerance.
  /// \param gap If the distance between two points are less than 'gap' they
  ///            are viewed as identical.
  /// \param neighbour Maximum distance between surfaces viewed as adjacent.
  /// \param kink If two adjacent surfaces meet with an angle less than 'kink',
  ///             they are seen as G1 continous. (angles in radians)
  /// \param bend If two surfaces meet along a common boundary and corresponding
  ///             surface normals form an angle which is larger than 'bend', there 
  ///             is an intentional sharp edge between the surfaces.(angles in radians)
  /// \param faces A vector of faces.
  /// \param adjacency_set True if the application knows that twin information between edges is set.
  SurfaceModel(double approxtol,
	       double gap,   // Gap between adjacent surfaces
	       double neighbour,  // Threshold for whether surfaces are adjacent
	       double kink,  // Kink between adjacent surfaces 
	       double bend, // Intended G1 discontinuity between adjacent surfaces
	       std::vector<shared_ptr<ftSurface> >& faces, // Input faces
	       bool adjacency_set = false);   // If the application knows that twin
                                              // information between edges is set, a more 
                                              // simple topology analysis may be performed

  /// Constructor taking a vector of faces.
  /// \param faces A vector of faces.
  /// \param space_epsilon 
  /// \param kink 
  /// \param adjacency_set
  // @@@jbt - In many STEP-files only one tolerance is given. A
  // contructor which reflects this would be useful.
  // @@@vsk - The surface model checks itself for gaps and kinks. Thus, some kink tolerance
  // is needed. I also included the adjacency_set since it will often be the case for STEP 
  // files
  SurfaceModel(std::vector<shared_ptr<ftSurface> >& faces,
	       double space_epsilon, double kink = 0.01,
	       bool adjacency_set = false);

  /// Constructor taking a vector of parametric surfaces
  /// \param approxtol Approximation error tolerance. Not used.
  /// \param gap If the distance between two points is less than 'gap' they
  ///            are viewed as identical.
  /// \param neighbour Maximum distance between surfaces or curves viewed as adjacent.
  /// \param kink If two adjacent surfaces meet with an angle less than 'kink',
  ///             they are seen as G1 continous. (angles in radians)
  /// \param bend If two surfaces meet along a common boundary and corresponding
  ///             surface normals form an angle which is larger than 'bend', there 
  ///             is an intentional sharp edge between the surfaces.(angles in radians)
  /// \param surfaces A vector of surfaces.
  SurfaceModel(double approxtol,
	       double gap,   // Gap between adjacent surfaces
	       double neighbour,  // Threshold for whether surfaces are adjacent
	       double kink,  // Kink between adjacent surfaces 
	       double bend, // Intended G1 discontinuity between adjacent surfaces
	       std::vector<shared_ptr<ParamSurface> >& surfaces); // Input surfaces

 protected:
  SurfaceModel(double approxtol,
	       double gap,   // Gap between adjacent surfaces
	       double neighbour,  // Threshold for whether surfaces are adjacent
	       double kink,  // Kink between adjacent surfaces 
	       double bend);  // Intended G1 discontinuity between adjacent surfaces

 public:
  /// Constructor
  /// Create shallow copy of another SurfaceModel
  /// \param sm The SurfaceModel to be copied 
  SurfaceModel(const SurfaceModel& sm);

  /// Destructor
  virtual ~SurfaceModel();

  /// Set or reset topology tolerances
  /// \param approxtol Approximation error tolerance.
  /// \param gap If the distance between two points is less than 'gap' they
  ///            are viewed as identical.
  /// \param kink If two adjacent surfaces meet with an angle less than 'kink',
  ///             they are seen as G1 continous. (angles in radians)
  void setTolerances(double approxtol, double gap, double kink);

  /// Set or reset topology tolerances
  /// \param approxtol Approximation error tolerance. Not used.
  /// \param gap If the distance between two points is less than 'gap' they
  ///            are viewed as identical.
  /// \param neighbour Maximum distance between surfaces or curves viewed as adjacent.
  /// \param kink If two adjacent surfaces meet with an angle less than 'kink',
  ///             they are seen as G1 continous. (angles in radians)
  /// \param bend If two surfaces meet along a common boundary and corresponding
  ///             surface normals form an angle which is larger than 'bend', there 
  ///             is an intentional sharp edge between the surfaces.(angles in radians)
  void setTolerances(double approxtol, double gap, double neighbour, double kink, double bend);

  /// Return surface model pointer
  /// \return Pointer to this SurfaceModel
  virtual SurfaceModel* asSurfaceModel()
      {
	  return this;
      }

  /// Make a copy of the current model
  /// \return Pointer to a copy of this SurfaceModel
  virtual SurfaceModel* clone() const
  {
    return new SurfaceModel(*this);
  }

  /// Number of simple entities
  /// \return Number of simple entities
  virtual int nmbEntities() const;

  /// Return one face
  /// \param idx Index of face
  /// \return Pointer to face with index idx
  shared_ptr<ftSurface> getFace(int idx) const;

  /// Return all faces
  /// \return Vector of pointer to all faces
  std::vector<shared_ptr<ftSurface> > allFaces() const;

  /// Return one surface
  /// \param idx Index of surface
  /// \return Pointer to ParamSurface
  shared_ptr<ParamSurface> getSurface(int idx) const;

  /// Return one surface as SplineSurface if possible
  /// \param idx Index of surface
  /// \return Pointer to SplineSurface
  shared_ptr<SplineSurface> getSplineSurface(int idx) const;

  /// Given a face in the surface model, return the index of this face
  /// \param face Shared pointer to face
  /// \return Index to face
  int getIndex(shared_ptr<ftSurface> face) const;

  /// Given a face in the surface model, return the index of this face
  /// \param face Pointer to face
  /// \return Index to face
  int getIndex(ftSurface* face) const;

  /// Given a surface in the surface model, return the index of this face
  /// \param surf Pointer to surface
  /// \return Index to face
  int getIndex(ParamSurface* surf) const;

  /// Return a specified face as a shared pointer
  /// \param face Pointer to face
  /// \return Shared pointer to face
  shared_ptr<ftSurface> fetchAsSharedPtr(ftFaceBase *face) const;

  /// Switches the ordering of the faces internal in this construction.
  /// To be used if the order in which the faces are handled has a
  /// significance
  void swapFaces(int idx1, int idx2);

  /// Creates the CellDivision object
  void initializeCelldiv();

  /// Return a cell in the cell division
  /// \param i Index of cell
  /// \return The cell
  const ftCell& getCell(int i) const;

  /// Set limiting volume to mark area of interest
  void limitVolume(double xmin, double xmax,
		   double ymin, double ymax,
		   double zmin, double zmax);

  /// Check if a point is within the volume of interest.
  /// \param point The point to check.
  /// \return Whether the point is within the limits.
  bool pointWithinLimits(const ftPoint& point);

  /// Check if a point is within the volume of interest.
  /// \param point The point to check.
  /// \return Whether the point is within the limits.
  bool pointWithinLimits(const Point& point);

  /// Evaluate position
  /// \param idx Index of surface
  /// \param par[] Parameter value
  /// \param pnt Result
  virtual void evaluate(int idx,      // Index of surface
  			double par[], // Parameter value
			Point& pnt) const;  // Result


  /// Evaluate position and a number of derivatives
  /// The sequence is position, first derivative in first parameter direction,
  /// first derivative in second parameter direction,
  /// second derivative in first parameter direction, mixed second derivative etc.
  /// \param idx Index
  /// \param par Parameter value
  /// \param nder Number of derivatives to compute, 0=only position
  /// \param der Result
  virtual void evaluate(int idx,      // Index
  			double par[], // Parameter value
  			int nder,     // Number of derivatives to compute, 0=only position
  			std::vector<Point>& der) const;  // Result


  /// Closest point between a given point and this surface model
  /// Returns one point
  /// \param pnt Input point
  /// \param clo_pnt Found closest point
  /// \param idx Index of surface where the closest point is found
  /// \param clo_par[] Parameter value corresponding to the closest point
  /// \param dist Distance between input point and found closest point
    virtual void
    closestPoint(Point& pnt,     // Input point
  		 Point& clo_pnt, // Found closest point
  		 int& idx,          // Index of surface where the closest point is found
  		 double clo_par[],  // Parameter value corresponding to the closest point
  		 double& dist);     // Distance between input point and found closest point

  /// Closest point between a given point and this surface model
  /// \param point Input point
  /// \return Closest point
  ftPoint closestPoint(const Point& point);

  /// Closest point between a given point and this surface model
  /// \param point Input point
  /// \return Closest point
  ftPoint closestPoint(const ftPoint& point) { return closestPoint(point.position()); }


  /// Extremal point(s) in a given direction
  /// Note that the found extremal point may be less accurate for trimmed surfaces
  /// \param dir Direction 
  /// \param ext_pnt Found extremal point
  /// \param idx Index of surface where the extremal point is found
  /// \param ext_par[] Parameter value of extremal point
  virtual void extremalPoint(Point& dir, Point& ext_pnt, int& idx, double ext_par[]);

  /// Bounding box of the entire surface model
  /// \return Bounding box
  virtual BoundingBox boundingBox();

  /// Bounding box corresponding to one surface
  /// \param idx Index of surface
  /// \return Bounding box of surface with index idx.
  virtual BoundingBox boundingBox(int idx) const;  // Index of surface

  /// Whether one particular surface is degenerated such that one boundary
  /// degenerates to a point
  /// \param idx Index of surface
  /// \return Whether the surface is degenerated 
  virtual bool isDegenerate(int idx) const;

  /// Whether one particular surface is degenerated
  /// More specific information, relevant only for surfaces
  /// \param idx Index of surface
  /// \retval b The bottom curve degenerates to a point
  /// \retval t The top curve degenerates to a point
  /// \retval l The left curve degenerates to a point
  /// \retval r The right curve degenerates to a point
  /// \return Whether the surface is degenerated
  bool isDegenerate(int idx,  // Index of surface
		    bool& b,  // The bottom curve degenerates to a point
		    bool& t,  // The top curve degenerates to a point
		    bool& l,  // The left curve degenerates to a point
		    bool& r) const;  // The right curve degenerates to a point

  // Intersection with another SurfaceModel
  // Do we also want to be able to perform intersections between curves and surfaces?
  // Output is points and/or curves. Curves are represented as a vector of ftCurve, 
  // and points as ftPoint.Currently, ftCurve expects two surfaces as input and no curves.
  // Thus, it has to be extended to become more flexible. This will be done by letting ftCurve
  // be a superclass and add new subclasses that keep the pointers to the curves and surfaces
  // being input to the intersection function.
  // ftPoint has a pointer to one surface, not two. Thus, also this class must be extended.
  // Moreover, it does not expect curves. This will be handled the same way as with ftCurve.
  // Not yet implemented
/*   void intersect(shared_ptr<SurfaceModel>, // The other surface model */
/* 		 double tol,  // Is this input or class content? */
/* 		 std::vector<ftCurve>& int_curves, // Intersection curves, one curve may */
/* 		 // cross several of the surfaces in the model, but each curve is connected */
/* 		 // and simple. */
/* 		 std::vector<ftPoint>& int_points) const;  // Found intersection points */

  /// Intersection with a plane.
  /// \param plane The plane.
  /// \return Pointer to an IntResultsModel. 
     virtual shared_ptr<IntResultsModel> intersect_plane(const ftPlane& plane);

  /// Intersection with a line. Expected output is points, probably one point. Curves 
  /// can occur in special configurations
  /// \param line The line.
  /// \return Pointer to an IntResultsModel.
     virtual shared_ptr<IntResultsModel> intersect(const ftLine& line);

  /// Intersection with a line. Expected output is points, probably one point. Curves 
  /// can occur in special configurations
  /// \param line Consist of one point and one direction represented by Point.
  /// \retval int_curves Intersection curves, one curve may cross several of the
  ///                    surfaces in the model, but each curve is connected and simple.
  /// \retval int_points Found intersection points.
  void 
    intersect(const ftLine& line, // Consist of one point and one direction
	      // represented by Point. Just storage, not much content (yet)
	      ftCurve& int_curves, // Intersection curves, one curve may
	      // cross several of the surfaces in the model, but each curve is connected
	      // and simple.
	      std::vector<ftPoint>& int_points);  // Found intersection points

  /// Test if a line with direction 'dir' through the point 'point' hits this surface
  /// model. If it hits, return 'true' and the intersection point closest to 'point'.
  /// \param point Point on the line.
  /// \param dir Line direction.
  /// \retval result Closest intersection point.
  /// \return Whether the line hits or not.
  bool hit(const Point& point, const Point& dir, ftPoint& result);

/*   /// The two surface models are intersected and this model is trimmed with respect to the  */
/*   /// intersection result.  */
/*   void booleanIntersect(shared_ptr<SurfaceModel>, // The other model */
/* 			double tol);    // Tolerance */

  /** Intersect the surface model with a plane.
      \param plane The plane.
      \return Intersection curve.
  */
  ftCurve intersect(const ftPlane& plane);

  /** Intersect the model with a plane and trim this model with respect to the
      plane, the part of the model at the positive side of the plane is removed.
      \param plane The plane.
  */
  void booleanIntersect(const ftPlane& plane);

  /** Intersect the model with a plane and return the surface model trimmed with respect
      to this plane, the part of the model at the positive side of the plane is removed.
      \param plane The plane.
      \return Pointer to trimmed surface model
  */
  shared_ptr<SurfaceModel> trimWithPlane(const ftPlane& plane);

   /** Intersect the model with a line.
       \param line The intersecting line
       \retval represent_segment ???
       \return Vector of intersection points */
  std::vector<ftPoint> intersect(const ftLine& line, 
				 std::vector<bool>& represent_segment);

   /** Intersect the model with a SplineCurve.
       \param crv The intersecting SplineCurve.
       \retval represent_segment ???
       \return Vector of pairs of intersection points and their curve parameter value
   */
  std::vector<std::pair<ftPoint, double> > 
    intersect(shared_ptr<SplineCurve> crv,
	      std::vector<bool>& represent_segment);

  /// Split two surface models according to intersections between them.
  /// \param model2 The other model.
  /// \return Vector of new surface models.
  std::vector<shared_ptr<SurfaceModel> > 
    splitSurfaceModels(shared_ptr<SurfaceModel>& model2);

  /// Split a surface models according to a set of faces belonging to
  /// a specified surface model
  void
  splitSurfaceModel(std::vector<shared_ptr<ftSurface> >& faces,
		    Body* model2,
		    std::vector<std::vector<shared_ptr<ParamSurface> > >& result);

  /// Check if a spline surface intersects the current surface model
  bool doIntersect(shared_ptr<SplineSurface> sf);

  // Gaussian curvature
  // Not yet implemented
  //  virtual double curvature(int idx, // Index of surface
  //			   double *par) const;  // Parameter value at which to compute curvature
  virtual double curvature(int idx, double *par) const { return 0.0; }

/*   // Gaussian curvature */
/*   // Not yet implemented */
/*   double gaussCurvature(int idx, // Index of surface */
/* 			double par[]) const;  // Parameter value at which to compute curvature */
/*   // Mean curvature */
/*   // Not yet implemented */
/*   double meanCurvature(int idx, // Index of surface */
/* 		       double par[]) const;  // Parameter value at which to compute curvature */


  /// Turn parameter directions of one surface Will turn the orientation of the
  /// surface normal An update in the topology structures is required.
  /// \param idx Index of surface
  virtual void turn(int idx);

  /// Turn parameter directions of all surfaces  An update in the topology structures
  /// is required.
  virtual void turn();

  /// Append a new face to the surface model. The face is included in the topological
  /// structure
  // Should it be possible to append a SurfaceModel entitiy?
  // It could also be possible to append a ParamSurface?
  // This function is probably not virtual. It does not make sense to append a curve to
  // a surface model or a surface to a composite curve.
  /// \param face The new face
  /// \param set_twin If true, set twin face info. 
  void append(shared_ptr<ftSurface> face, bool set_twin = true,
	      bool adjacency_set = false, bool remove_twins = false,
	      int idx = -1);

  /// Append a vector of faces to the surface model. The faces are included in the topological
  /// structure
  /// \param faces Vector of pointers to the new faces
  void append(std::vector<shared_ptr<ftSurface> > faces, 
	      bool adjacency_set = false, bool set_twin = true);

  /// Append all faces from another surface model. The faces are included in the topological
  /// structure
  /// \param anotherModel Pointer to the other surface model
  void append(shared_ptr<SurfaceModel> anotherModel);

  /// Remove one face from the face set
  /// \param face Pointer to the face to be removed
  /// \return Whether the face was removed
  bool removeFace(shared_ptr<ftSurface> face);

  /// Update neighbourhood information related to face
  /// \param face Pointer to the face
  void updateFaceTopology(shared_ptr<ftSurface> face);

/*   // Join a new surface to an existing surface in the surface set. */
/*   // An update of the topology structure is performed. */
/*   // This function makes sense only in some combination of input. */
/*   // Both surfaces involved must be non-trimmed spline surfaces. Rational surfaces */
/*   // may be joined only if the weights satisfy certain conditions along the common */
/*   // boundary */
/*   // The surfaces should already be adjacent. */
/*   // Not yet implemented */
/*   void join(int idx,      // The index of the surface in this surface set involved in */
/* 	    // the operation */
/* 	    shared_ptr<ftSurface> face, // A pointer to the other surface (face) */
/* 	    int continuity);   // Expected continuity at the joint (0 = position only) */
/*   // Maximum continuity possible is 2 */

/*   // Adapt one surface to a set of points */
/*   // Not yet implemented */
/*   void adapt(int idx,   // Index of face to change */
/* 	     int cont,  // Continuity to maintain to adjacent faces, limited by existing */
/* 	     // continuity. Could be specified independent for each boundary */
/* 	     std::vector<shared_ptr<Point> >& points, // Points to adapt to */
/* 	     double approx_tol);  // Required accuracy in approximation */
	       
/*   // Adapt the surface model to a set of points */
/*   // Not yet implemented */
/*   void adapt(int fix,  // Number of derivatives to keep fixed at the outer boundary */
/* 	     // Could be specified independent for each boundary */
/* 	     double approx_tol);  // Required accuracy in approximation */

/*   // Other adapt functions can be with respect to curves, or points with normal vector */
/*   // Points sorted with respect to the surfaces in the surfaces set or parameterized */
/*   // points will simplify the function */
/*   // The adapt functions will be of the last ones to be implemented */

/*   // Closest point between two surface models. Not to be implemented yet */
/*   // More points can be given as output as for the other closest point functions */
/*   // Not yet implemented */
/*   void closestPoint(shared_ptr<SurfaceModel> other,  // The other surface model */
/* 		    int& idx1,  // Index of surface where the point lies in this model */
/* 		    Point &pnt1, // Closest point in this model */
/* 		    double par1[], // Parameter value of closest point in this model */
/* 		    int& idx2,  // Index of surface where the point lies in the other model */
/* 		    Point &pnt2, // Closest point in the other model */
/* 		    double par2[]) const; // Parameter value of closest point in the other model */

/*   // Return the shortest distance curve between two points in the surface model. The curve */
/*   // will follow the shortest path independent of surface boundararies. */
/*   // More specification is required. */
/*   // Not yet implemented */
/*   ftCurve shortestDistCurve(int idx1,      // Index of first surface */
/* 			    double par1[], // Parameter value in first surface */
/* 			    int idx2,      // Index of second surface */
/* 			    double par2[], // Parameter value in second surface */
/* 			    double approx_tol) const; // The maximum distance between the */
/*   // returned curve and the true surface curve. */

/*   // Draw. This function does not draw itself, but produces information which openGl can use */
/*   // to draw this surface model */
/*   // Commentet for avoiding compilation probalems when testing */
/*   // virtual void draw(/\* Some appropriate parameter list *\/) const; */
/*   virtual void draw() const { } */

  /// Tesselate surface model
  /// Tesselate all surfaces with respect to a default resolution
  /// \retval meshes Tesselated model
  virtual void tesselate(std::vector<shared_ptr<GeneralMesh> >& meshes) const;

  /// Tesselate all surfaces with respect to a given total resolution.
  /// The resolution in each parameter direction is set from the method
  /// \param uv_res Tesselation resolution
  /// \retval meshes Tesselated model
  void tesselate(int uv_res,
		 std::vector<shared_ptr<GeneralMesh> >& meshes) const;

  /// Tesselate all surfaces with respect to given resolutions in each
  /// parameter direction.
  /// \param resolution[] Tesselation resolution
  /// \retval meshes Tesselated model
  virtual
  void tesselate(int resolution[],
		 std::vector<shared_ptr<GeneralMesh> >& meshes) const;

  /// Tesselate all surfaces with respect to a given tesselation density
  /// \param density Tesselation density
  /// \retval meshes Tesselated model
  virtual
  void tesselate(double density,
		 std::vector<shared_ptr<GeneralMesh> >& meshes) const;

  /// Tesselate specified surfaces with respect to a given total resolution.
  /// The resolution in each parameter direction is set from the method
  /// \param faces Specified surfaces
  /// \param uv_res Tesselation resolution
  /// \retval meshes Tesselated surfaces
  void tesselate(const std::vector<shared_ptr<ftFaceBase> >& faces,
		 int uv_res,
		 std::vector<shared_ptr<GeneralMesh> >& meshes) const;

  /// Tesselate specified surfaces with respect to given resolutions in each
  /// parameter direction.
  /// \param faces Specified surfaces
  /// \param resolution[] Tesselation resolutions
  /// \retval meshes Tesselated model
  void tesselate(const std::vector<shared_ptr<ftFaceBase> >& faces,
		 int resolution[],
		 std::vector<shared_ptr<GeneralMesh> >& meshes) const;

  /// Tesselate specified surfaces with respect to a given tesselation density
  /// \param faces Specified surfaces
  /// \param density Tesselation density
  /// \retval meshes Tesselated model
  void tesselate(const std::vector<shared_ptr<ftFaceBase> >& faces,
		 double density,
		 std::vector<shared_ptr<GeneralMesh> >& meshes) const;

  /// Return a tesselation of the control polygon of all surfaces
  /// \retval ctr_pol Tesselation of the control polygon of all surfaces.
  virtual 
    void tesselatedCtrPolygon(std::vector<shared_ptr<LineCloud> >& ctr_pol) const;

  /// Return a tesselation of the control polygon of specified surfaces
  /// \param faces Specified surfaces
  /// \retval ctr_pol Tesselation of the control polygon of the specified surfaces.
  void tesselatedCtrPolygon(const std::vector<shared_ptr<ftFaceBase> >& faces,
			    std::vector<shared_ptr<LineCloud> >& ctr_pol) const;  

  /// Fetch sample points from all faces including surface normal and curvature
  /// information. The distribution is computed from the input parameter density
  void fetchSamplePoints(double density,
			 std::vector<SamplePointData>& sample_points) const;

  /// Set the elements in boundary_curves_ (based on top_table_).
  void setBoundaryCurves();

  /** Construct the topology information regarding the input geometry.
      \return Messages */
  ftMessage buildTopology(int first_idx=0, bool set_twin_face_info=true);

  /** Fetch the topology information from twin information in the input geometry.
      \return Messages */
  ftMessage setTopology();

  /// Add information about faces at the boundary meeting only in vertices
  void setVertexIdentity();

  /// Add information about twin faces
  void setTwinFaceInfo();

  /// Check if the SurfaceModel is a closed shell
  /// \return \a true if closed, \a false if open
  bool isClosed() const;

  /** Return the number of boundaries of a surface set (including holes). */
  int nmbBoundaries() const ;

  /** Return a given boundary (may be a hole). 
  \param whichbound refers to element number in boundary_curves_.
  \return boundary curve */ 
  ftCurve getBoundary(int whichbound);

  /** Return information about all gaps. */
  ftCurve getGaps();

  /** Information about gaps on an alternative format */
  void getGaps(std::vector<ftEdge*>& gaps);

  /** Return information about all kinks. */
  ftCurve getKinks();

  /** Information about kinks on an alternative format */
  void getKinks(std::vector<ftEdge*>& kinks);

  /** Return information about all G1 discontinuities. */
  ftCurve getG1Disconts();

  /** Information about kinks on an alternative format */
  void getCorners(std::vector<ftEdge*>& corners);

  /** Return single surface */
  ftSurface* getSurface2(int index) const;

  /** Return all compact face sets */
  std::vector<shared_ptr<SurfaceModel> > getConnectedModels() const;

  /** Return pointers to pairs of faces that have a problem with the consistency
      of the face orientation */
  void getInconsistentFacePairs(std::vector<std::pair<ftFaceBase*, ftFaceBase*> >& faces)
  {
    faces = inconsistent_orientation_;
  }

  /// Return pointers to pairs of overlapping edges.
  /// \param tol Overlap tolerance
  /// \retval edges Vector of pairs of overlapping edges
  void getOverlappingEdges(double tol,
			   std::vector<std::pair<shared_ptr<ftEdgeBase>, 
			   shared_ptr<ftEdgeBase> > >& edges);

  /** Triangulate the complete surface set with a prescribed point density. NB! The density
      may be adjusted if the number of points become too high and it is used as an indicator,
      not an absolute measure. */
  shared_ptr<ftPointSet> triangulate(double density) const;

  /// Return pointers to pairs of overlapping faces.
  /// \param tol Overlap tolerance
  /// \retval faces Vector of pairs of overlapping faces
  void getOverlappingFaces(double tol,
			   std::vector<std::pair<ftSurface*, ftSurface*> >& faces);

  /// Return all vertices associated with this surface model
  /// \retval vertices Vector of pointers to all vertices.
  void getAllVertices(std::vector<shared_ptr<Vertex> >& vertices) const;

  /// Fetch vertices at the boundaries
  /// \retval vertices Vector of pointers to vertices at the boundaries.
  void 
    getBoundaryVertices(std::vector<shared_ptr<Vertex> >& vertices) const;

  /** Fetch all edges at all the boundaries of this model
      \retval edges Vector of pointers to edges at the boundaries. */
  std::vector<shared_ptr<ftEdge> > getBoundaryEdges() const;

  /** Fetch all edges at one of the boundaries of this model
      \param boundary_idx Index of one boundary
      \retval edges Vector of pointers to the edges at one boundary.*/
  std::vector<shared_ptr<ftEdge> > getBoundaryEdges(int boundary_idx) const;

  /** Fetch all interval unique inner edges in this model, i.e.
      a ftEdge for each boundary edge with a twin, only one of the
      edges in the pair is returned.
      \retval Vector of pointers unique inner edges*/
  std::vector<shared_ptr<ftEdge> > getUniqueInnerEdges() const;

  /// Return body (if any)
  Body* getBody();

  /// Return approximation tolerance.
  /// \return Approximation tolerance.
  double getApproximationTol() const
  {
    return approxtol_;
  }

  /// Simplify fragmented trimming loops
  /// \retval max_dist Maximum local distorsion around the transitions. ??? 
  /// \return Whether any faces were modified. 
  bool simplifyTrimLoops(double& max_dist);

  /// Check if a surface model represents a rotational object and extract
  /// the eventual rotational axis and angle
  bool isAxisRotational(Point& centre, Point& axis, Point& vec, 
			double& angle, double& min_ang);

  /// Check if a surface model represents a simple linearily swept object
  /// that can be generated by sweeping a set of planar surfaces along a
  /// vector parallel to the plane normal. Return plane information and
  /// distance of sweep
  bool isLinearSwept(Point& pnt, Point& axis, double& len);

  /// Fetch all planar faces lying in a specified plane
  std::vector<shared_ptr<ftSurface> >  facesInPlane(Point& pnt, Point& axis);

  /// Check if all entities are NURBS
  /// \return Whether all entities are NURBS
  bool allSplines() const;

  /// Check if the model has got a corner-to-corner configuaration
  /// \return Whether the model has got a corner-to-corner configuaration
  bool isCornerToCorner() const;

  /// Ensure that the blocks in the model meet in a corner-to-corner
  /// configuration. NB! This function must not be applied to closed
  /// surface sets.
  void makeCornerToCorner();

  /// Ensure that the blocks in the model has got common spline spaces
  void makeCommonSplineSpaces();

  /// Modify adjacent spline surfaces from having almost co-linear
  /// coefficients to exact co-linearity
  void enforceCoLinearCoefs();

  /// Regularize face to mimic the division of a twin surface
  /// \param face
  /// \retval twinset
  void regularizeTwin(ftSurface *face, 
		      std::vector<shared_ptr<ftSurface> >& twinset);

  /// Merge two faces
  /// \return Pointer to resulting face
  shared_ptr<ftSurface> 
    mergeFaces(ftSurface* face1, int pardir1, double parval1,
	       bool atstart1, ftSurface* face2, int pardir2, 
	       double parval2, bool atstart2,
	       std::pair<Point,Point> co_par1, 
	       std::pair<Point,Point> co_par2,
	       std::vector<Point>& seam_joints);

  /// Merge two faces
  /// \return Pointer to resulting face
  shared_ptr<ftSurface> 
    mergeSeamFaces(ftSurface* face1, ftSurface* face2, int pardir,
		   std::vector<Point>& seam_joints);

  /// Merge two faces
  /// \return Pointer to resulting face
  shared_ptr<ftSurface> 
    mergeSeamCrvFaces(ftSurface* face1, ftSurface* face2, 
		      std::vector<Point>& seam_joints);

  /// Approximate regular trimmed surfaces with spline
  /// surfaces and replace
  void replaceRegularSurfaces();

  shared_ptr<ftSurface> replaceRegularSurface(ftSurface* face, 
					      bool only_corner=false);

  /// Simplify current shell by merging surfaces with a smooth connection
  /// both with regard to the surfaces themselves and the associated boundary
  /// curves. The common boundary must be isoparametric in both surfaces
  void simplifyShell();

  /// Approximate surface sets with 4 boundaries with a spline surface
  /// If the set has more than 4 corners, no surface will be produced
  shared_ptr<SplineSurface> approxFaceSet(double& error, int degree=3);

  /// Debug. Check topology
  bool checkShellTopology();

 protected:

  double approxtol_;
  double tol2d_;  // Tolerance to use for decisions in the parameter domain

  // Engine for the topology analysis
  //tpTopologyTable<ftEdgeBase, ftFaceBase> top_table_;
    
  // Storage of faces.
  // ftSurface is inherited from ftFaceBase. It might be that the faces_
  // vector should point to the base class, but the other children are 
  // very specialized and seem to be irrelevant here.
  std::vector<shared_ptr<ftFaceBase> > faces_;

  // For each separate object, we store all boundary loops
  // First element is (what is supposed to be) the objects outer boundary.
  std::vector<std::vector<shared_ptr<Loop> > > boundary_curves_;

  shared_ptr<CellDivision> celldiv_ ;   // To gain speedup in closest point and intersections
  mutable std::vector<bool> face_checked_;
  mutable int highest_face_checked_;
  //  mutable BoundingBox big_box_;
  BoundingBox limit_box_;

  std::vector<std::pair<ftFaceBase*, ftFaceBase*> > inconsistent_orientation_;

  void addSegment(ftCurve& cv, ftEdgeBase* edge, ftCurveType ty);

 private:

  void getCurveofType(ftCurveType type, ftCurve& curve);

  std::vector<ftCurveSegment> intersect(const ftPlane& plane, ftSurface* sf);
  ftCurve localIntersect(const ftPlane& plane, ftSurface* sf);

  void localIntersect(const ftLine& line, ftSurface* sf, 
		      std::vector<ftPoint>& result,
		      std::vector<ftCurveSegment>& line_segments) const;

  void localIntersect(shared_ptr<SplineCurve> crv,
		      ftSurface* sf,
		      std::vector<std::pair<ftPoint, double> >& result,
		      std::vector<ftCurveSegment>& crv_segments,
		      std::vector<std::pair<double,double> >& crv_bound) const;

  ftPoint closestPointLocal(const ftPoint& point) const;

  void localExtreme(ftSurface *face, Point& dir, 
		    Point& ext_pnt, int& ext_id,
		    double ext_par[]);

  void tesselateOneSrf(shared_ptr<ParamSurface> surf,
		       shared_ptr<GeneralMesh>& mesh,
		       int n=20, int m=20) const;

  void meshToTriang(shared_ptr<ftSurface> face,
		    shared_ptr<GeneralMesh> mesh,
		    int n, int m, shared_ptr<ftPointSet> triang,
		    bool check_endpoint_identity = false) const;

  void setResolutionFromDensity(shared_ptr<ParamSurface> surf,
				double density,
				int min_nmb, int max_nmb,
				int& u_res, int& v_res) const;

  void 
    getCurrConnectedModel(shared_ptr<ftSurface>& face,
			  std::vector<shared_ptr<ftSurface> >& curr_set,
			  std::vector<shared_ptr<ftSurface> >& all_sets) const;

  bool isInside(const Point& pnt, double& dist);

  shared_ptr<ftSurface> 
    performMergeFace(shared_ptr<ParamSurface> base,
		     vector<CurveLoop>& loops1, 
		     vector<CurveLoop>& loops2,
		     Body* bd,
		     std::vector<Point>& seam_joints,
		     int reverse, int cont=1);

  void makeCoonsBdCvs(vector<shared_ptr<ParamCurve> >& cvs1,
		      double tol, int degree,
		      vector<shared_ptr<ParamCurve> >& cvs2);

  // Find candidate faces to be merged for simplifyShell
  bool
    mergeSituation(ftSurface* face1, int& dir1, double& val1, bool& atstart1,
		   ftSurface* face2, int& dir2, double& val2, bool& atstart2,
		   pair<Point, Point>& co_par1, pair<Point, Point>& co_par2);
};


} // namespace Go


#endif // _SURFACEMODEL_H
