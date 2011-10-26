//===========================================================================
//                                                                           
// File: CompositeCurve.h                                                    
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

#ifndef _COMPOSITECURVE_H
#define _COMPOSITECURVE_H

#include "GoTools/topology/tpJointType.h"
#include "GoTools/compositemodel/CompositeModel.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/compositemodel/ftPoint.h"
#include "GoTools/compositemodel/ftCurve.h"
#include "GoTools/compositemodel/ftPlane.h"
#include "GoTools/compositemodel/ftLine.h"
#include "GoTools/compositemodel/IntResultsCompCv.h"
#include <vector>


namespace Go
{
    

//===========================================================================
/** A composite curve including topological information
*/
// Note that the functions below may throw exceptions. More information will be
// added regarding the functions that may throw during implementation.
//
//===========================================================================

class CompositeCurve : public CompositeModel
  {
  public:

    /// Constructor.
    /// The sequence of the input curves will be kept, but a permutation array
    /// will sort the curves according to continuity
    /// \param gap If the distance between two points is less than 'gap' they
    ///            are viewed as identical.
    /// \param neighbour Maximum distance between curves viewed as adjacent.
    /// \param kink If two adjacent curves meet with an angle less than 'kink',
    ///             they are seen as G1 continous. (angles in radians)
    /// \param bend If two curves meet with an angle larger than 'bend', there 
    ///             is an intentional corner.(angles in radians)
    /// \param curves A vector of curves.
    CompositeCurve(double gap,   // Gap between adjacent curves
		   double neighbour,  // Threshold for whether curves are adjacent
		   double kink,  // Kink between adjacent curves
		   double bend, // Intended G1 discontinuity between adjacent curves
		   std::vector<std::shared_ptr<ParamCurve> >& curves); 

    /// Destructor
    ~CompositeCurve();

    /// Make a copy of the current model
    /// \return Pointer to the composite curve
      virtual CompositeCurve* clone() const;

  /// Number of simple entities
  /// \return Number of simple entities
      virtual int nmbEntities() const;

  /// Return one curve.
  /// Note that the index corresponds to the sequence of which the curves
  /// are added to the CompositeCurve, not the position in the composite
  /// curve
  /// \param idx Index of curve.
  /// \return Pointer to the curve.
  std::shared_ptr<ParamCurve> getCurve(int idx) const;

  /// Given a curve in the composite curve, return the index of this curve
  /// \param curve Pointer to curve
  /// \return Index to curve
  int getIndex(ParamCurve* curve) const;

  /// Get the parameter range.
  /// \retval start Parameter start. 
  /// \retval end Parameter end. 
  void parameterRange(double& start, double& end) const;

  /// Conversion between parameter settings
  /// \param idx Index of sub curve
  /// \param par Parameter in sub curve
  /// \return the global parameter corresponding to a sub curve parameter
  double getGlobalPar(int idx, double par) const;

  /// Return the sub curve parameter corresponding to a global parameter.
  /// At joints, the first sub curve parameter will be returned.
  /// \param global_par Parameter in composite curve
  /// \param idx Index of corresponding sub curve
  /// \retval local_par Parameter in sub curve
  void getLocalPar(double global_par, int& idx, double& local_par) const; 

  /// Evaluate position
  /// \param idx Index of curve
  /// \param par[] Parameter value
  /// \retval pnt Result
  virtual void evaluate(int idx, double par[], Point& pnt) const;

  /// Evaluate position and a number of derivatives
  /// The sequence is position, first derivative, second derivative, etc.
  /// \param idx Index of curve
  /// \param par Parameter value
  /// \param nder Number of derivatives to compute, 0=only position
  /// \retval der Result
  virtual void evaluate(int idx, double par[], int nder,
			std::vector<Point>& der) const;

  /// Evaluate with respect to the composite curve.
  /// Evaluate position
  /// \param par Parameter value
  /// \retval pnt Result
  void evaluateCurve(double par, Point& pnt) const;

  /// Evaluate position and a number of derivatives
  /// The sequence is position, first derivative, second derivative, etc.
  /// \param par Parameter value
  /// \param nder Number of derivatives to compute,  0=only position
  /// \retval der Result
  void evaluateCurve(double par, int nder, std::vector<Point>& der) const;

    /// Closest point between a given point and this composite curve
    /// Returns one point
    /// \param pnt Input point
    /// \retval clo_pnt Found closest point
    /// \retval idx Index of curve where the closest point is found
    /// \retval clo_par[] Parameter value corresponding to the closest point
    /// \retval dist Distance between input point and found closest point
    virtual void closestPoint(Point& pnt, Point& clo_pnt, int& idx,
			      double clo_par[], double& dist);

    /// Closest point between a given point and this composite curve
    /// Returns one point
    /// \param pnt Input point
    /// \return Closest point
    PointOnCurve closestPoint(Point& pnt);

    
    /// Extremal point(s) in a given direction
    /// \param dir Direction 
    /// \retval ext_pnt Found extremal point
    /// \retval idx Index of curve where the extremal point is found
    /// \retval ext_par[] Parameter value of extremal point
    virtual void
    extremalPoint(Point& dir, Point& ext_pnt, int& idx, double ext_par[]);

    /// Bounding box of the entire composite curve
    /// \return Bounding box
    virtual BoundingBox boundingBox();

    /// Bounding box corresponding to one curve
    /// \param idx Index of curve
    /// \return Bounding box of curve with index idx.
    virtual BoundingBox boundingBox(int idx) const;

    /// Whether one particular curve is degenerated
    /// \param idx Index of curve
    /// \return Whether the curve is degenerated
    virtual bool isDegenerate(int idx) const;

    /// Curvature of a curve 
    /// \param idx Index of curve
    /// \param par Parameter value at which to compute curvature
    /// \return The curvature.
    virtual double curvature(int idx, double *par) const;

    /// Turn parameter direction of one curve. An update in the topology structures is required.
    /// \param idx Index of curve
    virtual void turn(int idx);

   /// Turn parameter directions of all curves.  An update in the topology structures 
   /// is required.
   virtual void turn();  

/*     // Intersection with another composite curve */
/*     // Not implemented */
/*     void intersect(std::shared_ptr<CompositeCurve>, // The other composite curve */
/* 		  double tol,  // Is this input or class content? */
/* 		  std::vector<ftCurve>& int_curves, // Intersection curves, not the expected */
/* 		   // outcome, but it can occur */
/* 		  std::vector<ftPoint>& int_points) const;  // Found intersection points */

/*     /// Intersection with a plane */
/*      void intersect(ftPlane& plane, /// Intersection with a plane (3D) */
/* 		  double tol,  /// Tolerance */
/* 		  std::vector<ftCurve>& int_curves, // Intersection curves, not the expected */
/* 		   // outcome, but it can occur */
/* 		  std::vector<ftPoint>& int_points) const;  // Found intersection points */

  /// Intersection with a line. Expected output is points, probably one point. Curves 
  /// can occur in special configurations. The function makes most sense in 2D,
  /// otherwise it is very dependent on the tolerance.
  /// \param line The line.
  /// \return Pointer to an IntResultsModel. 
     virtual std::shared_ptr<IntResultsModel> intersect(const ftLine& line);

  /// Intersection with a plane.
  /// \param plane The plane.
  /// \return Pointer to an IntResultsModel.  
      virtual std::shared_ptr<IntResultsModel> intersect_plane(const ftPlane& plane);

  /// Test if a line with direction 'dir' through the point 'point' hits this composite
  /// curve. If it hits, return 'true' and the intersection point closest to 'point'.
  /// \param point Point on the line.
  /// \param dir Line direction.
  /// \retval result Closest intersection point.
  /// \return Whether the line hits or not.
      bool hit(const Point& point, const Point& dir, PointOnCurve& result);
 
   /// Append a new curve to the composite curve. The curve is included in the
   /// topological structure.
   /// \param curve The new curve.     
   void append(std::shared_ptr<ParamCurve> curve);

/*     // Join a new curve to an existing curve in the composite curve. */
/*     // An update of the topology structure is performed. */
/*     // This function makes sense only if the given parametric curves are of the same type. */
/*     // The curves should already be adjacent. */
/*     void join(int idx,      // The index of the curve in this composite curve involved in */
/*                 	      // the operation */
/* 	      std::shared_ptr<ParamCurve> curve, // A pointer to the other curve */
/* 	      int continuity);   // Expected continuity at the joint (0 = position only) */
/*                                  // Maximum continuity possible is 2 */

/*     // Adapt one curve to a set of points */
/*     void adapt(int idx,   // Index of curve to change */
/* 	       int cont,  // Continuity to maintain to adjacent curves, limited by existing */
/* 	                  // continuity. Could be specified independent for each boundary */
/* 	       std::vector<std::shared_ptr<Point> >& points, // Points to adapt to */
/* 	       double approx_tol);  // Required accuracy in approximation */
	       
/*    // Adapt the composite curve to a set of points */
/*     void adapt(int fix,  // Number of derivatives to keep fixed at the outer endpoints */
/* 	                  // Could be specified independent for each end */
/* 	       std::vector<std::shared_ptr<Point> >& points, // Points to adapt to */
/* 	       double approx_tol);  // Required accuracy in approximation */

    // Other adapt functions can be specified. Especially points parameterized with
    // respect to either the global composite curve or the segments can be of
    // interest

    // Closest point between two composite curves. 
/*     // More points can be given as output as for the other closest point functions */
/*     void closestPoint(std::shared_ptr<CompositeCurve> other,  // The other curve model */
/* 		      int& idx1,  // Index of curve where the point lies in this model */
/* 		      Point &pnt1, // Closest point in this model */
/* 		      double par1[], // Parameter value of closest point in this model */
/* 		      int& idx2,  // Index of curve where the point lies in the other model */
/* 		      Point &pnt2, // Closest point in the other model */
/* 		      double par2[], // Parameter value of closest point in the other model */
/* 		      closestPointLevel = GLOBAL_SEARCH) const; // Choice of computational method */

/*     // Advance a given distance along a composite curve. The distance parameter can be signed */
/*     // depending on whether we should advance in positive or negative parameter direction from */
/*     // the start point. */
/*     void advance(int idx1,  // Index of curve in start point */
/* 		 double par1,  // Parameter in start point */
/* 		 double advance_dist, // Distance to advance */
/* 		 int& idx2,  // Index of curve after advancing */
/* 		 double& par2) const;  // Parameter value of found point */

    //  Draw. These functions do not draw itself, but produces information which openGl
    //  can us to draw this model
    /// Tesselate with respect to a default parameter
    /// \retval meshes Tesselated model
     virtual 
     void tesselate(std::vector<std::shared_ptr<GeneralMesh> >& meshes) const;

    /// Tesselate with respect to a given resolution
    /// \param resolution[] All curves are tesselated with resolution[0].  
    /// \retval meshes Tesselated model
     virtual
     void tesselate(int resolution[],
		    std::vector<std::shared_ptr<GeneralMesh> >& meshes) const;

    /// Tesselate specified curves with respect to a given resolution
    /// \param curves Specified curves
    /// \param resolution[] Specified curves are tesselated with resolution[0].
    /// \retval meshes Tesselated model
     void tesselate(const std::vector<std::shared_ptr<ParamCurve> >& curves,
		    int resolution[],
		    std::vector<std::shared_ptr<GeneralMesh> >& meshes) const;

    /// Tesselate with respect to a given tesselation density
    /// \param density Tesselation density
    /// \retval meshes Tesselated model
     virtual
     void tesselate(double density,
		    std::vector<std::shared_ptr<GeneralMesh> >& meshes) const;

    /// Tesselate specified curves with respect to a given tesselation density
    /// \param curves Specified curves
    /// \param density Tesselation density
    /// \retval meshes Tesselated model
     void tesselate(const std::vector<std::shared_ptr<ParamCurve> >& curves,
		    double density,
		    std::vector<std::shared_ptr<GeneralMesh> >& meshes) const;

     /// Return a tesselation of the control polygon of this composite curve
     /// \retval ctr_pol Tesselated control polygon of this composite curve
  virtual 
    void tesselatedCtrPolygon(std::vector<std::shared_ptr<LineCloud> >& ctr_pol) const;

     /// Return a tesselation of the control polygon of specified curves
     /// \param curves Specified curves
     /// \retval ctr_pol Tesselated control polygon of specified curve
    void tesselatedCtrPolygon(const std::vector<std::shared_ptr<ParamCurve> >& curves,
			      std::vector<std::shared_ptr<LineCloud> >& ctr_pol) const;  

private:

  /// Storage of curves.
  std::vector<std::shared_ptr<ParamCurve> > curves_;

  /// Connectivity
  std::vector<tpJointType> continuity_;

  /// Permutation array
  std::vector<int> perm_;

  /// Whether or not a curve must be turned to make up a composite curve
  std::vector<bool> turned_;

  /// Representation of the parameterization of the complete composite curve. 
  /// This is a feature distinguishing the composite curve from a surface model.
  std::vector<double> start_parameters_;

  /// Find sequence of curves and set involved parameters
  void orderCurves();

  /// Intersect one curve with a line, 2D is assumed
  void 
    localIntersect(const ftLine& line, ParamCurve* cv, 
		   vector<double>& pt_par,
		   vector<double>& line_par);

  /// Intersect one curve with a plane
  void 
    localIntersect(const ftPlane& plane, ParamCurve* cv, 
		   vector<double>& pt_par,
		   vector<double>& line_par);

  /// Intersection between two curves
  void 
    localIntersect(ParamCurve *cv1, ParamCurve* cv2, 
		   vector<double>& pt_par,
		   vector<double>& line_par);

};

} // namespace Go



#endif // _COMPOSITECURVE_H
