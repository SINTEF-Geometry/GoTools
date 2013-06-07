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

#ifndef _CURVEMODEL_H
#define _CURVEMODEL_H

#include "GoTools/compositemodel/CompositeModel.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/compositemodel/ftPoint.h"
#include "GoTools/compositemodel/ftCurve.h"
#include "GoTools/compositemodel/ftPlane.h"
#include "GoTools/compositemodel/ftLine.h"
#include <vector>


namespace Go
{

 class CompositeCurve;

//===========================================================================
/** A curve model including topological information
*/
// Note that the functions below may throw exceptions. More information will be
// added regarding the functions that may throw during implementation.
//
//===========================================================================

class CurveModel : public CompositeModel
  {
  public:

    /// Constructor
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
    CurveModel(double gap,   // Gap between adjacent curves
	       double neighbour,  // Threshold for whether curves are adjacent
	       double kink,  // Kink between adjacent curves
	       double bend, // Intended G1 discontinuity between adjacent curves
	       std::vector<shared_ptr<ParamCurve> >& curves); 

    /// Destructor
    ~CurveModel();

    /// Make a copy of the current model
    /// \return Pointer to a copy of this CurveModel
    virtual CurveModel* clone() const;

  /// Number of simple entities
  /// \return Number of simple entities
    virtual int nmbEntities() const;

  /// Return one curve
  /// Note that the index corresponds to the sequence of which the curves
  /// are added to the CurveModel, not the position in the composite
  /// curve
  /// \param idx Index of curve
  /// \return Pointer to curve
  shared_ptr<ParamCurve> getCurve(int idx) const;

  /// Given a curve in the composite curve, return the index of this curve
  /// \param curve Pointer to curve
  /// \return Index to curve
  int getIndex(ParamCurve* curve) const;

  /// Evaluate position
  /// \param idx Index of curve
  /// \param par[] Parameter value
  /// \retval pnt Result
  virtual void evaluate(int idx,      // Index of curve
			double par[], // Parameter value
			Point& pnt) const;  // Result

  /// Evaluate position and a number of derivatives
  /// The sequence is position, first derivative, second derivative, etc.
  /// \param idx Index
  /// \param par[] Parameter value
  /// \param nder Number of derivatives to compute, 0=only position
  /// \retval der Result
  virtual void evaluate(int idx,      // Index
			double par[], // Parameter value
			int nder,     // Number of derivatives to compute, 0=only position
			std::vector<Point>& der) const;  // Result

  /// Closest point between a given point and this composite curve
  /// Returns one point
  /// We could think of specifiying a function that returns all global closest points,
  /// only for global search. It is also possible to return closest points within a given
  /// range from the best closest point, but only one point for each curve in the model.
  /// \param pnt Input point
  /// \retval clo_pnt Found closest point
  /// \retval idx Index of curve where the closest point is found
  /// \retval clo_par[] Parameter value corresponding to the closest point
  /// \retval dist Distance between input point and found closest point
  virtual void
  closestPoint(Point& pnt,     // Input point
	       Point& clo_pnt, // Found closest point
	       int& idx,           // Index of curve where the closest point is found
	       double clo_par[],   // Parameter value corresponding to the closest point
	       double& dist);       // Distance between input point and found closest point

  /// Intersection with a line. Expected output is points, probably one point. Curves 
  /// can occur in special configurations.
  /// \param line The line.
  /// \return Pointer to an IntResultsModel. 
     virtual shared_ptr<IntResultsModel> intersect(const ftLine& line);

  /// Intersection with a plane.
  /// \param plane The plane.
  /// \return Pointer to an IntResultsModel. 
     virtual shared_ptr<IntResultsModel> intersect_plane(const ftPlane& plane);

    /// Extremal point(s) in a given direction.
    // Vector of points? Should point, index, parameter value and possibly distance
    // be stored in a struct?
    ///
    /// \param dir Direction 
    /// \retval ext_pnt Found extremal point
    /// \retval idx Index of curve where the extremal point is found
    /// \retval ext_par[] Parameter value of extremal point
    virtual void
    extremalPoint(Point& dir,     // Direction
		  Point& ext_pnt, // Found extremal point
		  int& idx,       // Index of curve where the extremal point is found
		  double ext_par[]);   // Parameter value of extremal point

    /// Bounding box of the entire composite curve
    /// \return Bounding box
    virtual BoundingBox boundingBox();

    /// Bounding box corresponding to one curve
    /// \param idx Index of curve
    /// \return Bounding box
    virtual BoundingBox boundingBox(int idx) const;  // Index of curve

    /// Whether one particular curve is degenerated
    /// \param idx Index of curve
    /// \return Whether the curve is degenerated 
    virtual bool isDegenerate(int idx) const;

    /// Curvature of a curve 
    /// \param idx Index of curve
    /// \param par Parameter value at which to compute curvature
    /// \return The curvature.
    virtual double curvature(int idx, // Index of curve
			     double *par) const;  // Parameter value at which to compute curvature

    /// Turn parameter direction of one curve. An update in the topology structures is required.
    /// \param idx Index of curve
    virtual void turn(int idx);  

   /// Turn parameter directions of all curves.  An update in the topology structures 
   /// is required.
    virtual void turn();  

  /// Tesselate all curves with respect to a default resolution
  /// \retval meshes Tesselated model
  virtual 
    void tesselate(std::vector<shared_ptr<GeneralMesh> >& meshes) const;

  /// Tesselate all curves with respect to a given resolution
  /// \param resolution[] All curves are tesselated with resolution[0]
  /// \retval meshes Tesselated model
  virtual
  void tesselate(int resolution[],
		 std::vector<shared_ptr<GeneralMesh> >& meshes) const;

  /// Tesselate all curves with respect to a given tesselation density
  /// \param density Tesselation density
  /// \retval meshes Tesselated model
  virtual
  void tesselate(double density,
		 std::vector<shared_ptr<GeneralMesh> >& meshes) const;

  /// Tesselate the control polygon of all curves.
  /// \retval ctr_pol Tesselation of the control polygon of all curves.
  virtual 
    void tesselatedCtrPolygon(std::vector<shared_ptr<LineCloud> >& ctr_pol) const;

  /// Fetch all uniquely connected composite curves
  /// \return Vector of pointers to the composite curves
  std::vector<shared_ptr<CompositeCurve> > fetchCompositeCurves() const;

private:
  std::vector<shared_ptr<ftEdge> > edges_;

  // Define connectivity between curves
  void buildTopology();


};

} // namespace Go



#endif // _CURVEMODEL_H
