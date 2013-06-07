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

#ifndef _SPLINESURFACEINT_H
#define _SPLINESURFACEINT_H


#include "GoTools/intersections/ParamSurfaceInt.h"


namespace Go {


class SplineSurface;
class AlgObj3DInt;


/// This class represents the "intersection object" of a spline
/// surface.

class SplineSurfaceInt : public ParamSurfaceInt {
public:
    /// Constructor.

    /// \param surf the parametric curve defining the intersection
    /// object.
    explicit SplineSurfaceInt(shared_ptr<ParamSurface> surf);

    /// Constructor.
    /// \param surf the parametric curve defining the intersection
    /// object.
    /// \param parent the parent object to this object.
    explicit SplineSurfaceInt(shared_ptr<ParamSurface> surf,
			      ParamGeomInt *parent);

    /// Destructor.
    virtual ~SplineSurfaceInt(){};

    /// Return an intersection object for the input surface, using
    /// this object as parent.
    /// \param surf the parametric surface defining the intersection
    /// object.
    virtual shared_ptr<ParamSurfaceInt> 
    makeIntObject(shared_ptr<ParamSurface> surf);

    /// Return an intersection object for the input curve, using
    /// parameter parent as parent.
    /// \param crv the parametric curve defining the intersection
    /// object.
    /// \param parent the parent to the created intersection object.
    virtual shared_ptr<ParamCurveInt> 
    makeIntCurve(shared_ptr<ParamCurve> crv, ParamGeomInt* parent);

    /// Return true if the object has any inner knots in the specified
    /// parameter direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual bool hasInnerKnots(int pardir) const;

    /// Return true if the object has any critical parameter values or
    /// inner knots in the specified parameter direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual bool hasCriticalValsOrKnots(int pardir) const;

    /// Return the inner knot values in the specified direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    /// \param sort the returned values may be sorted by the function.
    virtual std::vector<double> getInnerKnotVals(int pardir,
						 bool sort=false) const;

    /// Return the critical parameter values and inner knots for the
    /// object.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual std::vector<double> getCriticalValsAndKnots(int pardir) const;

    /// Return the size of the geometric sample mesh in the specified
    /// direction.
    /// \param dir the parameter direction in question. Indexing
    /// starts at 0.
    virtual int getMeshSize(int dir);

    /// Return the corresponding mesh parameter.
    /// \param dir the parameter direction in question. Indexing
    /// starts at 0.
    /// \param idx the mesh idx in the specified direction. Indexing
    /// starts at 0.
    virtual double paramFromMesh(int dir, int idx);

    /// Return the geometric sample mesh for the parametric function.
    virtual std::vector<double>::iterator getMesh();

    /// A cone which contains all normals of the object.
    /// \return A cone which contains all normals of the object.
    virtual DirectionCone directionCone() const;    

    /// Make exact direction cone. Only applicable for Bezier case.
    virtual DirectionCone reducedDirectionCone(bool reduce_at_bd[4],
					       double epsge) const;

    /// Return the boundary objects of this object.
    /// \param bd_objs the boundary objects of this object.
    virtual void 
    getBoundaryObjects(std::vector<shared_ptr<BoundaryGeomInt> >& bd_objs);

    /// Check if the object is periodic in the specified direction.
    /// Analyze periodicity of surface based on number of repeating
    /// knots and control points. The return value is -1 if the
    /// surface ends are disjoint, otherwise k if cv is C^k
    /// continuous. These are sufficient but not necessary conditions
    /// for periodicity, so it is possible that a higher degree of
    /// periodicity exists.  Should not be called on this layer,
    /// should be overruled by inherited class.
    /// \param pardir the parameter direction in question.
    /// \return -1 if the curve ends are disjoint, or k if the surface
    /// is proven to be C^k continuous.
    virtual int checkPeriodicity(int pardir) const;

    /// Estimates if the current surface is simple enough for a
    /// singularity iteration. Checks the span of the normal cone and
    /// the size of the surface
    /// \return True if the surface is characterized as simple.
    virtual bool isSimple();

    /// Verify whether the object is a spline.
    /// \return True if the object is a spline.
    virtual bool isSpline();

    /// Check if a curve is an iso parametric curve in the current surface.
    /// NB! Only valid for splines
    virtual bool isIsoParametric(ParamCurveInt *curve, int dir, double par,
				 double ptol, double tol);

    /// We try to treat problems which will never result in a simple
    /// case by shrinking the domain slightly, resulting in smaller
    /// cones.  This is useful for scenarios where the normals are
    /// parallell in a boundary point.
    /// \param axis1 first vector defining a projection plane
    /// \param axis2 second vector defining a projection plane
    /// \return The optimized cone angle
    virtual double getOptimizedConeAngle(Point& axis1, Point& axis2);

    /// Find the knot intervals for which u and v lie inside, moving
    /// the value u or v if they lie close to a knot.
    /// \param u the \a u parameter value
    /// \param v the \a v parameter value
    /// \param utol the parametric tolerance deciding if the input
    /// parameter u should be moved.
    /// \param vtol the parametric tolerance deciding if the input
    /// parameter v should be moved.
    virtual void knotIntervalFuzzy(double& u, double&v,
				   double utol, double vtol) const; 

    /// Return the value of the knot next to the input parameter par.
    /// \param dir the parameter direction in question. Indexing
    /// starts at 0.
    /// \param par the parameter value
    /// \param forward if true we return the closest knot to the
    /// right, otherwise the closest knot to the left.
    /// \param tol the tolerance to determine if \a par is already
    /// located on the start of the next segment.
    /// \return The knot closest to the input parameter.
    virtual double nextSegmentVal(int dir, double par,
				  bool forward, double tol) const;

    /// Create a box containing the geometric sample mesh in the input
    /// coordinate system.
    /// \param axis the axis defining the coordinate system. The size
    /// of vector axis is 2, the last point is given as the cross
    /// product axis[0]%axis[1].
    /// \return The rotated box
    virtual RotatedBox getRotatedBox(std::vector<Point>& axis) const;

    /// Requests from selfintersection computation.  Split at G1
    /// discontinuities.
    /// \param angtol angular tolerance defining G1 discontinuity
    /// \param subG1 vector of subdivided patches that are G1
    virtual void splitAtG0(double angtol,
			   std::vector<shared_ptr<ParamSurfaceInt> >& subG1);

    /// Returns the normal surface corresponding to this surface.
    /// \return Pointer to the SplineSurface defining the normal
    /// surface.
    virtual shared_ptr<ParamSurfaceInt> getNormalSurface() const;

    /// Check whether the objects fulfills the requirements to
    /// implicitize.
    /// \return Return true if we may implicitize the object.
    virtual bool canImplicitize();

    /// Implicitize the object.
    /// \param tol the geoemtric tolerance for the implicitization.
    /// \return True if the implicitization was a success.
    virtual bool implicitize(double tol);

    /// Get the implicit representation of the object.
    /// Garbage is returned if we are not able to implicitize.
    /// \param tol geometric tolerance for the implicitization
    /// procedure.
    /// \param tol2 geometric estimate for the accuracy of the
    /// implicitized object.  Not yet in use!!!
    /// \param alg_obj_3d_int the algebraic object containing the
    /// implicitized surface.
    /// \return True if the implicitization was a success.
    // @@sbr Todo (tol2)!!!
    virtual bool getImplicit(double tol, double& tol2,
			     AlgObj3DInt& alg_obj_3d_int);

    /// Generate and return the specified isocurve.
    /// \param parameter value of the fixed parameter.
    /// \param pardir_is_u 'true' if the first parameter is the
    /// running parameter, 'false' otherwise.
    /// \return Pointer to the created isocurve.
    SplineCurve* constParamCurve(double parameter,
				 bool pardir_is_u) const;

    /// Return pointer to the spline surface defining this object.
    /// \return The spline surface defining this object.
    shared_ptr<SplineSurface> getSplineSurface()
    { return spsf_; }

    /// Return pointer to the spline surface defining this object.
    /// \return The spline surface defining this object.
    virtual shared_ptr<const SplineSurface> splineSurface() const
    { return spsf_; }

    /// Choose a degree for the implicit representation.
    void setImplicitDeg();

protected:
    // Data members
    shared_ptr<SplineSurface> spsf_;   // shared_ptr to
					      // this surface
    mutable shared_ptr<SplineSurface> normalsf_;  // Normal
							 // surface
							 // corresponding
							 // to this
							 // spline
							 // surface

private:

};


} // namespace Go


#endif // _SPLINESURFACEINT_H
