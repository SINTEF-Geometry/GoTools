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

#ifndef _SPLINECURVEINT_H
#define _SPLINECURVEINT_H


#include "GoTools/intersections/ParamCurveInt.h"


namespace Go {


class SplineCurve;


/// Class representing the "intersection object" of a spline curve.

class SplineCurveInt : public ParamCurveInt {
public:
    /// Constructor.
    /// \param curve the parametric curve defining the intersection
    /// object.
    explicit SplineCurveInt(shared_ptr<ParamCurve> curve);

    /// Constructor.
    /// \param curve the parametric curve defining the intersection
    /// object.
    /// \param parent the parent object to this object. Can be either
    /// a curve or a surface.
    explicit SplineCurveInt(shared_ptr<ParamCurve> curve, 
			    ParamGeomInt *parent);

    /// Destructor.
    virtual ~SplineCurveInt(){};

    /// Return an intersection object for the input curve.  This
    /// object is used as the parent for the intersection object.
    /// \param curve the parametric curve defining the intersection
    /// object.
    virtual shared_ptr<ParamCurveInt> 
    makeIntObject(shared_ptr<ParamCurve> curve);

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

    /// Return the critical parameter values and inner knots for
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

    /// Check if the object is periodic.  Analyze periodicity of curve
    /// based on number of repeating knots and control points. The
    /// return value is -1 if the curve ends are disjoint, otherwise k
    /// if cv is C^k continuous. These are sufficient but not
    /// necessary conditions for periodicity, so it is possible that a
    /// higher degree of periodicity exists.  Should not be called on
    /// this layer, should be overruled by inherited class.
    /// \param pardir the parameter direction in question. Does not
    /// pertain to for a curve.
    /// \return -1 if the curve ends are disjoint, or k if the curve
    /// is proven to be C^k continuous.
    virtual int checkPeriodicity(int pardir = 0) const;

    /// Find the knot interval for which t lies inside, moving the
    /// value t if it lies close to a knot.
    /// \param t the parameter value
    /// \param tol the parametric tolerance deciding if the input
    /// parameter t should be moved.
    virtual int knotIntervalFuzzy(double& t, double tol) const;

    /// Return the value of the knot next to the input parameter par.
    /// \param par the parameter value
    /// \param forward if true we return the closest knot to the
    /// right, otherwise the closest knot to the left.
    /// \param tol the tolerance to determine if \a par is already
    /// located on the start of the next segment.
    /// \return The knot closest to the input parameter.
    virtual double nextSegmentVal(double par, bool forward, double tol) const;

    /// Verify whether the object is a spline.
    /// \return True if the object is a spline.
    virtual bool isSpline();

    virtual const SplineCurve* getSpline()
	{ return spcv_.get(); }

    /// We try to treat problems which will never result in a simple
    /// case by shrinking the domain slightly, resulting in smaller
    /// cones.  This is useful for scenarios where the normals are
    /// parallell in a boundary point.
    /// \param axis1 first vector defining a projection plane
    /// \param axis2 second vector defining a projection plane
    /// \return The optimized cone angle
    virtual double getOptimizedConeAngle(Point& axis1, Point& axis2);

    // Functions used for coincidence testint

    /// Create a box containing the geometric sample mesh in the input
    /// coordinate system.
    /// \param axis the axis defining the coordinate system. The size
    /// of vector axis is 2, the last point is given as the cross
    /// product axis[0]%axis[1].
    /// \return The rotated box
    virtual RotatedBox getRotatedBox(std::vector<Point>& axis) const;

protected:
    // Data members
    shared_ptr<SplineCurve> spcv_; // shared_ptr to this curve
    
};


} // namespace Go


#endif // _SPLINECURVEINT_H
