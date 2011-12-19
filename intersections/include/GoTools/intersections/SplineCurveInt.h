//===========================================================================
//                                                                           
// File: SplineCurveInt.h 
//                                                                           
// Created: 
//                                                                           
// Author: vsk
//                                                                           
// Revision: $Id: SplineCurveInt.h,v 1.20 2007-01-15 10:12:38 vsk Exp $
//                                                                           
// Description:
//
// NB! : TODO: Fjerne testutskrifter.  Sette aepsge, REL_PAR_RES og
// andre konstanter.  
// Se paa klassehierakiet. I ParamCurveInt.h skal vi ha
// baade checkCoincidence(... ParamSurfaceInt*...) og 
//       checkCoincidence(... SplineSurfaceInt*...) ?
//  min_step,max_step. 
// Se etter @bsp, @@ og @.
//                                                                           
//===========================================================================

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
