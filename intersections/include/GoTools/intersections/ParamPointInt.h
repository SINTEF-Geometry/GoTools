//===========================================================================
//                                                                           
// File: ParamPointInt.h 
//                                                                           
// Created: 
//                                                                           
// Author: vsk
//                                                                           
// Revision: $Id: ParamPointInt.h,v 1.37 2006-11-03 14:15:06 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _PARAMPOINTINT_H
#define _PARAMPOINTINT_H


#include "GoTools/intersections/ParamGeomInt.h"
#include "GoTools/utils/Values.h"


namespace Go {


class RotatedBox;


/// Class representing the "intersection object" of a parametric
/// curve.

class ParamPointInt : public ParamGeomInt {
public:
    /// Constructor.
    /// \param point the point defining the intersection
    /// object.
    /// \param parent the parent object to this object. Can be either
    /// a curve or a surface.
    explicit ParamPointInt(shared_ptr<Point> point,
			   ParamGeomInt* parent = 0);

    /// Destructor.
    virtual ~ParamPointInt() {};

    /// Evaluate the object in the input parameter.
    /// \param pt the Point to be returned.
    /// \param tpar the parameter in which to evaluate. The parameter
    /// is not used as the parameter domain of a point is of dimension
    /// 0.
    virtual void point(Point& pt, const double* tpar) const
    { pt = *(point_.get()); }

    /// Evaluate the object in the input parameter, with the specified
    /// number of derivatives.
    /// \param pt the Point to be returned.
    /// \param tpar the parameter in which to evaluate. The parameter
    /// is not used as the parameter domain of a point is of dimension
    /// 0.
    /// \param derivs the number of derivatives to calculate. Should
    /// be 0.
    /// \param from_right if true the evaluation is to be performed
    /// from the right side of the parameter value. Not used.
    /// \param resolution tolerance used when determining whether
    /// parameters are located at special values of the parameter
    /// domain (in particualar: knot values in case of spline
    /// objects. Not used.
    virtual void point(std::vector<Point>& pt, 
		       const double* tpar, 
		       int derivs,  
		       const bool* from_right = 0,
		       double resolution = 1.0e-12) const
    {
	DEBUG_ERROR_IF(derivs != 0,
		       "Asked for derivatives of a standalone point.");
	pt.resize(1);
	pt[0] = *(point_.get());
    }

    /// Return a pointer to this object.
    /// \return Pointer to this object.
    virtual ParamPointInt* getParamPointInt()
    { return this; }

    /// The number of parameters in the object.
    virtual int numParams() const
    { return 0; }

	
    /// Return an estimate on the size and wiggle of the object.  Does
    /// not apply to a point.
    /// \param length the approximative length of the object in the
    /// corresponding parameter directions.  The size of the array
    /// should be equal to numParams().
    /// \param wiggle a scalar representing the wiggle of the object
    /// in the corresponding parameter directions. The size of the
    /// array should be equal to numParams().
    virtual void getLengthAndWiggle(double *length, double *wiggle) {}

    /// Return true if the object has any inner knots in the specified
    /// parameter direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual bool hasInnerKnots(int pardir) const
    { return false;  }

    /// Return true if the object has any critical parameter values in
    /// the specified parameter direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual bool hasCriticalVals(int pardir) const
    { return false; }

    /// Return true if the object has any critical parameter values or
    /// inner knots in the specified parameter direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual bool hasCriticalValsOrKnots(int pardir) const
    { return false; }

    /// Return true if we are allowed to divide in the specified
    /// parameter direction.
    /// \param pardir the parameter direction in question.
    virtual bool canDivide(int pardir)
    { return false; }

    /// Return the critical parameter values in the specified direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual std::vector<double> getCriticalVals(int pardir) const
    {
	std::vector<double> seg;
	return seg;
    }

    /// Return the inner knot values in the specified direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    /// \param sort the returned values may be sorted by the function.
    virtual std::vector<double> getInnerKnotVals(int pardir, 
						 bool sort=false) const
    {
        std::vector<double> knots;
	return knots;
    }

    /// Return the critical parameter values and inner knots for
    /// object.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual std::vector<double> getCriticalValsAndKnots(int pardir) const
    {
	std::vector<double> vals;
	return vals;
    }

    /// Return the size of the geometric sample mesh in the specified
    /// direction.
    /// \param dir the parameter direction in question. Indexing
    /// starts at 0.
    virtual int getMeshSize(int dir)
    { return 1; }  // The point itself

    /// Return the corresponding mesh parameter.  Does not apply to a
    /// point.
    /// \param dir the parameter direction in question. Indexing
    /// starts at 0.
    /// \param idx the mesh idx in the specified direction. Indexing
    /// starts at 0.
    virtual double paramFromMesh(int dir, int idx)
    { return 0.0; }  // Arbitrary

    /// Return the geometric sample mesh for the parametric function.
    virtual std::vector<double>::iterator getMesh()
    { return coords_.begin(); }

    /// Return the start parameter value in the specified direction.
    /// Does not apply to a point.
    /// \param pardir the parameter direction in question. Indexing starts at 0.
    virtual double startParam(int pardir) const
    { return 0.0; }

    /// Return the end parameter in the specified direction.  Does not
    /// apply to a point.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual double endParam(int pardir) const
    { return 0.0; }

    /// Return true if the specified point lies within eps from the
    /// boundary.
    /// \param par the parameter in which to evaluate. Size of array
    /// should be equal to numParams().
    /// \param eps the tolerance defining boundary neighbourhood.
    virtual bool boundaryPoint(const double* par, double eps) const 
    { return true; }

    /// Subdivide the object in the specified parameter direction and
    /// parameter value.  Does not apply to a point.
    /// \param pardir direction in which to subdive. Indexing starts
    /// at 0.
    /// \param par parameter in which to subdivide.
    /// \param subdiv_objs The subparts of this object. Of the same
    /// geometric dimension as this object.
    /// \param bd_objs the boundaries between the returned \a
    /// subdiv_objs. Of geometric dimension 1 less than this object.
    virtual void 
    subdivide(int pardir, double par, 
	      std::vector<shared_ptr<ParamGeomInt> >& subdiv_objs,
	      std::vector<shared_ptr<ParamGeomInt> >& bd_objs) {}

    /// Create the CompositeBox for the parametric object.
    /// \return The CompositeBox for the parametric object.
    virtual CompositeBox compositeBox() const;

    /// A cone which contains all normals of the object.  As a normal
    /// is not defined for a point, the zero vector is returned.
    /// \return A cone which contains all normals of the object.
    virtual DirectionCone directionCone() const;

    /// Return the boundary objects of this object.
    /// \param bd_objs the boundary objects of this object.
    virtual void 
    getBoundaryObjects(std::vector<shared_ptr<BoundaryGeomInt> >& bd_objs);
       
    /// The dimension of the geometric space.
    virtual int dimension() const
    { return dim_; }

    /// Verify whether or not the object is degenerate in the
    /// specified direction and parameter.
    /// \param epsge the geometric tolerance defining degeneracy.
    /// \param dir the parameter direction in question.
    /// \param par the parameter in which to evaluate. Size of array
    /// should be equal to numParams().
    virtual bool isDegenerate(double epsge, int dir, double *par)
    { return false; }

    /// Verify whether the object is a spline.
    /// \return True if the object is a spline.
    virtual bool isSpline()
    { return false; }

    /// We try to treat problems which will never result in a simple
    /// case by shrinking the domain slightly, resulting in smaller
    /// cones.  This is useful for scenarios where the normals are
    /// parallell in a boundary point.
    /// \param axis1 first vector defining a projection plane
    /// \param axis2 second vector defining a projection plane
    /// \return The optimized cone angle
    virtual double getOptimizedConeAngle(Point& axis1, Point& axis2)
    { return M_PI; }

    /// Create a box containing the geometric sample mesh in the input
    /// coordinate system.
    /// \param axis the axis defining the coordinate system. The size
    /// of vector axis is 2, the last point is given as the cross
    /// product axis[0]%axis[1].
    RotatedBox getRotatedBox(std::vector<Point>& axis) const;

protected:
    // Data members

    shared_ptr<Point> point_;

    int dim_; // Space dimension.

    ParamCurveInt *parentcurve_;

    mutable std::vector<double> coords_;

};


} // namespace Go


#endif // _PARAMPOINTINT_H
