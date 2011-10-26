//===========================================================================
//                                                                           
// File: ParamGeomInt.h                                                      
//                                                                           
// Created: Thu Sep 30 11:38:25 2004                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: ParamGeomInt.h,v 1.36 2007-06-22 16:33:21 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _PARAMGEOMINT_H
#define _PARAMGEOMINT_H


#include "GoTools/intersections/ParamObjectInt.h"
#include "GoTools/intersections/BoundaryGeomInt.h"
#include <memory>


namespace Go {


class ParamGeomInt;
class ParamPointInt;
class ParamCurveInt;
class ParamSurfaceInt;
class DirectionCone;


/// This class is a base class providing an interface to the
/// parametric "intersection objects".

class ParamGeomInt : public ParamObjectInt {
public:
    /// Constructor
    /// \param parent the parent to this intersection object.
    ParamGeomInt(ParamGeomInt* parent = 0) : ParamObjectInt(parent) {}

    /// Destructor
    virtual ~ParamGeomInt(){}

    /// The dimension of the geometric space.
    virtual int dimension() const = 0;

    /// Check if the object is periodic.  Analyze periodicity of curve
    /// based on number of repeating knots and control points. The
    /// return value is -1 if the curve ends are disjoint, otherwise k
    /// if cv is C^k continuous. These are sufficient but not
    /// necessary conditions for periodicity, so it is possible that a
    /// higher degree of periodicity exists.  Should not be called on
    /// this layer, should be overruled by inherited class.
    /// \param pardir the parameter direction in question.
    /// \return -1 if the curve ends are disjoint, or k if the curve
    /// is proven to be C^k continuous.
    virtual int checkPeriodicity(int pardir) const
    { return -1; } // Don't know about periodicity. Overrule when
		   // required

    /// Subdivide the object in the specified parameter direction and
    /// parameter value.
    /// \param pardir direction in which to subdive. Indexing starts
    /// at 0.
    /// \param par parameter in which to subdivide.
    /// \param subdiv_objs The subparts of this object. Of the same
    /// geometric dimension as this object.
    /// \param bd_objs the boundaries between the returned \a
    /// subdiv_objs. Of geometric dimension 1 less than this object.
    virtual void subdivide(int pardir, double par, 
			   std::vector<std::
			   shared_ptr<ParamGeomInt> >& subdiv_objs,
			   std::vector<std::
			   shared_ptr<ParamGeomInt> >& bd_objs) = 0;

    /// Return the CompositeBox for the parametric object.
    /// \return The compositeBox for the parametric object.
    virtual CompositeBox compositeBox() const = 0;

    /// A cone which contains all normals of the object.
    /// \return A cone which contains all normals of the object.
    virtual DirectionCone directionCone() const = 0;

    /// Check if the associated direction cone has an angle larger than pi
    bool coneLargerThanPi();

    /// Try to make a smaller direction cone.
    /// \param reduce_at_bd indicates which boundaries the reduction
    /// attempt should be made
    /// \param epsge the geometric tolerancedefining degeneracy
    /// \return A direction cone. This is either the same cone or a
    /// smaller one.
    virtual DirectionCone reducedDirectionCone(bool reduce_at_bd[4],
					       double epsge) const;

    /// Return the boundary objects of this object.
    /// \param bd_objs the boundary objects of this object.
    virtual void 
    getBoundaryObjects(std::vector<std::
		       shared_ptr<BoundaryGeomInt> >& bd_objs) = 0;

    /// Return the size of the geometric sample mesh in the specified
    /// direction.
    /// \param dir the parameter direction in question. Indexing
    /// starts at 0.
    virtual int getMeshSize(int dir) = 0;

    /// Return the corresponding mesh parameter.
    /// \param dir the parameter direction in question. Indexing
    /// starts at 0.
    /// \param idx the mesh idx in the specified direction. Indexing
    /// starts at 0.
    virtual double paramFromMesh(int dir, int idx) = 0;

    /// Return the geometric sample mesh for the parametric function.
    virtual std::vector<double>::iterator getMesh() = 0;

    /// Return the number of boundary objects.
    /// \return The number of boundary objects.
    int nmbBdObj() const
    { return (int)boundary_obj_.size(); }

    /// Return the specified boundary object.
    /// \param bd_idx index of the boundary object.
    /// \return The boundary object.
    BoundaryGeomInt* getBoundaryObject(int bd_idx) const
    {
	return (bd_idx < 0 || bd_idx >= int(boundary_obj_.size()))
	    ? 0 : boundary_obj_[bd_idx].get();
    }

    /// Return true if the object is degenerate in the specified
    /// direction.  Specialized for surface
    /// \param epsge the geometric tolerance defining degeneracy.
    /// \param dir the parameter direction in question.
    virtual int isDegenerate(double epsge, int dir)
    { return 0; }

    /// Return true if the object is degenerate in the specified
    /// direction and parameter.
    /// \param epsge the geometric tolerance defining degeneracy.
    /// \param dir the parameter direction in question.
    /// \param par the parameter in which to evaluate. Size of array
    /// should be equal to numParams().
    virtual bool isDegenerate(double epsge, int dir, double *par) = 0;

    /// Check the linearity of the object
    virtual bool isLinear(double epsge);

    /// Return whether or not the object is a spline.
    /// \return True if the object is a spline.
    virtual bool isSpline() = 0;

    /// We try to treat problems which will never result in a simple
    /// case by shrinking the domain slightly, resulting in smaller
    /// cones.  This is useful for scenarios where the normals are
    /// parallell in a boundary point.
    virtual double getOptimizedConeAngle(Point& axis1, Point& axis2) = 0;

protected:
    std::vector<std::shared_ptr<BoundaryGeomInt> > boundary_obj_;

};

    
} // namespace Go


#endif // _PARAMGEOMINT_H

