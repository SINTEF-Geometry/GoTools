//===========================================================================
//                                                                           
// File: IntersectorAlgPar.h                                                 
//                                                                           
// Created: Wed Jan 26 15:20:41 2005                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: IntersectorAlgPar.h,v 1.13 2006-11-03 14:43:14 jbt Exp $
//                                                                           
// Description: Class for the intersection between algebraic and
// parametric object.
//                                                                           
//===========================================================================

#ifndef _INTERSECTORALGPAR_H
#define _INTERSECTORALGPAR_H


#include "GoTools/intersections/Intersector.h"


namespace Go {


class AlgObjectInt;
class ParamObjectInt;
class IntersectorFuncConst;
class SplineCurve;
class SplineSurface;
class AlgObj2DInt;
class AlgObj3DInt;


/// This class is an interface class used to compute the intersection
/// between an algebraic object and a parametric object.  The
/// algebraic object is "plugged" into the parametric object,
/// resulting in a scalar function.  We then use the
/// IntersectorFuncConst class to find the intersection.

class IntersectorAlgPar : public Intersector {
public:

    /// Constructor.
    /// The last two variables are relevant only if the parent has one
    /// more parameter than the Intersector to be constructed.
    /// \param alg_obj the algebraic object.
    /// \param param_obj the parametric object.
    /// \param epsge the associated tolerance.
    /// \param prev the "parent" Intersector (0 if there is no
    /// parent).
    /// \param eliminated_parameter the index of the parameter that
    /// was removed from the parent \a prev.
    /// \param eliminated_value the value of the parameter that was
    /// removed from the parent \a prev.
    IntersectorAlgPar(shared_ptr<AlgObjectInt> alg_obj,
		      shared_ptr<ParamObjectInt> param_obj,
		      shared_ptr<GeoTol> epsge,
		      Intersector* prev = 0,
		      int eliminated_parameter = -1,
		      double eliminated_value = 0);

    /// Destructor.
    virtual ~IntersectorAlgPar(){};

    /// Compute the current intersections (topology)
    /// \param compute_at_boundary flag to indicate that we compute at
    /// the boundary
    virtual void compute(bool compute_at_boundary = true);

    /// Return intersection points and curves.
    /// \param int_points vector of shared pointers to
    /// IntersectionPoints
    /// \param int_curves vector of shared pointers to
    /// IntersectionCurves
    // Sends request to IntersectionPool.
    virtual void
    getResult(std::vector<shared_ptr<IntersectionPoint> >& int_points,
	      std::vector<shared_ptr<IntersectionCurve> >& int_curves);

    /// Return the number of parameter directions for the object.
    /// \return thenumber of parameter directions
    virtual int numParams() const;

protected:

    shared_ptr<ParamObjectInt> param_int_;

    shared_ptr<AlgObjectInt> algobj_int_; // To be replaced by
						 // general algebraic
						 // object.

    // We create an object IntersectorFuncConst from the param_int_ &
    // algobj_int_.
    shared_ptr<IntersectorFuncConst> int_func_const_;

    virtual void print_objs();

    virtual int getBoundaryIntersections();

    virtual int performInterception();

    virtual int simpleCase();

    virtual bool isLinear();

    virtual bool complexityReduced();

    virtual void handleComplexity();

    virtual int checkCoincidence();
    
    virtual void microCase();
    
    virtual int updateIntersections();

    virtual int repairIntersections()
    { return 0; }

    virtual int linearCase();

    virtual int doSubdivide();

    virtual void printDebugInfo();

private:

    shared_ptr<IntersectorFuncConst> 
    insertCurveInAlgobj(SplineCurve* cv,
			AlgObj2DInt* alg_obj2d_int,
			shared_ptr<GeoTol> epsge,
			Intersector* intersector,
			int eliminated_parameter,
			double eliminated_value);

    shared_ptr<IntersectorFuncConst>
    insertSurfaceInAlgobj(SplineSurface* sf,
			  AlgObj3DInt* alg_obj3d_int,
			  shared_ptr<GeoTol>
			  epsge,
			  Intersector* intersector,
			  int eliminated_parameter,
			  double eliminated_value);

};


} // namespace Go


#endif // _INTERSECTORALGPAR_H

