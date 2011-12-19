//===========================================================================
//                                                                           
// File: SfCvIntersector.h 
//                                                                           
// Created: 
//                                                                           
// Author: vsk
//                                                                           
// Revision: $Id: SfCvIntersector.h,v 1.24 2007-11-01 14:31:37 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _SFCVINTERSECTOR_H
#define _SFCVINTERSECTOR_H


#include "GoTools/intersections/Intersector2Obj.h"


namespace Go {


/// This class performs intersection between a parametric surface and
/// a parametric curve.

class SfCvIntersector : public Intersector2Obj {
public:

    /// Constructor.
    /// One of the objects should refer to a surface, the other a
    /// curve (this is not checked compile-time, so we rely on the
    /// user to obey this rule).  The last two variables are relevant
    /// only if the parent has one more parameter than the Intersector
    /// to be constructed.
    /// \param obj1 either of type ParamSurfaceInt or ParamCurveInt.
    /// \param obj2 either of type ParamCurveInt or ParamSurfaceInt
    /// (not the same type as \a obj1).
    /// \param epsge the associated tolerance.
    /// \param prev the "parent" Intersector (0 if there is no
    /// parent).
    /// \param eliminated_parameter the index of the parameter that
    /// was removed from the parent \a prev.
    /// \param eliminated_value the value of the parameter that was
    /// removed from the parent \a prev.
    SfCvIntersector(shared_ptr<ParamGeomInt> obj1, 
		    shared_ptr<ParamGeomInt> obj2,
		    double epsge,
		    Intersector* prev = 0,
		    int eliminated_parameter = -1,
		    double eliminated_value = 0);

    /// Constructor.
    /// One of the objects should refer to a surface, the other a
    /// curve (this is not checked compile-time, so we rely on the
    /// user to obey this rule).  The last two variables are relevant
    /// only if the parent has one more parameter than the Intersector
    /// to be constructed.
    /// \param obj1 either of type ParamSurfaceInt or ParamCurveInt.
    /// \param obj2 either of type ParamCurveInt or ParamSurfaceInt
    /// (not the same type as \a obj1).
    /// \param epsge the associated tolerance.
    /// \param prev the "parent" Intersector (0 if there is no
    /// parent).
    /// \param eliminated_parameter the index of the parameter that
    /// was removed from the parent \a prev.
    /// \param eliminated_value the value of the parameter that was
    /// removed from the parent \a prev.
    SfCvIntersector(shared_ptr<ParamGeomInt> obj1, 
		    shared_ptr<ParamGeomInt> obj2,
		    shared_ptr<GeoTol> epsge,
		    Intersector* prev = 0,
		    int eliminated_parameter = -1,
		    double eliminated_value = 0);

    /// Destructor
    virtual ~SfCvIntersector();

//     // Validation of given intersection results
//     virtual void validate(int level, ValidationStat status);

    /// Return the number of parameter directions for the object.
    /// \return the number of parameter directions
    virtual int numParams() const
    { return 3; }

    void postIterateBd();

protected:
    // Data members

    virtual shared_ptr<Intersector> 
    lowerOrderIntersector(shared_ptr<ParamGeomInt> obj1,
			  shared_ptr<ParamGeomInt> obj2, 
			  Intersector* prev = 0,
			  int eliminated_parameter = -1,
			  double eliminated_value = 0);

    virtual int performRotatedBoxTest(double eps1, double eps2);

    virtual bool foundIntersectionNearBoundary();

    virtual int simpleCase();

    virtual int simpleCase2(Point& axis1, Point& axis2);

    virtual int checkCoincidence();

    virtual void microCase();
  
    virtual int updateIntersections();

    virtual int repairIntersections()
    { return 0; }

    virtual int linearCase();

    virtual int doSubdivide();

    virtual void postIterate(int nmb_orig, int dir=-1, bool keep_endpt=true);

private:

    int cv_idx_, sf_idx_;  // Indices to the curve object and the
			   // surface object, respectivily

    void doIterate(double par[], double& dist, double *guess=0);

   void doIterate2(double par[], double& dist, double guess[]);

    double distInCandidatePar(double par, int dir, const double* seed);

    int sortParameterDirections(int perm[]);

    SubdivisionClassification getSubdivisionParameter(int dir, double& par);

};


} // namespace Go


#endif  // _CVCVINTERSECTOR_H
