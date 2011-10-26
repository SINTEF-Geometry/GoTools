//===========================================================================
//                                                                           
// File: CvPtIntersector.h 
//                                                                           
// Created: 
//                                                                           
// Author: vsk
//                                                                           
// Revision: $Id: CvPtIntersector.h,v 1.21 2006-11-03 14:43:14 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _CVPTINTERSECTOR_H
#define _CVPTINTERSECTOR_H


#include "GoTools/intersections/Intersector2Obj.h"


namespace Go {


/// This class performs intersection between a parametric curve and a
/// point.

class CvPtIntersector : public Intersector2Obj {
public:

    /// Constructor.
    /// One of the objects should refer to a curve, the other a point
    /// (this is not checked compile-time, so we rely on the user to
    /// obey this rule).  The last two variables are relevant only if
    /// the parent has one more parameter than the Intersector to be
    /// constructed.
    /// \param obj1 either of type ParamCurveInt or ParamPointInt.
    /// \param obj2 either of type ParamPointInt or ParamCurveInt (not
    /// the same type as \a obj1).
    /// \param epsge the associated tolerance.
    /// \param prev the "parent" Intersector (0 if there is no
    /// parent).
    /// \param eliminated_parameter the index (0) of the parameter
    /// that was removed from the parent \a prev.
    /// \param eliminated_value the value of the parameter that was
    /// removed from the parent \a prev.
    CvPtIntersector(std::shared_ptr<ParamGeomInt> obj1,
		    std::shared_ptr<ParamGeomInt> obj2,
		    std::shared_ptr<GeoTol> epsge, 
		    Intersector *prev = 0,
		    int eliminated_parameter = -1,
		    double eliminated_value = 0);

    /// Constructor.
    /// One of the objects should refer a curve, the other a point
    /// (this is not checked compile-time, so we rely on the user to
    /// obey this rule).  The last two variables are relevant only if
    /// the parent has one more parameter than the Intersector to be
    /// constructed.
    /// \param obj1 either of type ParamCurveInt or ParamPointInt.
    /// \param obj2 either of type ParamPointInt or ParamCurveInt (not
    /// the same type as \a obj1).
    /// \param epsge the associated tolerance.
    /// \param prev the "parent" Intersector (0 if there is no
    /// parent).
    /// \param eliminated_parameter the index (0) of the parameter
    /// that was removed from the parent \a prev.
    /// \param eliminated_value the value of the parameter that was
    /// removed from the parent \a prev.
    CvPtIntersector(std::shared_ptr<ParamGeomInt> obj1,
		    std::shared_ptr<ParamGeomInt> obj2,
		    double epsge, 
		    Intersector *prev = 0,
		    int eliminated_parameter = -1,
		    double eliminated_value = 0);

    /// Destructor.
    virtual ~CvPtIntersector();

//     // Validation of given intersection results
//     virtual void validate(int level, ValidationStat status);

    /// Return the number of parameter directions for the object.
    /// \return the number of parameter directions
    virtual int numParams() const
    { return 1; }
	
protected:
    // Data members

    virtual std::shared_ptr<Intersector> 
    lowerOrderIntersector(std::shared_ptr<ParamGeomInt> obj1,
			  std::shared_ptr<ParamGeomInt> obj2, 
			  Intersector* prev = 0,
			  int eliminated_parameter = -1,
			  double eliminated_value = 0);

    virtual int checkCoincidence();

    virtual void microCase();
    
    virtual int updateIntersections();

    virtual int repairIntersections()
    { return 0; }

    virtual int linearCase();

    virtual int doSubdivide();

private:
    int cv_idx_, pt_idx_;  // Indices to the curve object and the
			   // point object, respectivily

    int sortParameterDirections();

    SubdivisionClassification getSubdivisionParameter(double& par);

};


} // namespace Go


#endif  // _CVPTINTERSECTOR_H
