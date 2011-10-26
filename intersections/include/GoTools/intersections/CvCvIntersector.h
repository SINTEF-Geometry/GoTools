//===========================================================================
//                                                                           
// File: CvCvIntersector.h 
//                                                                           
// Created: 
//                                                                           
// Author: vsk
//                                                                           
// Revision: $Id: CvCvIntersector.h,v 1.26 2006-11-16 17:13:12 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _CVCVINTERSECTOR_H
#define _CVCVINTERSECTOR_H


#include "GoTools/intersections/Intersector2Obj.h"


namespace Go {


/// Class that performs intersection between two parametric curves.

class CvCvIntersector : public Intersector2Obj {
public:

    /// Constructor.
    /// The last two variables pertains only if the parent has one
    /// more parameter than the Intersector to be constructed.
    /// \param curve1 of type ParamCurveInt.
    /// \param curve2 of type ParamCurveInt.
    /// \param epsge the associated tolerance.
    /// \param prev the "parent" Intersector (0 if there is no
    /// parent).
    /// \param eliminated_parameter the index of the parameter that
    /// was removed from the parent \a prev.
    /// \param eliminated_value the value of the parameter that was
    /// removed from the parent \a prev.
    CvCvIntersector(std::shared_ptr<ParamGeomInt> curve1, 
		    std::shared_ptr<ParamGeomInt> curve2,
		    double epsge,
		    Intersector* prev = 0,
		    int eliminated_parameter = -1,
		    double eliminated_value = 0);

    /// Constructor.
    /// The last two variables pertains only if the parent has one
    /// more parameter than the Intersector to be constructed.
    /// \param curve1 of type ParamCurveInt.
    /// \param curve2 of type ParamCurveInt.
    /// \param epsge the associated tolerance.
    /// \param prev the "parent" Intersector (0 if there is no
    /// parent).
    /// \param eliminated_parameter the index of the parameter that
    /// was removed from the parent \a prev.
    /// \param eliminated_value the value of the parameter that was
    /// removed from the parent \a prev.
    CvCvIntersector(std::shared_ptr<ParamGeomInt> curve1, 
		    std::shared_ptr<ParamGeomInt> curve2,
		    std::shared_ptr<GeoTol> epsge, 
		    Intersector* prev = 0,
		    int eliminated_parameter = -1,
		    double eliminated_value = 0);

    /// Destructor.
    virtual ~CvCvIntersector();

//     // Validation of given intersection results
//     virtual void validate(int level, ValidationStat status);

    /// Return the number of parameter directions for the
    /// intersection.
    /// \return The number of parameter directions for the
    /// intersection.
    virtual int numParams() const
    { return 2; }
	
protected:
    // Data members

    virtual std::shared_ptr<Intersector> 
    lowerOrderIntersector(std::shared_ptr<ParamGeomInt> obj1,
			  std::shared_ptr<ParamGeomInt> obj2, 
			  Intersector* prev = 0,
			  int eliminated_parameter = -1,
			  double eliminated_value = 0);

    virtual int performRotatedBoxTest(double eps1, double eps2);

    virtual bool foundIntersectionNearBoundary();

    virtual int simpleCase2(Point& axis1, Point& axis2);

    virtual int checkCoincidence();

    virtual void microCase();
  
    virtual int updateIntersections();

    virtual int repairIntersections()
    { return 0; }

    virtual int linearCase();

    virtual int doSubdivide();

    virtual void writeOut();

private:
    void doIterate(double& par1, double& par2,
		   double& dist, double *guess=0);

    double distInCandidatePar(double par, int dir, const double* seed);

    int sortParameterDirections(int perm[]);

    SubdivisionClassification getSubdivisionParameter(int dir,
						      double& par);

};


} // namespace Go


#endif  // _CVCVINTERSECTOR_H
