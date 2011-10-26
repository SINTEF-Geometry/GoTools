//===========================================================================
//                                                                           
// File: SfPtIntersector.h 
//                                                                           
// Created: 
//                                                                           
// Author: vsk
//                                                                           
// Revision: $Id: SfPtIntersector.h,v 1.13 2006-11-03 14:43:14 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _SFPTINTERSECTOR_H
#define _SFPTINTERSECTOR_H


#include "GoTools/intersections/Intersector2Obj.h"


namespace Go {


/// This class performs intersection between a parametric surface and
/// a point

class SfPtIntersector : public Intersector2Obj {
public:

    /// Constructor.
    /// One of the objects should refer to a surface, the other a
    /// point (this is not checked compile-time, so we rely on the
    /// user to obey this rule).  The last two variables are relevant
    /// only if the parent has one more parameter than the Intersector
    /// to be constructed.
    /// \param obj1 either of type ParamSurfaceInt or ParamPointInt.
    /// \param obj2 either of type ParamPointInt or ParamSurfaceInt
    /// (not the same type as \a obj1).
    /// \param epsge the associated tolerance.
    /// \param prev the "parent" Intersector (0 if there is no
    /// parent).
    /// \param eliminated_parameter the index of the parameter that
    /// was removed from the parent \a prev.
    /// \param eliminated_value the value of the parameter that was
    /// removed from the parent \a prev.
    SfPtIntersector(std::shared_ptr<ParamGeomInt> obj1,
		    std::shared_ptr<ParamGeomInt> obj2,
		    std::shared_ptr<GeoTol> epsge, 
		    Intersector *prev = 0,
		    int eliminated_parameter = -1,
		    double eliminated_value = 0);

    /// Destructor.
    virtual ~SfPtIntersector();

//     // Validation of given intersection results
//     virtual void validate(int level, ValidationStat status);

    /// Return the number of parameter directions for the object.
    /// \return the number of parameter directions
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

    virtual int checkCoincidence();

    virtual void microCase();
  
    virtual int updateIntersections();

    virtual int repairIntersections()
    { return 0; }

    virtual int linearCase();

    virtual int doSubdivide();

    virtual int performRotatedBoxTest(double eps1, double eps2);

private:

    int sf_idx_, pt_idx_; // Indices to the surface object and the
			  // point object, respectivily

    void doIterate(double clo_par[2], double& clo_dist, double *guess = 0);

    int sortParameterDirections(int perm[]);

    SubdivisionClassification getSubdivisionParameter(int dir, double& par);

};


} // namespace Go


#endif  // _SFPTINTERSECTOR_H
