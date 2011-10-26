//===========================================================================
//                                                                           
// File: PtPtIntersector.h 
//                                                                           
// Created: 
//                                                                           
// Author: vsk
//                                                                           
// Revision: $Id: PtPtIntersector.h,v 1.14 2006-11-03 14:43:14 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _PTPTINTERSECTOR_H
#define _PTPTINTERSECTOR_H


#include "GoTools/intersections/Intersector2Obj.h"


namespace Go {


/// This class performs intersection between two points.

class PtPtIntersector : public Intersector2Obj {
public:

    /// Constructor.
    /// The last two variables are relevant only if the parent has one
    /// more parameter than the Intersector to be constructed.
    /// \param point1 of type ParamPointInt.
    /// \param point2 of type ParamPointInt.
    /// \param epsge the associated tolerance.
    /// \param prev the "parent" Intersector (0 if there is no
    /// parent).
    /// \param eliminated_parameter the index (0) of the parameter
    /// that was removed from the parent \a prev.
    /// \param eliminated_value the value of the parameter that was
    /// removed from the parent \a prev.
    PtPtIntersector(std::shared_ptr<ParamGeomInt> point1, 
		    std::shared_ptr<ParamGeomInt> point2,
		    std::shared_ptr<GeoTol> epsge, 
		    Intersector* prev = 0,
		    int eliminated_parameter = -1,
		    double eliminated_value = 0);

    /// Constructor.
    /// The last two variables are relevant only if the parent has one
    /// more parameter than the Intersector to be constructed.
    /// \param point1 of type ParamPointInt.
    /// \param point2 of type ParamPointInt.
    /// \param epsge the associated tolerance.
    /// \param prev the "parent" Intersector (0 if there is no
    /// parent).
    /// \param eliminated_parameter the index (0) of the parameter
    /// that was removed from the parent \a prev.
    /// \param eliminated_value the value of the parameter that was
    /// removed from the parent \a prev.
    PtPtIntersector(std::shared_ptr<ParamGeomInt> point1, 
		    std::shared_ptr<ParamGeomInt> point2,
		    double epsge, 
		    Intersector* prev = 0,
		    int eliminated_parameter = -1,
		    double eliminated_value = 0);


    /// Destructor
    virtual ~PtPtIntersector();

//     // Validation of given intersection results
//     virtual void validate(int level, ValidationStat status);

    /// Return the number of parameter directions for the object.
    /// \return the number of parameter directions
    virtual int numParams() const
    { return 0; }

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

};


} // namespace Go


#endif  // _PTPTINTERSECTOR_H
