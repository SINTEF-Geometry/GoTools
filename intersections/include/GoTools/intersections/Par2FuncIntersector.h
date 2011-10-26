//===========================================================================
//                                                                           
// File: Par2FuncIntersector.h                                               
//                                                                           
// Created: Mon Sep 27 14:26:50 2004                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: Par2FuncIntersector.h,v 1.16 2006-11-03 14:43:14 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _PAR2FUNCINTERSECTOR_H
#define _PAR2FUNCINTERSECTOR_H


#include "GoTools/intersections/IntersectorFuncConst.h"
#include "GoTools/intersections/IntersectionPoint.h"


namespace Go {


/// This class is performing intersections between a 1-dimensional
/// parametric surface and a constant.

class Par2FuncIntersector : public IntersectorFuncConst {
public:

//     Par2FuncIntersector(std::shared_ptr<ParamFunctionInt> func,
// 			std::shared_ptr<ParamFunctionInt> C,
// 			double epsge,
// 			Intersector* prev = 0);

    /// Constructor.
    /// One of the objects should refer to a 1D-surface, the other a
    /// constant (this is not checked compile-time, so we rely on the
    /// user to obey this rule).  The last two variables are relevant
    /// only if the parent has one more parameter than the Intersector
    /// to be constructed.
    /// \param func of type Param2FunctionInt.
    /// \param C of type Param0FunctionInt.
    /// \param epsge the associated tolerance.
    /// \param prev the "parent" Intersector (0 if there is no
    /// parent).
    /// \param eliminated_parameter the index of the parameter that
    /// was removed from the parent \a prev.
    /// \param eliminated_value the value of the parameter that was
    /// removed from the parent \a prev.
    Par2FuncIntersector(std::shared_ptr<ParamFunctionInt> func,
			std::shared_ptr<ParamFunctionInt> C,
			std::shared_ptr<GeoTol> epsge,
			Intersector *prev = 0,
			int eliminated_parameter = -1,
			double eliminated_value = 0);

    /// Destructor.
    virtual ~Par2FuncIntersector();

//     // Validation of given intersection results
//     virtual void validate(int level, ValidationStat status);

    /// Return the number of parameter directions for the object.
    /// \return the number of parameter directions
    virtual int numParams() const
    { return 2; }
	
protected:
    // Data members

    virtual std::shared_ptr<Intersector> 
    lowerOrderIntersector(std::shared_ptr<ParamFunctionInt> obj1,
			  std::shared_ptr<ParamFunctionInt> obj2,
			  Intersector* prev = 0,
			  int eliminated_parameter = -1,
			  double eliminated_value = 0);

    virtual int checkCoincidence();

    virtual void microCase();
    
    virtual int updateIntersections();

    virtual int repairIntersections()
    { return 0; }

    bool
    isConnected(std::vector<std::shared_ptr<IntersectionPoint> > bd_ints,
		int nmbbd);

    bool isConnected(std::vector<std::
		     pair<std::shared_ptr<IntersectionPoint>,
		     IntPtClassification> >& bd_ints, 
		     int nmb_nottouch);

    bool connectDirected(std::vector<std::
			 pair<std::shared_ptr<IntersectionPoint>,
			 IntPtClassification> >& bd_ints,
			 int nmbbd);

    bool canConnect(std::shared_ptr<IntersectionPoint> pt1,
		    std::shared_ptr<IntersectionPoint> pt2);

    virtual int doSubdivide();

private:

    int getSubdivisionParameter(int dir, double& par);

    // We need to decide in which direction to subdivide.
    int sortParameterDirections(int perm[]); //, int deg_edge[]);

    int checkSubdivParam(int dir, double par, double ta, double tb,
			 std::vector<std::
			 shared_ptr<IntersectionPoint> >& int_pts);

    int checkIsoCurve(int pdir, bool first, double par,
 		      std::vector<std::
		      shared_ptr<IntersectionPoint> > int_pts);

    bool getSubdivAtSing(int dir, double ta, double tb, double& par);

//     void splitIntResults(std::vector<std::
// 			 shared_ptr<IntersectionPoint> >& int_pts,
// 			 int pardir, double par,
// 			 double start, double end);

//     void doIterate(int pardir, double parval, double param[], double& dist,
// 		   double seed[]);

//     // Utility function for sorting input bd_int's.
//     IntPtClassification bdDir(const IntersectionPoint& int_pt,
// 			      Point sorting_dir);

    void writeDebugConnect(std::vector<std::
			   pair<std::shared_ptr<IntersectionPoint>,
			   IntPtClassification> >& bd_ints);

};


} // namespace Go


#endif // _PAR2FUNCINTERSECTOR_H

