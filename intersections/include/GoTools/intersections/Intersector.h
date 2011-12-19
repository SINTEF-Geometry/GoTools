//===========================================================================
//                                                                           
// File: Intersector.h
//                                                                           
// Created: 
//                                                                           
// Author: vsk
//                                                                           
// Revision: $Id: Intersector.h,v 1.44 2007-06-22 16:33:21 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _INTERSECTOR_H
#define _INTERSECTOR_H


#include "GoTools/intersections/SubdivisionClassification.h"
#include "GoTools/intersections/SingularityClassification.h"
#include "GoTools/intersections/SingularityInfo.h"
#include "GoTools/intersections/ComplexityInfo.h"
#include "GoTools/utils/Point.h"
#include "GoTools/geometry/RectDomain.h"
#include <vector>
#include <memory>


namespace Go {


class IntersectionPoint;
class IntersectionCurve;
class IntersectionPool;
class GeoTol;
struct BoundaryGeomInt;


/// This class is an abstract class providing an interface to the
/// intersection functionality in GoTools.

class Intersector {
public:

    /// Default constructor
    Intersector() : prev_intersector_(0) {}

    /// Constructor.
    /// \param epsge the geometric tolerance for the intersector.
    /// \param prev the previous intersector.
    Intersector(double epsge, Intersector *prev = 0);

    /// Constructor.
    /// \param epsge the geometric tolerance for the intersector.
    /// \param prev the previous intersector.
    Intersector(shared_ptr<GeoTol> epsge, Intersector *prev = 0);

    /// Destructor
    virtual ~Intersector(){};

    /// Compute the current intersections (topology).
    /// \param compute_at_boundary if true we will include computation
    /// of boundary intersections.
    virtual void compute(bool compute_at_boundary=true);

    /// Get intersection points and curves.  Sends request to
    /// IntersectionPool. The intersection points are isolated.
    /// \param int_points vector of intersection points
    /// \param int_curves vector of intersection curves
    virtual void
    getResult(std::vector<shared_ptr<IntersectionPoint> >& int_points,
	      std::vector<shared_ptr<IntersectionCurve> >& int_curves);

//     // Validation of given intersection results
//     virtual void validate(int level, ValidationStat status);

    /// Get the IntersectionPool for the intersector.
    /// \return The IntersectionPool.
    shared_ptr<IntersectionPool> getIntPool()
    { return int_results_; }

    /// Validate this pool and its siblings.
    /// \return \a true if all the sibling pools are valid, \a false
    /// otherwise
    bool validateSiblingPools();


    /// Get the tolerance object used by this Intersector
    /// \return The tolerance object for the intersector.
    shared_ptr<GeoTol> getTolerance()
    { return epsge_;}

    /// Verify whether singularities has been set.
    /// \return True if info on singularities has been set.
    bool hasSingularityInfo()
    { return (singularity_info_.get() != 0); } 

    /// Get info regarding singularities.
    /// \return The singularify info for the intersector.
    shared_ptr<SingularityInfo> getSingularityInfo()
    { return singularity_info_; } 

    /// Set info regarding singularities.
    /// \param previous the singularity info.
    /// \param missing_dir if the dimension of the intersection
    /// problem has been reduced from previous, the index tells us
    /// which parameter direction that was removed. Indexing starts at
    /// 0.
    void setSingularityInfo(shared_ptr<SingularityInfo> previous, 
			    int missing_dir)
    {
	if (missing_dir < 0)
	    singularity_info_ = (shared_ptr<SingularityInfo>)
		(new SingularityInfo(previous));
	else
	    singularity_info_ = (shared_ptr<SingularityInfo>)
		(new SingularityInfo(previous, missing_dir));
    }

    /// Instruct the intersector about known singular points.
    /// \param par the parameter value of the singularity. Size of
    /// array should be numParams().
    void setHighPriSing(double* par);

    /// Verify whether there has been created info regarding the
    /// complexity of the problem.
    /// \return True if complexity info has been created.
    bool hasComplexityInfo()
    { return (complexity_info_.get() != 0); } 

    /// Get info on the complexity of the problem.
    /// \return Info on the complexity of the problem.
    shared_ptr<ComplexityInfo> getComplexityInfo()
    { return complexity_info_; } 

    /// Get the number of parameter directions for the intersection.
    /// \return The number of parameter directions for the
    /// intersection.
    virtual int numParams() const = 0;

    /// Count the number of boundary objects belonging to the
    /// specified ParamGeomInt.
    /// \param idx refers to obj1 or obj2. Indexing starts at 0.
    /// \return The number of boundary objects.
    virtual int nmbBdObj(int idx) const
    { return 0; }

    /// Get the specified boundary object belonging to the specified
    /// ParamGeomInt.
    /// \param idx refers to obj1 or obj2. Indexing starts at 0.
    /// \param bd_idx index of the boundary object in the
    /// ParamGeomInt.
    /// \return The boundary object.
    virtual BoundaryGeomInt* getBoundaryObject(int idx, int bd_idx) const
    { return 0; }

    /// The current recursion level.
    /// \return The current recursion level.
    int nmbRecursions()
    {
	if (prev_intersector_ == 0)
	    return 0;
	else 
	    return prev_intersector_->nmbRecursions() + 1;
    }

    /// Verify whether the surface is self-intersecting.
    /// \return True if the surface is self-intersecting.
    virtual bool isSelfIntersection()
    { return false; }  // Default

    /// Check if this intersection algorithm is performed in a
    /// self-intersection context.
    virtual int  isSelfintCase()
    { return 0; }  // Default behaviour

    virtual void addComplexDomain(RectDomain dom)
    { ; }

    /// Write diagnostic information about the intersection points
    void writeIntersectionPoints() const;

    friend class SfSfIntersector;
    friend class IntersectionPool;

protected:
    // Data members

    // @ Logical problem here?  An IntersectionPool refers to 2 objects,
    // @ while an Intersector does not (need Intersector2Obj to do
    // @ that...)
    shared_ptr<IntersectionPool> int_results_;
    std::vector<shared_ptr<Intersector> > sub_intersectors_;
    Intersector *prev_intersector_;
    shared_ptr<GeoTol> epsge_;
    shared_ptr<SingularityInfo> singularity_info_;
    shared_ptr<ComplexityInfo> complexity_info_;

    //     virtual shared_ptr<Intersector> 
    //       lowerOrderIntersector(shared_ptr<ParamObjectInt> obj1,
    // 			    shared_ptr<ParamObjectInt> obj2, 
    // 			    Intersector* prev = 0,
    // 			    int eliminated_parameter = -1,
    // 			    double eliminated_value = 0) = 0;

    virtual void print_objs() = 0;

    virtual int getBoundaryIntersections() = 0;

    virtual int performInterception() = 0;

    virtual int simpleCase() = 0;

    virtual bool isLinear() = 0;

    virtual bool degTriangleSimple()
    { 
	// Default implementation that is OK for most sub classes
	return false;   
    }

    virtual bool complexityReduced() = 0;

    virtual void handleComplexity() = 0;

    virtual int checkCoincidence() = 0;

    virtual void microCase() = 0;

    virtual int updateIntersections() = 0;

    virtual int repairIntersections() = 0;

    virtual int linearCase() = 0;

    virtual int doSubdivide() = 0;

    virtual int complexIntercept()
	{
	    return 0;  // Overridden when required
	}

    virtual int complexSimpleCase()
	{
	    return 0;  // Overridden when required
	}

    virtual void doPostIterate()
	{
	    ;  // Overridden when required
	}

    virtual void printDebugInfo() = 0;
private:

};


} // namespace Go


#endif  // _INTERSECTOR_H
