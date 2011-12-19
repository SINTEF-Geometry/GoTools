//===========================================================================
//                                                                           
// File: SfSfIntersector.h 
//                                                                           
// Created: 
//                                                                           
// Author: oan
//                                                                           
// Revision: $Id: SfSfIntersector.h,v 1.50 2007-11-01 14:31:37 vsk Exp $
//                                                                           
// Description:
//
//===========================================================================

#ifndef _SFSFINTERSECTOR_H
#define _SFSFINTERSECTOR_H


#include "GoTools/intersections/Intersector2Obj.h"
#include "GoTools/intersections/IntersectionPoint.h"
#include "GoTools/intersections/IntersectionPool.h"


namespace Go {


class Param2FunctionInt;


/// This class performs intersection between two parametric surfaces.

class SfSfIntersector : public Intersector2Obj {
public:

    friend class SfSelfIntersector;

    /// Default constructor
    SfSfIntersector() {}

    /// Constructor.
    /// Both objects should refer to surfaces (this is not checked
    /// compile-time, so we rely on the user to obey this rule).
    /// \param obj1 of type ParamSurfaceInt.
    /// \param obj2 of type ParamSurfaceInt.
    /// \param epsge the associated tolerance.
    /// \param prev the "parent" Intersector (0 if there is no
    /// parent).
    SfSfIntersector(shared_ptr<ParamGeomInt> obj1, 
		    shared_ptr<ParamGeomInt> obj2,
		    shared_ptr<GeoTol> epsge, 
		    Intersector* prev = 0);

    /// Constructor.
    /// Both objects should refer to surfaces (this is not checked
    /// compile-time, so we rely on the user to obey this rule).
    /// \param obj1 of type ParamSurfaceInt.
    /// \param obj2 of type ParamSurfaceInt.
    /// \param epsge the associated tolerance.
    /// \param prev the "parent" Intersector (0 if there is no
    /// parent).
   SfSfIntersector(shared_ptr<ParamGeomInt> obj1, 
		    shared_ptr<ParamGeomInt> obj2,
		    double epsge, 
		    Intersector* prev = 0);
    /// Destructor
    virtual ~SfSfIntersector();

    /// Return the number of parameter directions for the object.
    /// \return the number of parameter directions
    virtual int numParams() const
    { return 4; }

    void postIterate3(int nmb_orig, int dir)
	{
	    postIterate(nmb_orig, dir, false);
	}

    friend class IntersectionPool;

protected:

    virtual shared_ptr<Intersector> 
      lowerOrderIntersector(shared_ptr<ParamGeomInt> obj1,
			    shared_ptr<ParamGeomInt> obj2, 
			    Intersector* prev = 0,
			    int eliminated_parameter = -1,
			    double eliminated_value = 0);

    virtual int performRotatedBoxTest(double eps1, double eps2);

    virtual bool foundIntersectionNearBoundary();

    virtual int performInterceptionByImplicitization();

    virtual int interceptionBySeparationSurface();

    virtual int simpleCase2(Point& axis1, Point& axis2);

    virtual int simpleCaseByImplicitization();

    virtual bool complexityReduced();

    virtual void handleComplexity();

    virtual int checkCoincidence();
    
    virtual void microCase();
    
    virtual bool degTriangleSimple();

    virtual int updateIntersections();

    virtual int repairIntersections();

    void repairSingularityBox();

    void repairFalseBranchPoints();

    void removeIsolatedPoints();

    void repairMissingLinks();

    bool connectIfPossible(shared_ptr<IntersectionPoint> pt1,
			   shared_ptr<IntersectionPoint> pt2);

    /// Fix crossing intersection links.
    void fixCrossingLinks();

    bool checkCloseEndpoint(shared_ptr<IntersectionPoint> pnt, 
			    shared_ptr<IntersectionLink> link);

    void iterateOnIntersectionPoints();

    void getSingularityBox(shared_ptr<IntersectionPoint> sing,
			   double frompar[], double topar[]);

    shared_ptr<SfSfIntersector>
    getSubIntersector(double frompar[], double topar[]);

    virtual int linearCase();

    virtual int doSubdivide();

    void getApproxImplicit(std::vector<std::vector<
			   shared_ptr<Param2FunctionInt> > >& approx_implicit,
			   std::vector<double>& approx_implicit_err,
			   std::vector<double>& approx_implicit_gradsize,
			   std::vector<double>& approx_implicit_gradvar);
    
private:

    // Implicitization related data members

    // This object contains the result of plugging in the two surfaces
    // into the two approximate implicitizations. The convention is
    // that
    //
    // approx_implicit_[implicit][surface] ,
    //
    // where implicit = 0 or 1, and surface = 0 or 1, means the result
    // of plugging surface 'surface' into implicit function
    // 'implicit'.
    std::vector<std::vector<shared_ptr<Param2FunctionInt> > >
    approx_implicit_;

    // The "error" of the implicitizations. Equals the smallest
    // singular values of the implicitization "D-matrix", and bounds
    // the algebraic distance between the exact object and the
    // implicit approximation. The index (0 or 1) refers to the
    // implicit surface.
    std::vector<double> approx_implicit_err_;

    // The size of the gradient of the implicit surface over the
    // inserted surface. The measure of the size is taken to be the
    // l2-norm of the control points of approx_implicit_[0][1] and
    // approx_implicit_[1][0]. The index (0 or 1) refers to the
    // inserted parametric surface.
    std::vector<double> approx_implicit_gradsize_;

    // The variation of the gradient of the implicit surface over the
    // inserted surface. The measure of the variation is taken to be
    // max length divided by min length of the control vectors of
    // approx_implicit_[0][1] and approx_implicit_[1][0]. The index (0
    // or 1) refers to the inserted parametric surface.
    std::vector<double> approx_implicit_gradvar_;

    // Indicates if the results is fetched from this or the previous
    // level of the recursion
    bool prev_implicit_[2];

    // VSK, 0806. Information about monotone trend in simple case
    int sorting_obj_[2];  // -1 if not set, otherwise 0 or 1
    Point sorting_dir_[2];  // Direction along which to sort intersection
                         // points in the parameter domain of object sorting_obj_
    bool sorting_pardir_[2];  // Indicates whether the parameter directions in
                              // the intersection points are consistent

    // Service funtionality for subdivision
    int sortParameterDirections(int perm[], int deg_edge[]);

    SubdivisionClassification 
    getSubdivisionParameter(int dir, int deg_edge, double& par);

    bool getSubdivParDegEdge(ParamSurfaceInt *surf, int dir, int pdir, 
			     int deg_edge, double treshhold, double& par);

    void setDegTriangle(std::vector<shared_ptr<ParamGeomInt> >&
			sub_objects,
			int deg_edge, int pdir, double par);

    int checkSubdivParam(int dir, double par, double ta, double tb,
			 std::vector<shared_ptr<IntersectionPoint> >&
			 int_pts);

    int checkIsoCurve(int pdir, bool first, double par,
		      std::vector<shared_ptr<IntersectionPoint> >&
		      int_pts);

    bool getSubdivAtSing(int dir, double ta, double tb, double& par);
    
    int splitIntResults(std::vector<shared_ptr<IntersectionPoint> >&
			 int_pts, int pardir, double par, double start, 
			 double end, bool large_move);

    void doIterate(int pardir, double parval, double param[], double& dist,
		   double seed[]);

    void doIterate(int pardir, double param[], double limit1[], double w1,
		   double limit2[], double w2, double& dist, double *seed = 0);

    virtual void postIterate(int nmb_orig, int dir, bool keep_endpt=true);

    void postIterate2(int nmb_orig, int dir, double par,
		      bool along = true, bool only_isolated = true);

    virtual void removeDegenerateConnections();

    // Service functionality for connecting intersection points at the 
    // boundaries in a simple case situation
    bool isConnected(std::vector<shared_ptr<IntersectionPoint> >&
		     bd_ints);

    bool isConnected(std::vector<std::
		     pair<shared_ptr<IntersectionPoint>,
		     IntPtClassification> >& bd_ints, 
		     int nmb_nottouch);

    bool connectDirected(std::vector<std::
			 pair<shared_ptr<IntersectionPoint>,
			 IntPtClassification> >& bd_ints, 
			 int nmb_nottouch);

    bool canConnect(shared_ptr<IntersectionPoint> pt1,
		    shared_ptr<IntersectionPoint> pt2);

    bool findMiddlePoint(shared_ptr<IntersectionPoint> pt1,
			 shared_ptr<IntersectionPoint> pt2,
			 double param[], double& dist);

    // Service funtionality for handling tiny, degenerate triangles
    void getPointsAtDegEdges(std::vector<shared_ptr<IntersectionPoint> >& result);

    void getPointsOppositeDegEdges(std::vector<std::pair<
                                   shared_ptr<IntersectionPoint>, int> >&
				   result);

    void getLinksAtDegEdges(std::vector<shared_ptr<IntersectionPoint> >& deg_pnts, 
			    int idx, shared_ptr<IntersectionPoint> curr,
			    std::vector<shared_ptr<IntersectionLink> >& links);

    // Service funtionality and structure for handling of sub problems
    // that are recognize as complex and not recommended for further
    // subdivision

    struct IntersectionGroup {
	std::vector<shared_ptr<IntersectionPoint> > points;
	std::vector<IntPtClassification> type;
	IntPtClassification main_type;
	bool tangential;
	bool iso[4];
	
	IntersectionGroup(std::vector<shared_ptr<IntersectionPoint> >& pts)
	{
	    points = pts;
	    type.resize(pts.size());
	    for (size_t ki=0; ki<type.size(); ki++)
		type[ki] = DIR_UNDEF;
	    main_type = DIR_UNDEF;
	    tangential = false;
	    for (int dir=0; dir<4; dir++)
		iso[dir] = false;
	}

	void classifyPoints(ParamSurface *srf1, ParamSurface *srf2);
	void classifyCurve();
	shared_ptr<IntersectionPoint> 
	closestToPoint(IntersectionPoint* pnt);
	shared_ptr<IntersectionPoint> getBranchPoint();
	shared_ptr<IntersectionPoint> 
	getMainPoint(IntPtClassification &type_curr, 
		     IntPtClassification type_other);
	bool isIso()
	{
	    for (int ki=0; ki<4; ki++)
		if (iso[ki])
		    return true;
	    return false;
	}
    };
    
    void makeIntersectionGroups(std::vector<
				shared_ptr<IntersectionPoint> >& pts,
				std::vector<
				shared_ptr<IntersectionGroup> >& groups);

    void classifyGroups(std::vector<
			shared_ptr<IntersectionGroup> >& groups);

    void tryConnectGroups(shared_ptr<IntersectionGroup> group1,
			  shared_ptr<IntersectionGroup> group2);

    void mergeGroups(std::vector<shared_ptr<IntersectionGroup> >& group,
		     IntersectionPoint *sing);

    void connectToSing(shared_ptr<IntersectionGroup> group,
		       IntersectionPoint *sing);

    void selfintComplex(IntersectionPoint *sing);
 
    void selfintComplex2();
 
    // Perform an approximate implicitization on both surfaces, then
    // plug both surfaces into the resulting algebraic
    // expressions. This results in (up to) four combinations of
    // plugged in objects, as well as the two implicitization errors.
    void createApproxImplicit(double tol, int recursion_limit = 10);

    // If object has a prev intersector (i.e. a parent), we initialize
    // appr_impl_tol_ based on that.
    void setApproxImplicitFromPrev();

    // Set the values of approx_implicit_gradsize_[i] and
    // approx_implicit_gradvar_[i], where 'i' is the index (0 or 1) of
    // the inserted parametric surface.
    void setGradValues(int i);

    // DEBUG
    void writeDebugConnect(std::vector<std::
			   pair<shared_ptr<IntersectionPoint>,
			   IntPtClassification> >& bd_ints);

    void writeDebugLinear(std::vector<
			  shared_ptr<IntersectionPoint> >&  bd_ints);

};


} // namespace Go


#endif  // _SFSFINTERSECTOR_H
