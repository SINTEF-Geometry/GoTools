//===========================================================================
//                                                                           
// File: SfSelfIntersector.h
//                                                                           
// Created: 25.05.05
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: $Id: SfSelfIntersector.h,v 1.20 2008-12-10 16:21:05 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _SFSELFINTERSECTOR_H
#define _SFSELFINTERSECTOR_H


#include "GoTools/intersections/Intersector.h"


namespace Go {


class ParamSurfaceInt;
class SurfaceAssembly;


struct SingBox
{
    std::vector<std::pair<double, int> > box_limit_;
    std::shared_ptr<IntersectionPoint> sing_;

    SingBox(std::vector<std::pair<double, int> > box, 
	    std::shared_ptr<IntersectionPoint> sing)
	{
	    box_limit_ = box;
	    sing_ = sing;
	}
};

struct SingUnion
{
    double limit_[4];  // Box limit of union
    std::vector<int> singbox_idx_;

    SingUnion(double box[], std::vector<int> idx)
	{
	    for (int ki=0; ki<4; ki++)
		limit_[ki] = box[ki];
	    singbox_idx_ = idx;
	}

    bool isInside(double *mima)
	{
	    if (limit_[0] < mima[0])
		return false;
	    if (limit_[1] > mima[1])
		return false;
	    if (limit_[2] < mima[2])
		return false;
	    if (limit_[3] > mima[3])
		return false;
	    return true;
	}
};


/// This class finds self-intersections for a parametric surface.

class SfSelfIntersector : public Intersector {
public:

    /// Constructor.
    /// \param surf the parametric surface.
    /// \param epsge the associated tolerance.
    /// \param prev the "parent" Intersector (0 if there is no parent).
    SfSelfIntersector(std::shared_ptr<ParamSurfaceInt> surf,
		      double epsge, Intersector *prev = NULL);

    /// Constructor.
    /// \param surf the parametric surface.
    /// \param epsge the associated tolerance.
    /// \param prev the "parent" Intersector (0 if there is no parent).
    SfSelfIntersector(std::shared_ptr<ParamSurfaceInt> surf,
		      std::shared_ptr<GeoTol> epsge, 
		      Intersector *prev = NULL);

    /// Destructor.
    virtual ~SfSelfIntersector();

    /// Compute topology of self-intersection.
    /// \param compute_at_boundary indicate if we want to compute at
    /// the boundary
    virtual void compute(bool compute_at_boundary = true);

    /// Get the number of parameter directions for the object.
    /// \return the number of parameter directions (i.e. 2)
    virtual int numParams() const
    { return 2; }

    /// Set the maximum number of recursion levels.
    /// \param max_rec the number of levels
    void setMaxRec(int max_rec)
    { max_rec = max_rec_; }

    /// Count the number of complex domains in the object.
    /// \return The number of complex domains in the object.
    int getNmbComplexDomain()
    { return (int)complex_domain_.size(); }

    /// Verifies that this computation is a self intersection problem
    virtual bool isSelfIntersection()
    { return true;  }

protected:

    // Compute topology of self-intersection when the surface is known
    // to be G1
    bool computeG1();

    // Get non self-intersecting subsurfaces corresponding to
    // the current surface
    std::vector<std::shared_ptr<ParamSurfaceInt> > 
    getNonSelfintersecting()
    { return non_selfint_; }

    virtual void print_objs() {}

    virtual int getBoundaryIntersections()
    { return 0; }

    virtual int performInterception()
    { return 0; }

    virtual int simpleCase()
    { return 0; }

    virtual bool isLinear()
    { return false; }

    virtual bool complexityReduced()
    { return false; }

    virtual void handleComplexity() {}

    virtual int checkCoincidence()
    { return 0; }
    
    virtual void microCase() {}
   
    virtual int updateIntersections()
    { return 0; }

    virtual int repairIntersections();

    virtual int linearCase()
    { return 0; }

    virtual int doSubdivide()
    { return 0; }

    virtual void printDebugInfo() {}

    virtual void addComplexDomain(RectDomain dom)
    { complex_domain_.push_back(dom); }

private:

    std::shared_ptr<ParamSurfaceInt> surf_;
//     std::shared_ptr<GeoTol> epsge_;
//     std::shared_ptr<IntersectionPool> int_results_;
//     Intersector *prev_intersector_;

    std::vector<std::shared_ptr<ParamSurfaceInt> > non_selfint_;
    int max_rec_;
    std::vector<RectDomain> complex_domain_;
    std::shared_ptr<SurfaceAssembly> div_sf_;
    std::vector<SingBox> sing_box_;
    std::vector<SingUnion> sing_union_;
   
    // Define parameters in which to subdivide the surface
    bool setSubdivision();

    bool getSingularities(std::vector<std::
			  shared_ptr<IntersectionPoint> >& sing);

    void 
	getPrevDivision(std::vector<std::pair<double,int> >& divpar_u,
			std::vector<std::pair<double,int> >& divpar_v);

    void sortSingularities(std::vector<std::
			   shared_ptr<IntersectionPoint> >& sing,
			   std::vector<std::pair<double,int> >& divpar_u,
			   std::vector<std::pair<double,int> >& divpar_v);

    bool defineSingTopology(std::vector<std::
			   shared_ptr<IntersectionPoint> >& sing,
			   std::vector<std::pair<double,int> >& divpar_u,
			   std::vector<std::pair<double,int> >& divpar_v,
			   std::vector<RectDomain>& sing_domain);

     void splitUnions(std::vector<std::pair<double,int> >& divpar_u,
		     std::vector<std::pair<double,int> >& divpar_v);

    void divideOneUnion(std::vector<std::pair<double,int> >& divpar_u,
			std::vector<std::pair<double,int> >& divpar_v);

    bool isSeparated(std::shared_ptr<IntersectionPoint> sing1,
		     std::shared_ptr<IntersectionPoint> sing2,
		     std::vector<std::pair<double,int> >& divpar_u,
		     std::vector<std::pair<double,int> >& divpar_v);

    std::vector<std::pair<double, std::pair<double, double> > >
	getOverlapBox(double minval, double maxval, int dir);

    double 
	getSplitValue(std::vector<std::pair<double, std::pair<double, double> > >& overlap,
		      double minval, double maxval, double divval);

    void validateSingularities(std::vector<std::
			       shared_ptr<IntersectionPoint> >& sing);

    void  refineSingularityLinks(std::vector<std::
				 shared_ptr<IntersectionPoint> >& sing,
				 std::vector<std::pair<double,int> >& divpar,
				 int dir);

    void makeDivisionLines(std::vector<std::
			   shared_ptr<IntersectionPoint> >& sing,
			   std::vector<double>& sing_par, int dir,
			   double start, double end,
			   std::vector<std::pair<double,int> >& divpar);    

    bool hasBoundaryIntersections(std::shared_ptr<ParamSurfaceInt> sub_sf);

    void computeBoundaryIntersections(std::shared_ptr<ParamSurfaceInt>
				      sub_sf,
				      IntersectionPoint *sing,
				      int is_handled[]);

    void getNmbDivision(int& div1, int& div2);

    void setDivisionValues(std::vector<std::pair<double,int> >& divpar,
			   int& nmbpar, int dir);
	
    void checkDivisionLines(std::vector<std::
			    shared_ptr<IntersectionPoint> >& sing,
			    std::vector<std::pair<double,int> >& divpar_u,
			    std::vector<std::pair<double,int> >& divpar_v);

    void checkOneDivLine(std::vector<std::
			 shared_ptr<IntersectionPoint> >& sing,
			 int dir, double frac,
			 std::vector<std::pair<double,int> >& div1,
			 std::vector<std::pair<double,int> >& div2);

     void getMaxCurvatures(std::shared_ptr<ParamSurfaceInt> surf, int nsample, 
			   std::vector<std::pair<double,double> >& max_curv1,
			   std::vector<std::pair<double,double> >& max_curv2);

    // Check if two sub surfaces share a singularity
    bool shareSingularity(std::shared_ptr<ParamSurfaceInt> sub1, 
			  std::shared_ptr<ParamSurfaceInt> sub2, 
			  double sing[]);

    bool hasSingularity(std::shared_ptr<ParamSurfaceInt> curr_sub);

    void modifyDivisionValues(std::vector<std::shared_ptr<IntersectionPoint> >& sing,
			      std::vector<std::pair<double,int> >& divpar, int dir);

    bool closeKnot(double par, std::vector<double>& knots, size_t& knot_idx,
		   double par_div, double& curr_knot);

    void estimateSingBox(std::shared_ptr<ParamSurfaceInt> normsf,
			 double param[], 
			 std::vector<std::pair<double, int> >& sing_box);

    void makeSingularityUnions();

    double getMinCurvatureRadAlongCurve(std::shared_ptr<ParamSurfaceInt> normsf, 
					int dir, double par, 
					double tmin, double tmax, double param[]);

    double getMinDistAlongCurve(double param[], int dir, double par, 
				double tmin, double tmax);

    int isInComplexDomain(std::shared_ptr<ParamSurfaceInt> subsf);

    void writeDebugComplex(int file, std::vector<RectDomain>& domain);

};


} // namespace Go


#endif  // _SFSELFINTERSECTOR_H


