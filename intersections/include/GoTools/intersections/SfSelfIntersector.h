/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

#ifndef _SFSELFINTERSECTOR_H
#define _SFSELFINTERSECTOR_H


#include "GoTools/intersections/Intersector.h"


namespace Go {


class ParamSurfaceInt;
class SurfaceAssembly;


struct SingBox
{
    std::vector<std::pair<double, int> > box_limit_;
    shared_ptr<IntersectionPoint> sing_;

    SingBox(std::vector<std::pair<double, int> > box, 
	    shared_ptr<IntersectionPoint> sing)
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
    SfSelfIntersector(shared_ptr<ParamSurfaceInt> surf,
		      double epsge, Intersector *prev = NULL);

    /// Constructor.
    /// \param surf the parametric surface.
    /// \param epsge the associated tolerance.
    /// \param prev the "parent" Intersector (0 if there is no parent).
    SfSelfIntersector(shared_ptr<ParamSurfaceInt> surf,
		      shared_ptr<GeoTol> epsge, 
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
    std::vector<shared_ptr<ParamSurfaceInt> > 
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

    shared_ptr<ParamSurfaceInt> surf_;
//     shared_ptr<GeoTol> epsge_;
//     shared_ptr<IntersectionPool> int_results_;
//     Intersector *prev_intersector_;

    std::vector<shared_ptr<ParamSurfaceInt> > non_selfint_;
    int max_rec_;
    std::vector<RectDomain> complex_domain_;
    shared_ptr<SurfaceAssembly> div_sf_;
    std::vector<SingBox> sing_box_;
    std::vector<SingUnion> sing_union_;
   
    // Define parameters in which to subdivide the surface
    bool setSubdivision();

    bool getSingularities(std::vector<
			  shared_ptr<IntersectionPoint> >& sing);

    void 
	getPrevDivision(std::vector<std::pair<double,int> >& divpar_u,
			std::vector<std::pair<double,int> >& divpar_v);

    void sortSingularities(std::vector<
			   shared_ptr<IntersectionPoint> >& sing,
			   std::vector<std::pair<double,int> >& divpar_u,
			   std::vector<std::pair<double,int> >& divpar_v);

    bool defineSingTopology(std::vector<
			   shared_ptr<IntersectionPoint> >& sing,
			   std::vector<std::pair<double,int> >& divpar_u,
			   std::vector<std::pair<double,int> >& divpar_v,
			   std::vector<RectDomain>& sing_domain);

     void splitUnions(std::vector<std::pair<double,int> >& divpar_u,
		     std::vector<std::pair<double,int> >& divpar_v);

    void divideOneUnion(std::vector<std::pair<double,int> >& divpar_u,
			std::vector<std::pair<double,int> >& divpar_v);

    bool isSeparated(shared_ptr<IntersectionPoint> sing1,
		     shared_ptr<IntersectionPoint> sing2,
		     std::vector<std::pair<double,int> >& divpar_u,
		     std::vector<std::pair<double,int> >& divpar_v);

    std::vector<std::pair<double, std::pair<double, double> > >
	getOverlapBox(double minval, double maxval, int dir);

    double 
	getSplitValue(std::vector<std::pair<double, std::pair<double, double> > >& overlap,
		      double minval, double maxval, double divval);

    void validateSingularities(std::vector<
			       shared_ptr<IntersectionPoint> >& sing);

    void  refineSingularityLinks(std::vector<
				 shared_ptr<IntersectionPoint> >& sing,
				 std::vector<std::pair<double,int> >& divpar,
				 int dir);

    void makeDivisionLines(std::vector<
			   shared_ptr<IntersectionPoint> >& sing,
			   std::vector<double>& sing_par, int dir,
			   double start, double end,
			   std::vector<std::pair<double,int> >& divpar);    

    bool hasBoundaryIntersections(shared_ptr<ParamSurfaceInt> sub_sf);

    void computeBoundaryIntersections(shared_ptr<ParamSurfaceInt>
				      sub_sf,
				      IntersectionPoint *sing,
				      int is_handled[]);

    void getNmbDivision(int& div1, int& div2);

    void setDivisionValues(std::vector<std::pair<double,int> >& divpar,
			   int& nmbpar, int dir);
	
    void checkDivisionLines(std::vector<
			    shared_ptr<IntersectionPoint> >& sing,
			    std::vector<std::pair<double,int> >& divpar_u,
			    std::vector<std::pair<double,int> >& divpar_v);

    void checkOneDivLine(std::vector<
			 shared_ptr<IntersectionPoint> >& sing,
			 int dir, double frac,
			 std::vector<std::pair<double,int> >& div1,
			 std::vector<std::pair<double,int> >& div2);

     void getMaxCurvatures(shared_ptr<ParamSurfaceInt> surf, int nsample, 
			   std::vector<std::pair<double,double> >& max_curv1,
			   std::vector<std::pair<double,double> >& max_curv2);

    // Check if two sub surfaces share a singularity
    bool shareSingularity(shared_ptr<ParamSurfaceInt> sub1, 
			  shared_ptr<ParamSurfaceInt> sub2, 
			  double sing[]);

    bool hasSingularity(shared_ptr<ParamSurfaceInt> curr_sub);

    void modifyDivisionValues(std::vector<shared_ptr<IntersectionPoint> >& sing,
			      std::vector<std::pair<double,int> >& divpar, int dir);

    bool closeKnot(double par, std::vector<double>& knots, size_t& knot_idx,
		   double par_div, double& curr_knot);

    void estimateSingBox(shared_ptr<ParamSurfaceInt> normsf,
			 double param[], 
			 std::vector<std::pair<double, int> >& sing_box);

    void makeSingularityUnions();

    double getMinCurvatureRadAlongCurve(shared_ptr<ParamSurfaceInt> normsf, 
					int dir, double par, 
					double tmin, double tmax, double param[]);

    double getMinDistAlongCurve(double param[], int dir, double par, 
				double tmin, double tmax);

    int isInComplexDomain(shared_ptr<ParamSurfaceInt> subsf);

    void writeDebugComplex(int file, std::vector<RectDomain>& domain);

};


} // namespace Go


#endif  // _SFSELFINTERSECTOR_H


