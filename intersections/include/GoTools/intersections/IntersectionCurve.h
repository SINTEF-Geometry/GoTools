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

#ifndef _INTERSECTIONCURVE_H
#define _INTERSECTIONCURVE_H


#include "GoTools/intersections/IntersectionPoint.h"
#include "GoTools/intersections/IntersectionLink.h"
#include "GoTools/geometry/ParamCurve.h"
#include <memory>
#include <list>
#include <vector>


namespace Go {


enum EvalKind {SPACECURVE, PARAMCURVE_1, PARAMCURVE_2};
enum TangentDomain {GEOM, PARAM1, PARAM2};
enum EstimateDirection {FORWARDS, BACKWARDS};
    
/// Error object used internally in IntersectionCurve
class Zero_Parameter_Span_Error {}; 
class ParamSurfaceInt;
class IntersectionCurve;


// NB: Before you start to use the IntersectionCurve class, please
// read the following.
//
// IntersectionCurve has now been split up in three curve-types:
//
// * InterpolatedIntersectionCurve - defined by Hermite interpolation
// of IntersectionPoints
//
// * DegeneratedIntersectionCurve - represents an IntersectionCurve
// that is collapsed into a single point.
//
// * IsoparametricIntersectionCurve - An IntersectionCurve that can be
// described by an isoparametric curve in one of the underlying
// objects.
// 
// The distinction between these curve types makes the implementation
// of member functions much clearer.  However, it should not be
// necessary for the user to know whether a given IntersectionCurve is
// actually interpolated, degenerated or isoparametric.  Therefore,
// the construction of IntersectionCurves is left to a special
// function, 'constructIntersectionCurve(...)', into which the user
// only has to pass the range of IntersectionPoints making up the
// curve.  This function is the only way to construct an
// IntersectionCurve, as the actual constructors are protected.
// Therefore, please use the 'constructIntersectionCurve(...) whenever
// you want to make an IntersectionCurve.

/// Object describing the curve that results from an intersection of
/// two geometrical objects.

//===========================================================================
class IntersectionCurve {
//===========================================================================
public:
    /// Virtual destructor
    virtual ~IntersectionCurve() {};

    /// Get out a (shared) pointer to a parametric curve that
    /// approximates this IntersectionCurve.
    virtual shared_ptr<ParamCurve>
    getCurve() const = 0;

    /// Get out a (shared) pointer to the curve in the parametric
    /// plane of the first or second object.  Should only be called
    /// when the concerned object is 2-parametric.
    /// \param obj_nmb should be either 1 or 2, depending on whether
    /// you want to get the parametric curve in object 1 or object 2.
    /// \return a shared pointer to a curve that approximates the
    /// intersection curve in the parametric domain of the specified
    /// object.
    virtual shared_ptr<ParamCurve>
    getParamCurve(int obj_nmb) const = 0;

    /// Get the start and end value for the parametric span of the
    /// IntersectionCurve.
    /// \retval start upon function return this variable will contain
    /// the start value of the parametric span.
    /// \retval end upon function return this variable will contain
    /// the end value of the parametric span.
    virtual void getParamSpan(double& start, double& end) const = 0;

    /// Evaluate the IntersectionCurve at parameter value \a pval.
    /// Its position at this parameter value will be returned in \a
    /// pos and its tangent direction will be returned in \a tan.
    /// \param pval the parameter value for which to evaluate the
    /// IntersectionCurve.  It should be within the parametric span
    /// (which can be obtained from GetParamSpan().
    /// \retval pos the calculated position of the point on the curve
    /// at parameter \a pval.
    /// \param tan the calculated tangent of the curve at parameter
    /// \a pval.
    virtual void evaluateAt(double pval, Point& pos, Point& tan) = 0;
    
    /// Refine the curve (which means to increase the number of
    /// IntersectionPoint s defining it), so that the curve matches a
    /// specific positional and angle tolerance.
    /// \param pos_tol the positional tolerance that the curve should
    /// match after refinement.
    /// \param angle_tol the angular tolerance that the tangent of the
    /// curve should match after refinement.
    virtual void refine(const double& pos_tol, const double& angle_tol) = 0;
    
    /// Determine if the IntersectionCurve is in fact representing an
    /// isoparametric intersection.
    /// \return 'true' if the IntersectionCurve is part of an isocurve
    /// for one of the objects, 'false' otherwise.
    virtual bool isIsocurve() const = 0;

    /// Determine if the IntersectionCurve is in fact a degenerated
    /// curve (a single point in space).
    /// \return 'true' if the IntersectionCurve is
    /// degenerated. 'false' otherwise.
    virtual bool isDegenerated() const = 0;

    /// Get one of the specific guidepoints (IntersectionPoint) that
    /// define the curve.
    /// \param index the index of the requested IntersectionPoint(from
    /// 0 and up to the number of defining IntersectionPoint s - 1).
    /// To find out how many IntersectionPoint s participate in the
    /// definition of the curve, use the numGuidePoints()
    /// function.
    /// \return a shared pointer to the requested IntersectionPoint.
    shared_ptr<IntersectionPoint> getGuidePoint(int index) const;

    /// Get the tangent of a guidepoint that we know lies on the curve
    /// (just using the getTangent() function of the guidepoint can be
    /// deceiving, as it is not necessarily unique).
    /// \param pt shared pointer to the guidepoint.  It is supposed to
    /// be part of the definition of the curve.
    /// \param type specifies what kind of tangent the user wants.
    /// '1' means the parameter plane tangent for the first
    /// object. '2' means the parameter plane tangent for the second
    /// object.  Other values of 'type' means the tangent in 3D space.
    /// \retval tan the requested tangent
    /// \return 'true' if a tangent was found (i.e. \a pt was indeed
    /// a guidepoint of the IntersectionCurve).  In the opposite case,
    /// 'false' is returned.
    virtual bool 
    getGuidePointTangent(shared_ptr<IntersectionPoint> pt,
			 Point& tan, int type = 0) const;
    

    /// Get number of guide points defining the curve.
    /// \return the number of guide points currently defining the
    /// curve (may increase if the 'refine()' function is run).
    int numGuidePoints() const { return (int)ipoints_.size();}
    
    /// Write the IntersectionPoints defining this IntersectionCurve
    /// to stream
    /// \param os output stream
    void  writeIPointsToStream(std::ostream& os) const;

protected:
    
    template<class iterator>
    IntersectionCurve(iterator  begin, iterator end)
	: ipoints_(begin, end) {}

    // The list of intersection points defining the curve
    std::list<shared_ptr<IntersectionPoint> > ipoints_;

    // This is the function that should be used when creating an
    // IntersectionCurve.  It will take care of the decision of which
    // kind (concrete class) of curve that is appropriate.
    template<class iterator> friend 
    shared_ptr<IntersectionCurve>
    constructIntersectionCurve(const iterator begin,
			       const iterator end);


};


/// IntersectionCurve that is degenerated into a single point.

//===========================================================================
class DegeneratedIntersectionCurve : public IntersectionCurve {
//===========================================================================
public:
    virtual ~DegeneratedIntersectionCurve();

    virtual shared_ptr<ParamCurve> getCurve() const;

    virtual shared_ptr<ParamCurve> getParamCurve(int obj_nmb) const;

    virtual bool isIsocurve() const
    { return false; } // Maybe this should be investigated

    virtual bool isDegenerated() const
    { return true; }

    virtual void getParamSpan(double& start, double& end) const
    {
	start = 0;
	end = 1;
    }

    virtual void refine(const double& pos_tol, const double& angle_tol)
    {
	MESSAGE("Tried to refine a degenerate curve.  Ignoring...");
    }

    virtual void evaluateAt(double pval, Point& pos, Point& tan) 
    {
	pos = ipoints_.front()->getPoint();
	tan = ipoints_.front()->getTangent();
    }

private:

    template<class iterator>
    DegeneratedIntersectionCurve(const iterator begin, const iterator end) 
	: IntersectionCurve(begin, end) 
    {
	MESSAGE("Created a degenerate IntersectionCurve object.");
    }
    
    template<class iterator>
    static bool degenerated_range(const iterator begin, const iterator end);

    template<class iterator> friend 
    shared_ptr<IntersectionCurve>
    constructIntersectionCurve(const iterator begin,
			       const iterator end);
};


/// IntersectionCurve that cannot be evaluated

//===========================================================================
class NonEvaluableIntersectionCurve : public IntersectionCurve
//===========================================================================
{
public:
    virtual ~NonEvaluableIntersectionCurve();

    virtual shared_ptr<ParamCurve> getCurve() const
    { return shared_ptr<ParamCurve>(); }

    virtual shared_ptr<ParamCurve> getParamCurve(int obj_nmb) const
    { return shared_ptr<ParamCurve>(); }

    virtual void getParamSpan(double& start, double& end) const
    { start = end = 0; }

    virtual void evaluateAt(double pval, Point& pos, Point& tan) 
    {
	MESSAGE("Warning! Tried to evaluate non-evaluable IntersectionCurve.\n"
		"Nothing done.");
    }

    virtual void refine(const double& pos_tol, const double& angle_tol)
    {    
	MESSAGE("Warning! Tried to refine non-evaluable IntersectionCurve.\n"
		"Nothing done.");
    }

    virtual bool isIsocurve() const
    { return false; }

    virtual bool isDegenerated() const
    { return false; }

private:
    template<class iterator>
    NonEvaluableIntersectionCurve(const iterator begin, const iterator end)
	: IntersectionCurve(begin, end)
    {
	MESSAGE("Warning! Created a non-evaluable IntersectionCurve.");
    }
    template<class iterator> friend 
    shared_ptr<IntersectionCurve>
    constructIntersectionCurve(const iterator begin,
			       const iterator end);
};


/// IntersectionCurve defined by an isoparametric curve that can be
/// picked from the underlying object

//===========================================================================
class IsoparametricIntersectionCurve : public IntersectionCurve
//===========================================================================
{
public:
    virtual ~IsoparametricIntersectionCurve() {};

    virtual shared_ptr<ParamCurve> getCurve() const;

    virtual shared_ptr<ParamCurve> getParamCurve(int obj_nmb) const;

    virtual bool isIsocurve() const
    { return true; }

    virtual bool isDegenerated() const
    { return false; } // Maybe this should be investigated

    virtual void refine(const double& pos_tol, const double& angle_tol)
    { 
	MESSAGE("Tried to refine an isoparametric curve.  Ignoring...");   
    }

    virtual void getParamSpan(double& start, double& end) const 
    { 
	start = isopar_geom_curve_->startparam(); 
	end = isopar_geom_curve_->endparam();
    }

    virtual void evaluateAt(double pval, Point& pos, Point& tan) 
    {
	temp_.resize(2);
	isopar_geom_curve_->point(temp_, pval, 1);
	pos = temp_[0];
	tan = temp_[1];
    }

protected:
    shared_ptr<ParamCurve> isopar_geom_curve_;
    shared_ptr<ParamCurve> isopar_param_curve_1_;
    shared_ptr<ParamCurve> isopar_param_curve_2_;
    mutable std::vector<Go::Point> temp_;

    template<class iterator>
    IsoparametricIntersectionCurve(const iterator begin, const iterator end,
				   int isopar)
	: IntersectionCurve(begin, end) 
    {
	precalculate_iso_curves(isopar);
    }

    void precalculate_iso_curves(int isopar); 

    template<class iterator>
    static std::vector<int>
    resolve_isoparametric_directions(const iterator begin, 
				     const iterator end);

    template<class iterator> friend 
    shared_ptr<IntersectionCurve>
    constructIntersectionCurve(const iterator begin,
			       const iterator end);
};


/// IntersectionCurve that is defined by Hermite interpolation of a
/// number of guidepoints.

//===========================================================================
class InterpolatedIntersectionCurve : public IntersectionCurve
//===========================================================================
{
public:
    virtual ~InterpolatedIntersectionCurve() {}

    virtual void refine(const double& pos_tol, const double& angle_tol);

    virtual shared_ptr<ParamCurve> getCurve() const;

    virtual shared_ptr<ParamCurve> getParamCurve(int obj_nmb) const;

    virtual bool isIsocurve() const
    { return false; }

    virtual bool isDegenerated() const
    { return false; }

    virtual void getParamSpan(double& start, double& end) const
    {
	start = 0;
	end = 1;
    }

    virtual void evaluateAt(double pval, Point& pos, Point& tan);

    virtual bool
    getGuidePointTangent(shared_ptr<IntersectionPoint> pt,
			 Point& tan, int type) const;


private:
    // Mapping between intersection points and their tangents
    // (geometrical, parametrical, parametrical), as interpreted by
    // this curve. (We precalculate because it is unsafe to bluntly
    // rely on the getTangent() member function of the
    // IntersectionPoint.  The tangent might be undefined, have an
    // undefined orientation, or depending on whether we differentiate
    // from the right or from the left).  The tangents are calculated
    // at construction time, and are expected to remain constant
    // throughout the lifetime of the IntersectionCurve.  On
    // refinement, when a new IntersectionPoint is inserted, the
    // tangent map should be updated too.
    std::map<shared_ptr<IntersectionPoint>,
	     std::pair<Go::Array<Go::Point, 3>, bool> > tangents_;
    double certified_pos_tol_;
    double certified_angle_tol_;
    mutable std::vector<Go::Point> temp_;
    mutable bool geom_cached_;
    mutable bool par1_cached_;
    mutable bool par2_cached_;
    mutable shared_ptr<ParamCurve> cached_geom_curve_;
    mutable shared_ptr<ParamCurve> cached_param_curve_1_;
    mutable shared_ptr<ParamCurve> cached_param_curve_2_;

    template<class iterator>
    InterpolatedIntersectionCurve(const iterator begin, const iterator end)
	: IntersectionCurve(begin, end),
/* 	  certified_pos_tol_(std::numeric_limits<double>::max()), */
/* 	  certified_angle_tol_(std::numeric_limits<double>::max()), */
	  certified_pos_tol_(1.0e8),
	  certified_angle_tol_(1.0e8),
	  geom_cached_(false), 
	  par1_cached_(false), 
	  par2_cached_(false)
    {
	resolve_tangents();
    }

    // Since not all participating IntersectionPoints may have
    // uniquely defined and oriented tangents, this routine makes
    // sure, upon construction of the IntersectionCurve, that the
    // tangents are treated consistently (they are oriented the same
    // way, branch-points and higher-order points treated without
    // ambiguity, etc.).  Note that ONLY start- and endpoints are
    // allowed to be branchpoints/higher-order points.
    void resolve_tangents(); 

    // Helper function for 'resolve_tangents' and
    // 'refine_interval_recursive' that decides from which side in
    // each parameter differentiation should be carried out for a
    // specific IntersectionPoint.
    void choose_differentiation_side(std::list<shared_ptr<IntersectionPoint> >::
				     const_iterator pt) const;

    // Helper function for 'resolve_tangents()', used to determine
    // whether the list containing the intersection points should be
    // inversed or not, depending on the tangents.
    void invert_sequence_if_necessary();

    // Helper function for 'resolve_tangents()', used to search for a
    // point with a clearly defined tangent direction and orientation,
    // which can be used as a reference for setting the orientation of
    // other, orientation-less, tangents in the curve.
    std::list<shared_ptr<IntersectionPoint> >::iterator
    choose_reference_direction();

    // Helper function for 'resolve_tangents()', used to orient the
    // tangents in a range or IntersectionPoints, in accordance with
    // the tangent of a reference point, supposedly inside the given
    // range.
    void make_consistent_orientation(std::list<shared_ptr<IntersectionPoint> >::
				     iterator ref_elem);

    // Helper function that flips the tangent of a point (explicitly
    // if the point has a uniquely defined tangent, or in internal map
    // if not).
    void flip_tangent(std::list<shared_ptr<IntersectionPoint> >::iterator pt,
		      bool flip);

    // Helper function used to determine whether the orientation of
    // the tangent vector to 'mid' should be flipped in order to
    // correspond with its neighbours, 'prec' and 'next'.
    bool determine_flip(std::list<shared_ptr<IntersectionPoint> >::const_iterator prec,
			std::list<shared_ptr<IntersectionPoint> >::const_iterator mid,
			std::list<shared_ptr<IntersectionPoint> >::const_iterator next);

    // Determine the tangent of an IntersectionPoint (assumed to be a
    // member point of this IntersectionCurve.  For most purposes, the
    // tangent would be equal to the one obtained by calling
    // pt->getTangent(), but not always, since the point can be a
    // branchpoint or (worse) a higher-order point.  In those cases,
    // the tangent is determined by the curve's "individual
    // interpretation".  The possibility of an IntersectionPoint being
    // a branch/higher-order point is only present for endpoints.  It
    // is always assumed that interior points have well-defined
    // tangents.
    Point tangent_of(std::list<shared_ptr<IntersectionPoint> >::const_iterator pt) const;

    // Similar to the function above, except that it returns the
    // _parametrical_ tangent, in the parameter domain of one of the
    // intersecting objects.  If the function is called with
    // 'second_obj' set to false, then the parametrical tangent in the
    // first object is returned, and vice versa.
    Point param_tangent_of(std::list<shared_ptr<IntersectionPoint> >::const_iterator pt,
			   bool second_obj) const;

    // Helper function for 'resolve_tangents()', used to determine the
    // tangent direction (not orientation) for an IntersectionPoint
    // where this information cannot be unambiguously found in the
    // usual way (branchpoints or higher-order points), as well as
    // creating curve spesific copies of tangents of points without a
    // determined direction.
    void establish_curve_spesific_tangent(std::list<shared_ptr<IntersectionPoint> >::
					  const_iterator pt);

    // Helper function to calculate the tangent of an
    // IntersectionPoint based on the position and tangents of
    // neighbours, and without using tangent information in the
    // conserned IntersectionPoint itself.
    bool context_tangent_estimate(std::list<shared_ptr<IntersectionPoint> >::
				  const_iterator pt, 
				  TangentDomain tdom, 
				  EstimateDirection, 
				  Point& result) const;

//     // Using neighbouring points to estimate tangent at pt (necessary
//     // if pt is a higher-order point or a branchpoint).
//     Point geometric_tangent_estimate(std::list<std::
// 				     shared_ptr<IntersectionPoint> >::
// 				     const_iterator pt) const;

    // Helper function for refine() function
    void refine_interval_recursive(std::list<shared_ptr<IntersectionPoint> >::
				   iterator start_point,
				   const double& pos_tol,
				   const double& angle_tol); // Recursive
							     // function

    // Evaluates surface point and tangent.  If tangent could not be
    // determined, return 'false', else return 'true'.
    bool eval_surf_point(const Point& midpoint_pos,
			 const Point& midpoint_tan,
			 const ParamSurfaceInt* psurf1,
			 const Point& midpoint_param_pos_1,
			 const ParamSurfaceInt* psurf2,
			 const Point& midpoint_param_pos_2,
			 Point& surface_point,
			 Point& surface_tangent,
			 Point& surface_1_param,
			 Point& surface_2_param,
			 int& jstat) const;

    void hermite_interpol(std::list<shared_ptr<IntersectionPoint> >::
			  iterator start_point,
			  std::list<shared_ptr<IntersectionPoint> >::
			  iterator end_point,
			  Point& mid_position,
			  Point& mid_tangent,
			  EvalKind kind) const;

    template<class iterator> friend 
    shared_ptr<IntersectionCurve>
    constructIntersectionCurve(const iterator begin,
			       const iterator end);
};


/// Intersection between two (partially) coincident curves.

//===========================================================================
class CoincCurveIntersectionCurve : public IntersectionCurve
//===========================================================================
{
public:
  virtual ~CoincCurveIntersectionCurve() {}

  virtual void refine(const double& pos_tol, const double& angle_tol)
  { // Not relevant
  }

  virtual shared_ptr<ParamCurve> getCurve() const;

  virtual shared_ptr<ParamCurve> getParamCurve(int obj_nmb) const;

  virtual bool isIsocurve() const
  { return false; }

  virtual bool isDegenerated() const
  { return false; }

  virtual void getParamSpan(double& start, double& end) const
  {
    start = startpar_;
    end = endpar_;
  }

  virtual void evaluateAt(double pval, Point& pos, Point& tan);

  virtual bool
  getGuidePointTangent(shared_ptr<IntersectionPoint> pt,
		       Point& tan, int type) const;
  
protected:
  double startpar_;
  double endpar_;
  mutable shared_ptr<ParamCurve> cached_geom_curve_;
  mutable shared_ptr<ParamCurve> cached_param_curve_1_;
  mutable shared_ptr<ParamCurve> cached_param_curve_2_;
  mutable bool geom_cached_;
  mutable bool par1_cached_;
  mutable bool par2_cached_;
 
  template<class iterator>
  CoincCurveIntersectionCurve(const iterator begin, const iterator end)
    : IntersectionCurve(begin, end),
      geom_cached_(false), 
      par1_cached_(false), 
      par2_cached_(false)
  {
    startpar_ = ipoints_.front()->getPar1()[0];
    endpar_ = ipoints_.back()->getPar1()[0];
    if (startpar_ > endpar_)
      std::swap(startpar_, endpar_);
  }
    template<class iterator> friend 
    shared_ptr<IntersectionCurve>
    constructIntersectionCurve(const iterator begin,
			       const iterator end);
};
// Use the below function to create an IntersectionCurve

//===========================================================================
template<class iterator> 
shared_ptr<IntersectionCurve>
constructIntersectionCurve(const iterator begin,
			   const iterator end)
//===========================================================================
{
    typedef IsoparametricIntersectionCurve IsoCurve;
    typedef DegeneratedIntersectionCurve DegenCurve;
    typedef InterpolatedIntersectionCurve InterpolCurve;
    typedef NonEvaluableIntersectionCurve NonEvalCurve;
    typedef CoincCurveIntersectionCurve CoincCurve;

    // Do NOT override the below code.  curves that are degenerated
    // should NOT be treated as InterpolatedIntersectionCurve!!!!
    // Causes trouble and crashes.
//     if (0) { // @@sbr Degenerated space-curve seems to be the trend
// 	     // for a function ...
    std::vector<int> isopar
	= IsoCurve::resolve_isoparametric_directions(begin, end);
    int npar = (*begin)->numParams1() + (*begin)->numParams2();
    if ((*begin)->numParams1() == 1 && (*begin)->numParams2() == 1)
      {
	// Curve-curve intersection. Pick sub part of first curve
	try {
	    return shared_ptr<IntersectionCurve>
		(new CoincCurve(begin, end));
	//} catch (Zero_Parameter_Span_Error& e) {
	} catch (...) {
	  return shared_ptr<IntersectionCurve>
	    (new NonEvalCurve(begin, end));
	}
      }
    else if (isopar.size() > 0 && (int)isopar.size() < npar) {
	
	// The curve is isoparametric
	try {
	    return shared_ptr<IntersectionCurve>
		(new IsoCurve(begin, end, isopar[0]));
	//} catch (Zero_Parameter_Span_Error& e) {
	} catch (Zero_Parameter_Span_Error) {
	    // Degenerated parameter span.  Making degenerated curve
	    // instead.
	    MESSAGE("Warning: tried to generate an isocurve but made"
		    " a degenerated curve instead.");
	    return shared_ptr<IntersectionCurve>
		(new DegenCurve(begin, end));
	}
    } else if (DegenCurve::degenerated_range(begin, end)) {
	//@@sbr Temp!!! 0) { Currently not handling functions (i.e. 1-dim).
	// the curve is degenerated
	return shared_ptr<IntersectionCurve>
	    (new DegenCurve(begin, end));
    }
    
    // If we got here, we will try to make a normal, interpolated
    // IntersectionCurve
    shared_ptr<IntersectionCurve> res;
    try {
	res = shared_ptr<IntersectionCurve>
	    (new InterpolCurve(begin, end));
    } catch(...) {
	// Could not make InterpolatedIntersectionCurve.  Probably
	// problem with tangents.  Making curve without evaluation
	// functionality instead.
	res = shared_ptr<IntersectionCurve>
	    (new NonEvalCurve(begin, end));
    }
    return res;
}


//===========================================================================
template<class iterator> std::vector<int> IsoparametricIntersectionCurve::  
resolve_isoparametric_directions(const iterator begin, const iterator end)
//===========================================================================
{
    // This is a static function

    std::vector<shared_ptr<IntersectionPoint> > points(begin, end);
    ASSERT(points.size() >= 2);

    int num_par = points.front()->numParams1()
	+ points.front()->numParams2();
    std::vector<bool> iso_candidates(num_par, true); // Candidates for
						     // isoparametric
    
    // Disqualify candidates
    int i;
    for (i = 0; i < int(points.size()) - 1; ++i) {
	// Keep only candidates that are isoparametric
	shared_ptr<IntersectionLink> link
	    = points[i]->getIntersectionLink(points[i+1].get());
	ASSERT(link.get() != 0);
	for (int par = 0; par < num_par; ++par) {
	    if (!link->isIsoparametricIn(par)) {
		iso_candidates[par] = false;
	    }
	}
    }
    
    // We have now eliminated all directions that are not
    // isoparametric for all links in the curve
    std::vector<int> result;
    for (i = 0; i < num_par; ++i) {
	if (iso_candidates[i]) {
	    result.push_back(i);
	}
    }
    return result;
}

//===========================================================================
template<class iterator> bool DegeneratedIntersectionCurve::
degenerated_range(const iterator begin, const iterator end)
//===========================================================================
{
    const double tol = (*begin)->getTolerance()->getEpsge();

    // If all points are within a geometric distance of 'tol' to the
    // first point, then the range is degenerated

    Point ref_pos = (*begin)->getPoint();
    bool found = false;
    for (iterator i = begin; i != end && !found; ++i) {
	Point cur_pos = (*i)->getPoint();
	found = (ref_pos.dist2(cur_pos) > tol * tol);
    }
    return !found;
}


}; // namespace Go;


#endif  // _INTERSECTIONCURVE_H

