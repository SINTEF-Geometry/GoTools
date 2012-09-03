//===========================================================================
//                                                                           
// File: CurveLoop.C                                                       
//                                                                           
// Created: Thu Mar 15 09:13:01 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: CurveLoop.C,v 1.28 2008-12-01 14:01:45 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#define DEBUG

#include <memory>
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/orientCurves.h"
#include <fstream>


using namespace Go;
using std::vector;
using std::max;
using std::min;
using std::endl;
using std::cerr;


namespace Go {


//===========================================================================
CurveLoop::CurveLoop()
  : space_epsilon_(-1.0), valid_state_(0)
//===========================================================================
{
}

//===========================================================================
CurveLoop::CurveLoop(const std::vector< shared_ptr<ParamCurve> >& curves,
		     double space_epsilon)
    : valid_state_(0)
//===========================================================================
{
    setSpaceEpsilon(space_epsilon);
    setCurves(curves);
}




//===========================================================================
CurveLoop::~CurveLoop()
//===========================================================================
{
}


//===========================================================================
void CurveLoop::swap(CurveLoop& other)
//===========================================================================
{
    curves_.swap(other.curves_);
    std::swap(space_epsilon_, other.space_epsilon_);
}



//===========================================================================
void
CurveLoop::setCurves(const std::vector<shared_ptr<ParamCurve> >& curves)
//===========================================================================
{
    if (curves.empty()) {
	THROW("Loop must contain at least one curve");
    }
    double maxdist = computeLoopGap(curves);
//     if (maxdist > space_epsilon_) {
// 	THROW("Distance between curve-ends is larger than given epsilon: " 
// 	      << maxdist << " > " << space_epsilon_ );
//     }
    if (maxdist > space_epsilon_)
      {
	valid_state_ = -1;
	MESSAGE("Distance between curve-ends is larger than given epsilon: " 
		<< maxdist << " > " << space_epsilon_ <<
		". Creating invalid CurveLoop.");
#ifdef DEBUG
	  std::ofstream out_file("curve_loop.g2");
	  for (size_t kj=0; kj<curves.size(); ++kj)
	    {
	      shared_ptr<CurveOnSurface> sf_cv =
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(curves[kj]);
	      shared_ptr<ParamCurve> cv;
	      if (sf_cv.get())
		cv = sf_cv->spaceCurve();
	      else
		cv = curves[kj];
	      cv->writeStandardHeader(out_file);
	      cv->write(out_file);
	    }
#endif
 
      }
    else
      valid_state_ = 1;

    curves_ = curves;

    // Try to fix
    if (valid_state_ < 0)
      {
	fixInvalidLoop(maxdist);
	if (maxdist <= space_epsilon_)
	  {
	    valid_state_ = 1;
	    MESSAGE("Loop fixed");
	  }
      }

}


//===========================================================================
void CurveLoop::setSpaceEpsilon(const double space_epsilon)
//===========================================================================
{
    ALWAYS_ERROR_IF(space_epsilon < 0.0, "Space epsilon smaller than 0");

    MESSAGE_IF(space_epsilon > 1.0,
	       "Rather large space epsilon... space_eps = " << space_epsilon);

    space_epsilon_ = space_epsilon;

    if (!curves_.empty())
	setCurves(curves_);    // check that this space epsilon works OK
}


//===========================================================================
double CurveLoop::getSpaceEpsilon() const
//===========================================================================
{
    return space_epsilon_;
}


//===========================================================================
void CurveLoop::turnOrientation()
//===========================================================================
{
  int ki;
  int nmb_curves = (int)curves_.size();
  for ( ki=0; ki<nmb_curves; ki++)
    curves_[ki]->reverseParameterDirection();
  // In order for the curves to form a continuous loop, we turn the vector.
  for ( ki = 0; ki < nmb_curves/2; ++ki)
      std::swap(curves_[ki], curves_[nmb_curves-1-ki]);
}

//===========================================================================
shared_ptr<ParamCurve> CurveLoop::operator[] (int index) const
//===========================================================================
{
    return curves_[index];
}

//===========================================================================
void CurveLoop::closestPoint(const Point& pt, int& clo_ind, double& clo_par, 
			       Point& clo_pt, double& clo_dist) const
//===========================================================================
{
    clo_ind = 0;
    double tmp_par, tmp_dist;
    Point tmp_pt;
    curves_[0]->closestPoint(pt, clo_par, clo_pt, clo_dist);
    size_t ki;
    for (ki=1; ki < curves_.size(); ki++) {
	curves_[ki]->closestPoint(pt, tmp_par, tmp_pt, tmp_dist);
	if (tmp_dist < clo_dist) {
	    clo_dist = tmp_dist;
	    clo_pt = tmp_pt;
	    clo_par = tmp_par;
	    clo_ind = (int)ki;
	}
    }
}

//===========================================================================
void CurveLoop::closestParPoint(const Point& pt, int& clo_ind, 
				  double& clo_par, Point& clo_pt, 
				  double& clo_dist) const
//===========================================================================
{
  clo_ind = -1;
  double tmp_par, tmp_dist;
  Point tmp_pt;
  size_t ki;
  shared_ptr<CurveOnSurface> curr_crv;
  shared_ptr<ParamCurve> par_crv;
  for (ki=0; ki<curves_.size(); ki++)
    {
      if (curves_[ki]->instanceType() != Class_CurveOnSurface)
	continue;

      curr_crv = dynamic_pointer_cast<CurveOnSurface, ParamCurve>(curves_[ki]);
      par_crv = curr_crv->parameterCurve();
      if (par_crv.get() == 0)
	continue;

      par_crv->closestPoint(pt, clo_par, clo_pt, clo_dist);
      clo_ind = (int)ki;
      break;
    }

  for (ki++; ki<curves_.size(); ki++)
    {
      if (curves_[ki]->instanceType() != Class_CurveOnSurface)
	continue;

      curr_crv = dynamic_pointer_cast<CurveOnSurface, ParamCurve>(curves_[ki]);
      par_crv = curr_crv->parameterCurve();
      if (par_crv.get() == 0)
	continue;

      par_crv->closestPoint(pt, tmp_par, tmp_pt, tmp_dist);
      if (tmp_dist < clo_dist)
	{
	  clo_dist = tmp_dist;
	  clo_pt = tmp_pt;
	  clo_par = tmp_par;
	  clo_ind = (int)ki;
	}
    }

}


//===========================================================================
    /// Return joint points between curves
  vector<Point> CurveLoop::getCorners() const
//===========================================================================
  {
    vector<Point> res;
    for (size_t ki=0; ki<curves_.size(); ++ki)
      {
	Point curr = curves_[ki]->point(curves_[ki]->startparam());
	res.push_back(curr);
      }

    return res;
  }

  
//===========================================================================
void
CurveLoop::getSmoothCurves(vector<vector<shared_ptr<ParamCurve> > >& curves,
			   double angtol)
//===========================================================================
{
  curves.clear();

  // Find first corner
  int ki, kj, kr;
  int idx = -1;
  int nmb = (int)curves_.size();
  for (ki=nmb-1, kj=0, kr=0; kr<nmb; ki=(ki+1)%nmb, kj=(kj+1)%nmb, ++kr)
    {
      vector<Point> pt1(2);
      vector<Point> pt2(2);
      curves_[ki]->point(pt1, curves_[ki]->endparam(), 1);
      curves_[kj]->point(pt2, curves_[kj]->startparam(), 1);
      if (pt1[1].angle(pt2[1]) > angtol)
	{
	  idx = kj;
	  break;
	}
    }
     
  if (idx < 0)
    {
      // Smooth loop
      curves.push_back(curves_);
    }
  else
    {
      for (ki=idx, kj=(idx+1)%nmb; ; ki=(ki+1)%nmb, kj=(kj+1)%nmb)
	{
	  vector<shared_ptr<ParamCurve> > curr_cvs;
	  curr_cvs.push_back(curves_[ki]);
	  for (; kj!=idx; ki=(ki+1)%nmb, kj=(kj+1)%nmb)
	    {
	      vector<Point> pt1(2);
	      vector<Point> pt2(2);
	      curves_[ki]->point(pt1, curves_[ki]->endparam(), 1);
	      curves_[kj]->point(pt2, curves_[kj]->startparam(), 1);
	      if (pt1[1].angle(pt2[1]) > angtol)
		break;
	      curr_cvs.push_back(curves_[kj]);
	    }
	  curves.push_back(curr_cvs);
	  if (kj == idx)
	    break;
	}
    }
}

//===========================================================================
bool CurveLoop::isValid() const
//===========================================================================
{
    return (valid_state_ == 1);
}


//===========================================================================
bool CurveLoop::fixInvalidLoop(double& max_gap)
//===========================================================================
{
    if (valid_state_ == 1) {
	return true; // Nothing to be done.
    }

#ifdef SBR_DBG
    std::cout << "valid_state_ = " << valid_state_ << std::endl;
#endif

    // Make copy
    vector<shared_ptr<ParamCurve> > curves(curves_.size());
    for (size_t ki=0; ki<curves_.size(); ++ki)
      curves[ki] = shared_ptr<ParamCurve>(curves_[ki]->clone());

    max_gap = computeLoopGap(curves);
    double max_gap2 = max_gap;

    vector<shared_ptr<ParamCurve> > par_cvs, space_cvs;
    for (size_t ki = 0; ki < curves.size(); ++ki)
	if (curves[ki]->instanceType() == Class_CurveOnSurface) {
	    shared_ptr<CurveOnSurface> cv_on_sf =
		dynamic_pointer_cast<CurveOnSurface, ParamCurve>(curves[ki]);
	    par_cvs.push_back(cv_on_sf->parameterCurve());
	    space_cvs.push_back(cv_on_sf->spaceCurve());
	}

    // We analyze the loop. Possibly the order of the curves is
    // incorrect. In that case we reorder the curves.  Or perhaps we
    // need to change the direction of one/some of the loops.  Whether
    // the loop is ccw or not must be handled on the outside.
    if (max_gap > space_epsilon_) {
	// We first check if there exists a space curve with
	// legal definition (or parametric of we prefer space).
	// boundaries is a vector of CurveOnSurface.
	double maxgap_par = computeLoopGap(par_cvs);
	double maxgap_space = computeLoopGap(space_cvs);
	if ((maxgap_par >= 0.0) &&
	    (maxgap_par < space_epsilon_)) {
#ifdef SBR_DBG
	    std::ofstream debug("tmp/debug_loops.g2");
	    std::cout << "Loop size: " << curves.size() << std::endl;
#endif
	    for (int k = 0; k < (int)curves.size(); ++k)
	    {
#ifdef SBR_DBG
		if (par_cvs[k] != NULL) {
		    par_cvs[k]->writeStandardHeader(debug);
		    par_cvs[k]->write(debug);
		}
		if (space_cvs[k] != NULL) {
		    space_cvs[k]->writeStandardHeader(debug);
		    space_cvs[k]->write(debug);
		}
#endif
		if (curves[k]->instanceType() == Class_CurveOnSurface) {
		    shared_ptr<CurveOnSurface> cv_on_sf =
			dynamic_pointer_cast<CurveOnSurface, ParamCurve>
			(curves[k]);
		    cv_on_sf->makeCurvesConsistent(true);
		    curves[k] =
		      shared_ptr<CurveOnSurface>
		      (new CurveOnSurface
		       (cv_on_sf->underlyingSurface(),
			cv_on_sf->parameterCurve(),
			cv_on_sf->spaceCurve(),
			true));
		}
	    }
	} else if ((maxgap_space >= 0.0) &&
		   (maxgap_space < space_epsilon_)) {
	    for (int k = 0; k < (int)curves.size(); ++k)
	    {
		if (curves[k]->instanceType() == Class_CurveOnSurface) {
		    shared_ptr<CurveOnSurface> cv_on_sf =
			dynamic_pointer_cast<CurveOnSurface, ParamCurve>
			(curves[k]);
		    //cv_on_sf->makeCurvesConsistent(true);
		    cv_on_sf->makeCurvesConsistent(false);
		    curves[k] =
		      shared_ptr<CurveOnSurface>
		      (new CurveOnSurface
		       (cv_on_sf->underlyingSurface(),
			cv_on_sf->parameterCurve(),
			cv_on_sf->spaceCurve(),
			false));
		}
	    }
	} else {
	    // Try to fix by rearranging the segments.
	    vector<int> perm;
	    vector<bool> flip;
	    orientCurves::orientCurves(curves, perm, flip,
				       space_epsilon_, false);
	    // Making the new boundary vector
	    vector< shared_ptr<ParamCurve> > new_boundary;
	    new_boundary.reserve(curves.size());
	    for (size_t bi = 0; bi < curves.size(); ++bi) {
		new_boundary.push_back(curves[perm[bi]]);
		if (flip[bi]) {
		    new_boundary[bi]->reverseParameterDirection();
		}
	    }
	    curves.swap(new_boundary);
	    // We check if that helped.
	    max_gap2 = Go::computeLoopGap(curves);
	    if (max_gap2 > space_epsilon_) {
		cerr << "Gap > space_epsilon_: " << max_gap2 << " > "
		     << space_epsilon_ << endl;
		MESSAGE("Cannot fix boundary that does not form a loop.");
	    }
	}
    }

    // If we're still outside we try another approach.
    for (size_t ki = 0; ki < curves.size(); ++ki)
	if (curves[ki]->instanceType() == Class_CurveOnSurface) {
	    shared_ptr<CurveOnSurface> cv_on_sf =
		dynamic_pointer_cast<CurveOnSurface, ParamCurve>(curves[ki]);
	    par_cvs[ki] = cv_on_sf->parameterCurve();
	    space_cvs[ki] = cv_on_sf->spaceCurve();
	}
    double maxgap_par = computeLoopGap(par_cvs);
    double maxgap_space = computeLoopGap(space_cvs);
    if ((maxgap_space >= 0.0) &&
	(maxgap_space < space_epsilon_)) {
	for (int k = 0; k < (int)curves.size(); ++k)
	    if (curves[k]->instanceType() == Class_CurveOnSurface) {
		shared_ptr<CurveOnSurface> cv_on_sf =
		    dynamic_pointer_cast<CurveOnSurface, ParamCurve>
		    (curves[k]);
		curves[k] =
		    shared_ptr<CurveOnSurface>
		    (new CurveOnSurface
		     (cv_on_sf->underlyingSurface(),
		      cv_on_sf->parameterCurve(),
		      cv_on_sf->spaceCurve(),
		      false));
	    }
    } else if ((maxgap_par >= 0.0) &&
	       (maxgap_par < space_epsilon_)) {
	for (int k = 0; k < (int)curves.size(); ++k)
	    if (curves[k]->instanceType() == Class_CurveOnSurface) {
		shared_ptr<CurveOnSurface> cv_on_sf =
		    dynamic_pointer_cast<CurveOnSurface, ParamCurve>
		    (curves[k]);
		curves[k] =
		    shared_ptr<CurveOnSurface>
		    (new CurveOnSurface
		     (cv_on_sf->underlyingSurface(),
		      cv_on_sf->parameterCurve(),
		      cv_on_sf->spaceCurve(),
		      true));
	    }
    }

    max_gap2 = computeLoopGap(curves);

    if (max_gap2 < max_gap)
      {
	curves_ = curves;
	max_gap = max_gap2;
      }
    return true;
}


//===========================================================================
bool CurveLoop::simplify(double tol, double ang_tol, double& max_dist)
//===========================================================================
{
  if (curves_.size() <= 2)
    return false;   // Already very few curves

  max_dist = 0.0;

  bool modified = false;
  size_t ki;
  vector<Point> der1(2);
  vector<Point> der2(2);
  double dist;
      
  for (ki=1; ki<curves_.size(); ++ ki)
    {
      // Check if the two curves may be joined
      curves_[ki-1]->point(der1, curves_[ki-1]->endparam(), 1);
      curves_[ki]->point(der2, curves_[ki]->startparam(), 1);

      if (der1[0].dist(der2[0]) > tol || der1[1].angle(der2[1]) > ang_tol)
	continue;  // Not smooth

      // Append curves 
      shared_ptr<ParamCurve> cv1 = shared_ptr<ParamCurve>(curves_[ki-1]->clone());
      shared_ptr<ParamCurve> cv2 = shared_ptr<ParamCurve>(curves_[ki]->clone());

      cv1->appendCurve(cv2.get(), 1, dist, true);
      if (dist > tol)
	continue;  // Error not within tolerance

      modified = true;  // Joining performed

      max_dist = std::max(max_dist, dist);

      // Replace curves in the curve loop
      curves_[ki-1] = cv1;
      curves_.erase(curves_.begin()+ki);

      ki--;
    }
      
  // Check whether the first and last curve may be joined
  curves_[curves_.size()-1]->point(der1, curves_[curves_.size()-1]->endparam(), 1);
  curves_[0]->point(der2, curves_[0]->startparam(), 1);

  if (der1[0].dist(der2[0]) <= tol && der1[1].angle(der2[1]) <= ang_tol)
    {
      // Append curves 
      shared_ptr<ParamCurve> cv1 = 
	shared_ptr<ParamCurve>(curves_[curves_.size()-1]->clone());
      shared_ptr<ParamCurve> cv2 = shared_ptr<ParamCurve>(curves_[0]->clone());

      cv1->appendCurve(cv2.get(), 1, dist, true);
      if (dist <= tol)
	{
	  modified = true;
	  curves_[0] = cv1;
	  curves_.erase(curves_.end()-1);
	}
    }

  return modified;
}

// //===========================================================================
// double CurveLoop::maxGap(int nmb_seg_samples)
// //===========================================================================
// {
//     double max_loop_gap = computeLoopGap(curves_);

//     double max_gap = max(max_loop_gap, max_sf_dist);
//     return max_gap;
// }


}; // end namespace Go

