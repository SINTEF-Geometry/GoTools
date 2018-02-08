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

#include <memory>
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/orientCurves.h"
#include "GoTools/geometry/SplineCurve.h"
#include <fstream>

//#define DEBUG

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
		     double space_epsilon, bool allow_fix)
    : valid_state_(0)
//===========================================================================
{
    setSpaceEpsilon(space_epsilon);
    setCurves(curves, allow_fix);
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
CurveLoop::setCurves(const std::vector<shared_ptr<ParamCurve> >& curves,
		     bool allow_fix)
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
	MESSAGE("Distance between curve-ends is larger than given epsilon: " 
		<< maxdist << " > " << space_epsilon_ <<
		". Creating invalid CurveLoop.");
#ifndef NDEBUG
	{
	    std::ofstream out_file("tmp/curve_loop.g2");
	    for (size_t kj=0; kj<curves.size(); ++kj)
	    {
		shared_ptr<CurveOnSurface> sf_cv =
		    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(curves[kj]);
		if (sf_cv.get() != NULL)
		{
		    std::ofstream out_file2("tmp/under_sf.g2");
		    shared_ptr<ParamSurface> under_sf = sf_cv->underlyingSurface();
		    if (under_sf.get() != NULL)
		    {
			under_sf->writeStandardHeader(out_file2);
			under_sf->write(out_file2);
		    }
		}
		shared_ptr<ParamCurve> cv;
		if (sf_cv.get())
		    cv = sf_cv->spaceCurve();
		else
		    cv = curves[kj];
		if (cv.get())
		{
		    cv->writeStandardHeader(out_file);
		    cv->write(out_file);
		}
	    }
	    double debug_val = 0.0;
	}
#endif
	  valid_state_ = -1;
      }
    else
      valid_state_ = 1;

    curves_ = curves;

    // Try to fix
    if (valid_state_ < 0 && allow_fix)
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

    // MESSAGE_IF(space_epsilon > 1.0,
    // 	       "Rather large space epsilon... space_eps = " << space_epsilon);

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
		try {
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
			// We check if that helped.
			max_gap2 = Go::computeLoopGap(curves);
                        if (max_gap2 < maxgap_space)
                        {
                            curves.swap(new_boundary);
                            if (max_gap2 > space_epsilon_) {
                                MESSAGE("Gap > space_epsilon_: " << max_gap2 << " > "
                                        << space_epsilon_ 
                                        << ". Cannot fix boundary that does not form a loop.");
                            }
                        }
                        else
                        {
                            // Most likely the error is due to bad input (gap larger than tolerance).
                            MESSAGE("Failed finding smaller dist by rearranging the segments.");

                        }
		} catch (...) {
		    MESSAGE("Method failed: orientCurves()");
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

    // If we got closer we replace the loop curves.
    max_gap2 = computeLoopGap(curves);
    if (max_gap2 < max_gap)
      {
	curves_ = curves;
	max_gap = max_gap2;
      }

#if 1
    return true;
#else
    // @@sbr201710 This version is correct. But if the gap is too large it is most likely an input problem.
    // The solution should then be to fix the tolerance (or the gap) prior to calling this routine.
    bool is_valid = (max_gap < space_epsilon_);

    return is_valid;
#endif
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

      // Check for rational spline curves. In that case only C0 continuity
      // should be requested
      SplineCurve *spline1 = cv1->getSplineCurve();
      SplineCurve *spline2 = cv2->getSplineCurve();
      bool rat1 = false, rat2 = false;
      if (spline1)
	rat1 = spline1->rational();
      if (spline2)
	rat2 = spline1->rational();
      try {
	cv1->appendCurve(cv2.get(), (rat1 || rat2) ? 0 : 1, dist, true);
      }
      catch (...)
	{
	  continue;
	}
      if (dist > tol || dist < 0.0)
	continue;  // Error not within tolerance or append not performed

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

      // Check for rational spline curves. In that case only C0 continuity
      // should be requested
      SplineCurve *spline1 = cv1->getSplineCurve();
      SplineCurve *spline2 = cv2->getSplineCurve();
      bool rat1 = false, rat2 = false;
      if (spline1)
	rat1 = spline1->rational();
      if (spline2)
	rat2 = spline1->rational();
      try {
	cv1->appendCurve(cv2.get(), (rat1 || rat2) ? 0 : 1, dist, true);
      }
      catch (...)
	{
	  dist = HUGE;
	}
      if (dist <= tol && dist >= 0.0)
	{
	  modified = true;
	  max_dist = std::max(max_dist, dist);
	  curves_[0] = cv1;
	  curves_.erase(curves_.end()-1);
	}
    }

  return modified;
}

//===========================================================================
  vector<shared_ptr<ParamCurve> > CurveLoop::split(int idx, double par)
//===========================================================================
{
  vector<shared_ptr<ParamCurve> > sub_cvs;
  if (idx < 0 || idx >= (int)curves_.size())
    return sub_cvs;

  sub_cvs = curves_[idx]->split(par);
  if (sub_cvs.size() > 0)
    {
      curves_.erase(curves_.begin()+idx);
      curves_.insert(curves_.begin()+idx, sub_cvs.begin(), sub_cvs.end());
    }

  return sub_cvs;
}

//===========================================================================
  int CurveLoop::removeCrvAndFix(shared_ptr<CurveOnSurface> cv)
//===========================================================================
{
  // Find curve in loop
  for (size_t ki=0; ki<curves_.size(); ++ki)
    {
      if (curves_[ki].get() == cv.get())
	{
	  shared_ptr<ParamCurve> adjcv1 = (ki > 0) ? curves_[ki-1] :
	    curves_[curves_.size()-1];
	  shared_ptr<ParamCurve> adjcv2 = (ki < curves_.size()-1) ? 
	    curves_[ki+1] : curves_[0];
	  curves_.erase(curves_.begin()+ki);
	  shared_ptr<CurveOnSurface> sf_cv1 = 
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(adjcv1);
	  shared_ptr<CurveOnSurface> sf_cv2 = 
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(adjcv2);
	  if (sf_cv1.get() && sf_cv2.get())
	    {
#ifdef DEBUG
	      std::ofstream of1("replace_space.g2");
	      sf_cv1->spaceCurve()->writeStandardHeader(of1);
	      sf_cv1->spaceCurve()->write(of1);
	      sf_cv2->spaceCurve()->writeStandardHeader(of1);
	      sf_cv2->spaceCurve()->write(of1);
	      std::ofstream of2("replace_par.g2");
	      sf_cv1->parameterCurve()->writeStandardHeader(of2);
	      sf_cv1->parameterCurve()->write(of2);
	      sf_cv2->parameterCurve()->writeStandardHeader(of2);
	      sf_cv2->parameterCurve()->write(of2);
#endif
	      Point pos1 = adjcv1->point(adjcv1->endparam());
	      Point pos2 = adjcv2->point(adjcv2->startparam());
	      Point pos = 0.5*(pos1+pos2);
	      bool replacefirst = sf_cv1->replaceEndPoint(pos, false);
	      bool replacesecond = sf_cv2->replaceEndPoint(pos, true);
#ifdef DEBUG
	      sf_cv1->spaceCurve()->writeStandardHeader(of1);
	      sf_cv1->spaceCurve()->write(of1);
	      sf_cv2->spaceCurve()->writeStandardHeader(of1);
	      sf_cv2->spaceCurve()->write(of1);
	      sf_cv1->parameterCurve()->writeStandardHeader(of2);
	      sf_cv1->parameterCurve()->write(of2);
	      sf_cv2->parameterCurve()->writeStandardHeader(of2);
	      sf_cv2->parameterCurve()->write(of2);
#endif
	      if (replacefirst && replacesecond)
		return 2;
	      else
		return 1;
	    }
	  else
	    return 1;
	}
    }
  return 0;
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

