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

#include "GoTools/geometry/BoundedSurface.h"

#include "GoTools/utils/Logger.h"
#include "GoTools/utils/Array.h"
#include "GoTools/utils/MatrixXD.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/SplineUtils.h"
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/SplineDebugUtils.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/ElementaryCurve.h"
#include "GoTools/geometry/GoIntersections.h"
#include <fstream>

//#define DEBUG

using namespace Go;
using std::vector;
using std::swap;
using std::max;
using std::min;
using std::pair;
using std::make_pair;
using std::streamsize;
using std::endl;


//#define CHECK_PARAM_LOOP_ORIENTATION

#ifndef NDEBUG
#define SBR_DBG
#include "GoTools/geometry/SplineDebugUtils.h"
#endif

//===========================================================================
BoundedSurface::BoundedSurface()
  : ParamSurface(), surface_(NULL), iso_trim_(false), iso_trim_tol_(-1.0), valid_state_(0)
//===========================================================================
{
}


//===========================================================================
BoundedSurface::BoundedSurface(shared_ptr<ParamSurface> surf,
			       vector<shared_ptr<CurveOnSurface> > loop,
			       double space_epsilon,
			       bool fix_trim_cvs)
  : ParamSurface(), surface_(surf), iso_trim_(false), iso_trim_tol_(-1.0), valid_state_(0)
//===========================================================================
{
    ALWAYS_ERROR_IF(loop.size() == 0, "Empty loop.");
    ALWAYS_ERROR_IF(surf.get() == 0, "Missing surface.");

        // Check if the surfaces in this bounded surface and in the
        // curves-on-surface instances describing the boundary are
        // the same. Exactly the same surface is required, copies
        // are not allowed.

    for (size_t i = 0; i < loop.size(); ++i) {
	DEBUG_ERROR_IF(!(surf == loop[i]->underlyingSurface()),
		 "Inconsistent surface pointers.");
    }


#ifdef CHECK_PARAM_LOOP_ORIENTATION
    // Check that the outer loop is CCW in the parameter plane, and that
    // the inner loops are CW.
//     double pareps = space_epsilon*1e-6; // @afr: Should we do something better here?
    double int_tol = 1e-6; // Used by the intersection algorithm.
    bool ccw = LoopUtils::paramIsCCW(loop, int_tol);
    if (!ccw) {
	THROW("Outer loop not CCW in the parameter plane.");
    }
#endif
    
    // Create the outer boundary loop. First make ParamCurve pointers.
    vector<shared_ptr<ParamCurve> > curves;
    for (size_t i=0; i< loop.size(); i++) 
      {
	if (fix_trim_cvs)
	  {
	    // Try to generate the parameter curve if it does not
	    // exist already
	    (void)loop[i]->ensureParCrvExistence(space_epsilon);

	    if (!loop[i]->sameCurve(space_epsilon))
	      {
		if (loop[i]->parPref())
		  {
		    shared_ptr<ParamCurve> tmp_cv = loop[i]->spaceCurve();
		    loop[i]->unsetSpaceCurve();
		    loop[i]->ensureSpaceCrvExistence(space_epsilon);
		    if (!loop[i]->spaceCurve().get())
		      loop[i]->setSpaceCurve(tmp_cv);
		  }
		else
		  {
		    shared_ptr<ParamCurve> tmp_cv = loop[i]->parameterCurve();
		    loop[i]->unsetParameterCurve();
		    bool found = loop[i]->ensureParCrvExistence(space_epsilon);
		    if (!found)
		      loop[i]->setParameterCurve(tmp_cv);
		  }
	      }
	  }
	curves.push_back(loop[i]);
      }

#ifndef NDEBUG
	{
	    std::ofstream debug("tmp/cvs_on_sf.g2");
	    SplineDebugUtils::writeCvsOnSf(curves, space_epsilon, debug);
	    double debug_val = 0.0;
	}
#endif //NDEBUG

    boundary_loops_.push_back(
      shared_ptr<CurveLoop>(new CurveLoop(curves, space_epsilon)));

    // Parameter curves may be placed on the wrong side of the seam
    // of closed surfaces. This cannot be distinguished locally during
    // creation. Make a check and repair if necessary
    if (fix_trim_cvs)
    {
    	(void)checkParCrvsAtSeam();
    }
    
    if (fix_trim_cvs)
    {
      // We then analyze the loops and set valid_state_.
      analyzeLoops();
    }
}

//===========================================================================
BoundedSurface::
BoundedSurface(shared_ptr<ParamSurface> surf,
	       vector<vector<shared_ptr<CurveOnSurface> > > loops,
	       double space_epsilon,
	       bool fix_trim_cvs)
    : ParamSurface(), surface_(surf), iso_trim_(false), iso_trim_tol_(-1.0), valid_state_(0)
//===========================================================================
{
    // This form of the constructor exists for backwards
    // compatibility. The preferred form is the one that uses a vector
    // of space_epsilons, thus treating each loop on its
    // own. Technically, in order to avoid code duplication, we call
    // contructor_implementation(). @jbt

    int nloops = (int)loops.size();
    vector<double> space_epsilons(nloops, space_epsilon);
    constructor_implementation(surf, loops, space_epsilons, fix_trim_cvs);
}

//===========================================================================
BoundedSurface::
BoundedSurface(shared_ptr<ParamSurface> surf,
	       vector<vector<shared_ptr<CurveOnSurface> > > loops,
	       vector<double> space_epsilons,
	       bool fix_trim_cvs)
    : ParamSurface(), surface_(surf), iso_trim_(false), iso_trim_tol_(-1.0), valid_state_(0)
//===========================================================================
{
    // The code in this constructor has been moved into
    // contructor_implementation() in order to avoid code
    // duplication. @jbt

    constructor_implementation(surf, loops, space_epsilons, fix_trim_cvs);
}

//===========================================================================
void BoundedSurface::
constructor_implementation(shared_ptr<ParamSurface> surf,
			   vector<vector<shared_ptr<CurveOnSurface> > > loops,
			   vector<double> space_epsilons,
			   bool fix_trim_cvs)
//===========================================================================
{
    // This function makes it possible to have overloading of two
    // nearly equal constructors. @jbt

    ALWAYS_ERROR_IF(loops.size() == 0, "Empty loop.");
    ALWAYS_ERROR_IF(surf.get() == 0, "Missing surface.");

        // Check if the surfaces in this bounded surface and in the
        // curves-on-surface instances describing the boundaries are
        // the same. Exactly the same surface is required, copies
        // are not allowed.
    // bool pref_par = true;

    for (size_t j=0; j<loops.size(); j++)
      for (size_t i=0; i<loops[j].size(); i++)
	{
	  shared_ptr<ParamSurface> sf = loops[j][i]->underlyingSurface();
	  if (!(surf.get() == sf.get()))
	      ALWAYS_ERROR_IF(!(surf.get() == sf.get()),
			      "Inconsistent surface pointers.");
	}

#ifdef CHECK_PARAM_LOOP_ORIENTATION
    // Check that the outer loop is CCW in the parameter plane, and that
    // the inner loops are CW.
    double int_tol = 1e-6; // Used by the intersection algorightm.
    bool ccw = LoopUtils::paramIsCCW(loops[0], int_tol);
    if (!ccw) {
	THROW("Outer loop not CCW in the parameter plane.");
    }
    for (int loop_index = 1; loop_index < loops.size(); ++loop_index) {
// 	double pareps = space_epsilon[loop_index] * 1.0e-6;
	ccw = LoopUtils::paramIsCCW(loops[loop_index], int_tol);
	if (ccw) {
	    THROW("Inner loop not CW in the parameter plane.");
	}
    }
#endif
    // Create the boundary loops
    for (size_t j=0; j<loops.size(); j++)
    {
	// Make ParamCurve pointers
	vector<shared_ptr<ParamCurve> > curves;
	for (size_t i=0; i< loops[j].size(); i++) {
	    if (fix_trim_cvs)
	    {
		// Try to generate the parameter curve if it does not
		// exist already
		(void)loops[j][i]->ensureParCrvExistence(space_epsilons[j]);

		if (!loops[j][i]->sameCurve(space_epsilons[j]))
		  {
		    if (loops[j][i]->parPref())
		      {
			shared_ptr<ParamCurve> tmp_cv = 
			  loops[j][i]->spaceCurve();
			loops[j][i]->unsetSpaceCurve();
			loops[j][i]->ensureSpaceCrvExistence(space_epsilons[j]);
			if (!loops[j][i]->spaceCurve().get())
			  loops[j][i]->setSpaceCurve(tmp_cv);
		      }
		    else if (surface_->dimension() == 3)
		      {
			shared_ptr<ParamCurve> tmp_cv = 
			  loops[j][i]->parameterCurve();
			loops[j][i]->unsetParameterCurve();
			bool found = 
			  loops[j][i]->ensureParCrvExistence(space_epsilons[j]);
			if (!found)
			  loops[j][i]->setParameterCurve(tmp_cv);
		      }
		  }
	    }
	    curves.push_back(loops[j][i]);
	}
	boundary_loops_.push_back(shared_ptr<CurveLoop>(new CurveLoop(curves, 
								      space_epsilons[j],
								      fix_trim_cvs)));
    }

    if (fix_trim_cvs)
    {
    	// Parameter curves may be placed on the wrong side of the seam
    	// of closed surfaces. This cannot be distinguished locally during
    	// creation. Make a check and repair if necessary
    	(void)checkParCrvsAtSeam();
    }
    
    if (fix_trim_cvs)
    {
      // We then analyze the loops and set valid_state_.
      analyzeLoops();
    }
}

//===========================================================================
BoundedSurface::
BoundedSurface(shared_ptr<ParamSurface> surf,
	       double space_epsilon)
  : ParamSurface(), iso_trim_(false), iso_trim_tol_(-1.0), valid_state_(0)
//===========================================================================
{
  shared_ptr<BoundedSurface> bd_sf = 
    dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);
  if (bd_sf.get())
    {
      // Already a bounded surface. Copy content
      surface_ = bd_sf->surface_;
      boundary_loops_ = bd_sf->boundary_loops_;
      loop_fixed_ = bd_sf->loop_fixed_;
      iso_trim_= bd_sf->iso_trim_;
      iso_trim_tol_ = bd_sf->iso_trim_tol_;
      valid_state_ = bd_sf->valid_state_;
    }
  else
    {
      surface_ = surf;
      vector<CurveLoop> loops = SurfaceTools::allBoundarySfLoops(surf, space_epsilon);
      for (size_t ki=0; ki<loops.size(); ++ki)
	{
	  shared_ptr<CurveLoop> curr_loop =
	    shared_ptr<CurveLoop>(new CurveLoop(loops[ki].getCurves(),
						loops[ki].getSpaceEpsilon()));
	    boundary_loops_.push_back(curr_loop);
	    loop_fixed_.push_back(0);
	}
      iso_trim_ = true;
      iso_trim_tol_ = space_epsilon;

      // We then analyze the loops and set valid_state_.
      analyzeLoops();
    }
}

//===========================================================================
BoundedSurface::
BoundedSurface(shared_ptr<ParamSurface> surf,
	       std::vector<CurveLoop>& loops)
  : ParamSurface(), surface_(surf), iso_trim_(false), iso_trim_tol_(-1.0), valid_state_(0)
//===========================================================================
{
  for (size_t ki=0; ki<loops.size(); ++ki)
    {
      shared_ptr<CurveLoop> curr_loop =
	shared_ptr<CurveLoop>(new CurveLoop(loops[ki].getCurves(),
					    loops[ki].getSpaceEpsilon()));
      boundary_loops_.push_back(curr_loop);
      loop_fixed_.push_back(0);
    }
    
    // We then analyze the loops and set valid_state_.
    analyzeLoops();
}

//===========================================================================
BoundedSurface::
BoundedSurface(shared_ptr<ParamSurface> surf,
	       std::vector<shared_ptr<CurveLoop> >& loops)
  : ParamSurface(), surface_(surf), iso_trim_(false), iso_trim_tol_(-1.0), valid_state_(0)
//===========================================================================
{
  for (size_t ki=0; ki<loops.size(); ++ki)
    {
      boundary_loops_.push_back(loops[ki]);
      loop_fixed_.push_back(0);
    }
    
    // We then analyze the loops and set valid_state_.
    analyzeLoops();
}

 //===========================================================================
BoundedSurface::~BoundedSurface()
//===========================================================================
{
}

//===========================================================================
void BoundedSurface::read(std::istream& is)
//===========================================================================
{
    bool fix_trim_cvs = false;
    read(is, fix_trim_cvs);
}

//===========================================================================
void BoundedSurface::read(std::istream& is,
			  bool fix_trim_cvs)
//===========================================================================
{
    // We verify that the object is valid.
    bool is_good = is.good();
    if (!is_good) {
	THROW("Invalid geometry file!");
    }
    ALWAYS_ERROR_IF(!boundary_loops_.empty(),
		    "This surface already exists");
    ALWAYS_ERROR_IF(surface_.get()!=NULL,
		    "This surface already exists");


    int instance_type;
    is >> instance_type;
    ClassType type = ClassType(instance_type); // Needs this conversion

    shared_ptr<GeomObject> goobject(Factory::createObject(type));
    shared_ptr<ParamSurface> tmp_srf 
	= dynamic_pointer_cast<ParamSurface, GeomObject>(goobject);
    ALWAYS_ERROR_IF(tmp_srf.get() == 0,
		    "Can not read this instance type");

    try
    {
	tmp_srf->read(is);
	surface_ = tmp_srf;
    }
    catch (...)
    { // We want the read routine to continue reading data, not a good strategy to throw before all object data is parsed.
	LOG_INFO("Failed reading the surface.");
    }

    int no_boundary_loops;
    is >> no_boundary_loops;
    // We verify that the object is valid.
    is_good = is.good();
    if (!is_good) {
	THROW("Invalid geometry file!");
    }
    for (int i=0; i<no_boundary_loops; ++i) {
	vector< shared_ptr<CurveOnSurface> > curves;
	int boundary_loops_i_size;
	double space_epsilon;
	is >> boundary_loops_i_size;
	is >> space_epsilon;
	is_good = is.good();
	if (!is_good) {
	    THROW("Invalid geometry file!");
	}
	for (int j=0; j<boundary_loops_i_size; ++j) {
	    shared_ptr<CurveOnSurface> curve(new CurveOnSurface);
	    curve->setUnderlyingSurface(surface_);
	    curve->read(is);

	    // // TEST
	    // if (curve->spaceCurve().get())
	    //   curve->setParPref(false);

	    // Try to generate the parameter curve if it does not
	    // exist already
	    if (fix_trim_cvs)
	    {
		(void)curve->ensureParCrvExistence(space_epsilon);
	    }
	    curves.push_back(curve);
	}

	#ifdef CHECK_PARAM_LOOP_ORIENTATION
	// We check direction of loop.
	bool ccw = (i == 0) ? true : false;
// 	double pareps = space_epsilon*1e-6;
 	double int_tol = 1e-6;
	if (ccw != LoopUtils::paramIsCCW(curves, int_tol)) {
	    if (i == 0) {
		LOG_INFO("Outer loop not CCW in the parameter plane.");
		//THROW("Outer loop not CCW in the parameter plane.");
	    } else {
		LOG_INFO("Inner loop not CW in the parameter plane.");
		//THROW("Inner loop not CW in the parameter plane.");
	    }
	}
	#endif

	vector<shared_ptr<ParamCurve> > dummy_vec;
	for (size_t j = 0; j < curves.size(); ++j)
	   dummy_vec.push_back(curves[j]);
	shared_ptr<CurveLoop>
	   loop(new CurveLoop(dummy_vec, space_epsilon));    // will check input
	boundary_loops_.push_back(loop);
    }

    iso_trim_ = false;
    iso_trim_tol_ = -1.0;
    valid_state_ = 0;

    // Parameter curves may be placed on the wrong side of the seam
    // of closed surfaces. This cannot be distinguished locally during
    // creation. Make a check and repair if necessary
    if (fix_trim_cvs)
    {
	(void)checkParCrvsAtSeam();
    }
    
    // TESTING
    analyzeLoops();
    // Do we need this? @jbt
 //   is_good = is.good();
 //   if (!is_good) {
	//THROW("Invalid geometry file!");
 //   }
}


//===========================================================================
void BoundedSurface::write(std::ostream& os) const
//===========================================================================
{
    streamsize prev = os.precision(15);

    os << surface_->instanceType() << endl;
    surface_->write(os);
    os << endl
        << boundary_loops_.size() << endl;
    for (size_t i=0; i<boundary_loops_.size(); ++i) {
        os << boundary_loops_[i]->size() << ' ';
        os << boundary_loops_[i]->getSpaceEpsilon() << endl;
        for (int j=0; j<boundary_loops_[i]->size(); ++j)
        {
            (*boundary_loops_[i])[j]->write(os);
            os << endl;
        }
    }
    os.precision(prev);   // Reset precision to it's previous value

}

//===========================================================================
BoundedSurface* BoundedSurface::clone() const
//===========================================================================
{
  shared_ptr<ParamSurface> surf(surface_->clone());
  vector<shared_ptr<CurveLoop> > loops;
  for (size_t ki=0; ki<boundary_loops_.size(); ++ki)
    {
      double eps = boundary_loops_[ki]->getSpaceEpsilon();
      vector<shared_ptr<ParamCurve> > crvs;
      int nmb_crvs = boundary_loops_[ki]->size();
      for (int kj=0; kj<nmb_crvs; ++kj)
	{
	  shared_ptr<ParamCurve> cv((*boundary_loops_[ki])[kj]->clone());
	  shared_ptr<CurveOnSurface> sf_cv = 
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv);
	  if (sf_cv.get())
	    sf_cv->setUnderlyingSurface(surf);
	  crvs.push_back(cv);
	}

      shared_ptr<CurveLoop> loop(new CurveLoop(crvs, eps));
      loops.push_back(loop);
    }
  return (new BoundedSurface(surf, loops));
}

//===========================================================================
bool BoundedSurface::checkParCrvsAtSeam()
//===========================================================================
{
  bool changed = false;

  // Check if the underlying surface is closed
  double eps = boundary_loops_[0]->getSpaceEpsilon();  // Get tolerance
  bool closed_u, closed_v;
  SurfaceTools::checkSurfaceClosed(*surface_, closed_u, closed_v, eps);
  if ((!closed_u) && (!closed_v))
    return false;

  // Check continuity of parameter curve. Only the outer loop is considered
  // First collect all CurveOnSurface curves and make sure that the two
  // representations are consistent
  int nmb = boundary_loops_[0]->size();
  vector<shared_ptr<CurveOnSurface> > cvs(nmb);
  int ki;
  for (ki=0; ki<nmb; ++ki)
    {
      cvs[ki] = 
	dynamic_pointer_cast<CurveOnSurface,ParamCurve>((*boundary_loops_[0])[ki]);
      if (!cvs[ki]->parameterCurve().get())
	return false;   // Missing parameter curve. Skip fix
      if (cvs[ki].get())
	cvs[ki]->makeCurvesConsistent(cvs[ki]->parPref());
    }

  RectDomain dom = surface_->containingDomain();
  double umin = dom.umin();
  double umax = dom.umax();
  double vmin = dom.vmin();
  double vmax = dom.vmax();
  double per_u = umax - umin;
  double per_v = vmax - vmin;
  
  bool b, r, t, l;
  surface_->isDegenerate(b, r, t, l, eps);

  int kj, kr=nmb-1;
  for (ki=0; ki<nmb; ++ki, kr++)
    {
      kj = (ki+1)%nmb;
      kr = kr%nmb;

      Point pt1 = cvs[kr]->parameterCurve()->point(cvs[kr]->endparam());
      Point pt2 = cvs[ki]->parameterCurve()->point(cvs[ki]->startparam());
      Point pt3 = cvs[ki]->parameterCurve()->point(cvs[ki]->endparam());
      Point pt4 = cvs[kj]->parameterCurve()->point(cvs[kj]->startparam());

      if (pt1.dist(pt2) < eps && pt3.dist(pt4) < eps)
	continue;  // Continuity OK

      // if (closed_u && fabs(fabs(pt1[0]-pt2[0])-per_u) < eps &&
      // 	  fabs(fabs(pt3[0]-pt4[0])-per_u) < eps)
      if (closed_u && fabs(fabs(pt1[0]-pt2[0])-per_u) < eps &&
	  fabs(pt2[0]-pt3[0]) < eps &&
	  (fabs(umax-pt3[0]) < eps || fabs(pt3[0]-umin) < eps))
	{
	  // Translate parameter curve of curve ki
	  Point vec(pt1[0]-pt2[0], 0.0);
	  cvs[ki]->translateParameterCurve(vec);
	  changed = true;
	}
      
      else if (closed_u && 
          ((fabs(pt1[1]-vmin) < eps && fabs(pt2[1]-vmin) < eps && b) ||
           (fabs(pt1[1]-vmax) < eps && fabs(pt2[1]-vmax) < eps && t)) &&
	  fabs(pt2[0]-pt3[0]) < eps &&
	  (fabs(umax-pt3[0]) < eps || fabs(pt3[0]-umin) < eps))
	{
	  // Translate parameter curve of curve ki
	  Point vec(pt4[0]-pt3[0], 0.0);
	  cvs[ki]->translateParameterCurve(vec);
	  changed = true;
	}
      
      // else if (closed_v && fabs(fabs(pt1[1]-pt2[1])-per_v) < eps &&
      // 	  fabs(fabs(pt3[1]-pt4[1])-per_v) < eps)
      else if (closed_v && fabs(fabs(pt1[1]-pt2[1])-per_v) < eps &&
	       fabs(pt2[1]-pt3[1]) < eps &&
	       (fabs(vmax-pt3[1]) < eps || fabs(pt3[1]-vmin) < eps))
	{
	  // Translate parameter curve of curve ki
	  Point vec(0.0, pt1[1]-pt2[1]);
	  cvs[ki]->translateParameterCurve(vec);
	  changed = true;
	}
    }
      
  // Check degeneracy
  bool bd[4];  // left, right, bottom, top
  bool deg = surface_->isDegenerate(bd[2], bd[1], bd[3], bd[0], eps);

  if (deg)
    {
      // Check for missing trimming curves
      for (ki=0; ki<(int)cvs.size(); ++ki)
	{
	  kj = (ki+1)%((int)cvs.size());
	  Point pt1 = cvs[ki]->ParamCurve::point(cvs[ki]->endparam());
	  Point pt2 = cvs[kj]->ParamCurve::point(cvs[kj]->startparam());
	  Point par1 = cvs[ki]->parameterCurve()->point(cvs[ki]->endparam());
	  Point par2 = cvs[kj]->parameterCurve()->point(cvs[kj]->startparam());
	  if (pt1.dist(pt2) < eps && par1.dist(par2) > eps)
	    {
	      // Gap in the parameter curve loop. Check if it coincides with
	      // a degenerate edge
	      int idx = -1;
	      if (fabs(par1[0]-par2[0])<eps && fabs(par1[0]-umin) < eps)
		idx = 0;
	      else if (fabs(par1[0]-par2[0])<eps && fabs(umax-par1[0]) < eps)
		idx = 1;
	      else if (fabs(par1[1]-par2[1])<eps && fabs(par1[1]-vmin) < eps)
		idx = 2;
	      else if (fabs(par1[1]-par2[1])<eps && fabs(vmax-par1[1]) < eps)
		idx = 3;
	      if (idx >= 0 && bd[idx])
		{
		  // Add degenerate edge
		  shared_ptr<SplineCurve> par_cv(new SplineCurve(par1, par2));
		  shared_ptr<SplineCurve> space_cv(new SplineCurve(pt1, 
								   par_cv->startparam(),
								   pt2,
								   par_cv->endparam()));
		  shared_ptr<CurveOnSurface> tmp_cv(new CurveOnSurface(surface_,
								       par_cv,
								       space_cv,
								       true));
		  cvs.insert(cvs.begin()+kj, tmp_cv);
		  changed = true;
		}
	    }
	}
     
      if ((int)cvs.size() > nmb)
	{
	  // New curves are added. Update curve loop
	  vector<shared_ptr<ParamCurve> > tmp_loop_cvs(cvs.begin(), cvs.end());
	  shared_ptr<CurveLoop> tmp_loop(new CurveLoop(tmp_loop_cvs, eps));
	  boundary_loops_[0] = tmp_loop;
	}
    }
  return changed;
}

//===========================================================================
BoundingBox BoundedSurface::boundingBox() const
//===========================================================================
{
  if (box_.valid())
    return box_;

  RectDomain dom = containingDomain();
  vector<shared_ptr<ParamSurface> > sub_sfs;

  RectDomain dom2 = surface_->containingDomain();
  double tol1 = std::min(1.0e-1, 0.001*(dom.umax()-dom.umin()));
  double tol2 = std::min(1.0e-1, 0.001*(dom.vmax()-dom.vmin()));
  if (dom.umin()-dom2.umin()<tol1 && dom2.umax()-dom.umax()<tol1 &&
      dom.vmin()-dom2.vmin()<tol2 && dom2.vmax()-dom.vmax()<tol2)
    return surface_->boundingBox();
  else
    {
      double umin = std::max(dom.umin(), dom2.umin());
      double umax = std::min(dom.umax(), dom2.umax());
      double vmin = std::max(dom.vmin(), dom2.vmin());
      double vmax = std::min(dom.vmax(), dom2.vmax());
      try {
	sub_sfs = surface_->subSurfaces(umin, vmin, umax, vmax);
      }
      catch (...)
	{
	  box_ = surface_->boundingBox();
	  return box_;
	}
    }

  box_ = (sub_sfs.size() == 1) ? sub_sfs[0]->boundingBox() : 
    surface_->boundingBox();
  return box_;
}


//===========================================================================
DirectionCone BoundedSurface::normalCone() const
//===========================================================================
{
  RectDomain dom = containingDomain();
  vector<shared_ptr<ParamSurface> > sub_sfs;
  try {
    sub_sfs = surface_->subSurfaces(dom.umin(), dom.vmin(), 
				    dom.umax(), dom.vmax());
  }
  catch (...)
    {
      return surface_->normalCone();
    }

  return (sub_sfs.size() == 1) ? sub_sfs[0]->normalCone() : 
    surface_->normalCone();
}


//===========================================================================
DirectionCone BoundedSurface::tangentCone(bool pardir_is_u) const
//===========================================================================
{
  RectDomain dom = containingDomain();
  vector<shared_ptr<ParamSurface> > sub_sfs;
  try {
    sub_sfs = surface_->subSurfaces(dom.umin(), dom.vmin(), 
				    dom.umax(), dom.vmax());
  }
  catch (...)
    {
      return surface_->tangentCone(pardir_is_u);
    }

  return (sub_sfs.size() == 1) ? sub_sfs[0]->tangentCone(pardir_is_u) : 
    surface_->tangentCone(pardir_is_u);
}


//===========================================================================
int BoundedSurface::dimension() const
//===========================================================================
{
    return surface_->dimension();
}


//===========================================================================
ClassType BoundedSurface::instanceType() const
//===========================================================================
{
    return classType();
}

//===========================================================================
const CurveBoundedDomain& BoundedSurface::parameterDomain() const
//===========================================================================
{
  domain_ = CurveBoundedDomain(boundary_loops_);
  return domain_;
}


//===========================================================================
RectDomain BoundedSurface::containingDomain() const
//===========================================================================
{
   RectDomain dom1 = parameterDomain().containingDomain();
   RectDomain dom2 = surface_->containingDomain();
   dom1.intersectWith(dom2);
   return dom1;
}


//===========================================================================
bool BoundedSurface::inDomain(double u, double v, double eps) const 
//===========================================================================
{
    Array<double, 2> pnt(u, v);
    return parameterDomain().isInDomain(pnt, eps);
}

//===========================================================================
int BoundedSurface::inDomain2(double u, double v, double eps) const 
//===========================================================================
{
    Array<double, 2> pnt(u, v);
    return parameterDomain().isInDomain2(pnt, eps);
}

//===========================================================================
bool BoundedSurface::onBoundary(double u, double v, double eps) const 
//===========================================================================
{
    Array<double, 2> pnt(u, v);
    return parameterDomain().isOnBoundary(pnt, eps);
}

//===========================================================================
Point BoundedSurface::closestInDomain(double u, double v) const 
//===========================================================================
{
    Array<double, 2> pnt(u, v);
    Array<double, 2> close(0.0, 0.0);
    double eps = 1.0e-8;  // A small number
    parameterDomain().closestInDomain(pnt, close, eps);
    return Point(close.x(), close.y());
}

//===========================================================================
CurveLoop BoundedSurface::outerBoundaryLoop(double degenerate_epsilon) const
//===========================================================================
{
    // As the outer boundary loop has been moved up front, we return the first loop.
    std::vector<shared_ptr<ParamCurve> > curves;
    CurveLoop& loop = *(boundary_loops_[0]);
    for (int i = 0; i < loop.size(); ++i) {
	if (!loop[i]->isDegenerate(degenerate_epsilon))
	    curves.push_back(loop[i]);
    }
    return CurveLoop(curves, loop.getSpaceEpsilon());
}


//===========================================================================
std::vector<CurveLoop> BoundedSurface::allBoundaryLoops(double degenerate_epsilon) const
//===========================================================================
{
    std::vector<CurveLoop> clvec;
    for (size_t j=0; j<boundary_loops_.size(); j++) {
	std::vector<shared_ptr<ParamCurve> > curves;
	CurveLoop& loop = *(boundary_loops_[j]);
	for (int i = 0; i < loop.size(); ++i) {
	    if (!loop[i]->isDegenerate(degenerate_epsilon))
		curves.push_back(loop[i]);
	}
	if (curves.size() > 0)
	  clvec.push_back(CurveLoop(curves, loop.getSpaceEpsilon(),
				    false));
    }
    return clvec;
}

//===========================================================================
std::vector<CurveLoop> BoundedSurface::absolutelyAllBoundaryLoops() const
//===========================================================================
{
    std::vector<CurveLoop> clvec;
    for (size_t j=0; j<boundary_loops_.size(); j++) {
	double loop_tol = boundary_loops_[j]->getSpaceEpsilon();
	std::vector<shared_ptr<ParamCurve> > curves;
	CurveLoop& loop = *(boundary_loops_[j]);
	for (int i = 0; i < loop.size(); ++i) {
	    curves.push_back(loop[i]);
	}
	clvec.push_back(CurveLoop(curves, loop_tol));
    }
    return clvec;
}

//===========================================================================
void BoundedSurface::getLoopCvInfo(int& nmb_loops, int& nmb_cvs, 
				   int& nmb_corners, double& min_cv_len, 
				   double& max_cv_len, double angtol) const
//===========================================================================
{
  nmb_loops = (int)boundary_loops_.size();
  nmb_cvs = 0;
  nmb_corners = 0;
  min_cv_len = std::numeric_limits<double>::max();
  max_cv_len = 0.0;
  for (int ki=0; ki<nmb_loops; ++ki)
    {
      int nmb = boundary_loops_[ki]->size();
      nmb_cvs += nmb;
      vector<Point> prev(2);
      shared_ptr<ParamCurve> cv0 = (*boundary_loops_[ki])[nmb-1];
      cv0->point(prev, cv0->endparam(), 1);
      for (int kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ParamCurve> cv = (*boundary_loops_[ki])[kj];
	  double len = cv->estimatedCurveLength();
	  min_cv_len = std::min(min_cv_len, len);
	  max_cv_len = std::max(max_cv_len, len);
	  vector<Point> curr(2);
	  cv->point(curr, cv->startparam(), 1);
	  double ang = prev[1].angle(curr[1]);
	  if (ang > angtol)
	    nmb_corners++;
	  cv->point(prev, cv->endparam(), 1);
	}
    }
}

//===========================================================================
void BoundedSurface::point(Point& pt, double upar, double vpar) const
//===========================================================================
{
    surface_->point(pt, upar, vpar);
}


//===========================================================================
void BoundedSurface::point(std::vector<Point>& pts, 
			     double upar, double vpar,
			     int derivs, bool u_from_right,
			     bool v_from_right, double resolution) const
//===========================================================================
{
    surface_->point(pts, upar, vpar, derivs, u_from_right, v_from_right,
		    resolution);
}


// //===========================================================================
// void BoundedSurface::point(std::vector<Point>& pts, 
// 			     double upar, double vpar,
// 			     int derivs) const
// //===========================================================================
// {
//     surface_->point(pts, upar, vpar, derivs);
// }


//===========================================================================
void BoundedSurface::evalGrid(int num_u, int num_v, 
			      double umin, double umax, 
			      double vmin, double vmax,
			      std::vector<double>& points,
			      double nodata_val) const
//===========================================================================
{
  int dim = dimension();
  double tol = 1.0e-6;  // A good tolerance for intersections
  CurveBoundedDomain dom = parameterDomain();
  
 #ifdef DEBUG
    std::ofstream of("tmp_grid_bd.g2");
    (void)of.precision(15);
    of << "400 1 0 4 255 0 0 255" << std::endl;

    vector<double> tmppt;
#endif

 // Evaluate underlying surface in grid.
  // This is done to be able to utilize structures for grid
  // evaluation and improve performance
  surface_->evalGrid(num_u, num_v, umin, umax, vmin, vmax,
		     points, nodata_val);

  // Modify the value of points lying outside the bounded surface
  int ki, kj, kr, kh;
  double udel = (umax - umin)/(double)(num_u-1);
  double vdel = (vmax - vmin)/(double)(num_v-1);
  double upar, vpar;
  double *pos;
  for (ki=0, vpar=vmin, pos=&points[0]; ki<num_v; ++ki, vpar+=vdel)
    {
      // Make horizontal parameter curve
      SplineCurve cv(Point(umin-2*udel,vpar), umin-2*udel, 
		     Point(umax+2*udel,vpar), umax+2*udel);

      // Find inside intervals
      vector<double> par_intervals;
      dom.findPcurveInsideSegments(cv, tol, par_intervals);

      if (par_intervals.size() == 0)
	{
	  // No boundary found. Extend with appropriate parameter bound
	  Vector2D param(0.5*(umin+umax), vpar);
	  if (dom.isInDomain(param, tol))
	    {
	      par_intervals.push_back(umin);
	      par_intervals.push_back(umax);
	    }
	}
      if (par_intervals.size() % 2 == 1)
	{
	  // Boundary touch. Extend with appropriate parameter bound
	  Vector2D param(0.5*(umin+par_intervals[0]), vpar);
	  if (dom.isInDomain(param, tol))
	    par_intervals.insert(par_intervals.begin(), umin);
	  else
	    par_intervals.push_back(umax);
	}

      for (kj=0, kr=0, upar=umin; kj<num_u; ++kj, upar+=udel, pos+=dim)
	{
	  // Check if the point is inside the trimmed surface
	  for(; kr<(int)par_intervals.size(); kr+=2)
	    if (upar <= par_intervals[kr+1])
	      break;
	  if (!(kr<(int)par_intervals.size() && 
		upar>=par_intervals[kr] && upar<=par_intervals[kr+1]))
	    {
	      for (kh=0; kh<dim; ++kh)
		pos[kh] = nodata_val;
	    }
#ifdef DEBUG
	  else
	    {
	      tmppt.push_back(upar);
	      tmppt.push_back(vpar);
	      tmppt.push_back(pos[0]);
	    }
#endif
	}
    }
#ifdef DEBUG
      of << tmppt.size()/3 << std::endl;
      for (size_t ka=0; ka<tmppt.size(); ka+=3)
	of << tmppt[ka] << " " << tmppt[ka+1] << " " << tmppt[ka+2] << std::endl;
#endif
}
  
//===========================================================================
Point BoundedSurface::getInternalPoint(double& upar, double& vpar) const
//===========================================================================
{
  CurveBoundedDomain dom = parameterDomain();
  dom.getInternalPoint(upar, vpar);
  return ParamSurface::point(upar, vpar);
}

//===========================================================================
void BoundedSurface::normal(Point& n, double upar, double vpar) const
//===========================================================================
{
    surface_->normal(n, upar, vpar);
}


//===========================================================================
vector<shared_ptr<ParamCurve> >
BoundedSurface::constParamCurves(double parameter, bool pardir_is_u) const
//===========================================================================
{
    // We first create a parametric cv to intersect with the sf.
    double epsgeo = getEpsGeo();
    // We are interested in the parametrization wrt the surface, hence we use the underlying surface.
    // Additionally using *this instead of *surface_ is a bad idea as getParEpsilon() will trigger a call
    // to this->constParamCurves().
    Point pareps = SurfaceTools::getParEpsilon(*surface_, epsgeo);
    double eps = 0.5*(pareps[0]+pareps[1]);
    eps = std::max(1.0e-6, eps);

    vector<shared_ptr<ParamCurve> > dummy;

    RectDomain rect_dom = containingDomain();
    Point from_pt(2);
    Point to_pt(2);
    double startpar, endpar;
    if (!pardir_is_u) {
	from_pt[0] = to_pt[0] = parameter;
	startpar = from_pt[1] = rect_dom.vmin();
	endpar = to_pt[1] = rect_dom.vmax();
    } else {
	startpar = from_pt[0] = rect_dom.umin();
	endpar = to_pt[0] = rect_dom.umax();
	from_pt[1] = to_pt[1] = parameter;
    }
    if (endpar - startpar < eps)
      return dummy;

    shared_ptr<SplineCurve> par_iso_cv(new SplineCurve(from_pt, startpar,
						       to_pt, endpar));
    shared_ptr<CurveOnSurface> iso_cv_on_sf
	(new CurveOnSurface(surface_, par_iso_cv, true));

    // We now intersect with the trimmed sf.
    // @@sbr Suspecting that this may cause trouble for intersection analysis ...
    BoundedSurface bd_sf(*this); // @@sbr To keep function const ...
    std::vector<shared_ptr<CurveOnSurface> > int_cvs =
	BoundedUtils::intersectWithSurface(*iso_cv_on_sf,
					     bd_sf, eps);

    vector<shared_ptr<ParamCurve> > return_cvs(int_cvs.begin(), int_cvs.end());

    return return_cvs;
}


//===========================================================================
vector<shared_ptr<ParamSurface> >
BoundedSurface::subSurfaces(double from_upar,
			      double from_vpar,
			      double to_upar,
			      double to_vpar,
			      double fuzzy) const
//===========================================================================
{
    if (surface_->instanceType() == Class_SplineSurface) {
	// If boundaries are close to existing knots, we snap.
	shared_ptr<SplineSurface> spline_sf =
	    dynamic_pointer_cast<SplineSurface, ParamSurface>(surface_);
	const BsplineBasis& basis_u = spline_sf->basis_u();
	const BsplineBasis& basis_v = spline_sf->basis_v();
	basis_u.knotIntervalFuzzy(from_upar, fuzzy);
	basis_u.knotIntervalFuzzy(to_upar, fuzzy);
	basis_v.knotIntervalFuzzy(from_vpar, fuzzy);
	basis_v.knotIntervalFuzzy(to_vpar, fuzzy);
	// @@sbr Suppose the fuzzy value could be used more.
    }

    double min_eps = boundary_loops_[0]->getSpaceEpsilon(); //= 0.0001; // Just a guess
    for (size_t ki = 1; ki < boundary_loops_.size(); ++ki) {
	double eps = boundary_loops_[ki]->getSpaceEpsilon();
	if (eps < min_eps)
	    min_eps = eps;
    }

    std::ofstream out_file("tmp0_mini_surf.g2");

    // First fetch the surrounding domain of the current parameter domain
    RectDomain domain = containingDomain();
    //   // Fetch underlying surface
      /*
    shared_ptr<SplineSurface> under_sf = 
	dynamic_pointer_cast<SplineSurface, ParamSurface>(surface_);
      */
    shared_ptr<ParamSurface> under_sf =
	surface_;
    if (under_sf.get() == NULL)
      THROW("did not expect this!");

    // Make a copy of the current surface
    vector<shared_ptr<BoundedSurface> > sub_sfs;
    sub_sfs.push_back(shared_ptr<BoundedSurface>(clone()));

#ifdef SBR_DBG
    size_t kkr;
    for (kkr=0; kkr<sub_sfs.size(); kkr++)
    {
	sub_sfs[kkr]->writeStandardHeader(out_file);
	sub_sfs[kkr]->write(out_file);
    }
#endif

    vector<shared_ptr<BoundedSurface> > new_sub_sfs;
    // Make boundary loops to the subsurface
    if (from_upar > domain.umin()) {
	for (size_t ki = 0; ki < sub_sfs.size(); ++ki) {
	    // First make new trimming curves
	    vector<shared_ptr<CurveOnSurface> > trim_crvs;
	    CurveBoundedDomain curve_domain = sub_sfs[ki]->parameterDomain();
	    curve_domain.clipWithDomain(2, from_upar, min_eps, 
					sub_sfs[ki]->underlyingSurface(), trim_crvs);
	    if (trim_crvs.size() == 0)
	    {
		continue;
	    }
	    for (size_t kj = 0; kj < trim_crvs.size(); ++kj) // Curve(s) picked in pos direction.
		trim_crvs[kj]->reverseParameterDirection();
	    vector<vector<shared_ptr<CurveOnSurface> > > loop_curves;
	    try {
		loop_curves = BoundedUtils::getBoundaryLoops(*sub_sfs[ki], 
							     trim_crvs, min_eps);
	    } catch (...) {
		//THROW("Failed extracting boundary loops.");
		continue;
	    }

	    vector<shared_ptr<BoundedSurface> > local_new_sub_sfs =
	      BoundedUtils::createTrimmedSurfs(loop_curves, 
					       sub_sfs[ki]->underlyingSurface(), //surface_, 
					       min_eps);
	    new_sub_sfs.insert(new_sub_sfs.end(), local_new_sub_sfs.begin(), local_new_sub_sfs.end());
	}
	sub_sfs = new_sub_sfs;
	new_sub_sfs.clear();
    }

#ifdef SBR_DBG
    for (kkr=0; kkr<sub_sfs.size(); kkr++)
    {
	sub_sfs[kkr]->writeStandardHeader(out_file);
	sub_sfs[kkr]->write(out_file);
    }
#endif

    if (to_upar < domain.umax()) {
	for (size_t ki = 0; ki < sub_sfs.size(); ++ki) {
	    // First make new trimming curves
	    vector<shared_ptr<CurveOnSurface> > trim_crvs;
	    CurveBoundedDomain curve_domain = sub_sfs[ki]->parameterDomain();
	    curve_domain.clipWithDomain(2, to_upar, min_eps, 
					sub_sfs[ki]->underlyingSurface(), trim_crvs);
	    if (trim_crvs.size() == 0)
	    {
		continue;
	    }
	    //      for (kj = 0; kj < trim_crvs.size(); ++kj) // Curve(s) picked in pos direction.
	    //	  trim_crvs[kj]->reverseParameterDirection();
	    vector<vector<shared_ptr<CurveOnSurface> > > loop_curves;
	    try {
		loop_curves = BoundedUtils::getBoundaryLoops(*sub_sfs[ki], 
							     trim_crvs, min_eps);
	    } catch (...) {
		//THROW("Failed extracting boundary loops.");
		continue;
	    }

	    vector<shared_ptr<BoundedSurface> > local_new_sub_sfs =
		BoundedUtils::createTrimmedSurfs(loop_curves, 
						 sub_sfs[ki]->underlyingSurface(),
	      min_eps);
	    new_sub_sfs.insert(new_sub_sfs.end(), local_new_sub_sfs.begin(), 
			       local_new_sub_sfs.end());
	}
	sub_sfs = new_sub_sfs;
	new_sub_sfs.clear();
    }

#ifdef SBR_DBG
    for (kkr=0; kkr<sub_sfs.size(); kkr++)
    {
	sub_sfs[kkr]->writeStandardHeader(out_file);
	sub_sfs[kkr]->write(out_file);
    }
#endif

    if (from_vpar > domain.vmin()) {
	for (size_t ki = 0; ki < sub_sfs.size(); ++ki) {
	    // First make new trimming curves
	    vector<shared_ptr<CurveOnSurface> > trim_crvs;
	    CurveBoundedDomain curve_domain = sub_sfs[ki]->parameterDomain();
	    curve_domain.clipWithDomain(1, from_vpar, min_eps, 
					sub_sfs[ki]->underlyingSurface(), trim_crvs);
	    if (trim_crvs.size() == 0)
	    {
		continue;
	    }
// 	    for (size_t kj = 0; kj < trim_crvs.size(); ++kj) // Curve(s) picked in pos direction.
// 		trim_crvs[kj]->reverseParameterDirection();
	    vector<vector<shared_ptr<CurveOnSurface> > > loop_curves;
	    try {
		loop_curves = BoundedUtils::getBoundaryLoops(*sub_sfs[ki], 
							     trim_crvs, min_eps);
	    } catch (...) {
		//THROW("Failed extracting boundary loops.");
		continue;
	    }

	    vector<shared_ptr<BoundedSurface> > local_new_sub_sfs =
		BoundedUtils::createTrimmedSurfs(loop_curves, 
						 sub_sfs[ki]->underlyingSurface(),min_eps);
	    new_sub_sfs.insert(new_sub_sfs.end(), local_new_sub_sfs.begin(), local_new_sub_sfs.end());
	}
	sub_sfs = new_sub_sfs;
	new_sub_sfs.clear();
    }

#ifdef SBR_DBG
    for (kkr=0; kkr<sub_sfs.size(); kkr++)
    {
	sub_sfs[kkr]->writeStandardHeader(out_file);
	sub_sfs[kkr]->write(out_file);
    }
#endif

    if (to_vpar < domain.vmax()) {
	for (size_t ki = 0; ki < sub_sfs.size(); ++ki) {
	    // First make new trimming curves
	    vector<shared_ptr<CurveOnSurface> > trim_crvs;
	    CurveBoundedDomain curve_domain = sub_sfs[ki]->parameterDomain();
	    curve_domain.clipWithDomain(1, to_vpar, min_eps, 
					sub_sfs[ki]->underlyingSurface(), trim_crvs);
	    if (trim_crvs.size() == 0)
	    {
		continue;
	    }
	    for (size_t kj = 0; kj < trim_crvs.size(); ++kj) // Curve(s) picked in pos direction.
		trim_crvs[kj]->reverseParameterDirection();
	    vector<vector<shared_ptr<CurveOnSurface> > > loop_curves;
	    try {
		loop_curves = BoundedUtils::getBoundaryLoops(*sub_sfs[ki], 
							     trim_crvs, min_eps);
	    } catch (...) {
		//THROW("Failed extracting boundary loops.");
		continue;
	    }

	    vector<shared_ptr<BoundedSurface> > local_new_sub_sfs =
	      BoundedUtils::createTrimmedSurfs(loop_curves, sub_sfs[ki]->underlyingSurface(), min_eps);
	    new_sub_sfs.insert(new_sub_sfs.end(), local_new_sub_sfs.begin(), local_new_sub_sfs.end());
	}
	sub_sfs = new_sub_sfs;
    }

#ifdef SBR_DBG
    for (kkr=0; kkr<sub_sfs.size(); kkr++)
    {
	sub_sfs[kkr]->writeStandardHeader(out_file);
	sub_sfs[kkr]->write(out_file);
    }
#endif

    vector<shared_ptr<ParamSurface> > return_sfs(sub_sfs.begin(), sub_sfs.end());

    return return_sfs;
}


//===========================================================================
double BoundedSurface::nextSegmentVal(int dir, double par, bool forward, double tol) const
//===========================================================================
{ // @@sbr Possibly check that par is inside domain?
  return surface_->nextSegmentVal(dir, par, forward, tol);
}


//===========================================================================
void BoundedSurface::closestPoint(const Point& pt,
				    double&  clo_u,
				    double&  clo_v, 
				    Point& clo_pt,
				    double&  clo_dist,
				    double   epsilon,
                                    const RectDomain* domain_of_interest,
				    // domain_of_interest is never used
				    double   *seed) const
//===========================================================================
{
    const CurveBoundedDomain& dom = parameterDomain();

    Vector2D new_seed_vec;
    double *new_seed = seed;
    double domain_tol = boundary_loops_[0]->getSpaceEpsilon();
    domain_tol = std::max(domain_tol, 1.0e-7);
    if (seed) {
        new_seed_vec[0] = seed[0];
        new_seed_vec[1] = seed[1];
	// Check that the seed is inside the domain, if not we
	// find the closest in domain as a new seed.
	Vector2D old_seed(seed);
	bool in_domain = false;
	try {
	  in_domain = dom.isInDomain(old_seed, domain_tol);
	    }
	catch (...)
	  {
	    LOG_INFO("Failed deciding whether point was inside domain. ");
	    in_domain = true; // Don't mess more about it
	  }

	if (!in_domain) {
	  dom.closestInDomain(old_seed, new_seed_vec, epsilon);
	}
	new_seed = new_seed_vec.begin();
    }

    RectDomain d2 = dom.containingDomain();
    surface_->closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, epsilon, &d2,
			   new_seed);
    try {
	if (dom.isInDomain(Vector2D(clo_u, clo_v), domain_tol)) {
	    return;
	}
    } catch (...) {
	LOG_INFO("Failed deciding whether point was inside domain. "
		   "At least close, let's say we found it!");
	return;
    }

    // It was not in the required domain. We find the closest boundary point.
    // We try two approaches, both in parameter space and in the real space.
    Vector2D new_clo_vec;
    try {
	dom.closestInDomain(Vector2D(clo_u, clo_v), new_clo_vec, domain_tol);
    } catch (...) {
	THROW("Failed finding closest point in domain.");
    }
    Point new_clo_pt = surface_->point(new_clo_vec[0], new_clo_vec[1]);
    closestBoundaryPoint(pt, clo_u, clo_v, clo_pt, clo_dist,
			 epsilon, domain_of_interest, seed);

    if (pt.dist(new_clo_pt) < pt.dist(clo_pt)) {
	clo_pt = new_clo_pt;
	clo_dist = pt.dist(clo_pt);
	clo_u = new_clo_vec[0];
	clo_v = new_clo_vec[1];
    }
}


//===========================================================================
void BoundedSurface::closestBoundaryPoint(const Point& pt,
					    double&        clo_u,
					    double&        clo_v, 
					    Point&       clo_pt,
					    double&        clo_dist,
					    double epsilon,
					    const RectDomain* domain_of_interest,
					    double *seed) const
//===========================================================================
{
    // Find closest point on every subcurve of spatial boundaryloop and compare.
    double clo_t;
    (*(boundary_loops_[0]))[0]->closestPoint(pt, clo_t, clo_pt, clo_dist);
    Point tmp_clopt(dimension());
    double tmp_clot;
    double tmp_cld;
    int n1 = (int)boundary_loops_.size();
    int i, j;
    vector<double> new_seed(2);
    for (j=0; j<n1; j++) {
	int n2 = boundary_loops_[j]->size();
	for (i = 0; i < n2; ++i) { // As we may try to create seed we iterate from 0.
	    shared_ptr<ParamCurve> bd_cv = (*boundary_loops_[j])[i];
	    double clo_par_t;
	    Point clo_par_pt;
	    double* cv_seed = NULL;
	    if ((seed != 0) && (bd_cv->instanceType() == Class_CurveOnSurface)) {
		shared_ptr<CurveOnSurface> cv_on_sf =
		    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(bd_cv);
		if ((cv_on_sf->parameterCurve().get() != 0) && (cv_on_sf->parPref())) {

		    Point par_pt(seed[0], seed[1]);
		    double clo_par_dist;
		    cv_on_sf->parameterCurve()->closestPoint(par_pt, clo_par_t, clo_par_pt, clo_par_dist);
		    cv_seed = &clo_par_t;
// 		    // We use parameter curve to construct a sensible seed.
// 		    Point par_pt = cv_on_sf->parameterCurve()->point(clo_t);
 		    new_seed[0] = par_pt[0];
 		    new_seed[1] = par_pt[1];
// 		    seed = &new_seed[0];
		}
	    }

	    bd_cv->closestPoint(pt, bd_cv->startparam(), bd_cv->endparam(),
				tmp_clot, tmp_clopt, tmp_cld, cv_seed);
	    if (tmp_cld < clo_dist) {
		clo_t = tmp_clot;
		clo_pt = tmp_clopt;
		clo_dist = tmp_cld;
		if (cv_seed != 0) {
		    seed = &new_seed[0];
		}
	    }
	}
    }

    // We need to get parameter values for clo_pt (suppose we could use seed if set in above routine...).
    double tmp_u, tmp_v;
    RectDomain local_rect_dom = parameterDomain().containingDomain();
    const RectDomain* rect_dom = (domain_of_interest) ? domain_of_interest : &local_rect_dom;
    surface_->closestPoint(clo_pt, tmp_u, tmp_v, tmp_clopt, tmp_cld,
			   epsilon, rect_dom, seed);
    // Now, the parameter point (tmp_u, tmp_v) should be in the parameter
    // domain of the surface. If so, we return happily.
    // Otherwise, we find the closest point in the domain.

    // VSK, 0902. First check if the point found on the boundary and the point in the surface
    // is the same. In that case, we are done
    if (clo_pt.dist(tmp_clopt) <= epsilon)
      {
	clo_u = tmp_u;
	clo_v = tmp_v;
	clo_pt = tmp_clopt;
	clo_dist = clo_pt.dist(pt); 
      }
    else
      {
	//	clo_dist = tmp_cld;
	const CurveBoundedDomain& dom = parameterDomain();
	bool is_in_domain = false;
	double domain_tol = std::max(epsilon, 1.0e-7);
	try {
	  // Test is rather unstable when point is on/near boundary.
	  is_in_domain = dom.isInDomain(Vector2D(tmp_u, tmp_v), domain_tol);
	} catch (...) {
	  // 	MESSAGE("Failed deciding whether point was in domain.");
	  is_in_domain = true;
	}

	if (is_in_domain) {
	  clo_u = tmp_u;
	  clo_v = tmp_v;
	  clo_pt = tmp_clopt;
	  clo_dist = clo_pt.dist(pt); //	clo_dist = tmp_cld;
	} else {
	  // 	MESSAGE_IF(domain_of_interest,
	  // 		      "Input parameter domain_of_interest may not be used!");
	  Vector2D par_pt;
	  // @afr: Again, I use (spatial) epsilon for domain comparisons...
	  dom.closestInDomain(Vector2D(tmp_u, tmp_v), par_pt, epsilon);
	  point(clo_pt, par_pt[0], par_pt[1]);
	  clo_u = par_pt[0];
	  clo_v = par_pt[1];
	  clo_dist = clo_pt.dist(pt);	
	}
      }
}


//===========================================================================
void
BoundedSurface::getBoundaryInfo(Point& pt1, Point& pt2, 
				 double epsilon, SplineCurve*& cv,
				 SplineCurve*& crosscv, double knot_tol) const
//===========================================================================
{
  vector<shared_ptr<CurveOnSurface> > bd_cvs;
  getBoundaryInfo(pt1, pt2, bd_cvs);

  crosscv = NULL;  // Don't know anything about this
  if (bd_cvs.size() == 0)
    cv = NULL;
  else
    {
      cv = bd_cvs[0]->spaceCurve()->geometryCurve();
      for (size_t ki=1; ki<bd_cvs.size(); ++ki)
	{
	  shared_ptr<SplineCurve> cv2 =  
	    shared_ptr<SplineCurve>(bd_cvs[ki]->spaceCurve()->geometryCurve());
	  double dist;
	  cv->appendCurve(cv2.get(), 0, dist, false);
	}
    }
 
}


//===========================================================================
void BoundedSurface::getBoundaryInfo(Point& pt1, Point& pt2, 
				       vector<shared_ptr<CurveOnSurface> >& bd_cvs) const
//===========================================================================
{
  double ptol = 1.0e-8;
    bd_cvs.clear();

    // We run through bd_cvs locating params for closest pts in space.
    int global_clo_ind1, global_clo_ind2;
    double global_clo_par1, global_clo_par2;
    double global_clo_dist1 = 1e10, global_clo_dist2 = 1e10; // @@sbr ...
    Point global_clo_pt1, global_clo_pt2;
    int clo_loop_ind1, clo_loop_ind2;
    for (size_t ki = 0; ki < boundary_loops_.size(); ++ki) {
	int clo_ind;
	double clo_par, clo_dist;
	Point clo_pt;
	boundary_loops_[ki]->closestPoint(pt1, clo_ind, clo_par, clo_pt, clo_dist);
	if (clo_dist < global_clo_dist1) {
	    global_clo_ind1 = clo_ind;
	    global_clo_par1 = clo_par;
	    global_clo_pt1 = clo_pt;
	    global_clo_dist1 = clo_dist;
	    clo_loop_ind1 = (int)ki;
	}
	boundary_loops_[ki]->closestPoint(pt2, clo_ind, clo_par, clo_pt, clo_dist);
	if (clo_dist < global_clo_dist2) {
	    global_clo_ind2 = clo_ind;
	    global_clo_par2 = clo_par;
	    global_clo_pt2 = clo_pt;
	    global_clo_dist2 = clo_dist;
	    clo_loop_ind2 = (int)ki;
	}
    }

    if (clo_loop_ind1 != clo_loop_ind2) {
	LOG_INFO("Input pts not member of the same loop!");
    }

    int nmb_cvs = boundary_loops_[clo_loop_ind1]->size();
    int nmb_segments = (global_clo_ind2 < global_clo_ind1) ?
	global_clo_ind2 + nmb_cvs + 1 - global_clo_ind1 : global_clo_ind2 - global_clo_ind1 + 1;
    for (int ki = 0; ki < nmb_segments; ++ki) {
	int ind = (global_clo_ind1 + ki)%nmb_cvs;
	shared_ptr<CurveOnSurface> bd_cv =
	    dynamic_pointer_cast<CurveOnSurface, ParamCurve>((*boundary_loops_[clo_loop_ind1])[ind]);
	if (bd_cv.get() == 0) {
	    LOG_INFO("Wasn't expecting this ...");
	}
	double tmin = bd_cv->startparam();
	double tmax = bd_cv->endparam();
	if (ki == 0) {
	    if (global_clo_par1 < tmax-ptol) {
		bd_cvs.push_back(shared_ptr<CurveOnSurface>(bd_cv->subCurve(global_clo_par1, tmax)));
	    }
	} else if (ki == nmb_segments - 1) {
	    if (global_clo_par2 > tmin+ptol) {
		bd_cvs.push_back(shared_ptr<CurveOnSurface>(bd_cv->subCurve(tmin, global_clo_par2)));
	    }
	} else {
	    bd_cvs.push_back(shared_ptr<CurveOnSurface>(bd_cv->clone()));
	}
    }
}


//===========================================================================
void BoundedSurface::turnOrientation()
//===========================================================================
{
    // Turn orientation is ambigous, could mean "swap parameter
    // directions or "reverse parameter direction u (or v)". Adding a
    // message for now, but this should be fixed at some point. @jbt
    LOG_INFO("Note: 'Turn orientation' is ambigous - did you \n"
	    "mean 'swap parameter directions'? Continuing...");

    box_.unset();
    surface_->turnOrientation();
    for (size_t ki=0; ki<boundary_loops_.size(); ki++) {
	boundary_loops_[ki]->turnOrientation();
	// @afr: Not sufficient. We must also swap U and V coordinates
	// on any parametric curves that used to be in CurveOnSurface
	// objects in all loops. This, we can only do for spline curves,
	// and the code will throw if we find a non-spline parameter curve.
	// Yes, it's ugly. Yes, it's a hack.
	// The problem is due to the fact that changes to the underlying surface
	// are not detected by the CurveOnSurface objects associated with it.

	// @afr: Why can't we just call swapParameterDirection()?
	// And from that implementation, it looks like it doesn't matter that the
	// underlying surf is a spline or that the curves are splines.
	for (int kj = 0; kj < boundary_loops_[ki]->size(); ++kj) {
	    shared_ptr<ParamCurve> cv = boundary_loops_[ki]->operator[](kj);
	    CurveOnSurface* scv
		= dynamic_cast<CurveOnSurface*>(cv.get());
	    if (scv) {
		if(!(scv->underlyingSurface() == surface_)) {
		    LOG_INFO("The boundary curves lie on the wrong surface!");
		}
		ParamCurve* pcv = scv->parameterCurve().get();
		if (pcv) {
		    SplineCurve* spcv = dynamic_cast<SplineCurve*>(pcv);
		    if (spcv) {
			for (int kk = 0; kk < spcv->numCoefs(); ++kk) {
			    swap(*(spcv->coefs_begin() + 2*kk),
				 *(spcv->coefs_begin() + 2*kk+1));
			}
		    } else {
			LOG_INFO("The parameter curve is not a spline curve, so we can't switch "
			      "U and V coordinates.");
		    }
		}
	    }
	}
    }
}


//===========================================================================
void BoundedSurface::reverseParameterDirection(bool direction_is_u)
//===========================================================================
{
  box_.unset();

  RectDomain dom = surface_->containingDomain();
  double u1 = dom.umin();
  double u2 = dom.umax();
  double v1 = dom.vmin();
  double v2 = dom.vmax();

  // Reverse parameter direction for underlying surface
  surface_->reverseParameterDirection(direction_is_u);

  // Handle parameter curves
  for (size_t ki = 0; ki < boundary_loops_.size(); ++ki)
    {
      boundary_loops_[ki]->turnOrientation();
      for (int kj = 0; kj < (*boundary_loops_[ki]).size(); ++kj) {
	shared_ptr<CurveOnSurface> cv
	  (dynamic_pointer_cast<CurveOnSurface, ParamCurve>
	   ((*boundary_loops_[ki])[kj]));
	LOG_INFO("Expecting a CurveOnSurface.");
	Point dir(u1+u2, v1+v2);
	int pdir = (direction_is_u) ? 1 : 2;
	bool done = cv->translateSwapParameterCurve(dir, -1, pdir);
	// shared_ptr<SplineCurve> trim_cv = 
	//   dynamic_pointer_cast<SplineCurve, ParamCurve>(cv->parameterCurve());
	// if (trim_cv.get())
	// 	{
	// 	  if (trim_cv->rational())
	// 	    {
	// 	      vector<double>::iterator iter = trim_cv->rcoefs_begin();
	// 	      while (iter != trim_cv->rcoefs_end()) {
	// 		double w1 = iter[2];
	// 		if (direction_is_u)
	// 		  iter[0] = (u1 + u2 - iter[0]/w1)*w1;
	// 		else
	// 		  iter[1] = (v1 + v2 - iter[1]/w1)*w1;
	// 		iter+=3;
	// 	      }
	// 	      trim_cv->updateCoefsFromRcoefs();
	// 	    }
	// 	  else
	// 	    {
	// 	      vector<double>::iterator iter = trim_cv->coefs_begin();
	// 	      while (iter != trim_cv->coefs_end()) {
	// 		if (direction_is_u)
	// 		  iter[0] = u1 + u2 - iter[0];
	// 		else
	// 		  iter[1] = v1 + v2 - iter[1];
	// 		iter+=2;
	// 	      }
	// 	    }
	// 	}
	// else
	// 	{
	if (!done)
	  {
	    // Regenerate parameter curve
	    double eps = 1.0e-4;  // Arbitrary
	    cv->unsetParameterCurve();
	    cv->ensureParCrvExistence(eps);
	  }
      }
    }
}


//===========================================================================
void BoundedSurface::makeBoundaryCurvesG1(double kink)
//===========================================================================
{
  box_.unset();

    for (size_t ki = 0; ki < boundary_loops_.size(); ++ki) {
	vector<shared_ptr<ParamCurve> > curves;
	double space_epsilon = boundary_loops_[ki]->getSpaceEpsilon();
	for (int kj = 0; kj < boundary_loops_[ki]->size(); ++kj) {
	    shared_ptr<CurveOnSurface> cv
		(dynamic_pointer_cast<CurveOnSurface, ParamCurve>
		 ((*boundary_loops_[ki])[kj]));
	    vector<shared_ptr<CurveOnSurface> > cvs =
	        splitIntoC1Curves(cv, space_epsilon, kink);
	    curves.insert(curves.end(), cvs.begin(), cvs.end());
	}

	boundary_loops_[ki] =
	    shared_ptr<CurveLoop>(new CurveLoop(curves, space_epsilon));
    }
}


//===========================================================================
void BoundedSurface::removeSmallBoundaryCurves(double gap, double neighbour,
					       double kink)
//===========================================================================
{
  box_.unset();

    for (size_t ki = 0; ki < boundary_loops_.size(); ++ki) {
	vector<shared_ptr<ParamCurve> > curves;

	int loop_size = boundary_loops_[ki]->size();
	for (int kk = 0; kk < loop_size; ++kk) {
	  double len = (*boundary_loops_[ki])[kk]->estimatedCurveLength();
	  if (len < neighbour)
	    {
	      int kk1 = kk - 1;
	      if (kk1 < 0)
		kk1 = loop_size - 1;
	      int kk2 = (kk + 1) % loop_size;

	      // Check distance angles
	      vector<Point> pt1(2), pt2(2), pt3(2), pt4(2);
	      (*boundary_loops_[ki])[kk1]->point(pt1, 
						 (*boundary_loops_[ki])[kk1]->endparam(), 1);
	      (*boundary_loops_[ki])[kk]->point(pt2, 
						(*boundary_loops_[ki])[kk]->startparam(), 1);
	      (*boundary_loops_[ki])[kk]->point(pt3, 
						(*boundary_loops_[ki])[kk]->endparam(), 1);
	      (*boundary_loops_[ki])[kk2]->point(pt4, 
						 (*boundary_loops_[ki])[kk2]->startparam(), 1);
	      double ang1 = pt1[1].angle(pt2[1]);
	      double ang2 = pt3[1].angle(pt4[1]);
	      double d1 = pt1[0].dist(pt2[0]);
	      double d2 = pt3[0].dist(pt4[0]);
	      bool dismiss_next = false;
	      if ((d2 < gap && ang2 < kink) || (d1 < gap && ang1 < kink))
		{
		  if (d2 < gap && ang2 < kink)
		    {
		      // Merge with next curve
		      double dist;
		      (*boundary_loops_[ki])[kk]->appendCurve((*boundary_loops_[ki])[kk2].get(),
							      1, dist);
		      if (!(d1 < gap && ang1 < kink))
			curves.push_back((*boundary_loops_[ki])[kk]);
		      dismiss_next = true;  // Do not consider next curve
		      if (kk2 == 0)
			{
			  // Next curve already registered. Remove
			  curves.erase(curves.begin());
			}
		    }
		  if (d1 < gap && ang1 < kink)
		    {
		      double dist;
		      (*boundary_loops_[ki])[kk1]->appendCurve((*boundary_loops_[ki])[kk].get(),
							       1, dist);
		      if (kk1 < kk)
			curves.pop_back();  // Remove last curve
		      curves.push_back((*boundary_loops_[ki])[kk1]);
		      if (kk1 == loop_size - 1)
			loop_size--;  // Do not consider last curve
		    }
		  if (dismiss_next)
		    kk++;
		}
	      else 
		curves.push_back((*boundary_loops_[ki])[kk]);
	    }
	  else 
	    curves.push_back((*boundary_loops_[ki])[kk]);
		  
	    }

	if (boundary_loops_[ki]->size() != (int)curves.size())
	  boundary_loops_[ki] =
	    shared_ptr<CurveLoop>(new CurveLoop(curves, gap));
    }
}

//===========================================================================
void BoundedSurface::swapParameterDirection()
//===========================================================================
{
  box_.unset();
//     shared_ptr<SplineSurface> under_surf
// 	= dynamic_pointer_cast<SplineSurface, ParamSurface>(surface_);
//     LOG_INFO("Expecting underlying surface to be a spline surface.");
//     under_surf->swapParameterDirection();
  surface_->swapParameterDirection();

    // Reverse parameter directions for each segment, flip x and y,
    // and then reverse the order of the segments in all loops. @jbt
    for (size_t ki = 0; ki < boundary_loops_.size(); ++ki) {
	for (int kj = 0; kj < boundary_loops_[ki]->size(); ++kj) {
	    (*boundary_loops_[ki])[kj]->reverseParameterDirection(true);
	}
	reverse(boundary_loops_[ki]->begin(), boundary_loops_[ki]->end());
    }

}


//===========================================================================
void BoundedSurface::setParameterDomain(double u1, double u2, double v1, double v2)
//===========================================================================
{
  RectDomain dom = surface_->containingDomain();
  double u1_prev = dom.umin();
  double u2_prev = dom.umax();
  double v1_prev = dom.vmin();
  double v2_prev = dom.vmax();

  surface_->setParameterDomain(u1, u2, v1, v2);

  for (size_t ki = 0; ki < boundary_loops_.size(); ++ki)
    for (int kj = 0; kj < (*boundary_loops_[ki]).size(); ++kj) {
      shared_ptr<CurveOnSurface> cv
	(dynamic_pointer_cast<CurveOnSurface, ParamCurve>
	 ((*boundary_loops_[ki])[kj]));
      LOG_INFO("Expecting a CurveOnSurface.");
      cv->setDomainParCrv(u1, u2, v1, v2, u1_prev, u2_prev, v1_prev, v2_prev);
    }
}

//===========================================================================
void BoundedSurface::setParameterDomainBdLoops(double u1, double u2, 
					       double v1, double v2)
//===========================================================================
{
  RectDomain dom = surface_->containingDomain();
  double u1_prev = dom.umin();
  double u2_prev = dom.umax();
  double v1_prev = dom.vmin();
  double v2_prev = dom.vmax();

  for (size_t ki = 0; ki < boundary_loops_.size(); ++ki)
    for (int kj = 0; kj < (*boundary_loops_[ki]).size(); ++kj) {
      shared_ptr<CurveOnSurface> cv
	(dynamic_pointer_cast<CurveOnSurface, ParamCurve>
	 ((*boundary_loops_[ki])[kj]));
      LOG_INFO("Expecting a CurveOnSurface.");
      cv->setDomainParCrv(u1, u2, v1, v2, u1_prev, u2_prev, v1_prev, v2_prev);
    }
}

//===========================================================================
bool BoundedSurface:: hasUnderlyingSpline(shared_ptr<SplineSurface>& srf)
//===========================================================================
{
  shared_ptr<ParamSurface> sf = surface_;
  shared_ptr<BoundedSurface> bd_sf = 
    dynamic_pointer_cast<BoundedSurface, ParamSurface>(sf);
  while (bd_sf.get())
    {
      sf = bd_sf->underlyingSurface();
      bd_sf = dynamic_pointer_cast<BoundedSurface, ParamSurface>(sf);
    }
  srf = dynamic_pointer_cast<SplineSurface, ParamSurface>(sf);

  return (srf.get() != NULL);
}

//===========================================================================
bool BoundedSurface::isDegenerate(bool& b, bool& r,
				    bool& t, bool& l, double tolerance) const
//===========================================================================
{
  // Check the underlying surface
  bool degen = surface_->isDegenerate(b, r, t, l, tolerance);

  if (degen)
    {
      // Check if the degenerate boundaries lies inside the bounded surface
      RectDomain dom = surface_->containingDomain();

      if (b)
	{
	  vector<shared_ptr<ParamCurve> > cvs = 
	    constParamCurves(dom.vmin(), true);
	  if (cvs.size() == 0)
	    b = false;
	}

      if (r)
	{
	  vector<shared_ptr<ParamCurve> > cvs = 
	    constParamCurves(dom.umax(), false);
	  if (cvs.size() == 0)
	    r = false;
	}

      if (t)
	{
	  vector<shared_ptr<ParamCurve> > cvs = 
	    constParamCurves(dom.vmax(), true);
	  if (cvs.size() == 0)
	    t = false;
	}

      if (l)
	{
	  vector<shared_ptr<ParamCurve> > cvs = 
	    constParamCurves(dom.umin(), false);
	  if (cvs.size() == 0)
	    l = false;
	}

      degen = (b || r || t || l);
    }

  return degen;
}

//===========================================================================
void BoundedSurface::getDegenerateCorners(vector<Point>& deg_corners, double tol) const
//===========================================================================
{
    // Fetch all degenerate corners from the underlying surface
    vector<Point> curr_deg;
    surface_->getDegenerateCorners(curr_deg, tol);
    const CurveBoundedDomain& dom = parameterDomain();
    for (size_t ki = 0; ki<curr_deg.size(); ++ki)
    {
	// Check if the current corner lies within the domain
	Array<double,2> pnt2d(curr_deg[ki][0],curr_deg[ki][1]);
	bool in_domain = dom.isInDomain(pnt2d, tol);
	if (in_domain)
	    deg_corners.push_back(curr_deg[ki]);
    }
}

//===========================================================================
void BoundedSurface::getCornerPoints(vector<pair<Point,Point> >& corners) const
//===========================================================================
{
  for (size_t ki=0; ki<boundary_loops_.size(); ++ki)
    {
      int size = boundary_loops_[ki]->size();
      for (int kj=0; kj<size; ++kj)
	{
	  shared_ptr<ParamCurve> crv = (*boundary_loops_[ki])[kj];
	  Point par = getSurfaceParameter((int)ki, (int)kj, crv->startparam());
	  Point pos = ParamSurface::point(par[0], par[1]);
	  corners.push_back(make_pair(pos, par));
	}
    }
}

//===========================================================================
void BoundedSurface::splitSingleLoops()
//===========================================================================
{
    // Single loop may be connected to identical loop, hence 2 is not a good idea.
    int nmb_new_segments = 3;

    for (size_t ki = 0; ki < boundary_loops_.size(); ++ki)
	if (boundary_loops_[ki]->size() == 1) { // we split
	    // Currently the topology analysis does not detect edges joining smoothly.
	    // This is the case for two loops, not too uncommon a scenario.
	    // We run through edges, checking whether any of them appear to be a loop.
	    // If so, edge is split in the middle.
	    shared_ptr<ParamCurve> loop = (*boundary_loops_[ki])[0];
	    vector<shared_ptr<ParamCurve> > new_curves;
	    double tmin = loop->startparam();
	    double tmax = loop->endparam();
	    double tstep = (tmax - tmin)/(nmb_new_segments);
	    for (int kj = 0; kj < nmb_new_segments; ++kj) {
		double split_t_low = tmin + kj*tstep;
		double split_t_high = tmin + (kj + 1)*tstep;
		new_curves.push_back(shared_ptr<ParamCurve>(loop->subCurve(split_t_low, split_t_high)));
	    }
	    (*boundary_loops_[ki]) = CurveLoop(new_curves, boundary_loops_[ki]->getSpaceEpsilon());
	}
}


//===========================================================================
SplineCurve* BoundedSurface::constParamCurve(double parameter,
					       bool direction_is_u) const
//===========================================================================
{
    // // We extract constParamCurve from the underlying surface (assuming it is a spline),
    // // then extract the part lying on the surface. If return curve is not continuous,
    // // function fails (GO_ERROR).
    // shared_ptr<SplineSurface> under_sf =
    // 	dynamic_pointer_cast<SplineSurface, ParamSurface>(surface_);
    // LOG_INFO("Expecting underlying surface to be a spline surface.");

    int pardir = direction_is_u ? 1 : 2;
    double tolerance = 1e-05;
    vector<shared_ptr<CurveOnSurface> > trim_pieces;
    const CurveBoundedDomain& dom = parameterDomain();

    dom.clipWithDomain(pardir, parameter, tolerance, surface_, trim_pieces);

    LOG_INFO("Expecting iso curve to be connected...");

    return dynamic_cast<SplineCurve*>
      (trim_pieces[0]->spaceCurve()->geometryCurve()); // Return cv is NEWed.
}



//===========================================================================
bool
BoundedSurface::isIsoTrimmed(double tol) const
//===========================================================================
{
    if (iso_trim_tol_ >= 0.0 && iso_trim_tol_ - tol < 1.0e-12)
    {
      return (iso_trim_ > 0);   // Already checked with "the same" or smaller tolerance
    }
    iso_trim_tol_ = tol;
    iso_trim_ = 1;  // Until the opposite is found

    if (boundary_loops_.size() == 0)
    {
	// Surface not trimmed?
      return (iso_trim_ > 0);
    }

    if (boundary_loops_.size() > 1)
    {
	// The surface has inner trimming curves. Not iso-trimmed
	iso_trim_ = 0;
	return (iso_trim_ > 0);
    }
    
    // Check boxes around the 2D curves. NB! Assumes that the 2D curves is OK
   int nmb_crvs = boundary_loops_[0]->size();
    vector<double> par_u, par_v;
    for (int ki=0; ki<nmb_crvs; ++ki)
    {
	shared_ptr<ParamCurve> crv = (*boundary_loops_[0])[ki];
	shared_ptr<ParamCurve> pcrv;
	if (crv->dimension() == 2)
	    pcrv = crv;
	else
	{
	    shared_ptr<CurveOnSurface> sf_cv =
	      dynamic_pointer_cast<CurveOnSurface,ParamCurve>(crv);
	    pcrv = sf_cv->parameterCurve();
	}
	if (!pcrv.get())
	{
	    // No 2D parameter curve
	    iso_trim_ = 0;
	    return (iso_trim_ > 0);
	}

	// Make box around curve
	BoundingBox box2d = pcrv->boundingBox();
	Point high = box2d.high();
	Point low = box2d.low();
	if (high[0] - low[0] > tol && high[1]-low[1] > tol)
	{
	    // 2D curve different from a iso-line
	  iso_trim_ = 0;
	  return (iso_trim_ > 0);
	}
	if (high[0] - low[0] <= tol)
	  par_u.push_back(0.5*(low[0]+high[0]));
	if (high[1] - low[1] <= tol)
	  par_v.push_back(0.5*(low[1]+high[1]));
    }

    // Count number of distinct constant parameter values in both directions
    std::sort(par_u.begin(), par_u.end());
    std::sort(par_v.begin(), par_v.end());
    int nmb = 0;
    size_t kj;
    for (kj=1; kj<par_u.size(); kj++)
      if (par_u[kj] - par_u[kj-1] > tol)
	nmb++;
    if (nmb > 1)
      iso_trim_ = 0;

    nmb = 0;
    for (kj=1; kj<par_v.size(); kj++)
      if (par_v[kj] - par_v[kj-1] > tol)
	nmb++;
    if (nmb > 1)
      iso_trim_ = 0;

    // Check if the iso parameter coincides with the boundary
    if (surface_->instanceType() != Class_BoundedSurface)
      {
	RectDomain dom = surface_->containingDomain();
	int nmb1=0, nmb2=0;
	for (kj=0; kj<par_u.size(); kj++)
	  if (fabs(par_u[kj]-dom.umin()) < tol || fabs(par_u[kj]-dom.umax()) < tol)
	    nmb1++;
	for (kj=0; kj<par_v.size(); kj++)
	  if (fabs(par_v[kj]-dom.vmin()) < tol || fabs(par_v[kj]-dom.vmax()) < tol)
	    nmb2++;
	if (nmb1 == (int)par_u.size() && nmb2 == (int)par_v.size())
	  iso_trim_ = 2;
      }
	
    return (iso_trim_ > 0);
}

//===========================================================================
bool
BoundedSurface::isBoundaryTrimmed(double tol) const
//===========================================================================
{
    if (iso_trim_tol_ >= 0.0 && iso_trim_tol_ - tol < 1.0e-12)
    {
      return (iso_trim_ > 1);   // Already checked with "the same" or smaller tolerance
    }
    
    bool iso_trim = isIsoTrimmed(tol);
    return (iso_trim_ > 1);
}

//===========================================================================
shared_ptr<ParamSurface> 
BoundedSurface::getIsoTrimSurface(double tol) const
//===========================================================================
{
  shared_ptr<ParamSurface> dummy;
  if (!isIsoTrimmed(tol))
    return dummy;

  // Fetch parameter domain
  RectDomain dom1 = containingDomain();
  RectDomain dom2 = surface_->containingDomain();
  double umin = std::max(dom1.umin(), dom2.umin());
  double umax = std::min(dom1.umax(), dom2.umax());
  double vmin = std::max(dom1.vmin(), dom2.vmin());
  double vmax = std::min(dom1.vmax(), dom2.vmax());

  vector<shared_ptr<ParamSurface> > sub_sfs = 
    surface_->subSurfaces(umin, vmin, umax, vmax);
  if (sub_sfs.size() == 1)
    return sub_sfs[0];
  else
    return dummy;
}

//===========================================================================
void
BoundedSurface::turnLoopOrientation(int idx)
//===========================================================================
{
  box_.unset();

    if (loop_fixed_.size() != boundary_loops_.size())
    {
	loop_fixed_.resize(boundary_loops_.size());
	std::fill(loop_fixed_.begin(), loop_fixed_.end(), 0);
    }
	    
    if (idx >= 0 && idx < (int)boundary_loops_.size())
    {
	boundary_loops_[idx]->turnOrientation();
	loop_fixed_[idx] = 1 - loop_fixed_[idx];
    }
}


//===========================================================================
bool BoundedSurface::isValid(int& valid_state) const
//===========================================================================
{
    valid_state = valid_state_;
    return (valid_state_ == 1);
}


//===========================================================================
void BoundedSurface::analyzeLoops()
//===========================================================================
{
    // We then analyze the boundary curves, starting with state = -1
    // etc.
    bool analyze = true;

    // Then we see if the par cv and the space cv match.
    double max_tol_mult = 1.0;
    int nmb_seg_samples = 100;//20;
    bool cv_match_ok = fixParSpaceMismatch(analyze, max_tol_mult,
					   nmb_seg_samples);

    // We then look for missing par cvs.
    bool par_cv_missing = parameterCurveMissing();

    // We first check if the loop is closed within loop tolerance.
    double max_loop_gap = -1.0; // Only set if fixLoopGaps returns false.
    bool loop_gaps_ok;
    try {
	loop_gaps_ok= fixLoopGaps(max_loop_gap, analyze);
    } catch (...) {
	loop_gaps_ok = false;
    }

    // There is no point in testing the order of the boundary loops if
    // any of the above failed! That routine needs valid parameter
    // loops.

    // Finally we see if the loop order and orientation is according
    // to requirements (one ccw outer loop as the first element, the
    // rest should be cw and lying inside ccw loop). Expecting sf to
    // be connected, i.e. only on ccw loop.
    double degenerate_epsilon = DEFAULT_SPACE_EPSILON;

    valid_state_ = 0;
    if (!cv_match_ok)
	valid_state_ += -1;
    if (par_cv_missing)
	valid_state_ += -2;
    if (!loop_gaps_ok)
	valid_state_ += -4;

    bool loop_order_ok = false;
    if ((!par_cv_missing) && (cv_match_ok) && loop_gaps_ok) {
	loop_order_ok = orderBoundaryLoops(analyze, degenerate_epsilon);
	if (!loop_order_ok)
	    valid_state_ += -8;
    }

    if ((!par_cv_missing) && cv_match_ok && loop_gaps_ok && loop_order_ok)
	valid_state_ = 1;

#ifdef SBR_DBG
    if (0)//valid_state_ != 1)
    {
	LOG_INFO("valid_state_: " + std::to_string(valid_state_) + ", par_cv_missing: " + std::to_string(par_cv_missing) + ", cv_match_ok: " + std::to_string(cv_match_ok) + ", loop_gaps_ok: " + std::to_string(loop_gaps_ok) + ", loop_order_ok: " + std::to_string(loop_order_ok));
    }
#endif
}


//===========================================================================
void BoundedSurface::removeMismatchCurves(double max_tol_mult)
//===========================================================================
{
    if (valid_state_ > 0)
	return;

    box_.unset();

    bool analyze = false;
    int nmb_seg_samples = 20;//100;
    bool cvs_match;
    cvs_match = fixParSpaceMismatch(analyze, max_tol_mult, nmb_seg_samples);
    analyzeLoops();
}

//===========================================================================
void BoundedSurface::fixMismatchCurves(double eps)
//===========================================================================
{
  for (size_t ki=0; ki<boundary_loops_.size(); ++ki)
    boundary_loops_[ki]->fixMismatchCurves(eps);
} 

//===========================================================================
bool BoundedSurface::fixInvalidSurface(double& max_loop_gap, double max_tol_mult)
//===========================================================================
{
    // If the the valid_state_ flag was not set we must analyze the loops.
    if (valid_state_ == 0)
    {
        analyzeLoops();
    }

    if (max_tol_mult < 1.0)
    {
	max_tol_mult = 1.0;
    }

    max_loop_gap = -1.0; // Just in case the user did not initialize the value.
    if (valid_state_ == 1) {
	return true; // Nothing to be done.
    }

    box_.unset();

#ifdef SBR_DBG
    LOG_DEBUG("Must fix invalid surface! valid_state_ = " + std::to_string(valid_state_));
#endif

    bool analyze = false;

    // We first try to fix mismatch between par and space cvs.
    if ((int)fabs(double(valid_state_))%2 > 0) {
	int nmb_seg_samples = 20;//100;
	bool cvs_match;
	cvs_match = fixParSpaceMismatch(analyze, max_tol_mult,
					     nmb_seg_samples);
	analyzeLoops();
    }

    // valid_state_ == -2 is not handled. Projection of geometry
    // curves should be handled on the outside of this class.

    // We then try to fix gaps in the loops.
    if ((int)fabs(double(valid_state_))%8 > 3) {
	bool loop_gaps_ok;
	try {
	    loop_gaps_ok = fixLoopGaps(max_loop_gap, analyze);
	} catch (...) {
	    loop_gaps_ok = false;
	}
	analyzeLoops();
    }

    // Finally we try to sort (and possibly reverse orientation) of
    // the loops such that the outer loop lies first and is ccw, and
    // the rest of the loops are cw and inside the outer loop.
    if ((int)fabs(double(valid_state_))%16 > 7) {
	double degenerate_epsilon = DEFAULT_SPACE_EPSILON;
	bool loop_order_ok;
	loop_order_ok = orderBoundaryLoops(analyze, degenerate_epsilon);
	analyzeLoops();
	// If reordering of loops was a success, the valid_state_ should
	// have been updated.
    }
    if (valid_state_ == 1)
	return true;
    else
	return false;
}


//===========================================================================
bool BoundedSurface::fixLoopGaps(double& max_loop_gap, bool analyze)
//===========================================================================
{
    if (valid_state_ == 1)
	return true; // No point in calling this routine.

    // Do not call this routine if either par and space cv do not match or a parameter curve is missing.
    if ((analyze == false) && (((-valid_state_)%4) > 0)) // Note that valid_state is either 0 or negative.
	return true;

    // If the gaps are ok there is nothing to be done.
    if ((analyze == false) && (((-valid_state_)/4) > 1)) // Note that valid_state is either 0 or negative.
	return true;

    box_.unset();

    max_loop_gap = -1.0;
    // We check if the loops are valid.
#ifdef SBR_DBG
    if (!analyze)
		LOG_DEBUG("Invalid loop!");
#endif
    bool all_loops_valid = true;
    // Loops are not closed within tolerance.
    for (size_t ki = 0; ki < boundary_loops_.size(); ++ki) {
	bool valid_loop = boundary_loops_[ki]->isValid();
	if (!valid_loop) {
	    double max_gap = -1.0;
	    bool success = boundary_loops_[ki]->fixInvalidLoop(max_gap);
	    //bool success = false;
            if (success)
            {
                if (max_gap > max_loop_gap)
                    max_loop_gap = max_gap;
            }
            else
            {
		all_loops_valid = false;
		LOG_INFO("Failed fixing invalid loop! max_gap = " + std::to_string(max_gap) + ", epsgeo = " + std::to_string(boundary_loops_[ki]->getSpaceEpsilon()));
		break;
	    }
	}
    }

    if (all_loops_valid) {
	return true;
    } else { // Ordering the loops requires closed and simple loops.
	LOG_INFO("Failed fixing loop gap(s). BoundedSurface is still invalid.");
	return false;
    }
}

 
//===========================================================================
double BoundedSurface::maxLoopSfDist(int loop_ind, int nmb_seg_samples) const
//===========================================================================
{
    ASSERT(loop_ind >= 0 && loop_ind < (int)boundary_loops_.size());
    // Assuming loop is made of CurveOnSurface objects.
    shared_ptr<CurveLoop> loop = boundary_loops_[loop_ind];
    int nmb_segments = loop->size();
    double epsgeo = loop->getSpaceEpsilon();
    double max_sf_dist = -1;
    bool verified_not_closed = false; // If needed we should check if the underlying surface is closed in either dirs.
    for (int ki = 0; ki < nmb_segments; ++ki) {
	if ((*loop)[ki]->instanceType() == Class_CurveOnSurface) {
	    shared_ptr<CurveOnSurface> cv_on_sf =
		dynamic_pointer_cast<CurveOnSurface, ParamCurve>((*loop)[ki]);
	    shared_ptr<ParamSurface> sf = cv_on_sf->underlyingSurface();
	    shared_ptr<ParamCurve> pcv = cv_on_sf->parameterCurve();
	    double tmin = cv_on_sf->startparam();
	    double tmax = cv_on_sf->endparam();
	    double tstep = (tmax-tmin)/(double)(nmb_seg_samples-1);
	    double max_dist = -1.0;
	    double tpar;
	    double seed[2];
	    for (int kj = 0; kj < nmb_seg_samples; ++kj)
	    {
		tpar = (kj == nmb_seg_samples - 1) ? tmax : tmin + kj*tstep;
		Point cv_pt = cv_on_sf->ParamCurve::point(tpar);
		double clo_u, clo_v, clo_dist;
		Point clo_pt;
		double* local_seed = NULL;
		Point par_pt;
		if (pcv)
		{
		    par_pt = pcv->point(tpar);
		    local_seed = &par_pt[0];
		}
		else if (verified_not_closed && (kj > 0))
		{   // If verified_closed we should check if the seed should switch side.
                    const RectDomain& rect_dom = containingDomain();
                    if (seed[0] < rect_dom.umin() || seed[0] > rect_dom.umax() ||
                        seed[1] < rect_dom.vmin() || seed[1] > rect_dom.vmax())
                    {
                        LOG_INFO("Seed is outside the domain!");
                    }
		    local_seed = seed;
		}
		// if (kj == 0)
		sf->closestPoint(cv_pt, clo_u, clo_v,
				 clo_pt, clo_dist, epsgeo, NULL, local_seed);
		// else
		//     sf->closestPoint(cv_pt, clo_u, clo_v,
		// 		     clo_pt, clo_dist, epsgeo, NULL, seed);
		// We also check towards the boundary, may be more stable for areas with high curvature.
		if ((clo_dist > max_sf_dist) && (closeToUnderlyingBoundary(clo_u, clo_v)))
		{
		    double bd_clo_u, bd_clo_v, bd_clo_dist;
		    Point bd_clo_pt;
		    seed[0] = clo_u;
		    seed[1] = clo_v;
		    sf->closestBoundaryPoint(cv_pt, bd_clo_u, bd_clo_v, 
					     bd_clo_pt, bd_clo_dist, epsgeo, NULL, seed);
		    if (bd_clo_dist < clo_dist)
		    {
			clo_dist = bd_clo_dist;
			clo_u = bd_clo_u;
			clo_v = bd_clo_v;
		    }
		}
		if (clo_dist > max_dist)
		    max_dist = clo_dist;
		seed[0] = clo_u;
		seed[1] = clo_v;
	    }

	    if (max_dist > max_sf_dist)
		max_sf_dist = max_dist;
	}
    }
    return max_sf_dist;
}


//===========================================================================
double BoundedSurface::maxLoopGap()
//===========================================================================
{
    double max_loop_gap = -1.0;
    const bool analyze = true;
    for (size_t ki = 0; ki < boundary_loops_.size(); ++ki)
    {
	double loop_gap = boundary_loops_[ki]->getMaxCurveDist();
	if (loop_gap > max_loop_gap)
	{
	    max_loop_gap = loop_gap;
	}
    }

    return max_loop_gap;
}

//===========================================================================
bool
BoundedSurface::orderBoundaryLoops(bool analyze, double degenerate_epsilon)
//===========================================================================
{
    if (valid_state_ > 0)
	return true;
    else if ((!analyze) && ((-valid_state_)%16 <= 7))
	LOG_INFO("Unexpected valid_state_ = " + std::to_string(valid_state_));

    // We're assuming that we are given only one outer loop. Hence
    // only one loop shall be given an ccw orientation, all other
    // loops should be cw.

    CurveBoundedDomain bounded_domain;
    shared_ptr<ParamCurve> pcurve;
//     vector<int> outer_loop(boundary_loops_.size(), 0);

    // We store the orientation of the loops.
    vector<int> loop_is_ccw(boundary_loops_.size());

    double tpar;
    Point par_point; // Point in parameter domain.
    // Method: We construct a bounded surface for each loop, then
    // check if a random point on each of the other curves is on
    // bounded surface.

    // @@sbr072009 What about method firstLoopInsideSecond()?

    // We store all loops the a specific loop lies inside. With
    // "inside" we mean part of the bounded part, regardless of loop
    // orientation.

    //const double int_tol = 1e-03;
//    const double int_tol = GoTools::spaceEpsilon();
    const double int_tol = GoTools::parameterEpsilon();

    vector<vector<int> > lies_inside_loop(boundary_loops_.size());
    for (size_t ki = 0; ki < boundary_loops_.size(); ++ki) {
      //const double int_tol = 1e-06;
	bool ccw;
	try {
	    ccw = LoopUtils::loopIsCCW(*boundary_loops_[ki], int_tol);
	} catch (...) {
	    LOG_INFO("Failed analyzing loop.");
	    return false;
	}
	loop_is_ccw[ki] = (ccw) ? 1 : 0;

	// @@sbr072009 But CurveBoundedDomain expects ccw loop ...
	// Which would indicate that cw loops will treat loops outside
	// as inside ...
	bounded_domain = CurveBoundedDomain(boundary_loops_[ki]);
	// We sample a point on each of the other loops in the
	// parameter domain. We then test if it is part of the bounded
	// part of the parameter plane.
	size_t kj;
	for (kj = 0; kj < boundary_loops_.size(); ++kj) {
	    if (kj == ki)
		continue;
	    pcurve = (dynamic_pointer_cast<CurveOnSurface, ParamCurve>(
			 boundary_loops_[kj]->operator[](0)))->parameterCurve();
	    if (pcurve.get() == 0) {
		LOG_INFO("Missing parameter curve!");
		return false; // Method demands parameter curves.
	    }

	    tpar = 0.5*(pcurve->startparam() + pcurve->endparam());
	    pcurve->point(par_point, tpar); // We find mid point on first curve.

	    Vector2D p_point(par_point[0], par_point[1]);
	    // We may need to handle an exception.
	    bool is_in_domain;
	    try {
		// Note that for a cw loop, the "in domain"-test
		// refers to whether the point is part of the bounded
		// region.
		is_in_domain =
		    bounded_domain.isInDomain(p_point, degenerate_epsilon);
		if (is_in_domain)
		    lies_inside_loop[ki].push_back((int)kj);
	    } catch (...) {
		LOG_INFO("Caught exception! Nothing more we can do ...");
		return false;
	    }
	}

// 	MESSAGE("Orientation is ccw? " << ccw << ", # loops inside: "
// 		<< lies_inside_loop[ki].size());
    }

    // We should be done computing, time to move the outer loop up
    // front.  The outer loop should be the one with all the other
    // loops inside (i.e. part of the bounded region as defined by
    // loop).
    int nmb_outer_loops = 0;
    int outer_index = -1;
    for (size_t ki = 0; ki < lies_inside_loop.size(); ++ki)
	if (lies_inside_loop[ki].size() == lies_inside_loop.size() -1) {
	    ++nmb_outer_loops;
	    outer_index = (int)ki;
	}

#ifdef SBR_DBG
    std::ofstream debug("tmp/debug.g2");
    surface_->writeStandardHeader(debug);
    surface_->write(debug);
    for (size_t ki = 0; ki < (int)boundary_loops_.size(); ++ki)
	for (size_t kj = 0; kj < boundary_loops_[ki]->size(); ++kj) {
	    shared_ptr<ParamCurve> cv = (*boundary_loops_[ki])[kj];
	    if (cv->instanceType() == Class_CurveOnSurface) {
		shared_ptr<CurveOnSurface> cv_on_sf =
		    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(cv);
		if (cv_on_sf->parameterCurve() != NULL) {
		    shared_ptr<SplineCurve> pcv =
			dynamic_pointer_cast<SplineCurve, ParamCurve>
			(cv_on_sf->parameterCurve());
		    if (pcv.get() != NULL)
			SplineDebugUtils::writeSpaceParamCurve(*pcv, debug, 0.0);
		    else {
			cv_on_sf->parameterCurve()->writeStandardHeader(debug);
			cv_on_sf->parameterCurve()->write(debug);
		    }
		}
		if (cv_on_sf->spaceCurve() != NULL) {
		    cv_on_sf->spaceCurve()->writeStandardHeader(debug);
		    cv_on_sf->spaceCurve()->write(debug);
		}
	    } else {
		cv->writeStandardHeader(debug);
		cv->write(debug);
	    }
	}
    double debug_val = 0.0;
#endif

    if (nmb_outer_loops != 1)
    {
	if (analyze)
	{
// 	    valid_state_ += -2;
	    LOG_INFO(std::to_string(boundary_loops_.size()) + " loops in total. Found " + std::to_string(nmb_outer_loops) + " outer loops! BoundedSurface invalid. surface_->classtype(): " + std::to_string(surface_->instanceType()));
	    return false;
	}
	else {
	    LOG_INFO(std::to_string(boundary_loops_.size()) + " loops in total. Failed finding one outer loop (" + std::to_string(nmb_outer_loops) + ")! BoundedSurface invalid.");
// 	    MESSAGE("No success in finding a single outer loop!");
	    return false;
	}
    }

    // If not already up front, we move it there.
    shared_ptr<CurveLoop> dummy_loop;
    if (outer_index != 0) {
	if (analyze) {
	    LOG_INFO("First loop is not the outer loop, BoundedSurface invalid.");
	    return false;
	} else {
	    std::swap(boundary_loops_[0], boundary_loops_[outer_index]);
	    std::swap(loop_is_ccw[0], loop_is_ccw[outer_index]);
	}
    }
//     else if (analyze)
// 	return true;

    // Finally we make sure that the orientation of the loops are as
    // required, i.e. the first (and outer) loops should be ccw, the
    // "neighbor" loops inside should be cw, their "neighbor" loops
    // inside should be ccw etc.  Not that we normally encounter do
    // not more than a ccw loop, with possibly a few cw loops inside.
    if (!loop_is_ccw[0]) {
	if (analyze) {
	    return false;
	}
	else {
	    boundary_loops_[0]->turnOrientation(); // Reverse direction of loop.
	}
    }
    for (size_t ki = 1; ki < boundary_loops_.size(); ++ki) {
	if (loop_is_ccw[ki]) {
	    if (analyze) {
		return false;
	    }
	    else {
		boundary_loops_[ki]->turnOrientation();
	    }
	}
    }

//     MESSAGE("Method under construction");

    return true;
}


//===========================================================================
bool BoundedSurface::parameterCurveMissing()
//===========================================================================
{
    // We run through all loops, checking whether all parameter curves
    // are present.
    for (size_t ki = 0; ki < boundary_loops_.size(); ++ki)
	for (int kj = 0; kj < boundary_loops_[ki]->size(); ++kj) {
	    shared_ptr<CurveOnSurface> cv_on_sf =
		dynamic_pointer_cast<CurveOnSurface, ParamCurve>
		((*boundary_loops_[ki])[kj]);
	    ASSERT(cv_on_sf.get() != NULL);
	    if (cv_on_sf->parameterCurve().get() == NULL)
		return true;
	}

    return false;
}


//===========================================================================
bool BoundedSurface::fixParSpaceMismatch(bool analyze, double max_tol_mult,
					 int nmb_seg_samples)
//===========================================================================
{
    // We run through all loop segments, checking whether the
    // direction and trace of the parameter curve matches that of the
    // space curve, as well as the corresponding parameter domains.
//     int nmb_samples = 100;
    if (max_tol_mult < 1.0)
    {
	max_tol_mult = 1.0;
    }

    for (size_t ki = 0; ki < boundary_loops_.size(); ++ki) {
      double loop_tol = boundary_loops_[ki]->getSpaceEpsilon();
      double space_eps = loop_tol;//1e03*loop_tol;
	bool cv_replaced = false;
	vector<shared_ptr<ParamCurve> > new_loop_cvs
	    (boundary_loops_[ki]->size());
	for (int kj = 0; kj < boundary_loops_[ki]->size(); ++kj) {
	    new_loop_cvs[kj] = (*boundary_loops_[ki])[kj];
	    shared_ptr<CurveOnSurface> cv_on_sf =
		dynamic_pointer_cast<CurveOnSurface, ParamCurve>
		((*boundary_loops_[ki])[kj]);
	    ASSERT(cv_on_sf.get() != NULL);
	    shared_ptr<ParamCurve> par_cv = cv_on_sf->parameterCurve();
	    shared_ptr<ParamCurve> space_cv = cv_on_sf->spaceCurve();
	    // Routine only checks for mismatch between existing curves.
	    if ((par_cv.get() == NULL) || (space_cv.get() == NULL))
		continue;

	    bool same_par_domain = cv_on_sf->sameParameterDomain(); // tol = 1e-12.
	    if ((!same_par_domain) && analyze)
		return false;

	    // We then check if orientation is correct, according to
	    // end pts.
	    bool same_orientation = cv_on_sf->sameOrientation();
	    if ((!same_orientation) && analyze)
		return false;

	    // There is no reason to use the loop tolerance as max
	    // allowed dist between par and space cv.
	    bool cv_consistent = true;
	    if ((!same_orientation) || (!same_par_domain)) {
// 		MESSAGE("Fixing cv orientation or domain!");
		if (analyze)
		    return false;
		bool pref_par = cv_on_sf->parPref();
		int ccm = cv_on_sf->curveCreationMethod();
		// How the curve was created is more to trust.
		bool prefer_parameter = (ccm == 1) ? false : pref_par;
// 		cv_consistent =
		cv_on_sf->makeCurvesConsistent(prefer_parameter);
	    }

	    // Finally we check if curves have matching image/trace.	    
	    bool same_trace = cv_on_sf->sameTrace(space_eps, nmb_seg_samples);
	    if (!same_trace) {
		double max_trace_diff = cv_on_sf->maxTraceDiff(nmb_seg_samples);
		if (max_trace_diff < max_tol_mult*space_eps) {
		    if (!analyze) {
			// MESSAGE("max_trace_diff = " << max_trace_diff <<
			// 	".Altering tolerance! From: " << space_eps <<
			// 	", to: " << 1.1*max_trace_diff);
			boundary_loops_[ki]->setSpaceEpsilon
			    (max_tol_mult*max_trace_diff);
			space_eps = boundary_loops_[ki]->getSpaceEpsilon();
		    }
		    same_trace = cv_on_sf->sameTrace(space_eps,
						     nmb_seg_samples);
		    if (!same_trace)
			LOG_INFO("Strange, this should not happen ...");
		} else {
		    if (max_trace_diff > 1.0) // Rather arbitrary. But we want to detect bugs, not bad input tolerance.
		    { 
			LOG_INFO("Deviation too large, epsgeo: " + std::to_string(space_eps) + ", max_trace_diff: " + std::to_string(max_trace_diff));
		    }
		}
	    }
	    if ((!cv_consistent) || (!same_trace)) {
		if (analyze) {
		    return false;
		}
		else {
		    // We have to remove one of the
		    // curves. @@sbr072009 Remains to check whether
		    // the remaining curve matches parameter domain /
		    // surface location.
		    cv_replaced = true;
 		    bool par_pref = cv_on_sf->parPref();
		    int ccm = cv_on_sf->curveCreationMethod();
		    bool remove_space = (ccm == 1) ? false : par_pref;
		    if (remove_space) {
			new_loop_cvs[kj] =
			    shared_ptr<ParamCurve>
			    (new CurveOnSurface(cv_on_sf->underlyingSurface(),
						cv_on_sf->parameterCurve(),
						true));
		    }
		    else {
			new_loop_cvs[kj] =
			    shared_ptr<ParamCurve>
			    (new CurveOnSurface(cv_on_sf->underlyingSurface(),
						cv_on_sf->spaceCurve(),
						false));
		    }
		}
	    }
	}
	if (cv_replaced)
	    boundary_loops_[ki]->setCurves(new_loop_cvs);
    }

    return true;
}


//===========================================================================
vector<shared_ptr<CurveOnSurface> >
BoundedSurface::splitIntoC1Curves(shared_ptr<CurveOnSurface>& curve,
				    double space_epsilon, double kink)
    //===========================================================================
{
    // First make curve k-regular
    shared_ptr<SplineCurve> pcv = dynamic_pointer_cast
	<SplineCurve, ParamCurve>(curve->parameterCurve());
    shared_ptr<SplineCurve> gcv = dynamic_pointer_cast
	<SplineCurve, ParamCurve>(curve->spaceCurve());
    if (pcv.get() != 0)
	{
	    pcv->makeKnotStartRegular();
	    pcv->makeKnotEndRegular();
	}
    if (gcv.get() != 0)
	{
	    gcv->makeKnotStartRegular();
	    gcv->makeKnotEndRegular();
	}
    shared_ptr<ParamCurve> prefered_cv;
    if (curve->parPref())
	prefered_cv = curve->parameterCurve();
    else
        prefered_cv = curve->spaceCurve();
    shared_ptr<ElementaryCurve> elem_cv = dynamic_pointer_cast<ElementaryCurve, ParamCurve>(prefered_cv);
    vector<shared_ptr<CurveOnSurface> > return_curves;
    if (elem_cv.get() != 0)
      {
	// Elementary curves have continuity C-infinity everywhere
	return_curves.push_back(curve);
	return return_curves;
      }

    shared_ptr<SplineCurve> cv = dynamic_pointer_cast<SplineCurve, ParamCurve>(prefered_cv);
    ALWAYS_ERROR_IF(cv.get() == 0,
		    "Curve did not exist!");
    vector<double> joints;
    vector<double> cont(2);
    cont[0] = space_epsilon;
    cont[1] = kink;
    GeometryTools::getGnJoints(*cv, cont, joints);

    shared_ptr<CurveOnSurface> temp_cv;
    int i;
    for (i = 0; i < int(joints.size()) - 1; ++i) {
	if (joints.size() == 2) {
	    return_curves.clear();
	    return_curves.push_back(curve);
	} else {
	    try {
		temp_cv = shared_ptr<CurveOnSurface>
		    (dynamic_cast<CurveOnSurface*>
		     (curve->subCurve(joints[i], joints[i+1])));
	    } catch (...) {
		LOG_INFO("Failed extracting subcurve on [" + std::to_string(joints[i]) + ", " + std::to_string(joints[i+1]) + "].");
		joints.erase(joints.begin() + i + 1);
		--i;
		continue; // @@sbr Handle this in a more stable manner?
	    }
	    return_curves.push_back(temp_cv);
	}
    }	
    return return_curves;
}


//===========================================================================
double BoundedSurface::getEpsGeo() const
//===========================================================================
{
    double min_space_eps = boundary_loops_[0]->getSpaceEpsilon();
    for (size_t ki = 1; ki < boundary_loops_.size(); ++ki) {
	double space_eps = boundary_loops_[ki]->getSpaceEpsilon();
	if (space_eps < min_space_eps)
	    min_space_eps = space_eps;
    }

    return min_space_eps;
}


//===========================================================================
bool BoundedSurface::closeToUnderlyingBoundary(double upar, double vpar,
					       double domain_fraction) const
//===========================================================================
{
    const RectDomain& rect_domain = surface_->containingDomain();
    double umin = rect_domain.umin();
    double umax = rect_domain.umax();
    double vmin = rect_domain.vmin();
    double vmax = rect_domain.vmax();
    double eps_u = (umax - umin)*domain_fraction;
    double eps_v = (vmax - vmin)*domain_fraction;
    if ((fabs(upar - umin) < eps_u) || (fabs(umax - upar) < eps_u)
	|| (fabs(vpar - vmin) < eps_v) || (fabs(vmax - vpar) < eps_v))
	return true;
    else
	return false;
}


//===========================================================================
int BoundedSurface::ElementOnBoundary(int elem_ix, double eps)
//===========================================================================
{
  if (surface_->instanceType() != Class_SplineSurface)
    return -1;
  
  // Fetch point internal to element
  // First fetch associated spline volume
  SplineSurface *sf = surface_->asSplineSurface();
  if (!sf)
    return -1;
  
  // Fetch curves in the parameter domain surrounding the specified element
  double elem_par[4];
  vector<shared_ptr<SplineCurve> > side_cvs = sf->getElementBdParCvs(elem_ix,
								     elem_par);
  if (side_cvs.size() == 0)
    return -1;

  // Check for intersections with the boundary curves
  for (int ki=0; ki<(int)boundary_loops_.size(); ++ki)
    {
      int nmb_crvs = boundary_loops_[ki]->size();
      shared_ptr<ParamCurve> pcrv;
      for (int kj=0; kj<nmb_crvs; ++kj)
	{
	  shared_ptr<ParamCurve> crv = (*boundary_loops_[ki])[kj];
	  shared_ptr<CurveOnSurface> sf_cv;
	  if (crv->instanceType() == Class_CurveOnSurface) 
	    sf_cv = dynamic_pointer_cast<CurveOnSurface, ParamCurve>(crv);
	  else
	    sf_cv = 
	      shared_ptr<CurveOnSurface>(new CurveOnSurface(surface_, crv, false));

	  // Check if the trimming curve is also a boundary curve
	  bool orient;
	  int bd = sf_cv->whichBoundary(eps, orient);
	  if (bd >= 0)
	    continue;   // Boundary curve, touch is not counted

	  sf_cv->ensureParCrvExistence(eps);
	  pcrv = sf_cv->parameterCurve();
	      
	  for (size_t kr=0; kr<side_cvs.size(); ++kr)
	    {
	      vector<pair<double,double> > int_pts;
	      vector<int> pretop;
	      vector<pair<pair<double,double>, pair<double,double> > > int_cvs;
	      intersect2Dcurves(pcrv.get(), side_cvs[kr].get(), eps,
				int_pts, pretop, int_cvs);
	      if (int_pts.size() > 0 || int_cvs.size() > 0)
		return 1;
	    }	
	}
      
      // Check if the trimming curve is completely inside the element
      // Only necessary for one curve in the loop. Use the last one
      if (pcrv.get())
	{
	  Point pt = pcrv->point(pcrv->startparam());
	  if (pt[0] >= elem_par[0] && pt[0] <= elem_par[1] &&
	      pt[1] >= elem_par[2] && pt[1] <= elem_par[3])
	    return 1;
	}
    }
  return 0;
}

//===========================================================================
int BoundedSurface::ElementBoundaryStatus(int elem_ix, double eps)
//===========================================================================
{
  // Result: -1 = not a spline surface, 0 = outside, 1 = on boundary, 2 = inside

  if (surface_->instanceType() != Class_SplineSurface)
    return -1;
  
  int bdstat = ElementOnBoundary(elem_ix, eps);
  if (bdstat != 0)
    return bdstat;

  // Fetch point internal to element
  // First fetch associated spline volume
  SplineSurface *sf = surface_->asSplineSurface();
  if (!sf)
    return -1;
  
  // Fetch number of patches in all parameter directions
  int nu = sf->numberOfPatches_u();
  int nv = sf->numberOfPatches_v();

  if (elem_ix < 0 || elem_ix >= nu*nv)
    return 0;

  // 2-variate index
  int iv = elem_ix/nu;
  int iu = elem_ix - iv*nu;

  // Parameter value
  vector<double> knots_u;
  vector<double> knots_v;
  const BsplineBasis basis_u = sf->basis_u();
  const BsplineBasis basis_v = sf->basis_v();
  basis_u.knotsSimple(knots_u);
  basis_v.knotsSimple(knots_v);

  double upar = 0.5*(knots_u[iu]+knots_u[iu+1]);
  double vpar = 0.5*(knots_v[iv]+knots_v[iv+1]);
  
  // Check if the element is inside the trimming loop(s)
  bool inside = inDomain(upar, vpar, eps);
  return (inside) ? 2 : 0;
}

//===========================================================================
Point BoundedSurface::getSurfaceParameter(int loop_idx, int cv_idx,  
					  double bd_par) const
//===========================================================================
{
  // Fetch boundary curve
  shared_ptr<ParamCurve> cv = (*boundary_loops_[loop_idx])[cv_idx];
  shared_ptr<CurveOnSurface> sf_cv = 
    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(cv);
  if (sf_cv.get() && sf_cv->parPref())
    {
      shared_ptr<ParamCurve> pcrv = sf_cv->parameterCurve();
      if (pcrv)
	{
	  Point ppar = pcrv->point(bd_par);
	  return ppar;
	}
    }

  // No parameter curve exist. Perform closest point computation
  Point pos = cv->point(bd_par);

  double clo_u, clo_v, clo_dist;
  double eps = 1.0e-6;
  Point clo_pt;
  double *seed = NULL;
  Point ppar;
  if (sf_cv.get() && sf_cv->hasParameterCurve())
    {
      ppar = sf_cv->parameterCurve()->point(bd_par);
      seed = ppar.begin();
    }
  closestPoint(pos, clo_u, clo_v, clo_pt, clo_dist, eps, NULL, seed);
  return Point(clo_u, clo_v);
}

//===========================================================================
bool BoundedSurface::simplifyBdLoops(double tol, double ang_tol, double& max_dist)
//===========================================================================
{
  box_.unset();

  max_dist = 0;
  double dist;
  bool modified = false;
  for (size_t ki=0; ki<boundary_loops_.size(); ++ki)
    {
      dist = 0;
      bool curr_modified = boundary_loops_[ki]->simplify(tol, ang_tol, dist);

      if (curr_modified)
	modified = true;
      max_dist = std::max(max_dist, dist);
    }

  if (modified)
    {
#ifdef DEBUG
  std::ofstream of("simplified_bd.g2");
  writeStandardHeader(of);
  write(of);
#endif
    }
  return modified;
}


//===========================================================================
bool BoundedSurface::makeUnderlyingSpline()
//===========================================================================
{
  shared_ptr<SplineSurface> spl_surf = dynamic_pointer_cast<SplineSurface, ParamSurface>(surface_);
  if (spl_surf.get() != 0)
    // Alredy spline
    return true;

  shared_ptr<ElementarySurface> elem_surf = dynamic_pointer_cast<ElementarySurface, ParamSurface>(surface_);
  if (elem_surf.get())
    surface_ = shared_ptr<ParamSurface>(elem_surf->geometrySurface());
  else
    {
      shared_ptr<ParamSurface> spl_surf2 = 
	shared_ptr<ParamSurface>(surface_->asSplineSurface());
      if (!spl_surf2.get())
	return false;
      surface_  = spl_surf2;
    }
  for (int i = 0; i < (int)boundary_loops_.size(); ++i)
    {
      shared_ptr<CurveLoop> cl = boundary_loops_[i];
      
      for (vector<shared_ptr<ParamCurve> >::iterator it = cl->begin();
	   it != cl->end();
	   ++it)
	{
	  shared_ptr<CurveOnSurface> curve = dynamic_pointer_cast<CurveOnSurface, ParamCurve>(*it);
	  if (curve.get() != 0)
	    curve->setUnderlyingSurface(surface_);
	}
    }
  return true;
}


//===========================================================================
bool BoundedSurface::allIsSpline() const
//===========================================================================
{
  shared_ptr<ParamSurface> underlying = surface_;
  shared_ptr<BoundedSurface> under_bound = dynamic_pointer_cast<BoundedSurface, ParamSurface>(underlying);
  while (under_bound.get() != NULL)
    {
      underlying = under_bound->underlyingSurface();
      under_bound = dynamic_pointer_cast<BoundedSurface, ParamSurface>(underlying);
    }

  if (underlying->instanceType() != Class_SplineSurface)
    return false;

  for (size_t i = 0; i < boundary_loops_.size(); ++i)
    {
      int j = 0;
      for (vector<shared_ptr<ParamCurve> >::const_iterator it = boundary_loops_[i]->begin();
	   it != boundary_loops_[i]->end();
	   ++it, ++j)
	{
	  shared_ptr<CurveOnSurface> cos = dynamic_pointer_cast<CurveOnSurface, ParamCurve>(*it);
	  if (cos.get() == NULL)
	    return false;
	  shared_ptr<ParamCurve> pcrv;
	  if (cos->parPref())
	    pcrv = cos->parameterCurve();
	  else
	    pcrv = cos->spaceCurve();
	  if (pcrv->instanceType() != Class_SplineCurve)
	    return false;
	}
    }
  return true;
}


//===========================================================================
BoundedSurface* BoundedSurface::allSplineCopy() const
//===========================================================================
{
  shared_ptr<ParamSurface> underlying = surface_;
  shared_ptr<BoundedSurface> under_bound = dynamic_pointer_cast<BoundedSurface, ParamSurface>(underlying);
  while (under_bound.get() != NULL)
    {
      underlying = under_bound->underlyingSurface();
      under_bound = dynamic_pointer_cast<BoundedSurface, ParamSurface>(underlying);
    }

  shared_ptr<SplineSurface> underlying_spline;
  if (underlying->instanceType() == Class_SplineSurface)
    {
      shared_ptr<SplineSurface> underlying_spl_case = dynamic_pointer_cast<SplineSurface, ParamSurface>(underlying);
      underlying_spline = shared_ptr<SplineSurface>(underlying_spl_case->clone());
    }
  else
    {
      shared_ptr<ElementarySurface> underlying_el = dynamic_pointer_cast<ElementarySurface, ParamSurface>(underlying);
      if (underlying_el.get() == NULL)
	return NULL;
      underlying_spline = shared_ptr<SplineSurface>(underlying_el->geometrySurface());
    }

  int nmb_loops = (int)boundary_loops_.size();
  vector<double> space_epsilons(nmb_loops);
  vector<vector<shared_ptr<CurveOnSurface> > > spline_loops(nmb_loops);

  for (int i = 0; i < nmb_loops; ++i)
    {
      space_epsilons[i] = boundary_loops_[i]->getSpaceEpsilon();
      spline_loops[i].resize(boundary_loops_[i]->size());
      int j = 0;
      for (vector<shared_ptr<ParamCurve> >::const_iterator it = boundary_loops_[i]->begin();
	   it != boundary_loops_[i]->end();
	   ++it, ++j)
	{
	  shared_ptr<CurveOnSurface> cos = dynamic_pointer_cast<CurveOnSurface, ParamCurve>(*it);
	  if (cos.get() == NULL)
        return NULL;
	  shared_ptr<ParamCurve> pcrv;
	  bool par_pref = cos->parPref();
	  if (par_pref)
	    pcrv = cos->parameterCurve();
	  else
	    pcrv = cos->spaceCurve();
	  shared_ptr<SplineCurve> curve_spline;
	  if (pcrv->instanceType() == Class_SplineCurve)
	    {
	      shared_ptr<SplineCurve> pcrv_spl = dynamic_pointer_cast<SplineCurve, ParamCurve>(pcrv);
	      curve_spline = shared_ptr<SplineCurve>(pcrv_spl->clone());
	    }
	  else
	    curve_spline = shared_ptr<SplineCurve>(pcrv->geometryCurve());
	  spline_loops[i][j] = shared_ptr<CurveOnSurface> (new CurveOnSurface(underlying_spline, curve_spline, par_pref));
	}
    }

  return new BoundedSurface(underlying_spline, spline_loops, space_epsilons);
}

//===========================================================================
bool BoundedSurface::isPlanar(Point& normal, double tol)
//===========================================================================
{
  return surface_->isPlanar(normal, tol);
}

//===========================================================================
void BoundedSurface::estimateSfSize(double& u_size, double& v_size, int u_nmb,
				    int v_nmb) const
//===========================================================================
{
  if (nmb_size_u_ >= u_nmb && nmb_size_v_ >= v_nmb)
    {
      u_size = est_sf_size_u_;
      v_size = est_sf_size_v_;
    }
  else
    {
      RectDomain dom = containingDomain();
      double del_u = (dom.umax() - dom.umin())/(double)(u_nmb);
      double del_v = (dom.vmax() - dom.vmin())/(double)(v_nmb);
      double u_par, v_par;

      int ki;
      u_size = v_size = 0.0;
      int nmb = 0;
      for (ki=0, v_par=dom.vmin()+0.5*del_v; ki<v_nmb; ++ki, v_par+=del_v)
	{
	  vector<shared_ptr<ParamCurve> > cvs;
	  try {
	    cvs = constParamCurves(v_par, true);
	  }
	  catch (...)
	    {
	      continue;
	    }
	  double len = 0.0;
	  for (size_t kj=0; kj<cvs.size(); ++kj)
	    len += cvs[kj]->estimatedCurveLength();
	  u_size += len;
	  nmb++;
	}
      u_size /= (double)nmb;

      nmb = 0;
       for (ki=0, u_par=dom.umin()+0.5*del_u; ki<u_nmb; ++ki, u_par+=del_u)
	{
	  vector<shared_ptr<ParamCurve> > cvs;
	  try {
	    cvs = constParamCurves(u_par, false);
	  }
	  catch (...)
	    {
	      continue;
	    }
	  double len = 0.0;
	  for (size_t kj=0; kj<cvs.size(); ++kj)
	    len += cvs[kj]->estimatedCurveLength();
	  v_size += len;
	  nmb++;
	}
      v_size /= (double)nmb;

      est_sf_size_u_ = u_size;
      est_sf_size_v_ = v_size;
      nmb_size_u_ = u_nmb;
      nmb_size_v_ = v_nmb;
   }
}

//===========================================================================
void BoundedSurface::estimateSfSize(double& u_size, double& min_u, 
				    double& max_u, double& v_size, 
				    double& min_v, double& max_v,
				    int u_nmb, int v_nmb) const
//===========================================================================
{
  RectDomain dom = containingDomain();
  double del_u = (dom.umax() - dom.umin())/(double)(u_nmb-1);
  double del_v = (dom.vmax() - dom.vmin())/(double)(v_nmb-1);
  double u_par, v_par;

  int ki;
  u_size = v_size = max_u = max_v = 0.0;
  min_u = min_v = std::numeric_limits<double>::max();
  int nmb = 0;
  for (ki=0, v_par=dom.vmin(); ki<v_nmb; ++ki, v_par+=del_v)
    {
      vector<shared_ptr<ParamCurve> > cvs;
      try {
	cvs = constParamCurves(v_par, true);
      }
      catch (...)
	{
	  continue;
	}
      double len = 0.0;
      for (size_t kj=0; kj<cvs.size(); ++kj)
	len += cvs[kj]->estimatedCurveLength();
      v_size += len;
      nmb++;
      min_v = std::min(min_v, len);
      max_v = std::max(max_v, len);
    }
  v_size /= (double)nmb;

  nmb = 0;
  for (ki=0, u_par=dom.umin(); ki<u_nmb; ++ki, u_par+=del_u)
    {
      vector<shared_ptr<ParamCurve> > cvs;
      try {
	cvs = constParamCurves(u_par, false);
      }
      catch (...)
	{
	  continue;
	}
      double len = 0.0;
      for (size_t kj=0; kj<cvs.size(); ++kj)
	len += cvs[kj]->estimatedCurveLength();
      u_size += len;
      nmb++;
      min_u = std::min(min_u, len);
      max_u = std::max(max_u, len);
    }
  u_size /= (double)nmb;
}

//===========================================================================
bool BoundedSurface::isLinear(Point& dir1, Point& dir2, double tol)
//===========================================================================
{
  // No holes are allowed
  if (boundary_loops_.size() > 1)
    return false;

  // Check the underlying surface
  bool sf_linear = surface_->isLinear(dir1, dir2, tol);
  if (!sf_linear)
    return false;

  // Check if the configuration of the trimming curves is consistent with
  // the linearity. Each trimming curve must be linear and either parallel or 
  // perpendicular to the direction of linearity of the surface, or it must lie
  // in a plane with normal equal to the direction of linearity 
  int nmb = boundary_loops_[0]->size();
  for (int ki=0; ki<nmb; ++ki)
    {
      shared_ptr<ParamCurve> cv = (*boundary_loops_[0])[ki];
      DirectionCone cone = cv->directionCone();

      Point pos;
      bool planar = cv->isInPlane(dir1, tol, pos);
      if (planar)
	continue;

      if (dir2.dimension() > 0)
	{
	  planar = cv->isInPlane(dir2, tol, pos);
	  if (planar)
	    {
	      dir1 = dir2;
	      dir2.resize(0);
	      continue;
	    }
	}
      
      if (cone.angle() > tol)
	return false;
      
      double ang = cone.centre().angle(dir1);
      double ang2 = 2.0*M_PI;
      if (dir2.dimension() > 0)
	ang2 = cone.centre().angle(dir2);
      if (ang2 > tol && fabs(M_PI-ang2) > tol && fabs(0.5*M_PI-ang2) > tol)
	dir2.resize(0);  // At most one direction of linearity
      if (ang > tol && fabs(M_PI-ang) > tol && fabs(0.5*M_PI-ang2) > tol)
	{
	  if (dir2.dimension() > 0)
	    dir1 = dir2;
	  else
	    return false;
	}
    }
  return true;
}

//===========================================================================
bool BoundedSurface::isAxisRotational(Point& centre, Point& axis, Point& vec,
				     double& angle)
//===========================================================================
{
  double eps = boundary_loops_[0]->getSpaceEpsilon();
  double teps = 1.0e-12;  // Just for equality testing

  // Check underlying surface
  bool rotational = surface_->isAxisRotational(centre, axis, vec, angle);
  Point normal;
  bool planar = surface_->isPlanar(normal, eps);
  if ((!rotational) && (!planar))
    return false;

  // Check trimming curves
  double curve_ang = -1.0;
  if (rotational)
    {
      // All trimming curves must either have the same centre and axis
      // as the underlying surface or lie in a plane going through the
      // axis. At most two such planes are allowed
      if (boundary_loops_.size() > 1)
	return false;

      vector<Point> plane_norm(2);

      vector<vector<shared_ptr<ParamCurve> > > smooth_cvs;
      boundary_loops_[0]->getSmoothCurves(smooth_cvs, eps);
      for (size_t kr=0; kr<smooth_cvs.size(); ++kr)
	{
	  Point curr_vec;
	  double curr_ang = 0.0;
	  for (size_t kj=0; kj<smooth_cvs[kr].size(); ++kj)
	    {
	      Point centre2, axis2, vec2, dir;
	      double angle2;
	      bool rot = smooth_cvs[kr][kj]->isAxisRotational(centre2, axis2, vec2, angle2);
	      Point tmp_norm;
	      bool planar_cv = smooth_cvs[kr][kj]->isInPlane(centre, axis, 
							     eps, tmp_norm);
	      if (planar_cv)
		{
		  // Check the number of planes identified
		  int kp;
		  for (kp=0; kp<2; ++kp)
		    {
		      if (plane_norm[kp].dimension() == 0)
			{
			  plane_norm[kp] = tmp_norm;
			  break;
			}
		      else
			{
			  double norm_ang = plane_norm[kp].angle(tmp_norm);
			  if (norm_ang < eps || fabs(M_PI-norm_ang) < eps)
			    break;  // Questionable use of tolerance
			}
		    }
		  if (kp == 2)
		    return false;  // More than two planes
		}
	      else if (rot)
		{
		  double tmp_ang = axis.angle(axis2);
		  if (tmp_ang > eps && fabs(M_PI-tmp_ang) > eps)
		    return false;

		  if (tmp_ang > eps)
		    {
		      // The axes are oppositely oriented. Adjust the start vector
		      Array<double,3> tmp_vec(vec2[0], vec2[1], vec2[2]);
		      MatrixXD<double, 3> mat;
		      mat.setToRotation(angle2, axis2[0], axis2[1], axis2[2]);  // Rotate the 
		      // start vector the angle angle2 around axis2
		      Array<double,3> tmp_vec2 = mat*tmp_vec;
		      vec2 = Point(tmp_vec2[0], tmp_vec2[1], tmp_vec2[2]);
		    }

		  if (curr_vec.dimension() == 0)
		    {
		      curr_vec = vec2;
		      curr_ang = angle2;
		    }
		  else
		    {
		      // @@@ VSK. This logic may need some more work
		      double vec_ang = curr_vec.angle(vec2);
		      double vec_ang2 = vec.angle(vec2);
		      if (fabs(vec_ang-angle2) < eps && vec_ang2 < eps)
		      	curr_vec = vec2;
		      curr_ang += angle2;
		    }

		  Point tmp_vec = centre2 - centre;
		  tmp_ang = axis.angle(tmp_vec);
		  if (tmp_vec.length() > eps && tmp_ang > eps &&
		      fabs(M_PI-tmp_ang) > eps)
		    return false;
		}
	      else
		return false;
	    }
	  if (curr_vec.dimension() == vec.dimension())
	    {
	      if (curve_ang < 0)
		curve_ang = curr_ang;
	      else if (fabs(curve_ang - curr_ang) > eps)
		return false;  // Rotational angle varies around the loop

	      if (curr_ang < angle+teps)
		{
		  vec = curr_vec;
		  angle = curr_ang;
		}
	    }
	}
    }     
  else if (planar)
    {
      // All trimming curves must have the same axis which is the
      // same as the plane normal and the same centre
      // Since the trimming loops are closed and all trimming curves must
      // be circles, the rotation must be complete
      axis = normal;
      for (size_t ki=0; ki<boundary_loops_.size(); ++ki)
	{
	  vector<vector<shared_ptr<ParamCurve> > > smooth_cvs;
	  boundary_loops_[ki]->getSmoothCurves(smooth_cvs, eps);
	  if (smooth_cvs.size() > 1)
	    return false;  // Not rotational
	  for (size_t kj=0; kj<smooth_cvs[0].size(); ++kj)
	    {
	      Point centre2, axis2, vec2;
	      double angle2;
	      bool rot = smooth_cvs[0][kj]->isAxisRotational(centre2, axis2, vec2, angle2);
	      if (!rot)
		return false;
	      double tmp_ang = axis.angle(axis2);
	      if (tmp_ang > eps && fabs(M_PI-tmp_ang) > eps)
		return false;
	      if (centre.dimension() == centre2.dimension())
		{
		  Point tmp_vec = centre2 - centre;
		  tmp_ang = axis.angle(tmp_vec);
		  if (tmp_vec.length() > eps && tmp_ang > eps && 
		      fabs(M_PI-tmp_ang) > eps)
		    return false;
		}
	      else
		centre = centre2;
	    }
	}
      if (centre.dimension() == 0)
	return false;

      Point pt = (*boundary_loops_[0])[0]->point((*boundary_loops_[0])[0]->startparam());
      vec = pt - centre;
      vec.normalize();
      angle = 2.0*M_PI;
    }
  return true;
}

//===========================================================================
void BoundedSurface:: replaceSurf(shared_ptr<ParamSurface> sf)
//===========================================================================
{
  // Update pointers to surface 
  for (size_t ki=0; ki<boundary_loops_.size(); ++ki)
    {
      int nmb = boundary_loops_[ki]->size();
      for (int kr=0; kr<nmb; ++kr)
	{
	  shared_ptr<ParamCurve> cv = (*boundary_loops_[ki])[kr];
	  shared_ptr<CurveOnSurface> sf_cv = 
	    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(cv);
	  if (sf_cv.get())
	    sf_cv->setUnderlyingSurface(sf);
	}
    }

  surface_ = sf;
}
