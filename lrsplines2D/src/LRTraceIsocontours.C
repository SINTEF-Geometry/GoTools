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

#include <fstream> // @@ debug purpose
#include <functional>
#include "GoTools/lrsplines2D/LRTraceIsocontours.h"
#include "GoTools/lrsplines2D/SSurfTraceIsocontours.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/CurveBoundedDomain.h"
#include "GoTools/lrsplines2D/LRSplinePlotUtils.h" // debug
#include "GoTools/geometry/RectDomain.h"
#include "GoTools/geometry/SplineDebugUtils.h"
#include "GoTools/lrsplines2D/SSurfTraceIsocontours.h"
#include "GoTools/creators/TrimCrvUtils.h"

//#define DEBUG0
//#define DEBUG
//#define DEBUG2

using namespace std;
using namespace Go;

using LRSurfPtr = shared_ptr<LRSplineSurface>;
using IsectCurve = pair<CurvePtr, CurvePtr> ; // @@
namespace {

// ----------------------------------------------------------------------------
// Wrapper to simplify syntax when the shared pointer itself is not required
// except to manage the allocated memory
inline shared_ptr<SplineSurface> as_spline_surf(const shared_ptr<LRSplineSurface> lrs)
// ----------------------------------------------------------------------------
{
  return shared_ptr<SplineSurface>(lrs->asSplineSurface());
}

// ----------------------------------------------------------------------------  
vector<CurveVec> 
merge_isocontours(vector<vector<CurveVec>>& curve_fragments,
		  const vector<pair<LRSurfPtr,LRSplineSurface::PatchStatus> > surf_frags,
		  const LRSplineSurface& lrs,
		  const std::vector<double>& isovals,
		  const double tol,
		  const CurveBoundedDomain* domain);
// ----------------------------------------------------------------------------

  
};// end anonymous namespace


namespace Go
{

// ============================================================================
  vector<CurveVec> LRTraceIsocontours(const shared_ptr<ParamSurface>& surf,
				      const std::vector<double>& isovals,
				      const int threshold_missing,
				      const double tol)
// ============================================================================
{
  shared_ptr<LRSplineSurface> lrsurf = 
    dynamic_pointer_cast<LRSplineSurface, ParamSurface>(surf);
  shared_ptr<BoundedSurface> bdsurf = 
    dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);
  CurveBoundedDomain bddom;

#ifdef DEBUG
  if (bdsurf.get())
    {
      std::ofstream of0("trans_trim_init.g2");
      vector<CurveLoop> bd_loops = bdsurf->allBoundaryLoops();
      for (size_t ki=0; ki<bd_loops.size(); ++ki)
	{
	  int nmb = bd_loops[ki].size();
	  for (int kj=0; kj<nmb; ++kj)
	    {
	      shared_ptr<ParamCurve> cv = bd_loops[ki][kj];
	      shared_ptr<CurveOnSurface> bdcv =
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv);
	      if (bdcv.get())
		cv = bdcv->parameterCurve();
	      cv->writeStandardHeader(of0);
	      cv->write(of0);
	    }
	}
    }
#endif

  // Translate to parameter domain to origo
  RectDomain dom = surf->containingDomain();
  Point vec(-0.5*(dom.umin()+dom.umax()), -0.5*(dom.vmin()+dom.vmax()));
  TrimCrvUtils::translateSurfaceDomain(surf.get(), vec);

#ifdef DEBUG
  std::ofstream of("trans_trim.g2");
#endif
  if (bdsurf.get())
    {
      lrsurf = dynamic_pointer_cast<LRSplineSurface, ParamSurface>(bdsurf->underlyingSurface());
      bddom = bdsurf->parameterDomain();

#ifdef DEBUG
      vector<CurveLoop> bd_loops = bdsurf->allBoundaryLoops();
      for (size_t ki=0; ki<bd_loops.size(); ++ki)
	{
	  int nmb = bd_loops[ki].size();
	  for (int kj=0; kj<nmb; ++kj)
	    {
	      shared_ptr<ParamCurve> cv = bd_loops[ki][kj];
	      shared_ptr<CurveOnSurface> bdcv =
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv);
	      if (bdcv.get())
		cv = bdcv->parameterCurve();
	      cv->writeStandardHeader(of);
	      cv->write(of);
	    }
	}

      std::ofstream ofbd("sf_loop.g2");
      SplineDebugUtils::writeBoundary(*bdsurf, ofbd);
#endif
    }

  if (!lrsurf.get())
    {
      THROW("No LR B-spline surface is found");
    }
    
  CurveBoundedDomain *bddomptr = NULL;
  if (bddom.nmbLoops() > 0)
    bddomptr = &bddom;

  const bool include_3D = true;
  bool use_sisl_marching = false;

  const vector<CurveVec> curves = LRTraceIsocontours(*lrsurf,
					       isovals,
					       threshold_missing,
					       tol,
					       include_3D,
					       use_sisl_marching,
					       bddomptr);

  // Translate resulting curves back to the original parameter domain
  vec *= -1;
  Point vec2(vec[0], vec[1], 0.0);
  vector<CurveVec> curves2(curves.size());
  for (size_t ki=0; ki<curves.size(); ++ki)
    {
      curves2[ki].resize(curves[ki].size());
      for (size_t kj=0; kj<curves[ki].size(); ++kj)
  	{
  	  auto cv = curves[ki][kj];
  	  if (cv.first.get())
  	    {
	      shared_ptr<SplineCurve> cv2_1(cv.first->clone());
	      shared_ptr<SplineCurve> cv2_2(cv.second->clone());
  	      cv2_1->translateCurve(vec);
  	      cv2_2->translateCurve(vec2);
	      curves2[ki][kj] = make_pair(cv2_1, cv2_2);
  	    }
  	}
    }

  return curves2;
}

// ============================================================================
vector<CurveVec> LRTraceIsocontours(const LRSplineSurface& lrs,
				    const std::vector<double>& isovals,
				    const int threshold_missing,
				    const double tol,
				    const bool include_3D_curves,
				    const bool use_sisl_marching,
				    const CurveBoundedDomain* domain)
// ============================================================================
{
#ifdef DEBUG2
  std::ofstream mesh0("mesh_init.eps");
  writePostscriptMesh(lrs, mesh0, true);
#endif
  
  //break LR surface up into individual patches
  const vector<pair<LRSurfPtr,LRSplineSurface::PatchStatus> > surf_fragments = 
    lrs.subdivideIntoSimpler(threshold_missing, tol, domain);

#ifdef DEBUG0
  std::ofstream of("lrsurf_fragments.g2");
  for (size_t ki=0; ki<surf_fragments.size(); ++ki)
    {
      shared_ptr<LRSplineSurface> tmp0 = surf_fragments[ki].first;
      tmp0->writeStandardHeader(of);
      tmp0->write(of);
      of << std::endl;
    }

  std::ofstream of2("lrsurf_fragments2.g2");
  for (size_t ki=0; ki<surf_fragments.size(); ++ki)
    {
      if (surf_fragments[ki].second == LRSplineSurface::INTERSECT)
	{
	  shared_ptr<LRSplineSurface> tmp0 = surf_fragments[ki].first;
	  tmp0->writeStandardHeader(of2);
	  tmp0->write(of2);
	  of2 << std::endl;
	}
    }
#endif

#ifdef DEBUG2
  std::cout << "Number of surface fragments: " << surf_fragments.size() << std::endl;
  vector<shared_ptr<LRSplineSurface> > lr1(surf_fragments.size());
  vector<shared_ptr<LRSplineSurface> > lr2(surf_fragments.size());
  for (size_t ki=0; ki<surf_fragments.size(); ++ki)
    {
      lr1[ki] = surf_fragments[ki].first;
      lr2[ki] = shared_ptr<LRSplineSurface>(surf_fragments[ki].first->clone());
      lr2[ki]->expandToFullTensorProduct();
    }

  std::ofstream mesh1("mesh_lr.eps");
  std::ofstream mesh2("mesh_tp.eps");
  writePostscriptMesh(lrs, lr1, mesh1, true);
  writePostscriptMesh(lrs, lr2, mesh2, true);
#endif


  // // Setting up function to compute isocontours for a surface fragment
  // const function<vector<CurveVec>(LRSurfPtr)> compute_isovals {
  //   [&] (pair<LRSurfPtr,LRSplineSurface::PatchStatus> l)
  //     {
  // 	LRSurfPtr lrsurf = l.first;
  // 	LRSplineSurface::PatchStatus stat = l.second;
  // 	if (stat == LRSplineSurface::OUTSIDE)
  // 	  {
  // 	    vector<CurveVec> dummy;
  // 	    return dummy;
  // 	  }
  // 	else
  // 	  return SSurfTraceIsocontours(*as_spline_surf(lrsurf), isovals,
  // 				       tol, include_3D_curves,
  // 				       use_sisl_marching);
  //     }
  // };

  // // computing isocurves for each surface fragment (vector<vector<CurveVec>>)
  // const auto curve_fragments = apply_transform(surf_fragments, compute_isovals);

  vector<vector<CurveVec> > curve_fragments(surf_fragments.size());
  for (size_t ki=0; ki<surf_fragments.size(); ++ki)
    {
      LRSplineSurface::PatchStatus stat = surf_fragments[ki].second;
      if (stat == LRSplineSurface::OUTSIDE)
	{
	  vector<CurveVec> dummy(isovals.size());
	  curve_fragments[ki] = dummy;
	}
      else
	try {
	  curve_fragments[ki] =
	    SSurfTraceIsocontours(*as_spline_surf(surf_fragments[ki].first), 
				  isovals,
				  tol, include_3D_curves,
				  use_sisl_marching);
	}
	catch (...)
	  {
	    std::cout << "Tracing of curve " << ki << "failed" << std::endl;
	  }
    }

#ifdef DEBUG0
  std::cout << "Ready to merge isocontours" << std::endl;
  ofstream os_surf("lrsurf.g2");
  for (auto f : surf_fragments) {
    f.first->writeStandardHeader(os_surf);
    f.first->write(os_surf);
  }
  os_surf.close();

  ofstream os_cv("curvefrags.g2");
  for (auto frag : curve_fragments)
    for (auto ival : frag)
      for (auto c : ival) {
	if (c.second.get())
	  {
	    c.second->writeStandardHeader(os_cv);
	    c.second->write(os_cv);
	  }
      }
  os_cv.close();

#endif
  
  // trim with trimming loop if existing, merge isocontours across 
  // patches and return result
  return merge_isocontours(curve_fragments, surf_fragments, lrs, isovals, 
			   tol, domain);
}

}; // end namespace Go;

namespace {

// ----------------------------------------------------------------------------
  array<double, 4> parameter_domain(const LRSurfPtr& patch)
// ----------------------------------------------------------------------------
{
  return array<double, 4>{ { patch->paramMin(XFIXED),
	                     patch->paramMax(XFIXED), 
	                     patch->paramMin(YFIXED),
	                     patch->paramMax(YFIXED)} };
}


// ----------------------------------------------------------------------------
  int find_exit_edge(const IsectCurve& c, const array<double, 4>& domain,
		     const bool start, const double tol, int edge)
// ----------------------------------------------------------------------------
{
  const double t = start ? c.first->startparam() : c.first->endparam();
  Point uv;  // represent u and v parameters
  c.first->point(uv, t);  // evaluate (u, v) parameter pair at curve start or end

  const vector<double> dists = {fabs(domain[0] - uv[0]),
				fabs(domain[1] - uv[0]),
				fabs(domain[2] - uv[1]),
				fabs(domain[3] - uv[1])};
  int first = std::max(0, edge + 1);
  if (first >= (int)dists.size())
    return -1;
  auto min_it = min_element(dists.begin()+first, dists.end());
  if (*min_it > tol)
    return -1; // curve doesn't end on edge, but inside domain (either closed, or ends at
	       // singularity)

  return int(min_it - dists.begin()); // index of boundary edge where curve starts/ends
		 
}

// ----------------------------------------------------------------------------
void map_curve(const IsectCurve& c, const array<double, 4>& domain, const double tol,
	       const bool at_start, map<double, CurveVec>& u_map, map<double, CurveVec>& v_map)
// ----------------------------------------------------------------------------
{
  int edge = -1;
  while (true)
    {
      edge = find_exit_edge(c, domain, at_start, tol, edge);

      if (edge >=0) {
	Point p1; c.first->point(p1, at_start ? c.first->startparam() : c.first->endparam());
	//assert( (fabs(domain[edge] - p1[0]) < tol) | (fabs(domain[edge] - p1[1]) < tol) );
	if ( !((fabs(domain[edge] - p1[0]) < tol) | (fabs(domain[edge] - p1[1]) < tol) ))
	  {
	    std::cout << "Error in map curve" << std::endl;
	    THROW("Error in map_curve");
	  }
      }

      if (edge < 0 || edge > 3)
	break;

      // if none of the cases below apply, curve does not terminate at edge.
      switch (edge) {
      case 0:
      case 1:
	u_map[domain[edge]].push_back(c);
	break;
      case 2:
      case 3:
	v_map[domain[edge]].push_back(c);
	break;
      }
    }
}    

  // ----------------------------------------------------------------------------
void prepare_curvemaps(const array<double, 4>& domain, const double tol, const CurveVec& cvec,
		       map<double, CurveVec>& u_map, map<double, CurveVec>& v_map)
// ----------------------------------------------------------------------------
{
  for_each(cvec.begin(), cvec.end(), [&] (const IsectCurve& c) {
      if (c.first.get())
	{
	  map_curve(c, domain, tol, true, u_map, v_map);  // map curve according to start point
	  map_curve(c, domain, tol, false, u_map, v_map); // map curve according to end point
	}
  });
}

// ----------------------------------------------------------------------------
IsectCurve join_isectcurves(const IsectCurve& c1, const IsectCurve& c2,
			    bool c1_at_start, bool c2_at_start)
// ----------------------------------------------------------------------------
{
  shared_ptr<SplineCurve> pcurve1 = shared_ptr<SplineCurve>(c1.first->clone());
  shared_ptr<SplineCurve> scurve1; 
  shared_ptr<SplineCurve> pcurve2 = shared_ptr<SplineCurve>(c2.first->clone());
  shared_ptr<SplineCurve> scurve2; 

  if (c1.second.get())
    {
      scurve1 = shared_ptr<SplineCurve>(c1.second->clone());
      scurve2 = shared_ptr<SplineCurve>(c2.second->clone());
      if (c1_at_start) {
	pcurve1->reverseParameterDirection();
	scurve1->reverseParameterDirection();
      }

      if (!c2_at_start) {
	pcurve2->reverseParameterDirection();
	scurve2->reverseParameterDirection();
      }
  
      pcurve1->appendCurve(pcurve2.get());
      scurve1->appendCurve(scurve2.get());
    }
  else
    {
      if (c1_at_start) {
	pcurve1->reverseParameterDirection();
      }

      if (!c2_at_start) {
	pcurve2->reverseParameterDirection();
      }
  
      pcurve1->appendCurve(pcurve2.get());
    }
  return IsectCurve { pcurve1, scurve1 };
}

// ----------------------------------------------------------------------------
void replace_segments(const IsectCurve& old1, const IsectCurve& old2,
		      const IsectCurve& updated, map<double, CurveVec>& target)
// ----------------------------------------------------------------------------
{
  for (auto& it_map : target) 
    for (auto& it_vec : it_map.second) 
      if ((it_vec.first == old1.first) || (it_vec.first == old2.first)) 
	it_vec = updated;
}


// ----------------------------------------------------------------------------
pair<double, double> identify_truncated_endpoints(const IsectCurve& c, Direction2D d,
						  double pval, const double tol)
// ----------------------------------------------------------------------------
{
  Point startpoint, endpoint;
  c.first->point(startpoint, c.first->startparam());
  c.first->point(endpoint, c.first->endparam());

  double len = c.first->estimatedCurveLength(3);
  const int ix = (d==YFIXED) ? 1 : 0;
  auto NaN = numeric_limits<double>::quiet_NaN();
#ifdef DEBUG
  double d1 = fabs(startpoint[ix] - pval);
  double d2 = fabs(endpoint[ix] - pval);
  if (d1 < tol && d2 < tol)
    {
      int stop_break = 1;
    }
#endif
  return pair<double, double> {
    (fabs(startpoint[ix] - pval) < tol && len >= tol) ? startpoint[(ix+1)%2] : NaN,
    (fabs(endpoint[ix] - pval) < tol && len >= tol) ? endpoint[(ix+1)%2] : NaN
  };
  
  // if (d==YFIXED) {
  //   swap(startpoint[0], startpoint[1]);
  //   swap(endpoint[0], endpoint[1]);
  // }
  // auto NaN = numeric_limits<double>::quiet_NaN();
  // return pair<double, double> {
  //   (fabs(startpoint[0] - pval) < tol ? startpoint[1] : NaN),
  //   (fabs(endpoint[0] - pval) < tol ? endpoint[1] : NaN)
  // };
}


// ----------------------------------------------------------------------------
void merge_segments(map<double, CurveVec>& mergemap, // map whose segments should be merged
		    map<double, CurveVec>& othermap, // map whose entries must also be updated
		    Direction2D d,                   // the fixed parameter direction
		    const LRSplineSurface& lrs,
		    const double isoval,
		    const double tol,               
		    CurveVec& bcurves,            
		    CurveVec& finished_curves) // insert newly merged finished curves here
// ----------------------------------------------------------------------------
{
  struct EndPoint {double pval; IsectCurve icurve;bool at_start;};

  array<double, 4> domain{{lrs.paramMin(XFIXED),
	lrs.paramMax(XFIXED), lrs.paramMin(YFIXED), lrs.paramMax(YFIXED)}};

  while (!mergemap.empty()) {
    auto it = *mergemap.begin();   mergemap.erase(mergemap.begin());
    double min_len = std::numeric_limits<double>::max();
    for (auto& ic : it.second) 
      { 
	min_len = std::min(min_len, ic.first->estimatedCurveLength(3));
      }
    double tol2 = std::max(1.0e-4, std::min(0.5*min_len, tol));

	// loop over intersection curve segments cut by this parameter line to set tolerance
    // sorting incident curve segments so that those that should be merged will lie right 
    // next to each other
    vector<CurvePtr> encountered;
    vector<EndPoint> tp_vec;
    for (auto& ic : it.second) { // loop over intersection curve segments cut by this parameter line
      const auto ends = identify_truncated_endpoints(ic, d, it.first, tol2);
      if (!std::isnan(ends.first + ends.second)) {
	// both endpoints of this curve are truncated by the parameter line.
	if (find(encountered.begin(), encountered.end(), ic.first) != encountered.end())
	  continue;
	encountered.push_back(ic.first); // ensure we will only register the curve once
      }
      if (!std::isnan(ends.first))
	{
	  if (find_exit_edge(ic, domain, true, tol, -1) < 0)
	    tp_vec.push_back({ends.first, ic, true});
	}
      if (!std::isnan(ends.second)) 
	{
	  if (find_exit_edge(ic, domain, false, tol, -1) < 0)
	    tp_vec.push_back({ends.second, ic, false});
	}
    }

#ifdef DEBUG
	std::ofstream of0("merge_curves.g2");
	for (size_t ka=0; ka<tp_vec.size(); ++ka)
	  {
	    tp_vec[ka].icurve.first->writeStandardHeader(of0);
	    tp_vec[ka].icurve.first->write(of0);
	  }
#endif
    // if (tp_vec.size() % 2 != 0)
    //   {
	// Check for duplicates
	for (size_t ki=0; ki<tp_vec.size(); ++ki)
	  for (size_t kj=ki+1; kj<tp_vec.size(); )
	    {
	      if (tp_vec[ki].icurve.first.get() == tp_vec[kj].icurve.first.get() &&
		  tp_vec[ki].at_start == tp_vec[kj].at_start)
		tp_vec.erase(tp_vec.begin()+kj);
	      else
		++kj;
	    }

	if (tp_vec.size() % 2 != 0)
	  {
	    // we suppose that for each curve going in, there is one going out
	    MESSAGE("Inconsistent number of curve endpoints");
	  }
      // }

    sort(tp_vec.begin(), tp_vec.end(), [](const EndPoint& t1, const EndPoint& t2) {return t1.pval < t2.pval;});

    //for (size_t i = 0; i != tp_vec.size(); i += 2) {
    size_t nmb = 2*(tp_vec.size()/2);
    for (size_t i = 0; i < nmb; ) 
      {
	if (i+1 >= tp_vec.size())
	  break;
	const auto entry1 = tp_vec[i];
	const auto entry2 = tp_vec[i+1];

#ifdef DEBUG
	std::ofstream of1("curves_to_merge.g2");
	entry1.icurve.first->writeStandardHeader(of1);
	entry1.icurve.first->write(of1);
	entry2.icurve.first->writeStandardHeader(of1);
	entry2.icurve.first->write(of1);
#endif
	//assert(fabs(entry1.pval - entry2.pval) < tol);
	if (fabs(entry1.pval - entry2.pval) > tol || tp_vec.size() % 2 != 0)
	  {
	    // Check if the midpoint between the curve endpoint instances is also equal to the
	    // isovalue
	    Point t1, t2;
	    Point midval;	
	    entry1.icurve.first->point(t1, entry1.at_start ? entry1.icurve.first->startparam() :
				       entry1.icurve.first->endparam());
	    entry2.icurve.first->point(t2, entry2.at_start ? entry2.icurve.first->startparam() :
				       entry2.icurve.first->endparam());
	    lrs.point(midval, 0.5*(t1[0]+t2[0]), 0.5*(t1[1]+t2[1])); 
	  
	    // if (fabs(midval[0] - isoval) > tol)
	    //   {
#ifdef DEBUG
	    std::cout << "Too large distance between pvals:" << fabs(entry1.pval - entry2.pval) << ", " << tol << std::endl;
	    std::cout << "midval: " << midval[0] << std::endl;
	    std::cout << "pval: ";
	    for (size_t ka=0; ka<tp_vec.size(); ++ka)
	      std::cout << tp_vec[ka].pval << ", ";
	    std::cout << std::endl;
	    std::ofstream ofpt("mismatch.g2");
	    ofpt << "400 1 0 4 255 0 0 255" << std::endl;
	    ofpt << tp_vec.size() << std::endl;
	    for (size_t ka=0; ka<tp_vec.size(); ++ka)
	      {
		double par = tp_vec[ka].at_start ? tp_vec[ka].icurve.first->startparam() :
		  tp_vec[ka].icurve.first->endparam();
		Point pt = tp_vec[ka].icurve.first->ParamCurve::point(par);
		ofpt << pt << " 0.0" << std::endl;
	      }
#endif
	    if (i+2 < tp_vec.size() && 
		fabs(tp_vec[i+2].pval - tp_vec[i+1].pval) < 
		fabs(entry2.pval - entry1.pval))
	      {
		++i;
		continue;
	      }
	    else if (fabs(midval[0] - isoval) > tol)
	      {
		i += 2;
		continue;
	      }
	  }

	if (entry1.icurve.first == entry2.icurve.first) {
	  // this is a single curve whose endpoints meet across this edge.  There are no more
	  // merges to be done.  The curve is finished, and can be returned.
	  // First check that it is really a closed curve and not a duplicate short curve
	  double len = entry1.icurve.first->estimatedCurveLength();
	  Point pt1, pt2;
	  entry1.icurve.first->point(pt1, entry1.icurve.first->startparam());
	  entry1.icurve.first->point(pt2, entry1.icurve.first->endparam());
	  if (len < tol && pt1.dist(pt2) > 0.5*len)
	    {
	      ++i;
	      continue;  // Short curve
	    }
	  else
	    finished_curves.push_back(entry1.icurve);
	} else {
	  IsectCurve new_curve;
	  try {
	    new_curve = join_isectcurves(entry1.icurve, entry2.icurve, entry1.at_start, entry2.at_start);
	  }
	  catch (...)
	    {
	      i += 2;
	      continue;
	    }
#ifdef DEBUG
	  new_curve.first->writeStandardHeader(of1);
	  new_curve.first->write(of1);
#endif

	  // replace references to the old curves with references to new_curve throughout
	  replace_segments(entry1.icurve, entry2.icurve, new_curve, mergemap);
	  replace_segments(entry1.icurve, entry2.icurve, new_curve, othermap);

	  for (size_t j = 0; j < tp_vec.size(); ++j) {
	  //for (size_t j = i+2; j < nmb; ++j) {
	    if ((tp_vec[j].icurve.first == entry1.icurve.first) |
		(tp_vec[j].icurve.first == entry2.icurve.first)) {
	      tp_vec[j].icurve = new_curve;
	      int ix = (d == XFIXED) ? 1 : 0;
	      double val1 = *(new_curve.first->coefs_begin()+ix);
	      double val2 = *(new_curve.first->coefs_end()-2+ix);
	      tp_vec[j].at_start = (fabs(val1-tp_vec[j].pval) < fabs(val2-tp_vec[j].pval));
	      //*(new_curve.first->coefs_begin() + (d==XFIXED ? 1 : 0)) == tp_vec[j].pval;
	    }
	  }

	  // updating boundary curve pointers if necessary
	  transform(bcurves.begin(), bcurves.end(), bcurves.begin(), [&](const IsectCurve& c) {
	      return ((c.first == entry1.icurve.first) | (c.first == entry2.icurve.first)) ? new_curve : c;});
#ifdef DEBUG3
	  std::ofstream of2("bcurves.g2");
	  for (size_t ka=0; ka<bcurves.size(); ++ka)
	    {
	      bcurves[ka].first->writeStandardHeader(of2);
	      bcurves[ka].first->write(of2);
	    }
	  int stop_break = 1;
#endif
	}
	i += 2;
      }
  }
}

// ----------------------------------------------------------------------------
array<double, 4> outer_boundary_pvals(const vector<LRSurfPtr>& patches)
// ----------------------------------------------------------------------------
{
  const auto start_u = [](LRSurfPtr l) {return l->startparam_u();};
  const auto end_u   = [](LRSurfPtr l) {return l->endparam_u();};
  const auto start_v = [](LRSurfPtr l) {return l->startparam_v();};
  const auto end_v   = [](LRSurfPtr l) {return l->endparam_v();};
  
  return array<double, 4>
    {(*min_element(patches.begin(), patches.end(),
	   [&start_u](LRSurfPtr l1, LRSurfPtr l2) {return start_u(l1) < start_u(l2);}))->startparam_u(),
     (*max_element(patches.begin(), patches.end(),
	   [&end_u](LRSurfPtr l1, LRSurfPtr l2) {return end_u(l1) < end_u(l2);}))->endparam_u(),
     (*min_element(patches.begin(), patches.end(),
	   [&start_v](LRSurfPtr l1, LRSurfPtr l2) {return start_v(l1) < start_v(l2);}))->startparam_v(),
     (*max_element(patches.begin(), patches.end(),
	   [&end_v](LRSurfPtr l1, LRSurfPtr l2) {return end_v(l1) < end_v(l2);}))->endparam_v()};
}

// ----------------------------------------------------------------------------
template<typename T>
void add_to_vec(vector<T>& target, const vector<T>& new_elements)
// ----------------------------------------------------------------------------
{
  target.insert(target.end(), new_elements.cbegin(), new_elements.cend());
}

// ----------------------------------------------------------------------------
template<typename T>
vector<T> expand_vec(const vector<vector<T>>& vv)
// ----------------------------------------------------------------------------
{
  vector<T> result;
  for (auto v : vv)
    add_to_vec(result, v);
  return result; 
}

// ----------------------------------------------------------------------------
template<typename K, typename T>
map<K, vector<T>> map_combine(const map<K, vector<T>>& m1,
			      const map<K, vector<T>>& m2)
// ----------------------------------------------------------------------------
{
  auto result = m1;
  result.insert(m2.begin(), m2.end());
  return result;
}
  
// ----------------------------------------------------------------------------
template<typename K, typename T>
vector<T> expand_map(const map<K, vector<T>>& m)
// ----------------------------------------------------------------------------
{
  vector<T> result;
  for (auto it : m)
    result.insert(result.end(), it.second.begin(), it.second.end());
  return result;
	    
}

// ----------------------------------------------------------------------------
template<typename C>
const C sort_container(C co)
// ----------------------------------------------------------------------------
{
  sort(co.begin(), co.end());
  return co;
}

// ----------------------------------------------------------------------------
template<typename T>
vector<T> remove_duplicates(vector<T>& c)
// ----------------------------------------------------------------------------
{
  auto tmp = sort_container(c);
  const auto it = unique(tmp.begin(), tmp.end());
  tmp.erase(it, tmp.end());
  return tmp;
}

// ----------------------------------------------------------------------------
  CurveVec trim_with_domain(CurveVec& curves, 
			    const CurveBoundedDomain* domain,
			    const double tol)
// ----------------------------------------------------------------------------
  {
    CurveVec trim_crvs;
    for (int ki=(int)curves.size()-1; ki>=0; --ki)
      {
	if (!curves[ki].first.get())
	  continue;
	vector<double> par_intervals;
	vector<pair<shared_ptr<SplineCurve>,shared_ptr<SplineCurve> > > sub_crvs;
	domain->findPcurveInsideSegments(*curves[ki].first, tol, par_intervals);
	for (size_t kj=0; kj<par_intervals.size(); kj+=2)
	  {
	    pair<shared_ptr<SplineCurve>,shared_ptr<SplineCurve> > sub;
	    shared_ptr<SplineCurve> sub_cv1(curves[ki].first->subCurve(par_intervals[kj],
								       par_intervals[kj+1]));
	    if (curves[ki].second.get())
	      {
		shared_ptr<SplineCurve> sub_cv2(curves[ki].second->subCurve(par_intervals[kj],
							     par_intervals[kj+1]));
		sub = make_pair(sub_cv1, sub_cv2);
	      }
	    else
	      {
		shared_ptr<SplineCurve> dummy;
		sub = make_pair(sub_cv1, dummy);
	      }

	    sub_crvs.push_back(sub);
	    // curves.push_back(sub);
	    // trim_crvs.push_back(sub);
	  }

	// Check if the sub curves belongs to a closed curve and should be
	// merged
	double mindist = std::numeric_limits<double>::max();
	int minix1 = -1, minix2 = -1;
	bool at_start1, at_start2;
	for (size_t kj=0; kj<sub_crvs.size(); ++kj)
	  {
	    Point p1 = 
	      sub_crvs[kj].first->ParamCurve::point(sub_crvs[kj].first->startparam());

	    Point p2 = 
	      sub_crvs[kj].first->ParamCurve::point(sub_crvs[kj].first->endparam());

	    for (size_t kr=kj+1; kr<sub_crvs.size(); ++kr)
	      {
		Point p3 = 
		  sub_crvs[kr].first->ParamCurve::point(sub_crvs[kr].first->startparam());

		Point p4 = 
		  sub_crvs[kr].first->ParamCurve::point(sub_crvs[kr].first->endparam());
		double d1 = p1.dist(p3);
		double d2 = p1.dist(p4);
		double d3 = p2.dist(p3);
		double d4 = p2.dist(p4);
		double dd = std::min(std::min(d1, d2), std::min(d3, d4));
		if (dd < mindist)
		  {
		    mindist = dd;
		    minix1 = (int)kj;
		    minix2 = (int)kr;
		    at_start1 = (std::min(d1, d2) < std::min(d3, d4));
		    at_start2 = (std::min(d1, d3) < std::min(d2, d4));
		  }
	      }
	  }
	if (mindist < tol)
	  {
	    // Merge curves
	    if (at_start1)
	      {
	    	sub_crvs[minix1].first->reverseParameterDirection();
	    	sub_crvs[minix1].second->reverseParameterDirection();
	      }
	    if (!at_start2)
	      {
	    	sub_crvs[minix2].first->reverseParameterDirection();
	    	sub_crvs[minix2].second->reverseParameterDirection();
	      }
	    sub_crvs[minix1].first->appendCurve(sub_crvs[minix2].first.get());
	    sub_crvs[minix1].second->appendCurve(sub_crvs[minix2].second.get());
	    sub_crvs.erase(sub_crvs.begin()+minix2);
	  }

	// Update contour curve information
	curves.erase(curves.begin()+ki);
	curves.insert(curves.end(), sub_crvs.begin(), sub_crvs.end());
	trim_crvs.insert(trim_crvs.end(), sub_crvs.begin(), sub_crvs.end());
      }
    return trim_crvs;
  }

// ----------------------------------------------------------------------------
CurveVec 
single_isocontour_merge(vector<CurveVec>& curves,
			const vector<pair<LRSurfPtr,LRSplineSurface::PatchStatus> > surf_patches,
			const LRSplineSurface& lrs,
			const double isoval,
			const double tol,
			const CurveBoundedDomain* domain)
// ----------------------------------------------------------------------------
{
#ifdef DEBUG
  std::ofstream of1("levelcurves.g2");
  for (size_t kb=0; kb<curves.size(); ++kb)
    {
      for (size_t ka=0; ka<curves[kb].size(); ++ka)
	{
	  if (curves[kb][ka].first.get())
	    {
	      curves[kb][ka].first->writeStandardHeader(of1);
	      curves[kb][ka].first->write(of1);
	    }
	}
    }
#endif
  // Trim candidate curves for intersection with trimming loop agains the loop
  // and remember the curves attached to the trimming loop
  CurveVec bcurves;
  if (domain)
    {
      double int_tol = std::min(tol, 1.0e-3);
      for (size_t ki=0; ki<curves.size(); ++ki)
	{
	  if (curves[ki].size() > 0 &&
	      surf_patches[ki].second == LRSplineSurface::INTERSECT)
	    {
	      // A trimming candidate
	      CurveVec trimcurves = trim_with_domain(curves[ki], domain, int_tol);
	      if (trimcurves.size() > 0)
		bcurves.insert(bcurves.end(), 
			       trimcurves.begin(), trimcurves.end());
	    }
	}
    }

  map<double, CurveVec> u_map, v_map; // map curves exiting patch domains along u and/or v
					// parameter direction
  for (int i = 0; i != (int)surf_patches.size(); ++i) {
    prepare_curvemaps(parameter_domain(surf_patches[i].first), tol, curves[i], u_map, v_map);
  }
#ifdef DEBUG
  vector<pair<double,CurveVec> > u_vec;
  vector<pair<double,CurveVec> > v_vec;
  for (auto iter=u_map.begin(); iter!=u_map.end(); ++iter)
    u_vec.push_back(std::make_pair((*iter).first,(*iter).second));
  for (auto iter=v_map.begin(); iter!=v_map.end(); ++iter)
    v_vec.push_back(std::make_pair((*iter).first,(*iter).second));
  int stop_break = 1;
#endif
  // identify intersection curves that do not need to be merged, and output them directly
  CurveVec result; // this will contain all isocontours (merged if necessary)
  const auto all_icurves = sort_container(expand_vec(curves));
  auto vecu = expand_map(u_map);  // VSK 0219. Do not combine different maps that may have the same key
  auto vecv = expand_map(v_map);
  add_to_vec(vecu, vecv);
  const auto mapped_icurves = sort_container(vecu);
  //const auto mapped_icurves = sort_container(expand_map(map_combine(u_map, v_map)));
  set_difference(all_icurves.begin(), all_icurves.end(),
		 mapped_icurves.begin(), mapped_icurves.end(), back_inserter(result));

  // remove mapped entries corresponding with outer boundaries (no merge will take place across
  // these)
  if (!domain)
    {
      RectDomain dom = lrs.containingDomain();
      const array<double, 4> outer_bnd = 
	{dom.umin(), dom.umax(), dom.vmin(), dom.vmax()};
      auto it = u_map.find(outer_bnd[0]); if (it != u_map.end()) {add_to_vec(bcurves, it->second); u_map.erase(it);}
      it      = u_map.find(outer_bnd[1]); if (it != u_map.end()) {add_to_vec(bcurves, it->second); u_map.erase(it);}
      it      = v_map.find(outer_bnd[2]); if (it != v_map.end()) {add_to_vec(bcurves, it->second); v_map.erase(it);}
      it      = v_map.find(outer_bnd[3]); if (it != v_map.end()) {add_to_vec(bcurves, it->second); v_map.erase(it);}
    }
  // clean up
  bcurves = remove_duplicates(bcurves);

#ifdef DEBUG
  std::ofstream of2("bcurves_init.g2");
  for (size_t ka=0; ka<bcurves.size(); ++ka)
    {
      bcurves[ka].first->writeStandardHeader(of2);
      bcurves[ka].first->write(of2);
    }
#endif

  // looping through each parameter line, merging curve segments, and spitting out finished
  // curves
  CurveVec merged_segs;
  merge_segments(u_map, v_map, XFIXED, lrs, isoval, tol, bcurves, merged_segs);
  merge_segments(v_map, u_map, YFIXED, lrs, isoval, tol, bcurves, merged_segs);

  add_to_vec(result, remove_duplicates(bcurves));
  add_to_vec(result, merged_segs);
  
  return remove_duplicates(result);
  //return result;
}

  
// ----------------------------------------------------------------------------
vector<CurveVec> 
merge_isocontours(vector<vector<CurveVec>>& curve_fragments,
		  const vector<pair<LRSurfPtr,LRSplineSurface::PatchStatus> > surf_frags,
		  const LRSplineSurface& lrs,
		  const std::vector<double>& isovals,
		  const double tol,
		  const CurveBoundedDomain* domain)
// ----------------------------------------------------------------------------
{
  // curve_fragments is indexed [surface patch][isovalue][set of curves]
  const int num_patches	    = (int)curve_fragments.size();  
  //assert(num_patches > 0);
  if (num_patches == 0)
    {
      std::cout << "No patches in merge_isocontours" << std::endl;
      THROW("No patches in merge_isocontours");
    }
  const int num_isocontours = (int)curve_fragments[0].size();
  
  vector<CurveVec> result; // should contain one entry per isovalue

  for (int i = 0; i != num_isocontours; ++i) {

    // collect all curve fragments that belong to the set of isocurves for a
    // particular isovalue.
    // auto isocurves_to_merge =
    //   apply_transform(curve_fragments, function<CurveVec(const vector<CurveVec>&)> {
    // 	  [i] (const vector<CurveVec>& v) 
    // 	    {
    // 	      return v[i];
    // 	    }});
    vector<CurveVec> isocurves_to_merge(curve_fragments.size());
    for (size_t kj=0; kj<curve_fragments.size(); ++kj)
      {
	if (curve_fragments[kj].size() <= i)
	  continue;
	if (curve_fragments[kj][i].size() > 0)
	  isocurves_to_merge[kj] = curve_fragments[kj][i];
      }

    // merge the fragments into complete curves, and store the resulting curves
    // in 'result'
    result.emplace_back(single_isocontour_merge(isocurves_to_merge, surf_frags, 
						lrs, isovals[i], tol, domain));
  }

  return result;
}


  
}; //end anonymous namespace
