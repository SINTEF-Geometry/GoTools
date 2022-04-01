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


#include "GoTools/lrsplines2D/TrimSurface.h"
#include "GoTools/lrsplines2D/TrimUtils.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/RectDomain.h"
#include "GoTools/geometry/LoopUtils.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/CurveBoundedDomain.h"
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/lrsplines2D/TrimCrvUtils.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include <iostream>
#include <fstream>

//#define DEBUG

using namespace Go;
using std::vector;
using std::pair;

//==============================================================================
bool TrimSurface::makeBoundedSurface(shared_ptr<ParamSurface>& surf,
				     bool isotrim[], 
				     vector<double>& points,
				     int tightness,
				     shared_ptr<BoundedSurface>& trim_surf,
				     bool only_outer)
//==============================================================================
{
  // Translate surface and points to origo
  RectDomain dom = surf->containingDomain();
  double umin = dom.umin();
  double umax = dom.umax();
  double vmin = dom.vmin();
  double vmax = dom.vmax(); 

  // Translate point cloud and remove points lying outside of the
  // surface domain
  Point vec(-0.5*(umin+umax), -0.5*(vmin+vmax), 0.0);
  int ptdim = 3;  // Point dimension
  int nmb_pts = (int)points.size()/ptdim;
  int ki, kj;
  for (ki=0, kj=0; kj<nmb_pts; )
    {
      if (points[ki]<umin || points[ki]>umax || 
  	  points[ki+1]<vmin || points[ki+1]>vmax)
	{
	  for (int kk=0; kk<ptdim; ++kk)
	    std::swap(points[ki+kk], points[ptdim*(nmb_pts-1)+kk]);
	  nmb_pts--;
	}
      else
	{
	  for (int kk=0; kk<ptdim; ++kk)
	    points[ki+kk] += vec[kk];
	  ++kj;
	  ki += ptdim;
	}
    }

  int min_nmb = 3;
  if (nmb_pts < min_nmb)
    return true;

#ifdef DEBUG
  std::ofstream of1("translated_cloud.g2");
  of1 << "400 1 0 0 " << std::endl;
  of1 << nmb_pts << std::endl;
  for (ki=0; ki<nmb_pts; ++ki)
    {
      for (kj=0; kj<ptdim; ++kj)
	of1 << points[ptdim*ki+kj] << " ";
      of1 << std::endl;
    }
#endif

  // Update parameter domain
  TrimCrvUtils::translateSurfaceDomain(surf.get(), vec);

  // if (surf->dimension() == 3)
  //   {
  //     surf->translate(vec);
  //   }
  // Set parameters for computations of trimming sequence
  int max_rec;
  int nmb_div;
  if (tightness <= 2)
    {
      max_rec = 1;
      nmb_div = (tightness == 2) ? 20 : 15;
    }
  else if (tightness <= 5)
    {
      max_rec = 2;
      nmb_div = (tightness == 3) ? 8 : ((tightness == 4) ? 12 : 15);
    }
  else
    {
      max_rec = 3;
      nmb_div = (tightness == 6) ? 10 : ((tightness == 7) ? 12 : 15);
    }


  // Compute trimming seqence
  RectDomain dom2 = surf->containingDomain();
  double domain[4];
  domain[0] = dom2.umin();
  domain[1] = dom2.umax();
  domain[2] = dom2.vmin();
  domain[3] = dom2.vmax();
  vector<vector<double> > seqs;
  //TrimUtils trimutil(points2, 1, domain);
  TrimUtils trimutil(&points[0], nmb_pts, 1, domain);
  trimutil.computeTrimSeqs(max_rec, nmb_div, seqs, only_outer);
  
  double udel, vdel;
  trimutil.getDomainLengths(udel, vdel);

#ifdef DEBUG
  std::cout << "Minimum sub domains, diag = "<< trimutil.getDomainDiag();
  std::cout << ", udel = " << udel << ", vdel = " << vdel << std::endl;

  std::ofstream of0("translated_seq.g2");
  for (size_t kr=0; kr<seqs.size(); ++kr)
    {
      if (seqs[kr].size() <= 2)
	continue;
      of0 << "410 1 0 0" << std::endl;
      of0 << seqs[kr].size()/2-2 << std::endl;
      size_t ki, kj;
      for (ki=0; ki<seqs[kr].size()-4; ki+=2)
	{
	  for (kj=0; kj<2; ++kj)
	    of0 << seqs[kr][ki+kj] << "  ";
	  of0 << 0 << " ";
	  for (; kj<4; ++kj)
	    of0 << seqs[kr][ki+kj] << "  ";
	  of0 << 0 << std::endl;
	}
    }
#endif 

  // @@@ VSK 0115. Use only the first sequence. If there are more, they
  // are likely to be incomplete.
  // if (only_outer && seqs.size() > 1)
  //   seqs.erase(seqs.begin()+1, seqs.end());
  int nmb_loops = only_outer ? 1 : (int)seqs.size();
  vector<vector<vector<double> > >   all_seqs(nmb_loops);
  for (int kh=0; kh<nmb_loops; ++kh)
    all_seqs[kh].push_back(seqs[kh]);

  double eps = std::max(udel, vdel);
  bool found = defineBdSurface(surf, domain, isotrim, eps, all_seqs, trim_surf);

  if (found)
    {
      // Translate back
      trim_surf->setParameterDomain(domain[0]-vec[0], domain[1]-vec[0],
				    domain[2]-vec[1], domain[3]-vec[1]);
      surf->setParameterDomain(umin, umax, vmin, vmax);
    }

  return found;
}


//==============================================================================
bool TrimSurface::defineBdSurface(shared_ptr<ParamSurface>& surf,
				  double domain[], bool isotrim[], double eps,
				  vector<vector<vector<double> > >& seqs,
				  shared_ptr<BoundedSurface>& trim_surf)
//==============================================================================
{
#ifdef DEBUG
  std::ofstream of01("trimming_cvs.g2");
#endif
 
  // Compute trimming loop
  // First extract parts of the trimming sequences following iso trim curves
  int nmb_match = 4;
  double tol = 1.0e-5;
  int nmb_loops = (int)seqs.size();
  vector<vector<shared_ptr<CurveOnSurface> > > loop(nmb_loops);
  int kh;
  for (kh=0; kh<nmb_loops; ++kh)
    {
      for (int kj=0; kj<4; ++kj)
      	{
      	  if (!isotrim[kj])
      	  	continue;

      	  int ix = (kj < 2) ? 0 : 1;
      	  for (int ki=0; ki<(int)seqs[kh].size();)
      	    {
      	      vector<vector<double> > split_seqs1 = 
      		TrimCrvUtils::extractConstParSeqs(seqs[kh][ki], ix, 
      						  domain[kj], nmb_match, 
      						  tol, eps);
      	      if (split_seqs1.size() > 1 || split_seqs1[0].size() == 4)
      		{
      		  seqs[kh].erase(seqs[kh].begin()+ki);
      		  seqs[kh].insert(seqs[kh].begin()+ki, 
      				      split_seqs1.begin(), split_seqs1.end());
      		  ki += (int)split_seqs1.size();
      		}
      	      else
      		++ki;
      	    }
      	}

      // Check if a bounded surface is required
      if (nmb_loops == 1)
	{
	  size_t ka;
	  for (ka=0; ka<seqs[kh].size(); ++ka)
	    if (seqs[kh][ka].size() > 4)
	      break;
	  if (ka == seqs[kh].size())
	    {
	      // No trimming required
	      return false;
	    }
	}

      // Split the remaining sequences from the outer loop in kinks
      double kink_tol = 5e-01; // 0.1 deg => 5.7 degrees.
      vector<vector<double> > split_seqs;
      for (size_t ki=0; ki<seqs[kh].size(); ++ki)
	{
	  vector<vector<double> > curr_seqs = 
	    TrimCrvUtils::splitCurvePointsInKinks(seqs[kh][ki], kink_tol);
	  split_seqs.insert(split_seqs.end(), curr_seqs.begin(), curr_seqs.end());
	}

      // Ensure a closed trimming loop
      TrimCrvUtils::makeConnectedLoop(split_seqs, tol);
  
      // Create trimming curves
      const int par_dim = 2;
      const int max_iter = 5;
      vector<shared_ptr<SplineCurve> > par_cvs;
      for (size_t ki = 0; ki < split_seqs.size(); ++ki)
	{
	  shared_ptr<SplineCurve> spline_cv_appr_2d
	    (TrimCrvUtils::approximateTrimPts(split_seqs[ki], par_dim, eps, 
					      max_iter));
	  par_cvs.push_back(spline_cv_appr_2d);
	}

#ifdef DEBUG
      for (size_t kr=0; kr<par_cvs.size(); ++kr)
	{
	  par_cvs[kr]->writeStandardHeader(of01);
	  par_cvs[kr]->write(of01);
	}
#endif

      // The curve should be CCW.
      // @@@ VSK 150310. Only outer curve loops
      // Assume one outer and the rest inner (this is not necessarily true)
      const double int_tol = 1e-06;
      vector<shared_ptr<ParamCurve> > par_cvs2(par_cvs.begin(), par_cvs.end());
      bool loop_is_ccw = LoopUtils::loopIsCCW(par_cvs2, eps, int_tol);
      if ((kh==0 && !loop_is_ccw) || (kh>0 && loop_is_ccw))
	{
	  //MESSAGE("We should change direction of the loop cv!");
	  for (size_t ki = 0; ki < par_cvs.size(); ++ki)
	    {
	      par_cvs[ki]->reverseParameterDirection();
	    }
	  reverse(par_cvs.begin(), par_cvs.end());
	}

      if (par_cvs.size() == 1 && par_cvs[0]->numCoefs() < 2*par_cvs[0]->order())
	{
	  // Refine curve to bring the coefficients closer to the curve
	  vector<double> knots, newknots;
	  par_cvs[0]->basis().knotsSimple(knots);
	  double lim = (knots[knots.size()-1] - knots[0])/(double)par_cvs[0]->order();
	  for (size_t kv=1; kv<knots.size(); ++kv)
	    {
	      if (knots[kv]-knots[kv-1] > lim)
		{
		  int num = (int)(knots[kv]-knots[kv-1]/lim);
		  num = std::max(num, 1);
		  double del = (knots[kv]-knots[kv-1])/(double)(num+1);
		  for (int ka=0; ka<num; ++ka)
		    newknots.push_back(knots[kv-1]+(ka+1)*del);
		}
	    }
	  par_cvs[0]->insertKnot(newknots);
	}
	  
      TrimCrvUtils::moveCurveCoefsInsideSurfDomain(surf.get(), par_cvs);

      vector<shared_ptr<ParamCurve> >  par_loop;
      par_loop.insert(par_loop.end(), par_cvs.begin(), par_cvs.end());
    
      // Create bounded surface
      for (size_t ki = 0; ki < par_loop.size(); ++ki)
	{
	  shared_ptr<CurveOnSurface> cv_on_sf(new CurveOnSurface(surf, par_loop[ki], true));
	  loop[kh].push_back(cv_on_sf);
	}
    }

  const bool fix_trim_cvs = false;
  const double epsgeo_bd_sf = 1e-03;
  trim_surf = 
    shared_ptr<BoundedSurface>(new BoundedSurface(surf, loop, epsgeo_bd_sf, 
						  fix_trim_cvs));
  // Parameter not set
  // int valid_state = 0;
  // bool is_valid = trim_surf->isValid(valid_state);
  // if (!is_valid)
  // {
  // 	MESSAGE("Created invalid BoundedSurface, valid_state = " << valid_state);
  // }

#ifdef DEBUG
  std::ofstream of("translated_trimmed.g2");
  trim_surf->writeStandardHeader(of);
  trim_surf->write(of);

  if (false)//trim_surf->dimension() == 1)
    {
      std::ofstream of2("translated_trimmed_3D.g2");
      shared_ptr<BoundedSurface> tmp_bd(trim_surf->clone());
      shared_ptr<ParamSurface> tmp_sf = tmp_bd->underlyingSurface();
      if (tmp_sf.get())
	{
	  shared_ptr<LRSplineSurface> tmp_lr = 
	    dynamic_pointer_cast<LRSplineSurface, ParamSurface>(tmp_sf);
	  tmp_lr->to3D();
	  tmp_bd->writeStandardHeader(of2);
	  tmp_bd->write(of2);
	}
    }
#endif

    return (trim_surf.get());
}

