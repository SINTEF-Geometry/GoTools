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

#include "GoTools/geometry/ExtremalPoint.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/utils/Point.h"
#include "GoTools/utils/Array.h"
#include "sislP.h"
#include <fstream>
#include <iostream>

#define DEBUG

using std::vector;
using std::pair;
using std::make_pair;
using namespace Go;

namespace {
void iterateExtremalPoints(shared_ptr<ParamSurface>& surface,
			   const CurveBoundedDomain* bddomain,
			   const Point& dir, double tol, double& maxval,
			   SISLSurf *surf1D, 
			   vector<pair<Point, double> >& ext_par);

  void iterateMax(SISLSurf *surf, double& val, double par[], double tol,
		  double start[], double end[], double par2[]);

}; // end anonymous namespace 

//namespace Go {
//===========================================================================
void 
ExtremalPoint::computeExtremalPoints(vector<shared_ptr<ParamSurface> >& surfs,
				     const Point& dir, double tol,
				     vector<pair<Point,Point> >& extremal_points)
//===========================================================================
{
  double maxval = std::numeric_limits<double>::lowest();
  for (size_t ki=0; ki<surfs.size(); ++ki)
    {
      extremalPoints(surfs[ki], dir, tol, maxval, extremal_points);
    }

  // Remove non-significant points
  for (size_t ki=0; ki<extremal_points.size(); )
    {
      double val = extremal_points[ki].first*dir;
      if (val < maxval - tol)
	extremal_points.erase(extremal_points.begin()+ki);
      else
	++ki;
    }
}

//===========================================================================
void 
ExtremalPoint::extremalPoints(shared_ptr<ParamSurface>& surface,
			      const Point& dir, double tol, double& maxval,
			      vector<pair<Point,Point> >& extremal_points)
//===========================================================================
{
  if (dir.dimension() != surface->dimension())
    THROW("Inconsistent dimension of geometric space");

  // Fetch associated spline surface
  shared_ptr<SplineSurface> tmp_spline;
  SplineSurface* surf = surface->getSplineSurface();
  if (!surf)
    {
      // Convert to spline surface
      tmp_spline = shared_ptr<SplineSurface>(surface->asSplineSurface());
      surf = tmp_spline.get();
    }
  if (!surf)
    return;   // No tools to compute extremal points

  const CurveBoundedDomain* bddomain = 0;
  shared_ptr<BoundedSurface> bsurf = 
    dynamic_pointer_cast<BoundedSurface, ParamSurface>(surface);
  if (bsurf.get())
    bddomain = &(bsurf->parameterDomain());

  // Create 1D surface
  shared_ptr<SplineSurface> surf2(surf->multCoefs(dir));

  // Create SISL surface and object, avoid copying of arrays
  SISLSurf* sislsf = GoSurf2SISL(*surf2, false);
  SISLObject *obj = newObject(SISLSURFACE);
  obj->s1 = sislsf;

  // Find maxima
  SISLIntdat *qintdat = 0;
  double currmax = maxval;
  int kstat = 0;
  s1161(obj, &currmax, tol, &qintdat, &kstat);
  if (kstat < 0)
    {
      if (obj) freeObject(obj);
      if (qintdat) freeIntdat(qintdat);
      THROW("SISL extremal point computation failed");
    }
  
  // Check result and transform to output format
  vector<pair<Point, double> > ext_par;
  if (qintdat)
    {
      for (int ki=0; ki<qintdat->ipoint; ++ki)
	{
	  // In this format, curves are represented as linked points.
	  // Thus, an interval extrema will end up as several independent points.
	  // Is this a problem?
	  Point par(qintdat->vpoint[ki]->epar[0], qintdat->vpoint[ki]->epar[1]);
	  if (bddomain)
	    {
	      Vector2D currpar(par[0], par[1]);
	      if (bddomain->isInDomain(currpar, tol))
		ext_par.push_back(make_pair(par, currmax));
	    }
	  else
	    ext_par.push_back(make_pair(par, currmax));
	}
    }
  
  if (obj) 
    {
      obj->s1 = 0;  // Keep surface for possible further use
      freeObject(obj);
    }

  if (qintdat) freeIntdat(qintdat);

  if (ext_par.size() > 0)
    {
      maxval = currmax;
    }
  else if (bddomain && kstat > 0)
    {
      // Maximum value found outside trimmed surface. Iterate for internal
      // maximum values
      iterateExtremalPoints(surface, bddomain, dir, tol, maxval, sislsf,
			    ext_par);
    }

  if (sislsf) freeSurf(sislsf);

  for (size_t kr=0; kr<ext_par.size(); ++kr)
    {
      Point pos = surface->point(ext_par[kr].first[0], ext_par[kr].first[1]);
      extremal_points.push_back(make_pair(pos, ext_par[kr].first));
    }
}

//} // end namespace Go

namespace {

  void iterateExtremalPoints(shared_ptr<ParamSurface>& surface,
			     const CurveBoundedDomain* bddomain,
			     const Point& dir, double tol, double& maxval,
			     SISLSurf *surf1D, vector<pair<Point,double> >& ext_par)
  {
    // Compute start points for the iteration in a regular grid. 
    // Compute density
      int n, m;
      double density = 1.0;
      int min_nmb = 4, max_nmb = 50;
      SurfaceTools::setResolutionFromDensity(surface, density, 
					     min_nmb, max_nmb, n, m);
    
      RectDomain dom = surface->containingDomain();
      
      // Indicates if an iterative search for a better maximum point is
      // worthwhile
      double fac = 0.2;

#ifdef DEBUG
  std::ofstream of1("cvs_ext.g2");
#endif
      // Fetch constant parameter curves in the 1. parameter direction
      int min_samples = 1;
      double u1 = dom.umin();
      double u2 = dom.umax();
      double udel = (u2 - u1)/(n+1);
      double par[2];
      par[0] = u1+0.5*udel;
      while (par[0] < u2)
	{
	  vector<shared_ptr<ParamCurve> > crvs = 
	    surface->constParamCurves(par[0], false);
	  if (crvs.size() == 0)
	    {
	      par[0] += udel;
	      continue;  // Outside trimmed surface
	    }

#ifdef DEBUG
	  for (size_t kr=0; kr<crvs.size(); ++kr)
	    {
	      shared_ptr<SplineCurve> tmpspl(crvs[kr]->geometryCurve());
	      if (tmpspl.get())
		{
		  tmpspl->writeStandardHeader(of1);
		  tmpspl->geometryCurve()->write(of1);
		}
	    }
#endif
	  // Distribute sampling points
	  double av_len = 0.0;
	  vector<double> cv_len(crvs.size());
	  double curr_len = 0.0;
	  for (size_t kr=0; kr<crvs.size(); ++kr)
	    {
	      double len = crvs[kr]->estimatedCurveLength();
	      av_len += len;
	      cv_len[kr] = len;
	      curr_len += (crvs[kr]->endparam()-crvs[kr]->startparam());
	    }
	  av_len /= (double)crvs.size();
	  
	  // Evaluate sampling points
	  int curr_nmb = (int)(m*(curr_len/(dom.vmax()-dom.vmin()))) + 1;
	  for (size_t kr=0; kr<crvs.size(); ++kr)
	    {
	      int nmb = (int)(curr_nmb*cv_len[kr]/av_len);
	      nmb = std::max(nmb, min_samples);
	      double v1 = crvs[kr]->startparam();
	      double v2 = crvs[kr]->endparam();
	      double vdel = (v2 - v1)/(double)(nmb+1);
	      par[1] = v1 + 0.5*vdel;

	      while (par[1] < v2)
		{
		  Point pos = crvs[kr]->point(par[1]);
#ifdef DEBUG
		  of1 << "400 1 0 4 0 255 0 255" << std::endl;
		  of1 << "1" << std::endl;
		  of1 << pos << std::endl;
#endif
		  double val = pos*dir;
		  if (val > maxval-fac*fabs(maxval))
		    {
		      // Iterative search
		      double start[2], end[2];
		      start[0] = std::max(u1, par[0] - 0.5*udel);
		      start[1] = std::max(v1, par[1] - 0.5*vdel);
		      end[0] = std::min(u2, par[0] + 0.5*udel);
		      end[1] = std::min(v2, par[1] + 0.5*vdel);
		      double par2[2];
		      iterateMax(surf1D, val, par, tol, start, end,
				 par2);
		      if (val >= maxval)
			{
			  if (bddomain)
			    {
			      // Check if the found point is inside the
			      Vector2D currpar(par2[0], par2[1]);
			      if (!bddomain->isInDomain(currpar, tol))
				{
				  // Check closest boundary point
				  Vector2D bdpar;
				  bddomain->closestOnBoundary(currpar,
							      bdpar, tol);
				  Point bd_pos = surface->point(bdpar[0], bdpar[1]);
				  val = bd_pos*dir;
				  if (val >= maxval)
				    {
				      par2[0] = bdpar[0];
				      par2[1] = bdpar[1];
				    }
				}
			    }
			}
		      if (val >= maxval)
			{
			  // Check if any previously found points must be
			  // removed
			  for (size_t kr=0; kr<ext_par.size(); )
			    {
			      if (ext_par[kr].second < val)
				ext_par.erase(ext_par.begin()+kr);
			      else
				++kr;
			    }

			  ext_par.push_back(make_pair(Point(par2[0],par2[1]), val));
			  maxval = val;
			}
		    }
		  par[1] += vdel;
		}
	    }
	  par[0] += udel;
	}
  }
				
  void iterateMax(SISLSurf *surf, double& val, double par[], double tol,
		  double start[], double end[], double par2[])
  {
    // Create a SISL point greater than the surface 
    double tmax = std::numeric_limits<double>::lowest();

    int nmb = surf->in1*surf->in2;
    for (int ki=0; ki<nmb; ++ki)
      tmax = std::max(tmax, surf->ecoef[ki]);

    SISLPoint *p1 = newPoint(&tmax, 1, 1);
    if (!p1)
      THROW("Error in scratch allocation");

    int kstat = 0;
    double seed[2];
    seed[0] = par[0];
    seed[1] = par[1];
    s1173(p1, surf, tol, start, end, seed, par2, &kstat);
    freePoint(p1);
    if (kstat < 0)
      THROW("Error in maximum point iteration");
    
    // Evaluate surface
    double val2[1];
    int left1 = 0, left2 = 0;
    s1424(surf, 0, 0, par2, &left1, &left2, val2, &kstat);
    if (kstat < 0)
      THROW("Error in evaluation of SISL surface");

    if (val2[0] >= val)
      val = val2[0];
    else
      {
	par2[0] = par[0];
	par2[1] = par[1];
      }
  }


}; // end anonymous namespace
