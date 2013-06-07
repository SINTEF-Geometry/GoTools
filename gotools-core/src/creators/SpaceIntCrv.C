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

#include "GoTools/creators/SpaceIntCrv.h"
#include "GoTools/creators/CoonsPatchGen.h"
//#include "GoTools/geometry/closestPtSurfSurfPlane.h"
#include <fstream>

using namespace Go;
using std::vector;
using std::max;
using std::min;

//===========================================================================
SpaceIntCrv::SpaceIntCrv(shared_ptr<ParamCurve> init_crv, int pardir,
			 vector<shared_ptr<CurveOnSurface> >& sfcv1,
			 vector<double> start1, vector<double> end1,
			 vector<shared_ptr<CurveOnSurface> >& sfcv2,
			 vector<double> start2, vector<double> end2,
			 vector<bool> opposite, bool same_orient)
  : init_crv_(init_crv), sfcv1_(sfcv1), sfcv2_(sfcv2),
    start1_(start1), end1_(end1), 
    start2_(start2), end2_(end2),
    opposite_(opposite), same_orient_(same_orient)
//===========================================================================
{
  segment_.reserve(sfcv1_.size()+1);
  if (same_orient_)
    {
      Point pt1 = sfcv1_[0]->parameterCurve()->point(sfcv1_[0]->startparam());
      segment_.push_back(pt1[pardir]);

      for (size_t ki=0; ki<start1_.size(); ++ki)
	{
	  Point pt2 = sfcv1_[ki]->parameterCurve()->point(sfcv1_[ki]->endparam());
	  segment_.push_back(pt2[pardir]);
	}
    }
  else
    {
      int nmb = (int)start1_.size();
      Point pt1 = sfcv1_[nmb-1]->parameterCurve()->point(sfcv1_[nmb-1]->endparam());
      segment_.push_back(pt1[pardir]);
      for (int ki=nmb-1; ki>=0; --ki)
	{
	  Point pt2 = sfcv1_[ki]->parameterCurve()->point(sfcv1_[ki]->startparam());
	  segment_.push_back(pt2[pardir]);
	}

    }
}

//===========================================================================
SpaceIntCrv::~SpaceIntCrv()
//===========================================================================
{
}
 
//===========================================================================
Point SpaceIntCrv::eval(double t) const
//===========================================================================
{
  Point result;
  evaluate(t, 0, &result);
  return result;
}

//===========================================================================
void SpaceIntCrv::eval(double t, int n, Point der[]) const
//===========================================================================
{
  if (n > 1)
    n = 1;
  if (n < 0)
    n = 0;
  evaluate(t, n, der);

}

//===========================================================================
double SpaceIntCrv::start() const
//===========================================================================
{
  return segment_[0];
}

//===========================================================================
double SpaceIntCrv::end() const
//===========================================================================
{
  return segment_[segment_.size()-1];;
}

//===========================================================================
int SpaceIntCrv::dim() const
//===========================================================================
{
// One geometry curve
  return sfcv1_[0]->dimension();  
}

//===========================================================================
bool SpaceIntCrv::approximationOK(double par, 
				  Point approxpos,
				  double tol1, double tol2) const
//===========================================================================
{
  // Check geometry curve
  Point pos = eval(par);
  if (approxpos.dist(pos) > tol1)
    return false;

  return true;  // The accuracy is good enough
}


//===========================================================================
void SpaceIntCrv::evaluate(double t, int n, Point result[]) const
//===========================================================================
{
  // double ang_tol = 0.01;
  double eps = 1.0e-6;

  // Evaluation is performed by iterating to the intersection point between
  // the two surfaces defined by the surface curve and a plane being
  // orhtogonal to the intersection curve in the given parameter

  // NB! At most 1. derivitive of curves are computed
  if (n > 1)
    n = 1;

  // Find the correct surface curve
  size_t idx;
  double tstart = start();
  double tend = end();
  double fac = (tend - tstart)/(init_crv_->endparam() - init_crv_->startparam());
  double t1 = tstart + (t - init_crv_->startparam())*fac;

  for (idx=0; idx<sfcv1_.size(); ++idx)
    if (segment_[idx] <= t1 && segment_[idx+1] > t1)
      break;
  if (idx == sfcv1_.size())
    {
      if (t1 >= segment_[0]-eps && t1 <= segment_[1]+eps)
	idx = 0;
      else
	idx = sfcv1_.size() - 1;
    }

  // int idx2 = idx;
  if (!same_orient_) 
    idx = sfcv1_.size()-idx-1;
  if (!same_orient_)
    t1 = tend - t1 + tstart;
  double t2;
  //if (same_orient_)
    t2 = start1_[idx] + 
      (t1-segment_[idx])*(end1_[idx]-start1_[idx])/(segment_[idx+1]-segment_[idx]);
//   else
//     t2 = end1_[idx] - 
//       (t1-segment_[idx2])*(end1_[idx]-start1_[idx])/(segment_[idx2+1]-segment_[idx2]);
    //t2 = end1_[idx] - (segment_[idx2+1] - t1 + start1_[idx]);

  // Evaluate initial curve
  vector<Point> der1(n+1);
  init_crv_->point(der1, t, n);

  // Testing
  Point tmp_pt;
  sfcv1_[idx]->point(tmp_pt, t2);

  if (sfcv2_[idx].get())
    {
      // Find the corresponding point in the other surface curve
      double rel = (end2_[idx] - start2_[idx])/(end1_[idx] - start1_[idx]);
      double guess = (!opposite_[idx]) ? start2_[idx] + (t2-start1_[idx])*rel : 
	end2_[idx] - (t2-start1_[idx])*rel;
      guess = std::max(start2_[idx], std::min(end2_[idx], guess)); // Numerics
  
      double par, clo_dist;
      Point clo_pt;
      sfcv2_[idx]->closestPoint(der1[0], start2_[idx], end2_[idx], 
				par, clo_pt, clo_dist, &guess);
//       std::cout << "Opposite: " << opposite_[idx] << "< orient: " << same_orient_;
//       std::cout << ", guess: " << guess << ", par: " << par << std::endl;

      // Compute the closest point to the point on the first curve onto the surface
      // corresponding to the second curve
      double clo_u, clo_v, clo_dist2;
      Point clo_pt2;
      Point sf_par2 = sfcv2_[idx]->faceParameter(par);
      sfcv2_[idx]->underlyingSurface()->closestPoint(clo_pt, clo_u, clo_v, clo_pt2,
						     clo_dist2, eps, NULL, sf_par2.begin());

      int write_srf = false;
      if (write_srf)
	{
	  std::ofstream out("SpaceIntCrv_sf.g2");
	  sfcv2_[idx]->underlyingSurface()->writeStandardHeader(out);
	  sfcv2_[idx]->underlyingSurface()->write(out);
	}

      result[0] = clo_pt2;
      if (n == 1)
	{
	  // Evaluate second surface with 1. derivatives
	  vector<Point> der2 = 
	    sfcv2_[idx]->underlyingSurface()->point(clo_u, clo_v, 1);
	
	  // Project average of input tangent onto the second surface
	  double u, v;
	  double len = der1[1].length();
	  der1[1].normalize();
	  CoonsPatchGen::blendcoef(&der2[1][0], &der2[2][0],
				   &der1[1][0], 3, 1, &u, &v);
	  Point proj1 = u*der2[1] + v*der2[2];
	    
	  proj1.normalize();
	  result[1] = proj1;
	  result[1] *= len;

	}
    }
  else
    {
      // Use information related to the first surface curve
      result[0] = der1[0];
      if (n == 1)
	{
	  result[1] = der1[1];
	}
    }
}

