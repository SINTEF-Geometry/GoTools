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

#include "GoTools/lrsplines2D/LRSurfApprox.h"
#include "GoTools/lrsplines2D/LRSurfSmoothLS.h"
#include "GoTools/lrsplines2D/LRBSpline2DUtils.h"
#include "GoTools/lrsplines2D/LinDepUtils.h"
#include "GoTools/creators/ApproxSurf.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/lrsplines2D/LRSplinePlotUtils.h"
#include <iostream>
#include <iomanip>
#include <fstream>

#define DEBUG

using std::vector;

using namespace Go;

//==============================================================================
LRSurfApprox::LRSurfApprox(vector<double>& points, 
			   int dim, double epsge, bool closest_dist,
			   bool repar)
  : points_(points), maxdist_(-10000.0), avdist_(0.0), outsideeps_(0), aepsge_(epsge),
    smoothweight_(1.0e-3), repar_(repar), check_close_(closest_dist), to3D_(-1)
//==============================================================================
{
  edge_derivs_[0] = edge_derivs_[1] = edge_derivs_[2] = edge_derivs_[3] = 0;

  increase_domain_ = false; //true;
  increase_fac_ = 0.05; //0.2;
  fix_boundary_ = false; //true;

  // Create an LR B-spline surface with the domain given by the 
  // parameter domain of the points. Only one element will be
  // created
  makeInitSurf(dim);
}

//==============================================================================
LRSurfApprox::LRSurfApprox(shared_ptr<SplineSurface>& srf,
			   vector<double>& points, 
			   double epsge, bool closest_dist,
			   bool repar)
  : points_(points), maxdist_(-10000.0), avdist_(0.0), outsideeps_(0), aepsge_(epsge),
    smoothweight_(1.0e-3), repar_(repar), check_close_(closest_dist), to3D_(-1)
//==============================================================================
{
  edge_derivs_[0] = edge_derivs_[1] = edge_derivs_[2] = edge_derivs_[3] = 0;

  // Create an LR B-spline surface based on the given spline surface
  makeInitSurf(srf);
}

//==============================================================================
LRSurfApprox::LRSurfApprox(shared_ptr<LRSplineSurface>& srf,
			   vector<double>& points, 
			   double epsge, bool closest_dist,
			   bool repar)
//==============================================================================
  : points_(points), maxdist_(-10000.0), avdist_(0.0), outsideeps_(0), aepsge_(epsge),
    smoothweight_(1.0e-3), repar_(repar), check_close_(closest_dist), to3D_(-1)
{
  edge_derivs_[0] = edge_derivs_[1] = edge_derivs_[2] = edge_derivs_[3] = 0;
  srf_ = srf;
}

//==============================================================================
LRSurfApprox::LRSurfApprox(int ncoef_u, int order_u, int ncoef_v, int order_v,
			   vector<double>& points, 
			   int dim, double epsge, bool closest_dist,
			   bool repar)
//==============================================================================
  : points_(points), maxdist_(-10000.0), avdist_(0.0), outsideeps_(0), aepsge_(epsge),
    smoothweight_(1.0e-3), repar_(repar), check_close_(closest_dist), to3D_(-1)
{
  edge_derivs_[0] = edge_derivs_[1] = edge_derivs_[2] = edge_derivs_[3] = 0;

  increase_domain_ = false; //true;
  increase_fac_ = 0.05; //0.2;
  fix_boundary_ = false; //true;

  // Create an LR B-spline surface with unset coefficients and the domain
  // given by the parameter domain of the points. The size of the spline
  // space is given. The knots will be equally spaced
  makeInitSurf(dim, ncoef_u, order_u, ncoef_v, order_v);
}

//==============================================================================
LRSurfApprox::~LRSurfApprox()
//==============================================================================
{
}

//==============================================================================
 shared_ptr<LRSplineSurface> LRSurfApprox::getApproxSurf(double& maxdist, 
							 double& avdist,
							 int& nmb_out_eps, 
							 int max_iter)
//==============================================================================
{
// #ifdef DEBUG
//   std::ofstream of0("init0_sf.g2");
//   shared_ptr<LRSplineSurface> tmp0(srf_->clone());
//   tmp0->to3D();
//   tmp0->writeStandardHeader(of0);
//   tmp0->write(of0);
//   of0 << std::endl;
//   LineCloud lines0 = tmp0->getElementBds();
//   lines0.writeStandardHeader(of0);
//   lines0.write(of0);
// #endif

  // Initiate approximation engine
  if (fix_boundary_)
    setFixBoundary(true);

  LRSurfSmoothLS LSapprox(srf_, coef_known_);

  // Initiate with data points
  LSapprox.addDataPoints(points_);

  // Initial smoothing of LR B-spline surface
  performSmooth(&LSapprox);

#ifdef DEBUG
  std::ofstream of1("init_sf.g2");
  shared_ptr<LRSplineSurface> tmp;
  if (srf_->dimension() == 1)
    {
      tmp = shared_ptr<LRSplineSurface>(srf_->clone());
      tmp->to3D();
    }
  else
    tmp = srf_;
  tmp->writeStandardHeader(of1);
  tmp->write(of1);
  of1 << std::endl;
  LineCloud lines = tmp->getElementBds();
  lines.writeStandardHeader(of1);
  lines.write(of1);
#endif

  // Compute accuracy in data points
  computeAccuracy();

  for (int ki=0; ki<max_iter; ++ki)
    {
      // Check if the requested accuracy is reached
      if (maxdist_ <= aepsge_ || outsideeps_ == 0)
	break;

      // Refine surface
      refineSurf();
#ifdef DEBUG
      std::ofstream of2("refined_sf.g2");
      shared_ptr<LRSplineSurface> tmp2;
      if (srf_->dimension() == 1)
	{
	  tmp2 = shared_ptr<LRSplineSurface>(srf_->clone());
	  tmp2->to3D();
	}
      else
	tmp2 = srf_;
      tmp2->writeStandardHeader(of2);
      tmp2->write(of2);
      of2 << std::endl;
      LineCloud lines2 = tmp2->getElementBds();
      lines2.writeStandardHeader(of2);
      lines2.write(of2);

      std::ofstream of3("point_clouds.g2");
      int del = srf_->dimension() + 2;
      for (LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
	   it != srf_->elementsEnd(); ++it)
	{
	  vector<double>& elem_data = it->second->getDataPoints();
	  int nmb = (int)elem_data.size()/del;
	  if (elem_data.size() > 0)
	    {
	      vector<double> tmppt;
	      if (srf_->dimension() == 1)
		tmppt = elem_data;
	      else
		{
		  tmppt.reserve(3*elem_data.size()/del);
		  for (int kr=0; kr<nmb; ++kr)
		    tmppt.insert(tmppt.end(), elem_data.begin()+kr*del+2, 
				 elem_data.begin()+(kr+1)*del);
		  
		}
	      PointCloud3D cloud(tmppt.begin(), nmb);
	      cloud.writeStandardHeader(of3);
	      cloud.write(of3);
	    }
	}
#endif

      // Update coef_known from information in LR B-splines
      //updateCoefKnown();
      unsetCoefKnown();
      if (fix_boundary_)
	setFixBoundary(true);
  
      // Check for linear independence (overloading)
      vector<LRBSpline2D*> funs = LinDepUtils::unpeelableBasisFunctions(*srf_);
      std::cout << "Number of unpeelable functions: " << funs.size() << std::endl;
     
      // Update surface
      LSapprox.updateLocals();
      performSmooth(&LSapprox);
  
#ifdef DEBUG
      std::ofstream of4("updated_sf.g2");
      shared_ptr<LRSplineSurface> tmp3;
      if (srf_->dimension() == 1)
	{
	  tmp3 = shared_ptr<LRSplineSurface>(srf_->clone());
	  tmp3->to3D();
	}
      else
	tmp3 = srf_;
      tmp3->writeStandardHeader(of4);
      tmp3->write(of4);
      of4 << std::endl;
      LineCloud lines3 = tmp3->getElementBds();
      lines3.writeStandardHeader(of4);
      lines3.write(of4);
#endif
  
      if (ki == to3D_ && srf_->dimension() == 1)
	{
	  // Turn the current function into a 3D surface
	  // before continuing the iteration
	  turnTo3D();
	}

      computeAccuracy();
    }

  // Set accuracy information
  maxdist = maxdist_;
  avdist = avdist_;
  nmb_out_eps = outsideeps_;
  
  return srf_;
}

//==============================================================================
void LRSurfApprox::performSmooth(LRSurfSmoothLS *LSapprox)
//==============================================================================
{
  double wgt1 = 0.1*smoothweight_;
  double wgt3 = 0.5*smoothweight_;
  double wgt2 = (1.0 - wgt3)*smoothweight_;
  LSapprox->setOptimize(wgt1, wgt2, wgt3);
  
  double approx_weight = 1.0-wgt1-wgt2-wgt3;  
  LSapprox->setLeastSquares(approx_weight);

  shared_ptr<LRSplineSurface> lrsf_out;
  int isOK = LSapprox->equationSolve(lrsf_out);
#ifdef DEBUG
  std::cout << "isOK: " << isOK << std::endl;
#endif

  srf_ = lrsf_out;
}

//==============================================================================
void LRSurfApprox::computeAccuracy()
//==============================================================================
{
  // Check the accuracy of all data points, element by element
  // Note that only points more distant from the surface than the tolerance
  // are considered in avdist_

  // Initiate accuracy information
  maxdist_ = 0.0;
  avdist_ = 0.0;
  outsideeps_ = 0;

  bool remove_distant = false;
  double threshhold = 0;
#ifdef DEBUG
  if (false)
    {
  std::cout << "Remove distant points? (0/1) " << std::endl;
  std::cin >> remove_distant;
  std::cout << "Give threshhold: " << std::endl;
  std::cin >> threshhold;
    }
  std::ofstream of("error_pnts.g2");
#endif

  RectDomain rd = srf_->containingDomain();
  int dim = srf_->dimension();
  int del = 2 + dim;
  LRSplineSurface::ElementMap::const_iterator it;
  int num = srf_->numElements();
  int maxiter = 4;
  int kj;
  for (it=srf_->elementsBegin(), kj=0; it != srf_->elementsEnd(); ++it, ++kj)
    {
      if (!it->second->hasDataPoints())
	{
	  // Reset accuracy information in element
	  it->second->resetAccuracyInfo();
	  continue;   // No points in which to check accuracy
	}
      
      double umin = it->second->umin();
      double umax = it->second->umax();
      double vmin = it->second->vmin();
      double vmax = it->second->vmax();
      vector<Point> error_pts;
      vector<double>& points = it->second->getDataPoints();
      int nmb = it->second->nmbDataPoints();

      // Local error information
      double max_err = 0.0;
      double av_err = 0.0;
      int outside = 0;

      int ki;
      double *curr;
      for (ki=0, curr=&points[0]; ki<nmb; )
	{
	  double dist;
	  Point curr_pt(curr+(dim==3)*2, curr+del);
	  if (check_close_ && dim == 3)
	    {
	      // Compute closest point
	      // VSK. 052013. Should be changed to make use of the fact that we know
	      // the element (at least initially)
	      double upar, vpar;
	      Point close_pt;
	      srf_->closestPoint(curr_pt, upar, vpar, close_pt,
				 dist, aepsge_, maxiter, &rd, curr);
	      if (to3D_ >= 0)
		dist = fabs(curr[del-1]-close_pt[2]);
	      if (repar_)
		{
		  curr[0] = upar;
		  curr[1] = vpar;

		  if (upar < umin || upar > umax || vpar < vmin || vpar > vmax)
		    {
		      //std::cout << "Closest parameter out of element " << std::endl;

		      // Find element
		      Element2D *elem = srf_->coveringElement(upar, vpar);
		      elem->addDataPoints(points.begin()+ki*del, 
					  points.begin()+(ki+1)*del);
		      it->second->eraseDataPoints(points.begin()+ki*del, 
						  points.begin()+(ki+1)*del);
		      nmb--;
		    }
		  else
		    {
		      curr += del;
		      ki++;
		    }
		}
	      else
		{
		  curr += del;
		  ki++;
		}
	    }
	  else
	    {
	      // Evaluate
	      Point pos;
	      srf_->point(pos, curr[0], curr[1]);
	      dist = pos.dist(Point(curr+2, curr+del));
	      curr += del;
	      ki++;
	    }

	  if (remove_distant && dist > threshhold)
	    {
	      points.erase(points.begin()+ki*del, points.begin()+(ki+1)*del);
	      nmb--;
	      ki--;
	      curr -= del;
	    }
	  else
	    {
	      // Accumulate approximation error
	      maxdist_ = std::max(maxdist_, dist);
	      max_err = std::max(max_err, dist);
	      if (dist > aepsge_)
		{
		  avdist_ += dist;
		  outsideeps_++;
		  av_err += dist;
		  outside++;
		  
		  // Accumulate error points
		  error_pts.push_back(curr_pt);
		}
	    }
	}
      if (outside > 0)
	av_err /= (double)outside;

      // Store accuracy information in the element
      it->second->setAccuracyInfo(av_err, max_err, outside);

      if (remove_distant)
	{
	  // Replace point information
	  it->second->eraseDataPoints();
	  it->second->addDataPoints(points.begin(), points.end());
	}
#ifdef DEBUG
      if (error_pts.size() > 0)
	{
	  of << "400 1 0 4 255 0 0 255" << std::endl;
	  of << error_pts.size() << std::endl; 
	  for (size_t kh=0; kh<error_pts.size(); ++kh)
	    of << error_pts[kh] << std::endl;
	}
#endif
    }
  if (outsideeps_ > 0)
    avdist_ /= (double)outsideeps_;
}

//==============================================================================
void LRSurfApprox::refineSurf()
//==============================================================================
{
#ifdef DEBUG
  int idx=0;
  for (LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
       it != srf_->elementsEnd(); ++it)
    {
      if (it->second->hasAccuracyInfo())
  	{
  	  double av_err, max_err;
  	  int nmb_out;
  	  it->second->getAccuracyInfo(av_err, max_err, nmb_out);
  	  std::cout << "Element " << idx << ", domain: [(" << it->second->umin() << ",";
  	  std::cout << it->second->umax() << ")x(" << it->second->vmin() << ",";
  	  std::cout << it->second->vmax() << ")]" << std::endl;
  	  std::cout << "Average error: " << av_err << std::endl;
  	  std::cout << "Maximum error: " << max_err << std::endl;
  	  std::cout << "Number of points: " << it->second->nmbDataPoints();
  	  std::cout << ", number of error points: " << nmb_out << std::endl;
  	}
      ++idx;
    }
#endif

  // Construct indexed element array
  int num_el = srf_->numElements();
  vector<Element2D*> elem(num_el);
  size_t kr = 0;
  for (LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
       it != srf_->elementsEnd(); ++it)
    elem[kr++] = it->second.get();
  
  // Sort elements according to average error
  vector<int> el_perm(num_el);
  int ki, kj;
  for (ki=0; ki<num_el; ++ki)
    el_perm[ki] = ki;

  // Do the sorting
  for (ki=0; ki<num_el; ++ki)
    for (kj=ki+1; kj<num_el; ++kj)
      {
	double av_err1 = elem[el_perm[ki]]->getAverageError();
 	double av_err2 = elem[el_perm[kj]]->getAverageError();

	if (av_err1 < av_err2)
	  std::swap(el_perm[ki], el_perm[kj]);
      }
  
  // Define threshhold for refinement. The average error includes only those points that
  // are outside of the resolution
  double max_av = 0;
  double mean_av = 0.0;
  int num_err = 0;
  for (ki=0; ki<num_el; ++ki)
    {
      double av_err = elem[el_perm[ki]]->getAverageError();
      max_av = std::max(max_av, av_err);
      mean_av += av_err;
      if (elem[el_perm[ki]]->getMaxError() > aepsge_)
	num_err++;
    }
  mean_av /= (double)num_el;

  double threshhold = 0.25*mean_av;  // Need to experiment with this

  vector<LRSplineSurface::Refinement2D> refs;
  int nmb_ref = std::max(std::min((int)el_perm.size(), 4), (int)(0.75*num_err));
  
  for (kr=0; kr<el_perm.size(); )
    {
      size_t nmb_perm = el_perm.size();
      if (elem[el_perm[kr]]->getAverageError() < threshhold && refs.size() > nmb_ref)
	break;  // No more refinements at the current stage
      if (elem[el_perm[kr]]->getMaxError() < aepsge_)
	break;

      // Check feasability of split
      size_t nmb_refs = refs.size();
      vector<Element2D*> elements;  // Elements affected by the refinement(s)
      checkFeasibleRef(elem[el_perm[kr]], refs, elements);
      if (elements.size() > 0)
	{
	  // Remove affected elements from pool
	  for (size_t kh=0; kh<elements.size(); ++kh)
	    {
	      for (ki=0; ki<(int)el_perm.size();)
		{
		  if (elements[kh] == elem[el_perm[ki]])
		    {
		      el_perm.erase(el_perm.begin()+ki);
		      break;
		    }
		  else
		    ki++;
		}
	    }
	  if (el_perm.size() == nmb_perm)
	    kr++;
	}
      else
	kr++;
    }
  
#ifdef DEBUG
  std::ofstream of("refine0.dat");
  std::streamsize prev = of.precision(15);
  for (kr=0; kr<refs.size(); ++kr)
    {
      of << refs[kr].kval << "  " << refs[kr].start << "  " << refs[kr].end;
      of << "  " << refs[kr].d << "  " << refs[kr].multiplicity << std::endl;
    }
#endif

  for (kr=0; kr<refs.size(); ++kr)
    {
#ifdef DEBUG
      std::cout << "Refine nr " << kr << ": " << refs[kr].kval << "  " << refs[kr].start << "  ";
      std::cout << refs[kr].end << "  " << refs[kr].d << "  " << refs[kr].multiplicity << std::endl;
#endif
      // Perform refinements, one at the time to keep information stored in the elements
      srf_->refine(refs[kr], true /*false*/);
    }

  // // Update coef_known from information in LR B-splines
  // //updateCoefKnown();
  // unsetCoefKnown();
  //  if (fix_boundary_)
  //   setFixBoundary(true);

  
  // // Prelimenary implementation. Write accuracy information to standard output
  // // and get refinement informatiation
  // // For each element, check if any data points are stored
  // int idx=0;
  // for (LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
  //      it != srf_->elementsEnd(); ++it)
  //   {
  //     if (it->second->hasAccuracyInfo())
  // 	{
  // 	  double av_err, max_err;
  // 	  int nmb_out;
  // 	  it->second->getAccuracyInfo(av_err, max_err, nmb_out);
  // 	  std::cout << "Element " << idx << ", domain: [(" << it->second->umin() << ",";
  // 	  std::cout << it->second->umax() << ")x(" << it->second->vmin() << ",";
  // 	  std::cout << it->second->vmax() << ")]" << std::endl;
  // 	  std::cout << "Average error: " << av_err << std::endl;
  // 	  std::cout << "Maximume error: " << max_err << std::endl;
  // 	  std::cout << "Number of points: " << it->second->nmbDataPoints();
  // 	  std::cout << ", number of error points: " << nmb_out << std::endl;
  // 	}
  //     ++idx;
  //   }

  // // Specify refinement
  // while (true)
  //   {
  //     int more;
  //     std::cout << "Refine more? " << std::endl;
  //     std::cin >> more;
  //     if (!more)
  // 	break;
  //     double parval, start, end;
  //     int dir;
  //     int mult = 1;
  //     std::cout << "New knot, give pardir (0,1), fixed value, start and end: ";
  //     std::cout << std::endl;
  //     std::cin >> dir;
  //     std::cin >> parval;
  //     std::cin >> start;
  //     std::cin >> end;

  //     // The refinement must be performed inserting one knot segment at the
  //     // time in order to maintain information stored in the element
  //     srf_->refine((dir==0) ? YFIXED : XFIXED, parval, start, end, mult);

  //     // Update coef_known from information in LR B-splines
  //     updateCoefKnown();
  //   }
  
}

//==============================================================================
void LRSurfApprox::makeInitSurf(int dim)
//==============================================================================
{
  // Create a cubic Bezier surface represented as an LR B-spline surface
  // Compute domain
  double umin, umax, vmin, vmax;
  computeParDomain(dim, umin, umax, vmin, vmax);

  vector<double> knots_u(8);
  vector<double> knots_v(8);

  if (increase_domain_)
    {
      umin -= increase_fac_*(umax-umin);
      umax += increase_fac_*(umax-umin);
      vmin -= increase_fac_*(vmax-vmin);
      vmax += increase_fac_*(vmax-vmin);
    }

  int kj;
  for (kj=0; kj<4; ++kj)
    {
      knots_u[kj] = umin;
      knots_v[kj] = vmin;
      knots_u[kj+4] = umax;
      knots_v[kj+4] = vmax;
    }

  makeInitSurf(dim, 4, 4, 4, 4, &knots_u[0], &knots_v[0]);
}

//==============================================================================
void LRSurfApprox::makeInitSurf(int dim, int ncoef_u, int order_u, int ncoef_v, 
				int order_v)
//==============================================================================
{
  // Compute domain
  double umin, umax, vmin, vmax;
  computeParDomain(dim, umin, umax, vmin, vmax);

  vector<double> knots_u(ncoef_u+order_u);
  vector<double> knots_v(ncoef_v+order_v);

  if (increase_domain_)
    {
      umin -= increase_fac_*(umax-umin);
      umax += increase_fac_*(umax-umin);
      vmin -= increase_fac_*(vmax-vmin);
      vmax += increase_fac_*(vmax-vmin);
    }

  int kj;
  double del_u = (umax - umin)/(double)(ncoef_u-order_u+1);
  double del_v = (vmax - vmin)/(double)(ncoef_v-order_v+1);
  for (kj=0; kj<order_u; ++kj)
    knots_u[kj] = umin;
   for (; kj<ncoef_u; ++kj)
    knots_u[kj] = umin + (kj-order_u+1)*del_u;
  for (; kj<ncoef_u+order_u; ++kj)
    knots_u[kj] = umax;
  for (kj=0; kj<order_v; ++kj)
    knots_v[kj] = vmin;
  for (; kj<ncoef_v; ++kj)
    knots_v[kj] = vmin + (kj-order_v+1)*del_v;
  for (; kj<ncoef_v+order_v; ++kj)
    knots_v[kj] = vmax;

  makeInitSurf(dim, ncoef_u, order_u, ncoef_v, order_v, &knots_u[0], &knots_v[0]);
}

//==============================================================================
void LRSurfApprox::makeInitSurf(int dim, int ncoef_u, int order_u, int ncoef_v, 
				int order_v, double *knots_u, double *knots_v)
//==============================================================================
{
  vector<double> coefs(ncoef_u*ncoef_v*dim, 0.0);
  shared_ptr<SplineSurface> sf1(new SplineSurface(ncoef_u, ncoef_v, 
						  order_u, order_v,  
						  knots_u,     // ptr to knots in u
						  knots_v,     // ptr to knots in v
						  &coefs[0], // ptr to coefs
						  dim));

  // Reformatting data to the format that ApproxSurf wants (separating parameter values and 
  // data points in two distinct vectors).
  int del = 2 + dim;  // Number of entities for each point
  int nmb_pts = (int)points_.size()/del;
  vector<double> pts, param;
  pts.reserve(dim*nmb_pts);
  param.reserve(2 * nmb_pts);
  int ki;
  for (vector<double>::iterator it = points_.begin(); it != points_.end(); ) 
    {
      for (ki=0; ki<2; ++ki)
	param.push_back(*it++);
      for (ki=0; ki<dim; ++ki)
	pts.push_back(*it++);
    }
    
  // Approximate data points to create initial surface
  ApproxSurf asurf(sf1,
		   pts, //  scattered data points
		   param, // parameter values of scattered data points
		   dim, // dimension
		   aepsge_, // geometric tolerance
		   0, // 'constdir' - doesn't matter as we are not going to reparameterize anyway
		   false, // 'approx_orig' 
		   false, // 'close_belt' 
		   0,     // 'nmb_stabil' 
		   false); // 'repar' 

  asurf.setSmoothingWeight(smoothweight_);
  asurf.setDoRefine(false);
  asurf.setFixBoundary(false);
  int max_iter = 1;
  shared_ptr<SplineSurface> result_surf = asurf.getApproxSurf(maxdist_, avdist_, outsideeps_, 
							      max_iter);

  // Make LR spline surface
  double knot_tol = 1.0e-6;

  srf_ = shared_ptr<LRSplineSurface>(new LRSplineSurface(result_surf.get(), knot_tol));
  coef_known_.assign(srf_->numBasisFunctions(), 0.0);  // Initially nothing is fixed
}

//==============================================================================
void LRSurfApprox::makeInitSurf(shared_ptr<SplineSurface> surf)
//==============================================================================
{
  // Make LR spline surface
  double knot_tol = 1.0e-6;

  srf_ = shared_ptr<LRSplineSurface>(new LRSplineSurface(surf.get(), knot_tol));
  coef_known_.assign(srf_->numBasisFunctions(), 0.0);  // Initially nothing is fixed
}

//==============================================================================
void LRSurfApprox::computeParDomain(int dim, double& umin, double& umax, 
				    double& vmin, double& vmax)
//==============================================================================
{
  // Compute domain
  umin = umax = points_[0];
  vmin = vmax = points_[1];
  int del = 2 + dim;
  for (size_t ki=del; ki<points_.size(); ki+=del)
    {
      umin = std::min(umin, points_[ki]);
      umax = std::max(umax, points_[ki]);
      vmin = std::min(vmin, points_[ki+1]);
      vmax = std::max(vmax, points_[ki+1]);
    }
}

//==============================================================================
void LRSurfApprox::setCoefKnown()
//==============================================================================
{
  // Set coef_known from boundary information
  // Note that k-multiple knots at boundaries are expected, otherwise no
  // coefficients will be fixed
  
  for (size_t ki=0; ki<4; ++ki)
    {
      if (edge_derivs_[ki] == 0)
	continue;  // No coefficients to fix
      
      // Set boundary characteristica
      Direction2D d = (ki == 1 || ki == 3) ? XFIXED : YFIXED;  // Orthogonal to the curve
      Direction2D d2 = (ki == 1 || ki == 3) ? YFIXED : XFIXED; // Along the curve
      bool atstart = (ki == 0 || ki == 3);  // Whether the constant
      // parameter curve is in the start related to the opposite parameter direction
      // of the surface

      // Define coefficients as fixed in the associated b-splines
      int fixcoef = 1;
      setCoefKnown(d, d2, atstart, fixcoef);
    }

  // Transfer information to the global vector
  updateCoefKnown();
}

//==============================================================================
void LRSurfApprox::setCoefKnown(Direction2D d, Direction2D d2, 
				bool atstart, int fixcoef)
//==============================================================================
{
#ifdef DEBUG
  std::ofstream of("mesh.eps");
  writePostscriptMesh(*srf_, of);
#endif

  // // Check if the basis has k-tupple knots in the current direction
  // const Mesh2D& mesh = srf_->mesh();
  // int ix = (atstart) ? mesh.firstMeshVecIx(d) : 
  //   mesh.lastMeshVecIx(d);  
  // int mult = mesh.minMultInLine(d, ix);
  // int deg = srf_->degree(d);  // Degree orthogonal to the constant parameter curve
  // if (mult < deg)
  //   return;  // No fixed coefficients are set

  // // Fetch LRBSplines. Traverse the knot vector to make keys 
  // // for the bspline map.
  // int dir = d;
  // int startmult[2], endmult[2];
  // double startval[2], endval[2];
  // // int num = mesh.numDistinctKnots(d);  // Number of knots orthogonal to the
  // // // constant parameter curve
  // int num = mesh.numDistinctKnots(d); //mesh.numDistinctKnots(d2);  // Number of knots orthogonal to the
  // // constant parameter curve

  // // Set parameter value at the constant parameter curve
  // if (atstart)
  //   {
  //     startval[dir] = mesh.kval(d, 0);
  //   }
  // else
  //   {
  //     endval[dir] = mesh.kval(d, mesh.numDistinctKnots(d)-1);
  //   }
  // int deg2 = srf_->degree(d2);  // Degree along the constant parameter curve

  // vector<int> knot_idx =  LRBSpline2DUtils::derive_knots(mesh, d2, 
  // 							 mesh.firstMeshVecIx(d2),
  // 							 mesh.lastMeshVecIx(d2),
  // 							 atstart ? ix : ix-1,
  // 							 atstart ? ix+1 : ix);
  // // vector<int> knot_idx =  LRBSpline2DUtils::derive_knots(mesh, d, 
  // // 							 mesh.firstMeshVecIx(d2),
  // // 							 mesh.lastMeshVecIx(d2),
  // // 							 atstart ? ix : ix-1,
  // // 							 atstart ? ix+1 : ix);

  // // Since we have multiple knots in the surface boundary and we are only interested
  // // in the basis functions being non-zero along the boundary, we know that these
  // // basis functions must have a multiplicity equal to the order at the boundary
  // // and one in the other end.
  // startmult[dir] = atstart ? deg+1 : 1;
  // endmult[dir] = atstart ? 1 : deg+1;
  // size_t k1, k2;
  // vector<double> coefs;
  // for (k1=0, k2=deg2+1; k2<knot_idx.size(); ++k1, ++k2)
  //   {
  //     // Fetch domain of first minimal LR B-spline along the constant
  //     // parameter curve
  //     // First orthogonal to the curve
  //     if (atstart)
  // 	{
  // 	  int k_idx = 
  // 	    Mesh2DUtils::search_upwards_for_nonzero_multiplicity(mesh, d, 1, 
  // 								 knot_idx[k1], knot_idx[k2]);
  // 	  endval[dir] = mesh.kval(d, k_idx);
  // 	}
  //     else
  // 	{
  // 	  int k_idx = 
  // 	    Mesh2DUtils::search_downwards_for_nonzero_multiplicity(mesh, d, num-2, 
  // 								   knot_idx[k1], knot_idx[k2]);
  // 	  startval[dir] = mesh.kval(d, k_idx);
  // 	}

  //     // Along the curve
  //     startval[1-dir] = mesh.kval(d2, knot_idx[k1]);
  //     endval[1-dir] = mesh.kval(d2, knot_idx[k2]);

  //     // Count multiplicitity along the curve
  //     int km = 1;
  //     for (; km<=deg2; ++km)
  // 	if (knot_idx[k1+km] > knot_idx[k1])
  // 	  break;
  //     startmult[1-dir] = km;

  //     km = 1;
  //     for (; km<=deg2; ++km)
  // 	if (knot_idx[k2-km] < knot_idx[k2])
  // 	  break;
  //     endmult[1-dir] = km;

  //     // Fetch the associated LR B-spline
  //     LRSplineSurface::BSplineMap::iterator bm = srf_->bsplineFromDomain(startval[0], startval[1], 
  // 									 endval[0], endval[1], 
  // 									 startmult[0], startmult[1], 
  // 									 endmult[0], endmult[1]);

  //     // Set fixed coefficient information in the b-spline
  //     bm->second->setFixCoef(fixcoef);
  //   }

  // Traverse all B-splines and check whether they have maximum multiplicity along
  // the given edge
  for (LRSplineSurface::BSplineMap::iterator it=srf_->basisFunctionsBeginNonconst(); 
       it != srf_->basisFunctionsEndNonconst(); ++it)
    {
      int deg = it->second->degree(d);
      int mult = (d == XFIXED) ? it->second->endmult_u(atstart) :
	it->second->endmult_v(atstart);
      if (mult == deg+1)
	it->second->setFixCoef(fixcoef);
    }
  
}

//==============================================================================
void LRSurfApprox::updateCoefKnown()
//==============================================================================
{
  coef_known_.resize(srf_->numBasisFunctions());

  LRSplineSurface::BSplineMap::const_iterator it_bs;			
  size_t ki;
  for (it_bs=srf_->basisFunctionsBegin(), ki=0; it_bs!=srf_->basisFunctionsEnd(); 
       ++it_bs, ++ki)
    coef_known_[ki] = it_bs->second->coefFixed();
}

//==============================================================================
void LRSurfApprox::unsetCoefKnown()
//==============================================================================
{
  coef_known_.resize(srf_->numBasisFunctions());

  LRSplineSurface::BSplineMap::const_iterator it_bs;			
  size_t ki;
  for (it_bs=srf_->basisFunctionsBegin(), ki=0; it_bs!=srf_->basisFunctionsEnd(); 
       ++it_bs, ++ki)
    {
      coef_known_[ki] = 0;
      it_bs->second->setFixCoef(0);
    }
}

//==============================================================================
void LRSurfApprox::checkFeasibleRef(Element2D* elem, 
				    vector<LRSplineSurface::Refinement2D>& refs,
				    vector<Element2D*>& affected)
//==============================================================================
{
  // Fetch B-splines
  const vector<LRBSpline2D*>& bsplines = elem->getSupport();
  size_t nmb = bsplines.size();
  
 // Refine one B-spline with support in the parent element in one or two
  // parameter directions depending on how many elements with a large error
  // that lies in its support
  // First check refinement in the u-direction
  double u_par = 0.5*(elem->umin() + elem->umax());
  int ixu = -1;
  double max_wgt = 0.0;
  size_t ki, kj;
  for (ki=0; ki<nmb; ++ki)
    {
      // Count the number of elements with large error affected
      double curr_wgt = 0.0;
      vector<Element2D*> curr_el = bsplines[ki]->supportedElements();
      for (kj=0; kj<curr_el.size(); ++kj)
	{
	  if (curr_el[kj]->umax() < u_par || curr_el[kj]->umin() > u_par)
	    continue;  // Element not affected

	  // Compute weight for importance of refinement
	  double max_err, av_err;
	  int nmb_outside;
	  curr_el[kj]->getAccuracyInfo(av_err, max_err, nmb_outside);
	  int nmb_pts = curr_el[kj]->nmbDataPoints();
	  double wgt = av_err*(int)nmb_outside/(int)nmb_pts;
	  curr_wgt += wgt;
	}
      if (curr_wgt > max_wgt)
	{
	  max_wgt = curr_wgt;
	  ixu = (int)ki;
	}
    }
	  
  // The v-direction
  double v_par = 0.5*(elem->vmin() + elem->vmax());
  int ixv = -1;
  max_wgt = 0.0;
  for (ki=0; ki<nmb; ++ki)
    {
      // Count the number of elements with large error affected
      double curr_wgt = 0.0;
      vector<Element2D*> curr_el = bsplines[ki]->supportedElements();
      for (kj=0; kj<curr_el.size(); ++kj)
	{
	  if (curr_el[kj]->vmax() < v_par || curr_el[kj]->vmin() > v_par)
	    continue;  // Element not affected

	  // Compute weight for importance of refinement
	  double max_err, av_err;
	  int nmb_outside;
	  curr_el[kj]->getAccuracyInfo(av_err, max_err, nmb_outside);
	  int nmb_pts = curr_el[kj]->nmbDataPoints();
	  double wgt = av_err*(int)nmb_outside/(int)nmb_pts;
	  curr_wgt += wgt;
	}
      if (curr_wgt > max_wgt)
	{
	  max_wgt = curr_wgt;
	  ixv = (int)ki;
	}
    }
  
  // Fetch the affected elements in both parameter directions
  double fac = 0.1;
  double fac3 = 0.95;
  vector<Element2D*> aff_u;
  int nmb_u = 0;
  if (ixu >= 0)
    {
      vector<Element2D*> curr_el_u = bsplines[ixu]->supportedElements();
      for (kj=0; kj<curr_el_u.size(); ++kj)
	{
	  if (curr_el_u[kj]->umax() > u_par && curr_el_u[kj]->umin() < u_par)
	    {
	      double max_err, av_err;
	      int nmb_outside;
	      curr_el_u[kj]->getAccuracyInfo(av_err, max_err, nmb_outside);
	      int nmb_pts = curr_el_u[kj]->nmbDataPoints();
	      aff_u.push_back(curr_el_u[kj]);
	      if (nmb_outside > fac*nmb_pts || av_err > fac3*avdist_)
		nmb_u++;
	    }
	}
    }

  vector<Element2D*> aff_v;
  int nmb_v = 0;
  if (ixv >= 0)
    {
      vector<Element2D*> curr_el_v = bsplines[ixv]->supportedElements();
      for (kj=0; kj<curr_el_v.size(); ++kj)
	{
	  if (curr_el_v[kj]->vmax() > v_par && curr_el_v[kj]->vmin() < v_par)
	    {
	      double max_err, av_err;
	      int nmb_outside;
	      curr_el_v[kj]->getAccuracyInfo(av_err, max_err, nmb_outside);
	      int nmb_pts = curr_el_v[kj]->nmbDataPoints();
	      aff_v.push_back(curr_el_v[kj]);
	      if (nmb_outside > fac*nmb_pts || av_err > fac3*avdist_)
		nmb_v++;
	    }
	}
    }

  // Assemble information
  double fac2 = 0.5;
  std::set<Element2D*> affected_combined;
  if (ixu >= 0 && nmb_u > (int)fac2*bsplines[ixu]->degree(XFIXED))
    {
      affected_combined.insert(aff_u.begin(), aff_u.end());
      LRSplineSurface::Refinement2D curr_ref;
      curr_ref.setVal(u_par, bsplines[ixu]->vmin(), bsplines[ixu]->vmax(),
		      XFIXED, 1);
      refs.push_back(curr_ref);
    }
			       
    if (ixv >= 0 && nmb_v > (int)fac2*bsplines[ixv]->degree(YFIXED))
    {
      affected_combined.insert(aff_v.begin(), aff_v.end());
      LRSplineSurface::Refinement2D curr_ref;
      curr_ref.setVal(v_par, bsplines[ixv]->umin(), bsplines[ixv]->umax(),
		      YFIXED, 1);
      refs.push_back(curr_ref);
    }

    affected.insert(affected.end(), affected_combined.begin(), affected_combined.end());
}

//==============================================================================
void LRSurfApprox::turnTo3D()
//==============================================================================
{
  if (!srf_->dimension() == 1)
    return;  // Not possible to make 3D surface

  // Make data points 3D
  for (LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
       it != srf_->elementsEnd(); ++it)
    {
      it->second->makeDataPoints3D();
    }
  
  // Turn surface into 3D
  srf_->to3D();
}
