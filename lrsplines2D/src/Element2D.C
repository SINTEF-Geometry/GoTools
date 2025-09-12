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

#include "GoTools/lrsplines2D/Element2D.h"
#include "GoTools/lrsplines2D/LRBSpline2D.h"
#include <set>

using std::vector;
using std::set;



namespace {
  int el_compare_u_par(const void* el1, const void* el2)
  {
    if (((double*)el1)[0] < ((double*)el2)[0])
      return -1;
    else if (((double*)el1)[0] > ((double*)el2)[0])
      return 1;
    else
      return 0;
  }

  int el_compare_v_par(const void* el1, const void* el2)
  {
    if (((double*)el1)[1] < ((double*)el2)[1])
      return -1;
    else if (((double*)el1)[1] > ((double*)el2)[1])
      return 1;
    else
      return 0;
  }
}



namespace Go {

Element2D::Element2D() {
	start_u_ =  0;
	start_v_ =  0;
	stop_u_  =  0;
	stop_v_  =  0;
	//overloadCount_ = 0;
	is_modified_ = false;
	overload_ = false;

}

Element2D::Element2D(double start_u, double start_v, double stop_u, double stop_v) {

	start_u_ = start_u;
	start_v_ = start_v;
	stop_u_  = stop_u ;
	stop_v_  = stop_v ;
	//overloadCount_ = 0;
	is_modified_ = false;
	overload_ = false;

}

Element2D::~Element2D()
{
}

bool Element2D::resetOverload()
{
  int nmb = 0;
  for (size_t ki=0; ki<support_.size(); ++ki)
    if (support_[ki]->getOverload())
      {
	++nmb;
	if (nmb == 2)
	  break;
      }

  overload_ = (nmb >= 2);
  if (!overload_)
    {
      for (size_t ki=0; ki<support_.size(); ++ki)
	support_[ki]->eraseOverload();
    }
  return overload_;
}


void Element2D::removeSupportFunction(LRBSpline2D *f) {
  for (size_t i=0; i<support_.size(); i++) {
#ifndef NDEBUG
//      std::cout << "DEBUG: support_ i: " << i << std::endl;
#endif
    //if((support_[i]) && (*f == *support_[i])) {
      if((support_[i]) && (f == support_[i])) {
			support_[i] = support_.back();
			//support_[support_.size()-1] = NULL;
			support_.pop_back();
			return;
		}
	}
  is_modified_ = true;
}

void Element2D::addSupportFunction(LRBSpline2D *f) 
{
  for (size_t i=0; i<support_.size(); i++) 
    {
      if(f == support_[i]) 
	{
	  return;
	}
      if (*f == *support_[i])
      	{ // @@sbr I guess this is the correct solution, since we may update the element with a newer basis function.
	  //	  MESSAGE("DEBUG: We should avoid adding basis functions with the exact same support ...");
      	  support_[i] = f;
      	  return;
      	}
    }
  support_.push_back(f);
  // f->addSupport(this);
  is_modified_ = true;
}

bool Element2D::hasSupportFunction(LRBSpline2D *f) 
{
  for (size_t i=0; i<support_.size(); i++) {
    if(f == support_[i]) 
      return true;
  }
  return false;
}

Element2D* Element2D::copy()
	  {
	    Element2D *returnvalue = new Element2D();
	    
	    returnvalue->start_u_ = this->start_u_;
	    returnvalue->start_v_ = this->start_v_;
	    returnvalue->stop_u_ = this->stop_u_;
	    returnvalue->stop_v_ = this->stop_v_;
	    
	    // std::vector<LRBSpline2D*>::const_iterator it;
	    //for(it=element_[iEl]->supportBegin(); it<element_[iEl]->supportEnd(); it++)
	    //	      {
	    //	returnvalue -> support_.push_back(this->support_[i]->copy());
	    // }
	    /*
	    for(int i=0;i<support_.size(); i++)
	      {
		returnvalue -> support_.push_back(this->support_[i]->copy());
	      }
	    */
	    
	    return returnvalue;
	  }




Element2D* Element2D::split(bool split_u, double par_value) {
	Element2D *newElement2D = NULL;
	if(split_u) {
		if(par_value >= stop_u_ || par_value <= start_u_)
			return NULL;
		newElement2D = new Element2D(par_value, start_v_, stop_u_, stop_v_);
		stop_u_ = par_value;
	} else {
		if(par_value >= stop_v_ || par_value <= start_v_)
			return NULL;
		newElement2D = new Element2D(start_u_, par_value, stop_u_, stop_v_);
		stop_v_ = par_value;
	}
	for(size_t i=0; i<support_.size(); i++) {
		if(support_[i]->addSupport(newElement2D)) // tests for overlapping as well
			newElement2D->addSupportFunction(support_[i]);
		if(!support_[i]->overlaps(this)) {
			support_[i]->removeSupport(this);
			support_[i] = support_.back();
			support_.pop_back();
			i--;
		}
	}
	is_modified_ = true;
	newElement2D->setModified();
	return newElement2D;
}

void Element2D::updateBasisPointers(std::vector<LRBSpline2D*> &basis) {
	for(size_t i=0; i<support_.size(); i++) {
		// add pointer from LRBSpline2D back to Element2D
		support_.back()->addSupport(this);
	}
	is_modified_ = true;
}

void Element2D::swapParameterDirection()
{
    std::swap(start_u_, start_v_);
    std::swap(stop_u_, stop_v_);
    is_modified_ = true;
}

bool Element2D::isOverloaded()  {
  int n = (int)support_.size();
	if(n > 0) {
	  int p1 = support_.front()->degree(YFIXED) + 1;
	  int p2 = support_.front()->degree(XFIXED) + 1;

	  if(n > p1*p2)
	    return true;
	}
	return false;
}

bool Element2D::isOverloaded(int lowest_nmb)  {
  int n = (int)support_.size();
  if(n > lowest_nmb)
    return true;
  else
    return false;
}


  void Element2D::fetchNeighbours(vector<Element2D*>& neighbours) const
  {
    neighbours.clear();
    set<Element2D*> ngh;
    for (size_t ki=0; ki<support_.size(); ++ki)
      {
	vector<Element2D*>::iterator el;
	vector<Element2D*>::iterator last = support_[ki]->supportedElementEnd();
	for (el=support_[ki]->supportedElementBegin(); el!=last; ++el)
	  {
	    if ((*el) == this)
	      continue;
	    double u1 = (*el)->umin();
	    double u2 = (*el)->umax();
	    double v1 = (*el)->vmin(); 
	    double v2 = (*el)->vmax();
	    if (start_u_ > u2 || stop_u_ < u1 || start_v_ > v2 ||
		stop_v_ < v1)
	      continue;
	    if ((start_u_ == u2 && (start_v_ == v2 || stop_v_ == v1)) ||
		(stop_u_ == u1 && (start_v_ == v2 || stop_v_ == v1)))
	      continue;
	    ngh.insert(*el);
	  }
      }
    neighbours.insert(neighbours.end(), ngh.begin(), ngh.end());
  }


  int Element2D::nmbDataPoints()
  {
    if (LSdata_.get())
      {
	int dim =  (support_.size() == 0) ? 1 : support_[0]->dimension();
	int nmb = LSdata_->getNmbValPrPoint();
	if (nmb == 0)
	  nmb = 3+dim;
	return (LSdata_->dataPointSize()/nmb);
      }
    else
      return 0;
  }

  int Element2D::nmbSignificantPoints()
  {
    if (LSdata_.get())
      {
	int dim =  (support_.size() == 0) ? 1 : support_[0]->dimension();
	int nmb = LSdata_->getNmbValPrPoint();
	if (nmb == 0)
	  nmb = 3+dim;
	return (LSdata_->significantPointSize()/nmb);
      }
    else
      return 0;
  }

  int Element2D::nmbGhostPoints()
  {
    if (LSdata_.get())
      {
	int dim =  (support_.size() == 0) ? 1 : support_[0]->dimension();
	int nmb = LSdata_->getNmbValPrPoint();
	if (nmb == 0)
	  nmb = 3+dim;
	return (LSdata_->ghostPointSize()/nmb);
      }
    else
      return 0;
  }

void Element2D::getOutsidePoints(vector<double>& points, Direction2D d,
				 bool& sort_in_u)
  {
    if (LSdata_)
      {
	int dim =  (support_.size() == 0) ? 1 : support_[0]->dimension();
	double start = (d == XFIXED) ? start_u_ : start_v_;
	double end = (d == XFIXED) ? stop_u_ : stop_v_;
	LSdata_->getOutsidePoints(points, dim, d, start, end, sort_in_u);
      }
  }

void Element2D::getOutsideSignificantPoints(vector<double>& points, 
					    Direction2D d, bool& sort_in_u)
  {
    if (LSdata_)
      {
	int dim =  (support_.size() == 0) ? 1 : support_[0]->dimension();
	double start = (d == XFIXED) ? start_u_ : start_v_;
	double end = (d == XFIXED) ? stop_u_ : stop_v_;
	LSdata_->getOutsideSignificantPoints(points, dim, d, start, end, 
					     sort_in_u);
      }
  }

  void Element2D::getOutsideGhostPoints(vector<double>& points, Direction2D d,
					bool& sort_in_u)
  {
    if (LSdata_)
      {
	int dim =  (support_.size() == 0) ? 1 : support_[0]->dimension();
	double start = (d == XFIXED) ? start_u_ : start_v_;
	double end = (d == XFIXED) ? stop_u_ : stop_v_;
	LSdata_->getOutsideGhostPoints(points, dim, d, start, end, sort_in_u);
      }
  }

  void Element2D::updateLSDataParDomain(double u1, double u2, 
					double v1, double v2, 
					double u1new, double u2new, 
					double v1new, double v2new)
  {
    if (LSdata_.get())
      {
	int dim =  (support_.size() == 0) ? 1 : support_[0]->dimension();
	LSdata_->updateLSDataParDomain(u1, u2, v1, v2, u1new, 
				       u2new, v1new, v2new, dim);
      }
  }

  void Element2D::setLSMatrix()
  {
    // Count the number of coefficients that are not fixed
    int nmb=0;
    for (size_t ki=0; ki<support_.size(); ++ki)
      if (!support_[ki]->coefFixed())
	nmb++;

    // Get dimension of geometry space
    int dim = (support_.size() == 0) ? 1 : support_[0]->dimension();

    // Make scratch
    if (!LSdata_)
      LSdata_ = shared_ptr<LSSmoothData>(new LSSmoothData());
    LSdata_->setLSMatrix(nmb, dim);
    
  }

  bool Element2D::getDataBoundingBox(double bb[])
  {
    if (LSdata_.get())
      return LSdata_->getDataBoundingBox(support_[0]->dimension(), bb);
    else
      return false;
  }

  void Element2D::makeDataPoints3D()
  {
    if (LSdata_.get())
      {
	int dim =  (support_.size() == 0) ? 1 : support_[0]->dimension();
	if (dim != 1)
	  return;
	LSdata_->makeDataPoints3D(dim);
	is_modified_ = true;
      }
  }

  void Element2D::updateAccuracyInfo()
  {
    if (LSdata_.get())
      {
	int dim =  (support_.size() == 0) ? 1 : support_[0]->dimension();
	LSdata_->updateAccuracyInfo(dim);
      }
  }

double Element2D::sumOfScaledBsplines(double upar, double vpar)
{
  double val = 0.0;
  for (size_t ki=0; ki<support_.size(); ++ki)
    {
      double curr = support_[ki]->evalBasisFunction(upar, vpar);
      val += curr*support_[ki]->gamma();
    }
  return val;
}


  void LSSmoothData::getOutsidePoints(vector<double>& points, int dim,
				      Direction2D d, double start, double end,
				      bool& sort_in_u)
  {
    // Sort the points in the indicated direction
    int del = (pt_del_ > 0) ? pt_del_ : dim+3;   // Number of entries for each point
    int nmb = (int)data_points_.size()/del;  // Number of data points
    if (nmb == 0)
      return;  // No points to sort

    int ix = (d == XFIXED) ? 0 : 1;
    if (true)
      //dim > 1 || (d == XFIXED && !sort_in_u_) || (d == YFIXED && sort_in_u_))
      {
	qsort(&data_points_[0], nmb, del*sizeof(double), 
	      (d == XFIXED) ? el_compare_u_par : el_compare_v_par);
	sort_in_u_ = !sort_in_u_;
      }
    
    // Traverse point set and extract inside and outside sub sets
    vector<double>::iterator first1;
    vector<double>::iterator last1;
    vector<double>::iterator first2;
    vector<double>::iterator last2;
    if (data_points_[ix] >= start)
      {
	first1 = data_points_.begin();
	int ki;
	for (ki=0; ki<(int)data_points_.size(); ki+=del)
	  if (data_points_[ki+ix] > end)
	    break;
	last1 = first2 = data_points_.begin() + ki;
	last2 = data_points_.end();
      }
    else
      {
	first2 = data_points_.begin();
	int ki;
	for (ki=0; ki<(int)data_points_.size(); ki+=del)
	  if (data_points_[ki+ix] > end)
	    break;
	last2 = first1 = data_points_.begin() + ki;
	last1 = data_points_.end();
      }
    
    // Split vector
    points.insert(points.end(), first2, last2);
    data_points_.erase(first2, last2);
    sort_in_u = sort_in_u_;
  }

  void LSSmoothData::getOutsideSignificantPoints(vector<double>& points, 
						 int dim, Direction2D d, 
						 double start, double end,
						 bool& sort_in_u)
  {
    // Sort the points in the indicated direction
    int del = (pt_del_ > 0) ? pt_del_ : dim+3;   // Number of entries for each point
    int nmb = (int)significant_points_.size()/del;  // Number of data points
    if (nmb == 0)
      return;  // No points to sort

    int ix = (d == XFIXED) ? 0 : 1;
    if (true)
      //dim > 1 || (d == XFIXED && !sort_in_u_) || (d == YFIXED && sort_in_u_))
      {
	qsort(&significant_points_[0], nmb, del*sizeof(double), 
	      (d == XFIXED) ? el_compare_u_par : el_compare_v_par);
	sort_in_u_ = !sort_in_u_;
      }
    
    // Traverse point set and extract inside and outside sub sets
    vector<double>::iterator first1;
    vector<double>::iterator last1;
    vector<double>::iterator first2;
    vector<double>::iterator last2;
    if (significant_points_[ix] >= start)
      {
	first1 = significant_points_.begin();
	int ki;
	for (ki=0; ki<(int)significant_points_.size(); ki+=del)
	  if (significant_points_[ki+ix] > end)
	    break;
	last1 = first2 = significant_points_.begin() + ki;
	last2 = significant_points_.end();
      }
    else
      {
	first2 = significant_points_.begin();
	int ki;
	for (ki=0; ki<(int)significant_points_.size(); ki+=del)
	  if (significant_points_[ki+ix] > end)
	    break;
	last2 = first1 = significant_points_.begin() + ki;
	last1 = significant_points_.end();
      }
    
    // Split vector
    points.insert(points.end(), first2, last2);
    significant_points_.erase(first2, last2);
    sort_in_u = sort_in_u_;
  }

  void LSSmoothData::getOutsideGhostPoints(vector<double>& points, int dim,
					   Direction2D d, double start, 
					   double end, bool& sort_in_u)
  {
    // Sort the points in the indicated direction
    int del = (pt_del_ > 0) ? pt_del_ : dim+3;   // Number of entries for each point
    int nmb = (int)ghost_points_.size()/del;  // Number of data points
    int ix = (d == XFIXED) ? 0 : 1;
    if (true)
      //dim > 1 || (d == XFIXED && !sort_in_u_ghost_) || 
      //(d == YFIXED && sort_in_u_ghost_))
      {
	qsort(&ghost_points_[0], nmb, del*sizeof(double), 
	      (d == XFIXED) ? el_compare_u_par : el_compare_v_par);
	sort_in_u_ghost_ = !sort_in_u_ghost_;
      }
    
    if (nmb == 0)
      return;  // No points to sort

    // Traverse point set and extract inside and outside sub sets
    vector<double>::iterator first1;
    vector<double>::iterator last1;
    vector<double>::iterator first2;
    vector<double>::iterator last2;
    if (ghost_points_[ix] >= start)
      {
	first1 = ghost_points_.begin();
	int ki;
	for (ki=0; ki<(int)ghost_points_.size(); ki+=del)
	  if (ghost_points_[ki+ix] > end)
	    break;
	last1 = first2 = ghost_points_.begin() + ki;
	last2 = ghost_points_.end();
      }
    else
      {
	first2 = ghost_points_.begin();
	int ki;
	for (ki=0; ki<(int)ghost_points_.size(); ki+=del)
	  if (ghost_points_[ki+ix] > end)
	    break;
	last2 = first1 = ghost_points_.begin() + ki;
	last1 = ghost_points_.end();
      }
    
    // Split vector
    points.insert(points.end(), first2, last2);
    ghost_points_.erase(first2, last2);
    sort_in_u = sort_in_u_ghost_;
  }

  void LSSmoothData::makeDataPoints3D(int dim)
  {
    int del1 = (pt_del_ > 0) ? pt_del_ : dim+3;   // Number of entries for each point
    int nmb = (int)data_points_.size()/del1;
    int del2 = 2+del1;
    vector<double> points(del2*nmb);  // Parameter value + point + distance
    for (int ki=0; ki<nmb; ++ki)
      {
	points[del2*ki] = data_points_[del1*ki];
	points[del2*ki+1] = data_points_[del1*ki+1];
	for (int kj=0; kj<del1; ++kj)
	  points[del2*ki+2+kj] = data_points_[del1*ki+kj];
      }
    std::swap(data_points_, points);

    nmb = (int)significant_points_.size()/del1;
    vector<double> spoints(del2*nmb);  // Parameter value + point
    for (int ki=0; ki<nmb; ++ki)
      {
	spoints[del2*ki] = significant_points_[del1*ki];
	spoints[del2*ki+1] = significant_points_[del1*ki+1];
	for (int kj=0; kj<del1; ++kj)
	  spoints[del2*ki+2+kj] = significant_points_[del1*ki+kj];
      }
    std::swap(significant_points_, spoints);
 
    nmb = (int)ghost_points_.size()/del1;
    vector<double> gpoints(del2*nmb);  // Parameter value + point
    for (int ki=0; ki<nmb; ++ki)
      {
	gpoints[del2*ki] = ghost_points_[del1*ki];
	gpoints[del2*ki+1] = ghost_points_[del1*ki+1];
	for (int kj=0; kj<del1; ++kj)
	  gpoints[del2*ki+2+kj] = ghost_points_[del1*ki+kj];
      }
    std::swap(ghost_points_, gpoints);
   }

  void LSSmoothData::updateAccuracyInfo(int dim)
  {
    accumulated_error_ = 0.0;
    average_error_ = 0.0;
    max_error_ = -1.0;

    int del = (pt_del_ > 0) ? pt_del_ : dim+3;   // Number of entries for each point
    int ix = (del > dim+3) ? del-2 : del-1;
    int nmb = (int)data_points_.size()/del;
    for (int ki=0; ki<nmb; ++ki)
      {
	if (del > dim+3 && data_points_[ki*del+ix+1] < 0.0)
	  continue;
	double dist = data_points_[ki*del+ix];
	double dist2 = fabs(dist);
	max_error_ = std::max(max_error_, dist2);
	accumulated_error_ += dist2;
      }

    nmb = (int)significant_points_.size()/del;
    for (int ki=0; ki<nmb; ++ki)
      {
	double dist = significant_points_[ki*del+ix];
	double dist2 = fabs(dist);
	max_error_ = std::max(max_error_, dist2);
	accumulated_error_ += dist2;
      }
    average_error_ = -1.0; // No longer valid
    nmb_outside_tol_ = -1;
  }

  bool LSSmoothData::getDataBoundingBox(int dim, double bb[])
  {
    int del = (pt_del_ > 0) ? pt_del_ : dim+3;   // Number of entries for each point
    int nmb = (int)data_points_.size()/del;
    int nmb_sign = (int)significant_points_.size()/del;
    if (nmb+nmb_sign == 0)
      return false;
    if (dim == 1 && minheight_ <= maxheight_)
      {
	bb[0] = minheight_;
	bb[1] = maxheight_;
      }
    else
      {
	int ki, kj;
	if (nmb > 0)
	  {
	    for (kj=0; kj<dim; ++kj)
	      bb[2*kj] = bb[2*kj+1] = data_points_[2+kj];
	  }
	else if (nmb_sign > 0)
	  {
	    for (kj=0; kj<dim; ++kj)
	      bb[2*kj] = bb[2*kj+1] = significant_points_[2+kj];
	  }

	for (ki=1; ki<nmb; ++ki)
	  {
	    for (kj=0; kj<dim; ++kj)
	      {
		double val = data_points_[ki*del+2+kj];
		bb[2*kj] = std::min(bb[2*kj], val);
		bb[2*kj+1] = std::max(bb[2*kj+1], val);
	      }
	  }
	for (ki=0; ki<nmb_sign; ++ki)
	  {
	    for (kj=0; kj<dim; ++kj)
	      {
		double val = significant_points_[ki*del+2+kj];
		bb[2*kj] = std::min(bb[2*kj], val);
		bb[2*kj+1] = std::max(bb[2*kj+1], val);
	      }
	  }
      }
    return true;
  }

  void LSSmoothData::updateLSDataParDomain(double u1, double u2, 
					   double v1, double v2, 
					   double u1new, double u2new, 
					   double v1new, double v2new,
					   int dim)
  {
    int del = (pt_del_ > 0) ? pt_del_ : dim+3;   // Number of entries for each point
    double d1u = u2 - u1;
    double d2u = u2new - u1new;
    double d1v = v2 - v1;
    double d2v = v2new - v1new;
    size_t ki;
    for (ki=0; ki<data_points_.size(); ki+=del)
      {
	data_points_[ki] = (data_points_[ki]-u1)*d2u/d1u + u1new;
	data_points_[ki+1] = (data_points_[ki+1]-v1)*d2v/d1v + v1new;
      }
    for (ki=0; ki<significant_points_.size(); ki+=del)
      {
	significant_points_[ki] = (significant_points_[ki]-u1)*d2u/d1u + u1new;
	significant_points_[ki+1] = (significant_points_[ki+1]-v1)*d2v/d1v + v1new;
      }
    for (ki=0; ki<ghost_points_.size(); ki+=del)
      {
	ghost_points_[ki] = (ghost_points_[ki]-u1)*d2u/d1u + u1new;
	ghost_points_[ki+1] = (ghost_points_[ki+1]-v1)*d2v/d1v + v1new;
      }
  }



  vector<double> Element2D::unitSquareBernsteinBasis() const
  {
    vector<double> result;
    if (support_.size() == 0)
      return result;   // This should not happen, some B-splines must have support in this element

    int dim = support_[0]->dimension();
    int deg_u = support_[0]->degree(XFIXED);
    int deg_v = support_[0]->degree(YFIXED);
    int result_size = dim * (deg_u+1) * (deg_v+1);
    result.resize(result_size,0.0);

    // double inv_size_u = 1.0 / (stop_u_ - start_u_);
    // double inv_size_v = 1.0 / (stop_v_ - start_v_);

    // Add up the coefficients for each B-spline
    for (vector<LRBSpline2D*>::const_iterator it = support_.begin(); it != support_.end(); ++it)
      {
	vector<double> local_basis = (*it)->unitSquareBernsteinBasis(start_u_, stop_u_, start_v_, stop_v_);
	for (int i = 0; i < result_size; ++i)
	  result[i] += local_basis[i];
      }

    return result;
  }



  void Element2D::bernsteinEvaluation(int degree, double value, vector<vector<double> >& result) const
  {
    result.resize(degree + 1);
    result[0].resize(1);
    result[0][0] = 1.0;

    for (int i = 1; i <= degree; ++i)
      {
	result[i].resize(i + 1);
	for (int j = 0; j < i; ++j)
	  {
	    result[i][j + 1] += value * result[i - 1][j];
	    result[i][j] += (1.0 - value) * result[i - 1][j];
	  }
      }
  }


  void Element2D::univariateBernsteinEvaluationInLine(int degree, double start, double end, vector<vector<double> >& result) const
  {
    result.resize(degree + 1);
    for (int i = 0; i <= degree; ++i)
      result[i].resize(degree + 1);

    vector<vector<double> > bernstein_start;
    vector<vector<double> > bernstein_end;

    // Some special cases can be treated more efficiently
    if (start == 0.0)
      {
	if (end == 1.0)
	  for (int i = 0; i <= degree; ++i)
	    result[i][i] = 1.0;

	else
	  {
	    bernsteinEvaluation(degree, end, bernstein_end);
	    for (int i = 0; i <= degree; ++i)
	      for (int j = i; j <= degree; ++j)
		result[i][j] = bernstein_end[j][i];
	  }
      }

    else if (start == 1.0)
      {
	if (end == 0.0)
	  for (int i = 0; i <= degree; ++i)
	    result[degree - i][i] = 1.0;

	else
	  {
	    bernsteinEvaluation(degree, end, bernstein_end);
	    for (int i = 0; i <= degree; ++i)
	      for (int j = degree - i; j <= degree; ++j)
		result[i][j] = bernstein_end[j][i + j - degree];
	  }
      }

    else if (end == 0.0)
      {
	bernsteinEvaluation(degree, start, bernstein_start);
	for (int i = 0; i <= degree; ++i)
	  for (int j = 0; j <= degree - i; ++j)
	    result[i][j] = bernstein_start[degree - j][i];
      }

    else if (end == 1.0)
      {
	bernsteinEvaluation(degree, start, bernstein_start);
	for (int i = 0; i <= degree; ++i)
	  for (int j = 0; j <= i; ++j)
	    result[i][j] = bernstein_start[degree - j][i - j];
      }

    else  // The general case where start and end both differ from 0.0 and 1.0
      {
	bernsteinEvaluation(degree, start, bernstein_start);
	bernsteinEvaluation(degree, end, bernstein_end);
	for (int i = 0; i <= degree; ++i)
	  for (int j = 0; j <= degree; ++j)
	    {
	      int k_start = std::max(0, i + j - degree);
	      int k_end = std::min(i, j);
	      for (int k = k_start; k <= k_end; ++k)
		result[i][j] += bernstein_start[degree - j][i - k] * bernstein_end[j][k];
	    }
      }

    vector<int> binomial(degree + 1);   // binomial[i] shall be the binomial (degree Choose i)
    binomial[0] = binomial[degree] = 1;
    for (int i = 1; i + i <= degree; ++i)
      binomial[degree - i] = binomial[i] = (binomial[i - 1] * (degree - i + 1)) / i;

    for (int j = 0; j <= degree; ++j)
      {
	double binom_d_j = (double)binomial[j];
	for (int i = 0; i <= degree; ++i)
	  result[i][j] *= binom_d_j;
      }
  }


/*
int Element2D::overloadedBasisCount() const {
	int ans = 0;
	for(size_t i=0; i<support_.size(); i++)
		if(support_[i]->isOverloaded())
			ans++;
	return ans;
}
*/

} // end namespace Go
