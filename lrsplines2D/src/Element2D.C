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

using std::vector;

namespace Go {

Element2D::Element2D() {
	start_u_ =  0;
	start_v_ =  0;
	stop_u_  =  0;
	stop_v_  =  0;
	overloadCount_ = 0;
	is_modified_ = false;
}

Element2D::Element2D(double start_u, double start_v, double stop_u, double stop_v) {

	start_u_ = start_u;
	start_v_ = start_v;
	stop_u_  = stop_u ;
	stop_v_  = stop_v ;
	overloadCount_ = 0;
	is_modified_ = false;
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

bool Element2D::isOverloaded()  const {
  int n = (int)support_.size();
	if(n > 0) {
		int p1 = support_.front()->degree(YFIXED) + 1;
		int p2 = support_.front()->degree(XFIXED) + 1;
		if(n > p1*p2)
			return true;
	}
	return false;
}

  int Element2D::nmbDataPoints()
  {
    if (LSdata_.get())
      {
	int dim =  (support_.size() == 0) ? 1 : support_[0]->dimension();
	return (LSdata_->dataPointSize()/(2+dim));
      }
    else
      return 0;
  }

  int Element2D::nmbGhostPoints()
  {
    if (LSdata_.get())
      {
	int dim =  (support_.size() == 0) ? 1 : support_[0]->dimension();
	return (LSdata_->ghostPointSize()/(2+dim));
      }
    else
      return 0;
  }

  void Element2D::getOutsidePoints(vector<double>& points, Direction2D d)
  {
    if (LSdata_)
      {
	int dim =  (support_.size() == 0) ? 1 : support_[0]->dimension();
	double start = (d == XFIXED) ? start_u_ : start_v_;
	double end = (d == XFIXED) ? stop_u_ : stop_v_;
	LSdata_->getOutsidePoints(points, dim, d, start, end);
      }
  }

  void Element2D::getOutsideGhostPoints(vector<double>& points, Direction2D d)
  {
    if (LSdata_)
      {
	int dim =  (support_.size() == 0) ? 1 : support_[0]->dimension();
	double start = (d == XFIXED) ? start_u_ : start_v_;
	double end = (d == XFIXED) ? stop_u_ : stop_v_;
	LSdata_->getOutsideGhostPoints(points, dim, d, start, end);
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

int compare_u_par(const void* el1, const void* el2)
{
  if (((double*)el1)[0] < ((double*)el2)[0])
    return -1;
  else if (((double*)el1)[0] > ((double*)el2)[0])
    return 1;
  else
    return 0;
}

int compare_v_par(const void* el1, const void* el2)
{
  if (((double*)el1)[1] < ((double*)el2)[1])
    return -1;
  else if (((double*)el1)[1] > ((double*)el2)[1])
    return 1;
  else
    return 0;
}

  void LSSmoothData::getOutsidePoints(vector<double>& points, int dim,
				      Direction2D d, double start, double end)
  {
    // Sort the points in the indicated direction
    int del = dim+2;                   // Number of entries for each point
    int nmb = (int)data_points_.size()/del;  // Number of data points
    int ix = (d == XFIXED) ? 0 : 1;
    qsort(&data_points_[0], nmb, del*sizeof(double), 
	  (d == XFIXED) ? compare_u_par : compare_v_par);
    
    if (nmb == 0)
      return;  // No points to sort

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
  }

  void LSSmoothData::getOutsideGhostPoints(vector<double>& points, int dim,
					   Direction2D d, double start, double end)
  {
    // Sort the points in the indicated direction
    int del = dim+2;                   // Number of entries for each point
    int nmb = (int)ghost_points_.size()/del;  // Number of data points
    int ix = (d == XFIXED) ? 0 : 1;
    qsort(&ghost_points_[0], nmb, del*sizeof(double), 
	  (d == XFIXED) ? compare_u_par : compare_v_par);
    
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
  }

  void LSSmoothData::makeDataPoints3D(int dim)
  {
    int nmb = data_points_.size()/(2+dim);
    int del1 = 2+dim;
    int del2 = 2+del1;
    vector<double> points(del2*nmb);  // Parameter value + point
    for (int ki=0; ki<nmb; ++ki)
      {
	points[del2*ki] = data_points_[del1*ki];
	points[del2*ki+1] = data_points_[del1*ki+1];
	for (int kj=0; kj<del1; ++kj)
	  points[del2*ki+2+kj] = data_points_[del1*ki+kj];
      }
    std::swap(data_points_, points);

    nmb = ghost_points_.size()/(2+dim);
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
