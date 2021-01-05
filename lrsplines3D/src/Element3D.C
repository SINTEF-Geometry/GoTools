//===========================================================================
//                                                                           
// File: Element3D.C                                                         
//                                                                           
// Created: Fri Mar  8 15:02:31 2013                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/lrsplines3D/Element3D.h"
#include "GoTools/lrsplines3D/LRBSpline3D.h"

namespace Go {

Element3D::Element3D() {
	start_u_ =  0;
	start_v_ =  0;
	start_w_ =  0;
	stop_u_  =  0;
	stop_v_  =  0;
	stop_w_  =  0;
	overloadCount_ = 0;
	is_modified_ = false;
}

Element3D::Element3D(double start_u, double start_v, double start_w,
		     double stop_u, double stop_v, double stop_w)
{
    start_u_ = start_u;
    start_v_ = start_v;
    start_w_ = start_w;
    stop_u_  = stop_u ;
    stop_v_  = stop_v ;
    stop_w_  = stop_w ;
    overloadCount_ = 0;
    is_modified_ = false;
}

void Element3D::removeSupportFunction(LRBSpline3D *f) {
  for (size_t i=0; i<support_.size(); i++) {
      if((support_[i]) && (f == support_[i])) {
			support_[i] = support_.back();
			//support_[support_.size()-1] = NULL;
			support_.pop_back();
			is_modified_ = true;
			return;
		}
	}
  int stop_break = 1;
}

void Element3D::addSupportFunction(LRBSpline3D *f)
{
  for (size_t i=0; i<support_.size(); i++)
    {
      if(f == support_[i])
        {
          return;
        }
      if (*f == *support_[i])
        { // @@sbr I guess this is the correct solution, since we may update the element with a newer basis function.
          //      MESSAGE("DEBUG: We should avoid adding basis functions with the exact same support ...");
          support_[i] = f;
          return;
        }
    }
  support_.push_back(f);
  // f->addSupport(this);
  is_modified_ = true;
}

Element3D* Element3D::copy()
	  {
	    Element3D *returnvalue = new Element3D();
	    
	    returnvalue->start_u_ = this->start_u_;
	    returnvalue->start_v_ = this->start_v_;
	    returnvalue->start_w_ = this->start_w_;
	    returnvalue->stop_u_ = this->stop_u_;
	    returnvalue->stop_v_ = this->stop_v_;
	    returnvalue->stop_w_ = this->stop_w_;
	    
	    // std::vector<LRBSpline3D*>::const_iterator it;
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




Element3D* Element3D::split(Direction3D split_dir, double par_value)
{
  Element3D *newElement3D = NULL;
  if (split_dir == XDIR)
    {
      if(par_value >= stop_u_ || par_value <= start_u_)
	return NULL;
      newElement3D = new Element3D(par_value, start_v_, start_w_, stop_u_, stop_v_, stop_w_);
      stop_u_ = par_value;
    } else if (split_dir == YDIR)
    {
      if(par_value >= stop_v_ || par_value <= start_v_)
	return NULL;
      newElement3D = new Element3D(start_u_, par_value, start_w_, stop_u_, stop_v_, stop_w_);
      stop_v_ = par_value;
    }
  else
    {
      if(par_value >= stop_w_ || par_value <= start_w_)
	return NULL;
      newElement3D = new Element3D(start_u_, start_v_, par_value, stop_u_, stop_v_, stop_w_);
      stop_w_ = par_value;
    }

  for(size_t i=0; i<support_.size(); i++)
    {
      if(support_[i]->addSupport(newElement3D)) // tests for overlapping as well
	newElement3D->addSupportFunction(support_[i]);
      if(!support_[i]->overlaps(this))
	{
	  support_[i]->removeSupport(this);
	  support_[i] = support_.back();
	  support_.pop_back();
	  i--;
	}
    }
  is_modified_ = true;
  newElement3D->setModified();
  return newElement3D;
}

void Element3D::updateBasisPointers(std::vector<LRBSpline3D*> &basis) {
	for(size_t i=0; i<support_.size(); i++) {
		// add pointer from LRBSpline3D back to Element3D
		support_.back()->addSupport(this);
	}
	is_modified_ = true;
}

// void Element3D::swapParameterDirection()
// {
//     std::swap(start_u_, start_v_);
//     std::swap(stop_u_, stop_v_);
//     is_modified_ = true;
// }

bool Element3D::isOverloaded()  const {
  int n = (int)support_.size();
	if(n > 0) {
		int p1 = support_.front()->degree(ZDIR) + 1;
		int p2 = support_.front()->degree(YDIR) + 1;
		int p3 = support_.front()->degree(XDIR) + 1;
		if(n > p1*p2*p3)
			return true;
	}
	return false;
}

int Element3D::nmbDataPoints()
{
  if (approx_data_.get())
    {
      int dim =  (support_.size() == 0) ? 1 : support_[0]->dimension();
      return (approx_data_->dataPointSize()/(4+dim)); // @obar: is this related to del?
    }
  else
    return 0;
}

void Element3D::getOutsidePoints(std::vector<double>& points, Direction3D d,
                                 bool& sort_in_u)
  {
    if (approx_data_)
      {
        int dim =  (support_.size() == 0) ? 1 : support_[0]->dimension();
        double start = (d == XDIR) ? start_u_ : ((d == YDIR) ? start_v_ : start_w_);
        double end   = (d == XDIR) ? stop_u_ : ((d == YDIR) ? stop_v_ : stop_w_);
        approx_data_->getOutsidePoints(points, dim, d, start, end, sort_in_u);
      }
  }

void Element3D::updateApproxDataParDomain(double u1, double u2,
					  double v1, double v2,
					  double w1, double w2,
					  double u1new, double u2new,
					  double v1new, double v2new,
					  double w1new, double w2new)
{
  if (approx_data_.get())
    {
      int dim =  (support_.size() == 0) ? 1 : support_[0]->dimension();
      approx_data_->updateApproxDataParDomain(u1, u2, v1, v2, w1, w2,
                                              u1new, u2new, v1new, v2new, w1new, w2new,
                                              dim);
    }
}

void Element3D::updateAccuracyInfo()
{
  if (approx_data_.get())
    {
      int dim =  (support_.size() == 0) ? 1 : support_[0]->dimension();
      approx_data_->updateAccuracyInfo(dim);
    }
}

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

int el_compare_w_par(const void* el1, const void* el2)
{
  if (((double*)el1)[2] < ((double*)el2)[2])
    return -1;
  else if (((double*)el1)[2] > ((double*)el2)[2])
    return 1;
  else
    return 0;
}
}

/*
int Element3D::overloadedBasisCount() const {
	int ans = 0;
	for(size_t i=0; i<support_.size(); i++)
		if(support_[i]->isOverloaded())
			ans++;
	return ans;
}
*/

void Approx3DData::updateApproxDataParDomain(double u1, double u2,
					     double v1, double v2,
					     double w1, double w2,
					     double u1new, double u2new,
					     double v1new, double v2new,
					     double w1new, double w2new,
					     int dim)
{
  int del = 4+dim;  // Parameter triple, position and distance
  double d1u = u2 - u1;
  double d2u = u2new - u1new;
  double d1v = v2 - v1;
  double d2v = v2new - v1new;
  double d1w = w2 - w1;
  double d2w = w2new - w1new;
  size_t ki;
  for (ki; ki<data_points_.size(); ki+=del)
    {
      data_points_[ki] = (data_points_[ki]-u1)*d2u/d1u + u1new;
      data_points_[ki+1] = (data_points_[ki+1]-v1)*d2v/d1v + v1new;
      data_points_[ki+2] = (data_points_[ki+2]-w1)*d2w/d1w + w1new;
    }
}

void Approx3DData::getOutsidePoints(std::vector<double>& points, int dim,
				    Direction3D d, double start, double end,
				    bool& sort_in_u)
{
  // Sort the points in the indicated direction
  int del = dim+4;                         // Number of entries for each point
  int nmb = (int)data_points_.size()/del;  // Number of data points
  if (nmb == 0)
    return;  // No points to sort

  int ix = (d == XDIR) ? 0 : ((d == YDIR) ? 1 : 2);

  if (true)
    //dim > 1 || (d == XDIR && !sort_in_u_) || (d == YDIR && sort_in_u_))
    {
      qsort(&data_points_[0], nmb, del*sizeof(double),
          (d == XDIR) ? el_compare_u_par : ((d == YDIR) ? el_compare_v_par : el_compare_w_par));
      sort_in_u_ = !sort_in_u_;
    }

  // Traverse point set and extract inside and outside sub sets
  std::vector<double>::iterator first1;
  std::vector<double>::iterator last1;
  std::vector<double>::iterator first2;
  std::vector<double>::iterator last2;
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

void Approx3DData::updateAccuracyInfo(int dim)
{
  accumulated_error_ = 0.0;
  average_error_ = 0.0;
  max_error_ = -1.0;

  int del = 4+dim;  // Parameter pair, position and distance
  int nmb = data_points_.size()/del;
  for (int ki=0; ki<nmb; ++ki)
    {
      double dist = data_points_[ki*del+del-1];
      double dist2 = fabs(dist);
      max_error_ = std::max(max_error_, dist2);
      accumulated_error_ += dist2;
    }
  average_error_ = -1.0; // No longer valid
  nmb_outside_tol_ = -1;
}


} // end namespace LR
