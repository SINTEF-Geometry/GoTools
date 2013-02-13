
#include "GoTools/lrsplines2D/Element2D.h"
#include "GoTools/lrsplines2D/LRBSpline2D.h"

namespace Go {

Element2D::Element2D() {
	start_u_ =  0;
	start_v_ =  0;
	stop_u_  =  0;
	stop_v_  =  0;
	overloadCount_ = 0;
}

Element2D::Element2D(double start_u, double start_v, double stop_u, double stop_v) {

	start_u_ = start_u;
	start_v_ = start_v;
	stop_u_  = stop_u ;
	stop_v_  = stop_v ;
	overloadCount_ = 0;
}

void Element2D::removeSupportFunction(LRBSpline2D *f) {
  for (size_t i=0; i<support_.size(); i++) {
#ifndef NDEBUG
//      std::cout << "DEBUG: support_ i: " << i << std::endl;
#endif
      if((support_[i]) && (*f == *support_[i])) {
			support_[i] = support_.back();
			//support_[support_.size()-1] = NULL;
			support_.pop_back();
			return;
		}
	}
}

void Element2D::addSupportFunction(LRBSpline2D *f) {
  for (size_t i=0; i<support_.size(); i++) {
    if(f == support_[i]) {
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
	return newElement2D;
}

void Element2D::updateBasisPointers(std::vector<LRBSpline2D*> &basis) {
	for(size_t i=0; i<support_.size(); i++) {
		// add pointer from LRBSpline2D back to Element2D
		support_.back()->addSupport(this);
	}
}

void Element2D::swapParameterDirection()
{
    std::swap(start_u_, start_v_);
    std::swap(stop_u_, stop_v_);
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

/*
int Element2D::overloadedBasisCount() const {
	int ans = 0;
	for(size_t i=0; i<support_.size(); i++)
		if(support_[i]->isOverloaded())
			ans++;
	return ans;
}
*/

} // end namespace LR
