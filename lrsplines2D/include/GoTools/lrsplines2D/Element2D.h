#ifndef ELEMENT2D_H
#define ELEMENT2D_H

#include <vector>

namespace Go {

class LRBSpline2D;

class Element2D  
 {
public:
	Element2D();
	Element2D(double start_u, double start_v, double stop_u, double stop_v);
	void removeSupportFunction(LRBSpline2D *f);
	void addSupportFunction(LRBSpline2D *f);
	Element2D *split(bool split_u, double par_value);
	Element2D* copy();
	// get/set methods
	double umin() const         { return start_u_; };
	double vmin() const         { return start_v_; };
	double umax() const         { return stop_u_;  };
	double vmax() const         { return stop_v_;  };
	double area() const         { return (stop_v_-start_v_)*(stop_u_-start_u_);  };
	std::vector<LRBSpline2D*>::iterator supportBegin() { return support_.begin(); };
	std::vector<LRBSpline2D*>::iterator supportEnd()   { return support_.end();   };
	std::vector<LRBSpline2D*>::const_iterator supportBegin()const { return support_.begin(); };
	std::vector<LRBSpline2D*>::const_iterator supportEnd() const  { return support_.end();   };
	const std::vector<LRBSpline2D*>& getSupport() const
	{
	  return support_;
	}

	LRBSpline2D* supportFunction(int i) { return support_[i];   };
	int nmbBasisFunctions() const       { return (int)support_.size(); };
	void setUmin(double u)                           { start_u_ = u; };
	void setVmin(double v)                           { start_v_ = v; };
	void setUmax(double u)                           { stop_u_  = u; };
	void setVmax(double v)                           { stop_v_  = v; };

	bool isOverloaded() const;
	void resetOverloadCount()    { overloadCount_ = 0;      }
	int incrementOverloadCount() { return overloadCount_++; }
	int getOverloadCount() const { return overloadCount_;   }


	void updateBasisPointers(std::vector<LRBSpline2D*> &basis) ;

        void swapParameterDirection();

private:
	double start_u_;
	double start_v_;
	double stop_u_;
	double stop_v_;

	std::vector<LRBSpline2D*> support_;

	int overloadCount_ ;
	
};

} // end namespace Go

#endif

