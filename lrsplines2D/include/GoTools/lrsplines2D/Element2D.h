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

#ifndef ELEMENT2D_H
#define ELEMENT2D_H

#include <vector>
#include "GoTools/utils/config.h"
#include "GoTools/lrsplines2D/Direction2D.h"
#include "GoTools/geometry/SplineCurve.h"

namespace Go {

  class LRBSpline2D;

struct LSSmoothData
{
  LSSmoothData()
  {
    ncond_ = 0;
    average_error_ = 0.0;
    max_error_ = -1.0;
    nmb_outside_tol_ = 0;
  }

  bool hasDataPoints()
  {
    return (data_points_.size() > 0);
  }

  void eraseDataPoints()
  {
    data_points_.clear();
  }

  void addDataPoints(std::vector<double>::iterator start, 
		     std::vector<double>::iterator end)
  {
    data_points_.insert(data_points_.end(), start, end);
  }

  std::vector<double>& getDataPoints()
  {
    return data_points_;
  }

  void getOutsidePoints(std::vector<double>& points, int dim,
			Direction2D d, double start, double end);
  
  int dataPointSize()
  {
    return (int)data_points_.size();
  }

  bool hasLSMatrix()
  {
    return (LSmat_.size() > 0);
  }

  void setLSMatrix(int nmb, int dim)
  {
    ncond_ = nmb;
    LSmat_.assign(nmb*nmb, 0.0);
    LSright_.assign(nmb*dim, 0.0);
  }

  void getLSMatrix(double*& LSmat, double*& LSright, int& ncond)
  {
    LSmat = &LSmat_[0];
    LSright = &LSright_[0];
    ncond = ncond_;
  }

  bool hasAccuracyInfo()
  {
    return (max_error_ >=  0.0);
  }

  void getAccuracyInfo(double& average_error, double& max_error,
		       int& nmb_outside_tol)
  {
    average_error = average_error_;
    max_error = max_error_;
    nmb_outside_tol = nmb_outside_tol_;
  }

  void setAccuracyInfo(double average_error, double max_error,
		       int nmb_outside_tol)
  {
    average_error_ = average_error;
    max_error_ = max_error;
    nmb_outside_tol_ = nmb_outside_tol;
  }

  std::vector<double> data_points_;
  std::vector<double> LSmat_;
  std::vector<double> LSright_;
  int ncond_;

  double average_error_;
  double max_error_;
  int nmb_outside_tol_;
};



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

	/// Check if the element is associated data points to be used in 
	/// least squares approximation
	bool hasDataPoints()
	{
	  if (LSdata_.get())
	    return LSdata_->hasDataPoints();
	  else
	    return false;
	}

	/// Number of scattered data points
	int nmbDataPoints();

	/// Remove data points associated with the element
	void eraseDataPoints()
	{
	  if (LSdata_.get())
	    LSdata_->eraseDataPoints();
	}

	/// Add data points to the element
	void addDataPoints(std::vector<double>::iterator start, 
			   std::vector<double>::iterator end)
	{
	  if (!LSdata_)
	    LSdata_ = shared_ptr<LSSmoothData>(new LSSmoothData());
	  LSdata_->addDataPoints(start, end);
	}

	/// Fetch data points
	std::vector<double>& getDataPoints()
	  {
	    if (!LSdata_)
	      LSdata_ = shared_ptr<LSSmoothData>(new LSSmoothData());
	    return LSdata_->getDataPoints();
	  }

	/// Split point set according to a modified size of the element
	/// and return the points lying outside the current element
	void getOutsidePoints(std::vector<double>& points, Direction2D d);

	/// Check if a submatrix for least squares approximation exists
	bool hasLSMatrix()
	{
	  if (LSdata_.get())
	    return LSdata_->hasLSMatrix();
	  else
	    return false;
	}

	/// Create scratch to store local least squares approximation
	/// arrays
	void setLSMatrix();

	/// Fetch local least squares arrays
	void getLSMatrix(double*& LSmat, double*& LSright, int& ncond)
	{
	  if (!LSdata_)
	    LSdata_ = shared_ptr<LSSmoothData>(new LSSmoothData());
	  LSdata_->getLSMatrix(LSmat, LSright, ncond);
	}

	/// Check if the element has accuracy information
	bool hasAccuracyInfo()
	{
	  if (LSdata_.get())
	    return LSdata_->hasAccuracyInfo();
	  else
	    return false;
	}

	/// Fetch accuracy information
	void getAccuracyInfo(double& average_error, double& max_error,
			    int& nmb_outside_tol)
	{
	  if (LSdata_.get())
	    LSdata_->getAccuracyInfo(average_error, max_error, nmb_outside_tol);
	  else
	    {
	      average_error = max_error = 0.0;
	      nmb_outside_tol = 0;
	    }
	}

	/// Store accuracy information
	void setAccuracyInfo(double average_error, double max_error,
			    int nmb_outside_tol)
	{
	  if (!LSdata_)
	    LSdata_ = shared_ptr<LSSmoothData>(new LSSmoothData());
	  LSdata_->setAccuracyInfo(average_error, max_error, nmb_outside_tol);
	}


	/// Check if the element has been modified lately
	bool isModified()
	{
	  return is_modified_;
	}

	/// Reset modified flag
	void resetModificationFlag()
	{
	  is_modified_ = false;
	}

	/// Set modified flag to true
	void setModified()
	{
	  is_modified_ = true;
	}

	/// Get the coefficients of the underlying lr-splinesurface on this element, expressed by the
	/// Bernstein basis after a linear transformation sending this element to the unit square
	/// \return          a vector of the coefficients of the control points p_ij in order p_00[0], p_00[1], ..., p_10[0], ..., p_01[0], ...
	std::vector<double> unitSquareBernsteinBasis() const;

	/// Get the Bezier curve (as a spline curve) given as the image of a line segment in the parameter space of the spline
	/// surface. The parameter domain of the curve will be [0, 1].
	SplineCurve* curveOnElement(double start_u, double start_v, double end_u, double end_v) const;


private:
	double start_u_;
	double start_v_;
	double stop_u_;
	double stop_v_;

	std::vector<LRBSpline2D*> support_;

	int overloadCount_ ;

	bool is_modified_;

	// Information used in the context of least squares approximation
	// with smoothing
	mutable shared_ptr<LSSmoothData> LSdata_;

	// Get the evaluations of the Bernstein functions up to given degree.
	// The evaluation of the j-th Bernstein function of degree i will be
	// stored as result[i][j] where 0 <= j <= i <= degree
	void bernsteinEvaluation(int degree, double value, std::vector<vector<double> >& result) const;

	// For the linear function L, where L(0)=start and L(1)=end, we can express B^i_d(L) as
	// sum_{j=0} ^d c_{ij} B^i_d(t) where B^i_d is the i-th Bernstain basis function of degree d.
	// This method returns the c_{ij}-values, multiplied by binomial(d,j)
	void univariateBernsteinEvaluationInLine(int degree, double start, double end, std::vector<vector<double> >& result) const;
};

} // end namespace Go

#endif

