//===========================================================================
//                                                                           
// File: Element3D.h                                                         
//                                                                           
// Created: Mon Feb 25 11:07:36 2013                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _ELEMENT3D_H
#define _ELEMENT3D_H

#include "GoTools/lrsplines3D/Direction3D.h"
#include "GoTools/utils/config.h"
#include <vector>

namespace Go {

class LRBSpline3D;


  /// Storage for data points situated in this element,  and accuracy information derived from 
  /// these points with relation to the LR B-spline volume where the element belongs.
  /// This volume is expected to approximate the point cloud where the points
  /// of this element make up a sub set
struct Approx3DData
{
  Approx3DData()
  {
    ncond_ = 0;
    average_error_ = 0.0;
    max_error_ = max_error_prev_ = -1.0;
    nmb_outside_tol_ = -1;
    acc_err_pos_ = acc_err_neg_ = 0.0;
    nmb_err_pos_ = nmb_err_neg_ = 0;
  }

  bool hasDataPoints()
  {
    return (data_points_.size() > 0);
  }

  void eraseDataPoints()
  {
    data_points_.clear();
  }

  void eraseDataPoints(std::vector<double>::iterator start, 
		       std::vector<double>::iterator end)
  {
    data_points_.erase(start, end);
  }

  void addDataPoints(std::vector<double>::iterator start, 
		     std::vector<double>::iterator end,
		     bool sort_in_u)
  {
    data_points_.insert(data_points_.end(), start, end);
    sort_in_u_ = sort_in_u;
  }

  void addDataPoints(std::vector<double>::iterator start, 
		     std::vector<double>::iterator end,
		     int del, bool sort_in_u)
  {
    for (std::vector<double>::iterator curr=start; curr!= end; curr+=del)
      {
	data_points_.insert(data_points_.end(), curr, curr+del);
	data_points_.push_back(0.0);
      }
    sort_in_u_ = sort_in_u;
  }

   std::vector<double>& getDataPoints()
  {
   return data_points_;
  }
  
  int dataPointSize()
  {
    return (int)data_points_.size();
  }

  void getOutsidePoints(std::vector<double>& points, int dim,
                        Direction3D d, double start, double end,
                        bool& sort_in_u);

  void updateApproxDataParDomain(double u1, double u2,
                                 double v1, double v2,
                                 double w1, double w2,
                                 double u1new, double u2new,
                                 double v1new, double v2new,
                                 double w1new, double w2new,
                                 int dim);

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

  void getSignedAccInfo(double& acc_err_pos, int& nmb_err_pos,
			double& acc_err_neg, int& nmb_err_neg)
  {
    acc_err_pos = acc_err_pos_;
    nmb_err_pos = nmb_err_pos_;
    acc_err_neg = acc_err_neg_;
    nmb_err_neg = nmb_err_neg_;
  }
  
  int getNmbOutsideTol()
  {
    return nmb_outside_tol_;
  }

  double getAverageError()
  {
    return average_error_;
  }

  double getAccumulatedError()
  {
    return accumulated_error_;
  }

  double getAccumulatedOutside()
  {
    return accumulated_out_;
  }

  double getMaxError()
  {
    return max_error_;
  }

  void setAccuracyInfo(double accumulated_error,
		       double average_error, double max_error,
		       int nmb_outside_tol, double accumulated_out)
  {
    accumulated_error_ = accumulated_error;
    average_error_ = average_error;
    max_error_prev_ = max_error_;
    max_error_ = max_error;
    nmb_outside_tol_ = nmb_outside_tol;
    accumulated_out_ = accumulated_out;
 }

  void setSignedAccInfo(double acc_err_pos, int nmb_err_pos,
			double acc_err_neg, int nmb_err_neg)
  {
    acc_err_pos_ = acc_err_pos;
    nmb_err_pos_ = nmb_err_pos;
    acc_err_neg_ = acc_err_neg;
    nmb_err_neg_ = nmb_err_neg;
  }
  
  void resetAccuracyInfo()
  {
    accumulated_error_ = 0.0;
    average_error_ = 0.0;
    max_error_ = max_error_prev_ = -1.0;
    nmb_outside_tol_ = -1;
    acc_err_pos_ = acc_err_neg_= 0.0;
    nmb_err_pos_ = nmb_err_neg_ = 0;
  }

  bool getDataBoundingBox(int dim, double bb[]);

  void makeDataPoints3D(int dim);

  void updateAccuracyInfo(int dim);

  std::vector<double> data_points_;

  int ncond_;
  bool sort_in_u_;

  double accumulated_error_;
  double average_error_;
  double max_error_;
  double max_error_prev_;
  int nmb_outside_tol_;
  double accumulated_out_;
  double acc_err_pos_;
  double acc_err_neg_;
  int nmb_err_pos_;
  int nmb_err_neg_;
};


  /// An element (or a mesh cell) in an LR B-spline volume description. The
  /// elements contain information about the LRBSpline3D having this element
  /// in its support and has the potential of storing data points corresponding
  /// to this element with derived accuracy information.


class Element3D  
{
 public:
   /// Constructor for an empty element
  Element3D();
   /// Constructor giving element boundaries
  Element3D(double start_u, double start_v, double start_w, double stop_u, double stop_v, double stop_w);
   /// Remove a B-spline with this element in its support from the vector maintained in the element
  void removeSupportFunction(LRBSpline3D *f);
   /// Add a B-spline with this element in its support from the vector maintained in the elemen
  void addSupportFunction(LRBSpline3D *f);
   /// Split element in the given direction and value. Called from volume refinement.
  Element3D *split(Direction3D split_dir, double par_value);
  Element3D* copy();
  // get/set methods
   /// Start parameter in the first parameter direction
  double umin() const         { return start_u_; }
   /// Start parameter in the second parameter direction
  double vmin() const         { return start_v_; }
   /// Start parameter in the third parameter direction
  double wmin() const         { return start_w_; }
   /// End parameter in the first parameter direction
  double umax() const         { return stop_u_;  }
   /// End parameter in the second parameter direction
  double vmax() const         { return stop_v_;  }
   /// End parameter in the third parameter direction
  double wmax() const         { return stop_w_;  }
  /// Volume of element domain
  double volume() const         { return (stop_w_-start_w_)*(stop_v_-start_v_)*(stop_u_-start_u_);  }
   /// Start iterator to B-splines with this element in their support
  std::vector<LRBSpline3D*>::iterator supportBegin() { return support_.begin(); }
   /// End iterator to B-splines with this element in their support
  std::vector<LRBSpline3D*>::iterator supportEnd()   { return support_.end();   }
   /// Start iterator to B-splines with this element in their support
  std::vector<LRBSpline3D*>::const_iterator supportBegin()const { return support_.begin(); }
   /// End iterator to B-splines with this element in their support
  std::vector<LRBSpline3D*>::const_iterator supportEnd() const  { return support_.end();   }
   /// Return a reference to the vector of B-splines with this element in their support
  const std::vector<LRBSpline3D*>& getSupport() const
  {
    return support_;
  }

   /// Check if the parameter pair is contained in the element domain
  bool contains(double upar, double vpar, double wpar)
  {
    return (upar >= start_u_ && upar <= stop_u_ &&
            vpar >= start_v_ && vpar <= stop_v_ &&
            wpar >= start_w_ && wpar <= stop_w_);
  }


   /// Accsess one specified B-spline with this element in the support
  LRBSpline3D* supportFunction(int i) { return support_[i];   }
   /// Number of B-splines having this element in their support
  int nmbBasisFunctions() const       { return (int)support_.size(); }
   /// Modify the start of the element domain in the first parameter direction
  void setUmin(double u)                           { start_u_ = u; }
   /// Modify the start of the element domain in the second parameter direction
  void setVmin(double v)                           { start_v_ = v; }
   /// Modify the start of the element domain in the third parameter direction
  void setWmin(double w)                           { start_w_ = w; }
   /// Modify the end of the element domain in the first parameter direction
  void setUmax(double u)                           { stop_u_  = u; }
   /// Modify the end of the element domain in the second parameter direction
  void setVmax(double v)                           { stop_v_  = v; }
   /// Modify the end of the element domain in the third parameter direction
  void setWmax(double w)                           { stop_w_  = w; }

  bool isOverloaded() const;
  void resetOverloadCount()    { overloadCount_ = 0;      }
  int incrementOverloadCount() { return overloadCount_++; }
  int getOverloadCount() const { return overloadCount_;   }


  void updateBasisPointers(std::vector<LRBSpline3D*> &basis) ;

  /* // @@sbr201302 This must be extended to 3D case. Swap two dirs or rotate? */
  /* void swapParameterDirection(); */

  /// Check if the element is associated data points to be used in 
  /// least squares approximation
  bool hasDataPoints()
  {
    if (approx_data_.get())
      return approx_data_->hasDataPoints();
    else
      return false;
  }

  /// Number of scattered data points
  int nmbDataPoints();
 
  /// Remove data points associated with the element
  void eraseDataPoints()
  {
    if (approx_data_.get())
      approx_data_->eraseDataPoints();
  }
  
  void eraseDataPoints(std::vector<double>::iterator start, 
                       std::vector<double>::iterator end)
  {
    if (approx_data_.get())
      approx_data_->eraseDataPoints(start, end);
  }
  
  
  /// Add data points to the element
  void addDataPoints(std::vector<double>::iterator start, 
                     std::vector<double>::iterator end,
                     bool sort_in_u)
  {
    if (!approx_data_) approx_data_ = shared_ptr<Approx3DData>(new Approx3DData());
    approx_data_->addDataPoints(start, end, sort_in_u);
  }

  void addDataPoints(std::vector<double>::iterator start, 
                     std::vector<double>::iterator end,
                     int del, bool sort_in_u)
  {
    if (!approx_data_) approx_data_ = shared_ptr<Approx3DData>(new Approx3DData());
    approx_data_->addDataPoints(start, end, del, sort_in_u);
  }
  
  /// Fetch data points
  std::vector<double>& getDataPoints() 
  {
    if (!approx_data_)
      approx_data_ = shared_ptr<Approx3DData>(new Approx3DData());
    return approx_data_->getDataPoints();
  }

  /// Split point set according to a modified size of the element
  /// and return the points lying outside the current element
  void getOutsidePoints(std::vector<double>& points, Direction3D d,
                        bool& sort_in_u);

  void updateApproxDataParDomain(double u1, double u2,
                                 double v1, double v2,
                                 double w1, double w2,
                                 double u1new, double u2new,
                                 double v1new, double v2new,
                                 double w1new, double w2new);

  /// Check if the element has accuracy information
  bool hasAccuracyInfo()
  {
    if (approx_data_.get())
      return approx_data_->hasAccuracyInfo();
    else
      return false;
  }

  /// Fetch accuracy information
  void getAccuracyInfo(double& average_error, double& max_error,
                      int& nmb_outside_tol)
  {
    if (approx_data_.get())
      approx_data_->getAccuracyInfo(average_error, max_error, nmb_outside_tol);
    else
      {
        average_error = max_error = 0.0;
        nmb_outside_tol = 0;
      }
  }

  void getSignedAccInfo(double& acc_err_pos, int& nmb_err_pos,
			double& acc_err_neg, int& nmb_err_neg)
  {
    if (approx_data_.get())
      approx_data_->getSignedAccInfo(acc_err_pos, nmb_err_pos,
				     acc_err_neg, nmb_err_neg);
    else
      {
	acc_err_pos = acc_err_neg = 0.0;
	nmb_err_pos = nmb_err_neg = 0;
      }
  }
  
  int getNmbOutsideTol()
  {
    if (approx_data_.get())
      return approx_data_->getNmbOutsideTol();
    else
      return 0;
  }
	  

  double getAverageError()
  {
    if (!approx_data_.get())
      return 0.0;
    else
      return 
        approx_data_->getAverageError();
  }

  double getAccumulatedError()
  {
    if (!approx_data_.get())
      return 0.0;
    else
      return 
        approx_data_->getAccumulatedError();
  }

  double getAccumulatedOutside()
  {
    if (!approx_data_.get())
      return 0.0;
    else
      return 
	approx_data_->getAccumulatedOutside();
  }
  
 double getMaxError()
  {
    if (!approx_data_.get())
      return 0.0;
    else
      return 
        approx_data_->getMaxError();
  }

  /// Store accuracy information
  void setAccuracyInfo(double accumulated_error,
                       double average_error, double max_error,
                       int nmb_outside_tol, double accumulated_out=0.0)
  {
    if (!approx_data_)
      approx_data_ = shared_ptr<Approx3DData>(new Approx3DData());
    approx_data_->setAccuracyInfo(accumulated_error, average_error,
				  max_error, nmb_outside_tol, accumulated_out);
  }


  void setSignedAccInfo(double acc_err_pos, int nmb_err_pos,
			double acc_err_neg, int nmb_err_neg)
  {
    if (!approx_data_)
      approx_data_ = shared_ptr<Approx3DData>(new Approx3DData());
    approx_data_->setSignedAccInfo(acc_err_pos, nmb_err_pos,
				   acc_err_neg, nmb_err_neg);
  }
  
  void resetAccuracyInfo()
  {
    if (approx_data_.get())
      approx_data_->resetAccuracyInfo();
  }

  /// Update accuracy statistics in points. 
  /// Number of outside points is NOT changed
  void updateAccuracyInfo();

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

 private:
  double start_u_;
  double start_v_;
  double start_w_;
  double stop_u_;
  double stop_v_;
  double stop_w_;

  std::vector<LRBSpline3D*> support_;

  int overloadCount_ ;

  bool is_modified_;

  // Information used in the context of approximation
  mutable shared_ptr<Approx3DData> approx_data_;

};

} // end namespace Go

#endif // _ELEMENT3D_H

