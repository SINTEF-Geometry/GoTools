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

  /// Storage for data points situated in this element, sub matrices used in least
  // squares approximation, and accuracy information derived from these points
  /// with relation to the LR B-spline surface where the element belongs.
  /// This surface is expected to approximate the point cloud where the points
  /// of this element make up a sub set
struct LSSmoothData
{
  LSSmoothData()
  {
    pt_del_ = 0;
    ncond_ = 0;
    average_error_ = accumulated_error_ = 0.0;
    max_error_ = max_error_prev_ = -1.0;
    nmb_outside_tol_ = -1;
    nmb_sign_outside_tol_ = -1;
    minheight_ = std::numeric_limits<double>::max();
    maxheight_ = std::numeric_limits<double>::lowest();
  }

  bool hasDataPoints()
  {
    return (data_points_.size() > 0);
  }

  bool hasSignificantPoints()
  {
    return (significant_points_.size() > 0);
  }

  void eraseDataPoints()
  {
    data_points_.clear();
  }

  void eraseSignificantPoints()
  {
    significant_points_.clear();
  }

  void eraseGhostPoints()
  {
    ghost_points_.clear();
  }

  void eraseDataPoints(std::vector<double>::iterator start, 
		       std::vector<double>::iterator end)
  {
    data_points_.erase(start, end);
  }

  void addDataPoints(std::vector<double>::iterator start, 
		     std::vector<double>::iterator end,
		     bool sort_in_u, int del=0)
  {
    data_points_.insert(data_points_.end(), start, end);
    sort_in_u_ = sort_in_u;
    if (pt_del_ == 0)
      pt_del_ = del;
  }

  void addDataPoints(std::vector<double>::iterator start, 
		     std::vector<double>::iterator end,
		     int del, bool sort_in_u, 
		     bool prepare_outlier_detection)
  {
    for (std::vector<double>::iterator curr=start; curr!= end; curr+=del)
      {
	data_points_.insert(data_points_.end(), curr, curr+del);
	data_points_.push_back(0.0);
	if (prepare_outlier_detection)
	  data_points_.push_back(1.0);
      }
    sort_in_u_ = sort_in_u;
    if (pt_del_ == 0)
      pt_del_ = del+1+(prepare_outlier_detection);
  }

  void addSignificantPoints(std::vector<double>::iterator start, 
			    std::vector<double>::iterator end,
			    bool sort_in_u, int del=0)
  {
    significant_points_.insert(significant_points_.end(), start, end);
    sort_in_u_ = sort_in_u;
    if (pt_del_ == 0)
      pt_del_ = del;
  }

  void addSignificantPoints(std::vector<double>::iterator start, 
			    std::vector<double>::iterator end,
			    int del, bool sort_in_u, 
			    bool prepare_outlier_detection)
  {
    for (std::vector<double>::iterator curr=start; curr!= end; curr+=del)
      {
	significant_points_.insert(significant_points_.end(), curr, curr+del);
	significant_points_.push_back(0.0);
	if (prepare_outlier_detection)
	  significant_points_.push_back(1.0);
      }
    sort_in_u_ = sort_in_u;
    if (pt_del_ == 0)
      pt_del_ = del+1+(prepare_outlier_detection);
  }

  void addGhostPoints(std::vector<double>::iterator start, 
		      std::vector<double>::iterator end,
		      bool sort_in_u, int del=0)
  {
    ghost_points_.insert(ghost_points_.end(), start, end);
    sort_in_u_ = sort_in_u;
    if (pt_del_ == 0)
      pt_del_ = del;
  }

  void addGhostPoints(std::vector<double>::iterator start, 
		      std::vector<double>::iterator end,
		      int del, bool sort_in_u, 
		      bool prepare_outlier_detection)
  {
    for (std::vector<double>::iterator curr=start; curr!= end; curr+=del)
      {
	ghost_points_.insert(ghost_points_.end(), curr, curr+del);
	ghost_points_.push_back(0.0);
	if (prepare_outlier_detection)
	  ghost_points_.push_back(1.0);
      }
    sort_in_u_ghost_ = sort_in_u;
    if (pt_del_ == 0)
      pt_del_ = del+1+(prepare_outlier_detection);
  }

  int getNmbValPrPoint()
  {
    return pt_del_;
  }

   std::vector<double>& getDataPoints()
  {
   return data_points_;
  }

   std::vector<double>& getSignificantPoints()
  {
   return significant_points_;
  }

  std::vector<double>& getGhostPoints()
  {
   return ghost_points_;
  }

  void getOutsidePoints(std::vector<double>& points, int dim,
			Direction2D d, double start, double end,
			bool& sort_in_u);
  
  void getOutsideSignificantPoints(std::vector<double>& points, int dim,
				   Direction2D d, double start, double end,
				   bool& sort_in_u);
  
  void getOutsideGhostPoints(std::vector<double>& ghost, int dim,
			     Direction2D d, double start, double end,
			     bool& sort_in_u);
  
  int dataPointSize()
  {
    return (int)data_points_.size();
  }

  int significantPointSize()
  {
    return (int)significant_points_.size();
  }

  int ghostPointSize()
  {
    return (int)ghost_points_.size();
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
		       int& nmb_outside_tol, int& nmb_out_sign)
  {
    average_error = average_error_;
    max_error = max_error_;
    nmb_outside_tol = nmb_outside_tol_;
    nmb_out_sign = nmb_sign_outside_tol_;
  }

  int getNmbOutsideTol()
  {
    return nmb_outside_tol_;
  }

  int getNmbSignOutsideTol()
  {
    return nmb_outside_tol_;
  }

  void getInfoSignificantPoints(int dim, double& maxdist, double& avdist,
				int& nmb_pts)
  {
    int ix = 2 + dim;  // Distance information
    nmb_pts = (int)significant_points_.size()/pt_del_;
    maxdist = 0.0;
    avdist = 0.0;
    for (std::vector<double>::iterator it=significant_points_.begin(); 
	 it != significant_points_.end(); it+=pt_del_)
      {
	double d2 = fabs(*(it+ix));
	maxdist = std::max(maxdist, d2);
	avdist += d2;
      }
    if (nmb_pts > 0)
      avdist /= (double)nmb_pts;
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
		       int nmb_outside_tol, int nmb_outside_sign,
		       double accumulated_out)
		       
  {
    accumulated_error_ = accumulated_error;
    average_error_ = average_error;
    max_error_prev_ = max_error_;
    max_error_ = max_error;
    nmb_outside_tol_ = nmb_outside_tol;
    nmb_sign_outside_tol_ = nmb_outside_sign;
    accumulated_out_ = accumulated_out;
  }

  void resetAccuracyInfo()
  {
    accumulated_error_ = 0.0;
    average_error_ = 0.0;
    max_error_ = max_error_prev_ = -1.0;
    nmb_outside_tol_ = -1;
    accumulated_out_ = 0.0;
    nmb_sign_outside_tol_ = -1;
  }

  void setHeightInfo(double minheight, double maxheight)
  {
    minheight_ = minheight;
    maxheight_ = maxheight;
  }

  void resetHeightInfo()
  {
    minheight_ = std::numeric_limits<double>::max();
    maxheight_ = std::numeric_limits<double>::lowest();
  }

  int getNmbOutliers()
  {
    if (pt_del_ != 5)
      return 0;  // No outlier information
    int nmb = 0;
    for (std::vector<double>::iterator it=data_points_.begin(); 
	 it != data_points_.end(); it+=pt_del_)
      {
	if (*(it+4) < 0.0)
	  ++nmb;
      }
    return nmb;
  }

  void getOutlierPts(std::vector<double>& outliers)
  {
    if (pt_del_ != 5)
      return;  // No outlier information
    int ix1 = pt_del_ - 5;
    int ix2 = pt_del_ - 2;
    for (std::vector<double>::iterator it=data_points_.begin(); 
	 it != data_points_.end(); it+=pt_del_)
      {
	if (*(it+4) < 0.0)
	  outliers.insert(outliers.end(), it+ix1, it+ix2);
      }
  }

  void getRegularPts(std::vector<double>& regular)
  {
    if (pt_del_ != 5)
      return;  // No outlier information
    int ix1 = pt_del_ - 5;
    int ix2 = pt_del_ - 2;
    for (std::vector<double>::iterator it=data_points_.begin(); 
	 it != data_points_.end(); it+=pt_del_)
      {
	if (*(it+4) > 0.0)
	  regular.insert(regular.end(), it+ix1, it+ix2);
      }
    for (std::vector<double>::iterator it=significant_points_.begin(); 
	 it != significant_points_.end(); it+=pt_del_)
      {
	if (*(it+4) > 0.0)
	  regular.insert(regular.end(), it+ix1, it+ix2);
      }
  }

  void getClassifiedPts(std::vector<double>& outliers,
			std::vector<double>& regular)
  {
    if (pt_del_ != 5)
      return;  // No outlier information
    int ix1 = pt_del_ - 5;
    int ix2 = pt_del_ - 2;
    for (std::vector<double>::iterator it=data_points_.begin(); 
	 it != data_points_.end(); it+=pt_del_)
      {
	if (*(it+4) < 0.0)
	  outliers.insert(outliers.end(), it+ix1, it+ix2);
	else if (*(it+4) > 0.0)
	  regular.insert(regular.end(), it+ix1, it+ix2);
      }
    for (std::vector<double>::iterator it=significant_points_.begin(); 
	 it != significant_points_.end(); it+=pt_del_)
      {
	if (*(it+4) < 0.0)
	  outliers.insert(outliers.end(), it+ix1, it+ix2);
	else if (*(it+4) > 0.0)
	  regular.insert(regular.end(), it+ix1, it+ix2);
      }
  }

  bool getDataBoundingBox(int dim, double bb[]);

  void makeDataPoints3D(int dim);

  void updateAccuracyInfo(int dim);

  void updateLSDataParDomain(double u1, double u2, 
			     double v1, double v2, 
			     double u1new, double u2new, 
			     double v1new, double v2new,
			     int dim);

  std::vector<double> data_points_;
  int pt_del_;
  std::vector<double> significant_points_;
  std::vector<double> ghost_points_;
  std::vector<double> LSmat_;
  std::vector<double> LSright_;
  int ncond_;
  bool sort_in_u_;
  bool sort_in_u_ghost_;

  double accumulated_error_;
  double average_error_;
  double max_error_;
  double max_error_prev_;
  int nmb_outside_tol_;
  int nmb_sign_outside_tol_;
  double accumulated_out_;
  double minheight_;
  double maxheight_;
};


  /// An element (or a mesh cell) in an LR B-spline surface description. The
  /// elements contain information about the B-splines having this element
  /// in its support and has the potential of storing data points corresponding
  /// to this element with derived accuracy information.
class Element2D  
 {
public:
   /// Constructor for an empty element
   Element2D();
   
   /// Constructor giving element boundaries
   Element2D(double start_u, double start_v, double stop_u, double stop_v);
   
   /// Destructor
   ~Element2D();
   
   /// Remove a B-spline with this element in its support from the vector maintained in the element
   void removeSupportFunction(LRBSpline2D *f);
   
   /// Add a B-spline with this element in its support from the vector maintained in the element
   void addSupportFunction(LRBSpline2D *f);
   
   /// Check if a B-spline has this element in its support
   bool hasSupportFunction(LRBSpline2D *f);
   
   /// Split element in the given direction and value. Called from surface refinement.
   Element2D *split(bool split_u, double par_value);

   Element2D* copy();
   
   // get/set methods
   /// Start parameter in the first parameter direction
   double umin() const         { return start_u_; };
   /// Start parameter in the second parameter direction
   double vmin() const         { return start_v_; };
   /// End parameter in the first parameter direction
	double umax() const         { return stop_u_;  };
   /// End parameter in the second parameter direction
   double vmax() const         { return stop_v_;  };
   /// Area of element domain
   double area() const         { return (stop_v_-start_v_)*(stop_u_-start_u_);  };

   /// Start iterator to B-splines with this element in their support
   std::vector<LRBSpline2D*>::iterator supportBegin() { return support_.begin(); };
   
   /// End iterator to B-splines with this element in their support
   std::vector<LRBSpline2D*>::iterator supportEnd()   { return support_.end();   };
   /// Start iterator to B-splines with this element in their support
   std::vector<LRBSpline2D*>::const_iterator supportBegin()const { return support_.begin(); };
   /// End iterator to B-splines with this element in their support
   std::vector<LRBSpline2D*>::const_iterator supportEnd() const  { return support_.end();   };

   /// Return a reference to the vector of B-splines with this element in their support
   const std::vector<LRBSpline2D*>& getSupport() const
   {
     return support_;
   }

   /// Number of B-splines with this element in their support
   int nmbSupport() const
   {
     return (int)support_.size();
   }

   /* std::vector<LRBSpline2D*> getSupport()  */
   /* { */
   /*   return support_; */
   /* } */

   /// Check if the parameter pair is contained in the element domain
   bool contains(double upar, double vpar)
   {
     return (upar >= start_u_ && upar <= stop_u_ && 
	     vpar >= start_v_ && vpar <= stop_v_);
   }

   /// Accsess one specified B-spline with this element in the
   LRBSpline2D* supportFunction(int i) { return support_[i];   };
   /// Number of B-splines having this element in their support
   int nmbBasisFunctions() const       { return (int)support_.size(); };
   /// Modify the start of the element domain in the first parameter direction
   void setUmin(double u)                           { start_u_ = u; };
   /// Modify the start of the element domain in the second parameter direction
   void setVmin(double v)                           { start_v_ = v; };
   /// Modify the end of the element domain in the first parameter direction
   void setUmax(double u)                           { stop_u_  = u; };
   /// Modify the end of the element domain in the second parameter direction
   void setVmax(double v)                           { stop_v_  = v; };

   bool isOverloaded() const;
   void resetOverloadCount()    { overloadCount_ = 0;      }
   int incrementOverloadCount() { return overloadCount_++; }
   int getOverloadCount() const { return overloadCount_;   }


   void updateBasisPointers(std::vector<LRBSpline2D*> &basis) ;

   void swapParameterDirection();

   /// Fetch neighbouring elements based on information on
   /// supporting B-splines. Corner touch excluded
   void fetchNeighbours(std::vector<Element2D*>& neighbours) const;

   /// Check if the element is associated data points to be used in 
   /// least squares approximation
   bool hasDataPoints()
   {
     if (LSdata_.get())
       return LSdata_->hasDataPoints();
     else
       return false;
   }

   bool hasSignificantPoints()
   {
     if (LSdata_.get())
       return LSdata_->hasSignificantPoints();
     else
       return false;
   }

   /// Number of double values for each point
   int getNmbValPrPoint()
   {
     if (LSdata_.get())
       return LSdata_->getNmbValPrPoint();
     else
       return 0;
   }
	  
   /// Number of scattered data points
   int nmbDataPoints();

   /// Number of significant data points
   int nmbSignificantPoints();

   /// Number of ghost points
   int nmbGhostPoints();

   /// Remove data points associated with the element
   void eraseDataPoints()
   {
     if (LSdata_.get())
       LSdata_->eraseDataPoints();
   }

   void eraseSignificantPoints()
   {
     if (LSdata_.get())
       LSdata_->eraseSignificantPoints();
   }

   void eraseDataPoints(std::vector<double>::iterator start, 
			std::vector<double>::iterator end)
   {
     if (LSdata_.get())
       LSdata_->eraseDataPoints(start, end);
   }

   void eraseGhostPoints()
   {
     if (LSdata_.get())
       LSdata_->eraseGhostPoints();
   }

   /// Add data points to the element
   void addDataPoints(std::vector<double>::iterator start, 
		      std::vector<double>::iterator end,
		      bool sort_in_u, int del=0)
   {
     if (!LSdata_)
       LSdata_ = shared_ptr<LSSmoothData>(new LSSmoothData());
     LSdata_->addDataPoints(start, end, sort_in_u, del);
   }

   void addDataPoints(std::vector<double>::iterator start, 
		      std::vector<double>::iterator end,
		      int del, bool sort_in_u, 
		      bool prepare_outlier_detection=false)
   {
     if (!LSdata_)
       LSdata_ = shared_ptr<LSSmoothData>(new LSSmoothData());
     LSdata_->addDataPoints(start, end, del, sort_in_u,
			    prepare_outlier_detection);
   }

   /// Add significante data points to the element
   void addSignificantPoints(std::vector<double>::iterator start, 
			     std::vector<double>::iterator end,
			     bool sort_in_u, int del=0)
   {
     if (!LSdata_)
       LSdata_ = shared_ptr<LSSmoothData>(new LSSmoothData());
     LSdata_->addSignificantPoints(start, end, sort_in_u, del);
   }

   void addSignificantPoints(std::vector<double>::iterator start, 
			     std::vector<double>::iterator end,
			     int del, bool sort_in_u, 
			     bool prepare_outlier_detection=false)
   {
     if (!LSdata_)
       LSdata_ = shared_ptr<LSSmoothData>(new LSSmoothData());
     LSdata_->addSignificantPoints(start, end, del, sort_in_u,
				   prepare_outlier_detection);
   }


   void addGhostPoints(std::vector<double>::iterator start, 
		       std::vector<double>::iterator end,
		       bool sort_in_u, int del=0)
   {
     if (!LSdata_)
       LSdata_ = shared_ptr<LSSmoothData>(new LSSmoothData());
     LSdata_->addGhostPoints(start, end, sort_in_u, del);
   }
   void addGhostPoints(std::vector<double>::iterator start, 
		       std::vector<double>::iterator end,
		       int del, bool sort_in_u,
		       bool prepare_outlier_detection)
   {
     if (!LSdata_)
       LSdata_ = shared_ptr<LSSmoothData>(new LSSmoothData());
     LSdata_->addGhostPoints(start, end, del, sort_in_u,
			     prepare_outlier_detection);
   }

   /// Fetch data points
   std::vector<double>& getDataPoints()
   {
     if (!LSdata_)
       LSdata_ = shared_ptr<LSSmoothData>(new LSSmoothData());
     return LSdata_->getDataPoints();
   }

   /// Fetch significant data points
   std::vector<double>& getSignificantPoints()
   {
     if (!LSdata_)
       LSdata_ = shared_ptr<LSSmoothData>(new LSSmoothData());
     return LSdata_->getSignificantPoints();
   }

   /// Fetch artificial data points intended for stabilization
   std::vector<double>& getGhostPoints()
   {
     if (!LSdata_)
       LSdata_ = shared_ptr<LSSmoothData>(new LSSmoothData());
     return LSdata_->getGhostPoints();
   }

   /// Split point set according to a modified size of the element
   /// and return the points lying outside the current element
   void getOutsidePoints(std::vector<double>& points, Direction2D d,
			 bool& sort_in_u);

   void getOutsideSignificantPoints(std::vector<double>& points, 
				    Direction2D d, bool& sort_in_u);

   void getOutsideGhostPoints(std::vector<double>& points, Direction2D d,
			      bool& sort_in_u);

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
			int& nmb_outside_tol, int& nmb_outside_sign)
   {
     if (LSdata_.get())
       LSdata_->getAccuracyInfo(average_error, max_error, 
				nmb_outside_tol, nmb_outside_sign);
     else
       {
	 average_error = max_error = 0.0;
	 nmb_outside_tol = nmb_outside_sign = 0;
       }
   }

   int getNmbOutsideTol()
   {
     if (LSdata_.get())
       return LSdata_->getNmbOutsideTol();
     else
       return 0;
   }
	  
   int getNmbSignOutsideTol()
   {
     if (LSdata_.get())
       return LSdata_->getNmbSignOutsideTol();
     else
       return 0;
   }
	  
   void getInfoSignificantPoints(int dim, double& maxdist, 
				 double& avdist, int& nmb_pts)
   {
     if (LSdata_.get())
       LSdata_->getInfoSignificantPoints(dim, maxdist, avdist, nmb_pts);
     else
       {
	 maxdist = avdist = 0.0;
	 nmb_pts = 0;
       }
   }

   double getAverageError()
   {
     if (!LSdata_.get())
       return 0.0;
     else
       return 
	 LSdata_->getAverageError();
   }

   double getAccumulatedError()
   {
     if (!LSdata_.get())
       return 0.0;
     else
       return 
	 LSdata_->getAccumulatedError();
   }


   double getAccumulatedOutside()
   {
     if (!LSdata_.get())
       return 0.0;
     else
       return 
	 LSdata_->getAccumulatedOutside();
   }

   double getMaxError()
   {
     if (!LSdata_.get())
       return 0.0;
     else
       return 
	 LSdata_->getMaxError();
   }

   /// Store accuracy information
   void setAccuracyInfo(double accumulated_error,
			double average_error, double max_error,
			int nmb_outside_tol, int nmb_outside_sign,
			double accumulated_out=0.0)
   {
     if (!LSdata_)
       LSdata_ = shared_ptr<LSSmoothData>(new LSSmoothData());
     LSdata_->setAccuracyInfo(accumulated_error, average_error, 
			      max_error, nmb_outside_tol, 
			      nmb_outside_sign, accumulated_out);
   }


   void resetAccuracyInfo()
   {
     if (LSdata_.get())
       LSdata_->resetAccuracyInfo();
   }

   void setHeightInfo(double minheight, double maxheight)
   {
     if (LSdata_.get())
       LSdata_->setHeightInfo(minheight, maxheight);
   }

   void resetHeightInfo()
   {
     if (LSdata_.get())
       LSdata_->resetHeightInfo();
   }

   int getNmbOutliers()
   {
     if (LSdata_.get())
       return LSdata_->getNmbOutliers();
     else
       return 0;
   }

   void getOutlierPts(std::vector<double>& outliers)
   {
     if (LSdata_.get())
       LSdata_->getOutlierPts(outliers);
   }

   void getRegularPts(std::vector<double>& regular)
   {
     if (LSdata_.get())
       LSdata_->getRegularPts(regular);
   }

   void getClassifiedPts(std::vector<double>& outliers,
			 std::vector<double>& regular)
   {
     if (LSdata_.get())
       LSdata_->getClassifiedPts(outliers, regular);
   }

   /// Get box bounding the data points: min, max for each coordinate
   bool getDataBoundingBox(double bb[]);

   /// Turn data points into 3D using the parameter values
   /// as the x- and y-coordinates
   void makeDataPoints3D();

   /// Update accuracy statistics in points. Number of outside
   /// points is NOT changed
   void updateAccuracyInfo();

   /// Reparameterize LS data informaiton
   void updateLSDataParDomain(double u1, double u2, 
			      double v1, double v2, 
			      double u1new, double u2new, 
			      double v1new, double v2new);

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
   /// \param start_u The u-value of the start point of the line segment
   /// \param start_v The v-value of the start point of the line segment
   /// \param end_u   The u-value of the start point of the line segment
   /// \param end_v   The v-value of the start point of the line segment
   /// \return    The spline curve
   SplineCurve* curveOnElement(double start_u, double start_v, double end_u, double end_v) const;

   // DEBUG
   double sumOfScaledBsplines(double upar, double vpar);

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
	void bernsteinEvaluation(int degree, double value, std::vector<std::vector<double> >& result) const;

	// For the linear function L, where L(0)=start and L(1)=end, we can express B^i_d(L) as
	// sum_{j=0} ^d c_{ij} B^i_d(t) where B^i_d is the i-th Bernstein basis function of degree d.
	// This method returns the c_{ij}-values, multiplied by binomial(d,j)
	void univariateBernsteinEvaluationInLine(int degree, double start, double end, std::vector<std::vector<double> >& result) const;
};

} // end namespace Go

#endif

