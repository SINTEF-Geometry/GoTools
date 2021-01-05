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

#include "GoTools/lrsplines3D/LRSpline3DEvalGrid.h"
#include "GoTools/lrsplines3D/LRSpline3DBezierCoefs.h"


inline int flat(int i, int j, int k, int l, int I, int J, int K)
{
  return (i + I*(j + J*(k + K*l)));
}

inline int flat(int i, int j, int k, int l, int m, int I, int J, int K, int L) { 
  return (i + I*(j + J*(k + K*(l + L*m))));
}

//==============================================================================
namespace Go
//==============================================================================
{


LRSpline3DBezierCoefs::LRSpline3DBezierCoefs() {}


LRSpline3DBezierCoefs::LRSpline3DBezierCoefs(LRSplineVolume& lr_spline)
    : dim_(lr_spline.dimension())
{
    lr_spline_ = lr_spline; 
    orig_dom_ = lr_spline.parameterSpan();
    order_u_ = 1 + lr_spline.degree(XDIR);
    order_v_ = 1 + lr_spline.degree(YDIR);
    order_w_ = 1 + lr_spline.degree(ZDIR);

    num_elements_ = lr_spline.numElements();

    // Currently only support non-rational LR-splines of degree 2 and 3 
    assert(!lr_spline.rational());
    assert(order_u_ <= 4 && order_u_ >= 3);
    assert(order_v_ <= 4 && order_v_ >= 3);
    assert(order_w_ <= 4 && order_w_ >= 3);
}

void LRSpline3DBezierCoefs::computeCoefsFromPts(const double *points, double *coefs) {
  std::vector<double> tmp_vec_x(dim_*order_u_*order_v_*order_w_);
  if (order_u_ == 3) {
    for (int jx=0; jx<order_v_; jx++) {
      for (int kx=0; kx<order_w_; kx++) {
        for (int d=0; d<dim_; d++) {
          tmp_vec_x[flat(d,0,jx,kx,dim_,order_u_,order_v_)] = points[flat(d,0,jx,kx,dim_,order_u_,order_v_)];
          tmp_vec_x[flat(d,1,jx,kx,dim_,order_u_,order_v_)] = M3_[3]*points[flat(d,0,jx,kx, dim_,order_u_,order_v_)]
                                                            + M3_[4]*points[flat(d,1,jx,kx, dim_,order_u_,order_v_)]
                                                            + M3_[5]*points[flat(d,2,jx,kx, dim_,order_u_,order_v_)];
          tmp_vec_x[flat(d,2,jx,kx,dim_,order_u_,order_v_)] = points[flat(d,2,jx,kx,dim_,order_u_,order_v_)];
	}
      } 
    }
  }
  else {
    for(int jx=0; jx<order_v_; jx++) {
      for(int kx=0; kx<order_w_; kx++) {
        for(int d=0; d<dim_; d++) {
          tmp_vec_x[flat(d,0,jx,kx,dim_,order_u_,order_v_)] =          points[flat(d,0,jx,kx,dim_,order_u_,order_v_)];
          tmp_vec_x[flat(d,1,jx,kx,dim_,order_u_,order_v_)] = M4_[ 4]*points[flat(d,0,jx,kx,dim_,order_u_,order_v_)]
                                                            + M4_[ 5]*points[flat(d,1,jx,kx,dim_,order_u_,order_v_)]
                                                            + M4_[ 6]*points[flat(d,2,jx,kx,dim_,order_u_,order_v_)]
                                                            + M4_[ 7]*points[flat(d,3,jx,kx,dim_,order_u_,order_v_)];
          tmp_vec_x[flat(d,2,jx,kx,dim_,order_u_,order_v_)] = M4_[ 8]*points[flat(d,0,jx,kx,dim_,order_u_,order_v_)]
                                                            + M4_[ 9]*points[flat(d,1,jx,kx,dim_,order_u_,order_v_)]
                                                            + M4_[10]*points[flat(d,2,jx,kx,dim_,order_u_,order_v_)]
                                                            + M4_[11]*points[flat(d,3,jx,kx,dim_,order_u_,order_v_)];
          tmp_vec_x[flat(d,3,jx,kx,dim_,order_u_,order_v_)] =          points[flat(d,3,jx,kx,dim_,order_u_,order_v_)];
        }
      }
    }
  }

  std::vector<double> tmp_vec_xy(dim_*order_u_*order_v_*order_w_);
  if (order_v_ == 3) {
    for(int ix=0; ix<order_u_; ix++) {
      for(int kx=0; kx<order_w_; kx++) {
        for(int d=0; d<dim_; d++) {
          tmp_vec_xy[flat(d,ix,0,kx, dim_,order_u_,order_v_)] =        tmp_vec_x[flat(d,ix,0,kx,dim_,order_u_,order_v_)];
          tmp_vec_xy[flat(d,ix,1,kx, dim_,order_u_,order_v_)] = M3_[3]*tmp_vec_x[flat(d,ix,0,kx,dim_,order_u_,order_v_)]
                                                              + M3_[4]*tmp_vec_x[flat(d,ix,1,kx,dim_,order_u_,order_v_)]
                                                              + M3_[5]*tmp_vec_x[flat(d,ix,2,kx,dim_,order_u_,order_v_)];
          tmp_vec_xy[flat(d,ix,2,kx, dim_,order_u_,order_v_)] =        tmp_vec_x[flat(d,ix,2,kx,dim_,order_u_,order_v_)];
        }
      }
    }
  }
  else {
    for(int ix=0; ix<order_u_; ix++) {
      for(int kx=0; kx<order_w_; kx++) {
        for(int d=0; d<dim_; d++) {
          tmp_vec_xy[flat(d,ix,0,kx, dim_,order_u_,order_v_)] =         tmp_vec_x[flat(d,ix,0,kx,dim_,order_u_,order_v_)];
          tmp_vec_xy[flat(d,ix,1,kx, dim_,order_u_,order_v_)] = M4_[ 4]*tmp_vec_x[flat(d,ix,0,kx,dim_,order_u_,order_v_)]
                                                              + M4_[ 5]*tmp_vec_x[flat(d,ix,1,kx,dim_,order_u_,order_v_)]
                                                              + M4_[ 6]*tmp_vec_x[flat(d,ix,2,kx,dim_,order_u_,order_v_)]
                                                              + M4_[ 7]*tmp_vec_x[flat(d,ix,3,kx,dim_,order_u_,order_v_)];
          tmp_vec_xy[flat(d,ix,2,kx, dim_,order_u_,order_v_)] = M4_[ 8]*tmp_vec_x[flat(d,ix,0,kx,dim_,order_u_,order_v_)]
                                                              + M4_[ 9]*tmp_vec_x[flat(d,ix,1,kx,dim_,order_u_,order_v_)]
                                                              + M4_[10]*tmp_vec_x[flat(d,ix,2,kx,dim_,order_u_,order_v_)]
                                                              + M4_[11]*tmp_vec_x[flat(d,ix,3,kx,dim_,order_u_,order_v_)];
          tmp_vec_xy[flat(d,ix,3,kx, dim_,order_u_,order_v_)] =         tmp_vec_x[flat(d,ix,3,kx,dim_,order_u_,order_v_)];
        }
      }
    }
  }

  if (order_w_ == 3) {
    for(int ix=0; ix<order_u_; ix++) {
      for(int jx=0; jx<order_v_; jx++) {
        for(int d=0; d<dim_; d++) {
          coefs[flat(d,ix,jx,0,dim_,order_u_,order_v_)] =        tmp_vec_xy[flat(d,ix,jx,0,dim_,order_u_,order_v_)];
          coefs[flat(d,ix,jx,1,dim_,order_u_,order_v_)] = M3_[3]*tmp_vec_xy[flat(d,ix,jx,0,dim_,order_u_,order_v_)]
                                                        + M3_[4]*tmp_vec_xy[flat(d,ix,jx,1,dim_,order_u_,order_v_)]
                                                        + M3_[5]*tmp_vec_xy[flat(d,ix,jx,2,dim_,order_u_,order_v_)];
          coefs[flat(d,ix,jx,2,dim_,order_u_,order_v_)] =        tmp_vec_xy[flat(d,ix,jx,2,dim_,order_u_,order_v_)];
        }
      }
    }
  }
  else {
    for(int ix=0; ix<order_u_; ix++) {
      for(int jx=0; jx<order_v_; jx++) {
        for(int d=0; d<dim_; d++) {
          coefs[flat(d,ix,jx,0,dim_,order_u_,order_v_)] =         tmp_vec_xy[flat(d,ix,jx,0,dim_,order_u_,order_v_)];
          coefs[flat(d,ix,jx,1,dim_,order_u_,order_v_)] = M4_[ 4]*tmp_vec_xy[flat(d,ix,jx,0,dim_,order_u_,order_v_)]
                                                        + M4_[ 5]*tmp_vec_xy[flat(d,ix,jx,1,dim_,order_u_,order_v_)]
                                                        + M4_[ 6]*tmp_vec_xy[flat(d,ix,jx,2,dim_,order_u_,order_v_)]
                                                        + M4_[ 7]*tmp_vec_xy[flat(d,ix,jx,3,dim_,order_u_,order_v_)];
          coefs[flat(d,ix,jx,2,dim_,order_u_,order_v_)] = M4_[ 8]*tmp_vec_xy[flat(d,ix,jx,0,dim_,order_u_,order_v_)]
                                                        + M4_[ 9]*tmp_vec_xy[flat(d,ix,jx,1,dim_,order_u_,order_v_)]
                                                        + M4_[10]*tmp_vec_xy[flat(d,ix,jx,2,dim_,order_u_,order_v_)]
                                                        + M4_[11]*tmp_vec_xy[flat(d,ix,jx,3,dim_,order_u_,order_v_)];
          coefs[flat(d,ix,jx,3,dim_,order_u_,order_v_)] =         tmp_vec_xy[flat(d,ix,jx,3,dim_,order_u_,order_v_)];
        }
      }
    }
  }
}


void LRSpline3DBezierCoefs::getBezierCoefs() {
  // Storage for points and coefs
  std::vector<double> points(dim_*order_u_*order_v_*order_w_);
  std::vector<double> coefs(dim_*order_u_*order_v_*order_w_);

  boxes_.clear();
  LRSpline3DEvalGrid eval_grid(lr_spline_);

  // Loop through elements, sample and compute Bezier coefficients from samples
  int i = 0;
  for(auto it=eval_grid.elements_begin(); it!=eval_grid.elements_end(); it++, i++) {
    
    // Print percentage
    static int printedPercentage = 0;
    int currPercentage = (int)(100 * (float)i / ((float)eval_grid.numElements() - 1));
    if (currPercentage > printedPercentage) {
      std::cout<< currPercentage <<"\% done    \r" << std::flush;
      printedPercentage = currPercentage;
    }

    // Get box
    //Go::Point low(it->umin(), it->vmin(), it->wmin()), high(it->umax(), it->vmax(), it->wmax());
    double ll_x, ll_y, ll_z, ur_x, ur_y, ur_z;
    eval_grid.low(*it, ll_x, ll_y, ll_z);
    eval_grid.high(*it, ur_x, ur_y, ur_z);
    Point low(ll_x, ll_y, ll_z);
    Point high(ur_x, ur_y, ur_z);
    boxes_.push_back(Go::BoundingBox(low, high));

    // determine domain
    double diag_len = (high-low).length();
    if (i==0) {
      min_box_diagonal_ = diag_len; 
      max_box_diagonal_ = diag_len;
    }
    else {
      min_box_diagonal_ = std::min( diag_len, min_box_diagonal_ );
      max_box_diagonal_ = std::max( diag_len, max_box_diagonal_ );
    }

    // Sample points in element
    eval_grid.evaluateGrid(*it, &points[0]);
    
    // Compute Bezier coefficients
    computeCoefsFromPts(&points[0], &coefs[0]);

    for (size_t ix = 0; ix<coefs.size(); ++ix) {
      bezier_coefs_.push_back(coefs[ix]);
    }
  }
  std::cout << std::endl; 
  calcMinMaxCoefficientsValue();
}

void LRSpline3DBezierCoefs::calcMinMaxCoefficientsValue() {
  min_coef_value_.resize(dim_+1);
  max_coef_value_.resize(dim_+1);
  int element=0;
  int ox=0;
  int oy=0;
  int oz=0;
  double len = 0;
  for (int cd=0; cd<(int)dim_; cd++) {
    const double tmp = bezier_coefs_[flat(cd,ox,oy,oz,element, (int)dim_, order_u_, order_v_, order_w_)];
    min_coef_value_[cd] = tmp;
    max_coef_value_[cd] = tmp;
    len += tmp*tmp;
  }
  min_coef_value_[(int)dim_] = std::sqrt(len);
  max_coef_value_[(int)dim_] = std::sqrt(len);
  for(int element=0; element<num_elements_; element++) {
    for(int ox=0; ox<order_u_; ox++) {
      for(int oy=0; oy<order_v_; oy++) { 
	for(int oz=0; oz<order_w_; oz++) {
	  double len = 0;
	  for (int cd=0; cd<(int)dim_; cd++) {
            const double tmp = bezier_coefs_[flat(cd,ox,oy,oz,element, (int)dim_, order_u_, order_v_, order_w_)];
	    min_coef_value_[cd] = std::min( min_coef_value_[cd], tmp );
	    max_coef_value_[cd] = std::max( max_coef_value_[cd], tmp );
	    len += tmp*tmp;
	  }
	  min_coef_value_[(int)dim_] = std::min(min_coef_value_[(int)dim_],std::sqrt(len));
	  max_coef_value_[(int)dim_] = std::max(max_coef_value_[(int)dim_],std::sqrt(len));
	}
      }
    }
  }
}


void LRSpline3DBezierCoefs::writeToFile(const std::string& filename) {
  std::ofstream ofs;
  ofs.open(filename, std::ios::out | std::ios::binary);
  writeToStream(ofs);
  ofs.close();
}

void LRSpline3DBezierCoefs::writeToStream(std::ostream& os) {
  const int dim_dom = 3;
  const int version = 3; // this might need updating later
  const int num_boxes = boxes_.size();
  float lowx = static_cast<float>(orig_dom_[0]);
  float lowy = static_cast<float>(orig_dom_[2]);
  float lowz = static_cast<float>(orig_dom_[4]);
  float highx = static_cast<float>(orig_dom_[1]);
  float highy = static_cast<float>(orig_dom_[3]);
  float highz = static_cast<float>(orig_dom_[5]);
  float min_box_diagonal = static_cast<float>(min_box_diagonal_);
  float max_box_diagonal = static_cast<float>(max_box_diagonal_);
  os.write((char*)&version, sizeof(int));
  os.write((char*)&dim_dom, sizeof(int));
  os.write((char*)&dim_, sizeof(int));
  os.write((char*)&order_u_, sizeof(int));
  os.write((char*)&order_v_, sizeof(int));
  os.write((char*)&order_w_, sizeof(int));
  os.write((char*)&num_boxes, sizeof(int));
  os.write((char*)&lowx, sizeof(float));
  os.write((char*)&lowy, sizeof(float));
  os.write((char*)&lowz, sizeof(float));
  os.write((char*)&highx, sizeof(float));
  os.write((char*)&highy, sizeof(float));
  os.write((char*)&highz, sizeof(float));
  os.write((char*)&min_box_diagonal, sizeof(float));
  os.write((char*)&max_box_diagonal, sizeof(float));
  
  std::cout << "version " << version << std::endl;
  std::cout << "dim_dom " << dim_dom << std::endl;
  std::cout << "dim_ " << dim_ << std::endl;
  std::cout << "order_u_ " << order_u_ << std::endl;
  std::cout << "order_v_ " << order_v_ << std::endl;
  std::cout << "order_w_ " << order_w_ << std::endl;
  std::cout << "num_boxes " << num_boxes << std::endl;
  std::cout << "low " << lowx << " " << lowy << " " << lowz << std::endl;
  std::cout << "high " << highx << " " << highy << " " << highz << std::endl;
  std::cout << "min_box_diagonal " << min_box_diagonal << std::endl;
  std::cout << "max_box_diagonal " << max_box_diagonal << std::endl;



  for (size_t ix = 0; ix < min_coef_value_.size(); ix++)
  {
    float tmp = static_cast<float>(min_coef_value_[ix]);
    os.write((char *)&tmp, sizeof(float));
  }
  for (size_t ix = 0; ix < max_coef_value_.size(); ix++)
  {
    float tmp = static_cast<float>(max_coef_value_[ix]);
    os.write((char *)&tmp, sizeof(float));
  }

  for (size_t ix = 0; ix < boxes_.size(); ix++)
  {
    Go::Point low = boxes_[ix].low();
    Go::Point high = boxes_[ix].high();
    lowx = static_cast<float>(low[0]);
    lowy = static_cast<float>(low[1]);
    lowz = static_cast<float>(low[2]);
    highx = static_cast<float>(high[0]);
    highy = static_cast<float>(high[1]);
    highz = static_cast<float>(high[2]);
    os.write((char *)&lowx, sizeof(float));
    os.write((char *)&lowy, sizeof(float));
    os.write((char *)&lowz, sizeof(float));
    os.write((char *)&highx, sizeof(float));
    os.write((char *)&highy, sizeof(float));
    os.write((char *)&highz, sizeof(float));
  }
  for (size_t ix = 0; ix < bezier_coefs_.size(); ix++)
  {
    float bezier_coef = static_cast<float>(bezier_coefs_[ix]);
    os.write((char *)&bezier_coef, sizeof(float));
  }

}

} // end namespace Go
