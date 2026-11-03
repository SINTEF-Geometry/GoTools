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

#ifndef _BASISGRIDRESULT_H
#define  _BASISGRIDRESULT_H

#include "GoTools/utils/config.h"

namespace Go
{
  /// Structure for storage of results of grid evaluation of the basis function of a spline surface.
  /// Positional evaluation information in one parameter value
  struct BasisPtsSf
  {
    /// Parameter tripple in which the basis functions are evaluated
    double param[2]; 
    /// Index of the knot interval where the parameter value is situated for all
    /// parameter directions. The indices of the non-zero basis functions are
    /// left_idx[i]-order[i]+1, ..., left_idx[i] for i=0,1
    int left_idx[2]; 
    /// The value of all basis functions, size equal to (degree_u+1)*(degree_v+1)*(degree_w+1)
    std::vector< double > basisValues; 

    void preparePts(double u, double v, int idx_u, int idx_v, int size)
    {
      param[0] = u;
      param[1] = v;
      left_idx[0] = idx_u;
      left_idx[1] = idx_v;
      basisValues.resize(size);
    }
  };

  /// Structure for storage of results of grid evaluation of the basis function of a spline surface.
  /// Position and first derivatives
  struct BasisDerivsSf
  {
    /// Parameter double in which the basis functions are evaluated
    double param[2];
    /// Index of the knot interval where the parameter value is situated for all
    /// parameter directions. The indices of the non-zero basis functions are
    /// left_idx[i]-order[i]+1, ..., left_idx[i] for i=0,1
    int left_idx[2];
    /// The value of all basis functions, size equal to (degree_u+1)*(degree_v+1)
    std::vector< double > basisValues; 
    /// the derivative of all basis functions in u direction, same size as previous
    std::vector< double > basisDerivs_u; 
    /// the derivative of all basis functions in v direction, same size as previous
    std::vector< double > basisDerivs_v;

    void prepareDerivs(double u, double v, int idx_u, int idx_v, int size)
    {
      param[0] = u;
      param[1] = v;
      left_idx[0] = idx_u;
      left_idx[1] = idx_v;
      basisValues.resize(size);
      basisDerivs_u.resize(size);
      basisDerivs_v.resize(size);
    }
  };

  /// Structure for storage of results of grid evaluation of the basis function of a spline surface.
  /// Position, first and second derivatives
  struct BasisDerivsSf2
  {
    /// Parameter double in which the basis functions are evaluated
    double param[2];   
    /// Index of the knot interval where the parameter value is situated for all
    /// parameter directions. The indices of the non-zero basis functions are
    /// left_idx[i]-order[i]+1, ..., left_idx[i] for i=0,1
    int left_idx[2];   
    /// The value of all basis functions, size equal to (degree_u+1)*(degree_v+1)
    std::vector< double > basisValues; 
  
    /// the derivative of all basis functions in u direction, same size as previous
    std::vector< double > basisDerivs_u;
    /// the derivative of all basis functions in v direction, same size as previous
    std::vector< double > basisDerivs_v;

    /// the second derivative of all basis functions twice in u direction, same size as previous
    std::vector< double > basisDerivs_uu;
    /// the second derivative of all basis functions in u and v direction, same size as previous
    std::vector< double > basisDerivs_uv;
    /// the second derivative of all basis functions twice in v direction, same size as previous
    std::vector< double > basisDerivs_vv;

    void prepareDerivs(double u, double v, int idx_u, int idx_v,
		       int size)
    {
      param[0] = u;
      param[1] = v;
      left_idx[0] = idx_u;
      left_idx[1] = idx_v;
      basisValues.resize(size);
      basisDerivs_u.resize(size);
      basisDerivs_v.resize(size);
      basisDerivs_uu.resize(size);
      basisDerivs_uv.resize(size);
      basisDerivs_vv.resize(size);
    }
  };

  /// Structure for storage of results of grid evaluation of the basis function of a spline surface.
  /// Position, first, second and third derivatives
  struct BasisDerivsSf3
  {
    /// Parameter double in which the basis functions are evaluated
    double param[2];
    /// Index of the knot interval where the parameter value is situated for all
    /// parameter directions. The indices of the non-zero basis functions are
    /// left_idx[i]-order[i]+1, ..., left_idx[i] for i=0,1
    int left_idx[2];
    /// The value of all basis functions, size equal to (degree_u+1)*(degree_v+1)
    std::vector< double > basisValues;

    /// the derivative of all basis functions in u direction, same size as previous
    std::vector< double > basisDerivs_u;
    /// the derivative of all basis functions in v direction, same size as previous
    std::vector< double > basisDerivs_v;

    /// the second derivative of all basis functions twice in u direction, same size as previous
    std::vector< double > basisDerivs_uu;
    /// the second derivative of all basis functions in u and v direction, same size as previous
    std::vector< double > basisDerivs_uv;
    /// the second derivative of all basis functions twice in v direction, same size as previous
    std::vector< double > basisDerivs_vv;

    std::vector< double > basisDerivs_uuu;
    std::vector< double > basisDerivs_uuv;
    std::vector< double > basisDerivs_uvv;
    std::vector< double > basisDerivs_vvv;

    void prepareDerivs(double u, double v, int idx_u, int idx_v, int size)
    {
      param[0] = u;
      param[1] = v;
      left_idx[0] = idx_u;
      left_idx[1] = idx_v;
      basisValues.resize(size);
      basisDerivs_u.resize(size);
      basisDerivs_v.resize(size);
      basisDerivs_uu.resize(size);
      basisDerivs_uv.resize(size);
      basisDerivs_vv.resize(size);
      basisDerivs_uuu.resize(size);
      basisDerivs_uuv.resize(size);
      basisDerivs_uvv.resize(size);
      basisDerivs_vvv.resize(size);
    }
  };

  /// Structure for storage of results of grid evaluation of the basis function of a spline surface.
  /// Position, and a given number of uni-directed derivatives
  struct BasisDerivsSfU
  {
    size_t derivs;
    /// Parameter double in which the basis functions are evaluated
    double param[2];
    /// Index of the knot interval where the parameter value is situated for all
    /// parameter directions. The indices of the non-zero basis functions are
    /// left_idx[i]-order[i]+1, ..., left_idx[i] for i=0,1
    int left_idx[2];

    /// The value of all basis functions and derivatives
    std::vector< std::vector<double> > values;

    /// Resize data structures
    /// \param u u parameter
    /// \param u v parameter
    /// \param idx_u u index
    /// \param idx_u v index
    /// \param deriv Number of derivatives
    /// \param size Number of functions
    void prepareDerivs(double u, double v, int idx_u, int idx_v,
		       size_t deriv, size_t size)
    {
      derivs = deriv;
      param[0] = u;
      param[1] = v;
      left_idx[0] = idx_u;
      left_idx[1] = idx_v;
      values.resize(size);
      for (size_t i=0;i<size;++i)
	values[i].resize(2*derivs+1);
    }
  };
} // end namespace
#endif

