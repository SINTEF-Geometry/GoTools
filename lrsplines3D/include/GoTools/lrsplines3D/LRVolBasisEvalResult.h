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

#ifndef _LRVOLBASISEVALRESULT_H
#define  _LRVOLBASISEVALRESULT_H

#include "GoTools/utils/config.h"
#include "GoTools/lrsplines3D/Element3D.h"

namespace Go
{
  /// Structure for storage of results of grid evaluation of the basis function of a spline volume.
  /// Positional evaluation information in one parameter value
  struct LRBasisPtsVol
  {
    /// Parameter tripple in which the basis functions are evaluated
    double param[3];

    /// Element associated with the parameter pair
    Element3D* element;

     /// Associated B-spline
    std::vector<LRBSpline3D*> bsplines;

     /// The value of all basis functions
    std::vector< double > basisValues; 

    void preparePts(double u, double v, double w, Element3D* elem, size_t size)
    {
      param[0] = u;
      param[1] = v;
      param[2] = w;
      element = elem;
      bsplines.resize(size);
      basisValues.resize(size);
    }
  };

  /// Structure for storage of results of grid evaluation of the basis function of a spline volume.
  /// Position and first derivatives
  struct LRBasisDerivsVol
  {
    /// Parameter double in which the basis functions are evaluated
    double param[3];

     /// Element associated with the parameter pair
    Element3D* element;

    /// Associated B-spline
    std::vector<LRBSpline3D*> bsplines;

   /// The value of all basis functions
    std::vector< double > basisValues; 
    /// the derivative of all basis functions in u direction, same size as previous
    std::vector< double > basisDerivs_u; 
    /// the derivative of all basis functions in v direction, same size as previous
    std::vector< double > basisDerivs_v;
    /// the derivative of all basis functions in w direction, same size as previous
    std::vector< double > basisDerivs_w;

    void prepareDerivs(double u, double v, double w, Element3D* elem, size_t size)
    {
      param[0] = u;
      param[1] = v;
      param[2] = w;
      element = elem;
      bsplines.resize(size);
      basisValues.resize(size);
      basisDerivs_u.resize(size);
      basisDerivs_v.resize(size);
      basisDerivs_w.resize(size);
    }
  };

  /// Structure for storage of results of grid evaluation of the basis function of a spline volume.
  /// Position, first and second derivatives
  struct LRBasisDerivsVol2
  {
    /// Parameter double in which the basis functions are evaluated
    double param[3];   

     /// Element associated with the parameter pair
    Element3D* element;

    /// Associated B-spline
    std::vector<LRBSpline3D*> bsplines;

     /// The value of all basis functions
    std::vector< double > basisValues; 
  
    /// the derivative of all basis functions in u direction, same size as previous
    std::vector< double > basisDerivs_u;
    /// the derivative of all basis functions in v direction, same size as previous
    std::vector< double > basisDerivs_v;
    /// the derivative of all basis functions in w direction, same size as previous
    std::vector< double > basisDerivs_w;

    /// the second derivative of all basis functions twice in u direction, same size as previous
    std::vector< double > basisDerivs_uu;
    /// the second derivative of all basis functions in u and v direction, same size as previous
    std::vector< double > basisDerivs_uv;
     /// the second derivative of all basis functions in u and w direction, same size as previous
    std::vector< double > basisDerivs_uw;
   /// the second derivative of all basis functions twice in v direction, same size as previous
    std::vector< double > basisDerivs_vv;
    /// the second derivative of all basis functions in v and w direction, same size as previous
    std::vector< double > basisDerivs_vw;
   /// the second derivative of all basis functions twice in w direction, same size as previous
    std::vector< double > basisDerivs_ww;

    void prepareDerivs(double u, double v, double w, Element3D* elem, size_t size)
    {
      param[0] = u;
      param[1] = v;
      param[2] = w;
      element = elem;
      bsplines.resize(size);
      basisValues.resize(size);
      basisDerivs_u.resize(size);
      basisDerivs_v.resize(size);
      basisDerivs_w.resize(size);
      basisDerivs_uu.resize(size);
      basisDerivs_uv.resize(size);
      basisDerivs_uw.resize(size);
      basisDerivs_vv.resize(size);
      basisDerivs_vw.resize(size);
      basisDerivs_ww.resize(size);
    }
  };

  /// Structure for storage of results of grid evaluation of the basis function of a spline volume.
  /// Position, first, second and third derivatives
  struct LRBasisDerivsVol3
  {
    /// Parameter double in which the basis functions are evaluated
    double param[3];

     /// Element associated with the parameter pair
    Element3D* element;

    /// Associated B-spline
    std::vector<LRBSpline3D*> bsplines;

     /// The value of all basis functions, size equal to (degree_u+1)*(degree_v+1)
    std::vector< double > basisValues;

    /// the derivative of all basis functions in u direction, same size as previous
    std::vector< double > basisDerivs_u;
    /// the derivative of all basis functions in v direction, same size as previous
    std::vector< double > basisDerivs_v;
    /// the derivative of all basis functions in w direction, same size as previous
    std::vector< double > basisDerivs_w;

    /// the second derivative of all basis functions twice in u direction, same size as previous
    std::vector< double > basisDerivs_uu;
    /// the second derivative of all basis functions in u and v direction, same size as previous
    std::vector< double > basisDerivs_uv;
     /// the second derivative of all basis functions in u and w direction, same size as previous
    std::vector< double > basisDerivs_uw;
    /// the second derivative of all basis functions twice in v direction, same size as previous
    std::vector< double > basisDerivs_vv;
    /// the second derivative of all basis functions in v and w direction, same size as previous
    std::vector< double > basisDerivs_vw;
   /// the second derivative of all basis functions twice in w direction, same size as previous
    std::vector< double > basisDerivs_ww;

    std::vector< double > basisDerivs_uuu;
    std::vector< double > basisDerivs_uuv;
    std::vector< double > basisDerivs_uuw;
    std::vector< double > basisDerivs_uvv;
    std::vector< double > basisDerivs_uvw;
    std::vector< double > basisDerivs_uww;
    std::vector< double > basisDerivs_vvv;
    std::vector< double > basisDerivs_vvw;
    std::vector< double > basisDerivs_vww;
    std::vector< double > basisDerivs_www;

    void prepareDerivs(double u, double v, double w, Element3D* elem, size_t size)
    {
      param[0] = u;
      param[1] = v;
      param[2] = w;
      element = elem;
      bsplines.resize(size);
      basisValues.resize(size);
      basisDerivs_u.resize(size);
      basisDerivs_v.resize(size);
      basisDerivs_w.resize(size);
      basisDerivs_uu.resize(size);
      basisDerivs_uv.resize(size);
      basisDerivs_uw.resize(size);
      basisDerivs_vv.resize(size);
      basisDerivs_vw.resize(size);
      basisDerivs_ww.resize(size);
      basisDerivs_uuu.resize(size);
      basisDerivs_uuv.resize(size);
      basisDerivs_uuw.resize(size);
      basisDerivs_uvv.resize(size);
      basisDerivs_uvw.resize(size);
      basisDerivs_uww.resize(size);
      basisDerivs_vvv.resize(size);
      basisDerivs_vvw.resize(size);
      basisDerivs_vww.resize(size);
      basisDerivs_www.resize(size);
    }
  };

  /// Structure for storage of results of grid evaluation of the basis function of a spline volume.
  /// Position, and a given number of uni-directed derivatives
  struct LRBasisDerivsVolU
  {
    size_t derivs;
    /// Parameter double in which the basis functions are evaluated
    double param[3];

    /// Element associated with the parameter pair
    Element3D* element;

   /// Associated B-spline
    std::vector<LRBSpline3D*> bsplines;

 
    /// The value of all basis functions and derivatives
    std::vector< std::vector<double> > values;

    /// Resize data structures
    /// \param u u parameter
    /// \param u v parameter
    /// \param deriv Number of derivatives
    /// \param size Number of functions
    void prepareDerivs(double u, double v, double w, size_t deriv, Element3D* elem, size_t size)
    {
      derivs = deriv;
      param[0] = u;
      param[1] = v;
      param[2] = w;
      element = elem;
      bsplines.resize(size);
      values.resize(size);
      for (size_t i=0;i<size;++i)
	values[i].resize(3*derivs+1);
    }
  };
} // end namespace
#endif

