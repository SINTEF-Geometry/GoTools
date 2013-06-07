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

#include "GoTools/geometry/HermiteInterpolator.h"
#include <vector>

using namespace Go;
using namespace std;

//===========================================================================
HermiteInterpolator::~HermiteInterpolator()
//===========================================================================
{
}

//===========================================================================
const BsplineBasis& HermiteInterpolator::basis()
//===========================================================================
{
    return basis_;
}

//===========================================================================
void HermiteInterpolator::interpolate(int num_points,
				       int dimension,
				       const double* param_start,
				       const double* data_start,
				       std::vector<double>& coefs)
//===========================================================================
{

  ALWAYS_ERROR_IF(num_points%2 != 0, "An odd number of points");
  for (int i = 1; i < num_points-1; i+=2) {
    ALWAYS_ERROR_IF(param_start[i] != param_start[i-1] ||
		param_start[i] >= param_start[i+1],
		    "Parameter sequence not fit for Hermite interpolation.");

  }
    
  // First we make a knot vector and define the spline space
  int order = 4;
  int num_coefs = num_points;
  ALWAYS_ERROR_IF(num_coefs < 4,
		  "Insufficient number of points.");


  std::vector<double> knots;
  knots.reserve(num_coefs + order);
  knots.insert(knots.end(), order, param_start[0]);
  for (int i = 2; i < num_points-3; i+=2)
    knots.insert(knots.end(), 2, param_start[i]);
  knots.insert(knots.end(), order, param_start[num_points-1]);
  basis_ = BsplineBasis(num_coefs, order, &knots[0]);

  if (int(coefs.size()) != dimension*num_coefs) {
    //	MESSAGE("Output coefficient vector had to be resized.");
    coefs.resize(dimension*num_coefs);
  }

  // Create the spline coefficients
  double scale;
  Point p1(dimension), p2(dimension), p3(dimension), p4(dimension);
  Point d1(dimension), d2(dimension);
  p1.setValue(data_start);
  d1.setValue(data_start+dimension);

  int h = 0;

  for (; h < dimension; ++h)
    coefs[h] = data_start[h];

  for (int i = 2; i < num_points; i+=2)
    {
      scale = (param_start[i] - param_start[i-1])/3.0;
      p4.setValue(data_start+i*dimension);
      d2.setValue(data_start+(i+1)*dimension);
      p2 = p1 + d1*scale;
      p3 = p4 - d2*scale;

      for (int j = 0;  j < dimension; ++j, ++h)
	coefs[h] = p2[j];
      for (int j = 0;  j < dimension; ++j, ++h)
	coefs[h] = p3[j];

      p1 = p4;
      d1 = d2;
    }

  for (int j = 0;  j < dimension; ++j, ++h)
    coefs[h] = p4[j];

}


//===========================================================================
void HermiteInterpolator::interpolate(const std::vector<Point>& data,
					const std::vector<double>& param,
					std::vector<double>& coefs)
//===========================================================================
{
  ALWAYS_ERROR_IF(2*param.size() != data.size(), 
	      "An odd number of points");
  for (size_t i1=1; i1<param.size(); i1++)
    ALWAYS_ERROR_IF(param[i1] <= param[i1-1],
		    "Parameter sequence must be strictly increasing.");


  // First we make a knot vector and define the spline space
  int order = 4;
  int num_coefs = (int)data.size();
  ALWAYS_ERROR_IF(num_coefs < 4,
		  "Insufficient number of points.");


  int nmb = (int)param.size()-1;
  std::vector<double> knots;
  knots.reserve(num_coefs + order);
  knots.insert(knots.end(), order, param[0]);
  for (int i1=1; i1<nmb; i1++)
    knots.insert(knots.end(), 2, param[i1]);
  knots.insert(knots.end(), order, param[nmb]);
  basis_ = BsplineBasis(num_coefs, order, &knots[0]);

  // Compute the coeficients
  int dim = data[0].dimension();
  double scale;
  Point p2(dim), p3(dim);
  if (int(coefs.size()) != dim*num_coefs) 
    coefs.resize(dim*num_coefs);
  int h = 0;
  for (int j=0; j<dim; j++)
    coefs[h++] = data[0][j];      
  for (int i1 = 0, i2=1; i2<=nmb; i1++, i2++)
    {
      scale = (param[i2] - param[i1])/3.0;
      p2 = data[2*i1] + data[2*i1+1]*scale;
      p3 = data[2*i2] - data[2*i2+1]*scale;

      for (int j=0; j<dim; j++)
	coefs[h++] = p2[j];      
      for (int j=0; j<dim; j++)
	coefs[h++] = p3[j];
    }
  for (int j=0; j<dim; j++)
    coefs[h++] = data[2*nmb][j]; 
}

