//===========================================================================
//                                                                           
// File: HermiteInterpolator.C                                              
//                                                                           
// Created: 01.08.30
//                                                                           
// Author: Vibeke Skytt, SINTEF
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================


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

