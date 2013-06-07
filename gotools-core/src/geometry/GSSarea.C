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

#include "GoTools/geometry/SplineSurface.h"


using namespace std;
using namespace Go;


namespace Go
{


//===========================================================================
double SplineSurface::area(double tol) const
//===========================================================================
{

  vector<double> c1Knots_u;
  vector<double> c1Knots_v;

  basis_u_.cNDiscontinuities(c1Knots_u, 2);
  basis_v_.cNDiscontinuities(c1Knots_v, 2);

  double result = 0.0;

  for (int i = 0; i < (int)c1Knots_u.size() - 1; ++i)
    {
      double start_u = c1Knots_u[i];
      double end_u = c1Knots_u[i+1];
      for (int j = 0; j < (int)c1Knots_v.size() - 1; ++j)
	{
	  double start_v = c1Knots_v[j];
	  double end_v = c1Knots_v[j+1];
	  result += areaDiagonalCross(tol, start_u, end_u, start_v, end_v);
	}
    }
  return result;
}




//===========================================================================
double SplineSurface::areaDiagonalCross(double tol,
					double start_u, double end_u,
					double start_v, double end_v) const
//===========================================================================
{
  vector <vector<double> > area;
  int termination_depth = 10;

  vector<double> points;
  vector<double> param_u;
  vector<double> param_v;

  int numPoints = 2;

  for (int i = 0; i < termination_depth; ++i)
    {
      // Get all points in the grids
      gridEvaluator(numPoints, numPoints, points, param_u, param_v, start_u, end_u, start_v, end_v);

      // Reset array with area for each subdomain
      vector <vector<double> > subArea;
      subArea.resize(numPoints);

      // Calculate subdomain areas
      for (int j = 0; j < numPoints-1; ++j)
	{
	  subArea[j].resize(numPoints);
	  for (int k = 0; k < numPoints-1; ++k)
	    {

	      // Get corner points
	      double* pointStart = &points[(j*numPoints + k)*3];
	      Point p00 = Point(pointStart, &pointStart[3], false);
	      Point p01 = Point(&pointStart[3], &pointStart[6], false);
	      pointStart = &pointStart[3*numPoints];
	      Point p10 = Point(pointStart, &pointStart[3], false);
	      Point p11 = Point(&pointStart[3], &pointStart[6], false);

	      // Calulate area
	      Point diagonalCross = (p11-p00) % (p01 - p10);
	      subArea[j][k] = diagonalCross.length() / 2.0;

	    }
	}

      // Sum up the areas of the subdomains. Sum up pairwise in each parameter direction to better avoid round-off errors
      for (int j = (numPoints-1) >> 1; j>0; j >>= 1)
	for (int k = 0; k < j; ++k)
	  for (int l = 0; l < j; ++l)
	    subArea[k][l] =
	      subArea   [k<<1]       [l<<1]
	      + subArea [(k<<1) | 1] [l<<1]
	      + subArea [k<<1]       [(l<<1) | 1]
	      + subArea [(k<<1) | 1] [(l<<1) | 1];

      // Insert first new entry in Romberg table, the calculated area
      area.resize(i+1);
      area[i].resize(i+1);
      area[i][0] = subArea[0][0];
      numPoints = 2 * numPoints - 1;
      if (i==0) continue;

      // Compute rest of Romberg table
      double fac = 0.0;
      for (int j = 1; j < i+1; ++j) {
	fac *= 4.0;
	fac += 3.0;
	area[i][j] = area[i][j-1] + (area[i][j-1] - area[i-1][j-1])/fac;
      }

      double err = abs(area[i][i] - area[i-1][i-1]);
      if (err < tol * area[i][i]) return area[i][i];

    }

  return area[termination_depth-1][termination_depth-1];

}



} // namspace Go
