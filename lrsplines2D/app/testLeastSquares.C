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

#include "GoTools/utils/config.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/creators/ApproxSurf.h" 
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSurfSmoothLS.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;
using std::vector;

int main(int argc, char *argv[])
{
  if (argc != 3) {
    std::cout << "Usage: point cloud (.g2) lrspline_out.g2 " << std::endl;
    return -1;
  }

  std::ifstream filein(argv[1]);
  std::ofstream fileout(argv[2]);
  
  ObjectHeader header;
  header.read(filein);
  PointCloud3D points;
  points.read(filein);

  BoundingBox box = points.boundingBox();
  Point low = box.low();
  Point high = box.high();

  // Create initial spline surface (cubic, one inner knot in each
  // parameter direction, dim=1)
  const double AEPSGE = 1.0e-3;
  vector<double> knots_u(9);
  vector<double> knots_v(9);
  vector<double> coefs(25, 0.0);

  int ki;
  for (ki=0; ki<4; ++ki)
    {
      knots_u[ki] = low[0];
      knots_v[ki] = low[1];
    }
  knots_u[4] = 0.5*(low[0]+high[0]);
  knots_v[4] = 0.5*(low[1]+high[1]);;
  
  for (ki=0; ki<4; ++ki)
    {
      knots_u[5+ki] = high[0];
      knots_v[5+ki] = high[1];
    }
      
  shared_ptr<SplineSurface> sf1(new SplineSurface(5, 5, 4, 4,  // order in v
						  &knots_u[0],     // ptr to knots in u
						  &knots_v[0],     // ptr to knots in v
						  &coefs[0], // ptr to coefs
						  1));       // one-dimensional

  // reformatting data to the format that ApproxSurf wants (separating parameter values and 
  // z-values in two distinct vectors).
  int nmb_pts = points.numPoints();
  double *datapoints = points.rawData();
  vector<double> zvals, param;
  zvals.reserve(nmb_pts);
  param.reserve(2 * nmb_pts);
  double* end = datapoints + 3 * nmb_pts;
  for (double* it = datapoints; it != end; ) {
    param.push_back(*it++);
    param.push_back(*it++);
    zvals.push_back(*it++);
  }

  // approximate initial surface
  double smoothing_term = 0.1;
  ApproxSurf asurf(sf1,
		   zvals, // z-values of scattered data points
		   param, // paramete values of scattered data points
		   1, // dimension
		   AEPSGE, // geometric tolerance
		   0, // 'constdir' - doesn't matter as we are not going to reparameterize anyway
		   false, // 'approx_orig' (not sure what this one does, better leave default value)
		   false, // 'close_belt' (leave at 'false', since we want all coefficients to be modified)
		   0,     // 'nmb_stabil' - not sure what it does, leave at default
		   false); // 'repar' - we set this to 'false' to avoid reparameterization of points (they
                           //           will be discarded anyway, so no reason to spend time on this).

  asurf.setSmoothingWeight(smoothing_term);
  asurf.setDoRefine(false);
  asurf.setFixBoundary(false);
  double maxdist, avdist; // will be set below
  int nmb_out_eps;        // will be set below
  shared_ptr<SplineSurface> result_surf = asurf.getApproxSurf(maxdist, avdist, nmb_out_eps, 1);

  std::cout << "Maxdist= " << maxdist << ", avdist= " << avdist;
  std::cout << ", nmb out= " << nmb_out_eps << std::endl;

  // Make 3D surface (not effective)
  vector<double> coefs2;
  int kj;
  vector<double>::iterator it;
  for (kj=0, it=result_surf->coefs_begin(); kj<result_surf->numCoefs_v(); ++kj)
    {
      double y = result_surf->basis_v().grevilleParameter((int)kj);
      for (ki=0; ki<result_surf->numCoefs_u(); ++ki, ++it)
	{
	  double x = result_surf->basis_u().grevilleParameter((int)ki);
	  coefs2.push_back(x);
	  coefs2.push_back(y);
	  coefs2.push_back(*it);
	}
    }
  shared_ptr<SplineSurface> sf2(new SplineSurface(result_surf->numCoefs_u(),
						  result_surf->numCoefs_v(),
						  result_surf->order_u(),
						  result_surf->order_v(),
						  result_surf->basis_u().begin(),
						  result_surf->basis_v().begin(),
						  &coefs2[0], // ptr to coefs
						  3));       // one-dimensional
  
  sf2->writeStandardHeader(fileout);
  sf2->write(fileout);

  // Make LR spline surface
  double knot_tol = 1.0e-6;

  shared_ptr<LRSplineSurface> lrsf(new LRSplineSurface(result_surf.get(), 
						       knot_tol));
  
  // Least squares approximation
  int nmb_basis = lrsf->numBasisFunctions();
  vector<int> coef_known(nmb_basis, 0);  // No coefficients are fixed
  
  LRSurfSmoothLS LSapprox(lrsf, coef_known);
  
  // Attach data points
  vector<double> datapts(datapoints, end);
  LSapprox.addDataPoints(datapts);
  
  double weight = 1.0;
  LSapprox.setLeastSquares(weight, 1.0);

  shared_ptr<LRSplineSurface> lrsf_out;
  int isOK = LSapprox.equationSolve(lrsf_out);
  std::cout << "isOK: " << isOK << std::endl;

  lrsf->to3D();
  lrsf->writeStandardHeader(fileout);
  lrsf->write(fileout);
}

