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

#ifndef _IMPLICITAPPROX_H_
#define _IMPLICITAPPROX_H_

#include "GoTools/implicitization/BernsteinTetrahedralPoly.h"
#include "GoTools/utils/BaryCoordSystem.h"
#include "GoTools/utils/Point.h"

namespace Go
{
  class RevEngPoint;

  /**
     ImlicitApprox - Interface class to Go::implicitization for reverse engineering
     purposes. Used only to approximate a point cloud with a plane (degree one). Higher
     degree approximations are unstable for noisy point clouds
   *
   */
  class ImplicitApprox
  {
  public:
    /// Constructor
    ImplicitApprox();

    /// Destructor
    ~ImplicitApprox();

    /// Approximate a group of RevEngPoints with a surface of degree degree.
    /// Recommended only for degree one
    void approx(std::vector<RevEngPoint*> points, int degree);

     /// Approximate several groups of RevEngPoints with a surface of degree degree.
    /// Recommended only for degree one
   void approx(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
			     std::vector<RevEngPoint*>::iterator> >& points,
		int degree);

    /// Approximate a group of points with a surface of degree degree.
    /// Recommended only for degree one
    void approxPoints(std::vector<Point> points, int degree);

    /// Project initial point and direction onto implicit surface to get a point
    /// in the surface and corresponding surface normal
    bool projectPoint(Point point, Point dir,
		      Point& projpos, Point& normal);

    /// \param pt Input point
    /// \param val Implicit distance between the point pt and the surface
    /// \param grad Gradient at point pt
    void evaluate(Point& pt, double& val, Point& grad);
    
    
  private:
    /// Degree of implicit surface
    int degree_;

    /// Bernstein polynomials on a tetrahedron
    BernsteinTetrahedralPoly implicit_;

    /// implicit_ differentiated
    BernsteinTetrahedralPoly deriv1_, deriv2_, deriv3_, deriv4_;

    // Encapsulates a barycentric coordinate system.
    BaryCoordSystem3D bc_;
    double sigma_min_;
    double eps_;
    
    double estimateDist(RevEngPoint* pt);

    void visualize(std::vector<RevEngPoint*> points, std::ostream& os);

    void visualize(std::vector<Point> points, Point& dir, std::ostream& os);

    void polynomialSurf(std::vector<Point>& pos_and_der, int degree,
			std::vector<double>& coefs);

    void polynomialSurfAccuracy(std::vector<Point>& pos_and_der, 
				int degree, std::vector<double>& coefs,
				double& maxfield, double& avfield,
				double& maxdist, double& avdist,
				int& ndiv, double& maxang,
				double& avang);
  };
}

#endif // _IMPLICITAPPROX_H_
