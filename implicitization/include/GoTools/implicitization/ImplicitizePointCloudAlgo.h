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

#ifndef _IMPLICITIZEPOINTCLOUDALGO_H
#define _IMPLICITIZEPOINTCLOUDALGO_H


#include "GoTools/implicitization/BernsteinTetrahedralPoly.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/utils/BaryCoordSystem.h"


namespace Go {


/**
 * Class that implements an algorithm for finding an implicit surface
 * from a 3D point cloud.
 *
 * Input: A point cloud of type PointCloud3D, and the degree that is
 * chosen for the implicit representation.
 *
 * Output: A barycentric coordinate system, and the Bernstein
 * polynomial that represents the implicit surface.
 *
 * The algorithm first defines a barycentric coordinate system in
 * terms of a tetrahedron that is slightly larger than the bounding box
 * of the point cloud. It then computes the Bernstein polynomial
 * representing the implicit surface. If the chosen degree is too low to
 * produce the exact result, an approximate implicit surface is
 * produced. If the degree is too high, the resulting polynomial will
 * be reducible and contain the exact implicitization as a factor.
 */

class ImplicitizePointCloudAlgo {
public:
    /// Default constructor
    ImplicitizePointCloudAlgo() : tol_(3.0e-15) { }
    /// Constructor.
    /// \param deg degree of the implicit representation
    explicit ImplicitizePointCloudAlgo(int deg)
	: deg_(deg), tol_(3.0e-15) { }
    /// Constructor.
    /// \param cloud point cloud to be implicitized
    /// \param deg degree of the implicit representation
    ImplicitizePointCloudAlgo(const PointCloud3D& cloud, int deg)
	: cloud_(cloud), deg_(deg), tol_(3.0e-15) { }

    /// Load the point cloud to be implicitized.
    /// \param cloud point cloud on PointCloud3D form
    void usePointCloud(const PointCloud3D& cloud)
    { cloud_ = cloud; }

    /// Choose the degree of the implicit representation.
    /// \param deg degree
    void setDegree(int deg)
    { deg_ = deg; }

    /// Set the tolerance.
    /// \param tol tolerance; default value is 3.0e-15
    void setTolerance(double tol)
    { tol_ = tol; }

    /// Perform the implicitization.
    /// This function runs the implicitization algorithm.
    void perform();

    /// Get the result of the implicitization.
    /// \retval implicit a BernsteinTetrahedralPoly representing the
    /// implicit surface
    /// \retval bc the barycentric coordinate system in which the
    /// BernsteinTetrahedralPoly is defined
    void getResultData(BernsteinTetrahedralPoly& implicit,
		       BaryCoordSystem3D& bc, double& sigma_min)
    {
	implicit = implicit_;
	bc = bc_;
	sigma_min = sigma_min_;
    }

private:
    PointCloud3D cloud_;
    BernsteinTetrahedralPoly implicit_;
    BaryCoordSystem3D bc_;
    int deg_;
    double tol_;
    double sigma_min_;

};


} // namespace Go


#endif // _IMPLICITIZEPOINTCLOUDALGO_H

