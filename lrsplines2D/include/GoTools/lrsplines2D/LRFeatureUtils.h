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


#ifndef _LFEATUREUTILS_H
#define _LRFEATUREUTILS_H


#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include <iostream>


namespace Go
{
  /// Given a current LR B-spline surface with associated point cloud,
  /// compute feature output in ncell x ncell grid.
  /// Features (number of column in associated grid cell):
    /// \param 0: Average slope in cell (9 samples); 
  ///  \param 1: Average value of surface in cell (9 samples); 
  /// \param 2: Maximum difference of surface values in cell (9 samples); 
  /// \param 3: Average distance between surface and points for each cell; 
  /// \param 4: Maximum distance between surface and points in cell;
  /// \param 5: Average intensity/height value of points in cell;
  /// \param 6: Maximum difference of intensity values in cell;
  /// \param 7: Standard deviation of distances between point cloud and surface in cell;
  /// \param 8: Standard deviation of intensity values in cell;
  /// \param 9: Average distance between surface and points in cell divided by maximum distance;
  ///  \param 10: Maximum difference between signed distances between points and surface in cell;
  /// \param 11: Average distance between points with higher intensity than the surface and surface in cell;
  /// \param 12: Average distance between points with lower intensity than the surface and surface in cell;
  /// \param 13: Number of point with lower intensity than the surface where the intensity difference is larger than threshold divided by the number of points in the cell;
  /// \param 14: Number of point with higher intensity than the surface where the intensity difference is larger than threshold divided by the number of points in the cell;
  /// \param 15: Number of surface elements in cell;
  /// \param 16: Average laplacian in cell (9 samples);
  /// The entries are scaled to represent a number in the range [0,10]. 

  namespace LRFeatureUtils
  {
    // Write accuracy features to file
    void writeCellInfo(const LRSplineSurface& srf, 
		       double tol, int ncell,
		       std::ostream &out);
  };

 }; // End namespace Go


#endif // _LRFEATUREUTILS_H
