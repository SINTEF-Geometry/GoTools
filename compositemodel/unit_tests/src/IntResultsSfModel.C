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

#include "GoTools/compositemodel/IntResultsSfModel.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/tesselator/LineStrip.h"

using namespace std;

namespace Go
{
  //===========================================================================
  // Constructor
  IntResultsSfModel::IntResultsSfModel(SurfaceModel* sfmodel, const ftLine& line)
  //===========================================================================
    : IntResultsModel(SurfaceModel_Line), sfmodel1_(sfmodel)
  {
    addLineInfo(line);
  }

  //===========================================================================
  // Constructor
  IntResultsSfModel::IntResultsSfModel(SurfaceModel* sfmodel, const ftPlane& plane)
  //===========================================================================
    : IntResultsModel(SurfaceModel_Plane), sfmodel1_(sfmodel)
  {
    addPlaneInfo(plane);
  }

  //===========================================================================
  // Destructor
  IntResultsSfModel::~IntResultsSfModel()
  //===========================================================================
  {
  }

  //===========================================================================
  void IntResultsSfModel::addIntPts(vector<ftPoint>& intpts)
  //===========================================================================
  {
    int_pts_.insert(int_pts_.end(), intpts.begin(), intpts.end());
  }


  //===========================================================================
  void IntResultsSfModel::addIntCvs(ftCurve& cvs)
  //===========================================================================
  {
    intcvs_ = cvs;
  }

  //===========================================================================
  void IntResultsSfModel::tesselate(std::vector<shared_ptr<LineStrip> >& meshes,
				    PointCloud3D& points) const
  //===========================================================================
  {
    if (hasIntCurves())
      intcvs_.tesselate(meshes);
    
    vector<double> coords;
    for (size_t ki=0; ki<int_pts_.size(); ++ki)
      {
	const Point pt = int_pts_[ki].position();
	coords.insert(coords.end(), pt.begin(), pt.end());
      }

    points = PointCloud3D(&coords[0], (int)int_pts_.size());
  }

  //===========================================================================
  void IntResultsSfModel::tesselate(int resolution, 
				    std::vector<shared_ptr<LineStrip> >& meshes,
				    PointCloud3D& points) const
  //===========================================================================
  {
    if (hasIntCurves())
      intcvs_.tesselate(resolution, meshes);
    
    vector<double> coords;
    for (size_t ki=0; ki<int_pts_.size(); ++ki)
      {
	const Point pt = int_pts_[ki].position();
	coords.insert(coords.end(), pt.begin(), pt.end());
      }

    points = PointCloud3D(&coords[0], (int)int_pts_.size());
  }

  //===========================================================================
  void IntResultsSfModel::tesselate(double density, 
				    std::vector<shared_ptr<LineStrip> >& meshes,
				    PointCloud3D& points) const
  //===========================================================================
  {
    if (hasIntCurves())
      intcvs_.tesselate(density, meshes);
    
    vector<double> coords;
    for (size_t ki=0; ki<int_pts_.size(); ++ki)
      {
	const Point pt = int_pts_[ki].position();
	coords.insert(coords.end(), pt.begin(), pt.end());
      }

    points = PointCloud3D(&coords[0], (int)int_pts_.size());
  }

} // namespace Go
