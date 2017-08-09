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

#ifndef _INTRESULTSSFMODEL_H
#define _INTRESULTSSFMODEL_H

#include "GoTools/compositemodel/IntResultsModel.h"
#include "GoTools/geometry/PointCloud.h"


namespace Go
{

  class SurfaceModel;
  class LineStrip;

//===========================================================================
/** Storage of results from intersection operations related to a SurfaceModel
*/
//
//===========================================================================

class IntResultsSfModel : public IntResultsModel
{
 public:

  /// Constructor
  // Note that only pointers are set to the initiating SurfaceModel
  // If this model is deleted, this pointer will be obsolete
  IntResultsSfModel(SurfaceModel *sfmodel, const ftLine& line);

  /// Constructor
  IntResultsSfModel(SurfaceModel *sfmodel, const ftPlane& plane);

  /// Destructor
  ~IntResultsSfModel();

  /// Add an intersection point
  void addIntPts(std::vector<ftPoint>& intpts);

  /// Add an intersection curve
  void addIntCvs(ftCurve& cvs);

  /// Check if any intersection curves are found
  virtual bool hasIntCurves() const
  {
    return (intcvs_.numSegments() > 0);
  }
  
  /// Return the number of intersection curve segments
  virtual int nmbCurveSegments() const
  {
    return intcvs_.numSegments();
  }

  /// Fetch all intersection curve segments joined in one curve
  const ftCurve& getIntersectionCurves()
  {
    return intcvs_;
  }

  /// Check if any intersection points are found
  virtual bool hasIntPoints() const
  {
    return (int_pts_.size() > 0);
  }
  
  /// Return the number of intersection points
  int nmbIntPoints() const 
  {
    return (int)int_pts_.size();
  }

  /// Fetch all intersection points
  std::vector<ftPoint>& getIntersectionPoints()
    {
      return int_pts_;
    }

  /// Check if any intersections are found
  virtual bool hasIntersections() const
  {
    return (intcvs_.numSegments() > 0 || int_pts_.size() > 0);
  }

  /// Tesselation
  /// Tesselate with respect to a default value
  virtual void tesselate(std::vector<shared_ptr<LineStrip> >& meshes,
			 PointCloud3D& points) const;

  /// Tesselate with respect to a given resolution, the same for each
  /// intersection curve
  virtual void tesselate(int resolution, 
			 std::vector<shared_ptr<LineStrip> >& meshes,
			 PointCloud3D& points) const;

  /// Tesselate with respect to a given density
  virtual void tesselate(double density, 
			 std::vector<shared_ptr<LineStrip> >& meshes,
			 PointCloud3D& points) const;

 private:
  SurfaceModel *sfmodel1_;
  SurfaceModel *sfmodel2_;

  std::vector<ftPoint> int_pts_;
  ftCurve intcvs_;
};

} // namespace Go



#endif // _INTRESULTSCOMPCV_H

