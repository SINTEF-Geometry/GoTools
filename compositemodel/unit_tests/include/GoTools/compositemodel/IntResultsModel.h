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

#ifndef _INTRESULTSMODEL_H
#define _INTRESULTSMODEL_H

#include "GoTools/compositemodel/ftPoint.h"
#include "GoTools/compositemodel/ftCurve.h"
#include "GoTools/compositemodel/ftPlane.h"
#include "GoTools/compositemodel/ftLine.h"
#include "GoTools/geometry/PointCloud.h"
#include <vector>


namespace Go
{

enum IntersectionType
{
  No_Type = 0,
  SurfaceModel_SurfaceModel,
  SurfaceModel_Plane,
  SurfaceModel_Line,
  CompositeCurve_CompositeCurve,
  CompositeCurve_Plane,
  CompositeCurve_Line
};

class LineStrip;

//===========================================================================
/** Storage of results from intersection operations related to a CompositeModel
*/
//
//===========================================================================
class IntResultsModel
{
 public:

  /// Constructor
  IntResultsModel(IntersectionType type);

  /// Destructor
  ~IntResultsModel()
    {}

  /// Plane involved in the intersection
  void addPlaneInfo(const ftPlane& plane)
  {
    plane_ = plane;
  }

  /// Line involved in the intersection
  void addLineInfo(const ftLine& line)
  {
    line_ = line;
  }

  /// Check if any intersection curves are found
  virtual bool hasIntCurves() const = 0;

  /// Return the number of intersection curve segments
  virtual int nmbCurveSegments() const = 0;

  /// Check if any intersection point are found
  virtual bool hasIntPoints() const = 0;
  
  /// Return the number of intersection points
  virtual int nmbIntPoints() const = 0;

  /// Check if any intersections are found
  virtual bool hasIntersections() const = 0;
  
  /// Tesselation
  /// Tesselate with respect to a default value
  virtual void tesselate(std::vector<shared_ptr<LineStrip> >& meshes,
			 PointCloud3D& points) const = 0;

  /// Tesselate with respect to a given resolution, the same for each
  /// intersection curve
  virtual void tesselate(int resolution, 
			 std::vector<shared_ptr<LineStrip> >& meshes,
			 PointCloud3D& points) const = 0;

  /// Tesselate with respect to a given density
  virtual void tesselate(double density, 
			 std::vector<shared_ptr<LineStrip> >& meshes,
			 PointCloud3D& points) const = 0;

 protected:
  IntersectionType type_;
  ftLine line_;   // Is set if a line is involved in the current intersection
  ftPlane plane_; // Is set if a plane is involved in the current intersection
  int numpar_;    // Number of parameter values associated with an intersection result
};

} // namespace Go



#endif // _INTRESULTSMODEL_H
