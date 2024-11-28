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

#ifndef _HEDGESURFACE_H
#define _HEDGESURFACE_H

#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/geometry/ClassType.h"
//#include "GoTools/compositemodel/RevEngRegion.h"

namespace Go {

class RevEngRegion;
  class CurveOnSurface;
  
  /// Additional information to the ClassType of a ParamSurface to distinguish
  /// between different rotational surfaces.  LINEARSWEPT_SURF is currently not
  /// active while ROTATIONALSWEPT_SURF is not implemented
  enum
    {
     SURF_TYPE_UNDEF, LINEARSWEPT_SURF, ROTATIONALSWEPT_SURF
    };
  
  /**  HedgeSurface - A topological surface associated with a RevEngRegion. 
       The class provides extra information and functionality compared to ftSurface, 
       mainly related to information about the RevEngRegion.
   *
   */

class HedgeSurface : public ftSurface
{
public:
  
  /// Constructor
  HedgeSurface();

  /// Constructor given geometry surface and associated region
  HedgeSurface(shared_ptr<ParamSurface> sf, RevEngRegion *region);

  /// Constructor given geometry surface and a number of associated regions. Not active
  HedgeSurface(shared_ptr<ParamSurface> sf, std::vector<RevEngRegion*>& region);

  /// Destructor
  ~HedgeSurface();

  /// Add information required to define a linear swept spline surface
  void setLinearSweepInfo(shared_ptr<SplineCurve> profile,
			  Point startpt, Point endpt)
  {
    surf_code_ = LINEARSWEPT_SURF;
    profile_ = profile;
    sweep1_ = startpt;
    sweep2_ = endpt;
  }
  
  /// Add information required to define a rotational swept spline surface. Not implemented
  void setRotationalSweepInfo(shared_ptr<SplineCurve> profile,
			      Point location, Point axis)
  {
    surf_code_ = ROTATIONALSWEPT_SURF;
    profile = profile_;
    sweep1_ = location;
    sweep2_ = axis;
  }

  /// Dimension of geometry space (should be 3)
  int dimension()
  {
    return surface()->dimension();
  }

  /// Number of points in associated region(s)
  int numPoints();

  /// Class type of geometry surface and possible swept surface code
  ClassType instanceType(int& code);

  /// Enquire if the surfae is a plane
  bool isPlane()
  {
    int code;
    return (instanceType(code) == Class_Plane);
  }

  /// Enquire if the surfae is a cylinder
  bool isCylinder()
  {
    int code;
    return (instanceType(code) == Class_Cylinder);
  }

  /// Enquire if the surfae is a sphere
  bool isSphere()
  {
    int code;
    return (instanceType(code) == Class_Sphere);
  }

  /// Enquire if the surfae is a torus
  bool isTorus()
  {
    int code;
    return (instanceType(code) == Class_Torus);
  }

  /// Enquire if the surfae is a cone
  bool isCone()
  {
    int code;
    return (instanceType(code) == Class_Cone);
  }

  /// Enquire if the surfae is a freeform surface
  bool isSpline()
  {
    int code;
    return (instanceType(code) == Class_SplineSurface);
  }

  /// Fetch associated region(s)
  std::vector<RevEngRegion*> getRegions()
  {
    return regions_;
  }

  /// Number of associated regions. Expected to be one
  int numRegions()
  {
    return (int)regions_.size();
  }

  /// Fetch specified region
  RevEngRegion* getRegion(int ix)
  {
    if (ix < 0 || ix >= (int)regions_.size())
      return 0;
    else
      return regions_[ix];
  }

  /// Add region to collection of associated regions
  void addRegion(RevEngRegion* reg);
  
  /// Remove region from collection of associated regions
  bool removeRegion(RevEngRegion* reg);

  /// Bounding box containing associated regions points
  BoundingBox regionsBox()
  {
    return bbox_;
  }

  /// Check if the geometry surfaces if entity and other is of the same type and has
  /// roughly the same characteristica. Is it a potential to merge surfaces?
  bool isCompatible(HedgeSurface* other, double angtol, double approx_tol,
		    ClassType& type, double& score);

  /// Ensure that the associated geometry surface is bounded (e.g. not an unlimited plane)
  void ensureSurfaceBounded();

  /// Bound unbounded primary surfaces 
  void limitSurf(double diag = -1.0);

  /// Make bounded surface when trimming edges are missing. Bound the associated region points
  /// in the parameter domain of the surface and transfer this information to this surface
  bool trimWithPoints(double aeps);

  /// Store current stage of hedge surface to file
  void store(std::ostream& os) const;

  /// Read hedge surface stage from file
  void read(std::istream& is);
    
private:
  /// Region(s) to which this surface is associated (only one)
  std::vector<RevEngRegion*> regions_;

  /// Bounding box of the associated region
  BoundingBox bbox_;

  /// Additional class type information to specify swept spline surfaces
  int surf_code_;

  /// The profile curve in a swept surface
  shared_ptr<SplineCurve> profile_;

  /// Sweep direction
  Point sweep1_;
  Point sweep2_;

  bool updateSurfaceWithAxis(Point axis[3], int ix, double tol, double angtol);
    
  bool updatePlaneWithAxis(Point axis[3], int ix, double tol, double angtol);
    
  bool updateCylinderWithAxis(Point axis[3], int ix, double tol, double angtol);
    
  bool checkAccuracyAndUpdate(shared_ptr<ParamSurface> surf, double tol,
			      double angtol);

  bool hasBaseSf();

  // Enquire if it is safe to intersect this surface and surf. Tangential intersections are
  // unstable and can produce infinite loops
  bool isTangential(HedgeSurface* surf);

  void doTrim(std::vector<shared_ptr<CurveOnSurface> >& int_cvs,
	      shared_ptr<BoundedSurface>& bdsf,
	      double tol,
	      std::vector<shared_ptr<HedgeSurface> >& added_sfs);

};
}

#endif // _HEDGESURFACE_H
