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

#include "GoTools/geometry/ElementaryUtils.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/ElementarySurface.h"
#include "GoTools/geometry/Plane.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/Cone.h"
#include "GoTools/geometry/Sphere.h"
#include "GoTools/geometry/Torus.h"

using namespace Go;

//===========================================================================
bool
 ElementaryUtils::sameElementarySurface(ParamSurface* under1,
					ParamSurface* under2,
					double tol, double angtol)

//===========================================================================
{
  bool same = false;

  // Check for equality of elementary surfaces
  ElementarySurface *elem1 = under1->elementarySurface();
  ElementarySurface *elem2 = under2->elementarySurface();
	      
  if (elem1 && elem2)
    {
      // Both surfaces are elementary. Check for equality
      if (elem1->instanceType() == Class_Plane &&
	  elem2->instanceType() == Class_Plane)
	{
	  Plane* plane1 = (Plane*)elem1;
	  Plane* plane2 = (Plane*)elem2;
	  Point pt1 = plane1->getPoint();
	  Point pt2 = plane2->getPoint();
	  Point norm1 = plane1->getNormal();
	  Point norm2 = plane2->getNormal();
	  double ang = norm1.angle(norm2);
	  if (ang < angtol || M_PI-ang < angtol)
	    {
	      double len = fabs((pt2 - pt1)*norm1);
	      if (len < tol)
		same = true;
	    }
	}
      else if (elem1->instanceType() == Class_Cylinder &&
	       elem2->instanceType() == Class_Cylinder)
	{
	  Cylinder* cyl1 = (Cylinder*)elem1;
	  Cylinder* cyl2 = (Cylinder*)elem2;
	  Point pt1 = cyl1->getLocation();
	  Point pt2 = cyl2->getLocation();
	  Point axis1 = cyl1->getAxis();
	  Point axis2 = cyl2->getAxis();
	  double rad1 = cyl1->getRadius();
	  double rad2 = cyl2->getRadius();
	  double ang = axis1.angle(axis2);
	  if (fabs(rad1-rad2) < tol && 
	      (ang < angtol || M_PI-ang < angtol))
	    {
	      double len = fabs((pt2 - pt1)*axis1);
	      if (len < tol)
		same = true;
	    }
	}
      else if (elem1->instanceType() == Class_Cone &&
	       elem2->instanceType() == Class_Cone)
	{
	  Cone* cone1 = (Cone*)elem1;
	  Cone* cone2 = (Cone*)elem2;
	  Point pt1 = cone1->getLocation();
	  Point pt2 = cone2->getLocation();
	  Point axis1 = cone1->getAxis();
	  Point axis2 = cone2->getAxis();
	  double rad1 = cone1->getRadius();
	  double rad2 = cone2->getRadius();
	  double angle1 = cone1->getConeAngle();
	  double angle2 = cone2->getConeAngle();
	  double ang = axis1.angle(axis2);
	  if (fabs(angle1-angle2) < angtol && 
	      (ang < angtol || M_PI-ang < angtol))
	    {
	      double len = fabs((pt2 - pt1)*axis1);
	      double d = pt1.dist(pt2);
	      double tanalpha = tan(angle1);
	      if (fabs(tanalpha) > tol)
		{
		  double d1 = rad1/tanalpha;
		  double d2 = rad2/tanalpha - d;
		  if (len < tol && fabs(d1-d2) < tol)
		    same = true;
		}
	      else if (fabs(rad1-rad2) < tol)
		same = true;
	    }
	}
      else if (elem1->instanceType() == Class_Sphere &&
	       elem2->instanceType() == Class_Sphere)
	{
	  Sphere* sphere1 = (Sphere*)elem1;
	  Sphere* sphere2 = (Sphere*)elem2;
	  Point pt1 = sphere1->getLocation();
	  Point pt2 = sphere2->getLocation();
	  double rad1 = sphere1->getRadius();
	  double rad2 = sphere2->getRadius();
	  double d = pt1.dist(pt2);
	  if (d < tol && fabs(rad1-rad2) < tol)
	    same = true;
	}
      else if (elem1->instanceType() == Class_Torus &&
	       elem2->instanceType() == Class_Torus)
	{
	  Torus* tor1 = (Torus*)elem1;
	  Torus* tor2 = (Torus*)elem2;
	  Point pt1 = tor1->getLocation();
	  Point pt2 = tor2->getLocation();
	  double radmin1 = tor1->getMinorRadius();
	  double radmin2 = tor2->getMinorRadius();
	  double radmax1 = tor1->getMajorRadius();
	  double radmax2 = tor2->getMajorRadius();
	  double d = pt1.dist(pt2);
	  Point x1, y1, z1, x2, y2, z2;
	  tor1->getCoordinateAxes(x1, y1, z1);
	  tor2->getCoordinateAxes(x2, y2, z2);
	  double ang1 = x1.angle(x2);
	  double ang2 = x1.angle(x2);
	  double ang3 = x1.angle(x2);
	  if (d < tol && fabs(radmin1-radmin2) < tol &&
	      fabs(radmax1-radmax2) < tol &&
	      (ang1 < angtol || M_PI-ang1 < angtol) &&
	      (ang2 < angtol || M_PI-ang2 < angtol) &&
	      (ang3 < angtol || M_PI-ang3 < angtol))
	    same = true;
	}
    }
  return same;
}
