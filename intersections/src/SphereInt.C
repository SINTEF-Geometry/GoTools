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

#include "GoTools/intersections/SphereInt.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SISLconversion.h"
#include "sislP.h"


using std::vector;
using std::string;


namespace Go {


//===========================================================================
SphereInt::SphereInt()
  : AlgObj3DInt(1)
//===========================================================================
{
}


//===========================================================================
SphereInt::SphereInt(Point center, double radius)
    : AlgObj3DInt(2), center_(center), radius_(radius)
//===========================================================================
{
    // All values in terms_ initialized to 0.0.
    terms_.push_back(Alg3DElem(1.0, 2, 0, 0));
    if (center_[0] != 0.0)
	terms_.push_back(Alg3DElem(-2.0*center_[0], 1, 0, 0));
    terms_.push_back(Alg3DElem(1.0, 0, 2, 0));
    if (center_[1] != 0.0)
	terms_.push_back(Alg3DElem(-2.0*center_[1], 0, 1, 0));
    terms_.push_back(Alg3DElem(1.0, 0, 0, 2));
    if (center_[2] != 0.0)
	terms_.push_back(Alg3DElem(-2.0*center_[2], 0, 0, 1));
    terms_.push_back(Alg3DElem(center_.length2() - radius_*radius_, 0, 0, 0));
}


//===========================================================================
SphereInt::~SphereInt()
//===========================================================================
{
}


//===========================================================================
void SphereInt::read(std::istream& is)
//===========================================================================
{
  // Expecting input to be on form
  // center_pt
  // radius
  // @@sbr Possibly separator characters ... As well as comments ...
  Point center_pt(3);
  double radius;
  char character;
  is >> character;
  center_pt.read(is);
  is >> character >> radius;

  center_ = center_pt;
  radius_ = radius;
}


//===========================================================================
Point SphereInt::center() const
//===========================================================================
{
    return center_;
}


//===========================================================================
double SphereInt::radius() const
//===========================================================================
{
    return radius_;
}


//===========================================================================
SplineSurface* SphereInt::surface() const
//===========================================================================
{
    SplineSurface* go_sphere = NULL;

    // The easy way is to create a SISL surface and then convert it to
    // the GO format.  We define the center axis to be (0.0, 0.0,
    // 1.0).
    Point axis(0.0, 0.0, 1.0);
    Point equator(radius_, 0.0, 0.0);
    Point center = center_;

    int latitude = 2;
    int longitude = 4;
    int kstat = 0;
    SISLSurf* sisl_sphere = NULL;
    s1023(&center[0], &axis[0], &equator[0], latitude, longitude,
	  &sisl_sphere, &kstat);
    if (kstat < 0) {
	THROW("Failed creating sphere!");
    } else {
	// We then convert the sisl sf to the go-format.
	go_sphere = SISLSurf2Go(sisl_sphere);
    }

    // As we're not using a smart_ptr we must free memory allocated.
    if (sisl_sphere) freeSurf(sisl_sphere);

    return go_sphere;
}


} // namespace Go
