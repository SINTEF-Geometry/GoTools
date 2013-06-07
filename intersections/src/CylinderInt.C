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

#include "GoTools/intersections/CylinderInt.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SISLconversion.h"
#include "sislP.h"
#include <string>


using std::vector;


namespace Go {


//===========================================================================
CylinderInt::CylinderInt()
  : AlgObj3DInt(1)
//===========================================================================
{
}


//===========================================================================
CylinderInt::CylinderInt(Point ax_pt, Point ax_dir, double radius)
    : AlgObj3DInt(2), ax_pt_(ax_pt), ax_dir_(ax_dir), radius_(radius)
//===========================================================================
{
    setImplicitValues();
}


//===========================================================================
CylinderInt::~CylinderInt()
//===========================================================================
{
}


//===========================================================================
void CylinderInt::read(std::istream& is)
//===========================================================================
{
    // Expecting input to be on form
    // ax_pt
    // ax_dir
    // radius
    // @@sbr Possibly separator characters ... As well as comments ...
    Point ax_pt(3), ax_dir(3);
    double radius;
    std::string character;
    is >> character;
    ax_pt.read(is);
    is >> character;
    ax_dir.read(is);
    is >> character >> radius;

    ax_pt_ = ax_pt;
    ax_dir_ = ax_dir;
    radius_ = radius;

    setImplicitValues();
}


//===========================================================================
Point CylinderInt::ax_pt() const
//===========================================================================
{
    return ax_pt_;
}


//===========================================================================
Point CylinderInt::ax_dir() const
//===========================================================================
{
    return ax_dir_;
}


//===========================================================================
double CylinderInt::radius() const
//===========================================================================
{
    return radius_;
}


//===========================================================================
SplineSurface* CylinderInt::surface(Point bottom_pos, double height) const
//===========================================================================
{
  // Being lazy we use a sisl routine for creating the surface.
  SplineSurface* go_cylinder = NULL;

  // The easy way is to create a SISL surface and then convert it to the GO format.
  // We define the center axis to be (0.0, 0.0, 1.0).
//   Point axis(0.0, 0.0, 1.0);
//   Point equator(axis
//   double equator[3];
//   for (int ki = 0; ki < 3; ++ki)
//     equator[ki] = axis[ki]*radius_;

  int kstat = 0;
  SISLSurf* sisl_cylinder = NULL;
  // We create a vector which is not spanned by ax_dir_.
  Point ax_dir = ax_dir_;
  Point random_vec = ax_dir_;
  if (random_vec[0] == random_vec[1])
    random_vec[0] += 1.0;
  else
    std::swap(random_vec[0], random_vec[1]);
  Point bottom_axis = ax_dir_%random_vec;
  bottom_axis.normalize();
  bottom_axis *= radius_;
  Point bottom_pt = ax_pt_;
  s1021(&bottom_pos[0], &bottom_axis[0], 1.0, &ax_dir[0], height, &sisl_cylinder, &kstat);
  if (kstat < 0)
    {
      THROW("Failed creating cylinder!");
    }
  else
    {
      // We then convert the sisl sf to the go-format.
      go_cylinder = SISLSurf2Go(sisl_cylinder);
    }

  // As we're not using a smart_ptr we must free memory allocated.
  if (sisl_cylinder) freeSurf(sisl_cylinder);

  return go_cylinder;
}


//===========================================================================
void CylinderInt::setImplicitValues()
//===========================================================================
{
   // For easy construction of the implicit function we assume the ax_dir_
    // is one of the coordinate axes.
    ASSERT((ax_dir_[0]*ax_dir_[1] == 0.0) || (ax_dir_[0]*ax_dir_[2] == 0.0) ||
	   (ax_dir_[1]*ax_dir_[2] == 0.0));
    // We may thus assume that only one of the elements in ax_dir differs from 0.0.
    double sum_origo = -radius_*radius_; //0.0;
    if (ax_dir_[0] == 0.0) {
	terms_.push_back(Alg3DElem(1.0, 2, 0, 0));
	terms_.push_back(Alg3DElem(-2.0*ax_pt_[0], 1, 0, 0));
	sum_origo += ax_pt_[0]*ax_pt_[0];
    }
    if (ax_dir_[1] == 0.0) {
	terms_.push_back(Alg3DElem(1.0, 0, 2, 0));
	terms_.push_back(Alg3DElem(-2.0*ax_pt_[1], 0, 1, 0));
	sum_origo += ax_pt_[1]*ax_pt_[1];
    }
    if (ax_dir_[2] == 0.0) {
	terms_.push_back(Alg3DElem(1.0, 0, 0, 2));
	terms_.push_back(Alg3DElem(-2.0*ax_pt_[2], 0, 0, 1));
	sum_origo += ax_pt_[2]*ax_pt_[2];
    }
//     if (ax_dir_[0] != 0.0) {
// 	terms_.push_back(Alg3DElem(1.0, 1, 0, 0));
//     } else if (ax_dir_[1] != 0.0) {
// 	terms_.push_back(Alg3DElem(1.0, 0, 1, 0));
//     } else if (ax_dir_[2] != 0.0) {
// 	terms_.push_back(Alg3DElem(1.0, 0, 0, 1));
//     }
    terms_.push_back(Alg3DElem(sum_origo, 0, 0, 0));
}


} // namespace Go
