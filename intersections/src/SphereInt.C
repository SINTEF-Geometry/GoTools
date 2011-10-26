//===========================================================================
//                                                                           
// File: SphereInt.C                                                         
//                                                                           
// Created: Mon Jan 31 13:26:57 2005                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: SphereInt.C,v 1.9 2006-03-08 09:31:21 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


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
