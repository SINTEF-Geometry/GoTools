//===========================================================================
//                                                                           
// File: make_sphere_octant.C                                                
//                                                                           
// Created: Tue Jul 26 11:04:50 2005                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: make_sphere_octant.C,v 1.2 2006-04-19 10:31:52 jbt Exp $
//                                                                           
//===========================================================================

#include "GoTools/geometry/SISLconversion.h"

int main()
{
/*
    double centre[] = { 0.0, 0.0, 0.0 };
    double axis[] = { 0.0, 0.0, 1.0 };
    double equator[] = { 1.0, 0.0, 0.0 };
    SISLSurf* sphere;
    int stat;
    s1023(centre, axis, equator, 1, 1, &sphere, &stat);
    shared_ptr<Go::SplineSurface> gsf(Go::SISLSurf2Go(sphere));
    double cutoff = 0.0;
    gsf.reset(gsf->subSurface(cutoff, gsf->startparam_v(),
			      gsf->endparam_u(), gsf->endparam_v()));
    
    gsf->writeStandardHeader(std::cout);
    gsf->write(std::cout);

    Go::Point p(3);
    gsf->point(p, cutoff, 0.0);
    double angle = std::acos(p[2]);
    std::cerr << "Point is " << p << std::endl;
    std::cerr << "Angle is " << angle << std::endl;
*/

    return 0;
}
