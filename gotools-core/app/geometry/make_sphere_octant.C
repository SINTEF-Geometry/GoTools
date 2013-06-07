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
