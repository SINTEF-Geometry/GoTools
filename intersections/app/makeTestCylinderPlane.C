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

#include <fstream>

#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/geometry/SplineSurface.h"
#include "sislP.h"

char cyl_file[] = "data/cylinder.g2";
char plane_file[] = "data/plane.g2";
char sphere_file[] = "data/sphere.g2";

using namespace std;
using namespace Go;

int main()
{
    // making plane spline surface, laying in 0xy
    double knots[] = {0, 0, 1, 1};
    double coefs[] = {0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0};

    SplineSurface planesurf(2, 2, 2, 2, knots, knots, coefs, 3, false);

    ofstream os_plane(plane_file);
    planesurf.writeStandardHeader(os_plane);
    planesurf.write(os_plane);
    os_plane.close();

    // making cylinder that is tangent to the plane
    double bottom_pos[] = {0, 0, -1};
    double bottom_axis[] = {1, 0, 0};
    double ellipse_ratio = 1;
    double axis_dir[] = {0, 0, 1};
    double height = 2;
    SISLSurf *cyl;
    int stat;

    s1021(bottom_pos, bottom_axis, ellipse_ratio, axis_dir, height, &cyl, &stat);

    SplineSurface* ssurf = SISLSurf2Go(cyl);
    
    ofstream os_cyl(cyl_file);
    ssurf->writeStandardHeader(os_cyl);
    ssurf->write(os_cyl);
    os_cyl.close();

    delete ssurf;
    freeSurf(cyl);

    // making sphere that touches the plane
    double centre[] = {0, 0, 0};
    //double axis[] = {0, 0, 1};
    double axis[] = {1, 0, 1};
    //double equator[] = {1, 0, 0};
    double equator[] = {0, 1, 0};
    SISLSurf *sphere;
    s1023(centre, axis, equator, 2, 4, &sphere, &stat);
    SplineSurface* sphere_surf = SISLSurf2Go(sphere);
    ofstream os_sphere(sphere_file);
    sphere_surf->writeStandardHeader(os_sphere);
    sphere_surf->write(os_sphere);
    os_sphere.close();
    delete sphere_surf;
    freeSurf(sphere);

    return 1;
}
