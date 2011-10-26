//===========================================================================
//                                                                           
// File: makeTestCylinderPlane.C                                             
//                                                                           
// Created: Fri Jul 30 15:38:31 2004                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: makeTestCylinderPlane.C,v 1.3 2005-09-22 15:01:29 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

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
