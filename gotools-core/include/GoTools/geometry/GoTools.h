/*
 * GoTools.h
 *
 *  Created on: Nov 1, 2010
 *      Author: jbt
 */

#ifndef GOTOOLS_H_
#define GOTOOLS_H_


#define GO_VERSION_MAJOR 4
#define GO_VERSION_MINOR 0
#define GO_VERSION_PATCH 1


#include <string>
#include "GoTools/geometry/ClassType.h"

#include <GoTools/geometry/SplineCurve.h>
#include <GoTools/geometry/CurveOnSurface.h>
#include <GoTools/geometry/Line.h>
#include <GoTools/geometry/Circle.h>
#include <GoTools/geometry/Ellipse.h>
#include <GoTools/geometry/BoundedCurve.h>
#include <GoTools/geometry/Hyperbola.h>
#include <GoTools/geometry/Parabola.h>

#include <GoTools/geometry/SplineSurface.h>
#include <GoTools/geometry/BoundedSurface.h>
#include <GoTools/geometry/CompositeSurface.h>
#include <GoTools/geometry/Plane.h>
#include <GoTools/geometry/Cylinder.h>
#include <GoTools/geometry/Sphere.h>
#include <GoTools/geometry/Cone.h>
#include <GoTools/geometry/Torus.h>
#include <GoTools/geometry/SurfaceOfRevolution.h>
#include <GoTools/geometry/Disc.h>

#include <GoTools/geometry/PointCloud.h>
#include <GoTools/geometry/LineCloud.h>

#include <GoTools/geometry/RectGrid.h>

namespace Go {

/// Class containing some service functions for the GoTools "system".
class GoTools
{
public:
    GoTools();
    virtual ~GoTools();

    /// This functions sets up the factories for constructing GeomObjects
    /// at run-time (i.e. Registrator instantiations). If you need to
    /// construct GeomObjects at run-time, such as when reading from a file,
    /// then call this function earøy in the code.
    static void init();

    /// Conversion function that returns a std::string with the class name
    /// corresponding to ClassType class_type.
    static std::string className(ClassType class_type);

};

} // namespace Go

#endif /* GOTOOLS_H_ */
