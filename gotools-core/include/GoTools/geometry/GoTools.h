/*
 * GoTools.h
 *
 *  Created on: Nov 1, 2010
 *      Author: jbt
 */

#ifndef GOTOOLS_H_
#define GOTOOLS_H_

#include <string>
#include "GoTools/geometry/GoTools_version.h"
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
    /// then call this function early in the code.
    static void init();

    /// Conversion function that returns a std::string with the class name
    /// corresponding to ClassType class_type.
    static std::string className(ClassType class_type);

    /// Get geometry tolerance
    static double spaceEpsilon();

    /// Set geometry tolerance
    static void setSpaceEpsilon(double space_epsilon);

    /// Get parameter tolerance
    static double parameterEpsilon();

    /// Set parameter tolerance
    static void setParameterEpsilon(double parameter_epsilon);

private:
    // Tolerances
    static double space_epsilon_;
    static double parameter_epsilon_;
};

} // namespace Go

#endif /* GOTOOLS_H_ */
