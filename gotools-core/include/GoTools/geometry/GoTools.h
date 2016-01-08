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
#include <GoTools/geometry/SurfaceOfLinearExtrusion.h>
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

    /// Get knot tolerance
    static double knotEpsilon();

    /// Set knot tolerance
    static void setKnotEpsilon(double knot_epsilon);

private:
    // Tolerances
    static double space_epsilon_;
    static double parameter_epsilon_;
    static double knot_epsilon_;
};

} // namespace Go

#endif /* GOTOOLS_H_ */
