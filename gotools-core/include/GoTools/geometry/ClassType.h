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

#ifndef _CLASSTYPE_H
#define _CLASSTYPE_H

namespace Go {

/** All concrete classes that inherit GeomObject should
 *  have a matching enum in ClassType. It is necessary to
 *  have a central location for the type identification numbers,
 *  so that two classes in different modules (for example) does
 *  not try to use the same number. At the same time, it is
 *  something we'd like to avoid, because it forces you to
 *  update this file every time you add a new GeomObject-
 *  inheriting class.
 *
 *  If you need to add a new type, you may also need to add a "Registrator",
 *  see the 'GoTools' class. If the new class is part of 'geometry' in
 *  'gotools-core', do this in the init() function of the GoTools class
 *  itself. If it's outside the core, then write another class with its own
 *  static init() function containing the Registrator.
 *
 *  Also, if you need the name of that added type, i.e. a ClassType-to-string
 *  conversion, then you should add that as well. See GoTools::className().
 */
enum ClassType
{
    Class_Unknown = 0,

    Class_SplineCurve = 100,
    Class_CurveOnSurface = 110,
    Class_CurveOnVolume = 111,
    Class_Line = 120,
    Class_Circle = 130,
    Class_Ellipse = 140,
    Class_BoundedCurve = 150,
    Class_Hyperbola = 160,
    Class_Parabola = 170,

    Class_SplineSurface = 200,
    Class_BoundedSurface = 210,
    Class_SurfaceOnVolume = 211,
    Class_GoBaryPolSurface = 220,
    Class_GoHBSplineParamSurface = 230,
    Class_CompositeSurface = 240,
    Class_Plane = 250,
    Class_Cylinder = 260,
    Class_SurfaceOfLinearExtrusion = 261,
    Class_Sphere = 270,
    Class_Cone = 280,
    Class_Torus = 290,
    Class_SurfaceOfRevolution = 291,
    Class_Disc = 292,
    Class_LRSplineSurface = 293,
    Class_TSplineSurface = 294,

    Class_Go3dsObject = 300,
    Class_GoHeTriang = 310,
    Class_GoSdTriang = 320,
    Class_GoQuadMesh = 330,
    Class_GoHybridMesh = 340,
    Class_ParamTriang = 350,
    Class_GoVrmlGeometry = 360,

    Class_PointCloud = 400,
    Class_LineCloud = 410,

    Class_GoTriangleSets = 500,
    Class_RectGrid = 510,

    Class_SplineVolume = 700,
    Class_BoundedVolume = 710,
    Class_Parallelepiped = 720,
    Class_SphereVolume = 721,
    Class_CylinderVolume = 722,
    Class_ConeVolume = 723,
    Class_TorusVolume = 724,
    Class_LRSplineVolume = 793
};

}

#endif // _CLASSTYPE_H
