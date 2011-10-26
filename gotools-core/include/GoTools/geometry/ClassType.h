//===========================================================================
//                                                                           
// File: ClassType.h                                                      
//                                                                           
// Created: Wed Nov  8 13:02:29 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: ClassType.h,v 1.31 2009-02-27 13:50:29 jbt Exp $
//                                                                           
//===========================================================================

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
    Class_TorusVolume = 724
};

}

#endif // _CLASSTYPE_H
