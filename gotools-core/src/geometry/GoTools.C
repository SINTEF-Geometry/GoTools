/*
 * GoTools.cpp
 *
 *  Created on: Nov 1, 2010
 *      Author: jbt
 */

#include <GoTools/geometry/GoTools.h>
#include "GoTools/geometry/Factory.h"

using std::string;

namespace Go {

GoTools::GoTools()
{
    // TODO Auto-generated constructor stub

}

GoTools::~GoTools()
{
    // TODO Auto-generated destructor stub
}

void GoTools::init()
{
    Registrator<SplineCurve> r100;
    Registrator<CurveOnSurface> r110;
    Registrator<Line> r120;
    Registrator<Circle> r130;
    Registrator<Ellipse> r140;
    Registrator<BoundedCurve> r150;
    Registrator<Hyperbola> r160;
    Registrator<Parabola> r170;

    Registrator<SplineSurface> r200;
    Registrator<BoundedSurface> r210;
    Registrator<CompositeSurface> r240;
    Registrator<Plane> r250;
    Registrator<Cylinder> r260;
    Registrator<Sphere> r270;
    Registrator<Cone> r280;
    Registrator<Torus> r290;
    Registrator<SurfaceOfRevolution> r291;
    Registrator<Disc> r292;

    Registrator<PointCloud<2> > r400_2d;
    Registrator<PointCloud<3> > r400_3d;
    Registrator<LineCloud> r410;

    Registrator<RectGrid> r510;

}

string GoTools::className(ClassType class_type)
{
    switch (class_type) {
    case Class_SplineCurve:
        return "SplineCurve";
    case Class_CurveOnSurface:
        return "CurveOnSurface";
    case Class_Line:
        return "Line";
    case Class_Circle:
        return "Circle";
    case Class_Ellipse:
        return "Ellipse";
    case Class_BoundedCurve:
        return "BoundedCurve";
    case Class_Hyperbola:
        return "Hyperbola";
    case Class_Parabola:
        return "Parabola";

    case Class_SplineSurface:
        return "SplineSurface";
    case Class_BoundedSurface:
        return "BoundedSurface";
    case Class_CompositeSurface:
        return "CompositeSurface";
    case Class_Plane:
        return "Plane";
    case Class_Cylinder:
        return "Cylinder";
    case Class_Sphere:
        return "Sphere";
    case Class_Cone:
        return "Cone";
    case Class_Torus:
        return "Torus";
    case Class_SurfaceOfRevolution:
        return "SurfaceOfRevolution";
    case Class_Disc:
        return "Disc";
    case Class_LRSplineSurface:
        return "LRSplineSurface";
    case Class_TSplineSurface:
        return "TSplineSurface";

    case Class_PointCloud:
        return "PointCloud";
    case Class_LineCloud:
        return "LineCloud";

    case Class_RectGrid:
        return "RectGrid";

    default:
        return "Unknown";
    }
}

} // namespace Go
