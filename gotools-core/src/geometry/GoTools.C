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
    Registrator<SurfaceOfLinearExtrusion> r261;
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
    case Class_SurfaceOfLinearExtrusion:
        return "SurfaceOfLinearExtrusion";
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


double GoTools::spaceEpsilon()
{
    return space_epsilon_;
}


void GoTools::setSpaceEpsilon(double space_epsilon)
{
    space_epsilon_ = space_epsilon;
}


double GoTools::parameterEpsilon()
{
    return parameter_epsilon_;
}


void GoTools::setParameterEpsilon(double parameter_epsilon)
{
    parameter_epsilon_ = parameter_epsilon;
}


double GoTools::knotEpsilon()
{
    return knot_epsilon_;
}


void GoTools::setKnotEpsilon(double knot_epsilon)
{
    knot_epsilon_ = knot_epsilon;
}


// Set default tolerances
double GoTools::space_epsilon_ = 0.001;
double GoTools::parameter_epsilon_ = 0.001;
double GoTools::knot_epsilon_ = 1.0e-5;


} // namespace Go
