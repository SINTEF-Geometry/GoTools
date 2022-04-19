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

#ifndef _PARAMETRIZATION_DOXYMAIN_H
#define _PARAMETRIZATION_DOXYMAIN_H

/**
\page streamable_doc The g2-format, GoTools file format for geometry entities

The following geometry entities have read and write functionality using the
g2-format:

\arg \c \em SplineCurve, entity code = 100,
\arg \c \em CurveOnSurface, entity code = 110,
\arg \c \em Line, entity code = 120,
\arg \c \em Circle, entity code = 130,
\arg \c \em Ellipse, entity code = 140,
\arg \c \em BoundedCurve, entity code = 150,
\arg \c \em Hyperbola, entity code = 160,
\arg \c \em Parabola, entity code = 170,
\arg \c \em SplineSurface, entity code = 200,
\arg \c \em BoundedSurface, entity code = 210,
\arg \c \em SurfaceOnVolume, entity code = 211,
\arg \c \em GoBaryPolSurface, entity code = 220,
\arg \c \em GoHBSplineParamSurface, entity code = 230,
\arg \c \em CompositeSurface, entity code = 240,
\arg \c \em Plane, entity code = 250,
\arg \c \em Cylinder, entity code = 260,
\arg \c \em Sphere, entity code = 270,
\arg \c \em Cone, entity code = 280,
\arg \c \em Torus, entity code = 290,
\arg \c \em SurfaceOfRevolution, entity code = 291,
\arg \c \em Disc, entity code = 292,
\arg \c \em LRSplineSurface, entity code = 293,
\arg \c \em TSplineSurface, entity code = 294,
\arg \c \em Go3dsObject, entity code = 300,
\arg \c \em GoHeTriang, entity code = 310,
\arg \c \em GoSdTriang, entity code = 320,
\arg \c \em GoQuadMesh, entity code = 330,
\arg \c \em GoHybridMesh, entity code = 340,
\arg \c \em ParamTriang, entity code = 350,
\arg \c \em GoVrmlGeometry, entity code = 360,
\arg \c \em PointCloud, entity code = 400,
\arg \c \em LineCloud, entity code = 410,
\arg \c \em GoTriangleSets, entity code = 500,
\arg \c \em RectGrid, entity code = 510,
\arg \c \em SplineVolume, entity code = 700,
\arg \c \em Parallelepiped, entity code = 720,
\arg \c \em SphereVolume, entity code = 721,
\arg \c \em CylinderVolume, entity code = 722,
\arg \c \em ConeVolume, entity code = 723,
\arg \c \em TorusVolume, entity code = 724
\arg \c \em LRSplineVolume, entity code = 793

In the following the file format corresponding to  most of the entities listed 
above will be described. Note that the g2-file format is a simple format that
does not support pointers. Thus, duplication of information may occur.

\section g2_sec_objectHeader The object header
All geometry entities are preceeded by the object 
\beginlink \link Go::ObjectHeader ObjectHeader\endlink in a stream.
ObjectHeader contains information about the class type of the
geometry object to folow and the version number of the file format. 
The format of the header
is as follows:
\arg \c <em> class type </em>. The entity code given in the list above.
\arg \c <em> major version </em>. Always 1.
\arg \c <em> minor version </em>. Always 0.
\arg \c <em> auxillerary data </em>. 0 if the default colour is chosen, 4 if the 
rgb-colour code is chosen. In the latter case, the entity is followed by 4 integer
digits between 0 and 255 giving the colours red, green and blue and the xxx.

\section g2_sec_splineCurve SplineCurve
The entity code 100 is given in the header. 
\beginlink \link Go::SplineCurve SplineCurve \endlink
continues with the following information
- dimension of geometry space, whether or not the curve is rational: 1=rational,
0=non-rational
- the number of coefficients, the polynomial order (i.e. degree+1)
- the knot vector, multiple knots are represented by giving the knot value
several times
- the curve coefficients

The curve coefficients are given continuously, normally divided by a carriage
return for each coefficient. The sequence is, in the non-rational case, 
\f$ (x_i, y_i, z_i), \quad
i=1, \ldots n \f$ where \f$ n \f$ is the number of coefficients. In the rational case,
the coefficients are given as \f$ (x_i h_i, y_i h_i, z_i h_i, h_i), \quad
i=1, \ldots n \f$. \f$h_i\f$ is the weight associated the current coefficient. 
Note that
the coefficients are represented with the weight multiplied with the \f$x\f$, 
\f$y\f$ and \f$x\f$-value.

A linear non-rational spline curve with no inner knots in the parameter interval 
\f$[0.1]\f$is including the header given as:
\verbatim
100 1 0 0 
3 0
2 2
0 0 1 1
0 0 0
1 0 0
\endverbatim

\section g2_sec_curveOnSurface CurveOnSurface
The entity enumeration of a \beginlink \link Go::CurveOnSurface CurveOnSurface \endlink
is 110. The remaining file format is as follows:
- whether or not the parameter space curve is the master information (1=parameter
curve is master, 0=space curve is master), the entity type of the parameter curve
(0 if no such curve exists), the entity type of the space curve (0 if no such
curve exists)
- the description of the parameter curve according to the rules for that curve
type
- the description of the space curve according to the rules for that curve type

Either the parameter curve, the space curve or both are given. Information
regarding constant parameter information for the curve on surface is not covered
by the file format.

Note that the CurveOnSurface entity does also contain information about a 
\beginlink \link Go::ParamSurface ParamSurface\endlink. Thus, the g2-representation
of a CurveOnSurface is not stand alone. It needs to be combined with a \ref g2_sec_boundedSurface. That section gives an example of how the CurveOnSurface entity is used
in this context.

\section g2_sec_line Line
The \beginlink \link Go::Line line \endlink is given by the following 
information:
- the dimension of the geometry space
- a point on the line 
- the direction of the line
- a flag, 0 or 1, for unbounded or bounded
- if the line is bounded, the start and end parameters

The entity enumeration for a line is 120.

\section g2_sec_circle Circle
The \beginlink \link Go::Circle circle\endlink is given by:
- the dimension of the geometry space
- the radius of the circle
- the centre of the circle
- the normal of the plane in which the circle lies
- the vector from the centre of the circle towards the "x-axis" - the
start point of the default parametrization of the circle
- the start and end parameters

The circle has a default implicit parametrization from 0 to \f$ 2 \pi \f$.

The entity enumeration for a circle is 130.

Example:
\verbatim
130 1 0 0
dim
rad
c_x c_y c_z
n_x n_y n_z
v_x v_y v_z
\endverbatim

\section g2_sec_ellipse Ellipse
The \beginlink \link Go::Ellipse ellipse \endlink is given by:
- the dimension of the geometry space
- the major radius of the ellipse
- the minor radius of the ellipse
- the centre of the ellipse
- the normal of the plane in which the ellipse lies
- the vector from the centre of the ellipse towards the "x-axis" - the
start point of the default parametrization

The entity enumeration for an ellipse is 130.

\section g2_sec_boundedCurve Bounded Curve
A \beginlink \link Go::BoundedCurve bounded curve\endlink are expected to be bounded 
both parametrically and geometrically. Thus, a parameter on whether or not 
the limitation in the
parameter space or in geometry space is preffered is included in the format. The
format is:
- Whether the limitation in the parameter space is the master (1 if true)
- The start parameter of the bounded curve, the end parameter of the bounded curve
- The start point of the bounded curve in geometry space
- The end point of the bounded curve in geometry space
- The underlying curve described according to its entity type

All the information listed above is expected to be given in a g2-file. The 
corresponding entity enumeration is 150.

A bounded line from the point (0,0,0) to (1,0,0) where the position in space
is the master regarding the bounding points is expressed as
\verbatim
150 1 0 0
120 0 
0 1
3
0 0 0
1 0 0

3
0 0 0
1 0 0
\endverbatim
120 is the entity number for a line and the dimension of the geometry
space is 3.

\section g2_sec_splineSurface SplineSurface
The entity enumeration for a \beginlink \link Go::SplineSurface spline surface \endlink
is 200. The body of a spline surface is as follows:
- dimension of geometry space, whether or not the surface is rational: 1=rational,
0=non-rational
- the number of coefficients in the first parameter direction, 
the polynomial order in this direction (i.e. degree+1)
- the knot vector in the first parameter direction, 
multiple knots are represented by giving the knot value several times
- the number of coefficients in the second parameter direction, 
the polynomial order in this direction (i.e. degree+1)
- the knot vector in the second parameter direction
- the surface coefficients

The surface coefficients are given continuously and the 1. parameter direction
runs fastest. The sequence is
\f$ (x_{(1,1)}, y_{(1,1)}, z_{(1,1)}), (x_{(2,1)}, y_{(2,1)}, z_{(2,1)}),
\ldots, (x_{(n,1)}, y_{(n,1)}, z_{(n,1)}), (x_{(1,2)}, y_{(1,2)}, z_{(1,2)}),
\ldots, (x_{(n,2)}, y_{(n,2)}, z_{(n,2)}), \ldots, (x_{(1,m)}, y_{(1,m)}, z_{(1,m)}),
\ldots,  (x_{(n,m)}, y_{(n,m)}, z_{(n,m)}) \f$. Here \f$n\f$ is the number of 
coefficients in the 1. parameter direction, and \f$m\f$ is the number of coefficients
in the 2. parameter direction. The weights are multiplied with the coefficients
similarily to the case for the spline curve and is given as an extra entity for
each coefficient, again similar to the curve case.

A simple non-rational, bilinear surface with no inner knots is represented as:
\verbatim
200 1 0 0 
3 0
2 2
0 0 1 1
2 2
0 0 1 1
0 0 0
1 0 0
0 1 0
1 1 0
\endverbatim

\section g2_sec_boundedSurface BoundedSurface
\beginlink \link Go::BoundedSurface BoundedSurface\endlink has entity enumeration 210.
A bounded surface consists of an underlying rectangular parametric surface and one
or more trimming loops. The file format body consists of the following information:
<ul>
<li> The underlying surface represented according to its type
<li> The number of trimming loops
<li> For each trimming loop
<ol>
<li> The number of curves in the current trimming loop
<li> The tolerance within which the loop is found to be continuous
<li> The curves in the trimming loop, one by one
</ol>
</ul>
The trimming curves are of type \beginlink \link Go::ParamCurve ParamCurve\endlink, 
but it most cases they will be represented as \ref g2_sec_curveOnSurface.

The following example illustrates the format a bounded surface. The underlying
surface is a plane as can be seen from the entity number in line two.
The 3 curves in the trimming loop are all of type CurveOnSurface. Note that
the entity number of CurveOnSurface does not appear. The space curve
corresponding to the curve on surface is the master, and the only curve
describing the trimming curve. This is given by the line <em> 0 0 100 </em> 
which occur for every curve and is to be interpreted: <em> the space curve is 
the master, the parameter curve does not exist, the space curve is a 
spline curve </em>. Following this line is the body of a spline curve in the
g2-format.
\verbatim
210 1 0 0
250 
3
0 0 0
0 0 1
1 0 0
1
-3 3
-3 3
0

1
3 0.0001

0 0 100
3 0
2 2
0 0 1 1
-1 0 0
0 -1 0

0 0 100
3 0
2 2
0 0 1 1
0 -1 0
-0.5 2 0

0 0 100
3 0
4 3
0 0 0 0.5 1 1 1
-0.5 2 0
-1 2 0
-2 1 0
-1 0 0
\endverbatim
In this case the two first trimming curves are linear spline curves while
the last one is of order 3 (degree 2) and has one inner knot. All the curves
are non-rational.


\section g2_sec_plane Plane
A \beginlink \link Go::Plane Plane\endlink has entity enumeration 250 and the body
related to the g2-formation contains:
- The dimension of the geometry space
- A point in the plane
- The normal of the plane
- One vector in the set of vectors spanning the plane
- A flag, 0 or 1, for unbounded or bounded
- If the plane is bounded, the four bounds: <em>umin, umax, vmin, vmax</em>
- A flag, 0 or 1, indicating whether or not the \a u and \a v parameter
directions are swapped.

The plane used as an underlying surface in the previous example has
the format:
\verbatim
250 1 0 0
3
0 0 0
0 0 1
1 0 0
1
-3 3
-3 3
0
\endverbatim


\section g2_sec_cylinder Cylinder
A \beginlink \link Go::Cylinder Cylinder\endlink has entity enumeration 260 and the body
related to the g2-formation contains:
- The dimension of the geometry space
- The cylinder radius
- A point on the cylinder axis
- The cylinder axis
- The vector from the cylinder axis towards the cylinder surface giving the
start point of the default parametrization
- A flag, 0 or 1, for unbounded or bounded in the linear direction
- If unbounded, the bounding parameters: <em>umin, umax</em>
- If bounded, the bounding parameters: <em>umin, umax, vmin, vmax</em>
- A flag, 0 or 1, indicating whether or not the \a u and \a v parameter
directions are swapped

\section g2_sec_sphere Sphere
A \beginlink \link Go::Sphere Sphere\endlink has entity enumeration 270 and the body
related to the g2-formation contains:
- The dimension of the geometry space
- The radius of the sphere
- The centre of the sphere
- The axis on which the degeneracies of a NURBS representation of the sphere will 
lie, i.e. the sphere axis
- The vector from the sphere axis towards the sphere surface giving the
start point of the default parametrization
- The bounding parameters: <em>umin, umax, vmin, vmax</em>
- A flag, 0 or 1, indicating whether or not the \a u and \a v parameter
directions are swapped

A sphere with centre in (1.5, 0, 0) and radius 1.5
represented in the g2-format is
\verbatim
270 1 0 0
3
1.5
1.5 0 0
0 0 1
1 0 0
0 6.283185307179586
-1.5707963267948966 1.5707963267948966
0
\endverbatim

\section g2_sec_cone Cone
A \beginlink \link Go::Cone Cone\endlink has entity enumeration 280 and the body
related to the g2-formation contains:
- The dimension of the geometry space
- The cone radius at the position of the point on the cone axis
- A point on the cone axis
- The cone axis
- The vector from the cone axis towards the cone surface giving the
start point of the default parametrization
- The opening angle of the cone
- A flag, 0 or 1, for unbounded or bounded in the linear direction
- If unbounded, the bounding parameters: <em>umin, umax</em>
- If bounded, the bounding parameters: <em>umin, umax, vmin, vmax</em>
- A flag, 0 or 1, indicating whether or not the \a u and \a v parameter
directions are swapped

\section g2_sec_torus Torus
A \beginlink \link Go::Torus Torus\endlink has entity enumeration 290 and the body
related to the g2-formation contains:
- The dimension of the geometry space
- The major radius of the torus
- The minor radius of the torus
- The torus centre
- The mid axis of the torus
- The vector from the torus centre to the torus surface giving the startpoint of
the default parametrization
- Select outer (1=yes, 0=no): Defines the parameter domain of the torus in degenerate
cases.
- A flag, 0 or 1, for unbounded or bounded in the linear direction
- If unbounded, the bounding parameters: <em>umin, umax</em>
- If bounded, the bounding parameters: <em>umin, umax, vmin, vmax</em>
- A flag, 0 or 1, indicating whether or not the \a u and \a v parameter
directions are swapped

\section g2_sec_surfaceOfRevolution SurfaceOfRevolution
A \beginlink \link Go::SurfaceOfRevolution SurfaceOfRevolution\endlink has entity 
enumeration 291 and the body contains:
- The dimension of the geometry space
- A point of the axis of revolution
- The direction of the axis of revolution
- The \em body of the <em> spline curve </em> to rotate around this axis

A surface of revolution is
\verbatim
291 1 0 0
3
0 0 0
0 0 1

3 0
2 2
0 0 1 1
1 0 0
3 0 0
\endverbatim

\section g2_sec_disc Disc
The entity number for \beginlink \link Go::Disc Disc\endlink is 292. The g2 body
contains:
- The dimension of the geometry space
- The centre of the disc
- The disc radius
- The normal to the disc surface
- A vector from the disc centre to the disc boundary to give a start point
for the parametrization
- Whether a NURBS representation should have a degeneracy in the centre (=1) or
have degenerate corners (=0)
- The angles giving the four degeneracy points at the boundary. Not used if the
flag for centre degeneracy is false.

\section g2_sec_pointCloud PointCloud
A \beginlink \link Go::PointCloud PointCloud\endlink has enumeration 400. It has the
following body:
- Number of points
- The coordinates of each point. The points are given one by one.

The dimension of a point in a point cloud is 3 by default. An example of
a point cloud with 3 points is
\verbatim
400 1 0 0
3
0 0 0 
1 0 0 
1 1 0 
\endverbatim

\section g2_sec_lineCloud LineCloud
A \beginlink \link Go::LineCloud LineCloud\endlink has enumeration 410. It has the
following body:
- Number of line segments
- The coordinates of the endpoints of each line segments. 
The lines are given one by one.

The dimension of a line in a line cloud is 3 by default. A line cloud
with 3 lines follow:
\verbatim
410 1 0 0
3
0 0 0  0 0 1
1 0 0  1 0 1
1 1 0  1 1 1
\endverbatim

\section g2_sec_splineVolume SplineVolume
The entity enumeration for a \beginlink \link Go::SplineVolume spline volume \endlink
is 700. The spline volume entity is placed in the module trivariate.
The body of a spline volume is as follows:
- dimension of geometry space, whether or not the volume is rational: 1=rational,
0=non-rational
- the number of coefficients in the first parameter direction, 
the polynomial order in this direction (i.e. degree+1)
- the knot vector in the first parameter direction, 
multiple knots are represented by giving the knot value several times
- the number of coefficients in the second parameter direction, 
the polynomial order in this direction (i.e. degree+1)
- the knot vector in the second parameter direction
- the number of coefficients in the third parameter direction, 
the polynomial order in this direction (i.e. degree+1)
- the knot vector in the third parameter direction
- the volume coefficients

The volume coefficients are given continuously and the 1. parameter direction
runs fastest, the 2. parameter direction is next and the 3. parameter
direction is the most slow. For a non-rational spline volume in the 3-dimensional space, the sequence is
\f$ (x_{(1,1,1)}, y_{(1,1,1)}, z_{(1,1,1)}), (x_{(2,1,1)}, y_{(2,1,1)}, z_{(2,1,1)}),
\ldots, (x_{(n,1,1)}, y_{(n,1,1)}, z_{(n,1,1)}), 
(x_{(1,2,1)}, y_{(1,2),1}, z_{(1,2,1)}),
\ldots, (x_{(n,m,1)}, y_{(n,m,1)}, z_{(n,m,1)}), \ldots, 
\ldots,  (x_{(n,m,p)}, y_{(n,m,p)}, z_{(n,m,p)}) \f$. Here \f$n\f$ is the number of 
coefficients in the 1. parameter direction, \f$m\f$ is the number of coefficients
in the 2. parameter direction and \f$p\f$ is the number of coefficients
in the 3. parameter direction. The weights are multiplied with the coefficients
similarily to the case for the spline curve and the spline surface
and is given as an extra entity for
each coefficient, again similar to the curve and surface case.

A simple non-rational, trilinear volume with no inner knots is represented as:
\verbatim
700 1 0 0 
3 0
2 2
0 0 1 1
2 2
0 0 1 1
2 2
0 0 1 1
0 0 0
1 0 0
0 1 0
1 1 1
0 0 1
1 0 1
0 1 1
1 1 1
\endverbatim

\section g2_sec_parallelepiped Parallelepiped
\beginlink \link Go::Parallelepiped Parallelepiped\endlink has enumeration 720 and
the following g2-format:
- Dimension of geometry space
- Lower right corner
- Unit vector in geometry space representing the first parameter direction
- Length of parallelepiped in the first direction
- Unit vector in geometry space representing the second parameter direction
- Length of parallelepiped in the second direction
- Unit vector in geometry space representing the third parameter direction
- Length of parallelepiped in the third direction

\section g2_sec_sphereVolume SphereVolume
\beginlink \link Go::SphereVolume SphereVolume\endlink has enumeration 721 and
the following g2-format:
- Dimension of geometry space
- Sphere radius
- Centre of sphere
- Axis of sphere, i.e. the vector along which the degeneracies of a NURBS
representation of the sphere will lie
- Vector from the centre to the outer surface giving input to the sphere
parametrization

\section g2_sec_cylinderVolume CylinderVolume
\beginlink \link Go::CylinderVolume CylinderVolume\endlink has enumeration 722 and
the following g2-format:
- Dimension of geometry space
- Cylinder centre
- Cylinder axis
- Vector from the centre to the outer surface giving input to the 
parametrization of the cylinder
- Minium cylinder radius. If this number is larger than 0, the cylinder will
have a hole along the axis
- Maximum radius, the radius of the outer cylinder surface
- Minimum hight, =1: the cylinder is infinite, =0: the cylinder is finite and the
distance from the centre to the bottom of the cylinder folllows (positive or
negative number)
- Maximum hight, =1: the cylinder is infinite, =0: the cylinder is finite and the
total hight of the cylinder follows
- Whether a NURBS representation should have a degeneracy in the centre (=1) or
have degenerate corners (=0)
- The angles giving the four degeneracy points at the boundary. Not used if the
flag for centre degeneracy is false.

\section g2_sec_coneVolume ConeVolume
\beginlink \link Go::ConeVolume ConeVolume\endlink has enumeration 723 and
the following g2-format:
- Dimension of geometry space
- Radius of the cone at the centre
- A point at the cone axis (the centre)
- Cone axis
- The vector from the cone axis towards the boundary surface of the cone giving the
start point of the cone parametrization
- The opening angle of the cone
- Minimum hight the distance from the centre to the smallest end disc of the cone
(positive or negative number)
- Maximum hight, =1: the cone is infinite, =0: the cone is finite and the cone
hight follows
- Whether a NURBS representation should have a degeneracy in the centre (=1) or
have degenerate corners (=0)
- The angles giving the four degeneracy points at the boundary. Not used if the
flag for centre degeneracy is false.

\section g2_sec_torusVolume TorusVolume
\beginlink \link Go::TorusVolume TorusVolume\endlink has enumeration 724 and
the following g2-format:
- Dimension of geometry space
- The torus centre
- The normal to the plane through the big circle of the torus
- A vector from the torus centre to the torus surface giving the parametrization
- The major radius
- The minor radius
- Angle of revolution for the small circle
- Anlge of revolution for the large circle
- Whether a NURBS representation should have a degeneracy in the centre (=1) or
have degenerate corners (=0)
- The angles giving the four degeneracy points at the boundary. Not used if the
flag for centre degeneracy is false.
*/

#endif // _PARAMETRIZATION_DOXYMAIN_H
