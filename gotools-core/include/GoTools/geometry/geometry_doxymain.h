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

#ifndef _GEOMETRY_DOXYMAIN_H
#define _GEOMETRY_DOXYMAIN_H

/**
\page geometry_doc GoTools Core

All GoTools modules depend on gotools-core. The module contains 
parametric curves and surfaces and functionality related to these entities.
For clarity it is divided into four parts: \em geometry, \em utils,
\em creators and \em tesselator.  The \em geometry part contains
classes and functionality for representing, storing and manipulating
parameterized geometrical objects, whereas the
\em utils part contains general, low-level functionality.
\em Creators contains
functionality for generating curves and surfaces by approximation,
blending, etc., and and \em tesselator for making tesselations of the
geometry entities.

\section geom_sec1 Geometric Data Structures
\image html data_structure2.gif "Simplified overview of the geometry class hierarchy"

The figure above shows the main geometry classes in GoTools. 

All geometric objects are of type \beginlink \link Go::GeomObject GeomObject\endlink. 
This class has a function called
instanceType, and by calling this function, it is possible to check the
concrete type of a given object. A GeomObject is 
\beginlink \link Go::Streamable Streamable \endlink, which means
that they can be written to and read from a stream (typically a
file) in a uniform way.  
See \beginlink \link streamable_doc the g2-format documentation\endlink for more
information on this topic. Similarly, they can be created in a
uniform way by the \beginlink \link Go::Factory Factory \endlink, which is
useful when you want to generate objects whose exact kind are
unknown at compile time

\beginlink \link Go::ParamCurve ParamCurve\endlink 
is the base class for all parametric curves in GoTools and defines
a common interface. Many operations do not need to distinguish between 
different types of parametric curves. A parametric curve has functionality
of type:
- operations on the parameter interval; fetch start and end parameter, 
change or reverse the parameter interval
- evaluation
- closest point
- compute the bounding box
- compute the directional cone surronding all tangent directions of
this curve
- length measures of the curve
- fetch a sub curve of this curve
- closest point computatations
- append another curve to the current curve
- check for degeneracy
- make a copy of itself 

\beginlink \link Go::ParamSurface ParamSurface\endlink 
is the base class for a parametric surface and defines the
common interface for all parametric surfaces. The functionality is roughly
- parameter domain operations; fetch the domain and the rectangular 
surrounding domain, check if a parameter pair lies in the domain of a
surface, turn the directions of the parameter domain
- evaluation
- fetch constant parameter curves
- get the boundary curves surrounding this surface
- compute bounding box
- compute the directional cone surrounding all surface normal directions
of the current surface
- area computations
- fetch a sub surface of this surface
- check for degeneracy
- make a copy of itself

\section geom_sec2 B-spline Curves

A B-spline curve is represented by 
\beginlink \link Go::SplineCurve SplineCurve\endlink.

The curve is defined by the formula

\f[ {\bf c}(t) = \sum_{i=1}^{n} {\bf p}_{i} B_{i,k,{\bf t}}(t). \f]

The dimension of the curve \f${\bf c}\f$ is equal to that of its
\em control \em points \f${\bf p}_i\f$. For example, if the dimension of the
control points
is one, the curve is a function, if the dimension is two,
the curve is planar, and if the dimension is three,
the curve is spatial.
Usually the dimension of the curve will be at most three,
but higher dimensions are allowed.

A B-spline curve is a linear combination of a sequence of 
\beginlink \link Go::BsplineBasis B-splines\endlink
\f$B_{i,k,{\bf t}}\f$ (called a B-basis)
uniquely determined by a knot vector \f${\bf t}\f$ and
the order \f$k\f$. Order is equivalent to polynomial degree plus one.
For example, if the order is two, the degree is one and the B-splines
and the curve \f$c\f$ they generate are (piecewise) linear.
If the order is three, the degree is two and the B-splines and
the curve are quadratic. Cubic B-splines and
curves have order 4 and degree 3, etc.

The parameter range of a B-spline curve \f${\bf c}\f$ is the interval
\f$[t_k, t_{n+1}], \f$
and so mathematically, the curve is a mapping
\f${\bf c}: [t_k, t_{n+1}] \to {\bf R}^{dim}\f$, where \f$dim\f$ is the Euclidean space
dimension of its control points.

The complete representation of a B-spline curve consists of
\arg \c \f$dim\f$: The dimension of the underlying Euclidean space,
              \f$1,2,3,\ldots\f$.
\arg \c \f$n\f$: The number of control points (also the number of B-splines)
\arg \c \f$k\f$: The order (polynomial degree + 1) of the B-splines.
\arg \c \f${\bf t}\f$: The knot vector of the B-splines.
            \f${\bf t} = (t_1, t_2, \ldots, t_{n+k})\f$.
\arg \c \f${\bf p}\f$: The control points of the B-spline curve.
           \f$p_{d,i}\;,\; d=1,\ldots,dim\;,\;
                i=1,\ldots,n.\;\;\f$
                e.g. when \f$dim = 3\f$, we have
                \f${\bf p} = (x_1,y_1,z_1,x_2,y_2,z_2,\ldots,x_n,y_n,z_n)\f$.

We note that arrays in \f$c\f$ start at index 0 which means,
for example, that if the array \f$t\f$ holds the knot vector,
then \f$t[0] = t_1,\ldots, t[n+k-1] = t_{n+k}\f$
and the parameter interval goes from
\f$t[k-1]\f$ to \f$t[n]\f$. Similar considerations apply to the other arrays.

The data in the curve representation must satisfy certain conditions:
- The knot vector must be non-decreasing: \f$t_i \le t_{i+1}\f$.
      Moreover, two knots \f$t_i\f$ and \f$t_{i+k}\f$ must be distinct:
      \f$t_i < t_{i+k}\f$.
- The number of control points should be greater than or equal
      to the order of the curve: \f$n \ge k\f$.

To understand the representation better, we will look at
three parts of the
representation: the B-splines (the basis functions),
the knot vector and the control polygon.

\subsection geom_sec2_1 B-spline Basis Functions

The spline space is represented by the class 
\beginlink \link Go::BsplineBasis BsplineBasis\endlink.
A set of B-spline basis functions is determined by the order \f$k\f$
and the knots. For example,
to define a single  basis function of degree one, we need three knots.

\image html linear_B_spline.gif "A linear B-spline (order 2) defined by three knots"
In the figure the three knots are marked as dots.

\image html linear2_B_spline.gif "Linear B-splines of with multiple knots at one end"
In the figure above the knots at the ends are equal.

By taking a linear combination of the three basis functions shown above,
we can generate a linear spline function
as shown in the next figure.

\image html linear_B_spline_curve.gif "A B-spline curve of dimension 1 as a linear combination of a sequence of B-splines.  Each B-spline (dashed) is scaled by a coefficient."
In this figure a linear spline function is generated by taking a linear 
combination of the three basis functions shown in the previous figures.

A quadratic B-spline basis function is a linear combination of two linear
basis functions. In the next figure, we will see
a quadratic B-spline defined by four knots.
A quadratic B-spline is
the sum of two products, the first product between the linear B-spline
on the left and a corresponding line from 0 to 1,
the second product between the linear B-spline
on the right and a corresponding line from 1 to 0;
For higher degree B-splines there is a similar definition.
A B-spline of order \f$k\f$ is the sum of two B-splines of
order \f$k-1\f$, each weighted with weights in the interval [0,1].
In fact we define B-splines of order 1 explicitly as box functions,

\f[
  B_{i,1}(t) =  
  1  \quad if \quad t_i \leq t < t_{i+1} \quad and \quad 0 \quad otherwise
\f]

and then the complete definition of a \f$k\f$-th order B-spline is

\f[ B_{i,k}(t) = {t - t_i \over t_{i+k-1} - t_i} B_{i,k-1}(t)
             +
             {t_{i+k} - t \over t_{i+k} - t_{i+1}} B_{i-1,k-1}(t). \f]

\image html quadratic_B_spline_curve.gif "A quadratic B-spline curve"
The figure above shows two linear B-splines and the corresponding lines 
(dashed) in the quadratic B-spline definition in addition to the quadratic
curve itself.

B-spline basis functions satisfy some important properties for curve 
and surface design.
Each basis function is non-negative and it can be shown that they sum to
one,

\f[ \sum_{i=1}^{n}B_{i,k,{\bf t}}(t) = 1. \f]

These properties combined mean that B-spline curves
satisfy the \em convex \em hull \em property: the curve lies in the convex
hull of its control points.
Furthermore, the support of the B-spline basis function \f$B_{i,k,{\bf t}}\f$ is
the interval \f$[t_i,t_{i+k}]\f$ which means that B-spline
curves has \em local \em control: moving one control point only
alters the curve locally.

With \f$k\f$ multiple knots at the ends of the knot
vector, B-spline curves also have the \em endpoint \em property:
the start point of the B-spline curve
equals the first control point and the end point equals the
last control point, in other words

\f$ {\bf c}(t_k) = {\bf p}_1
     \qquad \hbox{and} \qquad
   {\bf c}(t_{n+1}) = {\bf p}_n. \f$

\subsection geom_sec2_2 The Control Polygon

The \em control \em polygon of a B-spline curve is the polygonal
arc formed by its control points,
\f${\bf p}_0, {\bf p}_1, \ldots,{\bf p}_n\f$.
This means that
the control polygon, regarded as a parametric curve,
is itself a piecewise linear B-spline curve (order two).
If we increase the order, the distance between the control polygon
and the curve increases.
A higher order B-spline curve tends to smooth
the control polygon and at the same time mimic its shape.
For example, if the control polygon is convex, so is the B-spline curve.

Another property of the control polygon is that it will get closer
to the curve if it is redefined by inserting knots into the curve
and thereby increasing the number of control points. This can be seen in the next
figure.
If the refinement is infinite then the control polygon converges to the curve.
\image html B_spline_curves.gif "Linear, quadratic, and cubic B-spline curves sharing the same control polygon"

The control polygon in the figure above is equal to the linear B-spline curve. 
The curves are planar, i.e. the space dimension is two.

\image html cubic_B_spline_curve.gif "The cubic B-spline curve with a redefined knot vector"

\subsection geom_sec2_3 The Knot Vector

The knots of a B-spline curve describe the following properties of the curve:
- The parameterization of the B-spline curve
- The continuity at the joins between the adjacent polynomial
   segments of the B-spline curve.
In figure~\ref{curve7} we have two curves
with the same control polygon and order but
with different parameterization.

This example is not meant as an encouragement to use
parameterization for modelling, rather to make users
aware of the effect of parameterization. Something close to
curve length parameterization is in most cases
preferable. For interpolation, chord-length parameterization
is used in most cases.

The number of equal knots determines the degree
of continuity. If \f$k\f$ consecutive internal knots are equal,
the curve is discontinuous.
Similarly if \f$k-1\f$ consecutive internal knots are equal,
the curve is continuous but not in general differentiable.
A continuously differentiable curve
with a discontinuity in the second derivative
can be modelled using \f$k-2\f$ equal knots etc.
Normally, B-spline curves in GoTools are expected to be continuous.
For many algorithms, curves should
be continuously differentiable (\f$C^1\f$).

\image html two_curves.gif "Two quadratic B-spline curves" 
The curves in the figure above have the same control polygon but different knot 
vectors. The curves and the control polygons are two-dimensional

\image html quadratic_B_spline_curve2.gif "A quadratic B-spline curve with two inner knots"

\subsection geom_sec2_4 NURBS Curves

A NURBS (Non-Uniform Rational B-Spline) curve is a generalization
of a B-spline curve,

\f[ {\bf c}(t) = {\sum_{i=1}^{n} w_i {\bf p}_{i} B_{i,k,{\bf t}}(t)
                 \over
                 \sum_{i=1}^{n} w_i B_{i,k,{\bf t}}(t)} . \f]

In addition to the data of a B-spline curve, the NURBS curve
\f${\bf c}\f$ has a sequence of weights \f$w_1,\ldots,w_n\f$.
One of the advantages of NURBS curves over B-spline curves is that
they can be used to represent conic sections exactly (taking the
order \f$k\f$ to be three).
A disadvantage is that NURBS curves depend nonlinearly on their weights,
making some calculations, like the evaluation of derivatives,
more complicated and less efficient than with B-spline curves.

The representation of a NURBS curve is the same as for a B-spline
except that it also includes
\arg \c \f${\bf w}\f$: A sequence of weights
            \f${\bf w} = (w_1, w_2, \ldots, w_n)\f$.

We make the assumption that the weights are strictly positive: \f$w_i > 0\f$.

Under this condition, a NURBS curve, like its B-spline cousin,
enjoys the convex hull property. 
With \f$k\f$-fold knots at the ends of the knot vector,
also NURBS curves have the endpoint interpolation property.

A NURBS curve is in GoTools stored using the same entities as
non-rational spline curves. Note that the constructors of these entities assume that
the NURBS coefficients are given in the format: \f$w_i p_{i,1}, \ldots
w_i p_{i,dim},w_i\f$ for \f$i=1, \ldots n\f$, 
i.e. the coefficients are multiplied with the weigths.

\subsection geom_sec2_5 Spline Curve Functionality
In GoTools, \beginlink \link Go::SplineCurve SplineCurve \endlink 
inherites \beginlink \link Go::ParamCurve ParamCurve \endlink 
and has thereby the functionality
specified for ParamCurve objects. Functionality specific for a B-spline
curve can be found in the doxygen information. This functionality includes:
- Compute the derivative curve corresponding to a given curve
- Fetch information related to the spline space
- Fetch the control polygon of a spline curve
- Insert new knots into the knot vector of the curve and update the
curve accordingly
- Increase the polynomial degree of the curve
- Make sure that the curve has got an open knot vector, i.e. knot
multiplicity equal to the order in the ends of the curve

\section geom_sec3 B-spline Surfaces

A tensor product B-spline surface is represented by the class 
\beginlink \link Go::SplineSurface SplineSurface\endlink.

the B-spline surface is defined as

\f[ {\bf s}(u,v) = \sum_{i=1}^{n_1}\sum_{j=1}^{n_2}{\bf p}_{i,j} 
	B_{i,k_1,{\bf u}}(u) B_{j,k_2,{\bf v}}(v)  \f]

with control points \f${\bf p}_{i,j}\f$ and two variables
(or parameters) \f$u\f$ and \f$v\f$.
The formula shows that a basis function of a B-spline surface is a
product of two basis functions of B-spline curves (B-splines).
This is why a B-spline surface is called a tensor-product surface.
The following is a list of the components of the representation:
\arg \c \f$dim\f$: The dimension of the underlying Euclidean space.
\arg \c \f$n_1\f$: The number of control points with respect to the first parameter.
\arg \c \f$n_2\f$: The number of control points with respect to the second parameter.
\arg \c \f$k_1\f$: The order (polynomial degree + 1) of the B-splines in the first parameter.
\arg \c \f$k_2\f$: The order of the B-splines in the second parameter.
\arg \c \f${\bf u}\f$: The knot vector of the B-splines with respect to
                  the first parameter,
                  \f${\bf u} = (u_1,u_2,\ldots,u_{n_1+k_1})\f$.
\arg \c \f${\bf v}\f$: The knot vector of the B-splines with respect to
                  the second parameter,
                  \f${\bf v} = (v_1,v_2,\ldots,v_{n_2+k_2})\f$.
\arg \c \f${\bf p}\f$: The control points of the B-spline surface,
           \f$c_{d,i,j}\f$, \f$d=1,\ldots,dim\f$, \f$i=1,\ldots,n_1\f$,
		\f$j=1,\ldots,n_2\f$.
	When \f$dim = 3\f$, we have
          \f${\bf p} = (x_{1,1},y_{1,1},z_{1,1},x_{2,1},y_{2,1},z_{2,1},\ldots\f$,
                  \f$x_{n_1,1},y_{n_1,1},z_{n_1,1},\ldots\f$,
                     \f$x_{n_1,n_2},y_{n_1,n_2},z_{n_1,n_2})\f$.


The data of the B-spline surface must fulfill the following requirements:
- Both knot vectors must be non-decreasing.
- The number of control points must be greater than or equal to the order
with respect to both parameters: \f$n_1 \ge k_1\f$ and \f$n_2 \ge k_2\f$.

The properties of the representation of a B-spline surface are
similar to the properties of the representation of a B-spline curve.
The control points \f${\bf p}_{i,j}\f$ form a \em control \em net as shown in
the figure below.
The control net has similar properties to the control
polygon of a B-spline curve, described in section~\ref{contrlpoly}.
A B-spline surface has two knot vectors, one

\image html surf1.gif "A B-spline surface and its control net" width=10cm

\subsection  geom_sec3_2 The Basis Functions

A basis function of a B-spline surface is the product of two
basis functions corresponding to B-spline curves,

\f$ B_{i,k_1,{\bf u}}(u) B_{j,k_2,{\bf v}}(v). \f$

Its support is the rectangle \f$[u_i,u_{i+k_1}] \times [v_j,v_{j+k_2}]\f$.
If the basis functions in both directions are of degree one and all
knots have multiplicity one, then the surface basis functions are
pyramid-shaped.
For higher degrees, the surface basis functions are bell shaped.

\image html surface_basis_function.gif "A basis function of degree one in	both variables"
\subsection geom_sec3_3 NURBS Surfaces

A NURBS (Non-Uniform Rational B-Spline) surface is a generalization
of a B-spline surface,

\f[ {\bf s}(u,v) = {\sum_{i=1}^{n_1}\sum_{j=1}^{n_2} w_{i,j} {\bf p}_{i,j} 
	    B_{i,k_1,{\bf u}}(u) B_{j,k_2,{\bf v}}(v)  \over
             \sum_{i=1}^{n_1}\sum_{j=1}^{n_2} w_{i,j}
	         B_{i,k_1,{\bf u}}(u) B_{j,k_2,{\bf v}}(v)}. \f]

In addition to the data of a B-spline surface, the NURBS surface
has a weights \f$w_{i,j}\f$.
NURBS surfaces can be used to exactly represent several common
`analytic' surfaces such as spheres, cylinders, tori, and cones.
A disadvantage is that NURBS surfaces depend nonlinearly on their weights,
making some calculations less efficient.

The representation of a NURBS surface is the same as for a B-spline surface
except that it also includes

\arg \c \f${\bf w}\f$: The weights of the NURBS surface,
           \f$w_{i,j}\f$, \f$i=1,\ldots,n_1\f$, \f$j=1,\ldots,n_2\f$, so
          \f${\bf w} = (w_{1,1},w_{2,1},\ldots,w_{n_1,1},\ldots\f$,
                  \f$w_{1,2},\ldots,w_{n_1,n_2})\f$.

The weights are required to be strictly positive: \f$w_{i,j} > 0\f$.

The NURBS surface is represented by SISLSurf and SplineSurface in SISL and
GoTools, respectively. As for the curve case, the constructors of
SISLSurf and SplineSurface expect the
coefficients to be multiplied with the weights.

\subsection geom_sec3_4 Spline Surface Functionality
\beginlink \link Go::SplineSurface SplineSurface\endlink inherites 
\beginlink \link Go::ParamSurface ParamSurface \endlink 
in the GoTools data structure and 
implements the functionality described for ParamSurface. Important additional
functionality is listed below: 
- Compute the derivative surface and normal surface
 corresponding to a given surface
- Fetch information related to the spline spaces
- Fetch the control polygon of a spline surface
- Fetch all weights of a NURBS surface
- Grid evaluation of the basis functions related to a surface
- Grid evaluation of the surface
- Insert new knots into the knot vectors of the surface and adapt the
surface description accordingly
- Increase the polynomial degrees of the surface
- Fetch information related to the boundary curves of a surface

\section geom_sec4 Other Curve Types
As shown in the class hierarchy, a spline curve is not the only
parametric curve in GoTools although it is the most important one. There are
also other curves that share parts of the public interface:
\arg \c \beginlink \link Go::BoundedCurve BoundedCurve\endlink 
A bounded curve is a restriction of a parametric curve
to a given interval. The interval can be specified either in parameter
space, in geometry space or both. In the last case, it must be specified 
whether the parameter space bound or the geometry space bound is the
master. 
\arg \c \beginlink \link Go::ElementaryCurve ElementaryCurve\endlink 
An elementary curve is a container for a conic section 
or a line. This entity will be described in some detail later.
\arg \c \beginlink \link Go::CurveOnSurface CurveOnSurface\endlink 
This class represents a curve lying on a surface.

\subsection geom_sec4_2 Curve On Surface
The \beginlink \link Go::CurveOnSurface CurveOnSurface\endlink 
entity is often related to a bounded, or trimmed, surface. A
bounded surface is limited by a curve loop where each curve in the loop most
often is represented by a CurveOnSurface. However, a CurveOnSurface is an
entity in its own right. It simply represents a curve lying on a surface. It
can for instance originate from intersecting the surface with a plane.

A CurveOnSurface consists of a parametric surface and two 
corresponding curves, one curve in the parameter domain of the surface and one
spatial curve. The master representation is specified. Some operations may
required the existence of one particular of these curves, 
in particular the parameter curve
is requested. It exists functions that compute the missing curve
from the existing one.

\subsection geom_sec4_3 Elementary Curve
An \beginlink \link Go::ElementaryCurve elementary curve\endlink 
is a fairly new concept in GoTools, and until now 
such curves have
been entered as entities in an IGES or STEP file. The elementary curves 
may also
be represented as spline curves although non-bounded curves must be given a
finite extension. The elementary curves is of the types:
- \beginlink \link Go::Circle Circle\endlink 
The circle is defined by its center and radius and the plane
in which it lies if the dimension of the geometry space is larger than 2.
A circle is parameterized in terms of the angle. This is a bounded
parameterization.
- \beginlink \link Go::Ellipse Ellipse\endlink 
An ellips is represented by a centre, the direction of one
semi-axis and the length of the two semi-axis.  The plane
in which the ellipse lies is required
if the dimension of the geometry space is larger 
than 2. An ellipse is parameterized in terms of the angle which gives a 
bounded parameterization.
- \beginlink \link Go::Line Line\endlink 
A line is given by a point and a direction. It has an unbounded
parameterization, \f$t \in [- \infty , \infty ]\f$.
- \beginlink \link Go::Parabola Parabola\endlink 
A parabola is given by a position, \f$C\f$, a focal 
distance, \f$F\f$, and a coordinate system, \f$({\bf x}, {\bf y}, {\bf z})\f$. 
It is described as 
\f$f(t) = C + F(t^2 {\bf x} + 2t{\bf y})\f$ and has an unbounded parameterization,
\f$t \in [- \infty , \infty ]\f$.
- \beginlink \link Go::Hyperbola Hyperbola\endlink A hyperbola is given by a
position, a coordinate system and two distances. The parameterization is
unbounded.

\section geom_sec5 Other Surface Types
Similar to the curve case, GoTools supports also other types of parametric
surfaces than spline surfaces.
The additional surface types are bounded surface, elementary surface and 
composite surface. The elementary surfaces are conic surfaces, plane and 
torus. A composite surface consists of a number of spline surfaces where the
associated parameter domains are mapped into one composite domain. The class
is not complete in the sense that it does 
not support all functions specified in the
interface for a parametric surface.

\subsection geom_sec5_2 Bounded Surface
\beginlink \link Go::BoundedSurface A bounded surface\endlink 
is a trimmed surface defined by an underlying surface and a
trimming loop. The underlying surface is a parametric surface and the loop
consists of a closed sequence of ordered parametric curves, often of type
CurveOnSurface. The trimming loop must lie on the surface and be continuous
with respect to a given tolerance.

\subsection geom_sec5_3 Elementary Surface
Also \beginlink \link Go::ElementarySurface elementary surfaces\endlink 
are quite new in GoTools, and have
been entered as entities from IGES or STEP. The elementary surfaces 
may also
be represented as spline surfaces although non-bounded surfaces must be given a
finite extension. Operations 
involving this curve may, however, be made more efficiently by the knowledge
of the surface type.
The elementary surfaces is of the types:
- \beginlink \link Go::Cone Cone\endlink 
A cone is described by a location, the cone axis, the radius 
of the cross section corresponding to the location and the cone angle.
The parameterization is given by the angle and the distance along the axis,
thus it is unbounded in one parameter direction.
- \beginlink \link Go::Cylinder Cylinder\endlink 
A cylinder is given by its centre and radius and the
cylinder axis. It is parameterized by the angle around the circle and the
distance in the axis direction.
- \beginlink \link Go::Plane Plane\endlink 
A plane is given by a point and a normal. It has an unbounded
parameterization in both parameter directions.
- \beginlink \link Go::Sphere Sphere\endlink 
A sphere is given by its centre and radius. It has an angular
parameterization in both parameter directions and thereby bounded.
- \begnlink \link Go::Torus Torus\endlink 
A torus is given by the minor and major radius, a location and
an axis. The parameterization is angular in both parameter directions.
- \beginlink \link Go::Disc Disc \endlink
A disc is a part of a plane limited by a circle and given by its centre,
radius, plane normal and an additional vector to generate a coordinate 
system.

The elementary surfaces described above is placed in a coordinate system to
facilitate the parameterization. In addition to the geometric information
given for each entity, remaining coordinate system information must be 
specified.

\section geom_sec6 Geometry Construction
The combination of SISL and GoTools provide a lot of methods for geometry
construction. They partly overlap. The GoTools functionality is outlined
here. More details can be found in the appropriate classes and namespaces.
The SISL geometry construction methods are described in the SISL manual. 

The classes for parametric curves, surfaces and volumes in GoTools do not
contain construction facilites. They are solely for operating on these entites.
The construction is performed in separate classes.

\subsection geom_sec6_2 Curve Construction
Construction of elementary curves are being done by their definition or they are
read from an IGES file. The curve construction methods in GoTools are
concerned with spline curves. The methods are split between two sub modules
in gotools-core, namely geometry and creators. The relevant classes are
<ul>
<li> \beginlink \link Go::HermiteInterpolator HermiteInterpolator\endlink 
The class is placed in geometry. Interpolates a
set of points where each point is associated a tangent vector. The points
must be parameterized. This is a local curve construction method.
<li> \beginlink \link Go::SplineInterpolator SplineInterpolator\endlink 
A set of methods to interpolate a set of points. 
Some methods accept tangent vectors associated to some points. The 
points must be parameterized. The class gives the possibility to set 
particular conditions in the endpoints of the curve. This is a global method
for curve construction.
<li> \beginlink \link Go::SplineApproximator SplineApproximator\endlink 
Approximation in the least squares sense. 
Parameterized points and a spline space must be given. This class and the
two proceeding ones inherite the same class, 
\beginlink link Go::Interpolator Interpolator\endlink, and share parts
of the interface.
<li> \beginlink \link Go::ApproxCurve ApproxCurve\endlink 
Make a curve approximating a set of points with a given
accuracy. The points must be parameterized and it is possible to give an
initial spline space. The method iteratively applies the method in the class
SmoothCurve. This class and the following ones lie in the sub module
creators.
<li> \beginlink \link Go::SmoothCurve SmoothCurve\endlink 
Approximate a set of parameterized data points with a 
spline curve using a least squares approximation with a smoothing term. The
spline space must be given. The method can also be used to perform smoothing
on a given curve.
<li> \beginlink \link Go::AdaptCurve AdaptCurve\endlink 
Approximates an evaluator based curve with a given accuracy. The
method is based on least squares approxation with a smoothing term and spline
space refinement.
<li> \beginlink \link Go::ApproxCrvToSeqs ApproxCrvToSeqs\endlink 
Approximate a set of curve sequences by a set of B-spline
curves within a given accuracy. The curves are represented in the same spline space.
<li> \beginlink \link Go::HermiteAppC HermiteAppC\endlink 
Approximate an evaluator based curve represented by a
sub class of the class \beginlink \link Go::EvalCurve EvalCurve\endlink, 
by a spline curve. An hermite method is 
used, thus a \f$C^1\f$ continuous curve will be generated. The existing evaluator
based curves are
<ol>
      <li> \beginlink \link Go::ProjectCurve A curve in the parameter domain \endlink
of a surface is represented
as the projection of a 3D curve into the parameter domain of a given surface. The
output curve lies in the parameter domain of the surface.
      <li> Given data describing a 
\beginlink \link Go::ProjectIntersectionCurve surface-surface intersection curve\endlink, 
the curve
lying given offset distances from the two surface involved in the intersection
and parameterized similar to the intersection curve is represented. This
construction can be used to create a rolling ball blend.
<li> Given two approximations of the same intersection curves, reapproximate the
\beginlink \link Go::SpaceIntCrv intersection curve\endlink.
<li> Given a \beginlink \link Go::LiftCurve curve in the parameter domain \endlink
of a surface, the geometry space curve is evaluated.
<li> Represent \beginlink \link Go::CrossTangentOffset an offset curve\endlink 
from a given space curve. The direction of
the offset is given as an expression of corresponding curves. 
<li> An \beginlink \link Go::CrossTanOffDist offset curve \endlink to a 
given curve in the direction of a cross tangent curve constructed as a
blend between tangent curves.
<li> A \beginlink \link Go::EvalParamCurve parameter curve \endlink of any
type that is to be approximated by a spline curve.
</ol>
<li> \beginlink \link Go::HermiteAppS HermiteAppS\endlink 
Approximate an evaluator based curve set represented by a
sub class of the class \beginlink \link Go::EvalCurveSet EvalCurveSet\endlink. 
A number of curves defined in the same
spline space are generated. One concreate evaluator based curve set exists
currently:
<ol>
<li> Given a surface, a curve in geometry space and a corresponding cross
tangent curve, 
\beginlink \link Go::ProjectCurveAndCrossTan the parameter curve \endlink
representing the projection of the 
geometry curve into the surface and the corresponding projection of the cross
tangent in the parameter domain is represented.
<li> Given an \beginlink \link Go::IntCrvEvaluator intersection curve\endlink
represented as two CurveOnSurface instances,
reapproximate this intersection curve along with the two parameter domain
curves representing this curve.
<li> Improve the accuracy of an existing \beginlink \link Go::TrimCurve trimming
curve\endlink.
<li> Given an intersection curve between two surfaces, offset these surfaces
given distances and create the 
\beginlink \link Go::SmoothTransition intersections curve \endlink 
between the offset
surfaces and corresponding cross tangent curves that will give a smooth
transition between the offset surfaces.
</ol>
</ul>
\subsection geom_sec7 Surface Construction
Similar to the curve construction methods, these methods are placed in
the sub modules geometry and creators.
<ul>
<li> \beginlink \link Go::SweepSurfaceCreator SweepSurfaceCreator\endlink 
The class lies
in geometry. The class contains two methods; sweep a curve along another
curve to create a surface and rotational sweep.
<li> \beginlink \link Go::ApproxSurf ApproxSurf\endlink 
This class and the following lie in creators. A set of
parameterized scattered data are approximated by a spline surface. The boundary
curves to the surface can be defined a priori. The spline space must be
given. The method iteratively applies the method in SmoothSurface and performs
parameter iteration. It may also refine the spline space.
<li> \beginlink \link Go::SmoothSurf SmoothSurf\endlink 
Approximate a set of parameterized scattered data points
using a least squares approach with a smoothing term. The spline space must
be given. The method can also be used to smooth a given surface.
<li> \beginlink \link Go::LoftSurfaceCreator LoftSurfaceCreator \endlink
constructs a surface by interpolating a number of non-rational spline
curves.
<li> \beginlink \link Go::CoonsPatchGen CoonsPatchGen\endlink 
Construct a surface interpolating 4 boundary curves. The
curves may be attached cross tangent curves. To achieve interpolation, the
boundary conditions must be consistent.
<li> \beginlink \link Go::HahnsSurfaceGen HahnsSurfaceGen\endlink 
Fill an \f$n\f$-sided hole with \f$n\f$ surfaces. \f$n=3\f$, \f$5\f$ or
\f$6\f$.
</ul>

\section geom_sec8 The utils module
The 'utility' module was developed in parallel with the 'geometry' module,
and contains various useful classes and functions related to math. / 
computational geometry, while not being directly related to splines.

\subsection geom_sec8_1 Data structures
You can find templated arrays and optimised vectors.

\subsection geom_sec8_2 Geometrical objects
The module contain some objects that are much used in computational geometry,
like
<ul>
<li> \link Go::Point points \endlink </li>
<li> \link Go::BoundingBox bounding boxes \endlink, axis-aligned or 
     \link Go::RotatedBox rotated \endlink, eventually 
     \link Go::CompositeBox composite \endlink</li>
<li> \link Go::Rational rational \endlink numbers </li>
<li> \link Go::BaryCoordSystem barycentric \endlink coordinate systems </li>
<li> \link Go::DirectionCone direction cones \endlink </li>
</ul>
\subsection geom_sec8_3 Misc.
<ul>
<li> Compile-time computation of \link Go::Factorial factorials \endlink </li>
<li> A general, nonlinear \link Go::FunctionMinimizer function minimizer \endlink  
     in an arbitrary number of variables.  This can be useful in a great number
     of different settings! </li>
</ul>

\section geom_sec9 Tesselation
The geometry entities are tesselated in gotools-core in the sub module
tesselator. Depending on the type of the geometry object, the result of the
tesselation is stored in \beginlink \link Go::GenericTriMesh GenericTriMesh\endlink, 
\beginlink \link Go::LineStrip LineStrip\endlink, 
\beginlink \link Go::QuadMesh QuadMesh\endlink or
\beginlink \link Go::RegularMesh RegularMesh\endlink 
which all inherite the abstract class 
\beginlink \link Go::GeneralMesh GeneralMesh\endlink. The tesselators
in this module relate to a given resolution in each parameter direction of
the geometry object.

A curve is tesselated by the class 
\beginlink \link Go::CurveTesselator CurveTesselator\endlink 
and the tesselated curve
is returned by the function CurveTesselator::getMesh() as a shared pointer 
to a LineStrip object. 
A rectangular surface is tesselated in 
\beginlink \link Go::RectGridTesselator RectGridTesselator\endlink and the tesselator
returns a shared pointer to a RegularMesh. 
The triangles produced during tesselation is organized in a set of triangle
strips. This structure is taken advantage of to improve the drawing
efficiency.
A trimmed surface is tesselated by 
\beginlink \link Go::ParametricSurfaceTesselator ParametricSurfaceTesselator\endlink 
and the tesselation is represented as an unorganized bunch of triangles. 
The output from the functionality in the tesselation module is taken as
input to the GoTools viewers and the actual drawing is performed using
openGL.

 */

#endif // _GEOMETRY_DOXYMAIN_H
