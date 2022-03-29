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

#ifndef _TRIVARIATE_DOXYMAIN_H
#define _TRIVARIATE_DOXYMAIN_H

/**
\page trivariate GoTools Trivariate 

This module trivariate represents NURBS volumes and contains construction methods
and operations related to such volumes.

This module depends on:
- GoTools Core library
- SISL library

\section sec1 Data Structures

\image html data_structure.gif "Simplified overview of the geometry class hierarchy" 
\image latex data_structure.eps "Simplified overview of the geometry class hierarchy" 

The figure above shows the main geometric classes in GoTools
and how they are divided between modules. 

\section sec2 B-spline Volumes
A B-spline volume is represented in  \link Go::SplineVolume SplineVolume
\endlink in the GoTools module
trivariate.

The volume is defined by the formula

\f${\bf V}(u,v,w) = \sum_{i=1}^{n_1}\sum_{j=1}^{n_2}\sum_{h=1}^{n_3}{\bf p}_{i,j,h} 
	B_{i,k_1,{\bf u}}(u) B_{j,k_2,{\bf v}}(v) B_{h,k_3,{\bf w}}(w) \f$

with control points \f${\bf p}_{i,j,h}\f$ and three variables
(or parameters) \f$u\f$, \f$v\f$ and \f$w\f$. A basis function of a B-spline volume is a
product of three basis functions of B-spline curves (B-splines).

The following is a list of the components of the representation:

\arg \c \f$dim\f$: The dimension of the underlying Euclidean space.
\arg \c \f$n_1\f$: The number of control points with respect to the first parameter.
\arg \c \f$n_2\f$: The number of control points with respect to the second parameter.
\arg \c \f$n_3\f$: The number of control points with respect to the third parameter.
\arg \c \f$k_1\f$: The order (polynomial degree + 1) 
of the B-splines in the first parameter.
\arg \c \f$k_2\f$: The order of the B-splines in the second parameter.
\arg \c \f$k_3\f$: The order of the B-splines in the third parameter.
\arg \c \f${\bf u}\f$: The knot vector of the B-splines with respect to
                  the first parameter,
                  \f${\bf u} = (u_1,u_2,\ldots,u_{n_1+k_1})\f$.
\arg \c \f${\bf v}\f$: The knot vector of the B-splines with respect to
                  the second parameter,
                  \f${\bf v} = (v_1,v_2,\ldots,v_{n_2+k_2})\f$.
\arg \c \f${\bf w}\f$: The knot vector of the B-splines with respect to
                  the third parameter,
                  \f${\bf w} = (w_1,w_2,\ldots,w_{n_3+k_3})\f$.
\arg \c \f${\bf p}\f$: The control points of the B-spline volume,
           \f$c_{d,i,j,h}\f$, \f$d=1,\ldots,dim\f$, \f$i=1,\ldots,n_1\f$,
		\f$j=1,\ldots,n_2\f$, \f$h=1,\ldots,n_3\f$.
  When \f$dim = 3\f$, we have
  \f${\bf p} = (x_{1,1,1},y_{1,1,1},z_{1,1,1},x_{2,1,1},y_{2,1,1},z_{2,1,1},
  \ldots\f$,
  \f$x_{n_1,1,1},y_{n_1,1,1},z_{n_1,1,1},\ldots\f$,
  \f$x_{n_1,n_2,1},y_{n_1,n_2,1},z_{n_1,n_2,1},\ldots
  x_{n_1,n_2,n_3},y_{n_1,n_2,n_3},z_{n_1,n_2,n_3})\f$.

The data of the B-spline volume must fulfill the following requirements:

\li All knot vectors must be non-decreasing.
\li The number of control points must be greater than or equal to the order
with respect to all three parameters: \f$n_1 \ge k_1\f$, \f$n_2 \ge k_2\f$
and \f$n_3 \ge k_3\f$.

The properties of the representation of a B-spline volume are
similar to the properties of the representation of a B-spline curve or surface.
The control points \f${\bf p}_{i,j,h}\f$ form a \em control \em net.
The control net has similar properties to the control
polygon of a B-spline curve, described in the module gotools-core.
A B-spline volume has three knot vectors, one
for each parameter. 

\subsection sec2_1 The Basis Functions
A  \link Go::BsplineBasis basis function \endlink 
of a B-spline volume is the product of three
basis functions corresponding to B-spline curves,

\f$  B_{i,k_1,{\bf u}}(u) B_{j,k_2,{\bf v}}(v) B_{h,k_3,{\bf w}}(w)\f$ 

Its support is the box \f$[u_i,u_{i+k_1}] \times [v_j,v_{j+k_2}]
\times [w_h,w_{h+k_3}]\f$.


\subsection sec2_2 NURBS Volumes
A NURBS (Non-Uniform Rational B-Spline) volume is a generalization
of a B-spline volume,

\f$ {\bf V}(u,v,w) = {\sum_{i=1}^{n_1}\sum_{j=1}^{n_2} \sum_{r=1}^{n_3} 
h_{i,j,r} {\bf p}_{i,j,r} 
	    B_{i,k_1,{\bf u}}(u) B_{j,k_2,{\bf v}}(v)  B_{r,k_3,{\bf w}}(w)  \over
             \sum_{i=1}^{n_1}\sum_{j=1}^{n_2} \sum_{r=1}^{n_3} h_{i,j,r}
	         B_{i,k_1,{\bf u}}(u) B_{j,k_2,{\bf v}}(v)B_{r,k_3,{\bf w}}(w)}. \f$

In addition to the data of a B-spline surface, the NURBS surface
has a weights \f$h_{i,j,r}\f$.
NURBS volumes can be used to exactly represent volumes that have common
`analytic' surfaces such as spheres, cylinders, tori, and cones as
boundary surfaces.
A disadvantage is that NURBS volume depend nonlinearly on their weights,
making some calculations less efficient.

The representation of a NURBS volume is the same as for a B-spline volume
except that it also includes

\f${\bf h}\f$: The weights of the NURBS volume,
           \f$h_{i,j,r}\f$, \f$i=1,\ldots,n_1\f$, \f$j=1,\ldots,n_2\f$, \f$r=1,\ldots,n_3\f$, so
          \f${\bf h} = (h_{1,1,1},h_{2,1,1},\ldots,h_{n_1,1,1},
  h_{1,2,1},\ldots,h_{n_1,n_2,1}, \ldots h_{n_1,n_2,n_3})\f$.

The weights are required to be strictly positive: \f$h_{i,j,r} > 0\f$.

The NURBS volume is represented by SplineVolume. As for the curve and surface
cases, the constructor expects the
coefficients to be multiplied with the weights.

\subsection sec2_3 Spline Volume Functionality
The functionality of a  \link Go::SplineVolume spline volume \endlink 
to a large extend corresponds to the
functionality of a spline surface. Important functionality is:
- A NURBS volume is able to make a copy of itself
- Compute the bounding box of the volume
- Evaluation and grid evaluation
- Grid evaluation of basis functions
- Compute the derivative volume corresponding to a volume
- Closest point computation
- Fetch a sub volume of a given volume
- Fetch information related to the spline spaces
- Swap and reverse parameter directions in a volume
- Fetch the control polygon of the volume
- Fetch all weights of a NURBS volume
- Insert knots into the spline spaces of the volume and adapt the volume
description accordingly
- Increase the polynomial degree of the volume in one parameter direction
- Fetch a constant parameter surface from the volume
- Fetch all boundary surfaces surrounding a volume
- Check for periodicity and degeneracy

\subsection sec2_4 Construction Methods for SplineVolume
The following methods exist for construction of a spline volume. The
corresponding GoTools class names are given in brackets.
- Sweep a NURBS surface along a NURBS curve 
( \link Go::SweepVolumeCreator SweepVolumeCreator \endlink). 
- Rotational sweep of a NURBS surface 
( \link Go::SweepVolumeCreator SweepVolumeCreator \endlink).
- Lofting to interpolate a number of NURBS surfaces 
( \link Go::LoftVolumeCreator LoftVolumeCreator \endlink).
- Interpolate 6 boundary surface to create a volume using a Coons patch
approach ( \link Go::CoonsPatchVolumeGen CoonsPatchVolumeGen \endlink). 
This functionality applies only to non-rational spline surfaces.
- Represent an  \link Go::ElementaryVolume elementary volume \endlink 
as a spline volume. An elementary volume is a 
 \link Go::ParamVolume parametric volume \endlink similarily to
spline volumes. The elementary volumes are:
   -#  \link Go::SphereVolume Sphere  \endlink
   -#  \link Go::CylinderVolume Cylinder  \endlink
   -#  \link Go::ConeVolume Cone  \endlink
   -#  \link Go::Parallelepiped Parallelepiped  \endlink
   -#  \link Go::TorusVolume Torus  \endlink

A spline volume may have a well behaved outer boundary, but a bad
distribution of coefficients in the inner. This is in particular the
case if the volume is constructed by a Coons patch approach. The positioning
of the internal coefficients may be improved by smoothing. The coefficients
at the boundaries are kept fixed and the coefficients in the inner are
redistributed by solving a minimization problem. The smoothing is performed
in the class  \link Go::SmoothVolume SmoothVolume \endlink.

\section sec3 Surface on Volume
The module trivariate provides a surface which extends the class of parametric 
surfaces defined in gotools-core, namely  \link Go::SurfaceOnVolume
SurfaceOnVolume \endlink. This surface inherits most of the functionality defined
for parametric surfaces and takes the same role as CurveOnSurface for
parametric curves. The surface possesses information about
- The associated volume
- The geometric description of this surface and/or
- The description of this surface in the parameter domain of the given volume
- Constant parameter and volume boundary information. If this surface for
instance happens to 
be a constant parameter surface in the given volume, it knows the parameter
direction and the associated constant parameter value.

\section sec4 Curves Related to a Parametric Volume
Two evaluator based curves can be used in together with 
 \link Go::HermiteAppC HermiteAppC \endlink in gotools-core/creators
exist, namely 
 \link Go::VolumeParameterCurve VolumeParameterCurve \endlink and
 \link Go::VolumeSpaceCurve VolumeSpaceCurve \endlink. 
*/

#endif // _TRIVARIATE_DOXYMAIN_H
