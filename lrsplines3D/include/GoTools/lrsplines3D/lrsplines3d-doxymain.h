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

#ifndef _LRSPLINE3D_DOXYGEN_H
#define _LRSPLINE3D_DOXYGEN_H

/**
 \page lrsplines3d GoTools LRsplines3D

 This module contains classes and tools for working with LR spline volumes.


The module lrsplines3D represents LR-spline volumes.
This module depends on:
- GoTools Core library
- GoTools trivariate
- GoTools LRSplines2D
- SISL library

Example programs corresponding to this module are listed in 
\link examples_LRSplines3D examples_lrsplines3D \endlink

\section  LRvol LR spline volumes
In \link lrsplines2d \endlink it is described how the tensor-product 
structure that are the simplest way of generating a spline surface can lead to an
explosion in data size for surface with much local detail. This is even more so
for volumes.

LRSplineVolume inherits \link Go::ParamVolume ParamVolume \endlink and is a spline volume with the property of local refinement.
The lrsplines3D module in GoTools provides functionalities for working with LR (Locally Refined) spline volumes. These volumes are a powerful tool for representing complex 3D geometries with local refinement capabilities, offering flexibility in controlling the level of detail. 

An LR spline volume, \link Go::LRSplineVolume \endlink, is a piecewise polynomial or piecewise rational polynomial volume defined on an LR-mesh, \link Go::Mesh3D \endlink. An LR-mesh is a locally refined mesh made by applying a sequence of
refinements starting from a tensor-product mesh. LR spline volume are algorithmically defined throughout the refinement process of the mesh.
An LR spline volume is defined as

\f[ F(u,v,w) = \sum _{i=1}^L P_i s_i R_{i, p_1, p_2, p_3} (u,v,w) \f]

where P<SUB>i</SUB>, i=1, ... ,L are the surface coefficients, 
and s<SUB>i</SUB> are scaling factors introduced
to ensure partition of unity of the resulting collection of tensor product B-splines.
The tensor product B-splines R<SUB>i</SUB>(p<SUB>1</SUB>,p<SUB>2</SUB></SUB>,p<SUB>3</SUB></SUB>) 
are  of bi-degree (p<SUB>1</SUB>,p<SUB>2</SUB>,p<SUB>3</SUB>) 
defined on knot vectors of lengths p<SUB>1</SUB> + 2, p<SUB>2</SUB> + 2 and p<SUB>3</SUB> + 2 
on the parametric domain in the \[u\], \[v\] and \[w\] directions respectively. 

\image html vol_mesh.png "LR-mesh corresponding to volume. The mesh is shown before and after insertion of one mesh rectangle in the second parameter direction. The mesh is visualized by the mid point and corner curves of the elements."  width=600px 
A trivariate mesh corresponding to an LR volume of degree two is shown in the 
figure above. The first picture shows the mesh of a tensor-product volume. In the 
second picture the volume is refined by inserting one mesh rectangle. The example
program \link refine_lrvol refine_lrvol.C \endlink illustrates the refinement
procedure.

As for LR spline surfaces it is required that a new mesh rectangle splits the
support of at least one B-spline. The refinement procedure is similar to the
surface case.

\subsection Classes Classes involved in representing an LR spline volume

The core of lrsplines3D revolves around the LRSplineVolume class, which combines an underlying Mesh3D with a set of LRBSpline3D basis functions.

\link Go::LRVolApprox  LRVolApprox \endlink implements an adaptive and iterative algorithm for approximating a trivariate point 
cloud with format (x,y,z,f(x,y,z)) by local refinenement and approximation (\link Go::LRSpline3DMBA LRSpline3DMBA \endlink).

The key components and their relationships within the lrsplines3D module are as follows:



This class is the central component for defining, manipulating, and evaluating 3D LR spline volumes. An LRSplineVolume is essentially a sum of LRBSpline3D basis functions, each weighted by a control point, and defined over an adaptive Mesh3D.

\subsubsection Mesh3D
The LR-mesh is represented in \link Go::Mesh3D \endlink. It contains 
information about knot values in the three parameter directions and
mesh rectangles for each knot value. A mesh rectangle is described by the
indices of knots corresponding to the lower left and upper right corner of the
rectangle as well as multiplicity.

Mesh2D provides functionality to enquire properties of the mesh such as:
number of knots excluding multiplicity (numDistinctKnots), value of a given
knot (kval), iterators to knots (knotsBegin, knotsEnd), access to knots
(getKnots) and parameter domain (minParam, maxParam). See \link Go::Mesh3D \endlink
for the complete functionality.

\subsubsection LRBSpline3D
A \link Go::LRBSpline3D \endlink is constructed as a tensor product between
three univariate B-splines (\link Go::BSplineUniLR \endlink), but contains in 
addition the 
corresponding coefficient, the scaling factor and a possible rational weight.
The class contains information of the element in the support of the B-spline.

The class provides functionality to enquire the coefficient, scaling factor
and rational weight as well as geometry space dimension, associated knot 
vector and degree. The elements in the support are avaiable and the 
support limits can be requested. Functionality to evaluate position and derivatives
in a given parameter tripple is also available as well as the associated mesh.

\subsubsection Element3D
An element represents the domain of one polynomial patch in the LR spline
volume. It is limited by active mesh rectangles in the three parameter
direction. \link Go::Element3D \endlink contains information about the 
limits of this domain, the B-splines overlapping it and, in approximation
context, data points associated to this domain. Element3D provides access
to information about the domain properties, the associated B-eplines and 
neighbouring elements. The example program  
\link investigate_Element3D nvestigate_Element3D \endlink shows how to
obtain information related to the elements.

\subsubsection LRSplineVolume
\link Go::LRSplineVolume \endlink is a three-variate entity on which either a 
function (1D) or a 3D volume is represented.

The LRSplineVolume manages the geometric and topological information of the spline volume. It holds references to the underlying Mesh3D structure, which defines the knot intervals and element connectivity. It also maintains a collection of LRBSpline3D basis functions, which are the fundamental building blocks of the spline volume. The class provides methods for evaluating the volume at specific parameter points, evaluating on a grid, and performing refinement operations. The exmple program
\link evaluateLRVolume evaluateLRVolume.C \endlink demonstrates the various evaluation possibilities.

LRSplineVolume is the owner of all information required to represent the 
volume.
\par Key Data Members:

    Go::Mesh3D& mesh_: The underlying 3D mesh that defines the parametric domain and local refinement structure.
    BSplineMap bsplines_: An internal map of individual B-spline basis functions.
    ElementMap emap_: A map holding information about the individual mesh elements and their associated basis functions.
Univariate B-splines in the three parameter directions are contained in LRSplineVolume
to be referenced by LRBSpline3D.

Evaluation is performed through the functions point and elementGridEvaluate. The first
is implemented in several varieties: With and without computing
derivatives and with and without getting the relevant element as input. 
The latter evaluates a grid of points in a specified element. 

\subsection Classes2 Other classes, definitions and namespaces
- <b>LRVolApprox</b> Approximates a 4-dimensional point cloud by an LR spline function
(1D volume). The three first point coordinates are seen as the parameter tripple
corresponding to the points. The last coordinate is approximated.
 Refinement is performed according to the distance between the point
cloud and the volume and guided by a given tolerance. See \link Go::LRVolApprox \endlink

- <b>LRBSpline3DUtils</b> Provides utility functions for refinement related to LRBSpline3D objects.

- <b>Direction3D</b> Specifies the parameter direction of a volume, see
\link Go::Direction3D \endlink

- <b>LRSpline3DUtils</b>
Provides utility functions for spline volumes, mostly related to refinement,
but contains also functionality for evaluation of all B-splines in a specified
parameter tripple.

- <b>LRFeature3DUtils</b> Given a current LR spline volume with an
associated point cloud, compute feature output in a grid. Called from
LRVolApprox to visualize certain aspects of the approximation.

- <b>LRSpline3DMBA</b> Called from LRVolApprox. The name space provides 
functionality to update an LR spline volume using an adaptation to the  local
approximation method multi resolution B-spline approximation, see
\link Go::LRSpline3DMBA \endlink

- <b>LRSpline3DBezierCoefs</b> Bezier extraction. Coefficients of Bezier volumes 
are computed by interpolation of a set of sample points depending on the degrees of 
the LR B-spline volume. Only the quadratic and cubic cases are supported,
see also the example program \link Bezier_extraction Bezier_extraction.C \endlink

- <b>LRSpline3DEvalGrid</b> Grid evaluation of the elements of an LR spline volume.

- <b>LRVolStitch</b> Modifies a collection of trivariate LR spline functions 
organized in a regular pattern to obtain C<SUP>0</SUP> or C<SUP>1</SUP>
continuity between adjacent functions. The process involves an increase
in data size of the functions. See \link Go::LRVolStitch \endlink 
*/


#endif
