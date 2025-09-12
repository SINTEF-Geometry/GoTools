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

#ifndef _LRSPLINES2D_DOXYMAIN_H
#define _LRSPLINES2D_DOXYMAIN_H
/**
 \page lrsplines2d GoTools LRsplines2D

 This module contains classes and tools for working with LR spline surfaces.

This module depends on the following GoTools modules:
- GoTools Core
and the following modules external to GoTools:
- SISL (SINTEF)

Example programs corresponding to this module are listed in 
\link examples_LRSplines2D examples_lrsplines2D \endlink

\subsection  TP The problem with tensor-product grids

The simplest way to generate parametric surface patches is by applying the tensor-product construction to univariate parametric space curves. Examples of these are B-spline surfaces, B-spline volumes, NURBS surfaces and NURBS volumes, which are implemented in GoTools.

The number of control points of these tensor-product patches grows exponentially with the dimension of the geometric object. For instance, for curves with 10 control points, corresponding tensor-product surfaces will have 100 control points and tensor-product volumes will have 1000 control points. These control points form a rectlinear gridded structure, making it impossible to refine the model locally without refining across the entire domain. Hence, in practice, it becomes unfeasible to add sufficient detail to the model where it is needed most.

\subsection LR LR B-splines and LR spline surfaces

Locally Refined B-splines (LR B-splines) aim to solve this problem by providing a mathematical framework, generalizing the tensor-product construction, for refining the model locally. This framework extends to any dimension.

An LR spline surface, \link Go::LRSplineSurface \endlink, is a piecewise polynomial or piecewise rational polynomial surface defined on an LR-mesh, \link Go::Mesh2D \endlink. An LR-mesh is a locally refined mesh made by applying a sequence of
refinements starting from a tensor-product mesh. LR spline surfaces are algorithmically defined throughout the refinement process of the mesh.
An LR spline surface is defined as

\f[ F(u,v) = \sum _{i=1}^L P_i s_i R_{i, p_1, p_2} (u,v) \f]

where P<SUB>i</SUB>, i=1, ... ,L are the surface coefficients, 
and s<SUB>i</SUB> are scaling factors introduced
to ensure partition of unity of the resulting collection of tensor product B-splines.
The tensor product B-splines R<SUB>i</SUB>,p<SUB>1</SUB>,p<SUB>2</SUB></SUB> 
are  of bi-degree (p<SUB>1</SUB>,p<SUB>2</SUB>) 
defined on knot vectors of lengths p<SUB>1</SUB> + 2 and p<SUB>2</SUB> + 2 
on the parametric domain in the \[u\] and \[v\] directions respectively. 

\image html mesh.PNG "LR-mesh. The support of one tensor product B-spline visualized as a red pattern.  Initial knotlines are shown as black lines, the inserted knotline segments are blue."  width=600px 
An LR-mesh corresponding to an LR spline surface of bidegree two is shown in
the figure above. The mesh lines of the initial tensor-product surface are drawn
in black. The corresponding knots are: [u<SUB>1</SUB>, u<SUB>1</SUB>, u<SUB>1</SUB>, u<SUB>2</SUB>, u<SUB>4</SUB>, u<SUB>6</SUB>, u<SUB>7</SUB>, u<SUB>7</SUB>,u<SUB>7</SUB>] in the first parameter direction and [v<SUB>1</SUB>, v<SUB>1</SUB>, v<SUB>1</SUB>, v<SUB>3</SUB>, v<SUB>5</SUB>, v<SUB>6</SUB>, 4<SUB>6</SUB>, v<SUB>6</SUB>] in the second direction. 
The LR-mesh is constructed by first inserting knots at v<SUB>2</SUB> and 
v<SUB>4</SUB> covering parts of the parameter domain in the second parameter
direction follwed by knots at u<SUB>3</SUB> and u<SUB>5</SUB>. Note that a new
knot line must cover the support of at least one B-spline in the current 
surface. The parameter patches limited by knot lines are denoted elements,
\link Go::Element2D \endlink, and serve a role in manouvring in the parameter
domain of the surface. In approximation context (link to page) the elements
contain data points. The initial B-spline surface has 30 basis 
functions (B-splines) while the
constructed LR spline surface has 42  B-splines, \link Go::LRSpline2D \endlink. 
Some lines of the LR-mesh intersecting the support of such B-spline do not correspond to knotlines of its knot mesh as they do not traverse the support completely. The refinement process is performed in the example \link refine_lrsurf 
refine_lrsurf.C \endlink .


The procedure applied is the following:
<ol>
<li> Add a new line segment that triggers the refinement of at least one existing tensor-product B-spline.  It can be an extension of an existing line segment, can be independent of existing line segments, or increase the multiplicity of an existing the line segment. Thus, a line segment going from 
(u<SUB>5</SUB>,v<SUB>4</SUB>) to (u<SUB>5</SUB>,v<SUB>5</SUB>) is a legal choice.
<li> Subdivide all tensor product B-splines with support that is completely traversed by the new line/extended line.
<li> The tensor product B-splines are required to have minimal support. This means that all line segments traversing the support of a tensor product B-spline are required to be knotline segments of the B-spline taking the  multiplicity of line segments into account. After
inserting a new line segment and performing the subdivision in Step 2., there might still be tensor product B-splines that do not have minimal support with respect to the LR-mesh. 
Consequently all such B-splines  must be refined.
This process is continued until all tensor product B-splines have minimal support.
</ol>
If more than one new knotline segment is defined simultaneously, the refinement process is applied one segment at the time. Some details on how to choose
new knot line segments can be found in \link lrsplines2d_refine \endlink . 

LR spline surface possess most of the nice properties of spline surfaces, such as:
- Non-negative basis function
- Partition of unity
- The surface lies in the convex hull of its coefficients

However, depending on the refinement strategy, linear dependence situations may
occur. A case where such situations are identified and resolves is explained
in the example \link identify_and_resolve_linear_dependence
identify_and_resolve_linear_dependence.C \endlink.

\subsection Classes Classes involved in representing an LR spline surface

\subsubsection Mesh2D
The LR-mesh is represented in \link Go::Mesh2D \endlink. It contains 
information about knot values in the two parameter directions and
segments for each knot value. The latter has a compact encoding where a 
knot line segment is represented by the index of the start knot and the
knot multiplicity. Active segments have multiplicity one or higher while
segments with multiplicity 0 indicates that there is no knot in this part 
of the parameter domain.

Mesh2D provides functionality to enquire properties of the mesh such as:
number of knots excluding multiplicity (numDistinctKnots), value of a given
knot (kval), iterators to knots (knotsBegin, knotsEnd), access to knots
(getKnots), parameter domain (minParam, maxParam). Also mesh rectangle 
information such as knot multiplicity and active and not active knot line 
segments, is available.

\subsubsection BSplineUniLR
\link Go::BSplineUniLR \endlink represents a univariate B-spline by 
storing indices to the corresponding knot vector in Mesh2D to the active knots
of the current B-splines. Properties like degree, knot values and
parameter interval are deduced. Information about a possible overlap between
two univariate B-splines and the Greville parameter of this B-splines is
available.

\subsubsection LRBSpline2D
A \link Go::LRBSpline2D \endlink is constructed as a tensor product between
two univariate B-splines (BSplineUniLR), but contains in addition the 
corresponding coefficient, the scaling factor and a possible rational weight.
The class contains information of the element in the support of the B-spline.

The class provides functionality to enquire the coefficient, scaling factor
and rational weight as well as geometry space dimension, associated knot 
vector and degree. The elements in the support are avaiable and the 
support limits can be requested. The associated mesh is also available.

\subsubsection Element2D
An element represents the domain of one polynomial patch in the LR spline
surface. It is limited by active knot line segments in both parameter
direction. \link Go::Element2D \endlink contains information about the 
limits of this domain, the B-splines overlapping it and, in approximation
context, data points associated to this domain. Element2D provides access
to information about the domain properties, the associated B-eplines and 
neighbouring elements. The example program  
\link investigate_Element2D nvestigate_Element2D \endlink shows how to
obtain information related to the elements.

\subsubsection LRSplineSurface
\link Go::LRSplineSurface \endlink inherites \link Go::ParamSurface \endlink
and thus the functionality defined there. Some functionalities are not
 implemented, e.g. area. 

LRSplineSurface is the owner of all information required to represent the 
surface. Bivariate B-splines and elements are
contained in maps to facilitate rapid recognition of the element and B-splines 
having a given parameter pair in its support. A parameter remembering the last
element requested serves as a tool to speed up this search. The class also 
contains all individual univariate B-splines required to define the bivariate
B-splines.

Evaluation is performed through the functions point and evalGrid. The first
is implemented in several varieties: With and without computing
derivatives and with and without getting the relevant element as input. 
The latter
evaluates a regular grid of points, one by one, after connecting the elements
with a tensor-product grid. Connecting a parameter pair to an element is the
most time consuming part of the evaluation, and the grid ensures a rapid 
recognition. Evaluation is demonstrated in the example \link evaluateLRSurface.C
evaluateLRSurface \endlink . See also \link investigate_LRSplineSurface.C
investigate_LRSplineSurface \endlink for more functionality.

The LRSplineSurface is constructed through refinement given a tensor-product 
grid as input. The constructor can
receive a spline surface, \link Go::SplineSurface \endlink, or the information 
defining the grid as input. Refinement is triggered by LRSplineSurface, but
performed in \link Go::LRSplineUtils \endlink. Given an element or a B-spline
selected for refinement, a knot line segment that complies to the
requirements, see \link LR \endlink, must be specified. Several options
exists and are explained in some detail in \link lrsplines2d_refine \endlink .

\subsection Classes2 Other classes, definitions and namespaces
- <b>BSplineUniUtils</b> Utility functionality for keeping track of 
univariate B-splines. The univariate B-splines are stored as vector in
LRSplineSurface and these vectors need to be kept updated during refine 
operations.
- <b>Direction2D</b> Specifies the parameter direction of a surface, see
\link Go::Direction2D \endlink
- <b>LinDepUtils</b> The peeling algorithm, one step in identifiying linear
dependencies, see \link Go::LinDepUtils \endlink
- <b>LogLikelyhood</b> A statistical critierion for goodness of fit. 
Related to scattered data approximation.
- <b>LRApproxApp</b> Functionality related to the approximation of a point 
cloud by an LR spline surface: specific interfaces to the approximation and
computation of accuracy, see \link Go::LRApproxApp \endlink and 
example \link comparePointsLRSurf3D.C comparePointsLRSurf3D \endlink
- <b>LRBSpline2Dutils</b> LRBSpline2D related functionality used in refinement.
- <b>LRFeatureUtils</b> Given a current LR spline surface with an
associated point cloud, compute feature output in a grid. Called from
LRSurfApprox to visualize certain aspects of the approximation, see the
help documentation in \link PointCloud2LR.C \endlink
- <b>LRMinMax</b> Computes extrema of LR spline function. The functionality
requires associated contour curves as input, see \link Go::LRMinMax \endlink
- <b>LRSplineEvalGrid</b> Grid evaluation of the elements of an 
LR spline surface, \link Go::LRSplineEvalGrid \endlink
- <b>LRSplineMBA</b> Called from LRSurfApprox. The name space provides 
functionality to update an LR spline surface using an adaptation to the  local
approximation method multi resolution B-spline approximation, see
\link Go::LRSplineMBA \endlink
- <b>LRSplineUtils</b> Utilities, mostly related to refinement of an 
LR spline surface, but the namespace contains also some more functionality,
see \link Go::LRSplineUtils \endlink
- <b>LRSurfApprox</b> Approximate a scattered data point cloud by a 1D or 3D
LR spline surface. In the latter case, the point cloud must be parameterized. 
A large collection of parameters can be used to guide the approximation.
A short explanation can be found in the help text to the application
\link PointCloud2LR.C \endlink , which provides an interface to the
functionality. See also \link Go::LRSurfApprox \endlink and the example
programs \link approximateWithLRFunc approximateWithLRFunc.C \endlink and
\link approximateParPointsWithLRSurf approximateParPointsWithLRSurf.C \endlink .
A simplified interface can be found in \link Go::LRApproxApp \endlink.
- <b>LRSurfSmoothLS</b> Least squares approximation with a smoothing term.
Called from LRSurfApprox. This approximation approach is combined
with multi resolution B-spline approximation (MBA) to compute the
approximating surface. Least squares approximation is typically used in the
start of the process, then the process continues with MBA, 
\link Go::LRSurfSmoothLS \endlink
- <b>LRSurfStitch</b> Modifies a collection of bivariate LR spline functions 
organized in a regular pattern to obtain C<SUP>0</SUP> or C<SUP>1</SUP>
continuity between adjacent functions. The process involves an increase
in data size of the functions. See \link Go::LRSurfStitch \endlink 
- <b>LRTraceIsocontours</b> Compute the level-set curves of an LR spline 
function. The function is used to compute contour curves and the use is
illustrated in the app \link isoContours.C \endlink and the
example \link isoContoursLRFunc isoContoursLRFunc.C \endlink, see also
\link Go::LRTraceIsocontours \endlink
- <b>SSurfTraceIsocontours</b> Used by LRTraceIsocontours. 
- <b>TraceContoursTypedefs</b> Typedefs and tempate functionality used by 
LRTraceIsocontours. 
- <b>Mesh2DUtils</b> Utility functionality for manouvring in an LR mesh.
- <b>MeshLR</b> Base class for Mesh2D and \link Go::Mesh3D \endlink in 
lrsplines3D.
- <b>TrimSurface</b> Create bounded surface from an LR B-spline function 
and associated point cloud to make the function domain correspond to the
domain of the point cloud represented by its x- and y-values, see
\link Go::TrimSurface \endlink
- <b>LRSplinePlotUtils and LRSplinePlotUtils2 </b> Debug functionality

\subsection Sources

    T. Dokken, T. Lyche, K.F. Pettersen. Polynomial splines over locally refined box-partitions, Computer Aided Geometric Design, Volume 30, Issue 3, March 2013, Pages 331--356
**/
#endif // _LRSPLINES2D_DOXYMAIN_H
