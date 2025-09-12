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

#ifndef _COMPOSITEMODEL_DOXYMAIN_H
#define _COMPOSITEMODEL_DOXYMAIN_H

//===========================================================================
//                        DOCUMENTATION ON MAIN PAGE
//===========================================================================

/**

\page compositemodel GoTools CompositeModel

Contains tools for representing composite models.

This module depends on the following GoTools modules:
- GoTools Core
- GoTools Implicitization
- GoTools Intersection
- GoTools Parametrization
- GoTools Igeslib
- GoTools Topology

and the following modules external to GoTools:
- SISL (SINTEF)
- TTL (SINTEF)
- The newmat matrix library, which is freely available from <a
href ="http://www.robertnz.net/">www.robertnz.net</a>.

SISL, TTL and newmat is included for convenience.

\section comp_sec1 Brep topology
A CAD model can normally not be represented by one geometry entity like
the ones in gotools-core. The need for sets of entities, adjacency analysis and
representation of adjacency arises. There is a need for topological structures.

\link Go::Body Body \endlink is a boundary represented solid model. It is a virtual volume limited
by one or more shells, one outer and possibly a number of inner ones. The
shells are represented by the class SurfaceModel which is also the entity used
to represent a face set. In addition to being a topological entity,  \link Go::SurfaceModel SurfaceModel \endlink
provides an interface which views the set of surfaces as one entity. Thus, an
operation like closest point computation can be performed without caring
about which surface is closest to the point.

The surface model consists of a set of faces represented by the class 
 \link Go::ftSurface ftSurface \endlink. The face represents the abstract idea of a bounded surface, but it also
has a pointer to the geometric representation of this surface. The geometric
surface is a  \link Go::ParamSurface ParamSurface \endlink, it may be a NURBS surface, an elementary surface 
or a trimmed version thereof. The face is bounded by one or more loops,
one outer and possibly a number of inner ones. The inner loops represents
holes in the surface. Then the geometric surface will always be a  \link Go::BoundedSurface BoundedSurface \endlink. All
faces have an outer loop which either represents an outer trimming curve in
the case of a BoundedSurface or the surface boundary.

A loop consists of a number of edges represented by ftEdge. The edges
have a geometric representation by a ParamCurve. This is the level where
adjacency information related to faces or surfaces is stored. An ftEdge is a
half edge. It knows the face to which it belongs and the geometry of the curve
it corresponds to. This curve is normally a CurveOnSurface curve. Then,
both parameter and geometry information are available. The ftEdge has a
pointer to a twin edge. This edge belongs to the adjacent ftSurface meeting
the current one along the boundary represented by the current ftEdge. This
construction is a very flexible one, but it does not automatically ensure C 0
continuity as the curves pointed at by the two twin edges, are normally not
the same.
An edge is limited by two vertices which has a geometric representation
as a point. A vertex has knowledge about all edges meeting in this point.

\section comp_sec2 CompositeModel
SurfaceModel represents the shell in a boundary represented geometry model,
but it is also the entity that is used to view a surface set as one unity. In
this context, it inherits the abstract superclass \link Go::CompositeModel CompositeModel \endlink that defines
an interface for a unity of geometry entities of the same type. All composite
models can
Report how many entities it consist of
\li Evaluate a given entity in a given parameter or parameter pair
\li Compute the closest point from a given pont to the entity set
\li Compute the extremal point of the entity set in a given direction
\li Intersect itself with a plane
\li Intersect itself with a line
\li Compute the bounding box surrounding the entity set
\li Check whether one specified entity is degenerate
\li Turn parameter directions corresponding to one or all entities
\li Tesselate itself
\li Clone itself 

Tesselation is performed in the gotools-core submodule tesselator. The methods 
here expect information about the tesselation resulution. For each parameter 
direction in the object, the number of tringles or line strips must
be specified. The resolution is defined from the composite model. This can
either be done setting a resolution for all entities or by a density parameter
that governs the resolution. The composite model tesselation routines return
a vector of shared pointers to a GeneralMesh. The type of meshes used for the
concrete tesselation results are  \link Go::LineStrip LineStrip \endlink, 
 \link Go::RegularMesh RegularMesh \endlink and  \link Go::GenericTriMeshGenericTriMesh  GenericTriMeshGenericTriMesh  \endlink.

The control polygon of the composite model or a selected subset of the
model may be tesselated. In that case the output mesh is represented as a
 \link Go::LineCloud LineCloud \endlink. LineCloud is an entity in gotools-core/geometry.

\subsection comp_sec2_2 CompositeCurve
A  \link Go::CompositeCurve composite curve \endlink 
is an ordered collection of curves. The class expects a set
of curves and will orient them order them in a sequence. The curve set should
preferably be continuous, but also discontinuous set may be represented.

A composite curve will contain a number of  \link Go::ParamCurve ParamCurves \endlink, ordering 
information and continuity information. The curve sequence will be parameterized as a 
unity. Curve indices follow the indices for the curves given as
input to the object. The class has the following type specific functionality:
\li Fetch an individual curve
\li Fetch the parameter interval of the composite curve
\li Given a parameter value of one individual curve, get the corresponding parameter 
value in the composite curve, and vice versa. Given a
parameter value of the composite curve, get the index and parameter
value of the individual curve.
\li Evaluate the composite curve as one unity
Hit test. Compute the closest point in the composite curve to a given
point along a given line
\li Append a new curve to the composite curve

All the curves in the composite curve can be tesselated with respect to a
default value, or a value given by the application, or a given subset of the
composite curve can be tesselated with respect to the defined resolution. It
is also possible to set the number of line segments for each curve relative
to a given density parameter. Then the length of each line segment should
roughly approximate this density. The density should be set according to
the total size of the model. An upper bound of the tesselation size exists,
but a large model combined with a small density will still lead to a heavy
visualization model.

\subsection comp_sec2_3 CurveModel
A  \link Go::CurveModel curve model \endlink is an incomplete class 
inheriting a composite model. Not all
the functionality defined in CompositeModel is implemented in this class. Its
main purpose is to take a set of curves, for instance read from an IGES file,
compute the topology of the curve set, use the topology information to order
the curves into sequences, and return these sequences as composite curves.

Nevertheless, a curve model may be evaluated, tesselated and we can
fetch a particular curve in the set by index or find the associated index given
a curve. The index corresponds to the order of the input curves to this class.
It is also possible to compute the bounding box surronding the curve set and
the curve model may copy itself.

\subsection comp_sec2_4 SurfaceModel
A  \link Go::SurfaceModel surface model \endlink represents a surface set or a shell. 
It contains a number
of parametric surfaces and topology information representing adjacency 
between surfaces in the set. The input surfaces to a surface model may be
given as parametric surfaces or faces, i.e.  \link Go::ftSurface \endlink. 
In the latter case, the
adjacency information may be given or not. A surface model may consist of
several disjoint surface sets.
SurfaceModel has some specific functionality:
\li Reset topology tolerances and recompute adjacency information
Return a specified surface or face
\li Perform intersection with a line or a plane and represent the output in
a surface specific manner
\li Intersect the surface set with a plane, and trim the result with respect
to the orientation of the plane
\li Hit test. Compute the closest point in the surface set to a given point
along a given line
\li Append one or more new faces to the surface set
\li Append another surface set to the current one
\li Return information about the number of boundary loops in the surface
set. One boundary loop will correspond to either a continuous subset
of surface or a hole in the surface set.
\li Return informations of gaps between surfaces in the set
\li Return informations of kinks between surfaces in the set
\li Return informations of sharp edges between surfaces in the set
\li Triangulate the surface set with respect to a density parameter, and join
the individual triangulations for each surface into one triangulation.
\li Fetch information about neighbouring surfaces with inconsistent direc-
tion of the surface normals
\li Fetch all vertices in the surface set
\li Fetch vertices lying at the boundaries of the surface set

Some functionality is special for isogeometric models where each surface
is expected to be a non-trimmed spline surface and the surfaces should lie in
a corner-to-corner configuration

\li Check if all surfaces are splines
\li Check if the surface configuration is corner-to-corner
\li Ensure a corner-to-corner configuration if this is not the case. The
function applies only to spline surfaces
\li Ensure that the spline surfaces in the model shares the spline space at
common boundaries

All the surfaces in the surface model can be tesselated with respect to a
default value or a pair of values given by the application, or a given subset
of the surface model can be tesselated with respect to this resolution. The
application can also specify the approximate number of quads in the tessela-
tion. Each quad are divided into two triangles. Then a long and thin surface
will get more triangles in the long direction. It is also possible to set the
number of triangles in each parameter direction for each surface relative to
a given density parameter. Then the length of each triangle should roughly
approximate this density. The density should be set according to the total
size of the model. An upper bound of the tesselation size exists, but a large
model combined with a small density will still lead to a heavy visualization
model.

\subsection comp_sec2_5 Tolerances
The constructor of subclasses of \link Go::CompositeModel compositemodel \endlink requires a set of tolerances to be defined. With the exception of approxtol, these tolerances are used to compute adjacency between entities and to check the continuity between adjacent entities.

\subsection comp_sec2_6 Intersection results

Intersections performed in the composite model class either return the intersection results in a format particular to the concrete composite model subclass, or the intersection results are stored in a subclass of \link Go::IntResultsModel IntResultsModel \endlink. In the latter case, intersections can be performed without caring about the composite model type.

IntResultsModel stores the intersection results and the entities involved in the intersection and has the following functionality:

\li Report upon the existence and number of intersection points and intersection curve segments
\li Tesselate itself

The actual geometry of the intersection results must be fetched from the subclasses
\link Go::IntResultsSfModel IntResultsSfModel\endlink and \link Go::IntResultsCompCv IntResultsCompCv \endlink. The intersection results classes are very lean, but have the potential to get a rich set of functionality operating on intersection results.

SurfaceModel returns intersection results either as IntResultsSfModel or as
\link Go::ftCurve ftCurve \endlink and a vector of \link Go::ftPoints ftPoints \endlink. It is also possible to fetch intersection results from IntResultSfModel as ftCurve and ftPoint.

An ftCurve is composed of a number of ftCurveSegments. An ftCurveSegment represents a piece of an intersection curve. The segment is represented by a spatial curve and a curve in the parameter domain of the parametric entities involved in the intersection, which will be one or two. Thus, the curve segments distinguish between the individual surfaces in a surface set while ftCurve itself does not. ftCurve keeps track of the continuity between the segments it is composed of. An ftCurve is able to tesselate itself.

An ftPoint relates to one face, i.e., an \link Go::ftSurface ftSurface\endlink,
and has a geometric representation and a parameter domain representation. Each intersection point, for instance in an intersection between a surface model and a line, is represented as an ftPoint. The application has access to the position of this point, to the face with which it is associated and the parameter values in this face.

A composite curve represents its intersection results as an instance of
\link Go::IntResultsCompCv IntResultsCompCv\endlink. Intersection points, which are the expected result in this case, are represented as \link Go::PointOnCurve PointOnCurve\endlink. Intersection curve segments are represented as pairs of PointOnCurves. The two points limit the extension of the intersection curve segment. The class
PointOnCurve lies in the gotools-core/geometry submodule and holds a parametric curve, a parameter value and the geometric representation of the point. This content is accessible from the public user interface.

\subsection comp_sec2_7 The topological face
The topological face (\link Go::ftSurface ftSurface\endlink)
has access to the geometric representation of this face and information about neighbouring faces. In addition, it has knowledge about the solid (\link Go::Body Body\endlink) it belongs to, if any.

The geometric surface is of type ParamSurface and described in the documentation of the gotools-core/geometry submodule.

The face is limited by a number of \link Go::Loop Loops\endlink. The geometric surface corresponding to a face corresponds to the trimmed face, so the information present in the boundary loops of the face does also exist in the surface description.

ftSurface has a rich set of functionality of different types:

 - Access functionality
       - Fetch the corresponding surface
       - Fetch boundary loops
       - Fetch edges belonging to the boundary loops
       - Fetch neighbouring faces
       - Fetch all vertices or subsets of vertices
       - Fetch the body owning this face
       - Fetch all bodies meeting in this face, maximum 2
       - Fetch the coincident face in a volume model
       .
 - Quality checking
       - Check for discontinuities in the surface
       - Check the distance between a face and corresponding edges and vertices
       - Check the orientation of loops belonging to the face
       - Check if the face has a narrow region
       - Check for acute angles in the boundary loops
       .
 - Functionality related to an isogeometric model
       - Check if the corresponding surface is a spline
       - Check if this face and the neighbouring face have corresponding spline spaces
       - Ensure that this face and the neighbouring face have corresponding spline spaces
       - Check if this face and the neighbouring face satisfy a corner-to-corner condition
       - Ensure that this face and the neighbouring face satisfy a corner-to-corner condition
       - Fetch adjacency information between this face and a neighbouring face
       - Fetch information about boundaries where the face has no neighbouring face
       - Fetch information about coefficient enumerations for spline surface belonging to adjacent faces and free boundaries
       .
 - Other queries
       - Evaluations
       - Closest point computations
       - Perform smoothing on the current face
       - Compute bounding box

\subsection comp_sec2_8 The topological edge
The topological edge is implemented in the class \link Go::ftEdge ftEdge\endlink
and it is a half-edge. Thus, it contains information related to one face only. However, it has access to the corresponding half-edge when the edge represents a boundary between two adjacent faces. The topological edge has a geometric representation as a ParamCurve. The edge may represent a restriction of the curve. This restriction is implemented using limiting parameters in the parameter interval of the curve.

The main functionality is:

\li Fetch the curve representing the geometry description of the edge
\li Fetch the parameters limiting the edge compared to the corresponding curve
\li Fetch the face to which the edge belongs
\li Given an edge parameter, compute the corresponding face parameter
\li Access to geometry information like bounding box, results of evaluation, closest point computations and an estimate of the curve length
\li Information about the vertices limiting the edge, i.e., access to vertices and the edge parameter corresponding to a vertex
\li Access to the radial edge if this entity exist. A radial edge is defined in volume model and the entity is described in the [[documentation of the trivariatemodel module|Module TrivariateModel]].
\li Get all surfaces meeting in this edge, with a maximum of two


\section comp_sec3 Reverse engineering
The aim is to go from a triangulated point cloud to a boundary represented
CAD model. See \link reverseengineering_doc the reverse engineering page\endlink 
for more information.
*/
#endif // _COMPOSITEMODEL_DOXYMAIN_H
