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

\link Go::Body \endlink is a boundary represented solid model. It is a virtual volume limited
by one or more shells, one outer and possibly a number of inner ones. The
shells are represented by the class SurfaceModel which is also the entity used
to represent a face set. In addition to being a topological entity,  \link Go::SurfaceModel \endlink
provides an interface which views the set of surfaces as one entity. Thus, an
operation like closest point computation can be performed without caring
about which surface is closest to the point.

The surface model consists of a set of faces represented by the class 
 \link Go::ftSurface \endlink. The face represents the abstract idea of a bounded surface, but it also
has a pointer to the geometric representation of this surface. The geometric
surface is a  \link Go::ParamSurface \endlink, it may be a NURBS surface, an elementary surface 
or a trimmed version thereof. The face is bounded by one or more loops,
one outer and possibly a number of inner ones. The inner loops represents
holes in the surface. Then the geometric surface will always be a  \link Go::BoundedSurface \endlink. All
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
this context, it inherits the abstract superclass \link Go::CompositeModel \endlink that defines
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
concrete tesselation results are  \link Go::LineStrip \endlink, 
 \link Go::RegularMesh \endlink and  \link Go::GenericTriMesh \endlink.

The control polygon of the composite model or a selected subset of the
model may be tesselated. In that case the output mesh is represented as a
 \link Go::LineCloud \endlink. LineCloud is an entity in gotools-core/geometry.

\subsection comp_sec2_2 CompositeCurve
A  \link Go::CompositeCurve composite curve \endlink 
is an ordered collection of curves. The class expects a set
of curves and will orient them order them in a sequence. The curve set should
preferably be continuous, but also discontinuous set may be represented.

A composite curve will contain a number of  \link Go::ParamCurves \endlink, ordering 
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

*/
#endif // _COMPOSITEMODEL_DOXYMAIN_H
