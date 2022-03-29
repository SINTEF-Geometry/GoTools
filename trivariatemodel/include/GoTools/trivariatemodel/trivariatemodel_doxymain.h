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

#ifndef _TRIVARIATEMODEL_DOXYMAIN_H
#define _TRIVARIATEMODEL_DOXYMAIN_H

/**
\page trivariatemodel GoTools Trivariatemodel 

The module trivaratemodel contains the representation format for a 
set of volumes and
tools to act upon the volumemodel.

The module depends on:
- The GoTools Core library
- GoTools Composite Model
- GoTools Trivariate
- Modules that Composite Model depend upon 

\section volmod_sec1 Volume Topology
\image html vol_topology.gif "Topology structures for a volume model. The GoTools class names are given in brackets"

The volume topology in GoTools is an extension to the boundary represented
topology implemented in the module compositemodel. However, in contrast to
the boundary represented topology, the volume topology is not manifold. Thus,
we must expect more than two faces to meet in an edge. 

The top entity in the topology structure is the volume model 
( \link Go::VolumeModel VolumeModel \endlink) 
which consists
of a set of volumes. Each volume has a topological entity implemented in
 \link Go::ftVolume ftVolume \endlink 
and a geometrical representation implemented as 
 \link Go::ParamVolume ParamVolume \endlink. 
Information about the geometric representation of a volume can be found
in the module trivariate.

ftVolume inherits  \link Go::Body Body \endlink in compositemodel 
and is, as Body, surronded by one or more
shells represented as  \link Go::SurfaceModel SurfaceModels \endlink. 
A shell is a collection of faces
( \link Go::ftSurface ftSurface \endlink) which have a 
geometric representation as a 
 \link Go::ParamSurface ParamSurface \endlink. 
ParamSurface is described in the module gotools-core. 

A face is limited by a number of  \link Go::Loop loops \endlink 
that are sequences of edges 
( \link Go::ftEdge ftEdge \endlink).
In this context, the face is used to represent adjacency between volumes. 
An edge is limited by two  \link Go::Vertex vertices \endlink.

\section volmod_sec2 Volume Model
\image html composite_data_struct.gif "The inheritance tree for a volume model"
 \link Go::VolumeModel VolumeModel \endlink is a 
 \link Go::CompositeModel composite model \endlink as illustrated in 
the figure above. See also the description in the documentation of
composite model. As such it inherites a function interface, but some
methods are not implemented. Thus, the class is incomplete, but not in the
sense of being a topological entity.

VolumeModel has the following functionality:
- Fetch one entity in the set, either as a topological or geometric volume.
- Evaluate the volume model given information about which entity to
evaluate
- Compute  \link Go::BoundingBox bounding box \endlink
- Add a new volume to the volume model
- Remove one volume from the volume model
- Check if all geometrical volumes are NURBS
- Fetch all  \link Go::Vertex vertices \endlink in the model
- Fetch all  \link Go::EdgeVertex radial edges \endlink in the model
- Fetch non-radial edges, i.e. edges that does not belong to at least
two volumes.
- Fetch all inner faces in the model
- Fetch all boundary faces in the model
- Check if the model has a corner-to-corner configuration
- Ensure that the model has a corner-to-corner configuration
- Ensure that adjacent spline volumes share common spline spaces
at common boundaries
- Fetch the boundary of the volume model described as a
number of SurfaceModels
- Divide the volume model into connected volume models

The constructor of the volume model requires a number of \ref comp_sec2,
namely
gap, neighbour, kink and bend. These tolerances are the same as the ones
required for a composite model, and the tolerances are explained in the
documentation of the compositemodel module.

\section volmod_sec3 Topological Volume Entity
The topological volume entity is implemented in the class 
 \link Go::ftVolume ftVolume \endlink and it
inherits  \link Go::Body Body \endlink 
from the compositemodel module. The name ftVolume is choosen
to be in line with ftSurface and ftEdge. Those entities have got their names
for historical reasons. The prefix \em ft has no deeper meaning.

An ftVolume is a  \link Go::Body Body \endlink 
and has thus access to its shells, i.e. boundaries, and
can check whether two ftVolumes are adjacent. An ftVolume can also
- Fetch the corresponding geometry volume
- Fetch the bounding box of this volume
- Fetch all adjacent ftVolumes
- Fetch all radial edges belonging to this volume or being common
between this volume an another volume
- Fetch edges belonging only to the current volume
- Compute adjacency information between the current volume and another
volume
- Return information about outer boundaries
- Get local information about correspondence of coefficients of two
adjacent spline volumes
- Return the relation between the current volume and a given vertex

A ftVolume may be trimmed. That is, the boundary faces of the volume may limit the
extent of the underlying parametric volume. In this context the following functionality
apply:
- Check if the volume is hexahedral
- Approximate a hexahedral trimmed volume by one spline volume
- Split a trimmed volumes into a set of hexahedral trimmed volumes

<em> Note that the functionality related to trimmed volumes are under 
development and currently very unstable and limited. </em>

\section volmod_sec4 Extensions to Face
The initial implementation of a face as one entity in a boundary representation
solid was no longer
sufficient when the entity should serve as part of the boundary of a volume
belonging to a volume model. Some extensions turned out the be required.
The face  \link Go::ftSurface ftSurface \endlink 
has knowledge about the Body it belongs to, if any.
It has also a pointer to a boundary face belonging to an adjacent body. 
This pointer is used in the context 
of a volume model. During the assembly of a volume model, information about
coincident boundary faces of the ftVolumes are computed, and the pointer
representing adjacency is set accordingly.

An ftSurface provides access to its twin and to the, at most, two bodies sharing
the common face represented by this ftSurface.

When the ftSurface represents the boundary surface of an ftVolume, the
corresponding parametric surface will be of type 
 \link Go::SurfaceOnVolume SurfaceOnVolume \endlink. This class
is implemented in the module trivariate. A SurfaceOnVolume has knowledge
about the spatial representation of the surface, the ParamVolume on which it
belongs and the position of the surface in the parameter domain of the
volume. Currently, a SurfaceOnVolume is constructed only as a boundary surface
of a SplineVolume. It contains enough information to identify which 
volume boundary it represents and the orientation of the surface compared 
to this boundary.

\section volmod_sec5 Radial Edge
In a non-manifold model, more than two faces can meet in an edge, and
thus the half edge representation implemented in ftEdge is not sufficient
to hold the model. The radial edge, 
 \link Go::EdgeVertex EdgeVertex \endlink,
is an extension to the
topology structures of compositemodel and the class itself is placed in
compositemodel. An EdgeVertex contains information of
all pairs of half edges,  \link Go::ftEdge ftEdge \endlink, 
meeting in an edge, and each ftEdge
belonging to an EdgeVertex has access to this EdgeVertex.

An EdgeVertex instance is defined when adjacency between two ftSurfaces is
found. Thus, at edges in an ftVolume where no adjacent ftVolume exists, no
EdgeVertex instance will be defined. Then the edge is represented with a
pair of ftEdge instances.

The EdgeVertex class is placed in the module compositemodel, but is used
in connection with volume models. In addition to functionality related to
topology build, the class has some access functionality:
- Fetch edges meeting in the radial edge, either uniquely defined, i.e.
one instance for a pair or twins, or absolutely all edges.
- Fetch a specified edge
- Check if a ftEdge belongs to a radial edge
- Check if a ftEdge belongs to a radial edge and has no information about
a corresponding ftEdge.
In that case the corresponding surface does not belong to a volume.
- Fetch adjacent faces
- Fetch adjacent bodies

*/

#endif // _TRIVARIATEMODEL_DOXYMAIN_H
