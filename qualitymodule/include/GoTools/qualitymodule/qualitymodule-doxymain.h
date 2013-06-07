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

#ifndef _QUALITYMODULE_DOXYMAIN_H
#define _QUALITYMODULE_DOXYMAIN_H

/// \page qualitymodule GoTools QualityModule
///
/// This GoTools module consists of a set of tools to check the
/// quality of CAD models.
///
/// \section tol Tolerances
/// The tests make use of the tolerances
/// associated with a surface model, i.e.
/// 
/// <dl>
/// <dt><b>gap</b></dt>
/// <dd>CAD models are seen as continuous if they are continuous
/// within this tolerance. If the distance between two points are less
/// than this tolerance, they are viewed as identical. A reasonable
/// tolerance is in the specter \f$[10^{-6}, 10^{-2}]\f$, but it
/// depends on the expected quality of a model.</dd>
/// 
/// <dt><b>neighbour</b></dt>
/// <dd>This tolerance is used in adjacency analysis and
/// it represents the maximum distance between neighbouring surfaces
/// or curves. If two curves or surfaces lie more distant than this
/// tolerance, the entities are found not to be adjacent. The
/// neighbour tolerance should be larger than the gap tolerance. If
/// nothing specific is known, a factor of 10 makes sense, but if the
/// gap tolerance is really small, a larger factor should be
/// used. Surfaces that lie closer to each other than the neighbouring
/// tolerance is found to be adjacent, but if the distance between the
/// surfaces somewhere is larger than the gap tolerance, the surface
/// set contains gaps. This is an error in the model.</dd>
/// 
/// <dt><b>bend</b></dt>
/// <dd>If two surfaces meet along a common boundary and
/// corresponding surface normals form an angle which is larger than
/// this tolerance, it is assumed that there is an intended sharp edge
/// between the surfaces.  Similarly, if two curves in a composite
/// curve meet with an angle larger than this tolerance, there is an
/// intentional corner.</dd>
/// 
/// <dt><b>kink</b></dt>
/// <dd>If two adjacent curves or surfaces meet with an angle
/// less than this tolerance, they are seen as \f$G^1\f$
/// continuous. It the angle is larger than this tolerance, but less
/// than the bend tolerance, the intended \f$G^1\f$ continuity is
/// broken and it is an error in the model. The tolerance depends on
/// the continuity requirements of the application. One suggestion is
/// \f$10^{-2}\f$. The bend tolerance must be larger than the kink
/// tolerance, for instance by a factor of 10. Both angular tolerances
/// are given in radians.</dd>
/// </dl>
/// 
/// 
/// \section tests Tests
/// The available quality tests can be classified as follows:
/// 
/// <dl>
/// <dt><b>Face continuity</b></dt>
/// <dd>Checks for continuity between faces in a surface
/// model. Positional and tangential discontinuities with respect to
/// the gap and the kink tolerance, respective, are returned.
/// Discontinuities larger than the neighbour tolerance or bend
/// tolerance are assumed to be intentional, and thus no error is
/// reported.</dd>
/// 
/// <dt><b>Edge continuity</b></dt>
/// <dd>Checks for positional and tangential discontinuity between
/// edges in the model.  The same tolerances are used as for face
/// continuity.</dd>
/// 
/// <dt><b>Accuracy of bounding entities</b></dt>
/// <dd>Whether or not the distance between an edge and its associated
/// face, a vertex and its associated edges and a vertex and its
/// associated faces is less than the gap tolerance.</dd>
/// 
/// <dt><b>Acute angles</b></dt>
/// <dd>Check for acute angles between adjacent faces or edges. The
/// kink tolerance is applied.</dd>
/// 
/// <dt><b>Degeneracies</b></dt>
/// <dd>Identify surfaces with boundary curves degenerating to a point
/// and surfaces with degenerate corners. Identify vanishing surface
/// normals and curve tangents, intersections between boundary loops
/// belonging to a trimmed surface and self intersections of boundary
/// loops. The gap tolerance is used in the intersection and self
/// intersection tests while the neighbour and kink tolerance is used
/// in degeneracy tests. The gap tolerance is used in tests for
/// vanishing tangents and normals.</dd>
/// 
/// <dt><b>Small entities</b></dt>
/// <dd>Identify small edges, small faces, sliver surfaces, and narrow
/// regions in faces. The narrow region test relates to the neighbour
/// tolerance while the other tests use specific tolerances.</dd>
/// 
/// <dt><b>Consistency of orientation</b></dt>
/// <dd>Check for consistency in loops, i.e. the curves defining the
/// loop are head to tail oriented, or in surface models. For closed
/// face sets, all surface normals should point out of the model or
/// into it, preferable out. Also open face sets should have
/// consistent surface normal directions at face boundaries. However,
/// the face orientation test is currently not up to date and the
/// result can unfortunately not be trusted.</dd>
/// 
/// <dt><b>Identical entities</b></dt>
/// <dd>Check for identity of vertices, and identical or embedded
/// faces or edges. The neighbourhood tolerance is used.</dd>
/// 
/// <dt><b>Spline entity testing</b></dt>
/// <dd>Check spline curves and surfaces for \f$G^1\f$ and \f$C^1\f$
/// discontinuities. The test is localized to knots of high
/// multiplicity. The gap and the kink tolerances are used for these
/// tests. Check also for spline entities with close, but not
/// identical knots. Such entities can create problems for certain
/// operations. A specific tolerance is used in this test.</dd>
/// 
/// <dt><b>Compute curvature information with respect to single
/// entities</b></dt>
/// <dd>The minimum curvature radius and curvature radii less than a
/// given threshold can be obtained. The test is applied to all curve
/// or surface entities in a model.</dd>
/// </dl>
/// 
/// The tests described above are localized in the class \link
/// Go::FaceSetQuality FaceSetQuality\endlink.  See the documentation
/// of this class for more information.


#endif // _QUALITYMODULE_DOXYMAIN_H
