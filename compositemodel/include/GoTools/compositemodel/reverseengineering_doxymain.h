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

#ifndef _REVERSEENGINEERING_DOXYMAIN_H
#define _REVERSEENGINEERING_DOXYMAIN_H

/**
\page reverseengineering_doc Reverse engineering, from triangulated point cloud to CAD model.

The defined workflow constitues an early prototype for reverse engineering. The 
functionality is expected to be applied to mechanical objects. The mathematical 
description of these types of objects is characterized by predominantly use of
primary surfaces, e.g. planes, cylinders, cones, spheres and tori. Free form
surfaces are used mostly for blends and small details. The extent of sharp edges in
the  objects is small. Edges are in general blended. Relevant objects are often
made of casted iron or created by adaptive manufacturing giving rough surfaces. 
Only the parts of the surface that are critical for assembly is plastered, leaving
most of the part with small irregularities. Simple models and the main surfaces 
of more complex models is expected to be reconstructed.

\section re_sec1 Interface
The reverse engineering engine is the class
\link Go::RevEng RevEng\endlink. The class is instanciated with a triangulated surface 
represented as an \link Go::ftPointSet ftPointSet\endlink, which represents a triangulation
in GoTools. Assume that the triangulation
is represented as a collection of vertices and triangles in the format
\verbatim
  vector<double> vertices;
  vector<int> triangles;
\endverbatim
The points are stored consequetive in vertices as (x,y,z) while three indices in triangles
represents one triangle. The triangles are expected to have a consistent orientation. The
vertices and triangles are transferred to ftPointSet with the following code snippet:
\verbatim
  size_t nvertices = vertices.size()/3;
  size_t ntriangles = triangles.size()/3;
  shared_ptr<ftPointSet> tri_sf = shared_ptr<ftPointSet>(new ftPointSet());
  for (size_t ki=0; ki<nvertices; ++ki)
    {
      Vector3D xyz(vertices[3*ki], vertices[3*ki+1], vertices[3*ki+2]);
      shared_ptr<RevEngPoint> currpt(new RevEngPoint(xyz, -1));
      currpt->setIndex(ki);
      tri_sf->addEntry(currpt);
    }

  for (size_t ki=0; ki<ntriangles; ++ki)
    {
      ftSamplePoint* pt1 = (*tri_sf)[triangles[3*ki]];
      ftSamplePoint* pt2 = (*tri_sf)[triangles[3*ki+1]];
      ftSamplePoint* pt3 = (*tri_sf)[triangles[3*ki+2]];
      pt1->addTriangle(pt2, pt3);
      pt2->addTriangle(pt3, pt1);
      pt3->addTriangle(pt1, pt2);
    }
\endverbatim
\link Go::RevEngPoint RevEngPoint\endlink represents a triangulation vertex and is inherited
from \link Go::ftSamplePoint ftSamplePoint\endlink. RevEngPoint is enhanced with information
such as estimated surface normal and curvature as well as associated functionality.

\section re_sec2 Overview
The reverse engineering process is organized as a sequence of operations that together 
consitute a work flow. The process is as follows:

 * <ol>
 * <li> Enhance points
 * <li> Classify points according to Gauss and mean curvature
 * <li> Segment point cloud into regions
 * <li> Surface creation
 * <li> Compute global properties such as main axes and update surfaces accordingly
 * <li> Edge creation
 * <li> Define blend surfaces
 * <li> Trim surfaces with respect to identified edges, blend surfaces and adjacent regions
 * <li> Create CAD model 

Point 4 to 6 are repeated three times, each time differently. The process
is automated, but organized as a sequence of commands to RevEng. This allows for storing
the state at a number of locations to resume the computation at a convenient time. Note that
storing and reading the state can be time consuming. In the following, we will describe
the process in some detail. The figure below shows a triangulated surface with
400 thousand vertices. This data set will be used to illustrate the reverse
engineering process.

\image html Tarn_init.gif "Triangulated surface"  width=600px 

\section re_sec3 Workflow
The first function to call is RevEng::enhancePoints. The points are approximated by a surface
in a local neighbourhood. Surface normal and principal curvature estimates are 
computed from this surface. Approximation errors are registered and used to set an
approximation tolerance for the proceeding computations. An additional surface normal is
computed from the triangulation. The two versions of the surface normal have different
pros and cons, and both are used in the computations.

Classification is performed in RevEng::classifyPoints. It is based on the size and
sign of estimated Gauss and mean curvature in the points. Very small curvature values are
deciphered as zero. A small curvature radius compared to the average distance between
triangle vertices indicates that the point is a part of an edge. The expected typical 
measured objects has rounded edges. Thus, edge detection is not a prioritized topic in the current version of the 
functionality. As the initial triangulation may lack smoothness, the curvature
information is somewhat unstable, but still appropriate for recognizing significant 
regions suitable for being represented by one surface. In the image below, pink colour 
indicates that the area around a point is classified as flat, green and orange points
have Gauss curvature close to zero while the remaining colours indicate some 
curvature.

\image html Tarn_classify.gif "Vertices coloured according to curvature properties" width=600px

Next, the approximation tolerance is set by the call 
RevEng::setApproximationTolerance based on information from the preceeding computations. 
Alternatively, the application can use the
function RevEng::setApproxTol(double tol) if more control is preferred. As the given
point cloud is expected to be noisy, it is only required that a majority of the points associated 
to a surface will be fit by the surface within this tolerance. There are also requirements on the average approximation error and surface normal direction.

RevEng::segmentIntoRegions collects connected groups of points with the same classification.
Identified edge points are excluded. This is a two-stage process. First connected points
with the same classification are grouped, then small groups are merged with their
neighbour if feasible. Each group is stored in an instance of \link Go::RevEngRegion RevEngRegion\endlink. In the image below, the individual groups of some size are given 
separate colours. Very small groups are shown in black. 

\image html Tarn_segment.gif "Segmented point cloud" width=800px

The first surface creation is performed in RevEng::initialSurfaces. Regions with a 
significant number of points are selected and tentatively fitted by a plane, 
a cylinder or a cone. As the point classification can be misleading several
attempts are made and the best fit is selected if it satisfies the accuracy
requirements. Simultaneously, points that are found to belong to other surfaces are
dismissed from the region. For our test example, only regions with more than 1357
points are considered for surface creation. This number is estimated from the 
current region sizes. The surfaces are represented as 
\link Go::ParamSurface ParamSurface\endlink and embedded in 
\link Go::Hedgeurface HedgeSurface\endlink, which is inherited from
\link Go::ftEdge ftEdge\endlink. The collection of HedgeSurfaces are owned by 
the RevEng entity, but each surface is linked to the associated RevEngRegion. 
The result of the first surface creation for our test case is shown in the image below.

\image html Tarn_initsurf.gif "The first surfaces, associationed points and main regions" width=800px

The model is still incomplete. Some, even major, surfaces are missing. Transition zones
are missing and the surfaces don't have any matching directions. Another region growing is
performed in RevEng::growSurfaces, this time with the existence of some surfaces where
the distance between the points of a region and the surface of the adjacent region
can be measured. Next an identification of similar plane normals and rotational surface
axis is performed in RevEng::updateAxesAndSurfaces. Surfaces with almost identical 
axes are updated to be consistent. If this implies a major detoriation of the
surface accuracy, the update is dismissed. A coordinate system for the model is
defined.

The first edge recognition is performed in RevEng::firstEdges. Adjacent and almost
adjacent surfaces are intersected and the intersection curve is stored in 
\link Go::RevEngEdge RevEngEdge\endlink along with nearby regions. These regions is 
associated with the blend surface for which the edge is a place holder. The edge
also contains some context information. In the figure below, we see that some 
edges are still missing and that the edges don't join up. The main cylindrical surface
and the middle plane meets in five different edges. One is split at the seam of the
cylinder.

\image html Tarn_edges1.gif "The first edges and associationed points" width=600px

The first sequence of point 4 to 6 in the process overview is completed. Now 
RevEng::surfaceCreation is applied to continue the surface recognition. At this stage,
it is also possible to recognize spheres and tori. As in the first surface
recognition pass, more than one primary surface can be fitted to the points, and
the best fit is chosen if accurate enough. Some regions may be composed by several sub
groups of points that can be associated one surface. The identified model
coordinate system provides a tool to split these regions into consistent parts.
Adjacent surfaces with the same characteristics are merged and the region structure
is simplified further by including small regions into adjacent ones when possible.
The figure below shows the 16 surfaces identified in the example model and the
regions associated to these surfaces.

\image html Tarn_surfaces2.gif "Updated surface structure and associated regions" width=600px

RevEng::adaptToMainAxis updates the axis information obtained by updateAxesAndSurfaces
and complements information on axis direction with axis position. After harmonizing
the surfaces with respect to the updated global information, any potential missing edges are
computed.

At this stage the major surfaces and associated candidate blends are expected to
be indentified but there may still be unidentified smaller surfaces. Curvature  
information and point classification become unstable close to edges in the model and
the restriction on number of points to initiate surface fitting lead to a high
possibility of unrecognized surfaces. In RevEng::smallRegionSurfaces adjacent region 
fragments are collected with the aid of global information. Such merged regions are
fitted with planes, cylinders and cones if applicable. The region structure may
be updated by merging adjacent regions into those associated with the new surfaces and
by merging surfaces, and associated regions, representing the same feature. Additional
edges may be defined.

RevEng::manageBlends1 computed blend surface associated to the collected RevEngEdges.
The blend surfaces will either be cylinders or torii and some coordination of blend
radii is performed. Depending on the actual configuration of the adjacent surfaces,
some gaps or overlaps may occur. The current state for our example model is
visualized in the figure below, first all surfaces
including blends, then the large cylinder is removed for better visibility of the
blend surfaces and finally a detail is emphasized.

\image html Tarn_blend1.gif "Surfaces with edge blends" width=800px

RevEng::manageBlends2 bounds the blend surfaces in the blend direction and defines
corresponding trimming curves represented as \link Go::ftEdge ftEdge\endlink and
stored in the blend region. The
geometry representation of the trimming curves are stored in an instance of
\link Go::CurveOnSurface CurveOnSurface\endlink having
curves both in geometry and parameter space. The corresponding curves in the
surfaces being blended are also computed.

Edge blend surfaces will often meet in a corner blend. Configurations resulting in 
torus corners and some 4-sided free form (spline) corners are supported. Additional
trimming curves are added to the associated regions. The figure below shows, in the
first picture, all surfaces including blends and corner blends. In the following
pictures some surfaces are removed for visibility.

\image html Tarn_blend2.gif "Surfaces with edge and corner blends" width=800px

RevEng::trimSurfaces arranges the identified trimming curves for each surface in
loops and compute bounded surfaces. The figures below illustrate the procedure for
two surfaces in our example case. The points associated to the main cylinder surface
is shown to the left in the first figure along with trimming curves in geometry space.
The curves are extended beyond their real size and must be reduced. This computation
is performed in parameter space. Curves along the cylinder seam is added to the 
collection. Then the trimming curves are cut at intersections with the adjacent 
curves and arranged in a loop. The parameter curves before and after this modification
are shown in the middle two pictures. To the right, the final trimmed surface and the
modified trimming curves in geometry space are shown.

\image html Tarn_trim1.gif "Trimming of cylinder surface" width=800px
\image html Tarn_trim2.gif "Trimming of plane" width=800px

The trimming of the middle planar surface is shown in the figure above. The procedure
is similar to the cylinder case, but in this case we have an internal trimming curve.
This curve is computed as a boundary towards point regions without an associated
surface. The initial triangular surface lacks information in this area, thus the
expected cylindrical surface inside the model is not found. The shape of this
trimming curve is inferior to the ones found by the intersect and blend 
procedure and in newer versions it is dismissed during the trimming operation.



\image html Tarn_surface_collection.gif "All trimmed surfaces" width=600px

The final surface collection is shown in the image above. RevEng::createModel
finalizes the work flow. Adjacency analysis is performed using functionality in
\link Go::SurfaceModel SurfaceModel\endlink and \ref topology. The final model can be
exported in a \link Go::CompositeModelFileHandler g22-file\endlink.

\section re_sec4 Workflow summary - implementation
The workflow described above is implemented through a number of calls to functions in
the reverse engineering workflow RevEng. The function calls are listed below, and all 
must be executed.
- RevEng::Reveng(shared_ptr<ftPointSet> triangulation)
- RevEng::enhancePoints() *
- RevEng::classifyPoints() * 
- RevEng::setApproxTolerance() / RevEng::setApproxTol(double tolerance)
- RevEng::segmentIntoRegions() *
- RevEng::initialSurfaces() *
- RevEng::growSurfaces() *
- RevEng::updateAxesAndSurfaces() *
- RevEng::firstEdges() *
- RevEng::surfaceCreation() *
- RevEng::adaptToMainAxis() *
- RevEng::smallRegionSurfaces() * 
- RevEng::manageBlends1() *
- RevEng::manageBlends2()
- RevEng::trimSurfaces()
- shared_ptr<SurfaceModel> RevEng::createModel()

The workflow is automatic, the only possible current interaction by the application is
to set the tolerance. However, more user interaction is expected to be preferable.
Then user interaction can be used to perform quality control of the type of surfaces recognized and
to enable regularization of the model with respect to for instance parallelity, 
orthogonality and symmetry. This type of interference can in future versions of the
work flow be included between the calls to RevEng functions.

The workflow can be stopped and started at a number of stages. The possibilities
are marked with a star in the list above. To stop the execution after for instance
initialSurfaces and restart at a later stage. To store the required information,
the procedure is:
\verbatim
  RevEng reveng;
  std::ofstream of("initial_surface_stage.txt");
  of << SURFACEONE << std::endl;
  reveng.storeGrownRegions(of);
\endverbatim

To restart:
\verbatim
  std::ifstream is(char *data_file);
  int stage;
  is >> stage;
  reveng.readGrownRegions(is);
  if (stage <= SURFACEONE)
    {
      reveng.growSurfaces();

      // Continue with the remaining functions
    }
\endverbatim

Suggested flags for the stages are:
\verbatim 
  enum
    {
     ENHANCED, CLASSIFIED, SEGMENTED, SURFACEONE, GROWONE, UPDATEONE, EDGESONE, SURFACETWO, UPDATETWO, SURFACETHREE, BLENDONE
    };
\endverbatim
Note that storing and reading the execution stage may be time consuming.

Debug information can be fetched from RevEng after the definition of regions, that
means after executing segmentIntoRegions(), by calling
\verbatim
  void writeRegionStage(std::ostream& of, std::ostream& ofm, std::ostream& ofs) const;
  void writeRegionWithSurf(std::ostream& of) const;
  void writeEdgeStage(std::ostream& of) const;
\endverbatim
writeRegionStage will write points and surfaces organized in regions to three files.
The points of regions with a significant number of points are written to "of" 
together with possible associated surfaces. Regions with less points are written to 
"ofm" and all points from regions with few points are collected and written to "ofs".
writeRegionWithSurf outputs only those regions that has associated points.
writeEdgeStage writes all RevEngEdges and point groups associated to expected blend
surfaces. The debug information can be visualized in \ref viewlib.


\section re_sec5 Data structure
The reverse engineering functionality is integrated in GoTools. The main engine is
\link Go::RevEng RevEng\endlink where the workflow is run. RevEng is also the owner of the
data points, structures where data points are organized and computed surfaces
and edges. An overview is given in the figure below. 

\image html datastructure.png "Data structure" width=600px

The triangulation is represented by an \link Go::ftPointSet ftPointSet\endlink 
stored in RevEng. Segmented points are orgianized in \link Go::RevEngRegion RevEngRegion\endlink. RevEngRegion contains pointers to \link Go::RevEngPoint RevEngPoint\endlink
inherited from \link Go::ftSamplePoint ftSamplePoint\endlink. The points themselves
are still stored in ftPointSet while the point pointers are moved between regions as
the computation proceedes. Computed surfaces are stored in 
\link Go::HedgeSurface HedgeSurface\endlink, which again is stored in RevEng. A
surface is connected to one region, while a region may have a surface connected.
HedgeSurface is inherited from \link Go::ftSurface ftSurface\endlink, which is
input to the adjacency analysis, see \ref topology. A
\link Go::RevEngEdge RevEngEdge\endlink is computed from two surfaces but linked
to the associated regions. The instance itself is stored in RevEng. RevEngEdge also
contains pointers to regions associated to the foreseen blend surface and to the 
blend region collecting these associated regions when the blend surface is created.
\link Go::RevEngUtils RevEngUtils\endlink provides some utility functionality to the 
reverse engineering process.
*/
#endif // _REVERSEENGINEERING_DOXYMAIN_H
