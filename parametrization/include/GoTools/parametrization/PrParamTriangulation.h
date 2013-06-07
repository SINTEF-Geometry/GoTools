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

#ifndef PRPARAMTRIANGULATION_H
#define PRPARAMTRIANGULATION_H

#include "GoTools/parametrization/PrTriangulation_OP.h"
#include "GoTools/parametrization/PrSubTriangulation.h"
#include <vector>
using std::vector;
using std::istream;
using std::ostream;

/// This class represents a fine triangulation embedded in
/// a coarser one.  The points on the faces on the coarser 
/// triangulation can be seen as parametrizing the points on
/// the finer triangulation.
class PrParamTriangulation {

  PrTriangulation_OP*    mesh_;
  PrTriangulation_OP*    basemesh_;
  vector<int>            baseTriangle_;
  vector< vector<int> >  patch_triangles_;

public: 
  /// Empty default constructor
  PrParamTriangulation() {}

  /// Empty default destructor
  virtual ~PrParamTriangulation() {}

  /// Read from stream
  void open (istream& is);

  /// Write to stream
  void save (ostream& os);

  /// Get a pointer to the fine mesh
  PrTriangulation_OP* getFineMesh () { return mesh_; }

  /// Get a pointer to the coarse mesh
  PrTriangulation_OP* getCoarseMesh () { return basemesh_; }

  /// Reset this objects to handle the fine/coarse triangulation pair
  /// represented by 'mesh' and 'basemesh'.
  void attach(PrTriangulation_OP* mesh, PrTriangulation_OP* basemesh);

  /// Set all triangles contained in a sub-triangulation of the fine
  /// mesh to "correspond" to a certain triangle in the coarse mesh.
  /// \param idx index to the triangle in the coarse mesh
  /// \param sub_tri a sub-triangulation that represents those triangles
  ///                in the fine mesh that are to be associated with
  ///                the 'idx' triangle in the coarse mesh.
  void makeCorrespondences (int idx, PrSubTriangulation& sub_tri);

  /// (Print out debug information)
  void printBaseTriangles ();

  /// Get the 'parameter point' of a certain node on the fine mesh.  By
  /// this we mean to find its "corresponding" point on the coarse mesh.
  /// \param i index to the node on the fine mesh
  /// \return a 3D point on the coarse mesh that corresponds to this node
  Vector3D getParamPoint(int i);

  /// Given a triangle on the coarse mesh and the barycentric coordinates
  /// of a point in this triangle, find the triangle on the fine mesh that
  /// contains this point.
  /// \param coarseT index of the triangle on the coarse mesh
  /// \param coarseBC barycentric coordinates of the point with respect to
  ///                 the triangle on the coarse mesh
  /// \retval fineT gives the index to the triangle on the fine mesh
  ///               containing the specified point
  /// \retval fineBC gives the barycentric coordinates of the point with 
  ///                respect to the triangle on the fine mesh.
  void findTriangleContainingPoint(int coarseT, const Vector3D& coarseBC,
				   int& fineT, Vector3D& fineBC);

  /// Returns the barycentric coordinates of node with respect to a triangle
  /// on the fine mesh.   Usually, these are the uv-values obtained by the
  /// parameterization, but the handling becomes more complicated on edges.
  /// \param node index to a node on the fine mesh, for which we want to 
  ///             determine the barycentric coordinates
  /// \param tri index to a base triangle (ie. triangle on the coarse mesh),
  ///            for with we want to find the barycentric coordinates of the
  ///            node given by 'node'.
  /// \retval u the u-coordinate of the node
  /// \retval v the v-coordinate of the node (the w-coordinate is given by 1-u-v).
  void getUV (int node, int tri, double&u, double&v);

  /// Determine whether a given triangle on the fine mesh contains a point
  /// specified by its barycentric coordinates on the underlying coarse
  /// triangle.
  /// \param fineTri index of the triangle on the fine mesh
  /// \param coarseTri index of the triangle on the coarse mesh
  /// \param coarseBC barycentric coordinates of the point on the triangle
  ///                 in the coarse mesh
  /// \retval fineBC if the triangle on the fine mesh was found to contain
  ///                the point in question, its barycentric coordinates with
  ///                respect to this triangle will be returned here
  /// \return 'true' if the specified point is contained n the triangle on
  ///                the fine mesh, 'false' otherwise.
  bool triangleContainsPoint(int fineTri, int coarseTri,
			     const Vector3D& coarseBC, 
			     Vector3D& fineBC);

  /// Given a base triangle (triangle on the coarse mesh) and barycentric
  /// coordinates with respect to this triangle, compute the corresponding
  /// point on the fine mesh.
  /// \param coarseTri index of the base triangle in question
  /// \param coarseBC the barycentric coordinates
  /// \return a 3D point on the fine mesh corresponding to (parameterized by)
  ///         the point on the coarse mesh specified by 'coarseTri' and 
  ///         'coarseBC'.
  Vector3D getSurfPoint(int coarseTri, const Vector3D& coarseBC);

  /// Evaluate a point on a given mesh by specifying a triangle on the mesh
  /// and the barycentric coordinates within this triangle.  The point is 
  /// computed by linear interpolation of the triangle's corner nodes.
  static Vector3D evaluator(PrTriangulation_OP& mesh, int idx, Vector3D& bc);
};

#endif // PRPARAMTRIANGULATION_H

