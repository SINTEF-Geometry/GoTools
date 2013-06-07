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

#ifndef PRORGANIZEPOINTS_H
#define PRORGANIZEPOINTS_H

#include "GoTools/utils/Array.h"
#include <vector>
#include "GoTools/utils/config.h"

using Go::Vector3D;
using std::vector;


/*<PrOrganizedPoints-syntax: */

/** PrOrganizedPoints -  This class implements an interface to a planar graph
 * embedded in \f$R^2\f$ or \f$R^3\f$.
 * The two common examples are a triangulation and
 * a topologically rectangular grid, but it can also represent point clouds
 * with no explicit connectivity, as long as the boundary is identified.
 * The methods in this class are used by PrParametrize.
 */
class GO_API PrOrganizedPoints
{
public:
  /** @name Pure virtual functions... */
  //@{

  /// return the number of nodes in the graph. 
  virtual int         getNumNodes() const  = 0;
  /// return the i-th node in the graph if the nodes are three-dimensional
  virtual Vector3D    get3dNode(int i) const  = 0;
  /// set the coordinates of the i-th node in the graph if the nodes are
  /// three-dimensional.
  virtual void        set3dNode(int i, const Vector3D& p) = 0;

  /** Return the indices of the neighbours of the i-th node in:
  * 1. Any anticlockwise order if i is an interior node.
  * 2. The unique anticlockwise order if i is a boundary node.
  * This routine should be implemented in
  * any derived class by first calling clear()
  * followed successively by push_back() on the neighbours vector.
  */
  virtual void        getNeighbours(int i, std::vector<int>& neighbours) const= 0;

  /// is the i-th node a boundary node or interior node?
  virtual bool        isBoundary(int i) const = 0;

  /// return the u parameter value of the i-th node.
  virtual double      getU(int i) const = 0;
  /// return the v parameter value of the i-th node.
  virtual double      getV(int i) const = 0;

  /// reset the u parameter value of the i-th node.
  virtual void        setU(int i, double u) = 0;
  /// reset the v parameter value of the i-th node.
  virtual void        setV(int i, double v) = 0;
  /// Empty default destructor
  virtual ~PrOrganizedPoints(); // @afr: Moved empty def to .C file.
  //@}

  /** @name Other functions */
  //@{

  /// count the number of boundary nodes in the graph.
  virtual int  findNumBdyNodes() const;

  /// count the number of edges in the graph
  virtual int  findNumEdges() const;

  /// Check if there are any elements of the vector 'face' that are 
  /// inferior to the value of 'i'.
  static bool  isMinimum(int i, std::vector<int>& face);


  /// Compute the number of connected components of the graph
  int  findNumComponents() const;

  /// Locate all nodes connnected to node 'i', and mark their corresponding entry
  /// in the vector 'component' with the value 'ic'.
  /// \param i index specifying the node
  /// \param ic the 'label value' to write to the concerned entries in the 
  ///           'component' vector
  /// \retval component the vector where the labels will be written.  Its size should
  ///                   be equal to the number of nodes in the PrOrganizedPoints 
  ///                   object before calling this function.
  void labelNode(int i, int ic, std::vector<int>& component) const;
    
  /// Compute the number of separate boundaries of the graph
  int  findNumBdyComponents() const;

  /// Locate all nodes on the same boundary as the node 'i' (supposedly a boundary 
  /// node) and mark their corresponding entry in the vector 'component' with the
  /// value 'ic'.
  /// \param i index specifying a boundary node
  /// \param ic the 'label value' to write to the concerned entries in the 
  ///           'component' vector
  /// \retval component the vector where the labels will be written.  Its size
  ///                   should be equal to the number of nodes in the PrOrganizedPoints
  ///                   object before calling this function.
  void labelBdyNode(int i, int ic, std::vector<int>& component) const;

  /// Locate all connected components of the graph, give each of them an unique index
  /// (starting from 0), and giving all nodes of each component "local" indexes
  /// pertaining to that component.
  /// \retval component will be resized to number of nodes.  Each entry will contain
  ///                   the index of the component in which the corresponding node is
  ///                   located
  /// \retval newIndex  will be resized to number of nodes.  Each entry will contain
  ///                   the LOCAL index of the corresponding node in the component 
  ///                   it belongs to.
  int  indexComponents(std::vector<int>& component,
                          std::vector<int>& newIndex) const;
  /// Locate all nodes in the same connected component as the node 'i', and label
  /// their corresponding entries in the vector 'component' with the value 'ic'.
  /// These nodes are also given "local" indexes for that coponent.  These indexes
  /// will be written to the corresponding entries of the 'newIndex' vector.
  /// \param i the index of the node specifying the connected component to examine.
  /// \param ic the label to give to the entries in the 'component'-vector that 
  ///           correspond with nodes in the specified connected component.
  /// \param index start value for local indexing of nodes within the component.
  /// \retval component this vector is where the labeling is written.  Before calling
  ///                   the function, the user should ensure that its size is equal
  ///                   to the total number of nodes in the graph.
  /// \retval newIndex  this vector is where the new local coordinates are written.
  ///                   Before calling the function, the user should ensure that its
  ///                   size is equal to the total number of nodes in the graph.
  void labelNode(int i, int ic, int& index,
                         std::vector<int>& component,
                         std::vector<int>& newIndex) const;

  /** find the index of the node whose (x,y,z) point is exactly equal to
   * the given point. This can be useful for relocating "corner" points
   * after triangulating (if the graph is a triangulation).
   */
  int  findIndex(Vector3D& point) const;

  /// Compute for each node the number of edges that must be traversed in order to 
  /// reach the boundary (For boundary edges, this number is 0).
  /// \retval label this vector will contain the computed number for each node.  
  ///               Its size will therefore be equal to the number of nodes.
  void topologicalDistToBdy(std::vector<int>& label) const;

  // /// Compute for each node the shortest Euclidian distance that must be traversed 
  // /// along edges in order to reach the boundary.
  // void geometricalDistToBdy(std::vector<double>& label) const;

  /// Return the common neighbours of the j-th and k-th node in some arbitrary order
  void getCommonNeighbours(int j, int k, std::vector<int>& neighbours) const;
  /// Get the 2-neighbours of the node 'i'.
  void get2Neighbours(int i, std::vector<int>& neighbours) const;
 //@}

  /** @name Print routines */
  //@{

  /// Print out the XYZ nodes of the graph. If num = true, print first the number of nodes.
  virtual void printXYZNodes(std::ostream& os, bool num = 0) const;
  /// Print out the UV nodes of the graph. If num = 1, print first the number of nodes. 
  virtual void printUVNodes(std::ostream& os, bool num = 0) const;
  /// Print out the nodes of the graph. If num = 1, print first the number of nodes. 
  virtual void printUVXYZNodes(std::ostream& os, bool num = 0) const;
  /// Print out the edges of the graph.
  virtual void printXYZEdges(std::ostream& os) const;
  /// Print out the edges of the parametrization.
  virtual void printUVEdges(std::ostream& os) const;
  /// Print general information about this PrOrganizedPoints object (number of nodes, etc.)
  virtual void printInfo(std::ostream& os) const;
 //@}
};

/*>PrOrganizedPoints-syntax: */

/*Class:PrOrganizedPoints

Name:              "PrOrganizedPoints" - interface functions for a planar graph
Syntax:	           @PrOrganizedPoints-syntax
Keywords:
Description:       This class implements an interface to a planar graph
                   embedded in $R^2$ or $R^3$.
                   The two common examples are a triangulation and
                   a topologically rectangular grid.
                   The methods in this class are used by PrParametrize.
Member functions:
                   "getNumNodes()" --\\
                   return the number of nodes in the graph.
                 
                   "get3dNode(int i)" --\\
                   return the i-th node in the graph if the nodes are 
                   three-dimensional

                   "getNeighbours(int i, PrListInt& neighbours)" --\\
                   Return the indices of the neighbours of the i-th node in:
                   1. any anticlockwise order if i is an interior node
                   2. the unique anticlockwise order if i is a boundary node.
                   This routine should be implemented in
                   any derived class by first calling empty()
                   followed successively by append(); see PrIntList.h

                   "isBoundary(int i)" --\\
                   is the i-th node a boundary node or interior node?

                   "getU(int i)" --\\
                   return the u parameter value of the i-th node.

                   "getV(int i)" --\\Face
                   return the v parameter value of the i-th node.

                   "setU(int i, double u)" --\\
                   reset the u parameter value of the i-th node.

                   "setV(int i, double v)" --\\
                   reset the v parameter value of the i-th node.

                   "findNumBdyNodes()" --\\
                   count the number of boundary nodes in the graph.

                   "findIndex(Vector3D&)" --\\
                   find the index of the node whose (x,y,z) point is
                   exactly equal to the given point.
                   This can be useful for relocating "corner" points
                   after triangulating (if the graph is a triangulation).

Constructors:
Files:
Example:
See also:
Developed by:      SINTEF Applied Mathematics, Oslo, Norway
Author:	           Michael Floater, SINTEF
Date:              Mar. 97
*/

#endif // PRORGANIZEDPOINTS_H
