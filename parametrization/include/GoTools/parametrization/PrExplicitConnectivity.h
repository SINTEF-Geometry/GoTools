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

#ifndef _PREXPLICITCONNECTIVITY_H
#define _PREXPLICITCONNECTIVITY_H

#include "GoTools/parametrization/PrOrganizedPoints.h"

/*<PrExplicitConnectivity-syntax: */

/** PrExplicitConnectivity -  This class implements an interface to a planar graph
 * embedded in \f$R^2\f$ or \f$R^3\f$, where connectivity is explicitly defined.
 * The two common examples are a triangulation and
 * a topologically rectangular grid.  What differentiates this class from 
 * PrOrganizedPoints is that for PrExplicitConnectivity, it is possible to 
 * directly identify faces.
 * The methods in this class are used by PrParametrize.
 */
class PrExplicitConnectivity : public PrOrganizedPoints
{
public:

  /// count the number of faces in the graph
  virtual int  findNumFaces() const;

  /// Find the face of the graph which
  /// lies to the left of the directed edge whose end points
  /// are the node i and its j-th neighbour.
  /// The nodes of the face, starting with i, will be filled out in the
  /// list vector<int> face. The vector<int> neighbours is used for temporarily
  /// storing neighbours (for efficiency when calling this routine many times).
  /// \param i the index of the node specifying the directed edge
  /// \param j the 'j'th neighbour of 'i' will be the second node specifying the edge
  /// \retval neighbours a temporary vector holding the neighbours of the last node 
  ///         examined by the implementation (might not always be needed; its 
  ///         presence in the parameter list is mainly an efficiency issue).
  /// \retval face upon function return, will contain the indexes of the nodes
  ///         in the requested face (anticlockwise order).
  void findFace(int i, int j, std::vector<int>& neighbours, std::vector<int>& face) const;

  /// Compute the genus of the planar graph
  int  findGenus() const;

  /// Check if the planar graph is a triangulation
  bool isTriangulation() const;

  /// Print out the faces of the graph
  virtual void printXYZFaces(std::ostream& os) const;
  /// Print out the faces of the graph, ML format
  virtual void printXYZFacesML(std::ostream& os) const;
  /// Print out the faces of the graph, VRML format
  virtual void printXYZFacesVRML(std::ostream& os) const;
  /// Print out the faces of the parameterization
  virtual void printUVFaces(std::ostream& os) const;

  /// Find the next edge in the face which lies to the left of the directed
  /// edge whose end points are the node 'i' and its j-th neighbour.
  /// \param i the index of the node specifying the directed edge.  Upon function
  ///          return, it will contain the corresponding value for the next edge in the
  ///          face.
  /// \param j the 'j'th neighbour of 'i' will be the second node specifying the edge.
  ///          Upon function return, it will contain the contain the corresponding
  ///          value for the next edge in the face.
  /// \param neighbours when calling the function, this vector should contain the
  ///                   indexes of the neighbour nodes of node 'i'.  Upon function
  ///                   return, it will contain the corresponding value for the next
  ///                   node in the face.
  void findNextEdgeInFace(int& i, int& j, std::vector<int>& neighbours) const;


  virtual void printInfo(std::ostream& os) const;
  virtual void printTexture(std::ostream& os) const;
};



#endif // _PREXPLICITCONNECTIVITY_H

