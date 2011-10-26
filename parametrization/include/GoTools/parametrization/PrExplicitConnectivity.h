//===========================================================================
//                                                                           
// File: PrExplicitConnectivity.h                                            
//                                                                           
// Created: Wed Mar 29 16:40:21 2006                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: PrExplicitConnectivity.h,v 1.3 2007-03-02 16:18:19 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

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

