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

#ifndef PRTRIANGLE_H
#define PRTRIANGLE_H

#include <iostream>

/*<PrTriangle-syntax: */

/** PrTriangle -  This class represents a triangle for use in
 * PrTriangulation_OP. It has three indices of its
 * three vertices in some anticlockwise order.
 * Its three neighbouring triangles t1,t2,t3 are opposite
 * to its vertices. Thus t1 is the triangle opposite n1,
 * t2 is opposite n2 and t3 opposite n3. If any of the
 * neighbouring triangles do not exist in the triangulation
 * their index is -1 (which is out of range in the array
 * of triangles in the triangulation).   
 */
class PrTriangle
{
private:
  int n1_,n2_,n3_;
  int t1_,t2_,t3_;

public:

  /// Default constructor
  PrTriangle() {}
  /// Constructor
  PrTriangle(int n1, int n2, int n3, int t1, int t2, int t3)
    { n1_ = n1; n2_ = n2; n3_ = n3; t1_ = t1; t2_ = t2; t3_ = t3; }
 /// Constructor
  PrTriangle(const PrTriangle& t)
     {n1_ = t.n1_; n2_ = t.n2_; n3_ = t.n3_;
      t1_ = t.t1_; t2_ = t.t2_; t3_ = t.t3_; }
  /// Empty destructor
  ~PrTriangle() {}

  inline void init(int n1, int n2, int n3, int t1, int t2, int t3);

  inline const int& n1() const;
  inline const int& n2() const;
  inline const int& n3() const;
  inline const int& t1() const;
  inline const int& t2() const;
  inline const int& t3() const;

  inline int& n1();
  inline int& n2();
  inline int& n3();
  inline int& t1();
  inline int& t2();
  inline int& t3();

  /// Return the neighbouring triangle which lies opposite to the given node.
  int  getOppositeTriangle(int node) const;

  /// Return the neighbouring triangle which lies to the left of the given node. 
  int  getLeftTriangle(int node) const;

  /// Return the neighbouring triangle which lies to the right of the given node.
  int  getRightTriangle(int node) const;

  /// Return the node following "node" in an anticlockwise direction around
  /// the triangle.
  int  getAnticlockwiseNode(int node) const;

  /// Return the node following "node" in a clockwise direction around the triangle.
  int  getClockwiseNode(int node) const;

  /// Boolean routine: does "node" belong to the triangle?
  bool  isVertex(int node) const;

  /// If 'n1' is a node in the triangle, replace it with 'n2'.
  void replaceNode(int n1, int n2);

  /// If 't1' is a triangle linked to this triangle, replace it with 't2'.
  void replaceTriangle(int t1, int t2);

  /// if 'triangle' share an edge with this triangle, return the nodes on the 
  /// shared edge.
  int getEdge(int triangle, int &n1, int &n2) const;

  /// print contents to stream
  void print(std::ostream& os) const;

  /// read contents from stream
  void scan(std::istream& is);
};


/*>PrTriangle-syntax: */

//-----------------------------------------------------------------------------
inline void PrTriangle::init(int n1, int n2, int n3, int t1, int t2, int t3)
//-----------------------------------------------------------------------------
{
  n1_ = n1;
  n2_ = n2;
  n3_ = n3;
  t1_ = t1;
  t2_ = t2;
  t3_ = t3;
}

//-----------------------------------------------------------------------------
inline const int& PrTriangle::n1() const
//-----------------------------------------------------------------------------
{
  return n1_;
}

//-----------------------------------------------------------------------------
inline const int& PrTriangle::n2() const
//-----------------------------------------------------------------------------
{
  return n2_;
}

//-----------------------------------------------------------------------------
inline const int& PrTriangle::n3() const
//-----------------------------------------------------------------------------
{
  return n3_;
}

//-----------------------------------------------------------------------------
inline const int& PrTriangle::t1() const
//-----------------------------------------------------------------------------
{
  return t1_;
}

//-----------------------------------------------------------------------------
inline const int& PrTriangle::t2() const
//-----------------------------------------------------------------------------
{
  return t2_;
}

//-----------------------------------------------------------------------------
inline const int& PrTriangle::t3() const
//-----------------------------------------------------------------------------
{
  return t3_;
}

//-----------------------------------------------------------------------------
inline int& PrTriangle::n1()
//-----------------------------------------------------------------------------
{
  return n1_;
}

//-----------------------------------------------------------------------------
inline int& PrTriangle::n2()
//-----------------------------------------------------------------------------
{
  return n2_;
}

//-----------------------------------------------------------------------------
inline int& PrTriangle::n3()
//-----------------------------------------------------------------------------
{
  return n3_;
}

//-----------------------------------------------------------------------------
inline int& PrTriangle::t1()
//-----------------------------------------------------------------------------
{
  return t1_;
}

//-----------------------------------------------------------------------------
inline int& PrTriangle::t2()
//-----------------------------------------------------------------------------
{
  return t2_;
}

//-----------------------------------------------------------------------------
inline int& PrTriangle::t3()
//-----------------------------------------------------------------------------
{
  return t3_;
}

/*Class:PrTriangle

Name:              PrTriangle
Syntax:	           @PrTriangle-syntax
Keywords:
Description:       This class represents a triangle for use in
                   PrTriangulation_OP. It has three indices of its
                   three vertices in some anticlockwise order.
                   Its three neighbouring triangles t1,t2,t3 are opposite
                   to its vertices. Thus t1 is the triangle opposite n1,
                   t2 is opposite n2 and t3 opposite n3. If any of the
                   neighbouring triangles do not exist in the triangulation
                   their index is -1 (which is out of range in the array
                   of triangles in the triangulation).
Member functions:
                   "getOppositeTriangle(int node)" --\\
                   Return the neighbouring triangle which lies
                   opposite to the given node.
                 
                   "getLeftTriangle(int node)" --\\
                   Return the neighbouring triangle which lies
                   to the left of the given node.
                 
                   "getRightTriangle(int node)" --\\
                   Return the neighbouring triangle which lies
                   to the right of the given node.
                 
                   "getAnticlockwiseNode(int node)" --\\
                   Return the node following "node" in an anticlockwise
                   direction around the triangle.
                 
                   "getClockwiseNode(int node)" --\\
                   Return the node following "node" in a clockwise
                   direction around the triangle.
                 
                   "isVertex(int node)" --\\
                   Boolean routine: does "node" belong to the triangle?
                 
Constructors:
Files:
Example:
See also:
Developed by:      SINTEF Applied Mathematics, Oslo, Norway
Author:	           Michael Floater, SINTEF
Date:              August 97
*/

#endif // PRTRIANGLE_H
