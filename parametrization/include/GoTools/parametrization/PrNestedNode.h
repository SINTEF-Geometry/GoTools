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

#ifndef PRNESTEDNODE_H
#define PRNESTEDNODE_H

#include "GoTools/utils/Array.h"
using Go::Vector2D;
using Go::Vector3D;
#include <vector>
using std::vector;
using std::ostream;

/*<PrNestedNode-syntax: */

/** This class represents a node for use in
 * PrLevelTriangulation_OP and PrNestedTriangulation.
 * It has double x,y,z values and double u and v values (parameter
 * values). It also has the index of its "first" incident triangle 
 * at each triangulation level. If the node
 * is interior then tr[j] is any incident triangle. Otherwise
 * the node is a boundary node, in which case tr[j] is the
 * incident triangle which is the first one in the
 * anticlockwise direction around the node (equivalently
 * the only triangle containing the node whose right
 * neighbour from the point of view of the node is 0).
 */
class PrNestedNode // : public Vector3D
{
private:
  Vector3D pnt_;
  double u_,v_;
  int level_;      // 0 \le level_ \le k.
  vector<int> tr_; // array of pointers to triangles, one for each
                   // level: level_, ... up to k.
                   // Thus k  = level_ + tr_.size() - 1.

public:
  /// Empty default constructor
  PrNestedNode() {}
  /// Constructor
  PrNestedNode(double x, double y, double z,
	       double u, double v, int level)
            : pnt_(x, y, z) {u_ = u; v_ = v; level_ = level; }
//  PrNestedNode(const PrNestedNode& p)
//            : Vector3D(p) {u_ = p.u(); v_ = p.v(); level_ = p.level_; }
  /// Empty default destructor
  ~PrNestedNode() {}

  inline void init(double x, double y, double z,
		   double u, double v, int level);

  inline const Vector3D& pnt() const;
  inline const double& x() const;
  inline const double& y() const;
  inline const double& z() const;
  inline const double& u() const;
  inline const double& v() const;
  inline const int& level() const;
  inline const int& tr(int i) const;
  // Index access. Check for valid i (0, 1, 2)?
    // const double& operator [] (int i) const { return pnt_[i]; }

  inline       Vector3D& pnt();
  inline       double& x();
  inline       double& y();
  inline       double& z();
  inline       double& u();
  inline       double& v();
  inline       int& level();
  inline       int& tr(int i);
    //  double& operator [] (int i) { return pnt_[i]; }

  void addTrianglePtr(int t) {tr_.push_back(t); }

  // print and scan routines
  void print(ostream& os);
  void printXYZ(ostream& os);
  void printUV(ostream& os);
};

/*>PrNestedNode-syntax: */

//-----------------------------------------------------------------------------
inline void PrNestedNode::init(double x, double y, double z, 
			       double u, double v, int level)
//-----------------------------------------------------------------------------
{
  double p[3];
  p[0] = x;
  p[1] = y;
  p[2] = z;
  pnt_.setValue(p);
  u_ = u;
  v_ = v;
  level_ = level; 
}

//-----------------------------------------------------------------------------
inline const Vector3D& PrNestedNode::pnt() const
//-----------------------------------------------------------------------------
{
  return pnt_;
}

//-----------------------------------------------------------------------------
inline const double& PrNestedNode::x() const
//-----------------------------------------------------------------------------
{
  return pnt_.x();
}

//-----------------------------------------------------------------------------
inline const double& PrNestedNode::y() const
//-----------------------------------------------------------------------------
{
  return pnt_.y();
}

//-----------------------------------------------------------------------------
inline const double& PrNestedNode::z() const
//-----------------------------------------------------------------------------
{
  return pnt_.z();
}

//-----------------------------------------------------------------------------
inline const double& PrNestedNode::u() const
//-----------------------------------------------------------------------------
{
  return u_;
}

//-----------------------------------------------------------------------------
inline const double& PrNestedNode::v() const
//-----------------------------------------------------------------------------
{
  return v_;
}

//-----------------------------------------------------------------------------
inline const int& PrNestedNode::level() const
//-----------------------------------------------------------------------------
{
  return level_;
}

//-----------------------------------------------------------------------------
inline const int& PrNestedNode::tr(int i) const
//-----------------------------------------------------------------------------
{
   if (((unsigned int) i-level_)>=tr_.size())
      return (tr_.back());
   else
      return tr_[i-level_];
}

//-----------------------------------------------------------------------------
inline Vector3D& PrNestedNode::pnt() 
//-----------------------------------------------------------------------------
{
  return pnt_;
}

//-----------------------------------------------------------------------------
inline double& PrNestedNode::x() 
//-----------------------------------------------------------------------------
{
    return pnt_.x();
}

//-----------------------------------------------------------------------------
inline double& PrNestedNode::y() 
//-----------------------------------------------------------------------------
{
  return pnt_.y();
}

//-----------------------------------------------------------------------------
inline double& PrNestedNode::z() 
//-----------------------------------------------------------------------------
{
  return pnt_.z();
}

//-----------------------------------------------------------------------------
inline double& PrNestedNode::u() 
//-----------------------------------------------------------------------------
{
  return u_;
}

//-----------------------------------------------------------------------------
inline double& PrNestedNode::v() 
//-----------------------------------------------------------------------------
{
  return v_;
}

//-----------------------------------------------------------------------------
inline int& PrNestedNode::level() 
//-----------------------------------------------------------------------------
{
  return level_;
}

//-----------------------------------------------------------------------------
inline int& PrNestedNode::tr(int i) 
//-----------------------------------------------------------------------------
{
   if (((unsigned int)i-level_)>=tr_.size())
      return (tr_.back());
   else
      return tr_[i-level_];
}

/*Class:PrNestedNode

Name:              PrNestedNode
Syntax:	           @PrNestedNode-syntax
Keywords:
Description:       This class represents a node for use in
                   PrLevelTriangulation_OP and PrNestedTriangulation.
                   It has double x,y,z values contained
                   in its parent class Vector3D and double
                   u and v values (parameter values). It also has
                   the index of its "first" incident triangle at each
                   triangulation level. If the node
                   is interior then tr[j] is any incident triangle. Otherwise
                   the node is a boundary node, in which case tr[j] is the
                   incident triangle which is the first one in the
                   anticlockwise direction around the node (equivalently
                   the only triangle containing the node whose right
                   neighbour from the point of view of the node is 0).
Member functions:
                 
Constructors:
Files:
Example:
See also:
Developed by:      SINTEF Applied Mathematics, Oslo, Norway
Author:	           Michael Floater, SINTEF
Date:              Nov 2000
*/

#endif // PRNESTEDNODE_H
