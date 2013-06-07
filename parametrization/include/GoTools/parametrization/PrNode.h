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

#ifndef PRNODE_H
#define PRNODE_H


#include "GoTools/utils/Array.h"
using Go::Vector2D;
using Go::Vector3D;

/** This class represents a node for use in
 * PrTriangulation_OP. It has double x,y,z values contained
 * in its parent class Vector3D. It also has double
 * u and v values (parameter values) and the index
 * of its "first" incident triangle tr. If the node
 * is interior then tr is any incident triangle. Otherwise
 * the node is a boundary node, in which case tr is the
 * incident triangle which is the first one in the
 * anticlockwise direction around the node (equivalently
 * the only triangle containing the node whose right
 * neighbour from the point of view of the node is 0).
 */
class PrNode
{
private:
    Vector3D pt_;
    double u_,v_;
    int tr_;

public:
  PrNode()
      : pt_(), u_(0), v_(0), tr_(0)
	{}
  PrNode(double x, double y, double z, double u, double v, int tr)
      : pt_(x, y, z), u_(u), v_(v), tr_(tr)
	{}


  inline void init(double x, double y, double z, double u, double v, int tr);

  inline const Vector3D& point() const { return pt_; }
  inline const double& u() const;
  inline const double& v() const;
  inline const double& x() const;
  inline const double& y() const;
  inline const double& z() const;
  inline const int& tr() const;


  inline       double& u();
  inline       double& v();
  inline       double& x();
  inline       double& y();
  inline       double& z();
  inline       int& tr();

  // print and scan routines
  void print(std::ostream& os) const;
  void printXYZ(std::ostream& os) const;
  void printUV(std::ostream& os) const;
  void scan(std::istream& is);
};

/*>PrNode-syntax: */

//-----------------------------------------------------------------------------
inline void PrNode::init(double x, double y, double z, double u, double v, int tr)
//-----------------------------------------------------------------------------
{
    pt_[0] = x;
    pt_[1] = y;
    pt_[2] = z;
    u_ = u;
    v_ = v;
    tr_ = tr;
}

//-----------------------------------------------------------------------------
inline const double& PrNode::u() const
//-----------------------------------------------------------------------------
{
  return u_;
}

//-----------------------------------------------------------------------------
inline const double& PrNode::v() const
//-----------------------------------------------------------------------------
{
  return v_;
}

//-----------------------------------------------------------------------------
inline const int& PrNode::tr() const
//-----------------------------------------------------------------------------
{
  return tr_;
}

//-----------------------------------------------------------------------------
inline double& PrNode::u() 
//-----------------------------------------------------------------------------
{
  return u_;
}

//-----------------------------------------------------------------------------
inline double& PrNode::v() 
//-----------------------------------------------------------------------------
{
  return v_;
}

//-----------------------------------------------------------------------------
inline int& PrNode::tr() 
//-----------------------------------------------------------------------------
{
  return tr_;
}

//-----------------------------------------------------------------------------
inline const double& PrNode::x() const
//-----------------------------------------------------------------------------
{
  return pt_[0];
}

//-----------------------------------------------------------------------------
inline       double& PrNode::x()
//-----------------------------------------------------------------------------
{
  return pt_[0];
}

//-----------------------------------------------------------------------------
inline const double& PrNode::y() const
//-----------------------------------------------------------------------------
{
  return pt_[1];
}

//-----------------------------------------------------------------------------
inline       double& PrNode::y()
//-----------------------------------------------------------------------------
{
  return pt_[1];
}

//-----------------------------------------------------------------------------
inline const double& PrNode::z() const
//-----------------------------------------------------------------------------
{
  return pt_[2];
}

//-----------------------------------------------------------------------------
inline       double& PrNode::z()
//-----------------------------------------------------------------------------
{
  return pt_[2];
}

/*Class:PrNode

Name:              PrNode
Syntax:	           @PrNode-syntax
Keywords:
Description:       This class represents a node for use in
                   PrTriangulation_OP. It has double x,y,z values contained
                   in its parent class Vector3D. It also has double
                   u and v values (parameter values) and the index
                   of its "first" incident triangle tr. If the node
                   is interior then tr is any incident triangle. Otherwise
                   the node is a boundary node, in which case tr is the
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
Date:              August 97
*/

#endif // PRNODE_H
