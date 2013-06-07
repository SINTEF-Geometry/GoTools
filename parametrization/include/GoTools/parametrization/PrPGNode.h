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

#ifndef PRPGNODE_H
#define PRPGNODE_H

#include "GoTools/utils/Array.h"
using Go::Vector3D;

/** PrPGNode - This class represents a node for use in
 * PrPlanarGraph_OP. It has double x,y,z coordinate values. 
 * It also has double u and v values (parameter values) and the 
 * index of its last neighbour in the adj_ array of 
 * PrPlanarGraph_OP.
 */
class PrPGNode 
{
private:
  Vector3D pnt_;
  double u_,v_;
  int end_;

public:
  /// Empty default constructor
  PrPGNode() {}
  /// Constructor
  PrPGNode(double x, double y, double z, double u, double v, int end)
            : pnt_(x, y, z) {u_ = u; v_ = v; end_ = end; }
  /// Copy constructor
  PrPGNode(const PrPGNode& p)
	: pnt_(p.x(), p.y(), p.z()) {u_ = p.u(); v_ = p.v(); end_ = p.end(); }
  /// Empty destructor
  ~PrPGNode() {}
  
    //   Vector3D operator Vector3D() const
    //  { return Vector3D(begin()); }

  inline void init(double x, double y, double z, double u, double v, int end);

  inline const Vector3D& pnt() const;
  inline const double& x() const;
  inline const double& y() const;
  inline const double& z() const;
  inline const double& u() const;
  inline const double& v() const;
  inline const int& end() const;

  inline       Vector3D& pnt();
  inline       double& x();
  inline       double& y();
  inline       double& z();
  inline       double& u();
  inline       double& v();
  inline       int& end();

  // print and scan routines
  void print(std::ostream& os);
  void printXYZ(std::ostream& os);
  void printUV(std::ostream& os);
  void scan(std::istream& is);
};


/*>PrPGNode-syntax: */

//-----------------------------------------------------------------------------
inline void PrPGNode::init(double x, double y, double z, double u, double v, int end)
//-----------------------------------------------------------------------------
{
  double p[3];
  p[0] = x;
  p[1] = y;
  p[2] = z;
  pnt_.setValue(p);
  u_ = u;
  v_ = v;
  end_ = end;
}

//-----------------------------------------------------------------------------
inline const Vector3D& PrPGNode::pnt() const
//-----------------------------------------------------------------------------
{
  return pnt_;
}

//-----------------------------------------------------------------------------
inline const double& PrPGNode::x() const
//-----------------------------------------------------------------------------
{
    //return pnt_.x();
    return pnt_[0];
}

//-----------------------------------------------------------------------------
inline const double& PrPGNode::y() const
//-----------------------------------------------------------------------------
{
    //return pnt_.y();
    return pnt_[1];
}

//-----------------------------------------------------------------------------
inline const double& PrPGNode::z() const
//-----------------------------------------------------------------------------
{
    // return pnt_.z();
    return pnt_[2];
}

//-----------------------------------------------------------------------------
inline const double& PrPGNode::u() const
//-----------------------------------------------------------------------------
{
  return u_;
}

//-----------------------------------------------------------------------------
inline const double& PrPGNode::v() const
//-----------------------------------------------------------------------------
{
  return v_;
}

//-----------------------------------------------------------------------------
inline const int& PrPGNode::end() const
//-----------------------------------------------------------------------------
{
  return end_;
}

//-----------------------------------------------------------------------------
inline Vector3D& PrPGNode::pnt() 
//-----------------------------------------------------------------------------
{
  return pnt_;
}

//-----------------------------------------------------------------------------
inline double& PrPGNode::x() 
//-----------------------------------------------------------------------------
{
    //return pnt_.x();
    return pnt_[0];
}

//-----------------------------------------------------------------------------
inline double& PrPGNode::y() 
//-----------------------------------------------------------------------------
{
    //return pnt_.y();
    return pnt_[1];
}

//-----------------------------------------------------------------------------
inline double& PrPGNode::z() 
//-----------------------------------------------------------------------------
{
    //return pnt_.z();
    return pnt_[2];
}

//-----------------------------------------------------------------------------
inline double& PrPGNode::u() 
//-----------------------------------------------------------------------------
{
  return u_;
}

//-----------------------------------------------------------------------------
inline double& PrPGNode::v() 
//-----------------------------------------------------------------------------
{
  return v_;
}

//-----------------------------------------------------------------------------
inline int& PrPGNode::end() 
//-----------------------------------------------------------------------------
{
  return end_;
}

/*Class:PrPGNode

Name:              PrPGNode
Syntax:	           @PrPGNode-syntax
Keywords:
Description:       This class represents a node for use in
                   PrPlanarGraph_OP. It has double x,y,z values contained
                   in its parent class Vector3D. It also has double
                   u and v values (parameter values) and the index
                   of its last neighbour in the adj_ array of
                   PrPlanarGraph_OP.
Member functions:
                 
Constructors:
Files:
Example:
See also:
Developed by:      SINTEF Applied Mathematics, Oslo, Norway
Author:	           Michael Floater, SINTEF
Date:              August 97
*/

#endif // PRPGNODE_H
