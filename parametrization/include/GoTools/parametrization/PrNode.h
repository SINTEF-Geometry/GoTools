/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1997 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

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
