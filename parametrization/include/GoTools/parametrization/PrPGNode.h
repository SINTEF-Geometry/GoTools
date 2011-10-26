/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1997 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

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
