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

#include "GoTools/parametrization/PrParametrizeBdy.h"
#include "GoTools/utils/errormacros.h"
//#include "GoTools/utils/Values.h"
#ifdef _WIN32
#define M_PI 3.14159265358979
#endif

#ifdef __BORLANDC__
using std::cos;
using std::sin;
using std::sqrt;
#endif

// PRIVATE METHODS

//-----------------------------------------------------------------------------
int PrParametrizeBdy::findBdyNode()
//-----------------------------------------------------------------------------
// Find the bounday node with the least index.
{
  for(int i=0; i < g_->getNumNodes(); i++)
  {
    if(g_->isBoundary(i)) return i;
  }
  MESSAGE("Could not find boundary node, returning node -1");
  return -1;
}

//-----------------------------------------------------------------------------
int
PrParametrizeBdy::getNextBdyNode(int i)
//-----------------------------------------------------------------------------
// Given a boundary node i, return the next boundary node
// in an anticlockwise direction around the boundary.
{
  g_->getNeighbours(i,neighbours_);
  return neighbours_[0];
}

//-----------------------------------------------------------------------------
double
PrParametrizeBdy::chord(const Vector3D& a, const Vector3D& b)
//-----------------------------------------------------------------------------
//   Return "length" of chord between two points in R^3.
{
  switch(bdyparamtype_)
  {
    case PrCHORDLENGTHBDY: return a.dist(b);

    case PrCENTRIPETAL: return sqrt(a.dist(b));

    case PrUNIFBDY: return 1.0;
  }
  return 0.0;
}

//-----------------------------------------------------------------------------
double
PrParametrizeBdy::boundaryLength(int i1, int i2)
//-----------------------------------------------------------------------------
//   Calculate the "length" of the section of the xyz boundary
//   between nodes boundary i and j according to bdyparamtype_.
//   M.F. Mar. 97.
{
  int i,j;
  i = getNextBdyNode(i1);
  double tot_length = chord(g_->get3dNode(i1),g_->get3dNode(i));

  while(i != i2)
  {
    j = getNextBdyNode(i);
    tot_length += chord(g_->get3dNode(i),g_->get3dNode(j));
    i = j;
  }

  return tot_length;
}


// PUBLIC METHODS

//-----------------------------------------------------------------------------
PrParametrizeBdy::PrParametrizeBdy()
//-----------------------------------------------------------------------------
{
  bdyparamtype_ = PrCHORDLENGTHBDY;
}
//-----------------------------------------------------------------------------
PrParametrizeBdy::~PrParametrizeBdy()
//-----------------------------------------------------------------------------
{
}

//-----------------------------------------------------------------------------
void PrParametrizeBdy::attach(shared_ptr<PrOrganizedPoints> graph)
//-----------------------------------------------------------------------------
{
  g_ = graph;
}

//-----------------------------------------------------------------------------
bool
PrParametrizeBdy::parametrize()
//-----------------------------------------------------------------------------
//   Parametrize the boundary nodes of the given planar graph,
//   according to bdyparamtype_. The boundary parameter points
//   will be placed along the unit circle.
//   M.F. Mar. 97.
{
  int i1 = findBdyNode();
  if(i1 == -1) return false;
  double tot_length = boundaryLength(i1,i1);
  if(tot_length == 0.0) tot_length = 1.0;
    // This can only happen if all
    // the xyz boundary points are equal.

#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
    const double PI = 3.14159265358979323846264338327950288419716939937510;
    double factor = 2.0 * PI / tot_length;
#else
    double factor = 2.0 * M_PI / tot_length;
#endif

  g_->setU(i1,1.0);
  g_->setV(i1,0.0);
  int i,j;
  i = getNextBdyNode(i1);
  double acc_length = chord(g_->get3dNode(i1),g_->get3dNode(i));

  while(i != i1)
  {
    double theta = acc_length * factor;
    g_->setU(i,cos(theta));
    g_->setV(i,sin(theta));
    j = getNextBdyNode(i);
    acc_length += chord(g_->get3dNode(i),g_->get3dNode(j));
    i = j;
  }
  return true;
}

//-----------------------------------------------------------------------------
bool
PrParametrizeBdy::parametrizeSide(int i1, int i2)
//-----------------------------------------------------------------------------
//   Parametrize the boundary nodes between two boundary nodes i1,i2
//   of a given planar graph by chord length.
//   The two nodes i1 and i2 are assumed to be have u and v values.
//   The boundary parameter points between i1 and i2, starting from i1
//   in an anticlockwise direction around the boundary,
//   will be placed along the straight line between the two parameter
//   points.
//   M.F. Mar. 97.
{
  if (i1 == i2) {
    // Suppose edge is degenerate, must be represented by two separate points.
//     MESSAGE("Seems like input corner are identical!");
    THROW("Input corners are identical!");
    return false;
  }

  g_->getNeighbours(i1,neighbours_);
  if(neighbours_[0] == i2) return true;
       // i1 and i2 are adjacent, no nodes in between, nothing to do

  double tot_length = boundaryLength(i1,i2);
  if(tot_length == 0.0) return false;//tot_length = 1.0;
    // This can only happen if all
    // the xyz boundary points are equal between i1 and i2.

  int i,j;
  i = getNextBdyNode(i1);
  double acc_length = chord(g_->get3dNode(i1),g_->get3dNode(i));

  while(i != i2)
  {
    double lambda = acc_length / tot_length;
    double mu     = (tot_length - acc_length) / tot_length;
    g_->setU(i, mu * g_->getU(i1) + lambda * g_->getU(i2));
    g_->setV(i, mu * g_->getV(i1) + lambda * g_->getV(i2));
    j = getNextBdyNode(i);
    acc_length += chord(g_->get3dNode(i),g_->get3dNode(j));
    i = j;
  }
  return true;
}

//-----------------------------------------------------------------------------
bool
PrParametrizeBdy::parametrize(int c1, int c2, int c3, int c4,
                              double umin, double umax, double vmin, double vmax)
//-----------------------------------------------------------------------------
//   Parametrize the boundary nodes of the given planar graph,
//   according to bdyparamtype_. The boundary parameter points
//   will be placed along rectangle.
//   M.F. Mar. 97,
//   Revised MF, Nov. 98.
{
  if(!(g_->isBoundary(c1) && g_->isBoundary(c2) &&
       g_->isBoundary(c3) && g_->isBoundary(c4))) return false;

  g_->setU(c1, umin);
  g_->setV(c1, vmin);
  g_->setU(c2, umax);
  g_->setV(c2, vmin);
  g_->setU(c3, umax);
  g_->setV(c3, vmax);
  g_->setU(c4, umin);
  g_->setV(c4, vmax);

  bool b[4];
  b[0] = parametrizeSide(c1,c2);
  b[1] = parametrizeSide(c2,c3);
  b[2] = parametrizeSide(c3,c4);
  b[3] = parametrizeSide(c4,c1);
  int numnondegenerate = 0;
  for (int i = 0; i < 4; ++i) {
      if (b[i]) ++numnondegenerate;
  }
  return (numnondegenerate >= 2);
}

//-----------------------------------------------------------------------------
void
PrParametrizeBdy::findCornersFromXYZ(int* c)
//-----------------------------------------------------------------------------
//   Find four corner nodes c[0],c[1],c[2],c[3] by taking c[0] to be
//   an arbitrary boundary node and choosing c[1],c[2],c[3] so that
//   the lengths of four boundary curves delimited by c[0],c[1],c[2],c[3]
//   are as equal in "chord" length as possible.
//   The array c should be allocated outside with length 4.
//   M.F. Aug. 97.
{
  c[0] = findBdyNode();
  if(c[0] < 0) return; 
  double quarter_length = boundaryLength(c[0],c[0]) * 0.25;

  int i,j;
  i = getNextBdyNode(c[0]);
  double acc_length = chord(g_->get3dNode(c[0]),g_->get3dNode(i));

  for(int k=1; k<4; k++)
  {
    while(acc_length < quarter_length)
    {
      j = getNextBdyNode(i);
      acc_length += chord(g_->get3dNode(i),g_->get3dNode(j));
      i = j;
    }
    c[k] = i;
    acc_length -= quarter_length;
  }
  return;
}

//-----------------------------------------------------------------------------
void
PrParametrizeBdy::findCornersFromUV(int* c)
//-----------------------------------------------------------------------------
//   This is a simple routine which finds the indices
//   c[0],c[1],c[2],c[3] of the four vertices of the graph
//   whose (u,v) points are the furthest
//   SW, SE, NE, and NW in that order.
//   The indices can be used as corners for mapping to
//   the corners of the unit square when
//   reparametrising the boundary.
//   The array c should be allocated outside with length 4.
//   M.F. Aug. 97.
{
  for(int k=0; k<4; k++) c[k] = 0;
  for(int i=1; i< g_->getNumNodes(); i++)
  {
    if(g_->isBoundary(i))
    {
      if(g_->getU(i) + g_->getV(i) < g_->getU(c[0]) + g_->getV(c[0])) c[0] = i;
      if(g_->getU(i) - g_->getV(i) > g_->getU(c[1]) - g_->getV(c[1])) c[1] = i;
      if(g_->getU(i) + g_->getV(i) > g_->getU(c[2]) + g_->getV(c[2])) c[2] = i;
      if(g_->getU(i) - g_->getV(i) < g_->getU(c[3]) - g_->getV(c[3])) c[3] = i;
    }
  }
}



