//===========================================================================
//                                                                           
// File: PrPrmSymMeanValue.C                                                 
//                                                                           
// Created: Tue Jun 24 15:32:50 2003                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: PrPrmSymMeanValue.C,v 1.3 2007-03-02 16:35:52 jbt Exp $
//                                                                           
//===========================================================================


#include "GoTools/parametrization/PrPrmSymMeanValue.h"
#include "GoTools/parametrization/PrParamUtil.h"

// PRIVATE METHODS

//-----------------------------------------------------------------------------
bool
PrPrmSymMeanValue::makeWeights(int i)
//-----------------------------------------------------------------------------
//  Calculate mean value weights for the
//  interior node i of the graph
//  It is assumed here that the indices of the neighbours of i
//  have already been stored in neighbours_.
//  This is an attempt to improve the mean value method of
//  PrPrmMeanValue by making it symmetric.
{
    THROW("This method is bad, and suitable for experimentation only.");

  weights_.clear();
  int n = (int)neighbours_.size();

  //  std::cout << "Testing" << std::endl;


  int j;
  double tot_weight = 0.0, weight;
  for (j=0; j<n; j++)
  {
    int jj = neighbours_[j];
    int jprev = (j == 0 ? n-1 : j-1);
    int j1 = neighbours_[jprev];
    int jnext = (j == n-1 ? 0 : j+1);
    int j2 = neighbours_[jnext];
    Vector3D d = g_->get3dNode(j2);
    Vector3D a = g_->get3dNode(i);
    Vector3D b = g_->get3dNode(j1);
    Vector3D c = g_->get3dNode(jj);
    double t1 = tanThetaOverTwo(a,b,c);
    //std::cout << "t1 = " << t1 << endl;

    double t2 = tanThetaOverTwo(a,c,d);
    //std::cout << "t2 = " << t2 << endl;

    double t3 = tanThetaOverTwo(c,d,a);
    //std::cout << "t3 = " << t2 << endl;

    double t4 = tanThetaOverTwo(c,a,b);
    //std::cout << "t4 = " << t2 << endl;

    double len = a.dist(c);
    //std::cout << "len = " << len << endl;

    weight = (t1 + t2 + t3 + t4) / len;
    //std::cout << "weight = " << weight << endl;
    weights_.push_back(weight);
    tot_weight += weight;
  }

  // Scale the weights so that they sum to 1.

  double ratio = 1.0 / tot_weight;
  for(j=0; j<n; j++) weights_[j] *= ratio;

  return true;
}

//-----------------------------------------------------------------------------
double
PrPrmSymMeanValue::tanThetaOverTwo(Vector3D& a, Vector3D& b, Vector3D& c)
//-----------------------------------------------------------------------------
// Return tangent of half the angle between vectors b-a and c-a
// without using trig functions.
// Use fact that tan(alpha/2) = (1-cos(alpha)) / sin(alpha).
// and use scalar and dot products to get cos(alpha) and sin(alpha).
//  M.F. Apr. 2002.
{
  Vector3D bb = b - a;
  Vector3D cc = c - a;

  double cp = (bb.cross(cc)).length();
       // length of cross product of bb and cc
  double dp = bb * cc; // scalar product of bb and cc
  double bc = bb.length() * cc.length();
  return (bc - dp) / cp;

  return true;
}

// PUBLIC METHODS

//-----------------------------------------------------------------------------
PrPrmSymMeanValue::PrPrmSymMeanValue()
: PrParametrizeInt()
//-----------------------------------------------------------------------------
{
}

//-----------------------------------------------------------------------------
PrPrmSymMeanValue::~PrPrmSymMeanValue()
//-----------------------------------------------------------------------------
{
}

