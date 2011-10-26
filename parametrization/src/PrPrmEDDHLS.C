/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1997 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

/********************************************************************
 FILENAME    : PrPrmEDDHLS.C
 AUTHOR      : Michael Floater, SINTEF
 DATE        : Mar. 97
 DESCRIPTION : Implementation of methods in the class PrPrmEDDHLS.
 CHANGE LOG  :
*********************************************************************/

#include "GoTools/parametrization/PrPrmEDDHLS.h"
#include "GoTools/parametrization/PrParamUtil.h"

// PRIVATE METHODS

//-----------------------------------------------------------------------------
bool
PrPrmEDDHLS::makeWeights(int i)
//-----------------------------------------------------------------------------
//   Calculate Eck, DeRose, Duchamp, Hoppe, Lounsbery, Stuetzle weights
//   for the interior node i of the graph.
//   The formula for the weight of edge (i,j) is
//   w_jk = (L_{i,k_1}^2 + L_{j,k_1}^2 - L_{i,j}^2) / Area_{i,j,k_1}
//             +
//          (L_{i,k_2}^2 + L_{j,k_2}^2 - L_{i,j}^2) / Area_{i,j,k_2}.
//
//   These weights are not necessarily positive!!!
//   So it's better to use the class PrPrmShpPres.
//   M.F. Mar. 97.
{
  weights_.clear();
  int n = (int)neighbours_.size();

  int j;
  double tot_weight = 0.0, weight;
  for(j=0; j<n; j++)
  {
    int jj = neighbours_[j];
    int jprev = (j == 0 ? n-1 : j-1);
    int j1 = neighbours_[jprev];
    Vector3D a = g_->get3dNode(i);
    Vector3D b = g_->get3dNode(j1);
    Vector3D c = g_->get3dNode(jj);
    double area1 = area(a,b,c);

    int jnext = (j == n-1 ? 0 : j+1);
    int j2 = neighbours_[jnext];
    Vector3D d = g_->get3dNode(j2);
    double area2 = area(a,c,d);

    double l1 = a.dist2(c);
    double l2 = a.dist2(b);
    double l3 = c.dist2(b);
    double l4 = a.dist2(d);
    double l5 = c.dist2(d);

    weight = (l2+l3-l1) / area1 + (l4+l5-l1) / area2;
    weights_.push_back(weight);
#ifdef PRDEBUG
    if(weight < 0.0)
    {
      cout << "negative weight=" << weight;
      cout << " for edge(" << i << "," << jj << ")" << "\n";
    }
#endif
    tot_weight += weight;
  }

#ifdef PRDEBUG
    if(tot_weight < 0.0)
    {
      cout << "negative sum of weight=" << tot_weight;
      cout << " for node " << i << "\n";
    }
#endif

  // Scale the weights so that they sum to 1.

  double ratio = 1.0 / tot_weight;
  for(j=0; j<n; j++)
  {
    weights_[j] *= ratio;
#ifdef PRDEBUG
    if(weights_[j] < 0.0)
    {
      cout << "negative weight=" << weights_[j];
      cout << " for edge(" << i << "," << neighbours_[j] << ")" << "\n";
    }
#endif
  }

  return true;
}

// PUBLIC METHODS

//-----------------------------------------------------------------------------
PrPrmEDDHLS::PrPrmEDDHLS()
: PrParametrizeInt()
//-----------------------------------------------------------------------------
{
}

//-----------------------------------------------------------------------------
PrPrmEDDHLS::~PrPrmEDDHLS()
//-----------------------------------------------------------------------------
{
}

