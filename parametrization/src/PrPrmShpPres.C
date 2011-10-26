/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1997 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

/********************************************************************
 FILENAME    : PrPrmShpPres.C
 AUTHOR      : Michael Floater, SINTEF
 DATE        : Mar. 97
 DESCRIPTION : Implementation of methods in the class PrPrmShpPres.
 CHANGE LOG  :
*********************************************************************/

#include "GoTools/parametrization/PrPrmShpPres.h"
#include "GoTools/parametrization/PrParamUtil.h"
//#include "GoTools/utils/Values.h"
#ifdef _WIN32
#define M_PI 3.14159265358979
#endif

#ifdef __BORLANDC__
using std::cos;
using std::sin;
#endif

// PRIVATE METHODS


//-----------------------------------------------------------------------------
bool
PrPrmShpPres::makeWeights(int i)
//-----------------------------------------------------------------------------
//  Calculate shape-preserving weights for the
//  interior node i of the graph
//  It is assumed here that the indices of the neighbours of i
//  have already been stored in neighbours_.
//  This is the third parametrization in the article:
//  M. S. Floater, "Parametrization and smooth approximation of
//  surface triangulations", to appear in CAGD, 1997.
{
  weights_.clear();
  int n = (int)neighbours_.size();

  //  std::cout << "Testing" << std::endl;

  int j;
  for(j=0; j<n; j++) weights_.push_back(0.0);

  // Find local u,v's
  localParam(i);

  if (n == 2)
    {
      /* Added by VSK. Use chord length parametrization to get
	 the weights. */
      double length1 = g_->get3dNode(i).dist(g_->get3dNode(neighbours_[0]));
      double length2 = g_->get3dNode(i).dist(g_->get3dNode(neighbours_[1]));
      if (length1+length2 < 0.000000000001)
	  weights_[0] = weights_[1] = 0.5;
      else
	  {
	      weights_[0] = length1;
	      weights_[1] = length2;
 	  }
    }
  else
    {
      for(j=0; j<n; j++)
	{
	  /* Given the j-th neighbour of node i,
	     find the two neighbours by intersecting the
	     line through nodes i and j with all segments of the polygon
	     made by the neighbours. Take the two neighbours on
	     either side. Only one segment intersects this line. */

	  int k;
	  for(k=0; k<n; k++)
	    {
	      int kk = (k == n-1 ? 0 : k+1);
	      if(k == j || kk == j) continue;

	      double cross1 = det(u_[j],v_[j],u_[k],v_[k]);
	      double cross2 = det(u_[j],v_[j],u_[kk],v_[kk]);

	      if(cross1 * cross2 <= 0.0)
		{
		  double tau0,tau1,tau2;
		  baryCoords0(u_[j],v_[j],u_[k],v_[k],u_[kk],v_[kk],
			      tau0,tau1,tau2);
		  weights_[j]  += tau0;
		  weights_[k]  += tau1;
		  weights_[kk] += tau2;
		  break;
		}
	    }
	}
    }

  // Scale the weights so that they sum to 1.
  //  double ratio = 1.0 / (double)n;
  double sum = 0;
  for (j = 0; j < n; ++j)
      sum += weights_[j];
  double ratio = 1.0 / sum;
  for(j=0; j<n; j++) weights_[j] *= ratio;

  // Temporary code.
  /*
  ratio = 0.0;
  for(j=1; j<=n; j++)
  {
    weights_[j]  = 1.0 / weights_[j];
    ratio += weights_[j];
  }
  for(j=1; j<=n; j++) weights_[j] /= ratio;

  double tmp;

  tmp = weights_(1);
  weights_(1) = weights_(3);
  weights_(3) = tmp;

  tmp = weights_(2);
  weights_(2) = weights_(4);
  weights_(4) = tmp;
  */
  // End of temporary code.

  return true;
}

//-----------------------------------------------------------------------------
bool
PrPrmShpPres::localParam(int i)
//-----------------------------------------------------------------------------
//  Make a local parametrization for the interior
//  node i of the given graph.
//  This is based on a discretization of the geodesic polar map.
//  The distances of the neighbours from i are preserved
//  and the ratios of any two interior angles are also preserved
//  It is assumed here that the indices of the neighbours of i
//  have already been stored in neighbours_.
//  M.F. Mar. 97
{
  alpha_.clear();
  len_.clear();
  u_.clear();
  v_.clear();
  int n = (int)neighbours_.size();

  double alpha_sum = 0.0, alpha;

  int j;
  for (j=0; j<n; j++)
  {
    int j1 = neighbours_[j];
    Vector3D v1 = g_->get3dNode(j1) - g_->get3dNode(i);

    // get previous node
    int jprev = (j == 0 ? n-1 : j-1);
    int j2 = neighbours_[jprev];
    Vector3D v2 = g_->get3dNode(j2) - g_->get3dNode(i);

    len_.push_back(v1.length());
    //len_[j] = sqrt(v1.length()); // alternative
    //len_[j] = 1.0; // alternative
    alpha = v1.angle(v2);
    alpha_.push_back(alpha);
    //alpha_[j] = sqrt(v1.angle(v2)); // alternative
    //alpha_[j] = 1.0; // alternative
    alpha_sum += alpha;
  }

    double factor = 2.0 * M_PI / alpha_sum;

    for(j=0; j<n; j++) alpha_[j] *= factor;

    alpha_[0] = 0.0;
    for(j=1; j<n; j++) alpha_[j] += alpha_[j-1];

    for(j=0; j<n; j++)
    {
    //v_[j].init(len_[j] * cos(alpha_[j]), len_[j] * sin(alpha_[j]));
    //v_.push_back(CgVector2d(len_[j] * cos(alpha_[j]), len_[j] * sin(alpha_[j])));
      u_.push_back(len_[j] * cos(alpha_[j]));
      v_.push_back(len_[j] * sin(alpha_[j]));
    }

  return true;
}

// PUBLIC METHODS

//-----------------------------------------------------------------------------
PrPrmShpPres::PrPrmShpPres()
: PrParametrizeInt()
//-----------------------------------------------------------------------------
{
}

//-----------------------------------------------------------------------------
PrPrmShpPres::~PrPrmShpPres()
//-----------------------------------------------------------------------------
{
}

