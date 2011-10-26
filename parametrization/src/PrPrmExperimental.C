//===========================================================================
//                                                                           
// File: PrPrmExperimental.C                                                 
//                                                                           
// Created: Wed Jul 16 16:48:04 2003                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: PrPrmExperimental.C,v 1.6 2007-03-02 16:35:52 jbt Exp $
//                                                                           
//===========================================================================


#include "GoTools/parametrization/PrPrmExperimental.h"
#include "GoTools/parametrization/PrParamUtil.h"

using namespace std;

// PRIVATE METHODS
namespace {

    // @@ SLOW HACK!
    double power(double x, int n)
    {
	double val = 1.0;
	while(n-- > 0) val *= x;
	return val;
    }

}

//-----------------------------------------------------------------------------
bool
PrPrmExperimental::makeWeights(int i)
//-----------------------------------------------------------------------------
//  Calculate weights for the
//  interior node i of the graph
//  It is assumed here that the indices of the neighbours of i
//  have already been stored in neighbours_.
//  This is an experimental parametrization.
{
    weights_.clear();
    int n = (int)neighbours_.size();

    int j;
    double tot_weight = 0.0, weight;
    double anglesum, testsum;
    anglesum = 0.0;
    testsum = 0.0;
    for (j=0; j<n; j++) {
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

	Vector3D ab = b - a;
	Vector3D bc = c - b;
	Vector3D da = a - d;
	Vector3D cd = d - c;

	weight = 4*area1*area1*area1*(ab*bc)/(ab.length2()*bc.length2())
	    + 4*area2*area2*area2*(da*cd)/(da.length2()*cd.length2());

	weights_.push_back(weight);
	tot_weight += weight;

// 	testsum += cos(anglesum)*(t1+t2);
// 	anglesum += angle2;
    }
//     cout << "Test sum = " << testsum << endl;

  // Scale the weights so that they sum to 1.

  double ratio = 1.0 / tot_weight;
  for(j=0; j<n; j++) weights_[j] *= ratio;

  return true;
}

// PUBLIC METHODS

//-----------------------------------------------------------------------------
PrPrmExperimental::PrPrmExperimental()
: PrParametrizeInt()
//-----------------------------------------------------------------------------
{
}

//-----------------------------------------------------------------------------
PrPrmExperimental::~PrPrmExperimental()
//-----------------------------------------------------------------------------
{
}

