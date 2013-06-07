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

