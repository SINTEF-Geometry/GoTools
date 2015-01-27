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

#include <iostream>
#include <fstream>
#include <time.h>
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRSplineUtils.h"



using namespace Go;
using namespace std;

void test_surface(const char* msg, LRSplineSurface* lrs, int samples_u = 100, int samples_v = 100, double tol = 1.0e-16)
{
  clock_t t_start = clock();
  cout << endl << "STARTING NEW TEST" << endl << msg << endl;
  cout << "Grid size = " << samples_u << " x " << samples_v << endl;
  cout << "Number of basis functions = " << lrs->numBasisFunctions() << endl;
  int space_dim = lrs->dimension();
  if (space_dim > 1)
    {
      cout << "Lowering from dimension " << space_dim << " to 1" << endl;
      LRSplineSurface::BSplineMap::iterator bspl_end = lrs->basisFunctionsEndNonconst();
      for (LRSplineSurface::BSplineMap::iterator bspl_it = lrs->basisFunctionsBeginNonconst(); bspl_it != bspl_end; ++bspl_it)
      {
	Point new_point(1);
	new_point[0] = bspl_it->second->coefTimesGamma()[space_dim - 1];
	bspl_it->second->coefTimesGamma() = new_point;
      }
      cout << "Dimension lowering completed" << endl;
    }

  double min_u = lrs->paramMin(XFIXED);
  double max_u = lrs->paramMax(XFIXED);
  double min_v = lrs->paramMin(YFIXED);
  double max_v = lrs->paramMax(YFIXED);

  double step_u = (max_u - min_u) / (double)(samples_u - 1);
  double step_v = (max_v - min_v) / (double)(samples_v - 1);
  vector<double> pars_u(samples_u);
  vector<double> pars_v(samples_v);
  vector<vector<double> > fnc_evals(samples_v);

  double u_val = min_u;
  for (int i = 0; i < samples_u; ++i, u_val += step_u)
    pars_u[i] = min(u_val, max_u);

  double v_val = min_v;
  for (int i = 0; i < samples_v; ++i, v_val += step_v)
    pars_v[i] = min(v_val, max_v);

  for (int j = 0; j < samples_v; ++j)
    {
      fnc_evals[j].resize(samples_u);
      for (int i = 0; i < samples_u; ++i)
	{
	  Point pt;
	  lrs->point(pt, pars_u[i], pars_v[j]);
	  fnc_evals[j][i] = pt[0];
	}
    }

  cout << "Completed grid evaluations before calling insertParameterFunctions" << endl;
  clock_t t_before_insert = clock();
  // LRSplineUtils::insertParameterFunctions(lrs);
  lrs->to3D();
  clock_t t_after_insert = clock();
  cout << "InsertParameterFunctions completed, now testing against grid evaluation" << endl;

  for (int j = 0; j < samples_v; ++j)
    for (int i = 0; i < samples_u; ++i)
      {
	Point pt_evaluated;
	lrs->point(pt_evaluated, pars_u[i], pars_v[j]);
	Point pt_expected(pars_u[i], pars_v[j], fnc_evals[j][i]);
	Point pt_diff = pt_evaluated - pt_expected;
	double diff2 = pt_diff.length2();
	if (diff2 >= tol)
	  {
	    cerr << endl << "*** TEST FAILED!!" << endl;
	    cerr << "Parameters = (" << pars_u[i] << "," << pars_v[j] << ")" << endl;
	    cerr << "Parameters grid index = (" << i << "," << j << ")" << endl;
	    cerr << "Evaluated point = " << pt_evaluated << endl;
	    cerr << "Expected point  = " << pt_expected << endl;
	    cerr << "Vector diff = " << pt_diff << endl;
	    cerr << "Length2 = " << diff2 << " accepted tolerance is " << tol << endl;
	    exit(1);
	  }
      }

  cout << "Test completed succesfully" << endl;
  clock_t t_end = clock();
  double sec_insert = (double)(t_after_insert - t_before_insert) / CLOCKS_PER_SEC;
  double sec_total = (double)(t_end - t_start) / CLOCKS_PER_SEC;
  double sec_test_wo_insert = sec_total - sec_insert;
  cout << "Time spent (secs)    calling insert : " << sec_insert << "  rest during test : " << sec_test_wo_insert << "  total : " << sec_total << endl;
}


int main(int argc, char *argv[])
{
  if (argc != 1)
  {
      std::cout << "Usage: argv[0]" << std::endl;
      return -1;
  }

  // Test 1
  // Bilinear bezier without refinement
  double kuv_ss1[] = {0.0, 0.0, 1.0, 1.0};
  double ctr_ss1[] = {1.3, -1.1, -3.0, 2.2};
  SplineSurface* ssurf1 = new SplineSurface(2, 2, 2, 2, &kuv_ss1[0], &kuv_ss1[0], &ctr_ss1[0], 1);
  LRSplineSurface* lrs1 = new LRSplineSurface(ssurf1, 0.01);
  test_surface("Simple bilinear case", lrs1);

  // Test 2
  // Bidegree (2,3) with more random knots, no refinement
  double ord_u2 = 3;
  double ord_v2 = 4;
  double ncoefs_u2 = 6;
  double ncoefs_v2 = 8;
  double ku_ss2[] = {1.3, 1.3, 1.3, 1.5, 2.1, 2.8, 3.3, 3.3, 3.3};
  double kv_ss2[] = {-7.1, -7.1, -7.1, -7.1, -5.0, -3.5, -3.3, -0.4, 1.1, 1.1, 1.1, 1.1};
  double ctr_ss2[] = {
    1.3, -1.1, -3.0, 2.4, 1.0, 0.7,
    -5.6, -7.0, -2.5, 9.5, 2.2, -2.2,
    -4.6, 0.3, 7.9, 3.5, -5.7, 1.9,
    4.0, 1.4, 3.2, 3.5, 3.7, -4.4,
    0.2, -5.4, -9.7, -5.9, 6.7, 4.7,
    -1.9, -5.9, -2.7, -2.1, 3.1, 3.2,
    -2.1, 3.2, 3.3, -7.2, 6.0, 0.8,
    -9.2, 1.2, -3.2, -1.0, 5.0, 7.8 };
  SplineSurface* ssurf2 = new SplineSurface(ncoefs_u2, ncoefs_v2, ord_u2, ord_v2, &ku_ss2[0], &kv_ss2[0], &ctr_ss2[0], 1);
  LRSplineSurface* lrs2 = new LRSplineSurface(ssurf2, 0.01);
  test_surface("Tensor case of bidegree (2,3)", lrs2);

  // Test 3
  // Simple insertion, regular knots, bidegree (2,3)
  double ord_u3 = 3;
  double ord_v3 = 4;
  double ncoefs_u3 = 7;
  double ncoefs_v3 = 9;
  double ku_ss3[] = {0.0, 0.0, 0.0, 1.0, 2.0, 4.0, 5.0, 6.0, 6.0, 6.0};
  double kv_ss3[] = {0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 6.0, 6.0, 6.0};
  double ctr_ss3[] = {
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  SplineSurface* ssurf3 = new SplineSurface(ncoefs_u3, ncoefs_v3, ord_u3, ord_v3, &ku_ss3[0], &kv_ss3[0], &ctr_ss3[0], 1);
  LRSplineSurface* lrs3 = new LRSplineSurface(ssurf3, 0.01);
  lrs3->refine(XFIXED, 3.0, 1.0, 5.0);
  test_surface("Simple refinement on regular case, bidegree (2,3)", lrs3, 7, 7);

  // For further tests, create reusable spline surface with at least one b-spline with inner knots only
  double ord_u4 = 3;
  double ord_v4 = 4;
  double ncoefs_u4 = 7;
  double ncoefs_v4 = 9;
  double ku_ss4[] = {1.3, 1.3, 1.3, 1.5, 2.1, 2.6, 2.8, 3.3, 3.3, 3.3};
  double kv_ss4[] = {-7.1, -7.1, -7.1, -7.1, -5.5, -5.0, -3.5, -3.3, -0.4, 1.1, 1.1, 1.1, 1.1};
  double ctr_ss4[] = {
    1.3, -1.1, -3.0, 2.4, 1.0, 0.7, -3.3,
    -4.1, 2.0, 2.1, 6.6, -1.2, -7.2, 9.6,
    -5.6, -7.0, -2.5, 9.5, 2.2, -2.2, -1.1,
    -4.6, 0.3, 7.9, 3.5, -5.7, 1.9, 0.8,
    4.0, 1.4, 3.2, 3.5, 3.7, -4.4, 8.0,
    0.2, -5.4, -9.7, -5.9, 6.7, 4.7, 8.2,
    -1.9, -5.9, -2.7, -2.1, 3.1, 3.2, -2.7,
    -2.1, 3.2, 3.3, -7.2, 6.0, 0.8, -5.7,
    -9.2, 1.2, -3.2, -1.0, 5.0, 7.8, 1.2};
  SplineSurface* ssurf4 = new SplineSurface(ncoefs_u4, ncoefs_v4, ord_u4, ord_v4, &ku_ss4[0], &kv_ss4[0], &ctr_ss4[0], 1);

  // Test 4
  // Bidegree (2,3) with more random knots, one refinement
  LRSplineSurface* lrs4 = new LRSplineSurface(ssurf4, 0.01);
  lrs4->refine(XFIXED, 2.3, -5.5, -0.4);
  test_surface("Simple refinement, bidegree (2,3)", lrs4);

  // Test 5
  // Bidegree (2,3) with more random knots, more refinements
  LRSplineSurface* lrs5 = new LRSplineSurface(ssurf4, 0.01);

  lrs5->refine(XFIXED, 2.3, -5.5, -0.4);
  lrs5->refine(YFIXED, -3.4, 1.5, 2.6);
  lrs5->refine(YFIXED, -4.6, 1.5, 2.6);
  lrs5->refine(XFIXED, 1.6, -7.1, -3.5);
  lrs5->refine(XFIXED, 2.37, -5.0, -0.4);
  test_surface("Several refinements, bidegree (2,3)", lrs5);

  cout << endl << "All tests completed succesfully" << endl;
}
