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

#include <fstream>
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/errormacros.h"

using namespace Go;
using namespace std;


void test_pt(vector<Point> p, int p_pos, int left_u, int left_v, vector<double> bas, const SplineSurface& surf)
{
  int n0 = surf.numCoefs_u();
  vector<double>::const_iterator coefs = surf.ctrl_begin();

  double q0 = 0.0, q1 = 0.0, q2 = 0.0;
  int dim = surf.dimension();
  int kdim = dim + (surf.rational() ? 1 : 0);

  int start_u = left_u - surf.order_u() + 1;
  int start_v = left_v - surf.order_v() + 1;

  int pos_j = kdim * (start_u + n0*start_v);
  int bas_pos = 0;

  if (kdim == dim)
    for (int j = 0; j < surf.order_v(); ++j, pos_j += n0*kdim)
      {
	int pos_i = pos_j;
	for (int i = 0; i < surf.order_u(); ++i, pos_i += kdim, ++bas_pos)
	  {
	    q0 += coefs[pos_i] * bas[bas_pos];
	    q1 += coefs[pos_i+1] * bas[bas_pos];
	    q2 += coefs[pos_i+2] * bas[bas_pos];
	  }
      }
  else
    for (int j = 0; j < surf.order_v(); ++j, pos_j += n0*kdim)
      {
	int pos_i = pos_j;
	for (int i = 0; i < surf.order_u(); ++i, pos_i += kdim, ++bas_pos)
	  {
	    q0 += coefs[pos_i] * bas[bas_pos] / coefs[pos_i+3];
	    q1 += coefs[pos_i+1] * bas[bas_pos] / coefs[pos_i+3];
	    q2 += coefs[pos_i+2] * bas[bas_pos] / coefs[pos_i+3];
	  }
      }

  Point q(q0, q1, q2);

  if (q.dist2(p[p_pos]) > 1e-5)
    {
      cout << "Position " << p_pos << " was point-evaluated to (" << p[p_pos] << ") but base-evaluated to (" << q << ")" << endl;
      exit(-1);
    }

}



int main(int argc, char* argv[] )
{
    ALWAYS_ERROR_IF(argc != 5, "Usage: " << argv[0]
		    << " surfaceinfile numpts_u numpts_v left(0)/right(1)" << endl);

    // Open input surface file
    ifstream is(argv[1]);
    ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

    // Read surface from file
    ObjectHeader head;
    SplineSurface surf;
    is >> head >> surf;

    vector<vector<double> > pars(2);
    int nmb_grid[2];

    for (int i = 0; i < 2; ++i)
      {
	double current_par = (i == 0) ? surf.startparam_u() : surf.startparam_u();
	double end_par = (i == 0) ? surf.endparam_u() : surf.endparam_u();
	nmb_grid[i] = atoi(argv[2+i]);
	double step = (end_par - current_par)/(double)(nmb_grid[i]-1);
	for (int j = 0; j < nmb_grid[i]; ++j, current_par += step) pars[i].push_back(current_par);
      }

    bool from_right = atoi(argv[4]) ? true : false;

    vector<BasisDerivsSf2> pts_with_derivs;
    surf.computeBasisGrid(pars[0], pars[1],
			  pts_with_derivs, from_right);

    vector<Point> p;
    p.resize(10);

    for (int i = 0; i < nmb_grid[0]; ++i)
      {
	double par_u = pars[0][i];
	for (int j = 0; j < nmb_grid[1]; ++j)
	  {
	    double par_v = pars[1][j];
	    surf.point(p, par_u, par_v, 2, from_right, from_right);

	    BasisDerivsSf2 bd = pts_with_derivs[i+nmb_grid[0]*j];

	    test_pt(p, 0, bd.left_idx[0], bd.left_idx[1], bd.basisValues, surf);
	    test_pt(p, 1, bd.left_idx[0], bd.left_idx[1], bd.basisDerivs_u, surf);
	    test_pt(p, 2, bd.left_idx[0], bd.left_idx[1], bd.basisDerivs_v, surf);
	    test_pt(p, 3, bd.left_idx[0], bd.left_idx[1], bd.basisDerivs_uu, surf);
	    test_pt(p, 4, bd.left_idx[0], bd.left_idx[1], bd.basisDerivs_uv, surf);
	    test_pt(p, 5, bd.left_idx[0], bd.left_idx[1], bd.basisDerivs_vv, surf);
	  }
      }
}
