//===========================================================================
//
// File : test_combasGridDeriv2.C
//
// Created: Wed Jul  8 13:32:14 2009
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: test_combasGridDeriv2.C,v 1.1 2009/07/09 12:58:44 kfp Exp $
//
// Description:
//
//===========================================================================


#include <fstream>
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/errormacros.h"

using namespace Go;
using namespace std;


void test_pt(vector<Point> p, int p_pos, int left_u, int left_v, int left_w, vector<double> bas, const SplineVolume& vol)
{
  int n0 = vol.numCoefs(0);
  int n1 = vol.numCoefs(1);
  vector<double>::const_iterator coefs = vol.ctrl_begin();

  double q0 = 0.0, q1 = 0.0, q2 = 0.0;
  int dim = vol.dimension();
  int kdim = dim + (vol.rational() ? 1 : 0);

  int start_u = left_u - vol.order(0) + 1;
  int start_v = left_v - vol.order(1) + 1;
  int start_w = left_w - vol.order(2) + 1;

  int pos_k = kdim * (start_u + n0*(start_v + n1*start_w));
  int bas_pos = 0;

  if (kdim == dim)
    {
      for (int k = 0; k < vol.order(2); ++k, pos_k += n0*n1*kdim)
	{
	  int pos_j = pos_k;
	  for (int j = 0; j < vol.order(1); ++j, pos_j += n0*kdim)
	    {
	      int pos_i = pos_j;
	      for (int i = 0; i < vol.order(0); ++i, pos_i += kdim, ++bas_pos)
		{
		  q0 += coefs[pos_i] * bas[bas_pos];
		  q1 += coefs[pos_i+1] * bas[bas_pos];
		  q2 += coefs[pos_i+2] * bas[bas_pos];
		}
	    }
	}
    }
  else
    {
      for (int k = 0; k < vol.order(2); ++k, pos_k += n0*n1*kdim)
	{
	  int pos_j = pos_k;
	  for (int j = 0; j < vol.order(1); ++j, pos_j += n0*kdim)
	    {
	      int pos_i = pos_j;
	      for (int i = 0; i < vol.order(0); ++i, pos_i += kdim, ++bas_pos)
		{
		  q0 += coefs[pos_i] * bas[bas_pos] / coefs[pos_i+3];
		  q1 += coefs[pos_i+1] * bas[bas_pos] / coefs[pos_i+3];
		  q2 += coefs[pos_i+2] * bas[bas_pos] / coefs[pos_i+3];
		}
	    }
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
    ALWAYS_ERROR_IF(argc != 6, "Usage: " << argv[0]
		    << " volumeinfile numpts_u numpts_v numpts_w left(0)/right(1)" << endl);

    // Open input volume file
    ifstream is(argv[1]);
    ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

    // Read volume from file
    ObjectHeader head;
    SplineVolume vol;
    is >> head >> vol;

    vector<vector<double> > pars(3);
    int nmb_grid[3];

    for (int i = 0; i < 3; ++i)
      {
	double current_par = vol.startparam(i);
	nmb_grid[i] = atoi(argv[2+i]);
	double step = (vol.endparam(i) - current_par)/(double)(nmb_grid[i]-1);
	for (int j = 0; j < nmb_grid[i]; ++j, current_par += step) pars[i].push_back(current_par);
      }

   bool from_right = atoi(argv[5]) ? true : false;

    vector<BasisDerivs2> pts_with_derivs;
    vol.computeBasisGrid(pars[0], pars[1], pars[2],
			 pts_with_derivs, from_right);

    vector<Point> p;
    p.resize(10);

    for (int i = 0; i < nmb_grid[0]; ++i)
      {
	double par_u = pars[0][i];
	for (int j = 0; j < nmb_grid[1]; ++j)
	  {
	    double par_v = pars[1][j];
	    for (int k = 0; k < nmb_grid[2]; ++k)
	      {
		double par_w = pars[2][k];
		vol.point(p, par_u, par_v, par_w, 2, 
			  from_right, from_right, from_right);

		BasisDerivs2 bd = pts_with_derivs[i+nmb_grid[0]*(j+k*nmb_grid[1])];

		test_pt(p, 0, bd.left_idx[0], bd.left_idx[1], bd.left_idx[2], bd.basisValues, vol);
		test_pt(p, 1, bd.left_idx[0], bd.left_idx[1], bd.left_idx[2], bd.basisDerivs_u, vol);
		test_pt(p, 2, bd.left_idx[0], bd.left_idx[1], bd.left_idx[2], bd.basisDerivs_v, vol);
		test_pt(p, 3, bd.left_idx[0], bd.left_idx[1], bd.left_idx[2], bd.basisDerivs_w, vol);
		test_pt(p, 4, bd.left_idx[0], bd.left_idx[1], bd.left_idx[2], bd.basisDerivs_uu, vol);
		test_pt(p, 5, bd.left_idx[0], bd.left_idx[1], bd.left_idx[2], bd.basisDerivs_uv, vol);
		test_pt(p, 6, bd.left_idx[0], bd.left_idx[1], bd.left_idx[2], bd.basisDerivs_uw, vol);
		test_pt(p, 7, bd.left_idx[0], bd.left_idx[1], bd.left_idx[2], bd.basisDerivs_vv, vol);
		test_pt(p, 8, bd.left_idx[0], bd.left_idx[1], bd.left_idx[2], bd.basisDerivs_vw, vol);
		test_pt(p, 9, bd.left_idx[0], bd.left_idx[1], bd.left_idx[2], bd.basisDerivs_ww, vol);
	      }
	  }
      }
}
