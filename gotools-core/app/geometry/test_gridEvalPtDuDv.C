//===========================================================================
//
// File : test_gridEvalPtDuDv.C
//
// Created: Wed Dec  1 11:10:39 2010
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================



#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/Point.h"
#include <fstream>


using namespace Go;
using namespace std;


int main(int argc, char** argv)
{

  if (argc != 4)
    {
      cout << "Usage: " << argv[0] << " surfaceinfile num_u num_v" << endl;
      exit(-1);
    }

  ifstream filein(argv[1]);
  ALWAYS_ERROR_IF(filein.bad(), "Bad or no surface input filename");
  ObjectHeader head;
  filein >> head;
  if (head.classType() != SplineSurface::classType()) {
    THROW("Not a spline surface");
  }

  SplineSurface sf;
  filein >> sf;

  int num_u = atoi(argv[2]);
  int num_v = atoi(argv[3]);

  double start_u = sf.startparam_u();
  double end_u = sf.endparam_u();
  double start_v = sf.startparam_v();
  double end_v = sf.endparam_v();

  vector<double> params_u;
  vector<double> params_v;
  vector<double> points;
  vector<double> derivs_u;
  vector<double> derivs_v;

  double step, par;
  step = (end_u - start_u) / (double)(num_u);
  par = start_u + 0.5 * step;
  for (int i = 0; i < num_u; ++i, par += step)
    params_u.push_back(par);
  step = (end_v - start_v) / (double)(num_v);
  par = start_v + 0.5 * step;
  for (int i = 0; i < num_v; ++i, par += step)
    params_v.push_back(par);

  sf.gridEvaluator(params_u, params_v, points, derivs_u, derivs_v);

  int dim = sf.dimension();
  vector<double>::const_iterator pt_it = points.begin();
  vector<double>::const_iterator du_it = derivs_u.begin();
  vector<double>::const_iterator dv_it = derivs_v.begin();
  vector<Point> res;
  res.resize(3);
  for (int j = 0; j < num_v; ++j)
    for (int i = 0; i < num_u; ++i)
      {
	sf.point(res, params_u[i], params_v[j], 1);
	for (int k = 0; k < dim; ++k, ++pt_it, ++du_it, ++dv_it)
	  {
	    if (fabs((*pt_it)-res[0][k]) > 1.0e-9)
	      cout << "Point deviation at i=" << i << " j=" << j << " k=" << k << " girdeval=>" << (*pt_it) << " pteval=>" << res[0][k] << endl;
	    if (fabs((*du_it)-res[1][k]) > 1.0e-9)
	      cout << "Der_u deviation at i=" << i << " j=" << j << " k=" << k << " girdeval=>" << (*du_it) << " pteval=>" << res[1][k] << endl;
	    if (fabs((*dv_it)-res[2][k]) > 1.0e-9)
	      cout << "Der_v deviation at i=" << i << " j=" << j << " k=" << k << " girdeval=>" << (*dv_it) << " pteval=>" << res[2][k] << endl;
	  }
      }
}

