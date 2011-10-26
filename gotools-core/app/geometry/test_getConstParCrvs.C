//===========================================================================
//
// File : test_getConstParCrvs.C
//
// Created: Fri Jul 10 12:24:44 2009
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: test_getConstParCrvs.C,v 1.1 2009/07/10 10:38:16 kfp Exp $
//
// Description:
//
//===========================================================================



#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <memory>
#include <fstream>

using namespace Go;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::cout;
using std::endl;
using std::shared_ptr;;

int main(int argc, char** argv)
{

  if (argc != 5)
    {
      cout << "Usage: " << argv[0] << " surfaceinfile curvesoutfile nmb_crvs_u nmb_crvs_v" << endl;
      exit(-1);
    }

  ifstream filein(argv[1]);
  ALWAYS_ERROR_IF(filein.bad(), "Bad or no surface input filename");
  ObjectHeader head;
  filein >> head;
  if (head.classType() != SplineSurface::classType()) {
    THROW("Not a spline surface");
  }
  SplineSurface ss;
  filein >> ss;

  ofstream fileout(argv[2]);
  ALWAYS_ERROR_IF(fileout.bad(), "Bad curves output filename");

  int nmb_u = atoi(argv[3]);
  int nmb_v = atoi(argv[4]);
  vector<double> par_u, par_v;

  double step;
  double par;

  par = ss.startparam_u();
  if (nmb_u == 1)
    step = 0.0;
  else
    step = (ss.endparam_u()-par) / (double)(nmb_u-1);
  for (int i = 0; i < nmb_u; ++i, par += step)
    par_u.push_back(par);

  par = ss.startparam_v();
  if (nmb_v == 1)
    step = 0.0;
  else
    step = (ss.endparam_v()-par) / (double)(nmb_v-1);
  for (int i = 0; i < nmb_v; ++i, par += step)
    par_v.push_back(par);

  vector<shared_ptr<SplineCurve> > curves_u, curves_v;

  ss.getConstParamCurves(par_u, par_v, curves_u, curves_v);

  for (int i = 0; i < (int)curves_u.size(); ++i)
    {
      curves_u[i]->writeStandardHeader(fileout);
      curves_u[i]->write(fileout);
    }
  for (int i = 0; i < (int)curves_v.size(); ++i)
    {
      curves_v[i]->writeStandardHeader(fileout);
      curves_v[i]->write(fileout);
    }

  double eps = 1.0e-8;
  shared_ptr<BoundedSurface> bd_surf = 
    shared_ptr<BoundedSurface>(BoundedUtils::convertToBoundedSurface(ss, eps));

  size_t ki;
  for (ki=0; ki<par_u.size(); ++ki)
    {
      vector<shared_ptr<ParamCurve> > cvs = 
	bd_surf->constParamCurves(par_u[ki], false);
      for (size_t kj=0; kj<cvs.size(); ++kj)
	{
	  shared_ptr<SplineCurve> cc = 
	    shared_ptr<SplineCurve>(cvs[kj]->geometryCurve());
	  cc->writeStandardHeader(fileout);
	  cc->write(fileout);
	}
    }

  for (ki=0; ki<par_v.size(); ++ki)
    {
      vector<shared_ptr<ParamCurve> > cvs = 
	bd_surf->constParamCurves(par_v[ki], true);
      for (size_t kj=0; kj<cvs.size(); ++kj)
	{
	  shared_ptr<SplineCurve> cc = 
	    shared_ptr<SplineCurve>(cvs[kj]->geometryCurve());
	  cc->writeStandardHeader(fileout);
	  cc->write(fileout);
	}
    }
								      
}
