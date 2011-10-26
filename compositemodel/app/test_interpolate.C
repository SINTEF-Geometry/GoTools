//===========================================================================
//
// File : test_interpolate.C
//
// Created: 2009-08
//
// Author: Vibeke Skytt
//
// Revision: 
//
// Description:
//
//===========================================================================

#ifdef __BORLANDC__
#include <vcl.h>
#endif


#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/igeslib/IGESconverter.h"

#include <iostream>
#include <fstream>
#include <stdlib.h> // For atof()
using std::shared_ptr;
using std::vector;
using std::ofstream;
using namespace Go;




int main( int argc, char* argv[] )
{
#ifdef __BORLANDC__
  using Go::Point;
#endif

  // Test number of input arguments
  if (argc != 6)
    {
      std::cout << "Input arguments : Input file(g2), degree, open/closed(1/0), give parameterization (0/1), Give type (0/1)" << std::endl;
      exit(-1);
    }


  double gap = 0.001;
    double neighbour = 0.01;
  double kink = 0.01;
  double bend = 0.1;
  double approx = 0.001;


  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");
  int degree = atoi(argv[2]);
  int open = atoi(argv[3]);
  int give_par = atoi(argv[4]);
  int give_type = atoi(argv[5]);

  IGESconverter conv;
  conv.readgo(file1);
  // Read curves and parameter values
  vector<shared_ptr<GeomObject> > gogeom = conv.getGoGeom();
  vector<shared_ptr<SplineCurve> > crvs;
  size_t ki;
  for (ki=0; ki<gogeom.size(); ++ki)
    {
      if (gogeom[ki]->instanceType() == Class_SplineCurve)
	{
	  shared_ptr<GeomObject> lg = gogeom[ki];
	  shared_ptr<SplineCurve> gocv =
	    std::dynamic_pointer_cast<SplineCurve, GeomObject>(lg);
	  crvs.push_back(gocv);	  
	}
    }

  vector<int> type(crvs.size());
  if (give_type)
    {
      int tp;
      for (int ki=0; ki < (int)crvs.size(); ++ki)
	{
	  std::cout << "Type of curve nr " << ki << ": " << std::endl;
	  std::cin >> tp;
	  type[ki] = tp;
	}
    }

  vector<double> param;
  if (give_par)
    {
      double par;
      for (int ki=0; ki < (int)crvs.size(); ++ki)
	{
	  if (give_par && type[ki] != 1)
	    continue;
	  std::cout << "Parameter value for curve nr " << ki << ": " << std::endl;
	  std::cin >> par;
	  param.push_back(par);
	}
      if (open <= 0)
	{
	  std::cout << "Parameter value for last curve:" << std::endl;
	  std::cin >> par;
	  param.push_back(par);
	}
    }

  CompositeModelFactory factory(approx, gap, neighbour, kink, bend);
  SurfaceModel *model1;
  if (give_type)
    {
      if (give_par)
	model1 = factory.interpolateCurves2(crvs, type, param, open, degree);
      else
	model1 = factory.interpolateCurves(crvs, type, param, open, degree);
    }
  else
    {
      if (give_par)
	model1 = factory.interpolateCurves2(crvs, param, open, degree);
      else
	model1 = factory.interpolateCurves(crvs, param, open, degree);
    }

  ofstream out_file("interpolate_surf.g2");
  if (model1)
    {
      int nmb = model1->nmbEntities();
      for (int kj=0; kj<nmb; kj++)
      {
	  shared_ptr<ParamSurface> surf = model1->getSurface(kj);

	  // Just to be sure
	  shared_ptr<SplineSurface> splinesf = 
	    std::dynamic_pointer_cast<SplineSurface, ParamSurface>(surf);

	  SplineSurface* surf2 = 
	    splinesf->subSurface(splinesf->startparam_u(),
				 splinesf->startparam_v(),
				 splinesf->endparam_u(),
				 splinesf->endparam_v());
	  surf2->writeStandardHeader(out_file);
	  surf2->write(out_file);

	  delete surf2;
      }
    }
  delete model1;
}

      
  
