//===========================================================================
//
// File : testCompositeModel.C
//
// Created: Thu Feb 21 09:37:14 2008
//
// Author: Kjell Fredrik Pettersen
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
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"

#include <iostream>
#include <fstream>
#include <stdlib.h> // For atof()
using std::shared_ptr;
using std::cout;
using std::endl;
using namespace Go;




void dump(int modelNumber, SurfaceModel* sm)
{
  cout << "\n\n";
  cout << "Dumping model " << modelNumber << endl;

  printf("Model = %p\n",sm);

  cout << "#Faces = " << sm->nmbEntities() << endl;

  int top = sm->nmbEntities();

  for (int i=0; i<top; ++i)
    {
      shared_ptr<ftSurface> f;
      shared_ptr<ParamSurface> f2;
      f = sm -> getFace(i);
      f2 = sm -> getSurface(i);
      printf("Face %i = %p, param = %p\n", i, f.get(), f2.get());
    }
}


void dump2(int modelNumber, CompositeCurve* cv)
{
  cout << "\n\n";
  cout << "Dumping model " << modelNumber << endl;

  printf("Model = %p\n",cv);

  cout << "#Curves = " << cv->nmbEntities() << endl;

  int top = cv->nmbEntities();

  for (int i=0; i<top; ++i)
    {
      shared_ptr<ParamCurve> f2;
      f2 = cv -> getCurve(i);
      printf("Curve %i = %p\n", i, f2.get());
    }

  double start, end;
  cv->parameterRange(start, end);
  std::cout << "Parameter range: " << start << ", " << end << std::endl;
  
  double par = 0.2*start + 0.8*end;
  Point pt;
  cv->evaluateCurve(par, pt);
  std::cout << "Par: " << par << " Pos: " << pt[0] << " " << pt[1];
  std::cout << " " << pt[2] << std::endl;
}










int main( int argc, char* argv[] )
{
  // Test number of input arguments
  if (argc != 2)
    {
      std::cout << "Input arguments : Input file on IGES format" << std::endl;
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


  CompositeModelFactory factory(approx, gap, neighbour, kink, bend);
  CompositeModel* cm = factory.createFromIges(file1);

  CompositeModel *cm2 = cm->clone();
  SurfaceModel* sf1 = dynamic_cast<SurfaceModel*>(cm);
  SurfaceModel* sf2 = dynamic_cast<SurfaceModel*>(cm2);

  CompositeCurve *cv1 = dynamic_cast<CompositeCurve*>(cm);
  CompositeCurve *cv2 = dynamic_cast<CompositeCurve*>(cm2);
  if (sf1)
    {
      dump(1,sf1);
      dump(2,sf2);
    }
  else if (cv1)
    {
      dump2(1,cv1);
      dump2(2,cv2);
    }

}

















