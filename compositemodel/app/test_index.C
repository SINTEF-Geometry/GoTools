//===========================================================================
//
// File : test_index.C
//
// Created: Thu Apr  3 13:45:07 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: test_index.C,v 1.2 2009-05-13 07:29:53 vsk Exp $
//
// Description:
//
//===========================================================================


#ifdef __BORLANDC__
#include <vcl.h>
#endif



#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"

#include <iostream>
#include <fstream>
#include <stdlib.h> // For atof()
using std::shared_ptr;
using namespace std;
using namespace Go;




int main( int argc, char* argv[] )
{
  // Test number of input arguments
  if ((argc) != 2)
    {
      std::cout << "Input arguments : Input file on IGES format," << std::endl;
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
  SurfaceModel* sm = (SurfaceModel*) factory.createFromIges(file1);

  int nmb_faces = (int)sm->nmbEntities();
  int ers = 0;

  for (int i = 0; i < nmb_faces; ++i)
    {
      if (sm -> getIndex(sm -> getFace(i)) != i)
	{
	  ++ers;
	  std::cout << "Error for getFace(" << i << ")" << std::endl;
	}
      if (sm -> getIndex(sm -> getFace(i).get()) != i)
	{
	  ++ers;
	  std::cout << "Error for getFace(" << i << ").get()" << std::endl;
	}
    }
  
  std::cout << "Testing complete, found " << ers << " errors" << std::endl;

}
