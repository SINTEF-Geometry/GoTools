//===========================================================================
//
// File : testClosestPoint.C
//
// Created: Thu Feb 28 16:01:48 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: testClosestPoint.C,v 1.4 2009-05-13 07:29:53 vsk Exp $
//
// Description:
//
//===========================================================================

#ifdef __BORLANDC__
#include <vcl.h>
#endif

#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/ftMessage.h"
#include <istream>
#include <fstream>
#include <stdlib.h> // For atof()


using std::shared_ptr;
using std::vector;
using std::istream;


using namespace Go;



int main( int argc, char* argv[] )
{
  // Test number of input arguments
  if (argc < 5 || (argc % 3) != 2)
    {
      std::cout << "Input arguments : Input file on IGES format, ";
      std::cout << "p1X, p1Y, p1Z, ... , pnX, pxY, pnZ" << std::endl;

      exit(-1);
    }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");


  double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  // Removed from corresponding use of fairingToolbox
  double approx = 0.001;

  vector<double> coord;
  int nmb_points = (argc - 2) / 3;
  for (int i = 0; i < nmb_points * 3; ++i) coord.push_back (atof(argv[i+2]));

  vector<int> rgb;
  rgb.push_back(255);  rgb.push_back(0);  rgb.push_back(0);
  rgb.push_back(0);  rgb.push_back(255);  rgb.push_back(0);
  rgb.push_back(255);  rgb.push_back(255);  rgb.push_back(0);
  rgb.push_back(0);  rgb.push_back(255);  rgb.push_back(255);
  rgb.push_back(255);  rgb.push_back(128);  rgb.push_back(0);
  rgb.push_back(255);  rgb.push_back(64);  rgb.push_back(128);
  rgb.push_back(255);  rgb.push_back(0);  rgb.push_back(0);
  rgb.push_back(128);  rgb.push_back(255);  rgb.push_back(64);

  
  // Create the Surface model
  CompositeModelFactory factory(approx, gap, neighbour, kink, 10.0*kink);

  CompositeModel *model = factory.createFromIges(file1);


  // Left out from fairingToolbox
  // Make clean degenerate surfaces and remove surfaces degenerated to a line
  // int nmb_deg=0, nmb_bd=0, nmb_removed=0;
  // status = tool->ensureCleanDegeneracy(nmb_deg, nmb_bd, nmb_removed);
  // std::cout << "Clean degeneracy. Status message : " << status.getMessage() << std::endl;
  // std::cout << "Number of degenerate triangular or banana surfaces: " << nmb_deg << std::endl;
  // std::cout << "Number of modified degenerate boundaries: " << nmb_bd << std::endl;
  // std::cout << "Number of removed surfaces: " << nmb_removed << std::endl;
  
  std::ofstream baseStr("data/basePoint.g2");

  for (int i = 0; i < nmb_points; ++i)
    {

      Point p(coord[i*3], coord[i*3+1], coord[i*3+2]);
      Point clp;
      int idx;
      double clo_p[2];
      double dist;

      model->closestPoint(p, clp, idx, clo_p, dist);
      if (nmb_points == 1)
	{
	  //ftSurface* closestSurface;
	  //Point norm;
	  //Point normCross;

	  //closestSurface = model->getSurface2(idx);
	  //norm = closestSurface -> normal(clo_p[0], clo_p[1]);
	  //normCross = norm % (p - clp);

	  std::cout << "Closest point is " << clp << std::endl;
	  std::cout << "Surface index is " << idx << std::endl;
	  std::cout << "Parameters on surface are (" << clo_p[0] << ", " << clo_p[1] << ")" << std::endl;
	  std::cout << "Distance is " << dist << std::endl;
	  //std::cout << "Norm cross length is " << normCross.length() << std::endl;
	}

      baseStr << "400 1 0 4 \n";
      int colPos = (i & 7) * 3;
      baseStr << rgb[colPos] << " " << rgb[colPos + 1] << " " << rgb[colPos + 2] << " 255\n";
      baseStr << "1\n";
      baseStr << p << "\n";
      baseStr << "400 1 0 4 \n";
      baseStr << (rgb[colPos]>>1) << " " << (rgb[colPos + 1]>>1) << " " << (rgb[colPos + 2]>>1) << " 255\n";
      baseStr << "1\n";
      baseStr << clp << "\n";
  
    }
}
