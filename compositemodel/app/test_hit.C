//===========================================================================
//                                                                           
// File: test_intersect.C                                                     
//                                                                           
// Created: Tue Mar 21 16:44:59 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: test_lineIntersect.C,v 1.3 2009-05-13 07:29:54 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#ifdef __BORLANDC__
#include <vcl.h>
#endif

#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/ftPoint.h"
#include <fstream>
#include <stdlib.h> // For atof()


const double EPS_GAP = 1e-7;
const double EPS_NEIGHBOUR = 1e-4;
const double EPS_KINK = 1e-2;
const double EPS_BEND = 4e-2;

const double EPS_APPROX = 1e-4;


using namespace std;
using namespace Go;


int main( int argc, char* argv[] )
{
#ifdef __BORLANDC__
  using Go::Point;
#endif

  if (argc != 8) {
    std::cout << "Input parameters : Input file on iges format, p0_x, p0_y, .., p1_x" << std::endl;
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


  ifstream file(argv[1]);

  CompositeModelFactory factory(approx, gap, neighbour, kink, 10.0*kink);

  CompositeModel *model = factory.createFromIges(file);

  SurfaceModel *sfmodel = dynamic_cast<SurfaceModel*>(model);
  CompositeCurve *cvmodel = dynamic_cast<CompositeCurve*>(model);

  Point p0(atof(argv[2]), atof(argv[3]), atof(argv[4]));
  Point p1(atof(argv[5]), atof(argv[6]), atof(argv[7]));

  Point dir = p1-p0;

  std::cout << "Starting search" << std::endl;

  if (sfmodel)
    {
      ftPoint result(0.0, 0.0, 0.0);
      bool hit = sfmodel->hit(p0, dir, result);
      std::cout << "Hit found " << hit << std::endl;

      std::ofstream out_hit("hit.g2");
      out_hit << "400 1 0 4 255 0 0 255" << std::endl;
      out_hit << "1" << std::endl;
      out_hit << p0 << std::endl;
      out_hit << "410 1 0 4 255 0 0 255" << std::endl;
      out_hit << "1" << std::endl;
      out_hit << p0 << " " << p1 << std::endl;
      if (hit)
	{
	  out_hit << "400 1 0 4 0 255 0 255" << std::endl;
	  out_hit << 1 << std::endl;
	  out_hit << result.position() << std::endl;
	}
    }
  else if (cvmodel)
    {
      PointOnCurve result;
      bool hit = cvmodel->hit(p0, dir, result);
      std::cout << "Hit found " << hit << std::endl;

      std::ofstream out_hit("hit.g2");
      out_hit << "400 1 0 4 255 0 0 255" << std::endl;
      out_hit << "1" << std::endl;
      out_hit << p0 << std::endl;
      out_hit << "410 1 0 4 255 0 0 255" << std::endl;
      out_hit << "1" << std::endl;
      out_hit << p0 << " " << p1 << std::endl;
      if (hit)
	{
	  out_hit << "400 1 0 4 0 255 0 255" << std::endl;
	  out_hit << 1 << std::endl;
	  out_hit << result.getPos() << std::endl;
	}
    }
}
