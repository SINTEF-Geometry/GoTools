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
#include "GoTools/igeslib/IGESconverter.h"
#include "GoTools/compositemodel/ftMessage.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/IntResultsSfModel.h"
#include "GoTools/geometry/LineCloud.h"
#include "GoTools/tesselator/LineStrip.h"
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

  if (argc != 6 && argc != 8) {
    std::cout << "Input parameters : Input file on IGES format, p0_x, p0_y, .., p1_x" << std::endl;
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

  Point p0, p1;
  if (argc == 8)
    {
      p0 = Point(atof(argv[2]), atof(argv[3]), atof(argv[4]));
      p1 = Point(atof(argv[5]), atof(argv[6]), atof(argv[7]));
    }
  else
    {
      p0 = Point(atof(argv[2]), atof(argv[3]));
      p1 = Point(atof(argv[4]), atof(argv[5]));
    }


  ftLine line(p1-p0, p0);

  std::ofstream out("line_ints.g2");
  double density = 1.0;
  std::shared_ptr<IntResultsModel> results = model->intersect(line);

  std::vector<std::shared_ptr<LineStrip> > line_seg;
  PointCloud3D points;
  results->tesselate(density, line_seg, points);

  for (size_t ki=0; ki<line_seg.size(); ++ki)
    {
      LineCloud lines;
      vector<double> tmp_lines;
      int nmb_vert = line_seg[ki]->numVertices();
      double *vertices = line_seg[ki]->vertexArray();
      for (int kj=0; kj<nmb_vert-1; ++kj)
	tmp_lines.insert(tmp_lines.end(), vertices+kj*3, vertices+(kj+2)*3);

      lines.setCloud(&tmp_lines[0], nmb_vert-1);
      lines.writeStandardHeader(out);
      lines.write(out);
    }
  if (points.numPoints() > 0)
    {
      points.writeStandardHeader(out);
      points.write(out);
    }
      


  int stop_break;
  stop_break = 1;
}
