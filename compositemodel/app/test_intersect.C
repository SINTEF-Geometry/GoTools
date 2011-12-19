//===========================================================================
//                                                                           
// File: test_intersect.C                                                     
//                                                                           
// Created: Tue Mar 21 16:44:59 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: test_intersect.C,v 1.4 2009-05-13 07:29:53 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifdef __BORLANDC__
#include <vcl.h>
#endif

#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/igeslib/IGESconverter.h"
#include "GoTools/compositemodel/ftMessage.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/igeslib/ftTangPriority.h"
#include "GoTools/tesselator/GeneralMesh.h"
#include "GoTools/compositemodel/IntResultsModel.h"
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


using std::ifstream;
using std::ofstream;
using namespace Go;


int main( int argc, char* argv[] )
{
#ifdef __BORLANDC__
  using Go::Point;
#endif

  if (argc != 11) {
    std::cout << "Input parameters : Input file on IGES format, ";
    std::cout << "x-value of first point on plane, y-value f first point, ... , z-value of third point" << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");


  double gap = EPS_GAP; //0.0001;
  double neighbour = EPS_NEIGHBOUR; //0.001;
  double kink = EPS_KINK; //0.01;
  // Removed from corresponding use of fairingToolbox
  double approx = 0.001;


  ifstream file(argv[1]);

  CompositeModelFactory factory(approx, gap, neighbour, kink, 10.0*kink);

  shared_ptr<Go::CompositeModel> model(factory.createFromIges(file1));


  // Create the Surface model
  // Replaces:
  // ftFairingToolbox tool(approx, gap, neighbour, kink, 10.0*kink);

  ofstream os ("dump.g2");
  ofstream osPl ("dumpPlane.g2");

//   // Topology build
//   int ki;
//   ftMessage status;
//   status = tool->buildTopology();
//   std::cout << "Build topology. Status message : " << status.getMessage();
//   int nmbwarning = status.noOfWarnings();
//   std::cout << ", No of warnings : " << nmbwarning << std::endl;
//   for (ki=0; ki<nmbwarning; ki++)
//     std::cout << "Warning nr " << ki+1 << " : " << status.getWarning(ki) << std::endl;

  Point
    p_x (atof(argv[2]), atof(argv[3]), atof(argv[4])),
    p_y (atof(argv[5]), atof(argv[6]), atof(argv[7])),
    p_z (atof(argv[8]), atof(argv[9]), atof(argv[10]));

  ftPlane pl(p_x, p_y, p_z);

 osPl << "200 1 0 0" << std::endl;
  osPl << "3 0" << std:: endl;
  osPl << "2 2" << std::endl;
  osPl << "0 0 1 1" << std::endl;
  osPl << "2 2" << std::endl;
  osPl << "0 0 1 1" << std::endl;
  osPl << p_x << std::endl;
  osPl << p_x << std::endl;
  osPl << p_y << std::endl;
  osPl << p_z << std::endl;

  std::ofstream out("plane_ints.g2");
  double density = 1.0;
  std::vector<shared_ptr<LineStrip> > line_seg;
  PointCloud3D points;
  shared_ptr<IntResultsModel> results = model->intersect_plane(pl);

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

  shared_ptr<SurfaceModel> sfmodel = 
    dynamic_pointer_cast<SurfaceModel, CompositeModel>(model);
  if (sfmodel.get())
    {
//       shared_ptr<SurfaceModel> partmodel = sfmodel->trimWithPlane(pl);

//       std::ofstream out2("trimmed_model.g2");
//       int nmb = partmodel->nmbEntities();
//       for (int ki=0; ki<nmb; ++ki)
// 	{
// 	  shared_ptr<ParamSurface> sf = partmodel->getSurface(ki);
// 	  sf->writeStandardHeader(out2);
//  	  sf->write(out2);
// 	}

      sfmodel->booleanIntersect(pl);
      std::ofstream out3("trimmed_model2.g2");
      int nmb = sfmodel->nmbEntities();
      for (int ki=0; ki<nmb; ++ki)
	{
	  shared_ptr<ParamSurface> sf = sfmodel->getSurface(ki);
	  sf->writeStandardHeader(out3);
 	  sf->write(out3);
	}
    }
  shared_ptr<IntResultsCompCv> cvresults = 
    dynamic_pointer_cast<IntResultsCompCv, IntResultsModel>(results);

  if (cvresults.get())
    {
      std::cout << "Nmb int pt: " << cvresults->nmbIntPoints() << std::endl;
      vector<PointOnCurve> int_pt;
      cvresults->getIntersectionPoints(int_pt);
      for (size_t ki=0; ki<int_pt.size(); ++ki)
	std::cout << int_pt[ki].getPos() << std::endl;
    }
}
