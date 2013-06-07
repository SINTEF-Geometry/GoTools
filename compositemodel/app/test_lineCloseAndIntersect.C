/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

#include <vcl.h>
#endif

#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/igeslib/IGESconverter.h"
#include "GoTools/compositemodel/ftMessage.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/igeslib/ftTangPriority.h"
#include <fstream>
#include <stdlib.h> // For atof()
#include <cstdlib>
#include <time.h>


const double EPS_GAP = 1e-7;
const double EPS_NEIGHBOUR = 1e-4;
const double EPS_KINK = 1e-2;
const double EPS_BEND = 4e-2;

const double EPS_APPROX = 1e-4;


using std::vector;
using std::istream;
using namespace Go;


ftMessage readIges(istream& is, vector<shared_ptr<ftSurface> >& faces)
//===========================================================================
{

  ftMessage status;
  if (is.bad())
    {
      status.setError(FT_BAD_INPUT_FILE);
      return status;
    }

  IGESconverter conv;
  try
    {
      conv.readIGES(is);
    }
  catch (...)
    {
      status.setError(FT_ERROR_IN_READ_IGES);
      return status;
    }

  //      std::ofstream outfile("debug.out");
  //      conv.writedisp(outfile);

  vector<shared_ptr<GeomObject> > gogeom = conv.getGoGeom();

  int nmbgeom = (int)gogeom.size();
  faces.reserve(nmbgeom); // May be too much, but not really important
  int face_count = 0;


  // Remaining geometry.

  for (int i=0; i<nmbgeom; i++)
    {

      if (gogeom[i]->instanceType() == Class_SplineCurve)
	{
	  if (conv.getGroup().size() == 0)
	    {
	      // One mesh surface expected.

	    }
	  else
	    // Not expected. Ignore the current curve.
	    status.addWarning(FT_UNEXPECTED_INPUT_OBJECT_IGNORED);  
	}
      else if (gogeom[i]->instanceType() == Class_SplineSurface ||
	       gogeom[i]->instanceType() == Class_BoundedSurface)
	{
	  shared_ptr<GeomObject> lg = gogeom[i];
	  shared_ptr<ParamSurface> gosf =
	    dynamic_pointer_cast<ParamSurface, GeomObject>(lg);
	  shared_ptr<ftSurface> ftsf(new ftSurface(gosf, face_count++));
	  faces.push_back(ftsf);

	}
    }

  return status;

}



double randd()
{
  return (double)(rand()) / (double)(RAND_MAX);
}


int main( int argc, char* argv[] )
{
#ifdef __BORLANDC__
  using Go::Point;
#endif

//   if (argc != 5) {
//     std::cout << "Input parameters : Input file on IGES format, random seed, iterators, checkin?" << std::endl;
//     exit(-1);
//   }

//   // Read input arguments
//   std::ifstream file1(argv[1]);
//   ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");


//   double gap = 0.001;
//   double neighbour = 0.01;
//   double kink = 0.01;
//   // Removed from corresponding use of fairingToolbox
//   double approx = 0.001;


//   ifstream file(argv[1]);

//   vector<shared_ptr<ftSurface> > faces;
//   ftMessage status = readIges(file1, faces);
//   std::cout << "Read IGES file. Status message : " << status.getMessage();
//   int nmbwarning = status.noOfWarnings();
//   std::cout << ", No of warnings : " << nmbwarning << std::endl;
//   int ki;
//   for (ki=0; ki<nmbwarning; ki++)
//     std::cout << "Warning nr " << ki+1 << " : " << status.getWarning(ki) << std::endl;

//   // Create the Surface model
//   // Replaces:
//   // ftFairingToolbox tool(approx, gap, neighbour, kink, 10.0*kink);
//   SurfaceModel tool(approx, gap, neighbour, kink, 10.0*kink, faces);

//   // Topology build
//   status = tool.buildTopology();
//   std::cout << "Build topology. Status message : " << status.getMessage();
//   nmbwarning = status.noOfWarnings();
//   std::cout << ", No of warnings : " << nmbwarning << std::endl;
//   for (ki=0; ki<nmbwarning; ki++)
//     std::cout << "Warning nr " << ki+1 << " : " << status.getWarning(ki) << std::endl;

//   // Check input geometry
//   int nmbboundaries = tool.nmbBoundaries();
//   std::cout << "Number of boundary loops in input geometry: " << nmbboundaries << std::endl;
//   if (nmbboundaries > 0)
//   {
//       std::ofstream out_bd("boundaries.g2");
//       for (ki=0; ki<nmbboundaries; ki++)
//       {
// 	  ftCurve bdcv = tool.getBoundary(ki);
// 	  bdcv.writeSpaceCurve(out_bd);
//       }
//   }
      

//   ftCurve gapcv = tool.getGaps();
//   int nmbg = gapcv.numDisjointSubcurves();
//   std::cout << "Number of gaps in input geometry: " << nmbg << std::endl;
//   if (nmbg > 0)
//   {
//       std::ofstream out_gap("gap.g2");
//       gapcv.writeSpaceCurve(out_gap);
//   }



//   ftCurve kinkcv = tool.getKinks();
//   nmbg = kinkcv.numDisjointSubcurves();
//   std::cout << "Number of kinks in input geometry: " << nmbg << std::endl;
//   if (nmbg > 0)
//   {
//       std::ofstream out_kink("kink.g2");
//       kinkcv.writeSpaceCurve(out_kink);
//   }



//   ftCurve g1dcv = tool.getG1Disconts();
//   nmbg = g1dcv.numDisjointSubcurves();
//   std::cout << "Number of G1 discontinuities in input geometry: " << nmbg << std::endl;
//   if (nmbg > 0)
//   {
//       std::ofstream out_g1d("g1_disc.g2");
//       g1dcv.writeSpaceCurve(out_g1d);
//   }


//   srand(atoi(argv[2]));

//   int
//     far_out = 0,
//     far_in = 0,
//     close_out = 0,
//     close_in = 0;

//   BoundingBox bb = tool.boundingBox();
//   // Point diagonal = (bb.high() - bb.low()) * 2;
//   // Point p0 = bb.low() - diagonal / 4;
//   Point diagonal = bb.high() - bb.low();
//   Point p0 = bb.low();

//   int nmbLoops = atoi(argv[3]);
//   int t = -clock();
//   if (atoi(argv[4]) == 0)   // Do both sphere and exact
//     for (int i = 0; i<nmbLoops; ++i)
//       {

// 	Point q0 = Point(p0[0] + diagonal[0] * randd(),
// 			 p0[1] + diagonal[1] * randd(),
// 			 p0[2] + diagonal[2] * randd());
// 	Point q1 = Point(p0[0] + diagonal[0] * randd(),
// 			 p0[1] + diagonal[1] * randd(),
// 			 p0[2] + diagonal[2] * randd());

// 	ftLine aLine(q1-q0,q0);

// 	/*
// 	  std::cout << "Line direction " << aLine.direction() << " and point " << aLine.point() << std::endl;
// 	  std::cout << "Bounding box limits: " << bb.low() << " and " << bb.high() << std::endl;
// 	*/

// 	tool.testLineIntersect(aLine, far_out, far_in, close_out, close_in);
//       }

//   else if (atoi(argv[4]) == 1)   // Do sphere only
//     for (int i = 0; i<nmbLoops; ++i)
//       {

// 	Point q0 = Point(p0[0] + diagonal[0] * randd(),
// 			 p0[1] + diagonal[1] * randd(),
// 			 p0[2] + diagonal[2] * randd());
// 	Point q1 = Point(p0[0] + diagonal[0] * randd(),
// 			 p0[1] + diagonal[1] * randd(),
// 			 p0[2] + diagonal[2] * randd());

// 	ftLine aLine(q1-q0,q0);

// 	/*
// 	  std::cout << "Line direction " << aLine.direction() << " and point " << aLine.point() << std::endl;
// 	  std::cout << "Bounding box limits: " << bb.low() << " and " << bb.high() << std::endl;
// 	*/

// 	tool.testLineIntersect(aLine, far_in, close_in);

//       }

//   else if (atoi(argv[4]) == 2)   // Do both planes and exact
//     for (int i = 0; i<nmbLoops; ++i)
//       {

// 	Point q0 = Point(p0[0] + diagonal[0] * randd(),
// 			 p0[1] + diagonal[1] * randd(),
// 			 p0[2] + diagonal[2] * randd());
// 	Point q1 = Point(p0[0] + diagonal[0] * randd(),
// 			 p0[1] + diagonal[1] * randd(),
// 			 p0[2] + diagonal[2] * randd());

// 	ftLine aLine(q1-q0,q0);

// 	/*
// 	  std::cout << "Line direction " << aLine.direction() << " and point " << aLine.point() << std::endl;
// 	  std::cout << "Bounding box limits: " << bb.low() << " and " << bb.high() << std::endl;
// 	*/

// 	tool.testLineIntersect2(aLine, far_out, far_in, close_out, close_in);
//       }

//   else if (atoi(argv[4]) == 3)   // Do both planes only
//     for (int i = 0; i<nmbLoops; ++i)
//       {

// 	Point q0 = Point(p0[0] + diagonal[0] * randd(),
// 			 p0[1] + diagonal[1] * randd(),
// 			 p0[2] + diagonal[2] * randd());
// 	Point q1 = Point(p0[0] + diagonal[0] * randd(),
// 			 p0[1] + diagonal[1] * randd(),
// 			 p0[2] + diagonal[2] * randd());

// 	ftLine aLine(q1-q0,q0);

// 	/*
// 	  std::cout << "Line direction " << aLine.direction() << " and point " << aLine.point() << std::endl;
// 	  std::cout << "Bounding box limits: " << bb.low() << " and " << bb.high() << std::endl;
// 	*/

// 	tool.testLineIntersect2(aLine, far_in, close_in);
//       }

//   else if (atoi(argv[4]) == 4)   // Do exact only
//     for (int i = 0; i<nmbLoops; ++i)
//       {

// 	Point q0 = Point(p0[0] + diagonal[0] * randd(),
// 			 p0[1] + diagonal[1] * randd(),
// 			 p0[2] + diagonal[2] * randd());
// 	Point q1 = Point(p0[0] + diagonal[0] * randd(),
// 			 p0[1] + diagonal[1] * randd(),
// 			 p0[2] + diagonal[2] * randd());

// 	ftLine aLine(q1-q0,q0);

// 	/*
// 	  std::cout << "Line direction " << aLine.direction() << " and point " << aLine.point() << std::endl;
// 	  std::cout << "Bounding box limits: " << bb.low() << " and " << bb.high() << std::endl;
// 	*/

// 	tool.testLineIntersect3(aLine, far_in, close_in);
//       }

//   else if (atoi(argv[4]) == 5)   // Do nothing
//     for (int i = 0; i<nmbLoops; ++i)
//       {

// 	Point q0 = Point(p0[0] + diagonal[0] * randd(),
// 			 p0[1] + diagonal[1] * randd(),
// 			 p0[2] + diagonal[2] * randd());
// 	Point q1 = Point(p0[0] + diagonal[0] * randd(),
// 			 p0[1] + diagonal[1] * randd(),
// 			 p0[2] + diagonal[2] * randd());

// 	ftLine aLine(q1-q0,q0);

// 	/*
// 	  std::cout << "Line direction " << aLine.direction() << " and point " << aLine.point() << std::endl;
// 	  std::cout << "Bounding box limits: " << bb.low() << " and " << bb.high() << std::endl;
// 	*/

//       }

//   t += clock();
//   if (atoi(argv[4]) == 0)
//     {
//       std::cout << "Inside  and close(s): " << close_in << std::endl;
//       std::cout << "Outside but close(s): " << close_out << std::endl;
//       std::cout << "Outside and away(s) : " << far_out << std::endl;
//       std::cout << "BAD:inside + away(s): " << far_in << std::endl;
//     }
//   else if (atoi(argv[4]) == 1)
//     {
//       std::cout << "Close(s) : " << close_in << std::endl;
//       std::cout << "Far(s)   : " << far_in << std::endl;
//     }
//   else if (atoi(argv[4]) == 2)
//     {
//       std::cout << "Inside  and close(p): " << close_in << std::endl;
//       std::cout << "Outside but close(p): " << close_out << std::endl;
//       std::cout << "Outside and away(p) : " << far_out << std::endl;
//       std::cout << "BAD:inside + away(p): " << far_in << std::endl;
//     }
//   else if (atoi(argv[4]) == 3)
//     {
//       std::cout << "Close(p) : " << close_in << std::endl;
//       std::cout << "Far(p)   : " << far_in << std::endl;
//     }
//   else if (atoi(argv[4]) == 4)
//     {
//       std::cout << "Inside  : " << close_in << std::endl;
//       std::cout << "Outside : " << far_in << std::endl;
//     }

//   std::setprecision(3);
//   std::cout << "Time used: " << ((double)t/1000000.0) << std::endl;

    return 0;

}
