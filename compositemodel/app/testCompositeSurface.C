//===========================================================================
//
// File : testCompositeModel.C
//
// Created: Thu Feb 21 09:37:14 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: testCompositeSurface.C,v 1.3 2009-05-13 07:29:53 vsk Exp $
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
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/igeslib/ftTangPriority.h"
#include <istream>
#include <fstream>
#include <stdlib.h> // For atof()


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



int main( int argc, char* argv[] )
{
  // Test number of input arguments
  if (argc != 8)
    {
      std::cout << "Input arguments : Input file on IGES format, ";
      std::cout << "gap tolerance, neighbourhood tolerance. kink tolerance, approximation tolerance, ";
      std::cout << "boundary?,  output file" << std::endl;

      exit(-1);
    }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  double gap = atof(argv[2]);
  double neighbour = atof(argv[3]);
  double kink = atof(argv[4]);
  // Removed from corresponding use of fairingToolbox
  double approx = atof(argv[5]);
  int check_bd;
  check_bd = atoi(argv[6]);

  std::ofstream file2(argv[7]);

  vector<shared_ptr<ftSurface> > faces;
  ftMessage status = readIges(file1, faces);
  std::cout << "Read IGES file. Status message : " << status.getMessage();
  int nmbwarning = status.noOfWarnings();
  std::cout << ", No of warnings : " << nmbwarning << std::endl;
  int ki;
  for (ki=0; ki<nmbwarning; ki++)
    std::cout << "Warning nr " << ki+1 << " : " << status.getWarning(ki) << std::endl;

  // Create the Surface model
  // Replaces:
  // ftFairingToolbox tool(approx, gap, neighbour, kink, 10.0*kink);
  SurfaceModel tool(approx, gap, neighbour, kink, 10.0*kink, faces);

  // Left out from fairingToolbox
  // Make clean degenerate surfaces and remove surfaces degenerated to a line
  // int nmb_deg=0, nmb_bd=0, nmb_removed=0;
  // status = tool.ensureCleanDegeneracy(nmb_deg, nmb_bd, nmb_removed);
  // std::cout << "Clean degeneracy. Status message : " << status.getMessage() << std::endl;
  // std::cout << "Number of degenerate triangular or banana surfaces: " << nmb_deg << std::endl;
  // std::cout << "Number of modified degenerate boundaries: " << nmb_bd << std::endl;
  // std::cout << "Number of removed surfaces: " << nmb_removed << std::endl;
  
  // Topology build
  status = tool.buildTopology();
  std::cout << "Build topology. Status message : " << status.getMessage();
  nmbwarning = status.noOfWarnings();
  std::cout << ", No of warnings : " << nmbwarning << std::endl;
  for (ki=0; ki<nmbwarning; ki++)
    std::cout << "Warning nr " << ki+1 << " : " << status.getWarning(ki) << std::endl;

  // Check input geometry
  int nmbboundaries = tool.nmbBoundaries();
  std::cout << "Number of boundary loops in input geometry: " << nmbboundaries << std::endl;
  if (nmbboundaries > 0)
  {
      std::ofstream out_bd("boundaries.g2");
      for (ki=0; ki<nmbboundaries; ki++)
      {
	  ftCurve bdcv = tool.getBoundary(ki);
	  bdcv.writeSpaceCurve(out_bd);
      }
  }
      

  ftCurve gapcv = tool.getGaps();
  int nmbg = gapcv.numDisjointSubcurves();
  std::cout << "Number of gaps in input geometry: " << nmbg << std::endl;
  if (nmbg > 0)
  {
      std::ofstream out_gap("gap.g2");
      gapcv.writeSpaceCurve(out_gap);
  }

  /*for (ki=0; ki<nmbg; ki++)
    {
      double t1 = gapcv.startOfSegment(ki);
      double t2 = gapcv.endOfSegment(ki);
      double par, del;
      del = (t2 - t1)/(double)2;
      Point gappt;
      for (kj=0, par=t1; kj<3; kj++, par+=del)
	{
	  gapcv.point(par, ki, gappt);
	  std::cout << "Point nr " << ki+1 << " : ";
	  std::cout << gappt[0] << ", " << gappt[1];
	  std::cout << ", " << gappt[2] << std::endl;
	}
	}*/

  ftCurve kinkcv = tool.getKinks();
  nmbg = kinkcv.numDisjointSubcurves();
  std::cout << "Number of kinks in input geometry: " << nmbg << std::endl;
  if (nmbg > 0)
  {
      std::ofstream out_kink("kink.g2");
      kinkcv.writeSpaceCurve(out_kink);
  }

  /*for (ki=0; ki<nmbg; ki++)
    {
      double t1 = kinkcv.startOfSegment(ki);
      double t2 = kinkcv.endOfSegment(ki);
      double par, del;
      del = (t2 - t1)/(double)2;
      Point kinkpt;
      for (kj=0, par=t1; kj<3; kj++, par+=del)
	{
	  kinkcv.point(par, ki, kinkpt);
	  std::cout << "Point nr " << ki+1 << " : ";
	  std::cout << kinkpt[0] << ", " << kinkpt[1];
	  std::cout << ", " << kinkpt[2] << std::endl;
	}
	}*/

  ftCurve g1dcv = tool.getG1Disconts();
  nmbg = g1dcv.numDisjointSubcurves();
  std::cout << "Number of G1 discontinuities in input geometry: " << nmbg << std::endl;
  if (nmbg > 0)
  {
      std::ofstream out_g1d("g1_disc.g2");
      g1dcv.writeSpaceCurve(out_g1d);
  }

  /*for (ki=0; ki<nmbg; ki++)
    {
      double t1 = g1dcv.startOfSegment(ki);
      double t2 = g1dcv.endOfSegment(ki);
      double par, del;
      del = (t2 - t1)/(double)2;
      Point g1dpt;
      for (kj=0, par=t1; kj<3; kj++, par+=del)
	{
	  g1dcv.point(par, ki, g1dpt);
	  std::cout << "Point nr " << ki+1 << " : ";
	  std::cout << g1dpt[0] << ", " << g1dpt[1];
	  std::cout << ", " << g1dpt[2] << std::endl;
	}
	}*/


  
}

















