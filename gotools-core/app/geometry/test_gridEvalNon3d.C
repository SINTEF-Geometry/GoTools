//===========================================================================
//
// File : test_gridEvalNon3d.C
//
// Created: Tue Aug 11 11:30:51 2009
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================



#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/Point.h"
#include <fstream>


using namespace Go;
using namespace std;


int main(int argc, char** argv)
{

  if (argc != 6)
    {
      cout << "Usage: " << argv[0] << " surfaceinfile surface3doutfile points3doutfile num_u num_v" << endl;
      exit(-1);
    }

  ifstream filein(argv[1]);
  ALWAYS_ERROR_IF(filein.bad(), "Bad or no curvee input filename");
  ObjectHeader head;
  filein >> head;
  if (head.classType() != SplineSurface::classType()) {
    THROW("Not a spline surface");
  }

  SplineSurface sf;
  filein >> sf;

  ofstream fileoutsurf(argv[2]);
  ALWAYS_ERROR_IF(fileoutsurf.bad(), "Bad surface output filename");

  ofstream fileoutpts(argv[3]);
  ALWAYS_ERROR_IF(fileoutpts.bad(), "Bad points output filename");

  int num_u = atoi(argv[4]);
  int num_v = atoi(argv[5]);

  vector<double> pts, param_u, param_v;

  sf.gridEvaluator(num_u, num_v, pts, param_u, param_v);

  vector<double> coefs3d;
  vector<Point> pts3d;
  int dim = sf.dimension();
  bool rational = sf.rational();

  int ctrl_pts = sf.numCoefs_u() * sf.numCoefs_v();
  vector<double>::const_iterator it = sf.ctrl_begin();

  for (int i = 0; i < ctrl_pts; ++i)
    {
      if (dim <= 3)
	for (int j = 0; j < 3; ++j)
	  {
	    if (j>=dim)
	      coefs3d.push_back(0.0);
	    else
	      {
		coefs3d.push_back(*it);
		++it;
	      }
	  }
      else
	{
	  for (int j = 0; j < 3; ++j, ++it)
	    coefs3d.push_back(*it);
	  it += (dim-3);
	}
      if (rational)
	{
	  coefs3d.push_back(*it);
	  ++it;
	}
    }

  int pts_pos = 0;
  for (int i = 0; i < num_u*num_v; ++i)
    {
      double x, y, z;
      if (dim == 0)
	x = 0.0;
      else
	x = pts[pts_pos];
      if (dim <= 1)
	y = 0.0;
      else
	y = pts[pts_pos+1];
      if (dim <= 2)
	z = 0.0;
      else
	z = pts[pts_pos+2];
      pts_pos += dim;
      pts3d.push_back(Point(x, y, z));
    }

  SplineSurface sf3d(sf.basis_u(), sf.basis_v(), coefs3d.begin(), 3, rational);

  sf3d.writeStandardHeader(fileoutsurf);
  sf3d.write(fileoutsurf);

  fileoutpts << "400 1 0 4 255 255 0 255" << endl;
  fileoutpts << pts3d.size() << endl;
  for (int i = 0; i < (int)pts3d.size(); ++i)
    fileoutpts << pts3d[i] << endl;
}

