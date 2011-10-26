//===========================================================================
//                                                                           
// File: test_Par1FuncIntersector.C                                          
//                                                                           
// Created: Tue Oct  5 13:32:50 2004                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision:
// $Id: test_Par1FuncIntersector.C,v 1.4 2006-03-13 17:01:22 jbt Exp $
//                                                                           
// Description: Given input of 1 parameter spline function (R -> R) and a
//              double, compute intersections.
//                                                                           
//===========================================================================


#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/intersections/Par1FuncIntersector.h"
#include "GoTools/intersections/Param1FunctionInt.h"
#include "GoTools/intersections/Spline1FunctionInt.h"
#include "GoTools/intersections/Param0FunctionInt.h"
#include "GoTools/intersections/GeoTol.h"
#include "GoTools/geometry/SplineCurve.h"
#include <memory>
#include <fstream>
#include <iomanip>


using std::cout;
using std::endl;
using std::ifstream;
using std::cerr;
using namespace Go;
using std::shared_ptr;


int main(int argc, char** argv)
{
 
  if (argc != 3) {
    cout << "Usage: spline1function C" << endl;
    return 0;
  }

  ObjectHeader header;

  // Read the first curve from file
  ifstream filein(argv[1]);
  if (filein.bad()) {
    cerr << "File #1 error (no file or corrupt file specified)."
	 << std::endl;
    return 1;
  }
  double C(atof(argv[2]));

  header.read(filein);
  shared_ptr<ParamCurve> curve(new SplineCurve());
  curve->read(filein);
  filein.close();
    
  shared_ptr<Spline1FunctionInt> spline1_func_int(new Spline1FunctionInt(curve));
  shared_ptr<Param0FunctionInt> param0_func_int(new Param0FunctionInt(C));

  shared_ptr<GeoTol> aepsge(new GeoTol(1.e-6));
  //const double aepsge = 1.e-6;
  Par1FuncIntersector par1_func_intersect(spline1_func_int, param0_func_int, aepsge);

  par1_func_intersect.compute();
  std::vector<shared_ptr<IntersectionPoint> > int_points;
  std::vector<shared_ptr<IntersectionCurve> > int_curves;
  par1_func_intersect.getResult(int_points, int_curves);

  cout << "IntPoints found: " << int_points.size() << endl;
  cout << "IntCurves found: " << int_curves.size() << endl;

  return 0;
}
