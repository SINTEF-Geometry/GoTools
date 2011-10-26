//===========================================================================
//                                                                           
// File: test_Par2FuncIntersector.C                                          
//                                                                           
// Created: Wed Nov 17 08:59:05 2004                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision:
// $Id: test_Par2FuncIntersector.C,v 1.4 2006-03-29 14:33:33 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/intersections/Par2FuncIntersector.h"
#include "GoTools/intersections/Param2FunctionInt.h"
#include "GoTools/intersections/Spline2FunctionInt.h"
#include "GoTools/intersections/Param0FunctionInt.h"
#include "GoTools/geometry/RectDomain.h"
#include "GoTools/geometry/SplineUtils.h"
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
    cout << "Usage: spline2function C" << endl;
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
  shared_ptr<SplineSurface> surf(new SplineSurface());
  surf->read(filein);
  filein.close();
    
  shared_ptr<Spline2FunctionInt>
      spline2_func_int(new Spline2FunctionInt(surf));
  shared_ptr<Param0FunctionInt> param0_func_int(new Param0FunctionInt(C));

  shared_ptr<GeoTol> aepsge(new GeoTol(1.e-6));
  //const double aepsge = 1.e-6;
  Par2FuncIntersector par2_func_intersect(spline2_func_int,
					  param0_func_int, aepsge);

  par2_func_intersect.compute();
  std::vector<shared_ptr<IntersectionPoint> > int_points;
  std::vector<shared_ptr<IntersectionCurve> > int_curves;
  par2_func_intersect.getResult(int_points, int_curves);

  cout << "IntPoints found: " << int_points.size() << endl;
  cout << "IntCurves found: " << int_curves.size() << endl;

#ifdef INTERSECTIONS_DEBUG
  int ki;
  for (ki = 0; ki < int_curves.size(); ++ki) {
      int obj_nmb = 1;
      shared_ptr<SplineCurve> int_par_cv
	  = int_curves[ki]->getParamCurve(obj_nmb);
      std::ofstream debug("data/debug.g2");
      writeSpaceParamCurve(*int_par_cv, debug, C);
  }
#endif // INTERSECTIONS_DEBUG

  return 0;
}
