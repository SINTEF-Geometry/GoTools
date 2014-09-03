#include <fstream>
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/GoTools.h"


using namespace std;
using namespace Go;


int main(int argc, char** argv)
{
  GoTools::init();

  if (argc != 3) {
    cout << "Usage:  " << argv[0] << " infile outfile" << endl;
    return 1;
  }

  ifstream ins(argv[1]);
  ObjectHeader header;
  BoundedSurface* in_surf = new BoundedSurface();
  ins >> header;
  ins >> *in_surf;
  ins.close();

  /*
  shared_ptr<Sphere> sphere = dynamic_pointer_cast<Sphere>(in_surf->underlyingSurface());

  CurveLoop loop = in_surf->outerBoundaryLoop();
  for (int i = 0; i < loop.size(); ++i)
    {
      shared_ptr<ParamCurve> curve = loop[i];
      shared_ptr<CurveOnSurface> cos = dynamic_pointer_cast<CurveOnSurface>(curve);
      shared_ptr<ParamCurve> geo_curve = cos->spaceCurve();
      shared_ptr<Circle> curve_circ = dynamic_pointer_cast<Circle>(geo_curve);
      if (curve_circ.get())
	{
	  cout << "Got circle at i = " << i << endl;
	  shared_ptr<ParamCurve> par_curve = sphere->getElementaryParamCurve(curve_circ.get(), 1.0e-3);
	  cos->setParameterCurve(par_curve);
	}
    }
  */

  ofstream outs(argv[2]);
  in_surf->writeStandardHeader(outs);
  in_surf->write(outs);
  outs.close();
}
