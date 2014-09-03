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

  ofstream outs(argv[2]);

  shared_ptr<Cylinder> under_surf = dynamic_pointer_cast<Cylinder>(in_surf->underlyingSurface());
  RectDomain rd = under_surf->parameterDomain();

  for (int i = 0; i < 2; ++i)
    {
      double from_upar = (i == 0) ? rd.umin() : M_PI;
      double to_upar = (i == 0) ? M_PI : rd.umax();
      double from_vpar = rd.vmin();
      double to_vpar = rd.vmax();
      cout << "Splitting to [" << from_upar << ", " << to_upar << "]x[" << from_vpar << ", " << to_vpar << "]" << endl;
      vector<shared_ptr<ParamSurface> > sub_surfs = in_surf->subSurfaces(from_upar, from_vpar, to_upar, to_vpar);
      cout << "Splitting completed, number of surfaces is " << sub_surfs.size() << endl;
      sub_surfs[0]->writeStandardHeader(outs);
      sub_surfs[0]->write(outs);
    }
  outs.close();

}
