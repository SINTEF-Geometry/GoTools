#include <fstream>
#include <sstream>
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/GoTools.h"


using namespace std;
using namespace Go;


int main(int argc, char** argv)
{
  GoTools::init();

  if (argc != 3) {
    cout << "Usage:  " << argv[0] << " nmbSurfcs infile" << endl;
    return 1;
  }

  int nmb_surfs = atoi(argv[1]);
  vector<BoundedSurface*> surfs;

  ifstream ins(argv[2]);

  ObjectHeader header;
  for (int i = 0; i < nmb_surfs; ++i)
    {
      BoundedSurface* bs = new BoundedSurface();
      ins >> header >> *bs;
      surfs.push_back(bs);
    }
  ins.close();

  for (int i = 0; i < nmb_surfs; ++i)
    {
      CurveBoundedDomain cbd = surfs[i]->parameterDomain();
      RectDomain rd = cbd.containingDomain();
      cout << "Domain[" << i << "] can be limited to [" << rd.umin() << ", " << rd.umax() << "]x[" << rd.vmin() << ", " << rd.vmax() << "]" << endl;
    }
}

