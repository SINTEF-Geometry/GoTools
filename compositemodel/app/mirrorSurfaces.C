#include "GoTools/compositemodel/CompositeModelFactory.h"

#include "GoTools/geometry/ParamSurface.h"
#include <fstream>

//using namespace std;
using std::shared_ptr;
using std::dynamic_pointer_cast;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 9) {
    std::cout << "Input parameters : Input file on g2 format, output file, point in plane, normal" << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  std::ofstream file2(argv[2]);

  Point pnt(atof(argv[3]), atof(argv[4]), atof(argv[5]));
  Point norm(atof(argv[6]), atof(argv[7]), atof(argv[8]));

  double gap = 0.0001;
  double neighbour = 0.01;
  double kink = 0.01;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  shared_ptr<CompositeModel> model = shared_ptr<CompositeModel>(factory.createFromG2(file1));
  shared_ptr<SurfaceModel> sfmodel = dynamic_pointer_cast<SurfaceModel,CompositeModel>(model);

  int nmb_sfs = sfmodel->nmbEntities();
  std::vector<shared_ptr<ParamSurface> > mirrored;
  for (int ki=0; ki<nmb_sfs; ++ki)
    {
      shared_ptr<ParamSurface> curr = sfmodel->getSurface(ki);
      shared_ptr<ParamSurface> mod = 
	shared_ptr<ParamSurface>(curr->mirrorSurface(pnt, norm));
      mirrored.push_back(mod);
    }

  shared_ptr<SurfaceModel> mirrormodel = 
    shared_ptr<SurfaceModel>(new SurfaceModel(approxtol, gap, neighbour, 
					      kink, 10.0*kink, mirrored));

  for (int ki=0; ki<nmb_sfs; ++ki)
    {
      shared_ptr<ParamSurface> srf = mirrormodel->getSurface(ki);
       srf->writeStandardHeader(file2);
       srf->write(file2);
     }
}
  
