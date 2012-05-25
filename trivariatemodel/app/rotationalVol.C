#include <fstream>
#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/trivariatemodel/VolumeModelCreator.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"

using namespace Go;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::vector;

int main(int argc, char* argv[] )
{
  if (argc != 3)
      cout << "Usage: " << "<infile> <outfile>" << endl;

  ifstream infile(argv[1]);
  ALWAYS_ERROR_IF(infile.bad(), "Bad or no input filename");

  ofstream outfile(argv[2]);

  double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  CompositeModel *model = factory.createFromG2(infile);

  shared_ptr<SurfaceModel> sfmodel = 
    shared_ptr<SurfaceModel>(dynamic_cast<SurfaceModel*>(model));
  if (!sfmodel.get())
    exit(-1);
 
  shared_ptr<VolumeModel> volmod;
  bool rotational = VolumeModelCreator::createRotationalModel(sfmodel, volmod);
  std::cout << "Rotational: " << rotational << std::endl;

  if (rotational)
    {
        int nmb_vols = volmod->nmbEntities();
	for (int kr=0; kr<nmb_vols; ++kr)
	  {
	    shared_ptr<ParamVolume> curr_vol2 = volmod->getVolume(kr);
	    curr_vol2->writeStandardHeader(outfile);
	    curr_vol2->write(outfile);
	  }
    }

}

