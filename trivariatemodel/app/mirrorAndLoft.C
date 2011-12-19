#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/trivariate/LoftVolumeCreator.h"
#include <fstream>

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

  std::ofstream file3("loftvol.g2");

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
  vector<shared_ptr<ParamSurface> > mirrored;
  vector<shared_ptr<ftVolume> > blocks;
  for (int ki=0; ki<nmb_sfs; ++ki)
    {
      shared_ptr<ParamSurface> curr = sfmodel->getSurface(ki);
      shared_ptr<ParamSurface> mod = 
	shared_ptr<ParamSurface>(curr->mirrorSurface(pnt, norm));
      mirrored.push_back(mod);

      vector<shared_ptr<SplineSurface> > sfs(2);
      sfs[0] = dynamic_pointer_cast<SplineSurface,ParamSurface>(curr);
      sfs[1] = dynamic_pointer_cast<SplineSurface,ParamSurface>(mod);

      if (sfs[0].get() && sfs[1].get())
	{
	  shared_ptr<ParamVolume> vol = 
	    shared_ptr<ParamVolume>(LoftVolumeCreator::loftVolume(sfs.begin(), 2));

	  vol->writeStandardHeader(file3);
	  vol->write(file3);

	  shared_ptr<ftVolume> ftvol = 
	    shared_ptr<ftVolume>(new ftVolume(vol, gap, kink));
	  blocks.push_back(ftvol);
	}
    }

  shared_ptr<VolumeModel> volmodel = 
    shared_ptr<VolumeModel>(new VolumeModel(blocks, gap, neighbour, 
					    kink, 10.0*kink));


  int nmb_vol = volmodel->nmbEntities();
  for (int ki=0; ki<nmb_vol; ++ki)
    {
      shared_ptr<ParamVolume> vol = volmodel->getVolume(ki);
       vol->writeStandardHeader(file2);
       vol->write(file2);
     }
}
  
