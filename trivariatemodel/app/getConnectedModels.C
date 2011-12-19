#include <fstream>
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/SurfaceOnVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/geometry/Utils.h"

using namespace Go;


int main( int argc, char* argv[] )
{
  if (argc != 2) {
    std::cout << "Input parameters : Input file on g2 format" << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  vector<shared_ptr<ftVolume> > volumes;
  
  int ki;
  while (!file1.eof())
    {
      // Read volume from file
      ObjectHeader head;
      file1 >> head;
      shared_ptr<SplineVolume> vol2;
      vol2 = shared_ptr<SplineVolume>(new SplineVolume());
      vol2->read(file1);

      shared_ptr<ParamVolume> pvol
          = dynamic_pointer_cast<ParamVolume, SplineVolume>(vol2);
      volumes.push_back(shared_ptr<ftVolume>(new ftVolume(pvol)));
      // volumes.push_back(shared_ptr<ftVolume>(new ftVolume(vol2)));

      eatwhite(file1);
    }
    double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double bend = 0.05;

  shared_ptr<VolumeModel> model = 
    shared_ptr<VolumeModel>(new VolumeModel(volumes, gap, neighbour, 
					    kink, bend));
  std::ofstream of1("vol_boundary.g2");

  vector<shared_ptr<VolumeModel> > models2 =
    model->getConnectedModels();

  std::cout << "Number of connected models: " << models2.size() << std::endl;

  
  // Outer boundaries
  int kj;
  int nmb_bd = model->nmbBoundaries();
  std::cout << "Outer boundaries: " << nmb_bd << std::endl;

  for (ki=0; ki<(int)models2.size(); ++ki)
    {
      std::ofstream of2("connected.g2");
      int nmb = models2[ki]->nmbEntities();
      for (kj=0; kj<nmb; kj++)
      {
	  shared_ptr<SplineVolume> vol = models2[ki]->getSplineVolume(kj);

	  vol->writeStandardHeader(of2);
	  vol->write(of2);
      }
    }
}
