#include <fstream>
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/SurfaceOnVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/geometry/Utils.h"

using namespace Go;
using std::cout;
using std::endl;
using std::ifstream;
using std::vector;

int main(int argc, char* argv[] )
{
  if (argc != 2 && argc != 3)
      cout << "Usage: " << "infile2, (Insert knots)" << endl;

  ifstream is2(argv[1]);
  ALWAYS_ERROR_IF(is2.bad(), "Bad or no input filename");
  int insert = 0;
  if (argc == 3)
    insert = atoi(argv[2]);

  vector<shared_ptr<ftVolume> > volumes;
  
  //int ki;
  while (!is2.eof())
    {
      // Read volume from file
      ObjectHeader head;
      is2 >> head;
      shared_ptr<SplineVolume> vol2;
      vol2 = shared_ptr<SplineVolume>(new SplineVolume());
      vol2->read(is2);

      shared_ptr<ParamVolume> pvol
          = dynamic_pointer_cast<ParamVolume, SplineVolume>(vol2);
      volumes.push_back(shared_ptr<ftVolume>(new ftVolume(pvol)));
      // volumes.push_back(shared_ptr<ftVolume>(new ftVolume(vol2)));

      Utils::eatwhite(is2);
    }

  double gap = 0.0001;
  double neighbour = 0.01;
  double kink = 0.01;
  double bend = 0.05;
  shared_ptr<VolumeModel> model = 
    shared_ptr<VolumeModel>(new VolumeModel(volumes, gap, neighbour, 
					    kink, bend));

  if (model)
    {
      std::ofstream out_file("vol_bd.g2");
      int nmb = model->nmbBoundaries();
      std::cout << "Number of outer boundaries: " << nmb << std::endl;
      for (int ki=0; ki<nmb; ++ki)
	{
	  shared_ptr<SurfaceModel> curr = model->getOuterBoundary(ki);
	  int nmb_face = curr->nmbEntities();
	  for (int kj=0; kj<nmb_face; ++kj)
	    {
	      shared_ptr<ParamSurface> surf = curr->getSurface(kj);
	      shared_ptr<SurfaceOnVolume> volsf = 
		dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(surf);
	      if (volsf.get())
		surf = volsf->spaceSurface();

	      surf->writeStandardHeader(out_file);
	      surf->write(out_file);
	    }
	}
    }
}

