#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/trivariate/ParamVolume.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <fstream>

using namespace Go;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::shared_ptr;
using std::dynamic_pointer_cast;

int main(int argc, char* argv[] )
{
  if (argc != 3)
      cout << "Usage: " << "infile2 outfile" << endl;

  ifstream is2(argv[1]);
  ALWAYS_ERROR_IF(is2.bad(), "Bad or no input filename");

  ofstream of(argv[2]);

  vector<shared_ptr<ftVolume> > volumes;
  
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

      eatwhite(is2);
    }

  double gap = 0.000001;
  double neighbour = 0.00001;
  double kink = 0.01;
  double bend = 0.05;
  shared_ptr<VolumeModel> model = 
    shared_ptr<VolumeModel>(new VolumeModel(volumes, gap, neighbour, 
					    kink, bend));

  int nmb;
  std::cout << "Number of volumes: ";
  std::cin >> nmb;

  std::cout << "Give first index: ";
  for (int ki=0; ki<nmb; ++ki)
    {
      int idx;
      std::cin >> idx;
      shared_ptr<ParamVolume> curr = model->getVolume(idx);

      shared_ptr<SplineVolume> spl =
	dynamic_pointer_cast<SplineVolume, ParamVolume>(curr);
      if (spl.get())
	{
	  std::cout << "Swap parameter directions? ";
	  int swap;
	  std::cin >> swap;
	  if (swap)
	    {
	      int dir1, dir2;
	      std::cout << "The two directions to swap: ";
	      std::cin >> dir1;
	      std::cin >> dir2;
	      spl->swapParameterDirection(dir1, dir2);
	    }

	  std::cout << "Give parameter directions to turn: ";
	  int dir1;
	  std::cin >> dir1;
	  while (dir1 >= 0 && dir1 <= 2)
	    {
	      spl->reverseParameterDirection(dir1);
	      std::cin >> dir1;
	    }
	  std::cout << "Give next index: ";
	}
	
      curr->writeStandardHeader(of);
      curr->write(of);
    }
}
