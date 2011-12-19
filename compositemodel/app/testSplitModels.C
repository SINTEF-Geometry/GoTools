#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/ftCurve.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include <fstream>

using std::vector;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 3) {
    std::cout << "Input parameters : File 1 on g2 format, File 2 on g2 format" << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  std::ifstream file2(argv[2]);
  ALWAYS_ERROR_IF(file2.bad(), "Input file not found or file corrupt");

    double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  shared_ptr<CompositeModel> model1 = shared_ptr<CompositeModel>(factory.createFromG2(file1));
  shared_ptr<SurfaceModel> sfmodel1 = dynamic_pointer_cast<SurfaceModel,CompositeModel>(model1);

  shared_ptr<CompositeModel> model2 = shared_ptr<CompositeModel>(factory.createFromG2(file2));
  shared_ptr<SurfaceModel> sfmodel2 = dynamic_pointer_cast<SurfaceModel,CompositeModel>(model2);

  if (sfmodel1.get() && sfmodel2.get())
    {
      size_t ki;
      vector<shared_ptr<SurfaceModel> > res = 
	sfmodel1->splitSurfaceModels(sfmodel2);
      vector<shared_ptr<SurfaceModel> > res2;
      for (ki=0; ki<res.size(); ++ki)
	{
	  vector<shared_ptr<SurfaceModel> >res3 = res[ki]->getConnectedModels();
	  res2.insert(res2.end(), res3.begin(), res3.end());
	}

      std::cout << "Number of models: " << res2.size() << std::endl;
      for (ki=0; ki<res2.size(); ++ki)
	{
	  std::cout << "File nr " << ki << ": " << std::endl;
	  char outfile[80];
	  std::cin >> outfile;
	  std::ofstream of1(outfile);

	  int nmb = res2[ki]->nmbEntities();
	  for (int kj=0; kj<nmb; kj++)
	    {
	      shared_ptr<ParamSurface> surf = res2[ki]->getSurface(kj);
	      
	      surf->writeStandardHeader(of1);
	      surf->write(of1);
	    }
	}
    }
}

