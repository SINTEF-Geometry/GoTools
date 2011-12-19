#ifdef __BORLANDC__
#include <vcl.h>
#endif

#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeCurve.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include <fstream>

using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 2) {
    std::cout << "Input parameters : Input file on IGES format," << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

    double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  vector<shared_ptr<CompositeModel> > models = 
    factory.getModelsFromIges(file1);

  std::cout << "Number of models: " << models.size() << std::endl;
  
  std::ofstream out_file("iges2.g2");
  for (size_t kj=0; kj<models.size(); ++kj)
    {
      SurfaceModel *sfmodel = dynamic_cast<SurfaceModel*>(models[kj].get());

       CompositeCurve *cvmodel = 
	 dynamic_cast<CompositeCurve*>(models[kj].get());

       if (sfmodel)
	 {
	   int nmb = sfmodel->nmbEntities();
	   std::cout << kj << "sfmodel " << nmb << std::endl;
	   for (int ki=0; ki<nmb; ki++)
	     {
	       shared_ptr<ParamSurface> surf = sfmodel->getSurface(ki);

	       surf->writeStandardHeader(out_file);
	       surf->write(out_file);
	     }
	 }
	 if (cvmodel)
	   {
	     int nmb = cvmodel->nmbEntities();
	     std::cout << kj << "cvmodel " << nmb << std::endl;
	     for (int ki=0; ki<nmb; ki++)
	       {
		 shared_ptr<ParamCurve> crv = cvmodel->getCurve(ki);

		 crv->writeStandardHeader(out_file);
		 crv->write(out_file);
	       }
	   }

	 }
}
