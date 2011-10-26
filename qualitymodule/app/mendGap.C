#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/qualitymodule/FaceSetQuality.h"
#include "GoTools/qualitymodule/FaceSetRepair.h"
#include <fstream>

using std::cout;
using std::cin;
using std::shared_ptr;
using std::dynamic_pointer_cast;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 2 && argc != 3) {
    std::cout << "Input parameters : Input file on g2 format, (Insert knots)" << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

    double gap = 0.0001;
    double neighbour = 0.001;
  double kink = 0.01;
  double approxtol = 0.01;
  int insert = 0;
  if (argc == 3)
    insert = atoi(argv[2]);

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  shared_ptr<CompositeModel> model = shared_ptr<CompositeModel>(factory.createFromG2(file1));
  shared_ptr<SurfaceModel> sfmodel = dynamic_pointer_cast<SurfaceModel,CompositeModel>(model);

  if (insert)
    {
      cout << "Number of surfaces: " << sfmodel->nmbEntities() << std::endl;
      cout << "Surface to refine: ";
      int idx;
      cin >> idx;

      cout << "Parameter direction: ";
      int dir;
      cin >> dir;

      cout << "Number of knots: ";
      int nmb;
      cin >> nmb;

      cout << "Knots: ";
      vector<double> knots(nmb);
      for (int ki=0; ki<nmb; ++ki)
	cin >> knots[ki];

      shared_ptr<ParamSurface> srf = sfmodel->getSurface(idx);
      shared_ptr<SplineSurface> s1 = 
	dynamic_pointer_cast<SplineSurface,ParamSurface>(srf);
      if (s1.get())
	{
	  if (dir == 0)
	    s1->insertKnot_u(knots);
	  else
	    s1->insertKnot_v(knots);
	}
    }
      
  shared_ptr<FaceSetQuality> quality = 
    shared_ptr<FaceSetQuality>(new FaceSetQuality(gap, kink, approxtol));
  quality->attach(sfmodel);

  shared_ptr<FaceSetRepair> repair = 
    shared_ptr<FaceSetRepair>(new FaceSetRepair(quality));
 
  vector<pair<ftEdge*, ftEdge*> > gaps;
  quality->facePositionDiscontinuity(gaps);  

  std::cout << "Number of gaps: " << gaps.size() << std::endl;

  std::ofstream out_file("gaps_1.g2");
  size_t ki;
  for (ki=0; ki<gaps.size(); ++ki)
  {
    shared_ptr<ParamCurve> cv1 = gaps[ki].first->geomCurve();
    shared_ptr<ParamCurve> cv2 = gaps[ki].second->geomCurve();
    shared_ptr<SplineCurve> scv1 = shared_ptr<SplineCurve>(cv1->geometryCurve());
    shared_ptr<SplineCurve> scv2 = shared_ptr<SplineCurve>(cv2->geometryCurve());
    scv1->writeStandardHeader(out_file);
    scv1->write(out_file);
    scv2->writeStandardHeader(out_file);
    scv2->write(out_file);
  }

  repair->mendEdgeDistance();  

  std::ofstream out_model0("sfmodel0.g2");
  int nmb = sfmodel->nmbEntities();
  for (int ki=0; ki<nmb; ki++)
    {
      shared_ptr<ParamSurface> surf = sfmodel->getSurface(ki);

      surf->writeStandardHeader(out_model0);
      surf->write(out_model0);
    }

  repair->mendGaps();

  std::ofstream out_model("sfmodel.g2");
  nmb = sfmodel->nmbEntities();
  for (int ki=0; ki<nmb; ki++)
    {
      shared_ptr<ParamSurface> surf = sfmodel->getSurface(ki);

      surf->writeStandardHeader(out_model);
      surf->write(out_model);
    }

			      
  vector<pair<ftEdge*, ftEdge*> > gaps2;
  quality->facePositionDiscontinuity(gaps2);  

  std::cout << "Number of gaps after mending: " << gaps2.size() << std::endl;

  std::ofstream out_filen("gaps_2.g2");
  for (ki=0; ki<gaps2.size(); ++ki)
  {
    shared_ptr<ParamCurve> cv1 = gaps2[ki].first->geomCurve();
    shared_ptr<ParamCurve> cv2 = gaps2[ki].second->geomCurve();
    shared_ptr<SplineCurve> scv1 = shared_ptr<SplineCurve>(cv1->geometryCurve());
    shared_ptr<SplineCurve> scv2 = shared_ptr<SplineCurve>(cv2->geometryCurve());
    scv1->writeStandardHeader(out_filen);
    scv1->write(out_filen);
    scv2->writeStandardHeader(out_filen);
    scv2->write(out_filen);
  }

}

