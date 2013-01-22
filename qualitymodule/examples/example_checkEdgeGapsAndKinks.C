#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/qualitymodule/FaceSetQuality.h"
#include <fstream>

using std::vector;
using std::pair;
using namespace Go;

int main( int argc, char* argv[] )
{
  // Read input arguments
  std::ifstream file1("data/DemEx6.g2");
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  double gap = 0.001;
  double neighbour = 0.002;
  double kink = 0.01;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);
  //CompositeModelFactory factory(approxtol, gap, neighbour, kink, 100.0*kink);

  shared_ptr<CompositeModel> model = shared_ptr<CompositeModel>(factory.createFromG2(file1));
  shared_ptr<SurfaceModel> sfmodel = dynamic_pointer_cast<SurfaceModel,CompositeModel>(model);

  FaceSetQuality quality(gap, kink, approxtol);
  quality.attach(sfmodel);

  vector<pair<ftEdge*, ftEdge*> > gaps;
  vector<pair<ftEdge*, ftEdge*> > kinks;
  quality.edgePosAndTangDiscontinuity(gaps, kinks);  

  std::cout << "Number of gaps: " << gaps.size() << std::endl;
  std::cout << "Number of kinks: " << kinks.size() << std::endl;

  std::ofstream out_file("edge_gaps.g2");
  size_t ki;
  for (ki=0; ki<gaps.size(); ++ki)
  {
      shared_ptr<ParamCurve> crv1 = gaps[ki].first->geomCurve();
      shared_ptr<CurveOnSurface> sfcv = 
	dynamic_pointer_cast<CurveOnSurface,ParamCurve>(crv1);
      if (sfcv.get())
	{
	  sfcv->spaceCurve()->writeStandardHeader(out_file);
	  sfcv->spaceCurve()->write(out_file);
	}
      else
	{
	  crv1->writeStandardHeader(out_file);
	  crv1->write(out_file);
	}

      crv1 = gaps[ki].second->geomCurve();
      sfcv = dynamic_pointer_cast<CurveOnSurface,ParamCurve>(crv1);
      if (sfcv.get())
	{
	  sfcv->spaceCurve()->writeStandardHeader(out_file);
	  sfcv->spaceCurve()->write(out_file);
	}
      else
	{
	  crv1->writeStandardHeader(out_file);
	  crv1->write(out_file);
	}
  }

  std::ofstream out_file2("edge_kinks.g2");
  for (ki=0; ki<kinks.size(); ++ki)
  {
      shared_ptr<ParamCurve> crv1 = kinks[ki].first->geomCurve();
      shared_ptr<CurveOnSurface> sfcv = 
	dynamic_pointer_cast<CurveOnSurface,ParamCurve>(crv1);
      if (sfcv.get())
	{
	  sfcv->spaceCurve()->writeStandardHeader(out_file2);
	  sfcv->spaceCurve()->write(out_file2);
	}
      else
	{
	  crv1->writeStandardHeader(out_file2);
	  crv1->write(out_file2);
	}

      crv1 = kinks[ki].second->geomCurve();
      sfcv = dynamic_pointer_cast<CurveOnSurface,ParamCurve>(crv1);
      if (sfcv.get())
	{
	  sfcv->spaceCurve()->writeStandardHeader(out_file2);
	  sfcv->spaceCurve()->write(out_file2);
	}
      else
	{
	  crv1->writeStandardHeader(out_file2);
	  crv1->write(out_file2);
	}
  }
}

