#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/igeslib/IGESconverter.h"
#include <fstream>

using std::istream;
using std::vector;
using namespace Go;

ftMessage readg2(istream& is, vector<shared_ptr<ftSurface> >& faces)
//===========================================================================
{

  ftMessage status;
  if (is.bad())
    {
      status.setError(FT_BAD_INPUT_FILE);
      return status;
    }

  IGESconverter conv;
  try
    {
      conv.readgo(is);
    }
  catch (...)
    {
      status.setError(FT_ERROR_IN_READ_IGES);
      return status;
    }

  //      std::ofstream outfile("debug.out");
  //      conv.writedisp(outfile);

  vector<shared_ptr<GeomObject> > gogeom = conv.getGoGeom();

  int nmbgeom = (int)gogeom.size();
  faces.reserve(nmbgeom); // May be too much, but not really important
  int face_count = 0;


  // Remaining geometry.

  for (int i=0; i<nmbgeom; i++)
    {

      if (gogeom[i]->instanceType() == Class_SplineCurve)
	{
	  if (conv.getGroup().size() == 0)
	    {
	      // One mesh surface expected.

	    }
	  else
	    // Not expected. Ignore the current curve.
	    status.addWarning(FT_UNEXPECTED_INPUT_OBJECT_IGNORED);  
	}
      else if (gogeom[i]->instanceType() == Class_SplineSurface ||
	       gogeom[i]->instanceType() == Class_BoundedSurface)
	{
	  shared_ptr<GeomObject> lg = gogeom[i];
	  shared_ptr<ParamSurface> gosf =
	    dynamic_pointer_cast<ParamSurface, GeomObject>(lg);
	  shared_ptr<ftSurface> ftsf(new ftSurface(gosf, face_count++));
	  faces.push_back(ftsf);

	}
    }

  return status;

}



int main( int argc, char* argv[] )
{
  if (argc != 2) {
    std::cout << "Input parameters : Input file on g2 format," << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

    double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double approxtol = 0.01;

  vector<shared_ptr<ftSurface> > faces;
  ftMessage status = readg2(file1, faces);
  std::cout << "Read g2 file. Status message : " << status.getMessage() << std::endl;

  SurfaceModel model1(approxtol, gap, neighbour, kink, 10.0*kink, faces);
  int nmb_bd1 = model1.nmbBoundaries();
  std::cout << "No of boundaries (1): " << nmb_bd1 << std::endl;

  SurfaceModel model2(approxtol, gap, neighbour, kink, 10.0*kink, faces, true);
  int nmb_bd2 = model2.nmbBoundaries();
  std::cout << "No of boundaries (2): " << nmb_bd2 << std::endl;

  exit(0);
}



		      

