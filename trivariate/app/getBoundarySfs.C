#include <fstream>
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/Utils.h"
#include <memory>

using namespace Go;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::vector;


int main(int argc, char* argv[] )
{
  if (argc != 3)
      cout << "Usage: " << "infile outfile " << endl;

  // Open input surface file
  ifstream is(argv[1]);
  ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

  // Open outfile
  ofstream os(argv[2]);
  ALWAYS_ERROR_IF(os.bad(), "Bad or no output filename");

  while (!is.eof())
    {
      ObjectHeader head;
      is >> head;

      // Read volume from file
      SplineVolume vol;
      is >> vol;



      vector<shared_ptr<SplineSurface> > bd_sfs = vol.getBoundarySurfaces();
      vector<shared_ptr<SplineSurface> > bd_sfs2 = vol.getBoundarySurfaces();
      for (size_t ki=0; ki<bd_sfs2.size(); ++ki)
	{
	  bd_sfs2[ki]->writeStandardHeader(os);
	  bd_sfs2[ki]->write(os);
	}

      Utils::eatwhite(is);
    }

}
