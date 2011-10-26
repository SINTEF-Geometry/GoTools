#include <fstream>
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/SweepVolumeCreator.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/trivariatemodel/ftVolumeTools.h"
#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/compositemodel/ftSurface.h"

using namespace Go;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::shared_ptr;

int main(int argc, char* argv[] )
{
  if (argc != 13)
      cout << "Usage: " << "<infile1> <point> <axis> <angle> <bd idx> <bd idx2> <length> <outfile>" << endl;

  ifstream infile1(argv[1]);
  ALWAYS_ERROR_IF(infile1.bad(), "Bad or no input filename");

  ObjectHeader header;
  header.read(infile1);

  shared_ptr<SplineVolume> vol1(new SplineVolume());
  vol1->read(infile1);

  Point pos(atof(argv[2]),atof(argv[3]),atof(argv[4]));
  Point axis(atof(argv[5]),atof(argv[6]),atof(argv[7]));
  double angle = atof(argv[8]);
  int bd = atoi(argv[9]);
  int bd2 = atoi(argv[10]);
  double len = atof(argv[11]);

  ofstream outfile(argv[12]);

  shared_ptr<SplineSurface> bd_sf = vol1->getBoundarySurface(bd);

  SweepVolumeCreator creator;
  shared_ptr<SplineVolume> vol2 = 
    shared_ptr<SplineVolume>(creator.rotationalSweptVolume(*bd_sf,
							   angle,
							   pos, 
							   axis));

  shared_ptr<SplineSurface> bd_sf2 = vol2->getBoundarySurface(bd2);
  
  double u1 = 0.5*(bd_sf2->startparam_u() + bd_sf2->endparam_u());
  double v1 = 0.5*(bd_sf2->startparam_v() + bd_sf2->endparam_v());
  vector<Point> res(3);
  bd_sf2->point(res, u1, v1, 1);
  Point norm = res[1].cross(res[2]);
  norm.normalize();
  Point pnt2 = res[0] + len*norm;
  shared_ptr<SplineCurve> crv
    = shared_ptr<SplineCurve>(new SplineCurve(res[0], pnt2));
  shared_ptr<SplineVolume> vol3 = 
    shared_ptr<SplineVolume>(creator.linearSweptVolume(*bd_sf2,
						       *crv,
						       res[0]));
  


  vol1->writeStandardHeader(outfile);  
  vol1->write(outfile);
  vol2->writeStandardHeader(outfile);  
  vol2->write(outfile);
  vol3->writeStandardHeader(outfile);  
  vol3->write(outfile);
}
