#include <fstream>
#include <iomanip>
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/intersections/ParamSurfaceInt.h"
#include "GoTools/intersections/SplineSurfaceInt.h"
#include "GoTools/intersections/SfSfIntersector.h"
#include "GoTools/geometry/SISLconversion.h"
#include "sislP.h"


using namespace std;
using namespace Go;


int main(int argc, char** argv)
{
 


  ObjectHeader header;

  // Read the surface from file
  ifstream input1(argv[1]);
  if (input1.bad()) {
    return 1;
  }
  
  header.read(input1);
  shared_ptr<SplineSurface> surf(new SplineSurface());
  surf->read(input1);
  input1.close();

  DirectionCone cone = surf->normalCone();
  printf("Greater than pi : %d \n",cone.greaterThanPi());
  printf("Angle : %13.7f \n",cone.angle());

  int kstat = 0;
  SISLSurf *sf = GoSurf2SISL(*surf.get());
  s1990(sf, 0.000001, &kstat);
  printf("Sisl igtpi : %d, angle %13.7f \n",sf->pdir->igtpi,sf->pdir->aang);
  
  freeSurf(sf);
}
