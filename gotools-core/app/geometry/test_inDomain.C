#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/RectDomain.h"
#include <fstream>

using namespace Go;
using std::vector;

int main(int argc, char** argv)
{
  if (argc != 2 && argc != 4)
    {
      std::cout << "Usage; infile (nmb u pts) (nmb v pts)" << std::endl;
      return -1;
    }

  std::ifstream infile(argv[1]);
  ObjectHeader head;
  head.read(infile);
  BoundedSurface gosf;
  gosf.read(infile);

  const RectDomain& dom = gosf.containingDomain();
  double start_u = dom.umin();
  double end_u = dom.umax();
  double start_v = dom.vmin();
  double end_v = dom.vmax();

  int nmb_u = 20, nmb_v = 20;
  if (argc == 4)
    {
      nmb_u = atoi(argv[2]);
      nmb_v = atoi(argv[3]);
    }
  double del_u = (end_u - start_u)/(double)(nmb_u-1);
  double del_v = (end_v - start_v)/(double)(nmb_v-1);

  vector<Point> inside;
  int ki, kj;
  double u, v;
  for (v=start_v, ki=0; ki<nmb_v; v+=del_v, ki++)
    for (u=start_u, kj=0; kj<nmb_u; u+=del_u, kj++)
      {
	Point curr;

	if (!gosf.inDomain(u,v))
	  continue;

	gosf.point(curr, u, v);
	inside.push_back(curr);
      }

  std::ofstream out_file("inside_pts.g2");
  out_file << "400 1 0 4 255 0 0 255" << std::endl;
  out_file << inside.size() << std::endl;
  for (ki=0; ki < (int)inside.size(); ki++)
    {
      out_file << inside[ki][0] << " " << inside[ki][1] << " ";
      out_file << inside[ki][2] << std::endl;
    }
}


