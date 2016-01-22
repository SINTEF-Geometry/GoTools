#include "GoTools/utils/config.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/utils/Array.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;
using std::vector;



int main(int argc, char *argv[])
{
  if (argc != 4 && argc != 5) {
    std::cout << "Usage: surface in (.g2) point cloud (.g2) lrspline_out.g2 (grid (0/1))" << std::endl;
    return -1;
  }

  std::ifstream sfin(argv[1]);
  std::ifstream ptsin(argv[2]);
  std::ofstream fileout(argv[3]); 

  bool grid = false;
  if (argc == 5)
    grid = atoi(argv[4]);
  
  ObjectHeader header1;
  header1.read(sfin);
  shared_ptr<LRSplineSurface> sf1(new LRSplineSurface());
  sf1->read(sfin);

  ObjectHeader header2;
  header2.read(ptsin);
  PointCloud3D points;
  points.read(ptsin);

  int ki;

  double limit[2];
  double cell_del[2];
  if (grid)
    {
      std::cout << "Give domain start (umin, umax): " << std::endl;
      for (ki=0; ki<2; ++ki)
	std::cin >> limit[ki];
      std::cout << "Cell size (u, v): " << std::endl;
      for (ki=0; ki<2; ++ki)
	std::cin >> cell_del[ki];
    }
      

  int dim = sf1->dimension();
  RectDomain rd = sf1->containingDomain();
  int maxiter = 4;
  double aeps = 0.001;

  double umin = sf1->paramMin(XFIXED);
  double umax = sf1->paramMax(XFIXED);
  double vmin = sf1->paramMin(YFIXED);
  double vmax = sf1->paramMax(YFIXED);

  int nmb_pts = points.numPoints();
  vector<double> data(points.rawData(), points.rawData()+3*nmb_pts);

  double *curr;
  double dist;
  double maxdist = 0.0;
  double mindist = 0.0;
  double avdist = 0.0;
  double maxdist2 = 0.0;
  double mindist2 = 0.0;
  double avdist2 = 0.0;
  for (ki=0, curr=&data[0]; ki<nmb_pts; ++ki, curr+=3)
    {
      if (dim == 3)
	{
	  // Perform closest point
	  double upar, vpar;
	  Point close_pt;
	  Point curr_pt(curr, curr+dim);
	  sf1->closestPoint(curr_pt, upar, vpar, close_pt,
			    dist, aeps, maxiter, NULL, &rd, curr);
	  Point vec = curr_pt - close_pt;
	  Point norm;
	  sf1->normal(norm, upar, vpar);
	  if (vec*norm < 0.0)
	    dist *= -1;
	}
      else
	{
	  // Evaluate
	  Point pos;
	  sf1->point(pos, curr[0], curr[1]);
	  dist = curr[2]-pos[0];
	}

      maxdist = std::max(maxdist, dist);
      mindist = std::min(mindist, dist);
      avdist += fabs(dist);

      fileout << curr[0] << "  " << curr[1] << "  " << curr[2];
      fileout << "  " << dist;

      if (grid && dim == 1)
	{
	  // Compute cell distance
	  // First identify cell
	  int idx1 = (curr[0] - limit[0])/cell_del[0];
	  int idx2 = (curr[1] - limit[1])/cell_del[1];

	  // Evaluate corners
	  Point pos1, pos2, pos3, pos4;
	  double u1 = std::max(umin, limit[0]+idx1*cell_del[0]);
	  double u2 = std::min(umax, limit[0]+(idx1+1)*cell_del[0]);
	  double v1 = std::max(vmin, limit[1]+idx2*cell_del[1]);
	  double v2 = std::min(vmax, limit[1]+(idx2+1)*cell_del[1]);
	  sf1->point(pos1, u1, v1);
	  sf1->point(pos2, u2, v1);
	  sf1->point(pos3, u1, v2);
	  sf1->point(pos4, u2, v2);
	  double dist1 = curr[2]-pos1[0];
	  double dist2 = curr[2]-pos2[0];
	  double dist3 = curr[2]-pos3[0];
	  double dist4 = curr[2]-pos4[0];

	  if (dist1*dist2<0.0 || dist1*dist3<0.0 || dist1*dist4<0.0 || 
	      dist2*dist3<0.0 || dist2*dist4<0.0 || dist3*dist4<0.0)
	    dist = 0.0;
	  else
	    {
	      int sgn = (dist1 >= 0.0) ? 1 : -1;
	      dist = std::min(std::min(fabs(dist1),fabs(dist2)),
			      std::min(fabs(dist3),fabs(dist4)));
	      dist *= sgn;
	    }
	  maxdist2 = std::max(maxdist2, dist);
	  mindist2 = std::min(mindist2, dist);
	  fileout << " " << dist;
	  avdist2 += fabs(dist);
	}
     fileout << std::endl;

      //dist = fabs(dist);
     }

  avdist /= (double)nmb_pts;

  std::cout << "Max dist: " << maxdist << "Max dist below: " << mindist;
  std::cout << ", average dist: " << avdist << std::endl;

  if (grid)
    {
      avdist2 /= (double)nmb_pts;
      std::cout << "Max grid dist: " << maxdist2 << "Max dist below: " << mindist2;
      std::cout << ", average dist: " << avdist2 << std::endl;
    }
}


