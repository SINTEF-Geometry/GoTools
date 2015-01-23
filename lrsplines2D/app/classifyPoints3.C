#include "GoTools/utils/config.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/utils/Array.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSplineUtils.h"
#include "GoTools/lrsplines2D/LRBSpline2D.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;
using std::vector;



int colors[3][3] = {
  {0, 255, 0},
  {255, 255, 255},
  {255, 0, 0},
};

int compare_u_par(const void* el1, const void* el2)
{
  if (((double*)el1)[0] < ((double*)el2)[0])
    return -1;
  else if (((double*)el1)[0] > ((double*)el2)[0])
    return 1;
  else
    return 0;
}

int compare_v_par(const void* el1, const void* el2)
{
  if (((double*)el1)[1] < ((double*)el2)[1])
    return -1;
  else if (((double*)el1)[1] > ((double*)el2)[1])
    return 1;
  else
    return 0;
}


int main(int argc, char *argv[])
{
  if (argc != 6) {
    std::cout << "Usage: surface in (.g2), point cloud (.g2), points_out.g2, max level, nmb _levels" << std::endl;
    return -1;
  }

  std::ifstream sfin(argv[1]);
  std::ifstream ptsin(argv[2]);
  std::ofstream fileout(argv[3]); 
  
  double max_level = atof(argv[4]);
  int nmb_level = atoi(argv[5]);
  double min_level = -max_level;

  ObjectHeader header1;
  header1.read(sfin);
  shared_ptr<LRSplineSurface> sf1(new LRSplineSurface());
  sf1->read(sfin);

  ObjectHeader header2;
  header2.read(ptsin);
  PointCloud3D points;
  points.read(ptsin);

  int nmb_pts = points.numPoints();
  vector<double> data(points.rawData(), points.rawData()+3*nmb_pts);

  int dim = sf1->dimension();
  //RectDomain rd = sf1->containingDomain();

  int ki, kj, kr, ka;
  double *curr;
  double dist;
  vector<double> limits(2*nmb_level+1);
  vector<vector<double> > level_points(2*nmb_level+2);

  // Set distance levels 
  double del = max_level/(double)nmb_level;
  limits[nmb_level] = 0;
  for (ki=1; ki<=nmb_level; ++ki)
    {
      limits[nmb_level-ki] = -ki*del;
      limits[nmb_level+ki] = ki*del;
    }

  // Get all knot values in the u-direction
  const double* const uknots = sf1->mesh().knotsBegin(XFIXED);
  const double* const uknots_end = sf1->mesh().knotsEnd(XFIXED);
  int nmb_knots_u = sf1->mesh().numDistinctKnots(XFIXED);
  const double* knotu;

  // Get all knot values in the v-direction
  const double* const vknots = sf1->mesh().knotsBegin(YFIXED);
  const double* const vknots_end = sf1->mesh().knotsEnd(YFIXED);
  int nmb_knots_v = sf1->mesh().numDistinctKnots(YFIXED);
  const double* knotv;

  // Construct mesh of element pointers
  vector<Element2D*> elements;
  sf1->constructElementMesh(elements);

  double maxdist = 0.0;
  double mindist = 0.0;
  double avdist = 0.0;

  // For each point, classify according to distance
  // Sort points in v-direction
  qsort(&data[0], nmb_pts, 3*sizeof(double), compare_v_par);

  int pp0, pp1;
  Element2D* elem = NULL;
  for (pp0=0, knotv=vknots; pp0<(int)data.size() && data[pp0+1] < (*knotv); 
       pp0+=3);
  for (kj=0, ++knotv; knotv!= vknots_end; ++knotv, ++kj)
    {
      
      for (pp1=pp0; pp1<(int)data.size() && data[pp1+1] < (*knotv); pp1+=3);
      if (knotv+1 == vknots_end)
	for (; pp1<(int)data.size() && data[pp1+1] <= (*knotv); pp1+=3);
      // 	pp1 = (int)data.size();

      // Sort the current sub set of data according to the u-parameter
      qsort(&data[0]+pp0, (pp1-pp0)/3, 3*sizeof(double), compare_u_par);

      // Traverse the relevant data and identify the associated element
      int pp2, pp3;
      for (pp2=pp0, knotu=uknots; pp2<pp1 && data[pp2] < (*knotu); pp2+=3);
      for (ki=0, ++knotu; knotu!=uknots_end; ++knotu, ++ki)
	{
	  for (pp3=pp2; pp3<pp1 && data[pp3] < (*knotu); pp3 += 3);
	  if (knotu+1 == uknots_end)
	    for (; pp3<pp1 && data[pp3] <= (*knotu); pp3+=3);
	  //   pp3 = pp1;
	  
	  // Fetch associated element
	  elem = elements[kj*(nmb_knots_u-1)+ki];

	  int nump = (pp3 - pp2)/3;
	  for (kr=0, curr=&data[pp2]; kr<nump; ++kr, curr+=3)
	    {
	      // Evaluate
	      Point pos;
	      sf1->point(pos, curr[0], curr[1], elem);
	      dist = curr[2]-pos[0];

	      maxdist = std::max(maxdist, dist);
	      mindist = std::min(mindist, dist);
	      avdist += fabs(dist);
	  
	      // Find classification
	      for (ka=0; ka<limits.size(); ++ka)
		if (dist < limits[ka])
		  {
		    level_points[ka].push_back(curr[0]);
		    level_points[ka].push_back(curr[1]);
		    level_points[ka].push_back(curr[2]);
		    break;
		  }
	      if (ka == limits.size())
		{
		  level_points[ka].push_back(curr[0]);
		  level_points[ka].push_back(curr[1]);
		  level_points[ka].push_back(curr[2]);
		}
	    }
	  pp2 = pp3;
	}
      pp0 = pp1;
    }

  // Write to file
  for (ki=0; ki<level_points.size(); ++ki)
    {
      if (level_points[ki].size() == 0)
	continue;

      // Make point cloud
      PointCloud3D level_cloud(level_points[ki].begin(), level_points[ki].size()/3);

      double cc[3];
      if (ki <= nmb_level)
	{
	  cc[0] = ((nmb_level-ki)*colors[0][0] + ki*colors[1][0])/nmb_level;
	  cc[1] = ((nmb_level-ki)*colors[0][1] + ki*colors[1][1])/nmb_level;
	  cc[2] = ((nmb_level-ki)*colors[0][2] + ki*colors[1][2])/nmb_level;
	}
      else
	{
	  cc[0] = ((ki-nmb_level-1)*colors[2][0] + 
		   (2*nmb_level-ki+1)*colors[1][0])/nmb_level;
	  cc[1] = ((ki-nmb_level-1)*colors[2][1] + 
		   (2*nmb_level-ki+1)*colors[1][1])/nmb_level;
	  cc[2] = ((ki-nmb_level-1)*colors[2][2] + 
		   (2*nmb_level-ki+1)*colors[1][2])/nmb_level;
	}

      fileout << "400 1 0 4 " << cc[0] << " " << cc[1];
      fileout << " " << cc[2] << " 255" << std::endl;
      level_cloud.write(fileout);
    }

  avdist /= (double)nmb_pts;

  std::cout << "Max dist: " << maxdist << "Max dist below: " << mindist;
  std::cout << ", average dist: " << avdist << std::endl;
}

      




 
