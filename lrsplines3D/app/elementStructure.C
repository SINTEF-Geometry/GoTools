#include <iostream>
#include <fstream>

#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/lrsplines3D/Element3D.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/LineCloud.h"
#include "GoTools/geometry/GoTools.h"

using namespace std;
using namespace Go;

int colors[3][3] = {
  {0, 255, 0},
  {255, 255, 255},
  {255, 0, 0},
};


int main (int argc, char *argv[]) {

  if (argc != 4) {
    cout << "usage: ./elementStructure <input lrspline(.g2)> <output element mid points> <output element boundaries>" << endl;
    return -1;
  }

  ifstream ifs(argv[1]);
  ofstream ofs(argv[2]);
  ofstream ofs2(argv[3]);

  GoTools::init();

  ObjectHeader oh;
  oh.read(ifs);

  shared_ptr<LRSplineVolume> vol(new LRSplineVolume());
  vol->read(ifs);

  int num_el = vol->numElements();
  vector<double> mid_el(3*num_el);
  vector<double> bd_el(72*num_el);
  int ki=0, kj=0;
  Point pos, pos2, dir1, dir2, dir3;
  for (auto it=vol->elementsBegin(); it != vol->elementsEnd(); ++it)
    {
      const Element3D* elem = it->second.get();
      double umin = elem->umin();
      double umax = elem->umax();
      double vmin = elem->vmin();
      double vmax = elem->vmax();
      double wmin = elem->wmin();
      double wmax = elem->wmax();
      mid_el[ki++] = 0.5*(umin+umax);
      mid_el[ki++] = 0.5*(vmin+vmax);
      mid_el[ki++] = 0.5*(wmin+wmax);

      pos = Point(umin, vmin, wmin);
      dir1 = Point(umax-umin, 0.0, 0.0);
      dir2 = Point(0.0, vmax-vmin, 0.0);
      dir3 = Point(0.0, 0.0, wmax-wmin);

      pos2 = pos+dir1;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos2 = pos+dir2;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos2 = pos+dir3;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos += dir1;
      pos2 = pos+dir2;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos2 = pos+dir3;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos += dir2;
      pos2 = pos-dir1;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos2 = pos+dir3;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos += dir3;
      pos2 = pos-dir1;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos2 = pos-dir2;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos -= dir1;
      pos2 = pos-dir2;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos2 = pos-dir3;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos -= dir2;
      pos2 = pos+dir1;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];
     }

  // Make point cloud
  PointCloud3D elem_cloud(mid_el.begin(), mid_el.size()/3);

  elem_cloud.writeStandardHeader(ofs);
  elem_cloud.write(ofs);

  LineCloud elem_lines(bd_el.begin(), bd_el.size()/6);
  elem_lines.writeStandardHeader(ofs2);
  elem_lines.write(ofs2);
}
