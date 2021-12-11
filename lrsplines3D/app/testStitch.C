#include <iostream>
#include <fstream>

#include "GoTools/lrsplines3D/LRVolApprox.h"
#include "GoTools/lrsplines3D/LRVolStitch.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/PointCloud.h"

using namespace std;
using namespace Go;

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

int compare_w_par(const void* el1, const void* el2)
{
  if (((double*)el1)[2] < ((double*)el2)[2])
    return -1;
  else if (((double*)el1)[2] > ((double*)el2)[2])
    return 1;
  else
    return 0;
}

int main (int argc, char *argv[]) {

  if (argc != 9) {
    cout << "usage: ./testMBA <input 4d pt cloud(.g2)> <output lrspline(.g2)> <tolerance> <levels> <nmb split u> <nmb split v> <nmb split w> <cont>" << endl;
    return -1;
  }

  ifstream ifs(argv[1]);
  ofstream ofs(argv[2]);

  double epsge = atof(argv[3]);
  int levels = atoi(argv[4]);
  int split_u = atoi(argv[5]);
  int split_v = atoi(argv[6]);
  int split_w = atoi(argv[7]);
  int cont = atoi(argv[8]);
  double eps = 0.001;
  
  int num_pts;
  ifs >> num_pts;

  vector<double> pc4d;

  double domain[6];
  domain[0] = domain[2] = domain[4] = 1.0e8;
  domain[1] = domain[3] = domain[5] = -1.0e8;
  double mba_level = 0.0; //0.5*(minval+maxval); //0.0;
  double minval = std::numeric_limits<double>::max();
  double maxval = std::numeric_limits<double>::lowest();
  for (int ix=0; ix!=num_pts; ++ix)
    {
      double p0, p1, p2, q0;
      ifs >> p0 >> p1 >> p2 >> q0;
      pc4d.push_back(p0);
      pc4d.push_back(p1);
      pc4d.push_back(p2);
      pc4d.push_back(q0);
      
      domain[0] = std::min(domain[0], p0);
      domain[1] = std::max(domain[1], p0);
      domain[2] = std::min(domain[2], p1);
      domain[3] = std::max(domain[3], p1);
      domain[4] = std::min(domain[4], p2);
      domain[5] = std::max(domain[5], p2);

      mba_level += q0/(double)num_pts;
      minval = std::min(minval, q0);
      maxval = std::max(maxval, q0);
    }

  // Create sub point clouds
  double udel = (domain[1]-domain[0])/(double)(split_u+1);
  double vdel = (domain[3]-domain[2])/(double)(split_v+1);
  double wdel = (domain[5]-domain[4])/(double)(split_w+1);
  vector<vector<double> > pc4d_sub((split_u+1)*(split_v+1)*(split_w+1));
  vector<array<double, 6> > sub_domain((split_u+1)*(split_v+1)*(split_w+1));
  
  // Sort the points according to the w-parameter
  qsort(&pc4d[0], num_pts, 4*sizeof(double), compare_w_par);

  // Traverse points
  int pp0, pp1;
  double upar, vpar, wpar;
  int ki, kj, kr;
  for (pp0=0, wpar=domain[4]; pp0<(int)pc4d.size() && pc4d[pp0+2]<wpar; pp0+=4);
  for (kr=0, wpar+=wdel; kr<=split_w; ++kr, wpar+=wdel)
    {
      for (pp1=pp0; pp1<(int)pc4d.size() && pc4d[pp1+2]<=wpar; pp1+=4);

      // Sort the current sub set of points according to the v-parameter
      qsort(&pc4d[0]+pp0, (pp1-pp0)/4, 4*sizeof(double), compare_v_par);

      // Traverse sub set of points
      int pp2, pp3;
      for (pp2=pp0, vpar=domain[2]; pp2<pp1 && pc4d[pp2+1]<vpar; pp2+=4);
      for (kj=0, vpar+=vdel; kj<=split_v; ++kj, vpar+=vdel)
	{
	  for (pp3=pp2; pp3<pp1 && pc4d[pp3+1]<=vpar; pp3+=4);

	  // Sort the current sub set of points according to the u-parameter
	  qsort(&pc4d[0]+pp2, (pp3-pp2)/4, 4*sizeof(double), compare_u_par);

	  // Traverse sub set of points
	  int pp4, pp5;
	  for (pp4=pp2, upar=domain[0]; pp4<pp3 && pc4d[pp4]<upar; pp4+=4);
	  for (ki=0, upar+=udel; ki<=split_u; ++ki, upar+=udel)
	    {
	      for (pp5=pp4; pp5<pp3 && pc4d[pp5]<=upar; pp5+=4);

	      int ix = (kr*(split_v+1) + kj)*(split_u+1) + ki;
	      pc4d_sub[ix].insert(pc4d_sub[ix].begin(), pc4d.begin()+pp4, pc4d.begin()+pp5);
	      sub_domain[ix][0] = upar-udel;
	      sub_domain[ix][1] = upar;
	      sub_domain[ix][2] = vpar-vdel;
	      sub_domain[ix][3] = vpar;
	      sub_domain[ix][4] = wpar-wdel;
	      sub_domain[ix][5] = wpar;
	      
	      for (; pp5>pp4 && pc4d[pp5]==upar; pp5-=4);
	      pp4 = pp5;
	    }
	  qsort(&pc4d[0]+pp2, (pp3-pp2)/4, 4*sizeof(double), compare_v_par);
	  for (; pp3>pp2 && pc4d[pp3+1]==vpar; pp3-=4);
	  pp2 = pp3;
	}
      qsort(&pc4d[0]+pp0, (pp1-pp0)/4, 4*sizeof(double), compare_w_par);
      for (; pp1>pp0 && pc4d[pp1+2]==wpar; pp1-=4);
      pp0 = pp1;
    }

  // Check
  std::ofstream of0("subpc4d.g2");
  for (size_t kh=0; kh<pc4d_sub.size(); ++kh)
    {
      vector<double> pts(pc4d_sub[kh].size()*3/4);
      int kb = 0;
      for (size_t kw=0; kw<pc4d_sub[kh].size(); kw+=4)
	for (int ka=0; ka<3; ++ka)
	  pts[kb++] = pc4d_sub[kh][kw+ka];

      PointCloud3D cloud(pts.begin(), pts.size()/3);
      cloud.writeStandardHeader(of0);
      cloud.write(of0);
    }
  
  std::cout << "Domain: [" << domain[0] << "," << domain[1] << "]x[" << domain[2];
  std::cout << "," << domain[3] << "]x[" << domain[4] << "," << domain[5] << "]" << std::endl;
  std::cout << "Range: [" << minval << "," << maxval << "]" << std::endl;
  int dim = 1;
  int ncoef = 6;
  int order = 3;
  vector<shared_ptr<LRSplineVolume> > vols;
  for (size_t kh=0; kh<pc4d_sub.size(); ++kh)
    {
      LRVolApprox vol_approx(ncoef, order, ncoef, order, ncoef, order,
			     pc4d_sub[kh], dim, sub_domain[kh].data(), epsge, mba_level);
  
      double max, average, av_all;
      int num_out;
      cout << "Starting approximation..." << endl;

      shared_ptr<LRSplineVolume> result = vol_approx.getApproxVol(max,av_all,average,num_out,levels);
      vols.push_back(result);
      
      cout << "kh = " << kh << "\n"
	   << "max = " << max << " "
	   << "av_all = " << av_all << " "
	   << "average = " << average << " "
	   << "num_out ="  << num_out << " "
	   << "levels = " << levels << " "
	   << "num_elements = " << result->numElements() << " "
	   << endl;
   }

  std::ofstream of1("prestitch.g2");
  for (size_t ki=0; ki<vols.size(); ++ ki)
    {
      vols[ki]->writeStandardHeader(of1);
      vols[ki]->write(of1);
    }
  
  LRVolStitch stitch;
  stitch.stitchRegVols(vols, split_u+1, split_v+1, split_w+1, eps, cont);
  
  for (size_t ki=0; ki<vols.size(); ++ ki)
    {
      vols[ki]->writeStandardHeader(ofs);
      vols[ki]->write(ofs);
    }
}
