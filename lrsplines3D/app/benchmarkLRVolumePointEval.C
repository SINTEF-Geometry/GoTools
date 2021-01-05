//===========================================================================
//                                                                           
// File: benchmarkLRVolumePointEval.C                                       
//                                                                           
// Created: Tue Jan 29 13:57:16 2013                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description: Benchmark the point evaluation (with higher derivs)
//              for a LRSplineVolume.  If input is a regular
//              SplineVolume, we convert to LRSplineVolume and
//              compare. If input is a LRSplineVolume, we compare
//              with the corresponding (refined) SplineVolume.
//                                                                           
//===========================================================================


#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/lrsplines3D/Mesh3D.h"
#include "GoTools/lrsplines3D/Direction3D.h"
//#include "GoTools/lrsplines3D/LRBenchmarkUtils3D.h"
#include "GoTools/lrsplines3D/LRSpline3DUtils.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/ClassType.h"
#include "GoTools/utils/config.h"
#include "GoTools/geometry/PointCloud.h"


#include <iostream>
#include <assert.h>
#include <fstream>
#include <string>
#include "string.h"


using Go::ObjectHeader;
using Go::SplineVolume;
using Go::LRSplineVolume;
using std::vector;
using namespace Go;


// Assuming the domain is the same.
double maxDist(const SplineVolume& spline_vol,
	       const LRSplineVolume& lr_spline_vol,
	       int num_samples_u,
	       int num_samples_v,
	       int num_samples_w,
	       int der_u, int der_v, int der_w,
	       vector<double>& sampled_pts_lr);
#if 0
double maxDistNormals(const SplineVolume& spline_vol,
		      const LRSplineVolume& lr_spline_vol,
		      int num_samples_u,
		      int num_samples_v,
		      int num_samples_w);
#endif

int main(int argc, char *argv[])
{
  if (argc != 5)
  {
      std::cout << "Usage: (lr)spline_vol.g2 num_dir_samples sum_derivs num_iter" << std::endl;
      return -1;
  }

  std::cout.precision(15);

  char* filein_char = argv[1];
  std::ifstream filein(argv[1]); // Input lr spline.
  // If input surface is not fine enough we refine.
  int num_dir_samples = atoi(argv[2]);
  if (num_dir_samples < 2)
  {
      puts("num_dir_samples was less than the minimal value of 2, setting it to 2.");
      num_dir_samples = 2;
  }
  int sum_derivs = atoi(argv[3]);
  if (sum_derivs > 0)
  {
      puts("Derivative evaluation not supported yet, setting sum_derivs to 0.");
      sum_derivs = 0;
  }
  int num_iter = atoi(argv[4]);
  if (num_iter != 1)
  {
      puts("Not using num_iter yet.");
  }

  shared_ptr<SplineVolume> spline_vol;
  shared_ptr<LRSplineVolume> lr_spline_vol;
  double knot_tol = 1e-10;

  // If we refine the lr-surface in order to extract a regular
  // spline-version we must perform the operations on a copy.
  shared_ptr<LRSplineVolume> lr_spline_vol_copy;

  int order_u, order_v, num_coefs_u, num_coefs_v, dim, num_bases=-1;
  if (strstr(filein_char, ".g2")) // We do not expect an uppercase suffix (.G2).
    {
      ObjectHeader header;
      filein >> header;
      if (header.classType() == Class_SplineVolume)
	{
	  std::cout << "Input was a SplineVolume, creating a LRSplineVolume." << std::endl;
	  spline_vol = shared_ptr<SplineVolume>(new SplineVolume());
	  filein >> *spline_vol;
	  dim = spline_vol->dimension();
	  // We create the lr-version.
	  lr_spline_vol = shared_ptr<LRSplineVolume>(new LRSplineVolume(spline_vol.get(), knot_tol));
	}
      else if (header.classType() == Class_LRSplineVolume)
	{
	  std::cout << "Input was a LRSplineVolume, creating a SplineVolume (refining)." << std::endl;
	  lr_spline_vol = shared_ptr<LRSplineVolume>(new LRSplineVolume());
	  filein >> *lr_spline_vol;
	  dim = lr_spline_vol->dimension();

	  lr_spline_vol_copy = shared_ptr<LRSplineVolume>(lr_spline_vol->clone());
	  spline_vol = shared_ptr<SplineVolume>(LRSpline3DUtils::fullTensorProductVolume(*lr_spline_vol));
	  int num_coefs_u = spline_vol->numCoefs(0);
	  int num_coefs_v = spline_vol->numCoefs(1);
	  int num_coefs_w = spline_vol->numCoefs(2);
	  std::cout << "num_coefs_u = " << num_coefs_u << ", num_coefs_v = " << num_coefs_v <<
	      "num_coefs_w = " << num_coefs_w << std::endl;
	  std::ofstream spline_out("tmp/lr_spline_vol_version.g2");
	  spline_vol->writeStandardHeader(spline_out);
	  spline_vol->write(spline_out);
	}
      else
	{
	  std::cout << "Input was not a SplineVolume or a LRSplineVolume, exiting!" << std::endl;
	}
    }
  else
    {
      MESSAGE("Input was not a g2-file!");
      return -1;
    }

  // We include various testing of various functions.
  bool test_lr_functions = true;
  if (test_lr_functions)
    {
#if 0
      // We write to file the grid.
      std::ofstream grid_pre("tmp/lr_grid_pre.ps");
      writePostscriptMesh(*lr_spline_vol, grid_pre);
#endif

      bool use_unit_domain = false;
      if (use_unit_domain)
	{ // We rescale to the unit domain.
	    lr_spline_vol->setParameterDomain(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
	  std::ofstream lr_spline_out("tmp/lr_spline_vol_unit.g2");
	  lr_spline_vol->writeStandardHeader(lr_spline_out);
	  lr_spline_vol->write(lr_spline_out);
	}

      bool extract_subsurf = false;
      if (extract_subsurf)
	{ // We extract a subsurface. Set values manually.
	  double umin = 498850.13939999999 - 1;//0.0;
	  double umax = 498858.58250000002 + 1;//1.0;
	  double vmin = 3875136 - 1;//0.0;
	  double vmax = 3875144.2230000002 + 1;//1.0;
	  double wmin = 3875136 - 1;//0.0;
	  double wmax = 3875144.2230000002 + 1;//1.0;
	  shared_ptr<LRSplineVolume> sub_vol(lr_spline_vol->subVolume(umin, vmin, wmin, umax, vmax, wmax, knot_tol));
	  std::ofstream sub_vol_out("tmp/lr_spline_sub_vol.g2");
	  sub_vol->writeStandardHeader(sub_vol_out);
	  sub_vol->write(sub_vol_out);
	}

      bool rev_dir_u = false;
      if (rev_dir_u)
	{
	  lr_spline_vol->reverseParameterDirection(0);
	  spline_vol->reverseParameterDirection(0);
	}
      bool rev_dir_v = false;
      if (rev_dir_v)
	{
	  lr_spline_vol->reverseParameterDirection(1);
	  spline_vol->reverseParameterDirection(1);
	}
      bool rev_dir_w = false;
      if (rev_dir_w)
	{
	  lr_spline_vol->reverseParameterDirection(2);
	  spline_vol->reverseParameterDirection(2);
	}

      bool swap_par_dir = false;
      if (swap_par_dir)
	{
	  MESSAGE("Swapping parameter directions!");
	  lr_spline_vol->swapParameterDirection(0, 1);
	  spline_vol->swapParameterDirection(0, 1);
	  lr_spline_vol->swapParameterDirection(1, 2);
	  spline_vol->swapParameterDirection(1, 2);
	}

      bool test_constParamSurface = false;//true;
      if (test_constParamSurface)
	{
	  double wgt_u = 0.129;
	  double upar = wgt_u*lr_spline_vol->startparam_u() + (1.0 - wgt_u)*lr_spline_vol->endparam_u();
	  shared_ptr<LRSplineSurface> par_sf_u(lr_spline_vol->constParamSurface(upar, 0));
	  std::ofstream sf_out("tmp/const_par_sf.g2");
	  par_sf_u->writeStandardHeader(sf_out);
	  par_sf_u->write(sf_out);
	  double wgt_v = 0.279;
	  double vpar = wgt_v*lr_spline_vol->startparam_v() + (1.0 - wgt_v)*lr_spline_vol->endparam_v();
	  shared_ptr<LRSplineSurface> par_sf_v(lr_spline_vol->constParamSurface(vpar, 1));
	  par_sf_v->writeStandardHeader(sf_out);
	  par_sf_v->write(sf_out);
	  double wgt_w = 0.479;
	  double wpar = wgt_v*lr_spline_vol->startparam_w() + (1.0 - wgt_v)*lr_spline_vol->endparam_w();
	  shared_ptr<LRSplineSurface> par_sf_w(lr_spline_vol->constParamSurface(wpar, 2));
	  par_sf_w->writeStandardHeader(sf_out);
	  par_sf_w->write(sf_out);
	}

#if 0
      // We write to file the grid.
      std::ofstream grid_post("tmp/lr_grid_post.ps");
      writePostscriptMesh(*lr_spline_vol, grid_post);
#endif
    }

  vector<double> sampled_pts_lr;
  for (int kk = 0; kk < sum_derivs + 1; ++kk)
    for (int kj = 0; kj < sum_derivs + 1; ++kj)
      for (int ki = 0; ki < sum_derivs + 1 - kj; ++ki)
	{
	  double max_dist = maxDist(*spline_vol, *lr_spline_vol,
				    num_dir_samples, num_dir_samples, num_dir_samples,
				    ki, kj, kk,
				    sampled_pts_lr);

	  std::cout << "Max dist with der_u=" << ki << " & der_v=" << kj << " & der_w=" << kk << ": " << max_dist << std::endl;
	}

  if (sampled_pts_lr.size() > 0 && dim == 3)
    {
      std::ofstream fileout_pts("tmp/sampled_pts.g2");
      int num_pts = sampled_pts_lr.size()/dim;
      PointCloud3D pt_cl(sampled_pts_lr.begin(), num_pts);
      pt_cl.writeStandardHeader(fileout_pts);
      pt_cl.write(fileout_pts);
    }

#if 0
  double max_dist_normals = maxDistNormals(*spline_vol, *lr_spline_vol, num_dir_samples, num_dir_samples);
  std::cout << "Max dist normals = " << max_dist_normals << std::endl;
#endif

}


double maxDist(const SplineVolume& spline_vol,
	       const LRSplineVolume& lr_spline_vol,
	       int num_samples_u,
	       int num_samples_v,
	       int num_samples_w,
	       int der_u, int der_v, int der_w,
	       vector<double>& sampled_pts_lr)
{
  const int derivs = std::max(der_u, der_v);
  const int sum_derivs = der_u + der_v;
  const int num_pts = (sum_derivs + 1)*(sum_derivs + 2)/2;
  const int first_pos = (sum_derivs + 1)*sum_derivs/2;
  const int der_pos = first_pos + der_v;
  // Assuming the domain is the same.
  double umin = lr_spline_vol.startparam_u();
  double umax = lr_spline_vol.endparam_u();
  double vmin = lr_spline_vol.startparam_v();
  double vmax = lr_spline_vol.endparam_v();
  double wmin = lr_spline_vol.startparam_w();
  double wmax = lr_spline_vol.endparam_w();
  double ustep = (umax - umin)/((double)num_samples_u - 1);
  double vstep = (vmax - vmin)/((double)num_samples_v - 1);
  double wstep = (wmax - wmin)/((double)num_samples_w - 1);
  Point go_pt(3), lr_pt(3);
  double max_dist = -1.0;
// #ifndef NDEBUG
  double max_dist_u = 0.0;
  double max_dist_v = 0.0;
  double max_dist_w = 0.0;
  Point max_go_pt(3), max_lr_pt(3);
  vector<Point> go_pts(num_pts), lr_pts(num_pts);
// #endif
  for (int kk = 0; kk < num_samples_w; ++kk)
    {
      double wpar = wmin + kk*wstep;
      if (wpar > wmax)
	wpar = wmax;
      {
	for (int kj = 0; kj < num_samples_v; ++kj)
	  {
	    double vpar = vmin + kj*vstep;
	    if (vpar > vmax)
	      vpar = vmax;
	    for (int ki = 0; ki < num_samples_u; ++ki)
	      {
		double upar = umin + ki*ustep;
		if (upar > umax)
		  upar = umax;
                spline_vol.point(go_pts, upar, vpar, wpar, sum_derivs);
		go_pt = go_pts[der_pos];
//	  lr_spline_vol.point(lr_pts, upar, vpar, sum_derivs);
		Point lr_pt2 = lr_spline_vol(upar, vpar, wpar, der_u, der_v, der_w);
//lr_pts[der_pos];
		sampled_pts_lr.insert(sampled_pts_lr.end(), lr_pt2.begin(), lr_pt2.end());
		double dist = go_pt.dist(lr_pt2);
		if (dist > max_dist)
		  {
		    max_dist = dist;
// #ifndef NDEBUG
		    max_dist_u = upar;
		    max_dist_v = vpar;
		    max_dist_w = wpar;
		    max_go_pt = go_pt;
		    max_lr_pt = lr_pt2;
// #endif
		  }
	      }
	  }    
      }
    }
// #ifndef NDEBUG
  std::cout << "max_dist: u = " << max_dist_u << ", v = " << max_dist_v << ", w = " << max_dist_w << std::endl;
  std::cout << "max_go_pt = " << max_go_pt << std::endl;
  std::cout << "max_lr_pt = " << max_lr_pt << std::endl;
// #endif

  return max_dist;
}


#if 0
double maxDistNormals(const SplineVolume& spline_vol,
		      const LRSplineVolume& lr_spline_vol,
		      int num_samples_u,
		      int num_samples_v)
{
  // Assuming the domain is the same.
  double umin = lr_spline_vol.startparam_u();
  double umax = lr_spline_vol.endparam_u();
  double vmin = lr_spline_vol.startparam_v();
  double vmax = lr_spline_vol.endparam_v();
  double ustep = (umax - umin)/((double)num_samples_u - 1);
  double vstep = (vmax - vmin)/((double)num_samples_v - 1);
  Point go_normal(3), lr_normal(3);
  double max_dist_normals = -1.0;
// #ifndef NDEBUG
  double max_dist_u = 0.0;
  double max_dist_v = 0.0;
  Point max_go_normal(3), max_lr_normal(3);
// #endif
  for (int kj = 0; kj < num_samples_v; ++kj)
    {
      double vpar = vmin + kj*vstep;
      if (vpar > vmax)
	vpar = vmax;
      for (int ki = 0; ki < num_samples_u; ++ki)
	{
	  double upar = umin + ki*ustep;
	  if (upar > umax)
	    upar = umax;
	  spline_vol.normal(go_normal, upar, vpar);
//	  lr_spline_vol.point(lr_pts, upar, vpar, sum_derivs);
	  lr_spline_vol.normal(lr_normal, upar, vpar);
//lr_pts[der_pos];
	  double dist_normal = go_normal.dist(lr_normal);
	  if (dist_normal > max_dist_normals)
	    {
	      max_dist_normals = dist_normal;
// #ifndef NDEBUG
	      max_dist_u = upar;
	      max_dist_v = vpar;
	      max_go_normal = go_normal;
	      max_lr_normal = lr_normal;
// #endif
	    }
	}
    }    

// #ifndef NDEBUG
  std::cout << "max_dist_normal: u = " << max_dist_u << ", v = " << max_dist_v << std::endl;
  std::cout << "max_go_normal = " << max_go_normal << std::endl;
  std::cout << "max_lr_normal = " << max_lr_normal << std::endl;
// #endif

  return max_dist_normals;
}
#endif
