//===========================================================================
//                                                                           
// File: testLRsplineRefinement.cpp                                          
//                                                                           
// Created: Tue Nov 20 10:55:54 2012                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description: Compare result after refining different versions of lrsplines.
//                                                                           
//===========================================================================


//#include "GoTools/lr_splines/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSplinePlotUtils.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineSurface.h"


#include <iostream>
#include <assert.h>
#include <fstream>

// Assuming the domain is the same.
double maxDist(const Go::SplineSurface& spline_sf,
	       const Go::LRSplineSurface& lr_spline_sf,
	       int nmb_samples_u, int nmb_samples_v);


int main(int argc, char *argv[])
{
  if (argc != 3)
  {
      std::cout << "Usage: spline_sf.g2 ref_lr_spline_sf.g2" << std::endl;
      return -1;
  }

  std::ifstream filein(argv[1]); // Input (regular) spline surface.
  std::ofstream fileout(argv[2]); // Input (regular) spline surface.

  Go::ObjectHeader header;
  filein >> header;
  Go::SplineSurface spline_sf;
  if (header.classType() == Go::Class_SplineSurface)
  {
      filein >> spline_sf;
  }
  else
  {
      puts("Unexpected geometry input.");
      return -1;
  }

  puts("Now we convert from SplineSurface to LRSplineSurface!");
  double knot_tol = 1e-05;
  shared_ptr<Go::LRSplineSurface> lr_spline_sf(new Go::LRSplineSurface(spline_sf.clone(), knot_tol));
  puts("Done converting!");

  // We test to see if conversion was correct.
  int nmb_samples_u = 101; // Rather random.
  int nmb_samples_v = 36;
  double max_dist = maxDist(spline_sf, *lr_spline_sf, nmb_samples_u, nmb_samples_v);
  std::cout << "Max dist between input and converted surface: " << max_dist << std::endl;

  // We write to screen the number of element and basis functions.
  int num_basis_funcs = lr_spline_sf->numBasisFunctions();
  int num_elem = lr_spline_sf->numElements();
  std::cout << "num_basis_funcs: " << num_basis_funcs << ", num_elem: " << num_elem << std::endl;


  // We hardcode a refinement.
  // Hardcoded values for debugging.
  int dir = 0;
  // @@sbr Does it make sense that paramMin(XFIXED) yields xmin? And XFIXED <=> dir=0: Why?
  double umin = lr_spline_sf->paramMin(Go::XFIXED);
  double umax = lr_spline_sf->paramMax(Go::XFIXED);
  double vmin = lr_spline_sf->paramMin(Go::YFIXED);
  double vmax = lr_spline_sf->paramMax(Go::YFIXED);
  std::cout << "umin: " << umin << std::endl;
  std::cout << "umax: " << umax << std::endl;
  std::cout << "vmin: " << vmin << std::endl;
  std::cout << "vmax: " << vmax << std::endl;

  int num_coefs_u = spline_sf.numCoefs_u();
  int num_coefs_v = spline_sf.numCoefs_v();
  int order_u = spline_sf.order_u();
  int order_v = spline_sf.order_v();

  int mid_knot_ind_u = floor((num_coefs_u + order_u)/2);
  int mid_knot_ind_v = floor((num_coefs_u + order_u)/2);

  bool refine_at_line = false;//true;
  double parval = (refine_at_line) ? spline_sf.basis_v().begin()[mid_knot_ind_v] :
      0.5*(spline_sf.basis_v().begin()[mid_knot_ind_v] + spline_sf.basis_v().begin()[mid_knot_ind_v+1]);
  double start = spline_sf.basis_u().begin()[mid_knot_ind_u];
  double end = umax;
  int mult = 1;
  std::cout << "parval: " << parval << ", start: " << start << " end: " << end << std::endl;
  lr_spline_sf->refine((dir==0) ? Go::YFIXED : Go::XFIXED, parval, start, end, mult);
  lr_spline_sf->refine((dir==0) ? Go::YFIXED : Go::XFIXED, parval, start, end, mult);

  double max_dist_post_ref = maxDist(spline_sf, *lr_spline_sf, nmb_samples_u, nmb_samples_v);
  std::cout << "Max dist between input and converted surface: " << max_dist_post_ref << std::endl;

  // We write to screen the number of element and basis functions.
  num_basis_funcs = lr_spline_sf->numBasisFunctions();
  num_elem = lr_spline_sf->numElements();
  std::cout << "num_basis_funcs: " << num_basis_funcs << ", num_elem: " << num_elem << std::endl;

#ifndef NDEBUG
  std::ofstream lrsf_grid_ps("tmp/lrsf_grid.ps");
//  writePostscriptMesh(*lrsf);
  writePostscriptMesh(*lr_spline_sf, lrsf_grid_ps);
#endif NDEBUG

  lr_spline_sf->writeStandardHeader(fileout);
  lr_spline_sf->write(fileout);

}


double maxDist(const Go::SplineSurface& spline_sf,
	       const Go::LRSplineSurface& lr_spline_sf,
	       int nmb_samples_u, int nmb_samples_v)
{
    // Assuming the domain is the same.
    double umin = lr_spline_sf.startparam_u();
    double umax = lr_spline_sf.endparam_u();
    double vmin = lr_spline_sf.startparam_v();
    double vmax = lr_spline_sf.endparam_v();
    double ustep = (umax - umin)/((double)nmb_samples_u + 1);
    double vstep = (vmax - vmin)/((double)nmb_samples_v + 1);
    Go::Point go_pt, lr_pt;
    double max_dist = -1.0;
#ifndef NDEBUG
    double max_dist_u = 0.0;
    double max_dist_v = 0.0;
    Go::Point max_go_pt, max_lr_pt;
#endif
    for (int kj = 0; kj < nmb_samples_v; ++kj)
    {
	double vpar = vmin + kj*vstep;
	for (int ki = 0; ki < nmb_samples_u; ++ki)
	{
	    double upar = umin + ki*ustep;
	    spline_sf.point(go_pt, upar, vpar);
	    lr_spline_sf.point(lr_pt, upar, vpar);
	    double dist = go_pt.dist(lr_pt);
	    if (dist > max_dist)
	    {
		max_dist = dist;
#ifndef NDEBUG
		max_dist_u = upar;
		max_dist_v = vpar;
		max_go_pt = go_pt;
		max_lr_pt = lr_pt;
#endif
	    }
	}
    }    

#ifndef NDEBUG
    std::cout << "max_dist_u = " << max_dist_u << ", max_dist_v = " << max_dist_v << std::endl;
    std::cout << "max_go_pt = " << max_go_pt << std::endl;
    std::cout << "max_lr_pt = " << max_lr_pt << std::endl;
#endif
    return max_dist;


}
