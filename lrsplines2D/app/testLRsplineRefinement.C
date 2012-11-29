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


using namespace Go;

// Assuming the domain is the same.
double maxDist(const Go::ParamSurface* param_sf,
	       const Go::LRSplineSurface& lr_spline_sf,
	       int nmb_samples_u, int nmb_samples_v);


int main(int argc, char *argv[])
{
  if (argc != 3)
  {
      std::cout << "Usage: input_(lr_)spline_sf.g2 ref_lr_spline_sf.g2" << std::endl;
      return -1;
  }

  std::ifstream filein(argv[1]); // Input (regular) spline surface.
  std::ofstream fileout(argv[2]); // Input (regular) spline surface.
  shared_ptr<Go::ParamSurface> input_sf;
  Go::ObjectHeader header;
  filein >> header;
  shared_ptr<Go::LRSplineSurface> lr_spline_sf(new Go::LRSplineSurface());
  int num_coefs_u;
  int num_coefs_v;
  int order_u;
  int order_v;
  LRSplineSurface::Refinement2D ref_v, ref_u;
  // ref.d = Go::XFIXED; // Inserting a v/x knot.
  ref_u.multiplicity = ref_v.multiplicity = 1;
  LRSplineSurface lrsf;
  if (header.classType() == Go::Class_SplineSurface)
  {
      shared_ptr<Go::SplineSurface> spline_sf(new Go::SplineSurface());
      filein >> *spline_sf;
      num_coefs_u = spline_sf->numCoefs_u();
      num_coefs_v = spline_sf->numCoefs_v();
      order_u = spline_sf->order_u();
      order_v = spline_sf->order_v();

      input_sf = spline_sf;
      puts("Now we convert from SplineSurface to LRSplineSurface!");
      double knot_tol = 1e-05;
      lr_spline_sf = shared_ptr<Go::LRSplineSurface>(new Go::LRSplineSurface(spline_sf->clone(), knot_tol));
      puts("Done converting!");

#ifndef NDEBUG
      {
	std::vector<shared_ptr<LRBSpline2D> > bas_funcs;
	for (auto iter = lr_spline_sf->basisFunctionsBegin(); iter != lr_spline_sf->basisFunctionsEnd(); ++iter)
	  {
	    bas_funcs.push_back((*iter).second);
	  }
	puts("Remove when done debugging!");
      }
#endif

      double umax = lr_spline_sf->paramMax(Go::XFIXED);
      double vmax = lr_spline_sf->paramMax(Go::YFIXED);
      int mid_knot_ind_u = floor((num_coefs_u + order_u)/2);
      int mid_knot_ind_v = floor((num_coefs_u + order_u)/2);
      bool refine_at_line = false;//true;
//  std::cout << "parval: " << parval << ", start: " << start << " end: " << end << std::endl;
      ref_v.kval = (refine_at_line) ? spline_sf->basis_v().begin()[mid_knot_ind_v] :
	  0.5*(spline_sf->basis_v().begin()[mid_knot_ind_v] + spline_sf->basis_v().begin()[mid_knot_ind_v+1]);
      ref_v.start = spline_sf->basis_u().begin()[mid_knot_ind_u];
      ref_v.end = umax;
      ref_v.d = Go::YFIXED;

      ref_u.kval = (refine_at_line) ? spline_sf->basis_u().begin()[mid_knot_ind_u] :
	  0.5*(spline_sf->basis_u().begin()[mid_knot_ind_u] + spline_sf->basis_u().begin()[mid_knot_ind_u+1]);
      ref_u.start = spline_sf->basis_v().begin()[mid_knot_ind_v];
      ref_u.end = vmax;
      ref_u.d = Go::XFIXED;

  }
  else if (header.classType() == Go::Class_LRSplineSurface)
  {
      filein >> lrsf;
      lr_spline_sf = shared_ptr<LRSplineSurface>(lrsf.clone());
      puts("Refining a LRSplineSurface only working on specific cases, hardcoded values.");
      ref_v.kval = 0.625;
      ref_v.start = 0.5;
      ref_v.end = 1.0;
      ref_v.d = Go::YFIXED;
      ref_u = ref_v;
      ref_u.d = Go::XFIXED;
      input_sf = shared_ptr<ParamSurface>(lr_spline_sf->clone());

#ifndef NDEBUG
      {
	std::vector<shared_ptr<LRBSpline2D> > bas_funcs;
	for (auto iter = lr_spline_sf->basisFunctionsBegin(); iter != lr_spline_sf->basisFunctionsEnd(); ++iter)
	  {
	    bas_funcs.push_back((*iter).second);
	  }
	puts("Remove when done debugging!");
      }
#endif

  }
  else
  {
      puts("Unexpected geometry input.");
      return -1;
  }


#ifndef NDEBUG
  {
    std::vector<shared_ptr<LRBSpline2D> > bas_funcs;
    for (auto iter = lr_spline_sf->basisFunctionsBegin(); iter != lr_spline_sf->basisFunctionsEnd(); ++iter)
      {
	bas_funcs.push_back((*iter).second);
      }
    puts("Remove when done debugging!");
  }
#endif

  // We test to see if conversion was correct.
  int nmb_samples_u = 3;//101; // Rather random.
  int nmb_samples_v = 3;//36;
  double max_dist = maxDist(input_sf.get(), *lr_spline_sf, nmb_samples_u, nmb_samples_v);
  std::cout << "Max dist between input and converted surface: " << max_dist << std::endl;

  // We write to screen the number of element and basis functions.
  int num_basis_funcs = lr_spline_sf->numBasisFunctions();
  int num_elem = lr_spline_sf->numElements();
  std::cout << "Status prior to refinement: num_basis_funcs: " << num_basis_funcs << ", num_elem: " << num_elem << std::endl;

  // We hardcode a refinement.
  // Hardcoded values for debugging.

//  lr_spline_sf->refine((dir==0) ? Go::YFIXED : Go::XFIXED, parval, start, end, mult);
  std::vector<LRSplineSurface::Refinement2D> refs_single, refs_multi;
//  refs.push_back(ref);
  // LRSplineSurface::Refinement2D ref2 = ref;
  // ref2.kval = 0.4;
  // LRSplineSurface::Refinement2D ref3 = ref;
  // ref3.kval = 0.65;
//  LRSplineSurface::Refinement2D ref4 = ref;
  //ref4.d = Go::YFIXED;

//  refs.push_back(ref2);
//  refs.push_back(ref3);
  refs_single.push_back(ref_u);
  refs_single.push_back(ref_v);
  refs_multi.push_back(ref_u);
  refs_multi.push_back(ref_v);



#ifndef NDEBUG
  {
    std::vector<shared_ptr<LRBSpline2D> > bas_funcs;
    for (auto iter = lr_spline_sf->basisFunctionsBegin(); iter != lr_spline_sf->basisFunctionsEnd(); ++iter)
      {
	bas_funcs.push_back((*iter).second);
      }
    puts("Remove when done debugging!");
  }
#endif

  bool ref_single = true;

  if (ref_single)
    {

      for (uint ki = 0; ki < refs_single.size(); ++ki)
	{
	  lr_spline_sf->refine(refs_single[ki]);
	}
//  lr_spline_sf->refine(ref);
//  lr_spline_sf_multi->refine(refs);
//  lr_spline_sf->refine((dir==0) ? Go::YFIXED : Go::XFIXED, parval, start, end, mult);


#ifndef NDEBUG
      {
	std::vector<shared_ptr<LRBSpline2D> > bas_funcs;
	for (auto iter = lr_spline_sf->basisFunctionsBegin(); iter != lr_spline_sf->basisFunctionsEnd(); ++iter)
	  {
	    bas_funcs.push_back((*iter).second);
	  }
	puts("Remove when done debugging!");
      }
#endif

      double max_dist_post_ref = maxDist(input_sf.get(), *lr_spline_sf, nmb_samples_u, nmb_samples_v);
      std::cout << "Max dist input and ref surface (ref one at the time): " << max_dist_post_ref << std::endl;

      // We write to screen the number of element and basis functions.
      num_basis_funcs = lr_spline_sf->numBasisFunctions();
      num_elem = lr_spline_sf->numElements();
      std::cout << "Ref one at the time: num_basis_funcs: " << num_basis_funcs << ", num_elem: " << num_elem << std::endl;


// #ifndef NDEBUG
      std::ofstream lrsf_grid_ps("tmp/lrsf_grid.ps");
//  writePostscriptMesh(*lrsf);
// #endif NDEBUG


      writePostscriptMesh(*lr_spline_sf, lrsf_grid_ps);
      std::ofstream fileout2("tmp/ref_lr_single.g2");
      lr_spline_sf->writeStandardHeader(fileout2);
      lr_spline_sf->write(fileout2);

      lr_spline_sf->writeStandardHeader(fileout);
      lr_spline_sf->write(fileout);

    }
  else
    {
      shared_ptr<LRSplineSurface> lr_spline_sf_multi(new LRSplineSurface());
      *lr_spline_sf_multi = *lr_spline_sf;

#ifndef NDEBUG
      {
	std::vector<shared_ptr<LRBSpline2D> > bas_funcs;
	for (auto iter = lr_spline_sf_multi->basisFunctionsBegin(); iter != lr_spline_sf_multi->basisFunctionsEnd(); ++iter)
	  {
	    bas_funcs.push_back((*iter).second);
	  }
	puts("Remove when done debugging!");
      }
#endif

      lr_spline_sf_multi->refine(refs_multi);


#ifndef NDEBUG
      {
	std::vector<shared_ptr<LRBSpline2D> > bas_funcs;
	for (auto iter = lr_spline_sf_multi->basisFunctionsBegin(); iter != lr_spline_sf_multi->basisFunctionsEnd(); ++iter)
	  {
	    bas_funcs.push_back((*iter).second);
	  }
	puts("Remove when done debugging!");
      }
#endif

      double max_dist_post_ref_multi_ref = maxDist(input_sf.get(), *lr_spline_sf_multi, nmb_samples_u, nmb_samples_v);
      std::cout << "Max dist input and ref surface: " << max_dist_post_ref_multi_ref << std::endl;


      int num_basis_funcs_multi = lr_spline_sf_multi->numBasisFunctions();
      int num_elem_multi = lr_spline_sf_multi->numElements();
      std::cout << "Ref using a vector of refinements: num_basis_funcs: " << num_basis_funcs_multi <<
	", num_elem_multi: " << num_elem_multi << std::endl;


      std::ofstream fileout3("tmp/ref_lr_multi.g2");
      lr_spline_sf_multi->writeStandardHeader(fileout3);
      lr_spline_sf_multi->write(fileout3);

      std::ofstream lrsf_multi_grid_ps("tmp/lrsf_grid_multi.ps");
     writePostscriptMesh(*lr_spline_sf_multi, lrsf_multi_grid_ps);

    }

}


double maxDist(const Go::ParamSurface* param_sf,
	       const Go::LRSplineSurface& lr_spline_sf,
	       int nmb_samples_u, int nmb_samples_v)
{
    // Assuming the domain is the same.
    double umin = lr_spline_sf.startparam_u();
    double umax = lr_spline_sf.endparam_u();
    double vmin = lr_spline_sf.startparam_v();
    double vmax = lr_spline_sf.endparam_v();
    double ustep = (umax - umin)/((double)nmb_samples_u - 1);
    double vstep = (vmax - vmin)/((double)nmb_samples_v - 1);
    Go::Point go_pt(3), lr_pt(3);
    double max_dist = -1.0;
// #ifndef NDEBUG
    double max_dist_u = 0.0;
    double max_dist_v = 0.0;
    Go::Point max_go_pt(3), max_lr_pt(3);
// #endif
    for (int kj = 0; kj < nmb_samples_v; ++kj)
    {
	double vpar = vmin + kj*vstep;
	for (int ki = 0; ki < nmb_samples_u; ++ki)
	{
	    double upar = umin + ki*ustep;
	    param_sf->point(go_pt, upar, vpar);
	    lr_spline_sf.point(lr_pt, upar, vpar);
	    double dist = go_pt.dist(lr_pt);
	    if (dist > max_dist)
	    {
		max_dist = dist;
// #ifndef NDEBUG
		max_dist_u = upar;
		max_dist_v = vpar;
		max_go_pt = go_pt;
		max_lr_pt = lr_pt;
// #endif
	    }
	}
    }    

// #ifndef NDEBUG
    std::cout << "max_dist: u = " << max_dist_u << ", v = " << max_dist_v << std::endl;
    std::cout << "max_go_pt = " << max_go_pt << std::endl;
    std::cout << "max_lr_pt = " << max_lr_pt << std::endl;
// #endif

    return max_dist;


}
