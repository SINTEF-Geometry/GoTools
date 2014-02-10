/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

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
  std::ofstream fileout(argv[2]); // Output refined lr_spline_sf.
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
  double knot_tol = 1e-05;
  shared_ptr<Go::SplineSurface> spline_sf(new Go::SplineSurface());
  if (header.classType() == Go::Class_SplineSurface)
  {
      filein >> *spline_sf;

#if 0
      MESSAGE("Setting parameter domain to the unit square.");
      spline_sf->setParameterDomain(0.0, 1.0, 0.0, 1.0);
#endif


#if 1
      // Making the input k-regular results in an error for the refinement (max_dist is 3e-03).
      MESSAGE("Deactivated making the surface k-regular. Should test that for cone_sf.g2");
#else
      MESSAGE("Making sure the input surface is k-regular.");
      // Having some problems with a rational case (cylinder_sf.g2) ...

      // Example is already k-regular.
      spline_sf->makeSurfaceKRegular();
      // Swpping the parameter directions does not help.
//      spline_sf->swapParameterDirection();
#endif

#if 0
      // We refine the input spline surface prior to converting it to the lr format.
      const BsplineBasis& bas_u = spline_sf->basis_u();
      const auto knots = bas_u.begin();
      const int deg = bas_u.order() - 1;
      double w1 = 0.5;//34;
      double insert_u = w1*knots[deg] + (1.0 - w1)*knots[deg+1];
      spline_sf->insertKnot_u(insert_u);
#endif

#if 0
      order_u = spline_sf->order_u();
      int raise_u = 0;
      int raise_v = 0;
      if (order_u < 4)
	  raise_u = 4 - order_u;
      order_v = spline_sf->order_v();
      if (order_v < 4)
	  raise_v = 4 - order_v;
      if (raise_u > 0 || raise_v > 0)
	  spline_sf->raiseOrder(raise_u, raise_v);
#endif

      num_coefs_u = spline_sf->numCoefs_u();
      num_coefs_v = spline_sf->numCoefs_v();
      order_u = spline_sf->order_u();
      order_v = spline_sf->order_v();

      input_sf = spline_sf;
      puts("Now we convert from SplineSurface to LRSplineSurface!");
      lr_spline_sf = shared_ptr<Go::LRSplineSurface>(new Go::LRSplineSurface(spline_sf->clone(), knot_tol));
      puts("Done converting!");

      std::ofstream fileout_conv("tmp/file_conv.g2");
      lr_spline_sf->writeStandardHeader(fileout_conv);
      lr_spline_sf->write(fileout_conv);

#ifndef NDEBUG
      {
	std::vector<LRBSpline2D*> bas_funcs;
	for (auto iter = lr_spline_sf->basisFunctionsBegin(); iter != lr_spline_sf->basisFunctionsEnd(); ++iter)
	  {
	    bas_funcs.push_back(iter->second.get());
	  }
	puts("Remove when done debugging!");
      }
#endif

      const Mesh2D& mesh = lr_spline_sf->mesh();
      double umin = mesh.minParam(Go::XFIXED);
      int umin_ind = mesh.getKnotIdx(Go::XFIXED, umin, knot_tol);
      double w1 = 0.34; // Random value in the range (0.0, 1.0).
      ref_u.kval = w1*mesh.kval(Go::XFIXED, umin_ind) + (1.0 - w1)*mesh.kval(Go::XFIXED, umin_ind + 1);
      ref_u.start = mesh.minParam(Go::YFIXED);
      ref_u.end = mesh.maxParam(Go::YFIXED);
      ref_u.d = Go::XFIXED;
      double vmin = mesh.minParam(Go::YFIXED);
      int vmin_ind = mesh.getKnotIdx(Go::YFIXED, vmin, knot_tol);
      double w2 = 0.72;// Random value in the range (0.0, 1.0).
      ref_v.kval = w2*mesh.kval(Go::YFIXED, vmin_ind) + (1.0 - w2)*mesh.kval(Go::YFIXED, vmin_ind + 1);
      ref_v.start = mesh.minParam(Go::XFIXED);
      ref_v.end = mesh.maxParam(Go::XFIXED);
      ref_v.d = Go::YFIXED;

  }
  else if (header.classType() == Go::Class_LRSplineSurface)
  {
      filein >> lrsf;
      lr_spline_sf = shared_ptr<LRSplineSurface>(lrsf.clone());
      puts("Refining a LRSplineSurface only working on specific cases, hardcoded values.");

      const Mesh2D& mesh = lr_spline_sf->mesh();
      ref_u.kval = 0.5*(mesh.kval(Go::XFIXED, 0) + mesh.kval(Go::XFIXED, 1));
      ref_u.start = mesh.minParam(Go::YFIXED);
      ref_u.end = mesh.maxParam(Go::YFIXED);
      ref_u.d = Go::XFIXED;
      ref_v.kval = 0.5*(mesh.kval(Go::YFIXED, 0) + mesh.kval(Go::YFIXED, 1));
      ref_v.start = mesh.minParam(Go::XFIXED);
      ref_v.end = mesh.maxParam(Go::XFIXED);
      ref_v.d = Go::YFIXED;

      input_sf = shared_ptr<ParamSurface>(lr_spline_sf->clone());

#if 0//ndef NDEBUG
      {
	std::vector<LRBSpline2D*> bas_funcs;
	for (auto iter = lr_spline_sf->basisFunctionsBegin(); iter != lr_spline_sf->basisFunctionsEnd(); ++iter)
	  {
	    bas_funcs.push_back(iter->second.get());
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


#if 0//ndef NDEBUG
  {
    std::vector<LRBSpline2D*> bas_funcs;
    for (auto iter = lr_spline_sf->basisFunctionsBegin(); iter != lr_spline_sf->basisFunctionsEnd(); ++iter)
      {
	bas_funcs.push_back(iter->second.get());
      }
    puts("Remove when done debugging!");
  }
#endif

  // We test to see if conversion was correct.
  int nmb_samples_u = 101; // Rather random.
  int nmb_samples_v = 36;
  double max_dist = maxDist(input_sf.get(), *lr_spline_sf, nmb_samples_u, nmb_samples_v);
  std::cout << "Max dist between input and converted surface: " << max_dist << std::endl;

  // We write to screen the number of element and basis functions.
  int num_basis_funcs = lr_spline_sf->numBasisFunctions();
  int num_elem = lr_spline_sf->numElements();
  std::cout << "Status prior to refinement: num_basis_funcs: " << num_basis_funcs << ", num_elem: " << num_elem << std::endl;

  std::vector<LRSplineSurface::Refinement2D> refs_single, refs_multi;

#if 1
  refs_single.push_back(ref_u);
  refs_multi.push_back(ref_u);
#endif
#if 1
  refs_single.push_back(ref_v);
  refs_multi.push_back(ref_v);
#endif
#if 0
  LRSplineSurface::Refinement2D ref_v2 = ref_v;
  ref_v2.kval = -0.5;
  refs_single.push_back(ref_v2);
  refs_multi.push_back(ref_v2);
#endif

#if 0//ndef NDEBUG
  MESSAGE("Clearing the refinement vectors, for debugging!");
  refs_single.clear();
  refs_multi.clear();
#endif


#if 0//ndef NDEBUG
  {
    std::vector<LRBSpline2D*> bas_funcs;
    for (auto iter = lr_spline_sf->basisFunctionsBegin(); iter != lr_spline_sf->basisFunctionsEnd(); ++iter)
      {
	bas_funcs.push_back(iter->second.get());
      }
    puts("Remove when done debugging!");
  }
#endif

  shared_ptr<LRSplineSurface> lr_spline_sf_multi(new LRSplineSurface());
#if 1
  // @@sbr201301 Still having some problems with the copy constructor, still pointing to some shared data.
  *lr_spline_sf_multi = *lr_spline_sf;
//  lr_spline_sf_multi = shared_ptr<Go::LRSplineSurface>(new Go::LRSplineSurface(*lr_spline_sf));
#if 1//ndef NDEBUG
  {
    std::vector<LRBSpline2D*> bas_funcs;
    for (auto iter = lr_spline_sf_multi->basisFunctionsBegin(); iter != lr_spline_sf_multi->basisFunctionsEnd(); ++iter)
      {
	bas_funcs.push_back((*iter).second.get());
      }
    puts("Remove when done debugging!");
  }
#endif
#else
  MESSAGE("Deactivated copy constructor usage for LRSplineSurface due to problems.");
  lr_spline_sf_multi = shared_ptr<Go::LRSplineSurface>(new Go::LRSplineSurface(spline_sf->clone(), knot_tol));
#endif

  bool refine_single = true;
  if (refs_single.size() == 0)
      ;//refine_single = false;
  if (refine_single)
    {

#if 1//ndef NDEBUG
      {
	std::vector<LRBSpline2D*> bas_funcs;
	for (auto iter = lr_spline_sf_multi->basisFunctionsBegin(); iter != lr_spline_sf_multi->basisFunctionsEnd(); ++iter)
	  {
	    bas_funcs.push_back((*iter).second.get());
	  }
	puts("Remove when done debugging!");
      }
#endif

      for (unsigned int ki = 0; ki < refs_single.size(); ++ki)
	{
#ifndef NDEBUG
	  std::cout << "DEBUG: ki = " << ki << std::endl;
#endif
	  lr_spline_sf->refine(refs_single[ki]);
	}

#if 1//ndef NDEBUG
      {
	std::vector<LRBSpline2D*> bas_funcs;
	for (auto iter = lr_spline_sf_multi->basisFunctionsBegin(); iter != lr_spline_sf_multi->basisFunctionsEnd(); ++iter)
	  {
	    bas_funcs.push_back((*iter).second.get());
	  }
	puts("Remove when done debugging!");
      }
#endif

#if 0//ndef NDEBUG
      {
	std::vector<LRBSpline2D*> bas_funcs;
	for (auto iter = lr_spline_sf->basisFunctionsBegin(); iter != lr_spline_sf->basisFunctionsEnd(); ++iter)
	  {
	    bas_funcs.push_back(iter->second.get());
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

      std::ofstream lrsf_grid_ps("tmp/lrsf_grid.ps");

      writePostscriptMesh(*lr_spline_sf, lrsf_grid_ps);
      std::ofstream fileout2("tmp/ref_lr_single.g2");
      lr_spline_sf->writeStandardHeader(fileout2);
      lr_spline_sf->write(fileout2);

      lr_spline_sf->writeStandardHeader(fileout);
      lr_spline_sf->write(fileout);

      num_basis_funcs = lr_spline_sf->numBasisFunctions();
      num_elem = lr_spline_sf->numElements();
      std::cout << "Status after refinement (single): num_basis_funcs: " << num_basis_funcs << ", num_elem: " << num_elem << std::endl;

    }

  bool refine_multi = true;
  if (refs_multi.size() == 0)
      refine_multi = false;
  if (refine_multi)
    {

#if 1//ndef NDEBUG
      {
	std::vector<LRBSpline2D*> bas_funcs;
	for (auto iter = lr_spline_sf_multi->basisFunctionsBegin(); iter != lr_spline_sf_multi->basisFunctionsEnd(); ++iter)
	  {
	    bas_funcs.push_back(iter->second.get());
	  }
	puts("Remove when done debugging!");
      }
#endif

      lr_spline_sf_multi->refine(refs_multi);

      num_basis_funcs = lr_spline_sf_multi->numBasisFunctions();
      num_elem = lr_spline_sf_multi->numElements();
      std::cout << "Status after refinement (multi): num_basis_funcs: " << num_basis_funcs << ", num_elem: " << num_elem << std::endl;

#if 0//ndef NDEBUG
      {
	std::vector<LRBSpline2D*> bas_funcs;
	for (auto iter = lr_spline_sf_multi->basisFunctionsBegin(); iter != lr_spline_sf_multi->basisFunctionsEnd(); ++iter)
	  {
	    bas_funcs.push_back(iter->second.get());
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
