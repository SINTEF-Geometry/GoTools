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

#include <fstream>
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/SurfaceOnVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/trivariatemodel/ftVolumeTools.h"
#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/geometry/Utils.h"

using namespace Go;
using std::cout;
using std::endl;
using std::ifstream;
using std::vector;

// This program is intended to illustrate some GoTools functionality that
// can be useful in the context of performing isogeometric analysis on 
// trimmed volumes. 

int main(int argc, char* argv[] )
{
  if (argc != 2)
      cout << "Usage: " << "g2 Brep file" << endl;

  ifstream infile(argv[1]);
  ALWAYS_ERROR_IF(infile.bad(), "Bad or no input filename");

  // The tolerances must be set according to the properties of the model.
  // The neighbour tolerance must be smaller than the smallest entity in the
  // model, but larger than the largest gap.
  // The gap tolerance must be smaller than the neighbour tolerance
  double gap = 0.001; //0.001;
  double neighbour = 0.01; //0.01;
  double kink = 0.01;
  double approxtol = 0.001;
  int degree = 3;

  double eps = 1.0e-8;  // Computational tolerance

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  CompositeModel *model = factory.createFromG2(infile);

  shared_ptr<SurfaceModel> sfmodel = 
    shared_ptr<SurfaceModel>(dynamic_cast<SurfaceModel*>(model));
  if (!sfmodel.get())
    {
      std::cout << "No input model read" << std::endl;
      exit(-1);
    }
 
  if (sfmodel->nmbBoundaries() > 0)
    {
      std::cout << "Not a brep solid. Consider increasing the neighbour tolerance" << std::endl;
      exit(-1);
    }
      
  bool isOK = sfmodel->checkShellTopology();
  std::cout << "Shell topology: " << isOK << std::endl;

  // A spline volume will be created, but no attempt to orient it in space or adapt to
  // certain surfaces is performed
  shared_ptr<ftVolume> vol(new ftVolume(sfmodel));

  // Check if the volume can be represented without any trimming
  if (vol->isRegularized())
    {
      int degree = 3;  // cubic spline volume
      bool untrimmed = vol->untrimRegular(degree);  // An approximation is performed
      std::cout << "Untrimming performed: " << untrimmed << std::endl;
    }

  // Fetch spline volume
  bool is_spline = vol->isSpline();
  std::cout << "Spline? " << is_spline << std::endl;
  shared_ptr<ParamVolume> par_vol = vol->getVolume();
  SplineVolume *spline_vol = par_vol->asSplineVolume();
  if (spline_vol)
    {
      // Refine volume
      // First fetch parameter domain information
      // Sequence: umin, umax, vmin, vmax, wmin, wmax
      const Array<double, 6> pardom = spline_vol->parameterSpan();

      // Insert knots. The initial underlying volume will be without
      // inner knots. Normally, one would take existing knot vector into
      // considerations when inserting new knots
      vector<double> knots_u(8);
      double udel = (pardom[1] - pardom[0])/(double)9;
      double tt = pardom[0]+udel;
      for (int ki=0; ki<8; ++ki, tt+=udel)
	knots_u[ki] = tt;
      vector<double> knots_v(5);
      double vdel = (pardom[3] - pardom[2])/(double)6;
      tt = pardom[2]+vdel;
      for (int ki=0; ki<5; ++ki, tt+=vdel)
	knots_v[ki] = tt;
      vector<double> knots_w(4);
      double wdel = (pardom[5] - pardom[4])/(double)5;
      tt = pardom[4]+wdel;
      for (int ki=0; ki<4; ++ki, tt+=wdel)
	knots_w[ki] = tt;
      spline_vol->insertKnot(0, knots_u);
      spline_vol->insertKnot(1, knots_v);
      spline_vol->insertKnot(2, knots_w);

      // Insert one knot
      spline_vol->insertKnot(2, 0.5*(knots_w[2]+knots_w[3]));
 
      // Check element status
      // -1 = not a spline volume, 0 = outside trimmed volume, 1 = boundary
      // element, 2 = completely inside
      int elem_ix;
      std::cout << "Give element index" << std::endl;
      std::cin >> elem_ix;
      int element_status = vol->ElementBoundaryStatus(elem_ix);
      std::cout << "Element status: " << element_status << std::endl;

      // -1 = not spline, 0 = not on boundary, 1 = on boundary
      int on_bd = vol->ElementOnBoundary(elem_ix);
      std::cout << "On boundary: " << on_bd << std::endl;

      // Perform grid evaluation. In the case of a trimmed surface, the
      // grid will extend beyond the trimmed surface
      // Points and derivatives
      // Also other alternatives exist, check SplineVolume.h
     // Define grid in the internal of the volume (non-trimmed)
      int nmb_u = 5;
      int nmb_v = 5;
      int nmb_w = 5;
      vector<double> param_u(nmb_u);
      vector<double> param_v(nmb_v);
      vector<double> param_w(nmb_w);
      double del_u = (pardom[1] - pardom[0])/(double)(nmb_u+1);
      double del_v = (pardom[3] - pardom[2])/(double)(nmb_v+1);
      double del_w = (pardom[5] - pardom[4])/(double)(nmb_w+1);
      int kj, kr;
      double par;
      for (kj=0, par=pardom[0]+del_u; kj<nmb_u; ++kj, par+=del_u)
	param_u[kj] = par;
      for (kj=0, par=pardom[2]+del_v; kj<nmb_v; ++kj, par+=del_v)
	param_v[kj] = par;
      for (kj=0, par=pardom[4]+del_w; kj<nmb_w; ++kj, par+=del_w)
	param_w[kj] = par;

      vector<double> points;
      vector<double> derivs_u;
      vector<double> derivs_v;
      vector<double> derivs_w;
      spline_vol->gridEvaluator(param_u, param_v, param_w,
				points, derivs_u, derivs_v, derivs_w);

      // Basis functions including 1. derivatives
      // See SplineVolume.h for the struct BasisDerivs
      vector<BasisDerivs> basis_derivs1;
      spline_vol->computeBasisGrid(param_u, param_v, param_w,
				   basis_derivs1);

      // Including 2. derivatives
      // See SplineVolume.h for the struct BasisDerivs2
      vector<BasisDerivs2> basis_derivs2;
      spline_vol->computeBasisGrid(param_u, param_v, param_w,
				   basis_derivs2);

      int stop_debug_face = 1;
      
    }

  // Fetch parameter domain information. This function works also for 
  // other volume types
  const Array<double, 6> pardom = par_vol->parameterSpan();

  // Evaluate in midpoint
  Point pt;
  par_vol->point(pt, 0.5*(pardom[0]+pardom[1]),
		 0.5*(pardom[2]+pardom[3]), 0.5*(pardom[4]+pardom[5]));

  
  // Evaluate point and derivatives including 2nd order
  vector<Point> der(10);  // pos, der_u, der_v, der_w, der_uu, der_uv, der_uw, der_vv, der_vw, der_ww
  int nmb_derivs = 2;
  par_vol->point(der, 0.5*(pardom[0]+pardom[1]),
		 0.5*(pardom[2]+pardom[3]), 0.5*(pardom[4]+pardom[5]), 2);

  // Point in volume test
  bool inside = vol->isInside(pt);
  std::cout << "Point in volume test for midpoint: " << inside << std::endl;

  // Alternatively test parameter tripple
  inside = vol->ParamInVolume((2.0*pardom[0]+pardom[1])/3.0,
			      (2.0*pardom[2]+pardom[3])/3.0, 
			      (2.0*pardom[4]+pardom[5])/3.0);
  std::cout << "Point in volume test for point one third in in all directions: " << inside << std::endl;

  // Closest point (compute parameter value) of underlying volume
  double clo_u, clo_v, clo_w, clo_dist;
  Point clo_pt;
  par_vol->closestPoint(pt, clo_u, clo_v, clo_w, clo_pt, clo_dist, eps);
  std::cout << "Closest point parameters: (" << clo_u << ", " << clo_v << ", " << clo_w << ")" << std::endl;
  std::cout << "Distance: " << clo_dist << std::endl;

  // Closest point of trimmed volume
  double clo_u2, clo_v2, clo_w2, clo_dist2;
  Point clo_pt2;
  vol->closestPoint(pt, clo_u2, clo_v2, clo_w2, clo_pt2, clo_dist2, eps);
  std::cout << "Closest point parameters(trimmed): (" << clo_u2 << ", " << clo_v2 << ", " << clo_w2 << ")" << std::endl;
  std::cout << "Distance: " << clo_dist2 << std::endl;

  // Check surface category as volume boundary surfaces
  // First fetch outer boundary
  shared_ptr<SurfaceModel> shell = vol->getOuterShell();

  // Closest point of shell
  double clo_dist3;
  Point clo_pt3;
  double clo_par[2];
  int clo_idx;
  shell->closestPoint(pt, clo_pt3, clo_idx, clo_par, clo_dist3);
  std::cout << "Closest point face index: " << clo_idx << std::endl;
  std::cout << "Closest point parameters(shell): (" << clo_par[0] << ", " << clo_par[1] << ")" << std::endl;
  std::cout << "Distance: " << clo_dist3 << std::endl;

  // Number of surfaces in shell
  int nmb = shell->nmbEntities();
  for (int ki=0; ki<nmb; ++ki)
    {
      // Fetch face
      shared_ptr<ftSurface> face = shell->getFace(ki);

      // Check category
      // = -1 : not a volume boundary (trimming face)
      // = 0 : umin (the surface may be trimmed, it may not cover the entire volume boundar,
      // it may not have the same representation as the volume boundary surface, 
      // but it is coincident
      // = 1 : umax
      // = 2 : vmin
      // = 3 : vmax
      // = 4 : wmin
      // = 5 : wmax
      int bd_status = ftVolumeTools::boundaryStatus(vol.get(), face, gap);
      std::cout << "Boundary surface nr " << ki << ", boundary status: " << bd_status << std::endl;

      // Check element status
      // -1 = not a spline surface, 0 = outside trimmed surface, 1 = boundary
      // element, 2 = completely inside
      int sf_elem_ix = 0;  // Hardcoded just to show call
      // The functionality is also available through the associated surface
      int sf_elem_stat = face->ElementBoundaryStatus(sf_elem_ix, gap);
      std::cout << "Element status surface " << ki << ": " << sf_elem_stat << std::endl;

      // Fetch associated surface
      shared_ptr<ParamSurface> surf = face->surface();

      // Parameter domain containing the surface. In case of trimming, it might be larger
      RectDomain dom = surf->containingDomain();

      // Define grid in the internal of the surface (non-trimmed)
      int nmb_u = 5;
      int nmb_v = 5;
      vector<double> param_u(nmb_u);
      vector<double> param_v(nmb_v);
      double del_u = (dom.umax() - dom.umin())/(double)(nmb_u+1);
      double del_v = (dom.vmax() - dom.vmin())/(double)(nmb_v+1);
      int kj, kr;
      double par;
      for (kj=0, par=dom.umin()+del_u; kj<nmb_u; ++kj, par+=del_u)
	param_u[kj] = par;
      for (kj=0, par=dom.vmin()+del_v; kj<nmb_v; ++kj, par+=del_v)
	param_v[kj] = par;

       // Fetch associated spline surface. If the surface is trimmed, the underlying
      // non-trimmed surface is returned. If the surface is not a spline surface (this
      // can happen in theory), a null pointer is returned.
      SplineSurface* spline_surf = surf->asSplineSurface();

     // Check if the grid points are inside the trimmed surface. In that case
      // evaluate basis functions, point and surface normal. The surface is oriented
      // such that the normal points out of the volume.
      for (kj=0; kj<nmb_v; ++kj)
	for (kr=0; kr<nmb_u; ++kr)
	  {
	    if (face->pointInFace(param_u[kr], param_v[kj], gap))
	      {
		// Evaluate position
		Point pos = face->point(param_u[kr], param_v[kj]);

		// Evaluate normal
		Point normal = face->normal(param_u[kr], param_v[kj]);

		if (spline_surf)
		  {
		    // Compute value and 1. derivative with respect
		    // to the B-spline basis functions
		    vector<double> basisValues;
		    vector<double> basisDerivs_u;
		    vector<double> basisDerivs_v;
		    double param[2];
		    param[0] = param_u[kr];
		    param[1] = param_v[kj];
		    spline_surf->computeBasis(param, basisValues,
					      basisDerivs_u, basisDerivs_v);

		    // Alternative use struct for storage (see SplineSurface.h)
		    BasisDerivsSf result1;
		    spline_surf->computeBasis(param_u[kr], param_v[kj],
					      result1);

		    // Compute also 2. derivative
		    BasisDerivsSf2 result2;
		    spline_surf->computeBasis(param_u[kr], param_v[kj],
					      result2);
		  }

		// Given the position on the surface, find the corresponding
		// parameter value through a closest point iteration
		double upar, vpar, dist;  // Parameter of closest point and distance
		Point pos2;  // Closest point
		face->closestPoint(pos, upar, vpar, pos2, dist, eps);

		// Alternatively
		surf->closestPoint(pos, upar, vpar, pos2, dist, eps);

		int stop_debug_gridpoint = 1;
	      }
	  }

      // Perform grid evaluation. In the case of a trimmed surface, the
      // grid will extend beyond the trimmed surface
      // Points and derivatives
      vector<double> points;
      vector<double> derivs_u;
      vector<double> derivs_v;
      spline_surf->gridEvaluator(param_u, param_v, points, derivs_u, derivs_v);

      // Basis functions including 1. derivatives
      // See SplineSurface.h for the struct BasisDerivsSf
     vector<BasisDerivsSf> basis_derivs1;
      spline_surf->computeBasisGrid(param_u, param_v, basis_derivs1);

      // Including 2. derivatives
      // See SplineSurface.h for the struct BasisDerivsSf2
      vector<BasisDerivsSf2> basis_derivs2;
      spline_surf->computeBasisGrid(param_u, param_v, basis_derivs2);

      int stop_debug_face = 1;
    }
  
  
}

