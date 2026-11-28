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

#include "GoTools/lrsplines3D/SweepVolumeLR.h"
#include "GoTools/lrsplines3D/Mesh3D.h"
#include "GoTools/lrsplines2D/Mesh2D.h"
#include "GoTools/lrsplines2D/LRBSpline2D.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/BSplineUniUtils.h"
#include "GoTools/geometry/Circle.h"
#include "GoTools/utils/MatrixXD.h"
#include "GoTools/utils/errormacros.h"
#include <iostream>
#include <fstream>

using namespace Go;
using std::vector;

typedef std::map<LRSplineSurface::BSKey, std::unique_ptr<LRBSpline2D> > BSplineMap; // storage of basis functions

namespace {
void defineBsplinesUni(const LRSplineSurface &surface,
		       const SplineCurve& curve, int cvdir,
		       Mesh3D& mesh,
		       vector<std::unique_ptr<BSplineUniLR> >& uni1,
		       vector<std::unique_ptr<BSplineUniLR> >& uni2,
		       vector<std::unique_ptr<BSplineUniLR> >& unicv);
}; // end anonymous namespace 


//===========================================================================
LRSplineVolume* SweepVolumeLR::linearSweptVolume(const LRSplineSurface &surface,
						 const SplineCurve &curve,
						 const Point &pt)
//===========================================================================
{
  int dim = surface.dimension();
  if (curve.dimension() != dim || pt.dimension() != dim)
    THROW("linearSweepLR: dimension mismatch");
  bool cvrat = curve.rational();
  bool rat = cvrat || surface.rational();
  int kdim = (cvrat) ? dim+1 : dim;

  
  // Define trivariate mesh
  int cvdir = 3;
  Mesh2D mesh2d = surface.mesh();
  vector<double> cvknots(curve.knotsBegin(), curve.knotsEnd());
  Mesh3D mesh(mesh2d, cvknots, cvdir);

  // Define univariate B-splines
  vector<std::unique_ptr<BSplineUniLR> > bsplinesuni1, bsplinesuni2, bsplinesuni3;
  defineBsplinesUni(surface, curve, cvdir, mesh,
		    bsplinesuni1, bsplinesuni2, bsplinesuni3);
  
  // Create trivariate B-splines
  int ncoefcv = curve.numCoefs();
  vector<std::unique_ptr<LRBSpline3D> > bsplvol(ncoefcv*surface.numBasisFunctions());
  size_t b_ix = 0;
  for (BSplineMap::const_iterator bspl=surface.basisFunctionsBegin();
       bspl!=surface.basisFunctionsEnd(); ++bspl)
    {
      BSplineUniLR* uni1 = bspl->second->getUnivariate(XFIXED);
      int ix_uni1 = 0;
      bool found = BSplineUniUtils::identify_bsplineuni(uni1, bsplinesuni1, ix_uni1);
      if (!found)
	THROW("linearSweepLR: univariate B-spline not found");
	
      BSplineUniLR* uni2 = bspl->second->getUnivariate(YFIXED);
      int ix_uni2 = 0;
      found = BSplineUniUtils::identify_bsplineuni(uni2, bsplinesuni2, ix_uni2);
      if (!found)
	THROW("linearSweepLR: univariate B-spline not found");

      bool brat = bspl->second->rational();
      Point sfcoef = (brat) ? bspl->second->coefTimesGamma() :
	bspl->second->Coef();
      double gamma = bspl->second->gamma();
      double wgt1 = 1.0;
      if (brat)
	wgt1 = bspl->second->weight();
      
      int ix_uni3 = 0;
      vector<double>::const_iterator it = (cvrat) ? curve.rcoefs_begin() : curve.coefs_begin();
      vector<double>::const_iterator stop = (cvrat) ? curve.rcoefs_end() : curve.coefs_end();
      for (; it!=stop; it+=kdim, ++ix_uni3)
	{
	  Point cvcoef(it, it+dim);
	  double wgt2 = (cvrat) ? it[dim] : 1.0;
	  Point volcoef = (cvcoef/wgt2 + sfcoef - pt);
	  double wgt3 = (rat) ? wgt1*wgt2 : 1.0;
	  if (!brat)
	    volcoef *= gamma;
	  bsplvol[b_ix++] =
	    std::unique_ptr<LRBSpline3D>(new LRBSpline3D(volcoef, wgt3,
							 bsplinesuni1[ix_uni1].get(),
							 bsplinesuni2[ix_uni2].get(),
							 bsplinesuni3[ix_uni3].get(),
							 gamma, rat));
	}
    }
  
  // Define volume
  LRSplineVolume *sweepvol = new LRSplineVolume(mesh, bsplinesuni1, bsplinesuni2,
						bsplinesuni3, bsplvol,
						rat, surface.getKnotTol());
  return sweepvol;
}



//===========================================================================
LRSplineVolume* SweepVolumeLR::rotationalSweptVolume(const LRSplineSurface &surface,
						     double angle,
						     const Point &pt,
						     const Point &axis)
//===========================================================================
{
  ALWAYS_ERROR_IF(surface.dimension() != 3 || pt.dimension() != 3 || axis.dimension() != 3,
		  "Surface, point and axis must lie in three-dimensional space");
  int idim = 3;
  double TWOPI = 2.0*M_PI;
  double PIHALF = 0.5*M_PI;
  
  // First create a full unit circle in the xy-plane of order 3 and 9 control points
  Point centre(0.0, 0.0, 0.0);
  Point normal(0.0, 0.0, 1.0);
  Point x_axis(1.0, 0.0, 0.0);
  double radius = 1.0;
  shared_ptr<Circle> circ(new Circle(radius, centre, normal, x_axis));
  shared_ptr<SplineCurve> full_circle(circ->createSplineCurve());
  
  // Next create the circle segment to rotate along
  double abs_angle = angle >= 0.0 ? angle : -angle;
  if (abs_angle > TWOPI)
    abs_angle = TWOPI;
  double quadrant = floor (abs_angle / PIHALF);

  // Get difference between abs_angle and its closest diagonal line
  double diag_angle = abs_angle - (quadrant + 0.5) * PIHALF;
  // Get end parameter for circle segment
  double end_angle = PIHALF * (quadrant
			       + (1.0 
				  + (sqrt ((double) 2.0) + 1.0)
				    * tan (diag_angle / 2.0))
			         / 2.0);

  // Create circle segment
  shared_ptr<SplineCurve> circle_segment(full_circle->subCurve(0.0, end_angle));

  // Reverse (i.e. mirror across x-axis) if input angle was negative
  if (angle < 0.0)
    {
      vector<double>::iterator rcoefs = circle_segment->rcoefs_begin();
      for (int i=0; i<circle_segment->numCoefs(); ++i)
	rcoefs[i<<2 | 1] = -rcoefs[i<<2 | 1];
    }

  std::ofstream of1("circle_seg.g2");
  circle_segment->writeStandardHeader(of1);
  circle_segment->write(of1);

  // Create volume
  // Define trivariate mesh
  int cvdir = 1;
  Mesh2D mesh2d = surface.mesh();
  vector<double> cvknots(circle_segment->knotsBegin(), circle_segment->knotsEnd());
  Mesh3D mesh(mesh2d, cvknots, cvdir);

  // Define univariate B-splines
  vector<std::unique_ptr<BSplineUniLR> > uni_sf1, uni_sf2, uni_cv;
  defineBsplinesUni(surface, *circle_segment, cvdir,
		    mesh, uni_sf1, uni_sf2, uni_cv);

  // Create trivariate B-splines
  int ncoefcv = circle_segment->numCoefs();
  int ncoefsf = surface.numBasisFunctions();
  
  std::ofstream of2("circle_cfs.g2");
  of2 << "400 1 0 4 255 0 0 255" << std::endl;
  of2 << ncoefcv*ncoefsf << std::endl;
  
  vector<std::unique_ptr<LRBSpline3D> > bsplvol(ncoefcv*ncoefsf);
  const Point axis_norm = axis / axis.length();
  size_t b_ix = 0;
  for (BSplineMap::const_iterator bspl=surface.basisFunctionsBegin();
       bspl!=surface.basisFunctionsEnd(); ++bspl)
    {
      BSplineUniLR* uni1 = bspl->second->getUnivariate(XFIXED);
      int i_uni1 = 0;
      bool found = BSplineUniUtils::identify_bsplineuni(uni1, uni_sf1, i_uni1);
      if (!found)
	THROW("linearSweepLR: univariate B-spline not found");
	
      BSplineUniLR* uni2 = bspl->second->getUnivariate(YFIXED);
      int i_uni2 = 0;
      found = BSplineUniUtils::identify_bsplineuni(uni2, uni_sf2, i_uni2);
      if (!found)
	THROW("linearSweepLR: univariate B-spline not found");

      bool brat = bspl->second->rational();
      Point sfcoef = bspl->second->coefTimesGamma();
       // Point sfcoef = brat ? bspl->second->coefTimesGamma() :
       // 	bspl->second->Coef();
      double gamma = bspl->second->gamma();
      double wgt1 = 1.0;
      if (brat)
	wgt1 = bspl->second->weight();
      
       // Projection of surface coefficient to the axis
      Point surface_proj = pt + ((sfcoef-pt)*axis_norm) * axis_norm;

      // Vector from surface_proj to surface_pt
      Point xaxis = sfcoef - surface_proj;
      double rad = xaxis.length();

      // The vector yaxis will be the unique vector of same length as xaxsis,
      // making xaxis, yaxis, axis_norm an ortonormal right hand basis
      Point yaxis = axis_norm % xaxis;

      // Create the matrix M from xaxis, yaxis, axis_norm, scaled by the height of surface_pt
      // above the axis
      MatrixXD<double, 3> matr;
      for (int i = 0; i < 3; ++i)
	{
	  matr(i,0) = xaxis[i];
	  matr(i,1) = yaxis[i];
	  matr(i,2) = rad * axis_norm[i];
	}
      
     int i_uni3 = 0;
     int cv_kdim = 4;
      vector<double>::const_iterator cv_it = circle_segment->rcoefs_begin();
      vector<double>::const_iterator stop = circle_segment->rcoefs_end();
      for (; cv_it!=stop; cv_it+=cv_kdim, ++i_uni3)
	{
	  vector<double> vol_pt(3);
	  matr.mult(cv_it, vol_pt.begin());

	  // Rescale, as the curve segment coordinates where in a homogeneous system
	  for (int ki=0; ki<idim; ++ki)
	    vol_pt[ki] /= cv_it[idim];

	  of2 << vol_pt[0] << " " << vol_pt[1] << " " << vol_pt[2] << std::endl;

	  // Translate
	  //double gammafac = brat ? gamma : 1.0;
	  for (int ki=0; ki<idim; ++ki)
	    vol_pt[ki] += surface_proj[ki]; //*gammafac;

	  // Volume weight
	  double vol_wgt = cv_it[idim]*wgt1;

	  Point vol_coef(&vol_pt[0], &vol_pt[0]+idim);
	  // if (!brat)
	  //   vol_coef *= gamma;
	  bsplvol[b_ix++] =
	    std::unique_ptr<LRBSpline3D>(new LRBSpline3D(vol_coef, vol_wgt,
							 uni_cv[i_uni3].get(),
							 uni_sf1[i_uni1].get(),
							 uni_sf2[i_uni2].get(),
							 gamma, true));

	  
	}
    }
  
  // Define volume
  LRSplineVolume *sweepvol = new LRSplineVolume(mesh, uni_cv, uni_sf1, uni_sf2,
						bsplvol, true,
						surface.getKnotTol());
  return sweepvol;
}


namespace {
//===========================================================================
void defineBsplinesUni(const LRSplineSurface &surface,
		       const SplineCurve& curve, int cvdir,
		       Mesh3D& mesh,
		       vector<std::unique_ptr<BSplineUniLR> >& uni1,
		       vector<std::unique_ptr<BSplineUniLR> >& uni2,
		       vector<std::unique_ptr<BSplineUniLR> >& unicv)
//===========================================================================
{
  // Define univariate B-splines in the x-direction of the surface
  vector<double> cvknots(curve.knotsBegin(), curve.knotsEnd());
  int sfix1 = (cvdir % 3) + 1;
  int numuni1 = surface.numUnivariateBSplines(XFIXED);
  uni1.resize(numuni1);
  for (int u_ix=0; u_ix < numuni1; ++u_ix)
    {
      BSplineUniLR *curr = surface.getUnivariateBSpline(XFIXED, u_ix);
      uni1[u_ix] =
	std::unique_ptr<BSplineUniLR>(new BSplineUniLR(sfix1, curr->degree(),
						       &curr->kvec()[0],
						       &mesh));
    }
								 
  // Define univariate B-splines in the y-direction of the surface
  int sfix2 = ((cvdir+1) % 3) + 1;
  int numuni2 = surface.numUnivariateBSplines(YFIXED);
  uni2.resize(numuni2);
  for (int v_ix=0; v_ix < numuni2; ++v_ix)
    {
      BSplineUniLR *curr = surface.getUnivariateBSpline(YFIXED, v_ix);
      uni2[v_ix] = 
	std::unique_ptr<BSplineUniLR>(new BSplineUniLR(sfix2, curr->degree(),
						       &curr->kvec()[0],
						       &mesh));
    }
								 
  // Create univariate B-splines in the curve direction
  int ncoefcv = curve.numCoefs();
  int order = curve.order();
  unicv.resize(ncoefcv);
  for (int w_ix = 0; w_ix != ncoefcv; ++w_ix)
    {
      unicv[w_ix] = 
	std::move(std::unique_ptr<BSplineUniLR>(new BSplineUniLR(cvdir, order-1,
								 cvknots.begin() + w_ix,
								 &mesh)));
    }
}
} // end anonymous namespace

