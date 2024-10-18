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

#include "GoTools/geometry/SplineDebugUtils.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/Torus.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/Sphere.h"
#include "GoTools/geometry/Cone.h"

#include <fstream>


using std::vector;
using std::endl;
using std::setprecision;

namespace Go
{


//===========================================================================
void SplineDebugUtils::writeSpaceParamCurve(const SplineCurve& pcurve, std::ostream& os, double z)
//===========================================================================
{
    ALWAYS_ERROR_IF(pcurve.dimension() != 2,
		"Expecting input of 2D-curve.");

    std::vector<double> space_coefs;
    for (int i = 0; i < pcurve.numCoefs(); ++i) {
	space_coefs.insert(space_coefs.end(),
			   pcurve.coefs_begin() + i*2,
			   pcurve.coefs_begin() + (i + 1)*2);
	space_coefs.push_back(z); // Make param_curve live in plane parallell to the xy-plane.
    }

    SplineCurve space_pcurve =
	SplineCurve(pcurve.numCoefs(), pcurve.order(),
		      pcurve.basis().begin(), space_coefs.begin(), 3);
    space_pcurve.writeStandardHeader(os);
    space_pcurve.write(os);
}

//===========================================================================
void SplineDebugUtils::writeSpace1DCurve(const SplineCurve& curve, std::ostream& os, double z)
//===========================================================================
{
    ALWAYS_ERROR_IF(curve.dimension() != 1,
		"Expecting input of 1D-curve.");

    std::vector<double> coefs;
    std::vector<double> coefs_in(curve.coefs_begin(), curve.coefs_end());
    for (int i = 0; i < curve.numCoefs(); ++i) {
      double gp = curve.basis().grevilleParameter(i);
      coefs.push_back(gp);
      coefs.push_back(coefs_in[i]);
      coefs.push_back(z);
    }

    SplineCurve space_curve =
	SplineCurve(curve.numCoefs(), curve.order(),
		      curve.basis().begin(), coefs.begin(), 3);
    space_curve.writeStandardHeader(os);
    space_curve.write(os);
}



//===========================================================================
void SplineDebugUtils::writeSpaceParamCurve(const Line& pline, std::ostream& os, double z)
//===========================================================================
{
    ALWAYS_ERROR_IF(pline.dimension() != 2,
		"Expecting input of 2D-curve.");

    Point start_pt = pline.ParamCurve::point(pline.startparam());
    Point start_pt_3d(start_pt[0], start_pt[1], z);

    Point end_pt = pline.ParamCurve::point(pline.endparam());
    Point end_pt_3d(end_pt[0], end_pt[1], z);

    Line space_pline(start_pt_3d, end_pt_3d,
		     pline.startparam(), pline.endparam(),
		     pline.isReversed());

    space_pline.writeStandardHeader(os);
    space_pline.write(os);
}


//===========================================================================
void SplineDebugUtils::writeSpaceParamCurve(const Circle& pcirc, std::ostream& os, double z)
//===========================================================================
{
    ALWAYS_ERROR_IF(pcirc.dimension() != 2,
		"Expecting input of 2D-curve.");

    Point centre = pcirc.getCentre();
    Point vec1 = pcirc.getXAxis();
    Point normal(0.0, 0.0, 1.0);
    Point centre_3d(centre[0], centre[1], z);
    Point vec1_3d(vec1[0], vec1[1], 0.0);
    double radius = pcirc.getRadius();
    Circle circ_3d(radius, centre_3d, normal, vec1_3d);
    circ_3d.setParamBounds(pcirc.startparam(), pcirc.endparam());

    circ_3d.writeStandardHeader(os);
    circ_3d.write(os);
}


//===========================================================================
void SplineDebugUtils::writeSpaceParamCurve(shared_ptr<ParamCurve> pcurve,
					    std::ostream& os, double z)
//===========================================================================
{
  if (pcurve->instanceType() == Class_SplineCurve)
    {
      shared_ptr<SplineCurve> spline_cv =
	dynamic_pointer_cast<SplineCurve, ParamCurve>(pcurve);
      writeSpaceParamCurve(*spline_cv, os, z);
    }
  else if (pcurve->instanceType() == Class_Line)
    {
      shared_ptr<Line> line =
	dynamic_pointer_cast<Line, ParamCurve>(pcurve);
      writeSpaceParamCurve(*line, os, z);
    }
  else if (pcurve->instanceType() == Class_Circle)
    {
      shared_ptr<Circle> circle =
	dynamic_pointer_cast<Circle, ParamCurve>(pcurve);
      writeSpaceParamCurve(*circle, os, z);
    }
}

//===========================================================================
void SplineDebugUtils::writeTrimmedInfo(BoundedSurface& bd_sf,
		      std::ostream& os, double z)
//===========================================================================
{
    shared_ptr<ParamSurface> under_sf = bd_sf.underlyingSurface();
    under_sf->writeStandardHeader(os);
    under_sf->write(os);
    int nmb_loops = bd_sf.numberOfLoops();
    for (int kj = 0; kj < nmb_loops; ++kj) {
	shared_ptr<CurveLoop> loop = bd_sf.loop(kj);
	for (int kk = 0; kk < loop->size(); ++kk) {
	    shared_ptr<CurveOnSurface> cv_on_sf =
		dynamic_pointer_cast<CurveOnSurface, ParamCurve>
		((*loop)[kk]);
	    shared_ptr<ParamCurve> par_cv =
		cv_on_sf->parameterCurve();
	    if (par_cv.get() != NULL) {
		if (par_cv->instanceType() == Class_SplineCurve) {
		    shared_ptr<SplineCurve> spline_cv =
			dynamic_pointer_cast<SplineCurve, ParamCurve>
			(par_cv);
		    writeSpaceParamCurve(*spline_cv,
					 os);
		}
		else if (par_cv->instanceType() == Class_Line) {
		    shared_ptr<Line> line_cv =
			dynamic_pointer_cast<Line, ParamCurve>
			(par_cv);
		    writeSpaceParamCurve(*line_cv,
					 os);
		}
		else
		{
		    MESSAGE("Curve type not supported!");
		}
	    }
	    shared_ptr<ParamCurve> space_cv =
		cv_on_sf->spaceCurve();
	    if (space_cv.get() != NULL) {
		space_cv->writeStandardHeader(os);
		space_cv->write(os);
	    }
	}
    }
}


//===========================================================================
void SplineDebugUtils::writeBoundary(BoundedSurface& bd_sf,
				     std::ostream& os)
//===========================================================================
{
  int nmb_loops = bd_sf.numberOfLoops();
  for (int kj = 0; kj < nmb_loops; ++kj) {
    shared_ptr<CurveLoop> loop = bd_sf.loop(kj);
    for (int kk = 0; kk < loop->size(); ++kk) {
      shared_ptr<CurveOnSurface> cv_on_sf =
	dynamic_pointer_cast<CurveOnSurface, ParamCurve>((*loop)[kk]);
      shared_ptr<ParamCurve> par_cv = cv_on_sf->parameterCurve();
      shared_ptr<ParamCurve> space_cv = cv_on_sf->spaceCurve();
      if (!space_cv.get())
	    cv_on_sf->ensureSpaceCrvExistence(0.1);  // A rough tolerance
      if (space_cv.get() && space_cv->dimension() != 1)
	{
	  space_cv->writeStandardHeader(os);
	  space_cv->write(os);
	}
      else if (par_cv.get())
	{
	  shared_ptr<SplineCurve> par_spline = 
	    dynamic_pointer_cast<SplineCurve,ParamCurve>(par_cv);
	  shared_ptr<SplineCurve> space_spline = 
	    dynamic_pointer_cast<SplineCurve,ParamCurve>(space_cv);

	  if (par_spline.get() && space_spline.get() &&
	      space_spline->dimension() == 1 &&
	      par_spline->numCoefs() == space_spline->numCoefs() &&
	      par_spline->order() == space_spline->order())
	    {
	      std::vector<double> space_coefs;
	      for (int i = 0; i < par_spline->numCoefs(); ++i) {
		space_coefs.insert(space_coefs.end(),
				   par_spline->coefs_begin() + i*2,
				   par_spline->coefs_begin() + (i + 1)*2);
		space_coefs.insert(space_coefs.end(),
				   space_spline->coefs_begin() + i,
				   space_spline->coefs_begin() + (i + 1));
	      }
	      
	      SplineCurve space_pcurve =
		SplineCurve(par_spline->numCoefs(), par_spline->order(),
			    par_spline->basis().begin(), space_coefs.begin(), 3);
	      space_pcurve.writeStandardHeader(os);
	      space_pcurve.write(os);
	    }
	  else if (par_spline.get() && bd_sf.dimension() == 1)
	    {
	      std::vector<double> space_coefs;
	      for (int i = 0; i < par_spline->numCoefs(); ++i) {
		space_coefs.insert(space_coefs.end(),
				   par_spline->coefs_begin() + i*2,
				   par_spline->coefs_begin() + (i + 1)*2);
		double par = par_spline->basis().grevilleParameter(i);
		Point cvpar = par_spline->ParamCurve::point(par);
		Point pos = bd_sf.ParamSurface::point(cvpar[0],cvpar[1]);
		space_coefs.push_back(pos[0]);
	      }
	      
	      SplineCurve space_pcurve =
		SplineCurve(par_spline->numCoefs(), par_spline->order(),
			    par_spline->basis().begin(), space_coefs.begin(), 3);
	      space_pcurve.writeStandardHeader(os);
	      space_pcurve.write(os);
	    }
	  else if (par_spline.get())
	    writeSpaceParamCurve(*par_spline, os, 0.0);
	  else
	    std::cout << "Missing boundary curve" << std::endl;
	}
      else
	std::cout << "Missing boundary curve" << std::endl;
    }
  }
}


//===========================================================================
void SplineDebugUtils::writeOuterBoundaryLoop(ParamSurface& sf,
					      std::ostream& os)
//===========================================================================
{
    // We also write the boundary loops of the underlying surface.
    CurveLoop outer_loop = sf.outerBoundaryLoop();
    for (size_t ki = 0; ki < outer_loop.size(); ++ki)
    {
	outer_loop[ki]->writeStandardHeader(os);
	outer_loop[ki]->write(os);
    }
}


//===========================================================================
void SplineDebugUtils::objToFile(GeomObject* geom_obj, char *to_file)
//===========================================================================
{
    if (geom_obj) {
	std::ofstream debug(to_file);
	geom_obj->writeStandardHeader(debug);
	geom_obj->write(debug);
    }
}


//===========================================================================
void SplineDebugUtils::objsToFile(vector<shared_ptr<GeomObject> >& geom_objs, char *to_file)
//===========================================================================
{
    std::ofstream debug(to_file);
    for (size_t ki = 0; ki < geom_objs.size(); ++ki) {
	if (geom_objs[ki].get() != 0) {
	    geom_objs[ki]->writeStandardHeader(debug);
	    geom_objs[ki]->write(debug);
	}
    }
}

//===========================================================================
void SplineDebugUtils::writeSISLFormat(const SplineCurve& spline_cv, std::ostream& os)
//===========================================================================
{
  int i,j;
  int linenum;

  int dim = spline_cv.dimension();
  int num_coefs = spline_cv.numCoefs();
  int order = spline_cv.order();

  os << setprecision(15);

  os << "$ This is a B-Spline curve" << std::endl;
  os << "$ type: 0 is usual, 5 point, 6 analytic\n" << std::endl;
  os << "0" << std::endl;
  /* order */
  os << "$ order ik" << std::endl;
  os << order << std::endl;

  /* number of control vertices */
  os << "$ number of control vertices in" << std::endl;
  os << num_coefs << std::endl;

  /* dimension of geometry space  */
  os << "$ dimension" << std::endl;
  os << dim << std::endl;

  /* curve open/closed */
  os << "$ curve open/closed" << std::endl;
  os << 1 << std::endl;

  /* nonrational, i.e. polynomial */
  os << "$ rational or not" << std::endl;
  os << spline_cv.rational() << std::endl;

  /* knot vector */
  linenum = (num_coefs + order)/4;
  os << "$ knot vector" << std::endl;
  for (j=0; j < linenum; j++){
      for (i=0; i < 4; i++)
	  os << spline_cv.basis().begin()[j * 4 + i] << " ";
      os << "";
  }
  for (i = linenum * 4; i < (num_coefs + order); i++)
      os << spline_cv.basis().begin()[i] << " ";
  os << std::endl;

  /* control vertices */
  os << "$ control vertices" << std::endl;

  if (!spline_cv.rational())
  {
     for ( i = 0; i < num_coefs; i++ )
     {
	 for ( j = 0; j < dim; j++ )
	     os << spline_cv.coefs_begin()[i*dim+j] << " ";
	 os << 1.0 << std::endl;
     }
  }
  else
  {
     for ( i = 0; i < num_coefs; i++ )
     {
	for ( j = 0; j < dim+1; j++ )
	    os << spline_cv.rcoefs_begin()[i*(dim+1)+j] << " ";
	os << "" << std::endl;
     }
  }

  /* instance matrix */
  os << "$ instance matrix" << std::endl;
  if (dim == 3)
  {
     os << 1.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << std::endl;
     os << 0.0 << " " << 1.0 << " " << 0.0 << " " << 0.0 << std::endl;
     os << 0.0 << " " << 0.0 << " " << 1.0 << " " << 0.0 << std::endl;
     os << 0.0 << " " << 0.0 << " " << 0.0 << " " << 1.0 << std::endl;
  }
  else if (dim == 2)
  {
     os << 1.0 << " " << 0.0 << " " << 0.0 << " " << std::endl;
     os << 0.0 << " " << 1.0 << " " << 0.0 << " " << std::endl;
     os << 0.0 << " " << 0.0 << " " << 1.0 << " " << std::endl;
  }
  else
  {
      os << 1.0 << " " << 0.0 << std::endl;
      os << 0.0 << " " << 1.0 << std::endl;
  }

  os << "$ ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;

  return;
}


void SplineDebugUtils::writeCvsOnSf(const std::vector<shared_ptr<Go::ParamCurve> >& cvs,
				    double epsgeo,
				    std::ofstream& fileout)
{
    vector<shared_ptr<CurveOnSurface> > cvs_on_sf;
    for (size_t ki = 0; ki < cvs.size(); ++ki)
    {
	if (cvs[ki]->instanceType() == Class_CurveOnSurface)
	{
	    cvs_on_sf.push_back(dynamic_pointer_cast<CurveOnSurface>(cvs[ki]));
	}
    }

    writeCvsOnSf(cvs_on_sf, epsgeo, fileout);
}

void SplineDebugUtils::writeCvsOnSf(const std::vector<shared_ptr<Go::CurveOnSurface> >& cvs,
				    double epsgeo,
				    std::ofstream& fileout)
{
    for (size_t ki = 0; ki < cvs.size(); ++ki)
    {
	if (cvs[ki]->instanceType() != Class_CurveOnSurface)
	{
	    continue;
	}
	int prev_ind = (ki + cvs.size() - 1)%(cvs.size());

	shared_ptr<CurveOnSurface> cv_on_sf = cvs[ki];
	shared_ptr<ParamCurve> curr_par_cv = cv_on_sf->parameterCurve();
	if (curr_par_cv)
	{
	    if (curr_par_cv->instanceType() == Class_SplineCurve)
	    {
		SplineDebugUtils::writeSpaceParamCurve(*(dynamic_pointer_cast<SplineCurve>(curr_par_cv)), fileout);
	    }
	    else if (curr_par_cv->instanceType() == Class_Line)
	    {
		SplineDebugUtils::writeSpaceParamCurve(*(dynamic_pointer_cast<Line>(curr_par_cv)), fileout);
	    }
	    shared_ptr<ParamCurve> space_cv = cv_on_sf->spaceCurve();
	    if (space_cv)
	    {
		Point cv_space_pt_start = space_cv->point(space_cv->startparam());
		Point cv_space_pt_end = space_cv->point(space_cv->endparam());
		Point par_pt_start = curr_par_cv->point(curr_par_cv->startparam());
		Point par_pt_end = curr_par_cv->point(curr_par_cv->endparam());
		shared_ptr<ParamSurface> under_sf = cv_on_sf->underlyingSurface();
		Point lifted_par_pt_start = under_sf->point(par_pt_start[0], par_pt_start[1]);
		Point lifted_par_pt_end = under_sf->point(par_pt_end[0], par_pt_end[1]);
		double dist_start = cv_space_pt_start.dist(lifted_par_pt_start);
		double dist_end = cv_space_pt_end.dist(lifted_par_pt_end);
		if (dist_end > epsgeo || dist_start > epsgeo)
		{
		    MESSAGE("Mismatch between end points of par_cv & space_cv. ki = " << ki << ", dist_start: "
			    << dist_start << ", dist_end: " << dist_end);
		}
	    }
	}
	shared_ptr<CurveOnSurface> prev_cv_on_sf = cvs[prev_ind];
	shared_ptr<ParamCurve> prev_par_cv = prev_cv_on_sf->parameterCurve();
	if (prev_par_cv.get() != NULL && curr_par_cv.get() != NULL)
	{
	    Point prev_par_pt = prev_par_cv->point(prev_par_cv->endparam());
	    Point curr_par_pt = curr_par_cv->point(curr_par_cv->startparam());
	    double dist = prev_par_pt.dist(curr_par_pt);
	    double epspar = epsgeo;
	    if (dist > epspar   )
	    {
		MESSAGE("Dist from prev_curve end_pt is outside epspar! dist = " << dist << ", epspar: " << epspar);
	    }
	}
//	prev_par_cv = curr_par_cv;
	if (cv_on_sf.get() != NULL)
	{
	    shared_ptr<ParamCurve> space_cv = cv_on_sf->spaceCurve();
	    if (space_cv.get() != NULL)
	    {
		space_cv->writeStandardHeader(fileout);
		space_cv->write(fileout);
	    }
	}
    }
}

void SplineDebugUtils::writeSeamInfo(Go::BoundedSurface& bd_sf,
                                     std::ofstream& fileout)
{
    // We also write to file the boundary curves of the underlying surface and all end pts for the trim curves.
    shared_ptr<ParamSurface> under_sf = bd_sf.underlyingSurface();
    if (under_sf->instanceType() == Class_Cylinder ||
	under_sf->instanceType() == Class_Sphere ||
        under_sf->instanceType() == Class_Cone ||
	under_sf->instanceType() == Class_Torus)
    {
	shared_ptr<SplineSurface> spline_sf;
	if (under_sf->instanceType() == Class_Cylinder)
	{
	    shared_ptr<Cylinder> cyl = dynamic_pointer_cast<Cylinder>(under_sf);
	    spline_sf = shared_ptr<SplineSurface>(cyl->geometrySurface());
	}
	else if (under_sf->instanceType() == Class_Sphere)
	{
	    shared_ptr<Sphere> sphere = dynamic_pointer_cast<Sphere>(under_sf);
	    spline_sf = shared_ptr<SplineSurface>(sphere->geometrySurface());
	}
	else if (under_sf->instanceType() == Class_Cone)
	{
	    shared_ptr<Cone> cone = dynamic_pointer_cast<Cone>(under_sf);
	    spline_sf = shared_ptr<SplineSurface>(cone->geometrySurface());
	}
	else if (under_sf->instanceType() == Class_Torus)
	{
	    shared_ptr<Torus> torus = dynamic_pointer_cast<Torus>(under_sf);
	    spline_sf = shared_ptr<SplineSurface>(torus->geometrySurface());
	}

	// shared_ptr<Cylinder> cyl = dynamic_pointer_cast<Cylinder>(bd_sf.underlyingSurface());
//	 shared_ptr<SplineSurface> spline_sf(under_sf->asSplineSurface());
	// We then fetch the boundary curves from the spline_sf.
	CurveLoop loop = spline_sf->outerBoundaryLoop();
	for (size_t ki = 0; ki < loop.size(); ++ki)
	{
	    loop[ki]->writeStandardHeader(fileout);
	    loop[ki]->write(fileout);
	}
	// We also fetch all end pts from the trim cvs.
	vector<double> end_pts;
	vector<CurveLoop> bd_loops = bd_sf.absolutelyAllBoundaryLoops();
	bool missing_par_cv = false;
	for (int crv = 0; crv < int(bd_loops.size()); crv++)
	{
	    for (size_t ki = 0; ki < bd_loops[crv].size(); ++ki)
	    {
		shared_ptr<CurveOnSurface> cv_on_sf(dynamic_pointer_cast<
							CurveOnSurface, ParamCurve> (bd_loops[crv][ki]));
		assert(cv_on_sf.get() != 0);
		shared_ptr<ParamCurve> space_cv = cv_on_sf->spaceCurve();
		if (space_cv)
		{
		    space_cv->writeStandardHeader(fileout);
		    space_cv->write(fileout);
		    Point start_pt = space_cv->point(space_cv->startparam());
		    Point end_pt = space_cv->point(space_cv->startparam());
		    end_pts.insert(end_pts.end(), start_pt.begin(), start_pt.end());
		    end_pts.insert(end_pts.end(), end_pt.begin(), end_pt.end());
		}
	    }
	}
	const int dim = spline_sf->dimension();
	int num_pts = end_pts.size()/dim;
	PointCloud3D pt_cl(end_pts.begin(), num_pts);
	pt_cl.writeStandardHeader(fileout);
	pt_cl.write(fileout);
    }
    double break_pt_val = 0.0;
}



} // namespace Go
