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

#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/CurveModel.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/compositemodel/SurfaceModelUtils.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/creators/CreatorsUtils.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/geometry/Circle.h"
#include "GoTools/geometry/Ellipse.h"
#include "GoTools/geometry/ElementarySurface.h"
#include "GoTools/geometry/Utils.h"
#include "sislP.h"
#include <fstream>

using std::vector;

namespace Go
{



//===========================================================================
// Constructor
//===========================================================================
CompositeModelFactory::CompositeModelFactory(double approxtol,
			double gap,   // Gap between adjacent surfaces
			double neighbour,  // Threshold for whether surfaces are adjacent
			double kink,  // Kink between adjacent surfaces 
			double bend) // Intended G1 discontinuity between adjacent surfaces)
  : approxtol_(approxtol), gap_(gap), neighbour_(neighbour), kink_(kink), bend_(bend)
{
}

//===========================================================================
// Destructor
//===========================================================================
CompositeModelFactory::~CompositeModelFactory()
{
}

//===========================================================================
// Create empty surface model
//===========================================================================
SurfaceModel* CompositeModelFactory::createEmpty()
{
  vector<shared_ptr<ftSurface> > faces;
  SurfaceModel *model = new SurfaceModel(approxtol_, gap_, neighbour_, kink_, bend_, faces);;
  return model;
}

//===========================================================================
// Read IGES file
//===========================================================================
  CompositeModel* CompositeModelFactory::createFromIges(std::istream& is, 
							bool use_filetol,
							bool prefer_surfacemodel)
{
  CompositeModel *model = 0;

  IGESconverter conv;
  try
    {
      conv.readIGES(is);
    }
  catch (...)
    {
      return model;
    }

  return getGeometry(conv, use_filetol, prefer_surfacemodel);
}

//===========================================================================
  vector<shared_ptr<CompositeModel> > 
    CompositeModelFactory::getModelsFromIges(std::istream& is, 
					     bool use_filetol)
{
  vector<shared_ptr<CompositeModel> > models;

  IGESconverter conv;
  try
    {
      conv.readIGES(is);
    }
  catch (...)
    {
      return models;
    }

  // Get all geometric entities
  vector<shared_ptr<ftSurface> > faces;
  vector<shared_ptr<ParamCurve> > curves;
  getAllEntities(conv, faces, curves);

  if (faces.size() > 0)
    {
      CompositeModel *sfmodel = new SurfaceModel(approxtol_, gap_, neighbour_, 
						 kink_, bend_, faces);
      models.push_back(shared_ptr<CompositeModel>(sfmodel));
    }

  if (curves.size() > 0)
    {
      shared_ptr<CurveModel> cvmodel = 
	shared_ptr<CurveModel>(new CurveModel(gap_, neighbour_, 
					      kink_, bend_, curves));
      vector<shared_ptr<CompositeCurve> > comp_crvs = 
	cvmodel->fetchCompositeCurves();
      for (size_t ki=0; ki<comp_crvs.size(); ++ki)
	models.push_back(comp_crvs[ki]);
    }
  
  return models;
}

//===========================================================================
// Read g2 file
//===========================================================================
  CompositeModel* CompositeModelFactory::createFromG2(std::istream& is,
						      bool prefer_surfacemodel)
{
  CompositeModel *model = 0;

  IGESconverter conv;
  try
    {
      conv.readgo(is);
    }
  catch (...)
    {
      return model;
    }

  return getGeometry(conv, false, prefer_surfacemodel);
}


//===========================================================================
  vector<shared_ptr<CompositeModel> > 
    CompositeModelFactory::getModelsFromG2(std::istream& is, 
					     bool use_filetol)
{
  vector<shared_ptr<CompositeModel> > models;

  IGESconverter conv;
  try
    {
      conv.readgo(is);
    }
  catch (...)
    {
      return models;
    }

  // Get all geometric entities
  vector<shared_ptr<ftSurface> > faces;
  vector<shared_ptr<ParamCurve> > curves;
  getAllEntities(conv, faces, curves);

  if (faces.size() > 0)
    {
      CompositeModel *sfmodel = new SurfaceModel(approxtol_, gap_, neighbour_, 
						 kink_, bend_, faces);
      models.push_back(shared_ptr<CompositeModel>(sfmodel));
    }

  if (curves.size() > 0)
    {
      shared_ptr<CurveModel> cvmodel = 
	shared_ptr<CurveModel>(new CurveModel(gap_, neighbour_, 
					      kink_, bend_, curves));
      vector<shared_ptr<CompositeCurve> > comp_crvs = 
	cvmodel->fetchCompositeCurves();
      for (size_t ki=0; ki<comp_crvs.size(); ++ki)
	models.push_back(comp_crvs[ki]);
    }
  
  return models;
}

//===========================================================================
// Read geometry from file converter
//===========================================================================
  CompositeModel* CompositeModelFactory::getGeometry(IGESconverter& conv, 
						     bool use_filetol,
						     bool prefer_surfacemodel)
{
  CompositeModel *model = 0;
  vector<shared_ptr<ftSurface> > faces;
  vector<shared_ptr<ParamCurve> > curves;
  double filetol = gap_;
  if (use_filetol)
    {
      filetol = conv.minResolution();

      // Reset relevant tolerances
      gap_ = filetol;
      neighbour_ = 10.0*gap_;
    }

  getAllEntities(conv, faces, curves);

  // std::ofstream out("iges_sfs.g2");
  // for (size_t ki=0; ki<faces.size(); ++ki)
  //   {
  //     shared_ptr<ParamSurface> tmp = faces[ki]->surface();
  //     tmp->writeStandardHeader(out);
  //     tmp->write(out);
  //   }

  if (faces.size() > curves.size() ||
      (prefer_surfacemodel && faces.size() > 0))
    model = new SurfaceModel(approxtol_, gap_, neighbour_, kink_, bend_, faces);
  else
    model = new CompositeCurve(gap_, neighbour_, kink_, bend_, curves);

  return model;
}

//===========================================================================
// Read all geometry entities from file converter
//===========================================================================
  void 
  CompositeModelFactory::getAllEntities(IGESconverter& conv, 
					vector<shared_ptr<ftSurface> >& faces,
					vector<shared_ptr<ParamCurve> >& curves)
  {
  vector<shared_ptr<GeomObject> > gogeom = conv.getGoGeom();
  int nmbgeom = (int)gogeom.size();
  faces.reserve(nmbgeom); // May be too much, but not really important
  curves.reserve(nmbgeom);
  int face_count = 0;

    // std::ofstream out_file("failure.g2");
  for (int i=0; i<nmbgeom; i++)
    {
      if (gogeom[i].get() == 0)
	continue;
      shared_ptr<GeomObject> lg = gogeom[i];
      if (gogeom[i]->instanceType() == Class_SplineCurve ||
	  gogeom[i]->instanceType() == Class_CurveOnSurface)
	{
	    shared_ptr<ParamCurve> gocv =
	    dynamic_pointer_cast<ParamCurve, GeomObject>(lg);
	    curves.push_back(gocv);
	    if (gogeom[i]->instanceType() == Class_CurveOnSurface)
	      {
		shared_ptr<CurveOnSurface> tmp_cv = 
		  dynamic_pointer_cast<CurveOnSurface, ParamCurve>(gocv);
		replaceElementaryCurves(tmp_cv);
	      }
	}
      else if (gogeom[i]->instanceType() == Class_BoundedCurve)
	{
	  shared_ptr<BoundedCurve> bounded_cv =
	    dynamic_pointer_cast<BoundedCurve, GeomObject>(lg);
	  shared_ptr<ParamCurve> gocv = shared_ptr<ParamCurve>(bounded_cv->geometryCurve());
	  curves.push_back(gocv);
	}
      else if (gogeom[i]->instanceType() >= Class_Line &&
	       gogeom[i]->instanceType() <= Class_Parabola)
	{
	  shared_ptr<ElementaryCurve> elem_cv =
	    dynamic_pointer_cast<ElementaryCurve, GeomObject>(lg);
	  shared_ptr<SplineCurve> gocv = shared_ptr<SplineCurve>(elem_cv->geometryCurve());
	  gocv->setElementaryCurve(elem_cv);
	  curves.push_back(gocv);
	}
      else if (gogeom[i]->instanceType() == Class_SplineSurface)
 	{
	  shared_ptr<SplineSurface> gosf =
	    dynamic_pointer_cast<SplineSurface, GeomObject>(lg);

	  // Reparameterize
	  double usize, vsize;
	  gosf->estimateSfSize(usize, vsize);
	  
	  RectDomain dom = gosf->containingDomain();
	  gosf->setParameterDomain(dom.umin(), dom.umin()+usize,
	  			   dom.vmin(), dom.vmin()+vsize);

	  vector<shared_ptr<ParamSurface> > sfs = 
	    SurfaceModelUtils::checkClosedFaces(gosf, neighbour_);
	  for (size_t kr=0; kr<sfs.size(); ++kr)
	    {
	      shared_ptr<ftSurface> ftsf(new ftSurface(sfs[kr], face_count++));
	      faces.push_back(ftsf);
	    }

	}
      else if (gogeom[i]->instanceType() == Class_BoundedSurface)
	{
	  bool trim_failure = false;
	  shared_ptr<BoundedSurface> gosf =
	    dynamic_pointer_cast<BoundedSurface, GeomObject>(lg);

	  if (gosf->underlyingSurface()->instanceType() >= Class_Plane &&
	      gosf->underlyingSurface()->instanceType() <= Class_Torus)

	    {
	      // Replace elementary surface
	      shared_ptr<ElementarySurface> elem_sf = 
		dynamic_pointer_cast<ElementarySurface, ParamSurface>(gosf->underlyingSurface());

	      // Limit surface
	      RectDomain dom = gosf->containingDomain();
	      RectDomain dom2 = elem_sf->containingDomain();
	      double umin, umax, vmin, vmax;
	      // if (elem_sf->instanceType() == Class_Plane)
	      // 	{
		  umin = dom.umin()-0.1*(dom.umax()-dom.umin());
		  umax = dom.umax()+0.1*(dom.umax()-dom.umin());
		  vmin = dom.vmin()-0.1*(dom.vmax()-dom.vmin());
		  vmax = dom.vmax()+0.1*(dom.vmax()-dom.vmin());
	      // 	}
	      // else
	      // 	{
	      // 	  umin = dom.umin();
	      // 	  umax = dom.umax();
	      // 	  vmin = dom.vmin();
	      // 	  vmax = dom.vmax();
	      // 	}
	      umin = std::max(dom2.umin(), umin);
	      umax = std::min(dom2.umax(), umax);
	      vmin = std::max(dom2.vmin(), vmin);
	      vmax = std::min(dom2.vmax(), vmax);
	      elem_sf->setParameterBounds(umin, vmin, umax, vmax);
	      shared_ptr<SplineSurface> tmp_sf = 
		shared_ptr<SplineSurface>(elem_sf->geometrySurface());
	      tmp_sf->setElementarySurface(elem_sf);

	      vector<CurveLoop> bd_loops = gosf->allBoundaryLoops();
	      vector<vector<shared_ptr<CurveOnSurface> > > tmp_loops;
	      vector<double> space_eps;
	      for (size_t kr=0; kr<bd_loops.size(); ++kr)
		{
		  vector<shared_ptr<CurveOnSurface> > curr_loop;
		  double curr_eps = bd_loops[kr].getSpaceEpsilon();
		  space_eps.push_back(curr_eps);
		  int nmb_cvs = bd_loops[kr].size();
		  for (int kh=0; kh<nmb_cvs; ++kh)
		    {
		      shared_ptr<ParamCurve> tmp_crv = bd_loops[kr][kh];
		      shared_ptr<CurveOnSurface> tmp_sfcv = 
			dynamic_pointer_cast<CurveOnSurface, ParamCurve>(tmp_crv);
		      shared_ptr<CurveOnSurface> new_crv;
		      if (tmp_sfcv.get())
			{
			  replaceElementaryCurves(tmp_sfcv);
			  new_crv = 
			    shared_ptr<CurveOnSurface>(new CurveOnSurface(tmp_sf,
								      tmp_sfcv->parameterCurve(),
								      tmp_sfcv->spaceCurve(),
								      tmp_sfcv->parPref()));
			}
		      else if (tmp_crv->dimension() == tmp_sf->dimension())
			{
			  new_crv = shared_ptr<CurveOnSurface>(new CurveOnSurface(tmp_sf,
									      tmp_crv,
									      false));
			}
		      else
			{
			  new_crv = shared_ptr<CurveOnSurface>(new CurveOnSurface(tmp_sf,
									      tmp_crv,
									      true));
			}
			  
		      curr_loop.push_back(new_crv);
		    }
		  tmp_loops.push_back(curr_loop);
		}		  
	      gosf = shared_ptr<BoundedSurface>
		(new BoundedSurface(tmp_sf, tmp_loops, space_eps));

#ifdef DEBUG
	      int state;
	      bool valid = gosf->isValid(state);
	      if (!valid)
		{
		  std::cout << "Surface nr: " << i << ". Not valid. State:";
		  std::cout << state << std::endl;
		}
#endif
	    }
	  // else
	  //   {
	      // Reparameterize
	      double usize, vsize;
	      gosf->underlyingSurface()->estimateSfSize(usize, vsize);
	      
	      RectDomain dom3 = gosf->underlyingSurface()->containingDomain();
	      gosf->setParameterDomain(dom3.umin(), dom3.umin()+usize,
				       dom3.vmin(), dom3.vmin()+vsize);
	    // }
	  try {
	    CreatorsUtils::fixTrimCurves(gosf, 1.0, gap_, neighbour_, kink_);
	  }
	  catch(...)
	    {
	      // 	      std::cout << "Problems fixing trimming curve" << std::endl;
	      // 	      gosf->writeStandardHeader(out_file);
	      // 	      gosf->write(out_file);
	      trim_failure = true;
	    }
	  // if (trim_failure)
	  //   {
	  //     gosf->writeStandardHeader(out_file);
	  //     gosf->write(out_file);
	  //   }
	  // 	  if (!trim_failure)
	  // 	  {
	  // Test if this improves topology analysis
	  int fix = 0;
	  try {
	    fix = BoundedUtils::checkAndFixLoopOrientation(gosf);
	  }
	  catch(...)
	    {
	      MESSAGE("Problem with boundary loop");
	    }
	  if (fix == 2)
	    std::cout << "Turned boundary loop" << std::endl;
	      
	  vector<shared_ptr<ParamSurface> > sfs = 
	    SurfaceModelUtils::checkClosedFaces(gosf, neighbour_);
	  for (size_t kr=0; kr<sfs.size(); ++kr)
	    {
	      shared_ptr<ftSurface> ftsf(new ftSurface(sfs[kr], face_count++));
	      faces.push_back(ftsf);
	    }
	  //	  }

	}
      else if (gogeom[i]->instanceType() >= Class_Plane &&
	       gogeom[i]->instanceType() <= Class_Torus)
	{
	  shared_ptr<ElementarySurface> elem_sf = 
	    dynamic_pointer_cast<ElementarySurface,GeomObject>(lg);
	  shared_ptr<ParamSurface> gosf = shared_ptr<ParamSurface>(elem_sf->geometrySurface());
	  vector<shared_ptr<ParamSurface> > sfs = 
	    SurfaceModelUtils::checkClosedFaces(gosf, neighbour_);
	  for (size_t kr=0; kr<sfs.size(); ++kr)
	    {
	      shared_ptr<ftSurface> ftsf(new ftSurface(sfs[kr], face_count++));
	      faces.push_back(ftsf);
	    }
	}
    }
  }


// Read a vector of sisl surfaces
//===========================================================================
SurfaceModel* CompositeModelFactory::createFromSisl(vector<SISLSurf*>& surfaces)
//===========================================================================
{
  vector<shared_ptr<ParamSurface> > surf(surfaces.size());
  for (size_t ki=0; ki<surfaces.size(); ki++)
    {
      surf[ki] = shared_ptr<ParamSurface>(SISLSurf2Go(surfaces[ki]));
    }

  SurfaceModel *sfmodel = new SurfaceModel(approxtol_, gap_, neighbour_, kink_,
					     bend_, surf);

  return sfmodel;
}

// Make surface model from a box
// Input is one box corner, the vector from this corner towards another corner in the
// same box side, yet another vector in the plane defining this box side and the lengts
// of the box sides. The first length corresponds to the vector defining the box side,
// the second to the other side in the defined plane and the third defines the depth of
// the box. The lengths may be negative.
SurfaceModel* 
CompositeModelFactory::createFromBox(Point corner, Point side_vec, Point plane_vec,
				     double side1_length, double side2_length, 
				     double side3_length)
{
    SurfaceModel *sfmodel = 0;

    // Define all box corners
    side_vec.normalize();
    Point side3_vec = side_vec.cross(plane_vec);
    side3_vec.normalize();
    Point side2_vec = side3_vec.cross(side_vec);
    Point corner2 = corner + side1_length*side_vec;
    Point corner3 = corner + side2_length*side2_vec;
    Point corner4 = corner3 + side1_length*side_vec;
    Point corner5 = corner + side3_length*side3_vec;
    Point corner6 = corner2 + side3_length*side3_vec;
    Point corner7 = corner3 + side3_length*side3_vec;
    Point corner8 = corner4 + side3_length*side3_vec;

    // Make box sides
    vector<shared_ptr<ParamSurface> > box_sides(6);
    vector<Point> side_corners(4);
    vector<double> knots(4);
    knots[0] = knots[1] = 0.0;
    knots[2] = knots[3] = 1.0;
    
    // Front side
    side_corners[1] = corner;
    side_corners[0] = corner2;
    side_corners[3] = corner3;
    side_corners[2] = corner4;
    box_sides[0] = shared_ptr<ParamSurface>(fromKnotsAndCoefs(2, knots, 2, knots, 
							      side_corners));

    // Left side
    side_corners[1] = corner;
    side_corners[0] = corner3;
    side_corners[3] = corner5;
    side_corners[2] = corner7;
    box_sides[1] = shared_ptr<ParamSurface>(fromKnotsAndCoefs(2, knots, 2, knots, 
							      side_corners));

    // Bottom side
    side_corners[1] = corner;
    side_corners[0] = corner5;
    side_corners[3] = corner2;
    side_corners[2] = corner6;
    box_sides[2] = shared_ptr<ParamSurface>(fromKnotsAndCoefs(2, knots, 2, knots, 
							      side_corners));
    
    // Right side
    side_corners[1] = corner2;
    side_corners[0] = corner6;
    side_corners[3] = corner4;
    side_corners[2] = corner8;
    box_sides[3] = shared_ptr<ParamSurface>(fromKnotsAndCoefs(2, knots, 2, knots, 
							      side_corners));
    
    // Top side
    side_corners[1] = corner3;
    side_corners[0] = corner4;
    side_corners[3] = corner7;
    side_corners[2] = corner8;
    box_sides[4] = shared_ptr<ParamSurface>(fromKnotsAndCoefs(2, knots, 2, knots, 
							      side_corners));

    // Back side
    side_corners[1] = corner5;
    side_corners[0] = corner7;
    side_corners[3] = corner6;
    side_corners[2] = corner8;
    box_sides[5] = shared_ptr<ParamSurface>(fromKnotsAndCoefs(2, knots, 2, knots, 
							      side_corners));
    
    // Composite model
    sfmodel = new SurfaceModel(approxtol_, gap_, neighbour_, kink_, bend_, box_sides);
    return sfmodel;
}


// Make surface model from sphere
// Input is sphere center and radius
//===========================================================================
SurfaceModel* CompositeModelFactory::createFromSphere(Point centre, double radius)
//===========================================================================
{
    Point axis(0.0, 0.0, 1.0);
    Point equator(1.0, 0.0, 0.0);
    axis *= radius;
    equator *= radius;
    return createFromSphere(centre, axis, equator, 2, 4);
}

  // Make surface model from octants of a sphere
  // Input is sphere center, axis toward norht pole, vector from centre to
  // point on equator and the number of octants. The length on the equator vector
  // gives the sphere radius
  // latitude = 1 : Octants in northern hemisphere
  // latitude = 2 : Octants in both hemispheres
  // longitude = 1 : Octants in first quadrant
  // longitude = 2 : Octants in first and second quadrant
  // longitude = 3 : Octants in first, second and third quadrant
  // longitude = 4 : Octants in all quadrants
//===========================================================================
  SurfaceModel* CompositeModelFactory::createFromSphere(Point centre, 
							  Point axis, 
							  Point equator, 
							  int latitude, int longitude)
//===========================================================================
{
    SISLSurf *sphere = 0;
    int stat = 0;
    s1023(centre.begin(), axis.begin(), equator.begin(), latitude,longitude,
	  &sphere, &stat);

    if (stat < 0)
	return 0;

    vector<SISLSurf*> sfs;
    sfs.push_back(sphere);

    SurfaceModel* model = createFromSisl(sfs);
    freeSurf(sphere); 

    return model;
}




// Make surface model from truncated cylinder
//===========================================================================
SurfaceModel* 
CompositeModelFactory::createFromCylinder(Point bottom_pos, Point cylinder_axis,
					  Point major_axis, Point minor_axis)
//===========================================================================
{
    SISLSurf *cylinder = 0;
    int stat = 0;
    
    double ellipse_ratio = major_axis.length()/minor_axis.length();
    double height = cylinder_axis.length();
    cylinder_axis.normalize();

    s1021(bottom_pos.begin(), major_axis.begin(), ellipse_ratio,
	  cylinder_axis.begin(), height, &cylinder, &stat);
    if (stat < 0)
	return 0;

	
    vector<SISLSurf*> sfs;
    sfs.push_back(cylinder);

    SurfaceModel* model = createFromSisl(sfs);
    freeSurf(cylinder); 

    return model;
}


//  Interpolate a set of positional curves. Automatic parameterization
//===========================================================================
SurfaceModel* 
CompositeModelFactory::interpolateCurves(const vector<shared_ptr<SplineCurve> >& curves,
					 vector<double>& parvals, 
					 int open, int degree)
//===========================================================================
{
  // Check if any of the curves are rational
  bool rational = false;
  size_t nmb_cvs = curves.size();
  double startpar = 0.0;
  if (nmb_cvs < 2)
    return 0;

  size_t ki, kj;
  int kdim = curves[0]->numCoefs();

  for (ki=0; ki<nmb_cvs; ++ki)
    if (curves[ki]->rational())
      rational = true;

  int repeatcurve = 0;
  if ((open <= 0) && rational)
    {
	repeatcurve = 1;
    }

  // Convert to sisl curves
  vector<SISLCurve*> sisl_cvs(nmb_cvs+repeatcurve, 0);
  for (ki=0; ki<nmb_cvs; ++ki)
    sisl_cvs[ki] = (rational) ? Curve2SISL_rat(*curves[ki]) : 
      Curve2SISL(*curves[ki], false);
  if (repeatcurve)
    sisl_cvs[ki] = sisl_cvs[0];

  parvals.resize(sisl_cvs.size());

  int stat = 0;
  double *par_out = 0;           // Output parameters
  vector<int> type(nmb_cvs, 1);  // Only positional curves
  SISLSurf *sisl_surf;
  if (rational)
    {
      if (open < 0)
	open = 0;

      // Construct parameterization
      // Rational loft does not offer the possibility of creating a
      // closed curve. 
      // NB! The surface will be closed, but not smooth
      // Compute paremeter interval consequetive curves
      // First put the curves in the same spline space in order
      // to compute distance between coefficients
      double *knot = 0;
      double *scoef = 0;
      int order;
      int nmb_coef;

      s1930((int)sisl_cvs.size(), &sisl_cvs[0], &knot, &scoef, &nmb_coef, 
	    &order, &stat);
      for (ki=0; ki<nmb_cvs; ++ki)
	freeCurve(sisl_cvs[ki]);
      if (stat < 0)
	{
	  if (knot) free(knot);
	  if (scoef) free(scoef);
	  return 0;
	}

      // Remake curves
      for (ki=0; ki<nmb_cvs+1; ++ki)
	sisl_cvs[ki] = newCurve(nmb_coef, order, knot, 
				scoef+(kdim+1)*nmb_coef*ki, 2, kdim, 1);
	    
      if (knot) free(knot);
      if (scoef) free(scoef);

      double *sc1 = sisl_cvs[0]->ecoef;
      double *sc2;
      parvals[0] = startpar;
      for (kj=1; kj<sisl_cvs.size(); ++kj)
	{
	  sc2 = sisl_cvs[kj]->ecoef;
	  double dist = 0.0;
	  for (int ki=0; ki < nmb_coef; ++ki)
	    dist += Utils::distance_squared(sc1, sc1+kdim, sc2);
	  dist /= (double)nmb_coef;
	  startpar += dist;
	  parvals[kj] = startpar;
	  sc1 = sc2;
	      
	}
	      

      s1508((int)nmb_cvs+repeatcurve, &sisl_cvs[0], 
	    &parvals[0], &sisl_surf, &stat);
      if (stat >= 0 && degree > 3)
	{
	  // Raise degree as specified
	  SISLSurf *tmp_surf = 0;
	  int stat2;
	  s1387(sisl_surf, sisl_surf->in1, degree+1, &tmp_surf, &stat2);
	  if (stat2 >= 0)
	    {
	      freeSurf(sisl_surf);
	      sisl_surf = tmp_surf;
	      tmp_surf = 0;
	    }
	}
    }
  else
    {
      s1538((int)nmb_cvs, &sisl_cvs[0], &type[0], startpar, open, 
	    degree+1, 0, &sisl_surf, &par_out, &stat);
      for (ki=0; ki<sisl_cvs.size(); ++ki)
	parvals[ki] = par_out[ki];
    }
    if (stat < 0)
      {
	if (par_out)
	  free(par_out);
	for (ki=0; ki<nmb_cvs; ++ki)
	  freeCurve(sisl_cvs[ki]);
	return 0;
      }

    vector<SISLSurf*> sfs;
    sfs.push_back(sisl_surf);

    SurfaceModel* model = createFromSisl(sfs);
    freeSurf(sisl_surf); 
    for (ki=0; ki<nmb_cvs+repeatcurve; ++ki)
      freeCurve(sisl_cvs[ki]);
    if (par_out)
      free(par_out);

    return model;
  
}

 //  Interpolate a set of positional curves. Parameterization is given.
//===========================================================================
SurfaceModel* 
CompositeModelFactory::interpolateCurves2(const vector<shared_ptr<SplineCurve> >& curves,
					  vector<double>& param, 
					  int open, int degree)
//===========================================================================
{
  // Check if any of the curves are rational
  bool rational = false;
  size_t nmb_cvs = curves.size();
  if (nmb_cvs < 2)
    return 0;

  size_t ki;
  // int kdim = curves[0]->numCoefs();

  for (ki=0; ki<nmb_cvs; ++ki)
    if (curves[ki]->rational())
      rational = true;

  int repeatcurve = 0;
  if (open <= 0 && rational)
    {
	repeatcurve = 1;
    }

  // Convert to sisl curves
  vector<SISLCurve*> sisl_cvs(nmb_cvs+repeatcurve, 0);
  for (ki=0; ki<nmb_cvs; ++ki)
    sisl_cvs[ki] = (rational) ? Curve2SISL_rat(*curves[ki]) : 
      Curve2SISL(*curves[ki], false);
  if (repeatcurve)
    sisl_cvs[ki] = sisl_cvs[0];

  int stat = 0;
  double *par_out = 0;           // Output parameters
  vector<int> type(nmb_cvs, 1);  // Only positional curves
  SISLSurf *sisl_surf;
  if (rational)
    {
      if (open < 0)
	open = 0;

      s1508((int)nmb_cvs+repeatcurve, &sisl_cvs[0], 
	    &param[0], &sisl_surf, &stat);
      if (stat >= 0 && degree > 3)
	{
	  // Raise degree as specified
	  SISLSurf *tmp_surf = 0;
	  int stat2;
	  s1387(sisl_surf, sisl_surf->in1, degree+1, &tmp_surf, &stat2);
	  if (stat2 >= 0)
	    {
	      freeSurf(sisl_surf);
	      sisl_surf = tmp_surf;
	      tmp_surf = 0;
	    }
	}
    }
  else
    {
      s1539((int)nmb_cvs, &sisl_cvs[0], &type[0], &param[0], 0.0, open, 
	    degree+1, 0, &sisl_surf, &par_out, &stat);
    }
    if (stat < 0)
      {
	if (par_out)
	  free(par_out);
	for (ki=0; ki<nmb_cvs; ++ki)
	  freeCurve(sisl_cvs[ki]);
	return 0;
      }

    vector<SISLSurf*> sfs;
    sfs.push_back(sisl_surf);

    SurfaceModel* model = createFromSisl(sfs);
    freeSurf(sisl_surf); 
    for (ki=0; ki<nmb_cvs; ++ki)
      freeCurve(sisl_cvs[ki]);
    if (par_out)
      free(par_out);

    return model;
  
}

//  Interpolate a set of positional curves. Automatic parameterization
//===========================================================================
SurfaceModel* 
CompositeModelFactory::interpolateCurves(const vector<shared_ptr<SplineCurve> >& curves,
					 vector<int>& crv_type,
					 vector<double>& parvals, 
					 int open, int degree)
//===========================================================================
{
  // Check if any of the curves are rational
  bool rational = false;
  size_t nmb_cvs = curves.size();
  double startpar = 0.0;
  if (nmb_cvs < 2)
    return 0;

  size_t ki;
  // int kdim = curves[0]->numCoefs();

  for (ki=0; ki<nmb_cvs; ++ki)
    if (curves[ki]->rational())
      rational = true;

  if (rational)
    return 0;  // Not handled

  // Convert to sisl curves
  vector<SISLCurve*> sisl_cvs(nmb_cvs, 0);
  for (ki=0; ki<nmb_cvs; ++ki)
    sisl_cvs[ki] = Curve2SISL(*curves[ki], false);

  parvals.resize(sisl_cvs.size());

  int stat = 0;
  double *par_out = 0;           // Output parameters
  SISLSurf *sisl_surf;
  s1538((int)nmb_cvs, &sisl_cvs[0], &crv_type[0], startpar, open, 
	degree+1, 0, &sisl_surf, &par_out, &stat);
  if (stat < 0)
    {
      if (par_out)
	free(par_out);
      for (ki=0; ki<nmb_cvs; ++ki)
	freeCurve(sisl_cvs[ki]);
      return 0;
    }

  for (ki=0; ki<sisl_cvs.size(); ++ki)
    parvals[ki] = par_out[ki];

    vector<SISLSurf*> sfs;
    sfs.push_back(sisl_surf);

    SurfaceModel* model = createFromSisl(sfs);
    freeSurf(sisl_surf); 
    for (ki=0; ki<nmb_cvs; ++ki)
      freeCurve(sisl_cvs[ki]);
    if (par_out)
      free(par_out);

    return model;
  
 }

//  Interpolate a set of positional curves. Parameterization is given
//===========================================================================
SurfaceModel* 
CompositeModelFactory::interpolateCurves2(const vector<shared_ptr<SplineCurve> >& curves,
					 vector<int>& crv_type,
					 vector<double>& parvals, 
					 int open, int degree)
//===========================================================================
{
  // Check if any of the curves are rational
  bool rational = false;
  size_t nmb_cvs = curves.size();
  double startpar = 0.0;
  if (nmb_cvs < 2)
    return 0;

  size_t ki;
  // int kdim = curves[0]->numCoefs();

  for (ki=0; ki<nmb_cvs; ++ki)
    if (curves[ki]->rational())
      rational = true;

  if (rational)
    return 0;  // Not handled

  // Convert to sisl curves
  vector<SISLCurve*> sisl_cvs(nmb_cvs, 0);
  for (ki=0; ki<nmb_cvs; ++ki)
    sisl_cvs[ki] = Curve2SISL(*curves[ki], false);


  int stat = 0;
  double *par_out = 0;           // Output parameters
  SISLSurf *sisl_surf;
  s1539((int)nmb_cvs, &sisl_cvs[0], &crv_type[0], &parvals[0], startpar, open, 
	degree+1, 0, &sisl_surf, &par_out, &stat);
  if (stat < 0)
    {
      if (par_out)
	free(par_out);
      for (ki=0; ki<nmb_cvs; ++ki)
	freeCurve(sisl_cvs[ki]);
      return 0;
    }

    vector<SISLSurf*> sfs;
    sfs.push_back(sisl_surf);

    SurfaceModel* model = createFromSisl(sfs);
    freeSurf(sisl_surf); 
    for (ki=0; ki<nmb_cvs; ++ki)
      freeCurve(sisl_cvs[ki]);
    if (par_out)
      free(par_out);

    return model;
  
 }

//===========================================================================
CompositeCurve*
CompositeModelFactory::createLineSegment(Point startpt, Point endpt)
//===========================================================================
{
  SplineCurve *line_seg;
  double ptol = 1.0e-5;

  if (startpt.dist(endpt) < ptol)
    line_seg = new SplineCurve(startpt, 0.0, endpt, 1.0);
  else
    line_seg = new SplineCurve(startpt, endpt);

  vector<shared_ptr<ParamCurve> > curves;
  curves.push_back(shared_ptr<ParamCurve>(line_seg));

  CompositeCurve *cvmodel = new CompositeCurve(gap_, neighbour_, kink_, bend_,
					       curves);

  return cvmodel;
}



//===========================================================================
CompositeCurve*
CompositeModelFactory::createCircularArc(Point centre, Point start_pt, double angle,
					 Point axis)
//===========================================================================
{
  // Define full circle
  Point x_axis = start_pt - centre;
  double radius = (x_axis).length();
  x_axis.normalize();
  Circle circ(radius, centre, axis, x_axis);

  // Restrict circle
  shared_ptr<Circle> sub_circle = shared_ptr<Circle>(circ.subCurve(0.0, angle));

  // Make spline version
  shared_ptr<ParamCurve> circle_seg = shared_ptr<ParamCurve>(sub_circle->geometryCurve());

  vector<shared_ptr<ParamCurve> > curves;
  curves.push_back(shared_ptr<ParamCurve>(circle_seg));

  CompositeCurve *cvmodel = new CompositeCurve(gap_, neighbour_, kink_, bend_,
					       curves);

  return cvmodel;
}

//===========================================================================
CompositeCurve*
CompositeModelFactory::createEllipticArc(Point centre, Point direction, 
					 double r1, double r2,
					 double startpar, double angle,
					 Point axis)
//===========================================================================
{
  // Create full ellips
  Ellipse ell(centre, direction, axis, r1, r2);

  // Restrict the ellipsis
  double par2 = startpar + angle;
  if (par2 > 2.0*M_PI)
    par2 -= 2.0*M_PI;
  shared_ptr<Ellipse> sub_ellipse = 
    shared_ptr<Ellipse>(ell.subCurve(startpar, par2));

  // Make spline version
  shared_ptr<ParamCurve> ellipse_seg = shared_ptr<ParamCurve>(sub_ellipse->geometryCurve());

  vector<shared_ptr<ParamCurve> > curves;
  curves.push_back(shared_ptr<ParamCurve>(ellipse_seg));

  CompositeCurve *cvmodel = new CompositeCurve(gap_, neighbour_, kink_, bend_,
					       curves);

  return cvmodel;
  
}

  // Define full circle
// Make spline surface from control points
//===========================================================================
  SplineSurface* 
  CompositeModelFactory::fromKnotsAndCoefs(int order1, vector<double> knots1, int order2,
					   vector<double> knots2, vector<Point> coefs)
//===========================================================================
{
    int nmbcoef1 = (int)knots1.size() - order1;
    int nmbcoef2 = (int)knots2.size() - order2;

    if (nmbcoef1 < order1 ||
	nmbcoef2 < order2 ||
	coefs.size() == 0 ||
	nmbcoef1*nmbcoef2 != (int)coefs.size())
	return 0;  // Input does not define a surface

    int dim = coefs[0].size();
    vector<double> coefficients(dim*nmbcoef1*nmbcoef2);

    int ki, kr, kh;
    for (ki=0, kh=0; ki<(int)coefs.size(); ki++)
	for (kr=0; kr<dim; kr++, kh++)
	    coefficients[kh] = coefs[ki][kr];

    SplineSurface *surf = new SplineSurface(nmbcoef1, nmbcoef2, order1, order2,
					    &knots1[0], &knots2[0], &coefficients[0],
					    dim);

    return surf;
}

//===========================================================================
void CompositeModelFactory::replaceElementaryCurves(shared_ptr<CurveOnSurface> sf_cv)
//===========================================================================
{
  shared_ptr<ParamCurve> pcv = sf_cv->parameterCurve();
  shared_ptr<ElementaryCurve> ecv;
  double t1 = sf_cv->startparam();
  double t2 = sf_cv->endparam();
  if (pcv.get())
    {
      ecv = dynamic_pointer_cast<ElementaryCurve, ParamCurve>(pcv);
      if (!ecv.get())
	{
	  shared_ptr<BoundedCurve> bcv = dynamic_pointer_cast<BoundedCurve, ParamCurve>(pcv);
	  if (bcv.get())
	    {
	      ecv = dynamic_pointer_cast<ElementaryCurve, ParamCurve>(bcv->underlyingCurve());
	      t1 = bcv->startparam();
	      t2 = bcv->endparam();
	    }
	}
    }
  if (ecv.get())
    {
      ecv->setParamBounds(t1, t2);
      shared_ptr<SplineCurve> scv(ecv->createSplineCurve());
      scv->setElementaryCurve(ecv);
      sf_cv->setParameterCurve(scv);
    }

  shared_ptr<ParamCurve> spacecv = sf_cv->spaceCurve();
  shared_ptr<ElementaryCurve> ecv2;
  t1 = sf_cv->startparam();
  t2 = sf_cv->endparam();
  if (spacecv.get())
    {
      ecv2 = dynamic_pointer_cast<ElementaryCurve, ParamCurve>(spacecv);
      if (!ecv2.get())
	{
	  shared_ptr<BoundedCurve> bcv = dynamic_pointer_cast<BoundedCurve, ParamCurve>(spacecv);
	  if (bcv.get())
	    {
	      ecv2 = dynamic_pointer_cast<ElementaryCurve, ParamCurve>(bcv->underlyingCurve());
	      t1 = bcv->startparam();
	      t2 = bcv->endparam();
	    }
	}
    }
  if (ecv2.get())
    {
      ecv2->setParamBounds(t1, t2);
      shared_ptr<SplineCurve> scv(ecv2->createSplineCurve());
      scv->setElementaryCurve(ecv2);
      sf_cv->setSpaceCurve(scv);
    }
}


} // namespace Go

