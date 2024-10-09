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
#include "GoTools/compositemodel/SurfaceModelUtils.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/compositemodel/CompositeCurve.h"
#include "GoTools/compositemodel/ftPointSet.h"
#include "GoTools/compositemodel/AdaptSurface.h"
#include "GoTools/intersections/Identity.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/geometry/ElementaryUtils.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/geometry/RectDomain.h"
#include "GoTools/geometry/Plane.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/Cone.h"
#include "GoTools/geometry/Sphere.h"
#include "GoTools/geometry/Torus.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/compositemodel/Body.h"
#include "GoTools/topology/FaceConnectivityUtils.h"
#include "GoTools/tesselator/RectangularSurfaceTesselator.h"
#include "GoTools/tesselator/ParametricSurfaceTesselator.h"
#include "GoTools/tesselator/RegularMesh.h"
#include "GoTools/tesselator/GenericTriMesh.h"
#include "GoTools/tesselator/TesselatorUtils.h"
#include "GoTools/creators/CurveCreators.h"
#include "GoTools/geometry/SISLconversion.h"
#include "sislP.h"

#include <fstream>

//#define DEBUG

using std::vector;
using std::set;
using std::pair;
using std::make_pair;
using namespace Go;

//===========================================================================
vector<shared_ptr<ParamSurface> > 
SurfaceModelUtils::checkClosedFaces(shared_ptr<ParamSurface> surface, double tol)
//===========================================================================
{
  vector<shared_ptr<ParamSurface> > sfs;

#ifdef DEBUG
  std::ofstream of("close_sf.g2");
  surface->writeStandardHeader(of);
  surface->write(of);
#endif

  // Fetch non-trimmed surface to test
  shared_ptr<ParamSurface> sf;
  RectDomain dom = surface->containingDomain();
  shared_ptr<BoundedSurface> bd_sf = 
    dynamic_pointer_cast<BoundedSurface, ParamSurface>(surface);
  if (bd_sf.get())
    {
      shared_ptr<ParamSurface> tmp = bd_sf->underlyingSurface();
      vector<shared_ptr<ParamSurface> > sub_sfs;
      try {
	sub_sfs = tmp->subSurfaces(dom.umin(), dom.vmin(), 
				   dom.umax(), dom.vmax());
      }
      catch (...)
	{
	  sf = tmp;
	}

      if (sub_sfs.size() > 1)
	sf = tmp;
      else if (sub_sfs.size() == 1)
	sf = sub_sfs[0];
    }
  else
    sf = surface;
  
  // Fetch opposite boundary curves and check coincidence
  vector<shared_ptr<ParamCurve> > cvs_u1 = sf->constParamCurves(dom.umin(), false);
  vector<shared_ptr<ParamCurve> > cvs_u2 = sf->constParamCurves(dom.umax(), false);
  vector<shared_ptr<ParamCurve> > cvs_v1 = sf->constParamCurves(dom.vmin(), true);
  vector<shared_ptr<ParamCurve> > cvs_v2 = sf->constParamCurves(dom.vmax(), true);

  if ((cvs_u1.size() == 0) || (cvs_u2.size() == 0) || (cvs_v1.size() == 0) || (cvs_v2.size() == 0))
  {
	throw std::runtime_error("Seems like the object is missing the implementation of constParamCurves().");
  }

  // Compare also with middle curve
  vector<shared_ptr<ParamCurve> > cvs_u3 = 
    sf->constParamCurves(0.5*(dom.umin()+dom.umax()), false);
  vector<shared_ptr<ParamCurve> > cvs_v3 = 
    sf->constParamCurves(0.5*(dom.vmin()+dom.vmax()), true);
  
  Identity ident;
  int coinc1 = ident.identicalCvs(cvs_u1[0], cvs_u1[0]->startparam(), cvs_u1[0]->endparam(),
				  cvs_u2[0], cvs_u2[0]->startparam(), cvs_u2[0]->endparam(),
				  tol);
  int coinc2 = ident.identicalCvs(cvs_v1[0], cvs_v1[0]->startparam(), cvs_v1[0]->endparam(),
				  cvs_v2[0], cvs_v2[0]->startparam(), cvs_v2[0]->endparam(),
				  tol);
  int coinc3 = ident.identicalCvs(cvs_u1[0], cvs_u1[0]->startparam(), cvs_u1[0]->endparam(),
				  cvs_u3[0], cvs_u3[0]->startparam(), cvs_u3[0]->endparam(),
				  tol);
  int coinc4 = ident.identicalCvs(cvs_v1[0], cvs_v1[0]->startparam(), cvs_v1[0]->endparam(),
				  cvs_v3[0], cvs_v3[0]->startparam(), cvs_v3[0]->endparam(),
				  tol);  vector<shared_ptr<ParamSurface> > sub_sfs1;
  if (coinc3)
    coinc1 = false;   // More likely to be a sliver face
  if (coinc4)
    coinc2 = false;   // More likely to be a sliver face

  if (coinc1)
    {
      // Split in the first parameter direction
      vector<shared_ptr<ParamSurface> > sub_sfs2;
      double mid = 0.5*(dom.umin()+dom.umax());
      try {
	sub_sfs1 = surface->subSurfaces(dom.umin(), dom.vmin(), 
					mid, dom.vmax());
      }
      catch (...)
	{
	  sfs.push_back(surface);
	  return sfs;
	}
	
      try {
	sub_sfs2 = surface->subSurfaces(mid, dom.vmin(), 
					dom.umax(), dom.vmax());
      }
      catch(...)
	{
	  sfs.push_back(surface);
	  return sfs;
	}
      sub_sfs1.insert(sub_sfs1.end(), sub_sfs2.begin(), sub_sfs2.end());
    }
  else
    sub_sfs1.push_back(surface);
	
  if (coinc2)
    {
      for (size_t ki=0; ki<sub_sfs1.size(); ++ki)
	{
	  // Split in the second parameter direction
	  RectDomain dom2 = sub_sfs1[ki]->containingDomain();
	  vector<shared_ptr<ParamSurface> > sub_sfs2;
	  double mid = 0.5*(dom2.vmin()+dom2.vmax());
	  try {
	    sub_sfs2 = sub_sfs1[ki]->subSurfaces(dom2.umin(), dom2.vmin(), 
						 dom2.umax(), mid);
	  }
	  catch (...)
	    {
	      return sub_sfs1;
	    }
	    
	  sfs.insert(sfs.end(), sub_sfs2.begin(), sub_sfs2.end());
	  sub_sfs2.clear();
	  try {
	    sub_sfs2 = sub_sfs1[ki]->subSurfaces(dom2.umin(), mid,
						 dom2.umax(), dom2.vmax());
	  }
	  catch (...)
	    {
	      return sub_sfs1;
	    }
	  sfs.insert(sfs.end(), sub_sfs2.begin(), sub_sfs2.end());
	}
#ifdef DEBUG
      for (size_t kj=0; kj<sfs.size(); ++kj)
	{
	  sfs[kj]->writeStandardHeader(of);
	  sfs[kj]->write(of);
	}
#endif
      return sfs;
    }
  else
    {
#ifdef DEBUG
      for (size_t kj=0; kj<sub_sfs1.size(); ++kj)
	{
	  sub_sfs1[kj]->writeStandardHeader(of);
	  sub_sfs1[kj]->write(of);
	}
#endif
      return sub_sfs1;
    }
}

//===========================================================================
void
SurfaceModelUtils::sameUnderlyingSurf(vector<shared_ptr<ftSurface> >& sf_set,
				      double tol, double angtol,
				      vector<vector<shared_ptr<ftSurface> > >& faces,
				      vector<shared_ptr<ParamSurface> >& under_sfs)
//===========================================================================
{
  // Extract bounded surfaces
  vector<shared_ptr<ParamSurface> > cand_sfs;
  for (size_t ki=0; ki<sf_set.size(); ++ki)
    {
      shared_ptr<ParamSurface> sf = sf_set[ki]->surface();
      if (sf->instanceType() == Class_BoundedSurface)
	cand_sfs.push_back(sf);
    }

  // Make pairwise check of candidate surfaces for identical underlying surfaces
  for (size_t ki=0; ki<cand_sfs.size(); )
    {
      size_t incr = 1;
      shared_ptr<BoundedSurface> bd_sf1 =
	dynamic_pointer_cast<BoundedSurface, ParamSurface>(cand_sfs[ki]);
      shared_ptr<ParamSurface> under1 = bd_sf1->underlyingSurface();
      for (size_t kj=ki+1; kj<cand_sfs.size(); ++kj)
	{
	  shared_ptr<BoundedSurface> bd_sf2 =
	    dynamic_pointer_cast<BoundedSurface, ParamSurface>(cand_sfs[kj]);
	  shared_ptr<ParamSurface> under2 = bd_sf2->underlyingSurface();
	  bool same = false;
	  if (under1.get() == under2.get() || 
	      ElementaryUtils::sameElementarySurface(under1.get(), 
						     under2.get(), 
						     tol, angtol))
	    // Same underlying surface
	    same = true;
	  else
	    {
	      // Check coincidence
	      // Check first for bounded underlying surfaces
	      ElementarySurface *elem1 = under1->elementarySurface();
	      ElementarySurface *elem2 = under2->elementarySurface();
	      if (elem1 && (!elem1->isBounded()))
		same = false;
	      else if (elem2 && (!elem2->isBounded()))
		same = false;
	      else
		{
		  Identity ident;
		  int coinc = ident.identicalSfs(under1, under2, tol);
		  if (coinc >= 1)
		    same = true;
		}
	    }
	  if (same)
	    {
	      std::swap(cand_sfs[ki+incr], cand_sfs[kj]);
	      incr++;
	    }
	}
      if (incr > 1)
	{
	  // Identical underlying surfaces are found
	  vector<shared_ptr<ftSurface> > curr_faces;
	  for (size_t kj=ki; kj<ki+incr; ++kj)
	    {
	      for (size_t kr=0; kr<sf_set.size(); ++kr)
		{
		  if (sf_set[kr]->surface().get() == cand_sfs[kj].get())
		    {
		      curr_faces.push_back(sf_set[kr]);
		      break;
		    }
		}
	    }
	  shared_ptr<BoundedSurface> bd_sf =
	    dynamic_pointer_cast<BoundedSurface, ParamSurface>(cand_sfs[ki]);
	  if (bd_sf.get())
	    {
	      faces.push_back(curr_faces);
	      shared_ptr<ParamSurface> curr_under = bd_sf->underlyingSurface();
	      ElementarySurface *curr_elem = curr_under->elementarySurface();
	      if (curr_elem)
		{
		  // Construct underlying surface with a possibly extend 
		  // parameter domain covering all associated faces
		  // Case dependent
		  curr_under = extendedUnderlyingSurface(curr_faces, tol, angtol);
		}
	      if (curr_under.get())
		under_sfs.push_back(curr_under);
	      else
		under_sfs.push_back(bd_sf->underlyingSurface());
	    }
	}
      ki += incr;
    }
}


//===========================================================================
shared_ptr<ParamSurface>
SurfaceModelUtils::extendedUnderlyingSurface(vector<shared_ptr<ftSurface> >& sf_set,
					     double tol, double angtol)
//===========================================================================
{
  shared_ptr<ParamSurface> surf;  // Resulting surface

  ElementarySurface *elem1 = sf_set[0]->surface()->elementarySurface();
  if (!elem1)
    return surf; // Function only applicable for elementary surfaces. No result

  // Check that all surfaces is of the same type
  // Check also equality of surface descriptions
  RectDomain dom = elem1->getParameterBounds();
  RectDomain dom1 = dom;
  RectDomain pardom = elem1->containingDomain();
  Point loc1 = elem1->location();
  Point dir1 = elem1->direction();
  for (size_t ki=1; ki<sf_set.size(); ++ki)
    {
      ElementarySurface *elem2 = sf_set[ki]->surface()->elementarySurface();
      if (!elem2 || elem1->instanceType() != elem2->instanceType())
	return surf; 

      RectDomain dom2 = elem2->getParameterBounds();
      Point loc2 = elem2->location();
      Point dir2 = elem2->direction();
      double ang = (dir1.dimension() == 0) ? 0.0 : dir1.angle(dir2);
      if (ang > angtol && M_PI-ang > angtol)
	continue;   // Not the same surface

      bool swapped = ((elem1->isSwapped() && !elem2->isSwapped()) || 
		      (!elem1->isSwapped() && elem2->isSwapped()));
      if (swapped)
	continue;  // Different orientation

      // Case distinction
      if (elem1->instanceType() == Class_Plane)
	{
	  double len1 = fabs((loc1-loc2)*dir1);
	  if (len1 > tol)
	    continue;   // Not the same surface

	  Plane *plane1 = dynamic_cast<Plane*>(elem1);
	  Plane *plane2 = dynamic_cast<Plane*>(elem2);
	  if (plane1 == NULL || plane2 == NULL)
	    return surf;  // Something is wrong

	  double len2 = loc1.dist(loc2);
	  Point axis1_1, axis1_2, axis2_1, axis2_2;
	  plane1->getSpanningVectors(axis1_1, axis1_2);
	  plane2->getSpanningVectors(axis2_2, axis2_2);
	  if (plane1->isSwapped())
	      std::swap(axis1_1, axis1_2);
	  if (plane2->isSwapped())
	      std::swap(axis2_1, axis2_2);
	  double ang2 = axis1_1.angle(axis1_2);
	  if (len2 > tol || ang2 > angtol)
	    {
	      // Parameterization of the two planes differ.
	      // Modify parameter domain of surface 2 to match that of surface 1
	      Vector2D low = dom2.lowerLeft();
	      Vector2D high = dom2.upperRight();
	      if (len2 > tol)
		{
		  // Move domain
		  Point diff = loc2 - loc1;
		  low[0] += diff[0];
		  high[0] += diff[0];
		  low[1] += diff[1];
		  high[1] += diff[1];
		}
	      if (ang2 > angtol)
		{
		  // Extend domain to make sure to cover the rotated domain
		  double fac = 2.0/sqrt(2.0);
		  Vector2D mid = 0.5*(high + low);
		  low = mid - fac*(mid - low);
		  high += high + fac*(high - mid);
		}
	      RectDomain dom3(low, high);
	      dom1.addUnionWith(dom3);
	    }
	  else
	    {
	      // Simply extend the initial domain
	      dom1.addUnionWith(dom2);
	    }
	}
      else if (elem1->instanceType() == Class_Cylinder)
	{
	  double len1 = fabs((loc1-loc2)*dir1);
	  if (len1 > tol)
	    continue;   // Not the same surface

	  Cylinder *cyl1 = dynamic_cast<Cylinder*>(elem1);
	  Cylinder *cyl2 = dynamic_cast<Cylinder*>(elem2);
	  if (cyl1 == NULL || cyl2 == NULL)
	    return surf;  // Something is wrong

	  double rad1 = cyl1->getRadius();
	  double rad2 = cyl2->getRadius();
	  if (fabs(rad2 - rad1) > tol)
	    continue;   // Not the same surface

	  Point axis1_1, axis1_2, axis1_3, axis2_1, axis2_2, axis2_3;;
	  cyl1->getCoordinateAxes(axis1_1, axis1_2, axis1_3);
	  cyl2->getCoordinateAxes(axis2_1, axis2_2, axis2_3);
	  double len2 = loc1.dist(loc2);
	  double ang2 = axis1_1.angle(axis2_1);
	  if (len2 > tol || ang2 > angtol)
	    {
	      // Parameterization of the two cylinders differ.
	      // Modify parameter domain of surface 2 to match that of surface 1
	      Vector2D low = dom2.lowerLeft();
	      Vector2D high = dom2.upperRight();
	      if (len2 > tol)
		{
		  // Move domain
		  Point diff = loc2 - loc1;
		  int sgn = (diff * dir1 > 0.0);
		  low[1] += sgn*len2;
		  high[1] += sgn*len2;
		}
	      if (ang2 > angtol)
		{
		  Point v = axis1_1.cross(axis2_1);
#ifdef DEBUG
		  std::cout << "Reparameterization of cylinder. To be continued" << std::endl;
#endif
		}
	      RectDomain dom3(low, high);
	      dom1.addUnionWith(dom3);
	    }
	  else
	    {
	      dom1.addUnionWith(dom2);
	    }
	}
      else if (elem1->instanceType() == Class_Cone)
	{
	  Cone *cone1 = dynamic_cast<Cone*>(elem1);
	  Cone *cone2 = dynamic_cast<Cone*>(elem2);
	  if (cone1 == NULL || cone2 == NULL)
	    return surf;  // Something is wrong
	  double angle1 = cone1->getConeAngle();
	  double angle2 = cone2->getConeAngle();
	  double ang = dir1.angle(dir2);
	  if (fabs(angle1-angle2) > angtol || 
	      (ang > angtol && M_PI-ang > angtol))
	    continue; // Not the same surface

	  Point axis1_1, axis1_2, axis1_3, axis2_1, axis2_2, axis2_3;;
	  cone1->getCoordinateAxes(axis1_1, axis1_2, axis1_3);
	  cone2->getCoordinateAxes(axis2_1, axis2_2, axis2_3);
	  double len2 = loc1.dist(loc2);
	  double ang2 = axis1_1.angle(axis2_1);
	  if (len2 > tol || ang2 > angtol)
	    {
	      // The parameterization differ
	      // Modify parameter domain of surface 2 to match that of surface 1
	      Vector2D low = dom2.lowerLeft();
	      Vector2D high = dom2.upperRight();
	      if (len2 > tol)
		{
		  // Move domain
		  Point diff = loc2 - loc1;
		  int sgn = (diff * dir1 > 0.0);
		  low[1] += sgn*len2;
		  high[1] += sgn*len2;
		}
	      if (ang2 > angtol)
		{
#ifdef DEBUG
		  std::cout << "Reparameterization of cone. To be continued" << std::endl;
#endif

		}
	      RectDomain dom3(low, high);
	      dom1.addUnionWith(dom3);
	    }
	  else
	    {
	      dom1.addUnionWith(dom2);
	    }
	}
      else if (elem1->instanceType() == Class_Sphere)
	{
#ifdef DEBUG
	  std::cout << "Sphere" << std::endl;
#endif
	}
      else if (elem1->instanceType() == Class_Torus)
	{
#ifdef DEBUG
	  std::cout << "Torus" << std::endl;
#endif
	}
      else
	return surf;  // Surface type not supported
    }

  // Create surface. Case distinction
  double fac1 = (dom1.umax()-dom1.umin())/(dom.umax()-dom.umin());
  double fac2 = (dom1.vmax()-dom1.vmin())/(dom.vmax()-dom.vmin());
  if (elem1->instanceType() == Class_Plane)
    {
      Plane *plane1 = dynamic_cast<Plane*>(elem1);
      Point axis1, axis2;
      plane1->getSpanningVectors(axis1, axis2);
      shared_ptr<Plane> plane2(new Plane(loc1, dir1, axis1, 
					 plane1->isSwapped()));
      if (plane1->isSwapped())
	plane2->setParameterBounds(dom1.vmin(), dom1.umin(), 
				   dom1.vmax(), dom1.umax());
      else
	plane2->setParameterBounds(dom1.umin(), dom1.vmin(), 
				   dom1.umax(), dom1.vmax());
      if (plane1->isSwapped())
	std::swap(fac1, fac2);
      plane2->setParameterDomain(pardom.umin(), 
				 pardom.umin()+fac1*(pardom.umax()-pardom.umin()),
				 pardom.vmin(),
				 pardom.vmin()+fac2*(pardom.vmax()-pardom.vmin()));
      surf = plane2;
    }
  else if (elem1->instanceType() == Class_Cylinder)
    {
      Cylinder *cyl1 = dynamic_cast<Cylinder*>(elem1);
      double rad1 = cyl1->getRadius();
      Point axis1, axis2, axis3;
      cyl1->getCoordinateAxes(axis1, axis2, axis3);
      shared_ptr<Cylinder> cyl2(new Cylinder(rad1, loc1, dir1, axis1, 
					     cyl1->isSwapped()));
      if (cyl1->isSwapped())
	cyl2->setParameterBounds(dom1.vmin(), dom1.umin(), 
				 dom1.vmax(), dom1.umax());
      else
	cyl2->setParameterBounds(dom1.umin(), dom1.vmin(), 
				 dom1.umax(), dom1.vmax());
      if (cyl1->isSwapped())
	std::swap(fac1, fac2);
      cyl2->setParameterDomain(pardom.umin(), 
			       pardom.umin()+fac1*(pardom.umax()-pardom.umin()),
			       pardom.vmin(),
			       pardom.vmin()+fac2*(pardom.vmax()-pardom.vmin()));
      surf = cyl2;
    }
  else if (elem1->instanceType() == Class_Cone)
    {
      Cone *cone1 = dynamic_cast<Cone*>(elem1);
      double rad = cone1->getRadius();
      double cone_angle = cone1->getConeAngle();
      Point axis1, axis2, axis3;
      cone1->getCoordinateAxes(axis1, axis2, axis3);
      shared_ptr<Cone> cone2(new Cone(rad, loc1, dir1, axis1, 
				      cone_angle, cone1->isSwapped()));
      if (cone1->isSwapped())
	cone2->setParameterBounds(dom1.vmin(), dom1.umin(), 
				  dom1.vmax(), dom1.umax());
      else
	cone2->setParameterBounds(dom1.umin(), dom1.vmin(), 
				  dom1.umax(), dom1.vmax());

      if (cone1->isSwapped())
	std::swap(fac1, fac2);
      cone2->setParameterDomain(pardom.umin(), 
				pardom.umin()+fac1*(pardom.umax()-pardom.umin()),
				pardom.vmin(),
				pardom.vmin()+fac2*(pardom.vmax()-pardom.vmin()));
      surf = cone2;
    }

#ifdef DEBUG
  std::ofstream of("elem_faces.g2");
  for (size_t ka=0; ka<sf_set.size(); ++ka)
    {
      sf_set[ka]->surface()->writeStandardHeader(of);
      sf_set[ka]->surface()->write(of);
    }
  surf->writeStandardHeader(of);
  surf->write(of);
#endif
  return surf;
}

//===========================================================================
void SurfaceModelUtils::simplifySurfaceModel(shared_ptr<SurfaceModel>& model,
					    int degree)
//===========================================================================
{
  double eps = model->getTolerances().gap;
  double angtol = model->getTolerances().kink;
  double neighbour = model->getTolerances().neighbour;
  double bend = model->getTolerances().bend;
  int nmb = model->nmbEntities();

  vector<shared_ptr<ftSurface> > merge_master;
  vector<vector<ftSurface*> > merge_slave;
  vector<vector<pair<shared_ptr<Vertex>, shared_ptr<Vertex> > > > merge_vxs;
  for (int ki=0; ki<nmb; ++ki)
    {
      // For all faces, find the candidate neighbours for merging with
      // the current face
      // Fetch all neighbours
      shared_ptr<ftSurface> curr = model->getFace(ki);
      vector<pair<shared_ptr<Vertex>, shared_ptr<Vertex> > > common_vxs;
      vector<ftSurface*> mergecand = getMergeCandFaces(curr, common_vxs,
						       bend /*angtol*/);
      if (mergecand.size() == 0)
	continue;  // No candidate for merge

      merge_master.push_back(curr);
      merge_slave.push_back(mergecand);
      merge_vxs.push_back(common_vxs);
    }
  if (merge_master.size() == 0)
    return;  // No merge is possible

  // Sort master faces for merge with respect to the number of possible
  // merge configurations
  size_t kj, kr;
  for (kj=0; kj<merge_master.size(); ++kj)
    for (kr=kj+1; kr<merge_master.size(); ++kr)
      if (merge_slave[kj].size() < merge_slave[kr].size())
	{
	  std::swap(merge_master[kj], merge_master[kr]);
	  std::swap(merge_slave[kj], merge_slave[kr]);
	  std::swap(merge_vxs[kj], merge_vxs[kr]);
	}
   
  for (size_t ki=0; ki<merge_master.size(); ++ki)
    {
      shared_ptr<ftSurface> curr = merge_master[ki];
      vector<ftSurface*> mergecand = merge_slave[ki];

      // Prioritize according to shape of an eventual merged surface
     vector<double> sf_len_frac(mergecand.size());
     vector<double> other_len_frac(mergecand.size());
     vector<double> sf_reg(mergecand.size());
     for (kj=0; kj<mergecand.size();)
	{
	  estMergedSfSize(curr.get(), mergecand[kj], merge_vxs[ki][kj].first,
			  merge_vxs[ki][kj].second, sf_len_frac[kj], 
			  other_len_frac[kj], sf_reg[kj],
			  neighbour, bend);
	  if (sf_len_frac[kj] < eps)
	    {
	      // Not a proper configuration. Remove
	      mergecand.erase(mergecand.begin()+kj);
	      sf_len_frac.erase(sf_len_frac.begin()+kj);
	      other_len_frac.erase(other_len_frac.begin()+kj);
	      sf_reg.erase(sf_reg.begin()+kj);
	      merge_vxs[ki].erase(merge_vxs[ki].begin()+kj);
	    }
	  else
	    kj++;
	}
      
      // Sort candidates with respect to regularity of merged surface
     for (kj=0; kj<mergecand.size(); ++kj)
	{
	  for (kr=kj+1; kr<mergecand.size(); ++kr)
	    {
	      double frac1 = other_len_frac[kj]/sf_len_frac[kj];
	      double tmp1 = sf_len_frac[kj] - other_len_frac[kj];
	      double frac2 = other_len_frac[kr]/sf_len_frac[kr];
	      double tmp2 = sf_len_frac[kr] - other_len_frac[kr];
#ifdef DEBUG
	      std::ofstream of0("merge_sort.g2");
	      curr->surface()->writeStandardHeader(of0);
	      curr->surface()->write(of0);
	      mergecand[kj]->surface()->writeStandardHeader(of0);
	      mergecand[kj]->surface()->write(of0);
	      mergecand[kr]->surface()->writeStandardHeader(of0);
	      mergecand[kr]->surface()->write(of0);
#endif
	      // @@@ VSK, 190413
	      // This test is tuned towards a specific case (The TERRIFIC
	      // part with blends in the indentations) and it is likely
	      // that it must be refined when we have gained more experience
	      // with similar cases.
	      //if (frac2 < frac1)
	      if ((((tmp2 > tmp1 && tmp2 > 0.0 && frac2 < frac1) ||
		   (frac2 < frac1 && frac2 < 1.0)) &&
		  !(std::max(frac1,frac2) < 1.0 && std::min(tmp1,tmp2) > 0.0 &&
		    sf_reg[kj] > sf_reg[kr])) ||
		  (std::max(frac1,frac2) < 1.0 && std::min(tmp1,tmp2) > 0.0 &&
		   sf_reg[kr] > sf_reg[kj]))
		//		 &&  !(sf_len_frac[kr] < other_len_frac[kj]))
		{
		  std::swap(mergecand[kj], mergecand[kr]);
		  std::swap(sf_len_frac[kj], sf_len_frac[kr]);
		  std::swap(other_len_frac[kj], other_len_frac[kr]);
		  std::swap(sf_reg[kj], sf_reg[kr]);
		  std::swap(merge_vxs[ki][kj], merge_vxs[ki][kr]);
		}
	    }
	}
     merge_slave[ki] = mergecand;
    }

  // Try to merge surfaces
  for (size_t ki=0; ki<merge_master.size(); ++ki)
    {
      shared_ptr<ftSurface> curr = merge_master[ki];
      vector<ftSurface*> mergecand = merge_slave[ki];
      Body *bd = curr->getBody();

     for (kj=0; kj<mergecand.size(); ++kj)
	{
	  int pardir1, pardir2;
	  double parval1, parval2;
	  bool atstart1, atstart2;
	  pair<Point, Point> co_par1, co_par2;
	  shared_ptr<ftSurface> merged_face;
 
#ifdef DEBUG
	  std::ofstream of("merge_cand.g2");
	  curr->surface()->writeStandardHeader(of);
	  curr->surface()->write(of);
	  mergecand[kj]->surface()->writeStandardHeader(of);
	  mergecand[kj]->surface()->write(of);
#endif

	  shared_ptr<Vertex> vx1, vx2;
	  shared_ptr<ftEdge> edge1, edge2;
	  int adj_idx = 0;
	  while (curr->areNeighbours(mergecand[kj], edge1, edge2, adj_idx))
	    {
	      if (adj_idx == 0)
		vx1 = edge1->getVertex(true);
	      vx2 = edge1->getVertex(false);
	      adj_idx++;
	    }
	  int found_merge = mergeSituation(curr.get(), mergecand[kj], vx1, vx2,
					   pardir1, parval1, atstart1, pardir2,
					   parval2, atstart2, co_par1, co_par2,
					   eps);
	  if (!found_merge)
	    continue;  // Try next configuration
	  if (found_merge == 1)
	    {
	      // Try first a combined approximation
	      // Make intermediate surface model
	      vector<shared_ptr<ParamSurface> > sfs(2);
	      sfs[0] = shared_ptr<ParamSurface>(curr->surface()->clone());
	      sfs[1] = shared_ptr<ParamSurface>(mergecand[kj]->surface()->clone());
	      shared_ptr<SurfaceModel> 
		tmp_model(new SurfaceModel(model->getApproximationTol(),
					   model->getTolerances().gap,
					   model->getTolerances().neighbour,
					   model->getTolerances().kink,
					   model->getTolerances().bend,
					   sfs));

	      // Perform approximation
	      double error;
	      shared_ptr<ParamSurface> approx_surf = tmp_model->approxFaceSet(error, degree);
	      if (approx_surf.get())
		{
		  // Replace the two original faces in the model with the new
		  // surface
		  merged_face = shared_ptr<ftSurface>(new ftSurface(approx_surf, -1));
		  (void)merged_face->createInitialEdges(model->getTolerances().gap, 
							model->getTolerances().kink);
		  merged_face->setBody(bd);

		  shared_ptr<ftSurface> face2 = model->fetchAsSharedPtr(mergecand[kj]);
#ifdef DEBUG
		  std::ofstream of1_0("tmp_merged_model0.g2");
		  int nmb5 = model->nmbEntities();
		  for (int kv=0; kv<nmb5; ++kv)
		    {
		      shared_ptr<ParamSurface> tmp_sf = model->getSurface(kv);
		      tmp_sf->writeStandardHeader(of1_0);
		      tmp_sf->write(of1_0);
		    }
		  bool OK1 = model->checkShellTopology();
#endif
		  vector<shared_ptr<Vertex> > tmp_vx1 = curr->vertices();
		  vector<shared_ptr<Vertex> > tmp_vx2 = face2->vertices();
		  model->removeFace(curr);
		  model->removeFace(face2);
#ifdef DEBUG
		  std::ofstream of1_1("tmp_merged_model1.g2");
		  int nmb3 = model->nmbEntities();
		  for (int kv=0; kv<nmb3; ++kv)
		    {
		      shared_ptr<ParamSurface> tmp_sf = model->getSurface(kv);
		      tmp_sf->writeStandardHeader(of1_1);
		      tmp_sf->write(of1_1);
		    }
		  bool OK2 = model->checkShellTopology();
#endif
		  model->append(merged_face);
#ifdef DEBUG
		  bool OK3 = model->checkShellTopology();

		  vector<shared_ptr<Vertex> > tmp_vx = merged_face->vertices();

		  approx_surf->writeStandardHeader(of);
		  approx_surf->write(of);
		  std::ofstream of1_2("tmp_merged_model2.g2");
		  int nmb4 = model->nmbEntities();
		  for (int kv=0; kv<nmb4; ++kv)
		    {
		      shared_ptr<ParamSurface> tmp_sf = model->getSurface(kv);
		      tmp_sf->writeStandardHeader(of1_2);
		      tmp_sf->write(of1_2);
		    }
		  
		  model->setBoundaryCurves();
		  int nmb_bd = model->nmbBoundaries();
		  if (nmb_bd > 0)
		    {
		      std::cout << "Model boundaries exists " << std::endl;
		      vector<shared_ptr<ftEdge> > bd_edges = model->getBoundaryEdges();
		      std::ofstream of2_0("bd_edges.g2");
		      for (size_t b1=0; b1<bd_edges.size(); ++b1)
			{
			  shared_ptr<SplineCurve> tmp_cv(bd_edges[b1]->geomCurve()->geometryCurve());
			  tmp_cv->writeStandardHeader(of2_0);
			  tmp_cv->write(of2_0);
			}
		      int stop_break = 1;
		    }
#endif
		}
	      else
		{
		  // The faces meet with requested continuity, but the common boundary
		  // is not a constant parameter curve in both faces. If possible,
		  // approximate the face with a non-trimmed spline surface
		  int nmb_bd1 = curr->nmbOuterBdCrvs(model->getTolerances().gap,
						     model->getTolerances().neighbour,
						     model->getTolerances().bend); 
		  int nmb_bd2 = mergecand[kj]->nmbOuterBdCrvs(model->getTolerances().gap,
							      model->getTolerances().neighbour,
							      model->getTolerances().bend); 
		  if (nmb_bd1 == 4 && nmb_bd2 == 4)
		    {
		      // Both faces have 4 corners
		      shared_ptr<ftSurface> m1 = 
			model->replaceRegularSurface(curr.get(), true);
		      shared_ptr<ftSurface> m2 = 
			model->replaceRegularSurface(mergecand[kj], true);
		      if (m1.get() || m2.get())
			{
			  // Update candidate information with new surfaces
			  size_t kh, kh2;
			  for (kh=ki+1; kh<merge_master.size(); ++kh)
			    {
			      if (m2.get() && merge_master[kh].get() == mergecand[kj])
				merge_master[kh] = m2;
			      for (kh2=0; kh2<merge_slave[kh].size(); ++kh2)
				{
				  if (m1.get() && merge_slave[kh][kh2] == curr.get())
				    merge_slave[kh][kh2] = m1.get();
				  else if (m2.get() && merge_slave[kh][kh2] == mergecand[kj])
				    merge_slave[kh][kh2] = m2.get();
				}
			    }

			  // Success. Continue merging
			  if (m1.get())
			    curr = m1;
			  if (m2.get())
			    mergecand[kj] = m2.get();
			  shared_ptr<Vertex> vx1, vx2;
			  shared_ptr<ftEdge> edge1, edge2;
			  int adj_idx = 0;
			  while (curr->areNeighbours(mergecand[kj], edge1, edge2, adj_idx))
			    {
			      if (adj_idx == 0)
				vx1 = edge1->getVertex(true);
			      vx2 = edge1->getVertex(false);
			      adj_idx++;
			    }
			  found_merge = 
			    mergeSituation(curr.get(), mergecand[kj], vx1, vx2,
					   pardir1, parval1, atstart1, pardir2,
					   parval2, atstart2, co_par1, co_par2,
					   eps);
#ifdef DEBUG
			  of << "400 1 0 4 255 0 0 255" << std::endl;
			  of << "1" << std::endl;
			  of << vx1->getVertexPoint() << std::endl;
			  of << "400 1 0 4 255 0 0 255" << std::endl;
			  of << "1" << std::endl;
			  of << vx2->getVertexPoint() << std::endl;
			  std::cout << "found_merge = " << found_merge << std::endl;
#endif
			}
		      else
			continue;
		    }
		  else
		    continue;
		}
	    }
	  if (found_merge == 2)
	    {
#ifdef DEBUG
	      curr->surface()->writeStandardHeader(of);
	      curr->surface()->write(of);
	      mergecand[kj]->surface()->writeStandardHeader(of);
	      mergecand[kj]->surface()->write(of);
#endif
	      vector<Point> seam_joints;
	      merged_face = 
		model->mergeFaces(curr.get(), pardir1, parval1, atstart1, 
				  mergecand[kj], pardir2, parval2, atstart2, 
				  co_par1, co_par2, seam_joints);
	    }
	  if (merged_face.get())
	    {
	      // Remove input faces from candidate information
	      size_t kh, kh2;
	      for (kh=ki+1; kh<merge_master.size();)
		{
		  if (merge_master[kh].get() == mergecand[kj])
		    {
		      merge_master.erase(merge_master.begin()+kh);
		      merge_slave.erase(merge_slave.begin()+kh);
		    }
		  else
		    {
		      for (kh2=0; kh2<merge_slave[kh].size();)
			{
			  if (merge_slave[kh][kh2] == curr.get() ||
			      merge_slave[kh][kh2] == mergecand[kj])
			    {
			      merge_slave[kh].erase(merge_slave[kh].begin()+kh2);
			    }
			  else
			    kh2++;
			}
		      kh++;
		    }
		}
	      break;
	    }
	}
    }
#ifdef DEBUG
  std::ofstream of2("merged_model.g2");
  int nmb2 = model->nmbEntities();
  for (int ki=0; ki<nmb2; ++ki)
    {
      shared_ptr<ParamSurface> tmp_sf = model->getSurface(ki);
      tmp_sf->writeStandardHeader(of2);
      tmp_sf->write(of2);
    }

  vector<shared_ptr<Vertex> > vertices;
  model->getAllVertices(vertices);
  of2 << "400 1 0 4 255 0 0 255" << std::endl;
  of2 << vertices.size() << std::endl;
  for (size_t kv=0; kv<vertices.size(); ++kv)
    of2 << vertices[kv]->getVertexPoint() << std::endl;
  of2 << std::endl;
#endif
}

//===========================================================================
void SurfaceModelUtils::simplifySurfaceModel2(shared_ptr<SurfaceModel>& model,
					      int degree, bool remove_joints)
//===========================================================================
{
  Body *bd = model->getBody();

  // Collect edges across which the continuity is g1
  double bend = model->getTolerances().bend;
  FaceConnectivityUtils<ftEdgeBase,ftSurface> connectivity;
  vector<ftEdgeBase*> vec;
  vector<shared_ptr<ftSurface> > faces = model->allFaces();
  connectivity.smoothEdges(faces, vec, bend);
  
#ifdef DEBUG
      if (vec.size() > 0)
	{
	  std::ofstream edgof0("sortedvec0.g2");
	  for (size_t ki=0; ki<vec.size(); ++ki)
	    {
	      shared_ptr<ParamCurve> cv = vec[ki]->geomEdge()->geomCurve();
	      shared_ptr<CurveOnSurface> sf_cv =
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv);
	      if (sf_cv.get())
		cv = sf_cv->spaceCurve();
	      cv->writeStandardHeader(edgof0);
	      cv->write(edgof0);
	    }
	}
#endif

  // Group connected edges
  vector<int> grp_ix;
  int ki, kj, kr, kh;
  grp_ix.push_back(0);
  for (ki=0, kr=0; ki<(int)vec.size(); kr=(int)grp_ix.size()-1, ki=grp_ix[kr])
    {
      grp_ix.push_back(ki+1);
      for (kj=ki+1; kj<(int)vec.size(); ++kj)
	{
	  ftEdge *edg2 = vec[kj]->geomEdge();
	  for (kh=grp_ix[kr]; kh<grp_ix[kr+1]; ++kh)
	    {
	      ftEdge* edg1 = vec[kh]->geomEdge();
	      if (edg1->commonVertex(edg2))
		{
		  int ka = grp_ix[grp_ix.size()-1];
		  if (kj > ka)
		    std::swap(vec[ka], vec[kj]);
		  grp_ix[grp_ix.size()-1]++;
		  break;
		}
	    }
	}
    }
  if (grp_ix[grp_ix.size()-1] < (int)vec.size())
    grp_ix.push_back((int)vec.size());

#ifdef DEBUG
      if (vec.size() > 0)
	{
	  std::ofstream edgof("sortedvec.g2");
	  for (size_t ki=0; ki<vec.size(); ++ki)
	    {
	      shared_ptr<ParamCurve> cv = vec[ki]->geomEdge()->geomCurve();
	      shared_ptr<CurveOnSurface> sf_cv =
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv);
	      if (sf_cv.get())
		cv = sf_cv->spaceCurve();
	      cv->writeStandardHeader(edgof);
	      cv->write(edgof);
	    }
	}
#endif

  // For each edge group, assemble adjacent faces
  tpTolerances toptol = model->getTolerances();
  double tol = model->getApproximationTol();
  vector<set<shared_ptr<ParamSurface> > > all_sfs;
  vector<set<shared_ptr<ftSurface> > > all_faces;
  for (ki=1; ki<(int)grp_ix.size(); ++ki)
    {
      set<shared_ptr<ParamSurface> > curr_sfs;
      set<shared_ptr<ftSurface> > curr_faces;
      DirectionCone union_cone;
      for (kj=grp_ix[ki-1]; kj<grp_ix[ki]; ++kj)
	{
	  vector<ftSurface*> adj_faces = vec[kj]->geomEdge()->getAdjacentFaces();


	  for (size_t kr=0; kr<adj_faces.size(); ++kr)
	    {
	      int ix = model->getIndex(adj_faces[kr]);
	      shared_ptr<ParamSurface> adj = model->getSurface(ix);
	      DirectionCone cone = adj->normalCone();
	      if (union_cone.dimension() == 0)
		union_cone = cone;
	      else
		union_cone.addUnionWith(cone);
	      curr_faces.insert(model->getFace(ix));
	      curr_sfs.insert(adj);
	    }
	}
      if (!union_cone.greaterThanPi())
	{
	  all_sfs.push_back(curr_sfs);
	  all_faces.push_back(curr_faces);
	}
    }

  // Check if the face groups overlap
  for (ki=0; ki<(int)all_faces.size(); /*++ki*/)
    {
      bool overlap = false;
      for (auto it=all_faces[ki].begin(); it!=all_faces[ki].end(); ++it)
      	{
      	  for (kj=ki+1; kj<(int)all_faces.size(); )
      	    {
      	      if (all_faces[kj].find(*it) != all_faces[kj].end())
      		{
      		  all_faces[ki].insert(all_faces[kj].begin(), all_faces[kj].end());
      		  all_faces.erase(all_faces.begin()+kj);
      		  all_sfs[ki].insert(all_sfs[kj].begin(), all_sfs[kj].end());
      		  all_sfs.erase(all_sfs.begin()+kj);
		  overlap = true;
      		}
      	      else
      		++kj;
      	    }
	}
      if (!overlap)
	++ki;
      if (ki == (int)all_faces.size()-1)
	break;
    }
  
  for (ki=0; ki<(int)all_sfs.size(); ++ki)
    {
      vector<shared_ptr<ParamSurface> > grp_sfs(all_sfs[ki].begin(), 
						all_sfs[ki].end());
      vector<shared_ptr<ftSurface> > grp_faces(all_faces[ki].begin(), 
					       all_faces[ki].end());

      shared_ptr<SurfaceModel> grp_model(new SurfaceModel(tol, toptol.gap, 
							  toptol.neighbour,
							  toptol.kink, toptol.bend,
							  grp_sfs));

#ifdef DEBUG
      std::ofstream of("grp_sfs.g2");
      for (size_t kr=0; kr<grp_sfs.size(); ++kr)
	{
	  grp_sfs[kr]->writeStandardHeader(of);
	  grp_sfs[kr]->write(of);
	}
#endif
      // Check for sharp edges
      vector<ftEdge*> corners;
      grp_model->getCorners(corners);
      if (corners.size() == 0)
	{
	  // Continue with the merge

	  double dist;
	  shared_ptr<ParamSurface> approx_surf;
	  
	  // First check if an append situation is feasible
	  approx_surf = grp_model->mergeAllSurfaces();

	  if (!approx_surf.get())
	    {
	      try {
		approx_surf = grp_model->representAsOneSurface(dist);
	      }
	      catch (...)
		{
		  if (approx_surf.get())
		    approx_surf.reset();
		}
	    }
	  if (approx_surf.get())
	    {
	      // Replace surface
	      shared_ptr<ftSurface> new_face(new ftSurface(approx_surf,-1));
	      (void)new_face->createInitialEdges(toptol.gap, toptol.kink);
	      for (size_t kr=0; kr<grp_faces.size(); ++kr)
		model->removeFace(grp_faces[kr]);

	      if (remove_joints)
		{
		  // Try to simplify edge structure of adjacent faces
		  model->setBoundaryCurves();  // Get updated information
		  vector<shared_ptr<ftEdge> > bd_edgs = model->getBoundaryEdges();
		  vector<ftSurface*> bd_faces;
		  for (size_t kr=0; kr<bd_edgs.size(); ++kr)
		    {
		      ftSurface *curr_face = bd_edgs[kr]->face()->asFtSurface();
		      size_t kh=0;
		      if (curr_face)
			{
			  for (kh=0; kh<bd_faces.size(); ++kh)
			    if (bd_faces[kh] == curr_face)
			      break;
			}
		      if (curr_face && kh==bd_faces.size())
			bd_faces.push_back(curr_face);
		    }

		  for (size_t kr=0; kr<bd_faces.size(); ++kr)
		    bd_faces[kr]->joinFreeEdges();
		  model->setBoundaryCurves();  // Get updated information
		}
	      new_face->setBody(bd);
	      model->append(new_face);
	    }
	}
    }
  int stop_break = 1;
}

//===========================================================================
vector<ftSurface*> 
SurfaceModelUtils::getMergeCandFaces(shared_ptr<ftSurface> curr,
				     vector<pair<shared_ptr<Vertex>,
				     shared_ptr<Vertex> > >& common_vxs,
				     double angtol)
//===========================================================================
{
  vector<ftSurface*> neighbours;

  if (curr->twin())
    return neighbours;  // Not a merge candidate

  curr->getAdjacentFaces(neighbours);

  // Keep only those neighbours having sufficient continuity towards
  // the current face (face normals and edges)
  int kj;
  bool radial_edges = curr->hasRealRadialEdges();
 for (kj=0; kj<(int)neighbours.size(); ++kj)
    {
#ifdef DEBUG
      std::ofstream of("adj_sfs.g2");
      curr->surface()->writeStandardHeader(of);
      curr->surface()->write(of);
      neighbours[kj]->surface()->writeStandardHeader(of);
      neighbours[kj]->surface()->write(of);
#endif

      // Both faces cannot be a part of a volumetric adjacency situation
      // (including being adjacent to a twin face situation). Thus, both
      // faces cannot have radial edges
      if (neighbours[kj]->twin() || 
	  (radial_edges && neighbours[kj]->hasRealRadialEdges()))
	{
	  neighbours.erase(neighbours.begin()+kj);
	  kj--;
	  continue;
	}
      int status = 0;
      int adj_idx = 0;
      shared_ptr<ftEdge> edge1, edge2;
      vector<shared_ptr<ftEdge> > edg;
      while (curr->areNeighbours(neighbours[kj], edge1, edge2, adj_idx))
	{
	  int status2 = 0;
	  if (!edge1->hasConnectivityInfo())
	    status2 = 4;
	  else
	    status2 = edge1->getConnectivityInfo()->WorstStatus();
	  status = std::max(status, status2);
	  edg.push_back(edge1);
	  adj_idx++;
	}
      if (status >= 2)
	{
	  // Not a smooth transition
	  neighbours.erase(neighbours.begin()+kj);
	  kj--;
	}
      else
	{
	  // Check continuity of edges
	  // First fetch the vertices bounding the common edges
	  int kr;
	  vector<shared_ptr<Vertex> > vxs(2);
	  if (edg.size() == 1)
	    {
	      vxs[0] = edg[0]->getVertex(true);
	      vxs[1] = edg[0]->getVertex(false);
	    }
	  else
	    {
	      // Check that the edges represents a continuous chain
	      for (kr=1; kr<(int)edg.size(); ++kr)
		if (!edg[kr-1]->commonVertex(edg[kr].get()))
		  break;

	      if (kr < (int)edg.size())
		{
		  // Not a smooth transition
		  neighbours.erase(neighbours.begin()+kj);
		  kj--;
		  continue;
		}

	      shared_ptr<Vertex> tmp_vx = edg[0]->getCommonVertex(edg[1].get());
	      vxs[0] = edg[0]->getOtherVertex(tmp_vx.get());
	      tmp_vx = edg[edg.size()-1]->getCommonVertex(edg[edg.size()-2].get());
	      vxs[1] = edg[edg.size()-1]->getOtherVertex(tmp_vx.get());
	    }
	  
#ifdef DEBUG
	  of << "400 1 0 4 255 0 0 255 " << std::endl;
	  of << "1" << std::endl;
	  of << vxs[0]->getVertexPoint() << std::endl;
	  of << "400 1 0 4 255 0 0 255 " << std::endl;
	  of << "1" << std::endl;
	  of << vxs[1]->getVertexPoint() << std::endl;
#endif

	  for (kr=0; kr<2; ++kr)
	    {
	      vector<ftEdge*> edgf1 = vxs[kr]->getFaceEdges(curr.get());
	      for (size_t kf=0; kf<edgf1.size(); )
		{
		  if (edgf1[kf]->twin() && 
		      edgf1[kf]->twin()->geomEdge()->face() == neighbours[kj])
		    edgf1.erase(edgf1.begin()+kf);
		  else
		    kf++;
		}

	      vector<ftEdge*> edgf2 = vxs[kr]->getFaceEdges(neighbours[kj]);
	      for (size_t kf=0; kf<edgf2.size(); )
		{
		  if (edgf2[kf]->twin() && 
		      edgf2[kf]->twin()->geomEdge()->face() == curr.get())
		    edgf2.erase(edgf2.begin()+kf);
		  else
		    kf++;
		}

	      if (edgf1.size() != 1 || edgf2.size() != 1)
		{
		  neighbours.erase(neighbours.begin()+kj);
		  kj--;
		  break;
		}

	      // Check tangency
	      double t1 = edgf1[0]->parAtVertex(vxs[kr].get());
	      double t2 = edgf2[0]->parAtVertex(vxs[kr].get());
	      Point tan1 = edgf1[0]->tangent(t1);
	      Point tan2 = edgf2[0]->tangent(t2);
	      if (tan1.angle(tan2) > angtol)
		{
		  neighbours.erase(neighbours.begin()+kj);
		  kj--;
		  break;
		}
	    }
	  if (kr == 2)
	    common_vxs.push_back(make_pair(vxs[0], vxs[1]));
	}
    }

  return neighbours;
}

//===========================================================================
void SurfaceModelUtils::estMergedSfSize(ftSurface* face1, ftSurface* face2,
					shared_ptr<Vertex> vx1, 
					shared_ptr<Vertex> vx2,
					double& len_frac, double& other_frac,
					double& sf_reg, double neighbour,
					double bend)
//===========================================================================
{
  len_frac = other_frac = 0.0;

  // Group edges into smooth sequences
  vector<vector<shared_ptr<ftEdgeBase> > > edgegroup1, edgegroup2;
  face1->getBoundaryLoop(0)->groupSmoothEdges(neighbour,
					      bend,
					      edgegroup1);
  face2->getBoundaryLoop(0)->groupSmoothEdges(neighbour,
					      bend,
					      edgegroup2);

  if (edgegroup1.size() != 4 || edgegroup2.size() != 4)
    return;

  // Compute boundary lengths estimates
  size_t ki, kj;
  vector<double> bd_len1(edgegroup1.size(), 0.0);
  vector<double> bd_len2(edgegroup2.size(), 0.0);
  for (ki=0; ki<edgegroup1.size(); ++ki)
    for (kj=0; kj<edgegroup1[ki].size(); ++kj)
      bd_len1[ki] += edgegroup1[ki][kj]->geomEdge()->estimatedCurveLength();

  for (ki=0; ki<edgegroup2.size(); ++ki)
    for (kj=0; kj<edgegroup2[ki].size(); ++kj)
      bd_len2[ki] += edgegroup2[ki][kj]->geomEdge()->estimatedCurveLength();

  // Find the edge group indices corresponding to the edges limited by
  // the vertices vx1 and vx2
  int idx1, idx2;
  for (idx1=0; idx1<(int)edgegroup1.size(); ++idx1)
    {
      shared_ptr<Vertex> v1 = edgegroup1[idx1][0]->geomEdge()->getVertex(true);
      shared_ptr<Vertex> v2 = 
	edgegroup1[idx1][edgegroup1[idx1].size()-1]->geomEdge()->getVertex(false);
      if ((v1.get() == vx1.get() && v2.get() == vx2.get()) ||
	  (v1.get() == vx2.get() && v2.get() == vx1.get()))
	break;
    }
  if (idx1 == (int)edgegroup1.size())
    return;

  for (idx2=0; idx2<(int)edgegroup2.size(); ++idx2)
    {
      shared_ptr<Vertex> v1 = edgegroup2[idx2][0]->geomEdge()->getVertex(true);
      shared_ptr<Vertex> v2 = 
	edgegroup2[idx2][edgegroup2[idx2].size()-1]->geomEdge()->getVertex(false);
      if ((v1.get() == vx1.get() && v2.get() == vx2.get()) ||
	  (v1.get() == vx2.get() && v2.get() == vx1.get()))
	break;
    }
  if (idx2 == (int)edgegroup2.size())
    return;

  double size1, size2, size3, size4;
  size1 = 0.5*(bd_len1[idx1] + bd_len1[(idx1+2)%4]);
  size2 = 0.5*(bd_len2[idx2] + bd_len2[(idx2+2)%4]);
  size3 = 0.5*(bd_len1[(idx1+1)%4] + bd_len1[(idx1+3)%4]);
  size4 = 0.5*(bd_len2[(idx2+1)%4] + bd_len2[(idx2+3)%4]);

  len_frac = 
    std::min(0.5*(size1+size2), size3+size4)/std::max(0.5*(size1+size2), size3+size4);
  other_frac = std::min(size2,size4)/std::max(size2,size4);
  sf_reg = std::min(bd_len1[(idx1+2)%4],bd_len2[(idx2+2)%4])/
    std::max(bd_len1[(idx1+2)%4],bd_len2[(idx2+2)%4]); 
}

//===========================================================================
int SurfaceModelUtils::mergeSituation(ftSurface* face1, ftSurface* face2,
				      shared_ptr<Vertex> vx1, 
				      shared_ptr<Vertex> vx2,
				      int& dir1, double& val1, bool& atstart1, 
				      int& dir2, double& val2, bool& atstart2, 
				      pair<Point, Point>& co_par1, 
				      pair<Point, Point>& co_par2,
				      double eps)
//===========================================================================
{
  // Compute information required to specify a merge operation
  // Traverse edges between faces
  double ptol = 1.0e-8;
  vector<ftEdge*> edges = vx1->uniqueEdges();
  size_t ki;
  dir1 = dir2 = -1;
  val1 = val2 = 0.0;
  for (ki=0; ki<edges.size(); ++ki)
    {
      vector<ftSurface*> adj_faces = edges[ki]->getAdjacentFaces();
      if (adj_faces.size() != 2)
	continue;
      if ((adj_faces[0] == face1 && adj_faces[1] == face2) ||
	  (adj_faces[0] == face2 && adj_faces[1] == face1))
	break;
    }
  if (ki == edges.size())
    return 0;   // Not neighbours

  shared_ptr<Vertex> vx;
  shared_ptr<Vertex> vx0 = vx1;
  ftEdge* edg = edges[ki];
  if (edg->face() == face2)
    edg = edg->twin()->geomEdge();
  while (true)
    {
      shared_ptr<ParamCurve> cv1 = edg->geomCurve();
      shared_ptr<ParamCurve> cv2 = edg->twin()->geomEdge()->geomCurve();
      shared_ptr<CurveOnSurface> sf_cv1 = 
	dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv1);
      shared_ptr<CurveOnSurface> sf_cv2 = 
	dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv2);

      int pardir1, pardir2;
      double parval1, parval2;
       if (!sf_cv1.get() || 
	  !sf_cv1->isConstantCurve(ptol, pardir1, parval1))
	return 1;
      if (!sf_cv2.get() || 
	  !sf_cv2->isConstantCurve(ptol, pardir2, parval2))
	return 1;

      pardir1--;
      pardir2--;
      if (dir1 < 0)
	{
	  dir1 = pardir1;
	  val1 = parval1;
	  dir2 = pardir2;
	  val2 = parval2;
	}
      else if ((pardir1 != dir1 && pardir2 != dir1) || 
	       (pardir2 != dir2 && pardir1 != dir2))
	return 1;   // Not a consistent parameter direction troughout the faces
      // Skip testing of consistent parameter values for the time being. It should be
      // done, but it is a very special situation if the edge sequence is smooth, it is a
      // constant parameter curve and the parameter value in the underlying surface
      // changes.
     
      vx = edg->getOtherVertex(vx0.get());
      if (vx.get() == vx2.get())
	break;

      vector<ftEdge*> edges2 = vx->uniqueEdges();
      size_t kj;
      for (kj=0; kj<edges2.size(); ++kj)
	{
	  vector<ftSurface*> adj_faces = edges2[kj]->getAdjacentFaces();
	  if (adj_faces.size() != 2)
	    continue;
	  if (edges2[kj] == edg || edges2[kj] == edg->twin())
	    continue;
	  if ((adj_faces[0] == face1 && adj_faces[1] == face2) ||
	      (adj_faces[0] == face2 && adj_faces[1] == face1))
	    break;
	}
      if (kj == edges2.size())
	return 0;   // Not neighbours

      edg = edges2[kj];
      vx0 = vx;
    }
      
  Point par1_1 = vx1->getFacePar(face1);
  Point par1_2 = vx1->getFacePar(face2);
  Point par2_1 = vx2->getFacePar(face1);
  Point par2_2 = vx2->getFacePar(face2);
  co_par1 = make_pair(par1_1, par1_2);
  co_par2 = make_pair(par2_1, par2_2);

  // Check whether the parameter at which to merge is at start- or end of the surfaces  
  double u1, u2, v1, v2;
  face1->surface()->getInternalPoint(u1, v1);
  face2->surface()->getInternalPoint(u2, v2);
  atstart1 = (dir1 == 0) ? (u1 > val1) : (v1 > val1);
  atstart2 = (dir2 == 0) ? (u2 > val2) : (v2 > val2);

  return 2;  // Merge information computed, merge across a constant parameter curve
}

//===========================================================================
void 
SurfaceModelUtils::sortTrimmedSurfaces(vector<vector<shared_ptr<CurveOnSurface> > >& cvs1,
				       vector<shared_ptr<ParamSurface> >& sfs1,
				       vector<int>& at_bd1,
				       Body *model1,
				       vector<vector<shared_ptr<CurveOnSurface> > >& cvs2,
				       vector<shared_ptr<ParamSurface> >& sfs2,
				       vector<int>& at_bd2,
				       Body *model2, double eps, double angtol,
				       vector<vector<pair<shared_ptr<ParamSurface>, int> > >& groups,
				       SurfaceModel *shell1, 
				       SurfaceModel *shell2)
//===========================================================================
{
  if ((model1 == NULL && shell1 == NULL) ||
      (model2 == NULL && shell2 == NULL))
    return;   // Not possible to sort trimmed surfaces

  // Make trimmed surfaces and sort trimmed and non-trimmed surface according
  // to whether they are inside or outside the other surface model
  // Sequence: from first group and inside model2, from first group and outside model2,
  // from second group and inside model1, from second group and outside model1
  groups.resize(4);
  for (size_t ki=0; ki<cvs1.size(); ki++)
    {
      if (cvs1[ki].size() == 0)
	{
	  // The surface is not involved in any intersections. Check if
	  // it lies inside or outside the other surface model
	  // Fetch a point in the surface
	  double u, v;
	  Point pnt = sfs1[ki]->getInternalPoint(u,v);

#ifdef DEBUG
	  int state;
	  shared_ptr<BoundedSurface> bd_sf = 
	    dynamic_pointer_cast<BoundedSurface,ParamSurface>(sfs1[ki]);
	  if (bd_sf.get())
	    {
	      bd_sf->analyzeLoops();
	      bool valid = bd_sf->isValid(state);
	      // if (!valid)
	      // 	std::cout << "Surface not valid: " << state << std::endl;
	      std::ofstream of1("curr1.g2");
	      bd_sf->writeStandardHeader(of1);
	      bd_sf->write(of1);
	    }
#endif

	  double pt_dist, ang=0.0;
 	  bool inside = (model2 != NULL) ? 
	    model2->isInside(pnt, pt_dist, ang) : 
	    shell2->isInside(pnt, pt_dist);
	  shared_ptr<ParamSurface> tmp_surf = 
	    shared_ptr<ParamSurface>(sfs1[ki]->clone());
	  if (inside)
	    {
	      groups[0].push_back(make_pair(tmp_surf, ki));
	      // 17102017 A missing symmetry regarding surfaces
	      // coincident with the input model boundaries needs to
	      // be resolved or verified
	      if (fabs(pt_dist) < eps && M_PI-ang < angtol)
		groups[1].push_back(make_pair(shared_ptr<ParamSurface>(tmp_surf->clone()), ki));
	    }
	  else
	    {
	      groups[1].push_back(make_pair(tmp_surf, ki));
	      if (fabs(pt_dist) < eps && ang < angtol)
		groups[0].push_back(make_pair(shared_ptr<ParamSurface>(tmp_surf->clone()), ki));
	    }
	}
      else
	{
	  vector<shared_ptr<BoundedSurface> > trim_sfs;
	  shared_ptr<BoundedSurface> bd_sf1 = 
	    dynamic_pointer_cast<BoundedSurface,ParamSurface>(sfs1[ki]);
	  if (bd_sf1.get())
	    {
	      // First check size of underlying domain, large domains may
	      // make the intersections unstable
	      reduceUnderlyingSurface(bd_sf1, cvs1[ki]);

	      // Make trimmed surfaces
	      try {
		trim_sfs = 
		  BoundedUtils::splitWithTrimSegments(bd_sf1, cvs1[ki], eps);
	      }
	      catch(...)
		{
#ifdef DEBUG
		  std::cout << "Trimmed surfaces missing" << std::endl;
#endif
		}
	    }
	  for (size_t kr=0; kr<trim_sfs.size(); ++kr)
	    {
#ifdef DEBUG
	      int state;
	      trim_sfs[kr]->analyzeLoops();
 	      bool valid = trim_sfs[kr]->isValid(state);
// 	      if (!valid)
// 		std::cout << "Surface not valid: " << state << std::endl;
#endif

	      // // First check size of underlying domain
	      // vector<shared_ptr<CurveOnSurface> > dummy_vec;
	      // reduceUnderlyingSurface(trim_sfs[kr], dummy_vec);

	      // Check if the trimmed surface lies inside or outside the 
	      // other surface model.
	      double u, v;
	      Point pnt =  trim_sfs[kr]->getInternalPoint(u,v);

#ifdef DEBUG
	      std::ofstream of1("curr1.g2");
	      trim_sfs[kr]->writeStandardHeader(of1);
	      trim_sfs[kr]->write(of1);
#endif

	      double pt_dist, ang=0.0;
	      bool inside = (model2 != NULL) ? 
		model2->isInside(pnt, pt_dist, ang) : 
		shell2->isInside(pnt, pt_dist);
	      // double tol2 = model2->getTolerances().neighbour;
	      // if (pt_dist < tol2 /*eps*/)
	      // 	{
	      // 	  // Coincidence. Make extra test
	      // 	  double clo_par[2];
	      // 	  double clo_dist;
	      // 	  Point clo_pt;
	      // 	  int clo_ix;
	      // 	  if (model2 != NULL)
	      // 	    model2->getShell(0)->closestPoint(pnt, clo_pt, clo_ix, clo_par, clo_dist);
	      // 	  else
	      // 	    shell2->closestPoint(pnt, clo_pt, clo_ix, clo_par, clo_dist);
	      // 	  if (clo_ix != at_bd1[ki] && clo_dist < tol2 && clo_dist > eps/*eps*/)
	      // 	    {
	      // 	      shared_ptr<ftSurface> clo_face = (model2 != NULL) ?
	      // 		model2->getShell(0)->getFace(clo_ix) : shell2->getFace(clo_ix);
	      // 	      Point norm = clo_face->normal(clo_par[0], clo_par[1]);
	      // 	      if (norm.length() > eps)
	      // 		norm.normalize();
	      // 	      if ((pnt - clo_pt)*norm > 0.0)
	      // 		inside = false;
	      // 	      // Point pnt2 = pnt + 2.0*tol2*norm;
	      // 	      // double pt_dist2, ang2;
	      // 	      // inside = (model2 != NULL) ? model2->isInside(pnt2, pt_dist2, ang2) : 
	      // 	      // 	shell2->isInside(pnt2, pt_dist2);
	      // 	      int stbr = 1;
	      // 	    }
		      
	      	// }
	      if (inside)
		{
		  groups[0].push_back(make_pair(trim_sfs[kr], ki));
		  // 17102017 A missing symmetry regarding surfaces
		  // coindident with the input model boundaries needs to
		  // be resolved or verified
		  // if (fabs(pt_dist) < eps && M_PI-ang < angtol)
		  //   groups[1].push_back(shared_ptr<ParamSurface>(trim_sfs[kr]->clone()));
		}
			      
	      else
		{
		  groups[1].push_back(make_pair(trim_sfs[kr], ki));
		  if (at_bd1[ki] >= 0 && fabs(pt_dist) < eps && ang < angtol)
		    groups[0].push_back(make_pair(shared_ptr<ParamSurface>(trim_sfs[kr]->clone()), ki));
		}
	    }
	}
    }
  
  for (size_t ki=0; ki<cvs2.size(); ki++)
    {
     if (cvs2[ki].size() == 0)
	{
	  // The surface is not involved in any intersections. Check if
	  // it lies inside or outside the other surface model
	  // Fetch a point in the surface
	  double u, v;
	  Point pnt = sfs2[ki]->getInternalPoint(u,v);

#ifdef DEBUG
	  int state;
	  shared_ptr<BoundedSurface> bd_sf = 
	    dynamic_pointer_cast<BoundedSurface,ParamSurface>(sfs2[ki]);
	  if (bd_sf.get())
	    {
	      bd_sf->analyzeLoops();
	      bool valid = bd_sf->isValid(state);
	      // if (!valid)
	      // 	std::cout << "Surface not valid: " << state << std::endl;

	      std::ofstream of1("curr2.g2");
	      bd_sf->writeStandardHeader(of1);
	      bd_sf->write(of1);
	    }
#endif

	  double pt_dist, ang=0.0;
 	  bool inside = (model1 != NULL) ? 
	    model1->isInside(pnt, pt_dist, ang) : 
	    shell1->isInside(pnt, pt_dist);
	  shared_ptr<ParamSurface> tmp_surf = 
	    shared_ptr<ParamSurface>(sfs2[ki]->clone());
	  if (inside)
	    {
	      groups[2].push_back(make_pair(tmp_surf, ki));
	      // 17102017 A missing symmetry regarding surfaces
	      // coindident with the input model boundaries needs to
	      // be resolved or verified
	      if (fabs(pt_dist) < eps && M_PI-ang < angtol)
		groups[3].push_back(make_pair(shared_ptr<ParamSurface>(tmp_surf->clone()), ki));
	    }
	  else
	    {
	      groups[3].push_back(make_pair(tmp_surf, ki));
	      if (fabs(pt_dist) < eps && ang < angtol)
		groups[2].push_back(make_pair(shared_ptr<ParamSurface>(tmp_surf->clone()), ki));
	    }
	}
      else
	{
	  // Make bounded surfaces
	  vector<shared_ptr<BoundedSurface> > trim_sfs;
	  shared_ptr<BoundedSurface> bd_sf2 = 
	    dynamic_pointer_cast<BoundedSurface,ParamSurface>(sfs2[ki]);
	  if (bd_sf2.get())
	    {
	      // First check size of underlying domain, large domains may
	      // make the intersections unstable
	      reduceUnderlyingSurface(bd_sf2, cvs2[ki]);

	      try {
		trim_sfs = 
		  BoundedUtils::splitWithTrimSegments(bd_sf2, cvs2[ki], eps);
	      }
	      catch(...)
		{
#ifdef DEBUG
		  std::cout << "Trimmed surfaces missing" << std::endl;
#endif
		}
	    }
	  for (size_t kr=0; kr<trim_sfs.size(); ++kr)
	    {
#ifdef DEBUG
	      int state;
	      trim_sfs[kr]->analyzeLoops();
// 	      bool valid = trim_sfs[kr]->isValid(state);
// 	      if (!valid)
// 		std::cout << "Surface not valid: " << state << std::endl;
#endif

	      // // First check size of underlying domain
	      // vector<shared_ptr<CurveOnSurface> > dummy_vec;
	      // reduceUnderlyingSurface(trim_sfs[kr], dummy_vec);

	  // Check if the trimmed surface lies inside or outside the 
	  // other surface model.
	      double u, v;
	      Point pnt =  trim_sfs[kr]->getInternalPoint(u,v);

#ifdef DEBUG
	      std::ofstream of1("curr2.g2");
	      trim_sfs[kr]->writeStandardHeader(of1);
	      trim_sfs[kr]->write(of1);
#endif

	      double pt_dist, ang=0.0;
	      bool inside = (model1 != NULL) ? 
		model1->isInside(pnt, pt_dist, ang) : 
		shell1->isInside(pnt, pt_dist);
	      // double tol1 = model1->getTolerances().neighbour;
	      // if (pt_dist < tol1 /*eps*/)
	      // 	{
	      // 	  // Coincidence. Make extra test
	      // 	  double clo_par[2];
	      // 	  double clo_dist;
	      // 	  Point clo_pt;
	      // 	  int clo_ix;
	      // 	  if (model1 != NULL)
	      // 	    model1->getShell(0)->closestPoint(pnt, clo_pt, clo_ix, clo_par, clo_dist);
	      // 	  else
	      // 	    shell1->closestPoint(pnt, clo_pt, clo_ix, clo_par, clo_dist);
	      // 	  if (clo_dist < tol1 && clo_dist > eps/*eps*/)
	      // 	    {
	      // 	      shared_ptr<ftSurface> clo_face = (model1 != NULL) ?
	      // 		model1->getShell(0)->getFace(clo_ix) : shell2->getFace(clo_ix);
	      // 	      Point norm = clo_face->normal(clo_par[0], clo_par[1]);
	      // 	      if (norm.length() > eps)
	      // 		norm.normalize();
	      // 	      if ((pnt - clo_pt)*norm > 0.0)
	      // 		inside = false;
	      // 	      // Point pnt2 = pnt + 2.0*tol1*norm;
	      // 	      // double pt_dist2, ang2;
	      // 	      // inside = (model1 != NULL) ? model2->isInside(pnt2, pt_dist2, ang2) : 
	      // 	      // 	shell1->isInside(pnt2, pt_dist2);
	      // 	      int stpr = 1;
	      // 	    }
		      
	      	// }
	      if (inside)
		{
		  groups[2].push_back(make_pair(trim_sfs[kr], ki));
		  // 17102017 A missing symmetry regarding surfaces
		  // coindident with the input model boundaries needs to
		  // be resolved or verified
		  // if (fabs(pt_dist) < eps && M_PI-ang < angtol)
		  //   groups[3].push_back(shared_ptr<ParamSurface>(trim_sfs[kr]->clone(), ki));
		}			      
	      else
		{
		  groups[3].push_back(make_pair(trim_sfs[kr], ki));
		  // 17102017 A missing symmetry regarding surfaces
		  // coindident with the input model boundaries needs to
		  // be resolved or verified
		  // if (fabs(pt_dist) < eps && ang < angtol)
		  //   groups[2].push_back(shared_ptr<ParamSurface>(trim_sfs[kr]->clone(),ki));
		}
	    }
	}
    }
#ifdef DEBUG
  std::ofstream of1("inside1.g2");
  std::ofstream of2("outside1.g2");
  std::ofstream of3("inside2.g2");
  std::ofstream of4("outside2.g2");
  for (size_t ki=0; ki<groups[0].size(); ++ki)
    {
      groups[0][ki].first->writeStandardHeader(of1);
      groups[0][ki].first->write(of1);
    }
  for (size_t ki=0; ki<groups[1].size(); ++ki)
    {
      groups[1][ki].first->writeStandardHeader(of2);
      groups[1][ki].first->write(of2);
    }
  for (size_t ki=0; ki<groups[2].size(); ++ki)
    {
      groups[2][ki].first->writeStandardHeader(of3);
      groups[2][ki].first->write(of3);
    }
  for (size_t ki=0; ki<groups[3].size(); ++ki)
    {
      groups[3][ki].first->writeStandardHeader(of4);
      groups[3][ki].first->write(of4);
    }
#endif
}

//===========================================================================
void SurfaceModelUtils::intersectLine(shared_ptr<ParamSurface>& surface,
				      Point pnt, Point dir, double tol,
				      vector<pair<Point,Point> >& result,
				      vector<pair<shared_ptr<ParamCurve>, 
				      shared_ptr<ParamCurve> > >& line_seg)
//===========================================================================
{
  // Convert the surface to a SISLSurf in order to use SISL functions
  // on it. The "false" argument dictates that the SISLSurf will only    
  // copy pointers to arrays, not the arrays themselves.
  shared_ptr<SplineSurface> tmp_spline;
  SplineSurface* splinesf = surface->getSplineSurface();
  if (!splinesf)
    {
      // Convert to spline surface
      tmp_spline = shared_ptr<SplineSurface>(surface->asSplineSurface());
      splinesf = tmp_spline.get();
    }
  const CurveBoundedDomain* bdomain = 0;
  shared_ptr<BoundedSurface> bsurf = 
    dynamic_pointer_cast<BoundedSurface, ParamSurface>(surface);
    if (bsurf.get())
      bdomain = &(bsurf->parameterDomain());

  ASSERT(splinesf != 0);
  
  SISLSurf* sislsf = GoSurf2SISL(*splinesf, false);
  int dim = 3;
  double epsco = 1e-15; // Not used
  double epsge = 1e-6;
  int numintpt;  // number of single intersection points
  double* pointpar = 0; // array containing the parameter values of single intersect. pt.
  int numintcr; // number of intersection curves
  SISLIntcurve** intcurves = 0;
  int stat;

  // Find the intersection points
  s1856(sislsf, pnt.begin(), dir.begin(), dim, epsco, epsge,
	&numintpt, &pointpar, &numintcr, &intcurves, &stat);
  MESSAGE_IF(stat!=0, "s1856 returned code: " << stat);

  int i;
  for (i = 0; i < numintpt; ++i)
    {
      double u = pointpar [i<<1];
      double v = pointpar [i<<1 | 1];
      Point pt = surface->point(u, v);

      bool in_domain = true;
      if (bdomain != 0)
	{
	  // Check if the point is inside the trimmed surface
	  Array<double,2> tmp_pt(u,v);
	  in_domain = bdomain->isInDomain(tmp_pt, epsge);
	  
	}
      if (in_domain)
	result.push_back(std::make_pair(pt, Point(u, v)));
    }

  for (i=0; i<numintcr; i++)
  {
      // Evaluate endpoints of line segment and make geometry curve
      int npt = intcurves[i]->ipoint;
      Point pt1 = surface->point(intcurves[i]->epar1[0],intcurves[i]->epar1[1]);
      Point pt2 = surface->point(intcurves[i]->epar1[2*(npt-1)],
				 intcurves[i]->epar1[2*(npt-1)+1]);
      SplineCurve *gcv = new SplineCurve(pt1, pt2);

      // Project the curve into the parameter space of the surface
      shared_ptr<Point> pt1_2D = shared_ptr<Point>(new Point(intcurves[i]->epar1[0],intcurves[i]->epar1[1]));
      shared_ptr<Point> pt2_2D = shared_ptr<Point>(new Point(intcurves[i]->epar1[2*(npt-1)],intcurves[i]->epar1[2*(npt-1)+1]));
      shared_ptr<ParamCurve> gcv2 = shared_ptr<ParamCurve>(gcv->clone());
      SplineCurve *pcv = CurveCreators::projectSpaceCurve(gcv2, surface, 
							  pt1_2D, pt2_2D, tol);
      
 	vector<SplineCurve*> final_param_curves;
	vector<SplineCurve*> final_space_curves;
	if (bdomain != 0) {
	    // the surface was trimmed.  We must check for intersections with
	    // trimming curves
	    vector<double> params_start_end;
	    bdomain->findPcurveInsideSegments(*pcv, 
					      tol,
					      params_start_end);
	    int num_segments = (int)params_start_end.size() / 2;
	    //cout << "Num segments found: " << num_segments << endl;

	    for (int j = 0; j < num_segments; ++j) {
		SplineCurve* pcv_sub = pcv->subCurve(params_start_end[2 * j],
						     params_start_end[2 * j + 1]);
		final_param_curves.push_back(pcv_sub->clone());
		final_space_curves.push_back(0); //@ change this? (not necessary)
		delete pcv_sub;
	    }
	    // deleting curves that will not be directly used later
	    delete(gcv);
	    delete(pcv);
	} else {
	    final_param_curves.push_back(pcv);
	    final_space_curves.push_back(gcv);
	}

	// pushing back segments
	for (size_t j = 0; j < final_space_curves.size(); ++j) {
	  line_seg.push_back(std::make_pair(shared_ptr<ParamCurve>(final_param_curves[j]),
					    shared_ptr<ParamCurve>(final_space_curves[j])));
	}
  }

      free(pointpar);
      freeIntcrvlist(intcurves, numintcr);
      freeSurf(sislsf);
}

//===========================================================================
bool
SurfaceModelUtils::extremalPoint(shared_ptr<ParamSurface>& surface,
				 Point dir, tpTolerances& toptol,
				 Point& ext_pnt, double ext_par[])
//===========================================================================
{
  bool modified = false;
  double tol2d = 1.0e-4;

  dir.normalize();
#ifdef DEBUG
  std::ofstream of1("sf_ext.g2");
  surface->writeStandardHeader(of1);
  surface->write(of1);
#endif

  // First check bounding box
  BoundingBox box = surface->boundingBox();
  Point vec = box.high() - box.low();
  vec = dir*(vec*dir);
  if (vec.length() < toptol.gap)
    {
      Point tmp_pt = 0.5*(box.low() + box.high());
      double upar, vpar, dist;
      Point clo_pt;
      surface->closestPoint(tmp_pt, upar, vpar, clo_pt, dist, toptol.gap);
      ext_pnt = clo_pt;
      ext_par[0] = upar;
      ext_par[1] = vpar;
      return true;
    }

  // Convert the surface to a SISLSurf in order to use SISL functions
  // on it. The "false" argument dictates that the SISLSurf will only    
  // copy pointers to arrays, not the arrays themselves.
  shared_ptr<SplineSurface> tmp_spline;
  SplineSurface* surf = surface->getSplineSurface();
  if (!surf)
    {
      // Convert to spline surface
      tmp_spline = shared_ptr<SplineSurface>(surface->asSplineSurface());
      surf = tmp_spline.get();
    }
    const CurveBoundedDomain* bddomain = 0;
    shared_ptr<BoundedSurface> bsurf = 
      dynamic_pointer_cast<BoundedSurface, ParamSurface>(surface);
    if (bsurf.get())
      bddomain = &(bsurf->parameterDomain());

  ASSERT(surf != 0);

  SISLSurf* sislsf = GoSurf2SISL(*surf, false);
  double epsge = 1.0e-6;
  int numintpt;  // number of single extremal points
  double* pointpar = 0; // array containing the parameter values of single extremal. pt.
  int numintcr; // number of extremal curves
  SISLIntcurve** intcurves = 0;
  int stat = 0;

  s1921(sislsf, dir.begin(), dir.dimension(), 0.0, epsge, 
	&numintpt, &pointpar, &numintcr, &intcurves, &stat);
  MESSAGE_IF(stat!=0, "s1921 returned code: " << stat); 


  // Check if any of the found extremal points are better than the
  // current most extreme point
  vector<Point> curr_pnt;
  vector<double> curr_par;
  int ki;
  for (ki=0; ki<numintpt; ++ki)
    {
      // Evaluate surface
      Point pos = surf->ParamSurface::point(pointpar[2*ki],pointpar[2*ki+1]);
      if (ext_pnt.dimension() == 0 || pos*dir > ext_pnt*dir)
	{
	  curr_pnt.push_back(pos);
	  curr_par.insert(curr_par.end(), pointpar+2*ki, pointpar+2*(ki+1));
	}
    }

  for (ki=0; ki<numintcr; ++ki)
    {
      Point pos = surf->ParamSurface::point(intcurves[ki]->epar1[0],
					    intcurves[ki]->epar1[1]);
      if (ext_pnt.dimension() == 0 || pos*dir > ext_pnt*dir)
	{
	  curr_pnt.push_back(pos);
	  curr_par.insert(curr_par.end(), intcurves[ki]->epar1, 
			  intcurves[ki]->epar1+2);
	}

      double *pp = intcurves[ki]->epar1 + 2*(intcurves[ki]->ipoint-1);
      pos = surf->ParamSurface::point(pp[0], pp[1]);
      if (ext_pnt.dimension() == 0 || pos*dir > ext_pnt*dir)
	{
	  curr_pnt.push_back(pos);
	  curr_par.insert(curr_par.end(), pp, pp+2);
	}
    }

  if (sislsf)
    freeSurf(sislsf);
  if (pointpar)
    free(pointpar);
  if (intcurves)
    freeIntcrvlist(intcurves, numintcr);

  if (curr_pnt.size() == 0)
    {
      // No better extremal point is found. 
      return false;
    }

  if (bddomain != 0)
    {
      // The surface was trimmed
      // Remove extremal points which are outside the domain
      for (ki=0; ki<(int)curr_pnt.size();)
	{
	  Vector2D param(curr_par[2*ki], curr_par[2*ki+1]);
	  if (!bddomain->isInDomain(param, epsge))
	    {
	      curr_pnt.erase(curr_pnt.begin()+ki);
	      curr_par.erase(curr_par.begin()+2*ki, curr_par.begin()+2*(ki+1));
	    }
	  else
	    ++ki;
	}
    }
  if (curr_pnt.size() > 0)
    {
      //  Fetch the best extremal point
      modified = true;
      ext_pnt = curr_pnt[0];
      ext_par[0] = curr_par[0];
      ext_par[1] = curr_par[1];
      for (ki=1; ki<(int)curr_pnt.size(); ++ki)
	{
	  if (curr_pnt[ki]*dir > ext_pnt*dir)
	    {
	      ext_pnt = curr_pnt[ki];
	      ext_par[0] = curr_par[2*ki];
	      ext_par[1] = curr_par[2*ki+1];
 	    }
	}
    }  

  else if (bddomain != 0)
    {
      // It can be an extremal point inside the face that is better than
      // the previous one.
      // First get the extreme points on the boundary
      vector<CurveLoop> bd_loops = bsurf->allBoundaryLoops();
      for (ki=0; ki<(int)bd_loops.size(); ++ki)
	{
	  vector<shared_ptr<ParamCurve> > curr_loop(bd_loops[ki].begin(),
						    bd_loops[ki].end());
	  shared_ptr<CompositeCurve> comp_cv = 
	    shared_ptr<CompositeCurve>(new CompositeCurve(toptol.gap,
							  toptol.neighbour,
							  toptol.kink,
							  toptol.bend,
							  curr_loop));
				       
	  int idx;
	  Point bd_ext;
	  double bd_par;
	  comp_cv->extremalPoint(dir, bd_ext, idx, &bd_par);
	  if (ext_pnt.dimension() == 0 || bd_ext*dir > ext_pnt*dir)
	    {
	      modified = true;
	      ext_pnt = bd_ext;
	      Point param = bsurf->getSurfaceParameter(ki, idx, bd_par);
	      ext_par[0] = param[0];
	      ext_par[1] = param[1];
	    }
	}
      
      // Triangulate trimmed surface
      int n, m;
      double density = 1.0;
      int min_nmb = 4, max_nmb = 50;
      setResolutionFromDensity(surface, density, min_nmb, max_nmb, tol2d, 
      			       n, m);

      RectDomain dom = surface->containingDomain();
      
      // Fetch constant parameter curves in the 1. parameter direction
      int min_samples = 1;
      double u1 = dom.umin();
      double u2 = dom.umax();
      double udel = (u2 - u1)/(n+1);
      double par[2];
      par[0] = u1+udel;
      while (par[0] < u2)
	{
	  vector<shared_ptr<ParamCurve> > crvs = 
	    surface->constParamCurves(par[0], false);
	  if (crvs.size() == 0)
	    {
	      par[0] += udel;
	      continue;  // Outside trimmed surface
	    }
#ifdef DEBUG
	  for (size_t kr=0; kr<crvs.size(); ++kr)
	    {
	      shared_ptr<SplineCurve> tmpspl(crvs[kr]->geometryCurve());
	      tmpspl->writeStandardHeader(of1);
	      tmpspl->geometryCurve()->write(of1);
	    }
#endif
	  // Distribute sampling points
	  double av_len = 0.0;
	  vector<double> cv_len(crvs.size());
	  double curr_len = 0.0;
	  for (size_t kr=0; kr<crvs.size(); ++kr)
	    {
	      double len = crvs[kr]->estimatedCurveLength();
	      av_len += len;
	      cv_len[kr] = len;
	      curr_len += (crvs[kr]->endparam()-crvs[kr]->startparam());
	    }
	  av_len /= (double)crvs.size();
	  
	  // Evaluate sampling points
	  int curr_nmb = (int)(m*(curr_len/(dom.vmax()-dom.vmin()))) + 1;
	  for (size_t kr=0; kr<crvs.size(); ++kr)
	    {
	      int nmb = (int)(curr_nmb*cv_len[kr]/av_len);
	      nmb = std::max(nmb, min_samples);
	      double v1 = crvs[kr]->startparam();
	      double v2 = crvs[kr]->endparam();
	      double vdel = (v2 - v1)/(double)(nmb+1);
	      par[1] = v1 + vdel;
	    //for (kj=0, par[pt_dir]=v1+vdel; kj<nmb; ++kj, par[pt_dir]+=vdel)
	      while (par[1] < v2)
		{
		  Point pos = crvs[kr]->point(par[1]);
		  if (pos*dir > ext_pnt*dir)
		    {
#ifdef DEBUG
		      of1 << "400 1 0 4 255 0 0 255" << std::endl;
		      of1 << "1" << std::endl;
		      of1 << pos << std::endl;
#endif
		      modified = true;
		      ext_pnt = pos;
		      ext_par[0] = par[0];
		      ext_par[1] = par[1];
		    }
		  par[1] += vdel;
		}
	    }
	  par[0] += udel;
	}


      // Too time consuming
      // shared_ptr<GeneralMesh> mesh;
      // tesselateOneSrf(surface, mesh, tol2d, n, m);

      // // Get the most extreme triangulation nodes
      // double *nodes = mesh->vertexArray();
      // // int nmb_nodes = mesh->numVertices();
      // int num_triang = mesh->numTriangles();
      // double *par_nodes = mesh->paramArray();
      // unsigned int *triang_idx = mesh->triangleIndexArray();
      // for (ki=0; ki<num_triang; ++ki)
      // 	{
      // 	  // Due to the structure of the tesselation, the points must
      // 	  // be handled more than once
      // 	  for (int kj=0; kj<3; ++kj)
      // 	    {
      // 	      Point node_ext(nodes+3*triang_idx[ki+kj], 
      // 			     nodes+3*triang_idx[ki+kj]+3, false);
      // 	      if (ext_pnt.dimension() == 0 || node_ext*dir > ext_pnt*dir)
      // 		{
      // 		  modified = true;
      // 		  ext_pnt = node_ext;
      // 		  ext_par[0] = par_nodes[2*triang_idx[ki+kj]];
      // 		  ext_par[1] = par_nodes[2*triang_idx[ki+kj]+1];
      // 		}
      // 	    }
      // 	}

      // Use this value as a start point for an extreme point iteration

    }
  return modified;
}

  //===========================================================================
  void 
  SurfaceModelUtils::setResolutionFromDensity(shared_ptr<ParamSurface> surf,
					      double density, double tol2d,
					      int min_nmb, int max_nmb,
					      int& u_res, int& v_res)
  //===========================================================================
  {
	// Estimate size of surface/underlying surface
    shared_ptr<ParamSurface> sf;
    shared_ptr<BoundedSurface> bd_surf = 
      dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);
    if (bd_surf.get())
      {
	// A trimmed surface is found
	// Get underlying surface 
	sf = bd_surf->underlyingSurface();
	if (bd_surf->isIsoTrimmed(tol2d))
	  {
	    RectDomain domain = bd_surf->containingDomain();
	    RectDomain dom2 = sf->containingDomain();
	    double umin = std::max(domain.umin(), dom2.umin());
	    double umax = std::min(domain.umax(), dom2.umax());
	    double vmin = std::max(domain.vmin(), dom2.vmin());
	    double vmax = std::min(domain.vmax(), dom2.vmax());
    
	    vector<shared_ptr<ParamSurface> > sfs = sf->subSurfaces(umin, vmin, umax, vmax);
	    sf = sfs[0];
	  }
      }
    else 
      sf = surf;
	
    double len_u, len_v;
    GeometryTools::estimateSurfaceSize(*sf, len_u, len_v);

    u_res = (int)(len_u/density);
    v_res = (int)(len_v/density);
    double fac = len_u/len_v;
    u_res = std::max(min_nmb, std::min(u_res, (int)(fac*max_nmb)));
    v_res = std::max(min_nmb, std::min(v_res, (int)(max_nmb/fac)));

  }

  //===========================================================================
void SurfaceModelUtils::tesselateOneSrf(shared_ptr<ParamSurface> surf,
					shared_ptr<GeneralMesh>& mesh,
					double tol2d, int n, int m)
//===========================================================================
  {
      ClassType type = surf->instanceType();
      if (type == Class_SplineSurface)
      {
	  RectangularSurfaceTesselator tesselator(*surf.get(), n, m);
	  tesselator.tesselate();
	  mesh = tesselator.getMesh();
      }
      else if (type == Class_BoundedSurface)
      {
	  shared_ptr<BoundedSurface> bd_surf = 
	      dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);
	  if (bd_surf->isIsoTrimmed(tol2d))
	  {
	      // Get surrounding domain
	      RectDomain domain = bd_surf->containingDomain();
    
	      // Get smallest surrounding surface
	      shared_ptr<ParamSurface> base_sf = bd_surf->underlyingSurface();
	      while (base_sf->instanceType() == Class_BoundedSurface)
		  base_sf = dynamic_pointer_cast<BoundedSurface, ParamSurface>(base_sf)->underlyingSurface();
	      RectDomain dom2 = base_sf->containingDomain();  // To avoid problems due to numerics
	      double umin = std::max(domain.umin(), dom2.umin());
	      double umax = std::min(domain.umax(), dom2.umax());
	      double vmin = std::max(domain.vmin(), dom2.vmin());
	      double vmax = std::min(domain.vmax(), dom2.vmax());
    
	      vector<shared_ptr<ParamSurface> > sfs = base_sf->subSurfaces(umin, vmin, umax, vmax);
	      RectangularSurfaceTesselator tesselator(*(sfs[0].get()), n, m, false);
	      tesselator.tesselate();
	      mesh = tesselator.getMesh();
	  }
	  else
	  {
	      ParametricSurfaceTesselator tesselator(*surf.get());
	      int n2, m2;
	      tesselator.getRes(n2, m2);
	      if (n == n2 && m == m2)
		tesselator.tesselate();
	      else
		tesselator.changeRes(n, m);
	      mesh = tesselator.getMesh();
	  }
      }
  }

//===========================================================================
void SurfaceModelUtils::triangulateFaces(vector<shared_ptr<ftSurface> >& faces,
					 shared_ptr<ftPointSet>& triang,
					 double tol)
//===========================================================================
{
  vector<pair<int, int> > pnt_range;
  int nmb_pnt = 0;
  int nmb_sample = 35;
  vector<shared_ptr<ftSurface> > faces2;
  for (size_t ki=0; ki<faces.size(); ++ki)
    {
      shared_ptr<ftSurface> curr_face = faces[ki];
      shared_ptr<ParamSurface> surf = curr_face->surface();
      shared_ptr<ftPointSet> local_triang = shared_ptr<ftPointSet>(new ftPointSet());
      vector<int> local_corner;
      RectDomain dom = surf->containingDomain();
      AdaptSurface::createTriangulation(surf, dom, local_triang, local_corner,
					false, nmb_sample);
      triang->append(local_triang);

      // Handle common boundaries
      // First find pairs of faces meeting at a common boundary
      vector<ftSurface*> neighbours;
      curr_face->getAdjacentFaces(neighbours);
      for (size_t kj=0; kj<neighbours.size(); ++kj)
	{
	  // Check if this face is meshed already
	  size_t kr;
	  for (kr=0; kr<faces2.size(); ++kr)
	    if (faces2[kr].get() == neighbours[kj])
	      {
		// A common boundary is found
		triang->markLocalBoundary(faces2[kr], pnt_range[kr].first, 
					  pnt_range[kr].second, curr_face,
					  nmb_pnt, triang->size(), tol);
	      }
	}
      // Set range information
      faces2.push_back(curr_face);
      pnt_range.push_back(make_pair(nmb_pnt, triang->size()));
      nmb_pnt = triang->size();
    }

#ifdef DEBUG
  std::ofstream of0("triang.g2");
  triang->write(of0);
  std::ofstream pointsout("pointsdump.g2");
  vector<Vector3D> bd_nodes;
  vector<Vector3D> inner_nodes;
  int k2;
  for (k2=0; k2<(int)triang->size(); ++k2)
    {
      if ((*triang)[k2]->isOnBoundary())
	bd_nodes.push_back((*triang)[k2]->getPoint());
      else
	inner_nodes.push_back((*triang)[k2]->getPoint());
    }
		
  pointsout << "400 1 0 4 255 0 0 255" << std::endl;
  pointsout << bd_nodes.size() << std::endl;
  for (k2=0; k2<(int)bd_nodes.size(); ++k2)
    pointsout << bd_nodes[k2][0] << " " << bd_nodes[k2][1] << " " << bd_nodes[k2][2] << std::endl;
  pointsout << "400 1 0 4 0 255 0 255" << std::endl;
  pointsout << inner_nodes.size() << std::endl;
  for (k2=0; k2<(int)inner_nodes.size(); ++k2)
    pointsout << inner_nodes[k2][0] << " " << inner_nodes[k2][1] << " " << inner_nodes[k2][2] << std::endl;
#endif

}

//===========================================================================
void 
SurfaceModelUtils::reduceUnderlyingSurface(shared_ptr<BoundedSurface>& bd_sf,
					   vector<shared_ptr<CurveOnSurface> >& cvs)
//===========================================================================
{
  double domainfac = 10.0;
  double redfac = 0.1;

  // Check if the domain size should be reduced
  RectDomain dom1 = bd_sf->containingDomain();
  RectDomain dom2 = bd_sf->underlyingSurface()->containingDomain();
  double ll1 = dom1.diagLength();
  double ll2 = dom2.diagLength();
  if (ll2 > domainfac*ll1)
    {
      double u1 = std::max(dom1.umin()-redfac*ll1, dom2.umin());
      double u2 = std::min(dom1.umax()+redfac*ll1, dom2.umax());
      double v1 = std::max(dom1.vmin()-redfac*ll1, dom2.vmin());
      double v2 = std::min(dom1.vmax()+redfac*ll1, dom2.vmax());
      vector<shared_ptr<ParamSurface> > sub =
	bd_sf->underlyingSurface()->subSurfaces(u1, v1, u2, v2);
      if (sub.size() == 1)
	{
	  bd_sf->replaceSurf(sub[0]);
	  for (size_t kr=0; kr<cvs.size(); ++kr)
	    cvs[kr]->setUnderlyingSurface(sub[0]);
	}
    }
}	  
