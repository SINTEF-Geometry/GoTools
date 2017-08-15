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


#include "GoTools/viewlib/vol_and_lr/gvApplicationVolAndLR.h"
#include "GoTools/viewlib/gvView.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/tesselator/GenericTriMesh.h"
#include "GoTools/tesselator/CurveTesselator.h"
#include "GoTools/tesselator/RectangularSurfaceTesselator.h"
#include "GoTools/tesselator/ParametricSurfaceTesselator.h"
#include "GoTools/trivariate/RectangularVolumeTesselator.h"
#include <assert.h>

#include <QMenuBar>

using namespace Go;
using std::vector;


//===========================================================================
gvApplicationVolAndLR::gvApplicationVolAndLR(std::auto_ptr<DataHandler> dh,
					     QWidget * parent,
					     const char * name,
					     Qt::WindowFlags f)
    : gvApplication(dh, parent, name, f)
//===========================================================================
{
    buildExtraGUI();
}


//===========================================================================
gvApplicationVolAndLR::~gvApplicationVolAndLR()
//===========================================================================
{
    std::cout << "Hmmmm" << std::endl;
}

//===========================================================================
void gvApplicationVolAndLR::view_reset()
//===========================================================================
{
    vector<double> translate_vec(3, 0.0);
    for (int i = 0; i < data_.numObjects(); ++i)
    {
	GeneralMesh* gen_mesh = getMesh(i);
	if (gen_mesh)
	{
	    gen_mesh->translate(translate_vec);
	}
	else
	{
	    MESSAGE("Object does not support translation!");
	}
    }

    view_->setOriginFocalPoint(false);

    view_->focusOnBox();
    view_->updateGL();
}

//===========================================================================
void gvApplicationVolAndLR::translate_to_origin()
//===========================================================================
{
  data_.disableUpdates();
  bool was_translated = false;
  // We first locate the boundingbox for the whole model.
  const BoundingBox& mod_bd_box = data_.boundingBox();
  Point low = mod_bd_box.low();
  Point high = mod_bd_box.high();
  Point transl_pt = -0.5*(low + high);
  vector<double> translate_vec(transl_pt.begin(), transl_pt.end());
  for (int i = 0; i < data_.numObjects(); ++i)
    {
      if (data_.getSelectedStateObject(i))
	{
	  if (was_translated == true)
	    {
	      std::cout << "Currently only supporting one selected object." << std::endl;
	    }
	  else
	    {
#if 1 // We now translate the mesh, not the geometry.
		  Tesselator* tess = data_.tesselator(i).get();
		  ParametricSurfaceTesselator* param_surf_tess = dynamic_cast<ParametricSurfaceTesselator*>(tess);
		  if (param_surf_tess)
		  {
		      shared_ptr<GenericTriMesh> gen_mesh = param_surf_tess->getMesh();
		      gen_mesh->translate(translate_vec);
#ifndef NDEBUG
		      {
			  double debug_val = 1.0;
		      }
#endif
		  }
//		  tess->tesselate();
#else
	      was_translated = true;
	      shared_ptr<GeomObject> obj(data_.object(i));
	      // BoundingBox bd_box = obj->boundingBox();
	      if (obj->instanceType() == Class_SplineSurface)
		{
		  std::cout << "Translating SplineSurface!" << std::endl;
		  SplineSurface* spline_sf =
		    dynamic_cast<SplineSurface*>(obj.get());
		  int dim = spline_sf->dimension();
		  assert(!spline_sf->rational());
		  vector<double>::iterator iter = spline_sf->coefs_begin();
		  while (iter != spline_sf->coefs_end())
		    {
		      for (int kj = 0; kj < dim; ++kj)
			iter[kj] -= center_pt[kj];
		      iter += dim;
		    }
		  // We must retesselate.
		  Tesselator* tess = data_.tesselator(i).get();
		  tess->tesselate();
		}
	      else if (obj->instanceType() == Class_LRSplineSurface)
		{
		  std::cout << "Translating LRSplineSurface!" << std::endl;
		  LRSplineSurface* lrspline_sf =
		    dynamic_cast<LRSplineSurface*>(obj.get());
		  int dim = lrspline_sf->dimension();
		  assert(!lrspline_sf->rational());

		  auto iter = lrspline_sf->basisFunctionsBegin();
		  while (iter != lrspline_sf->basisFunctionsEnd())
		    {
		      LRBSpline2D* bsb = iter->second.get();
		      Point coef = bsb->Coef();
		      double gamma = bsb->gamma();
		      for (int kj = 0; kj < dim; ++kj)
			coef[kj] -= center_pt[kj];
		      bsb->setCoefAndGamma(coef, gamma);
		      ++iter;
		    }
		  // We must retesselate.
		  Tesselator* tess = data_.tesselator(i).get();
		  tess->tesselate();
		  
		}
	      else if (obj->instanceType() == Class_BoundedSurface)
		{
		  std::cout << "Translating BoundedSurface!" << std::endl;
		  BoundedSurface* bd_sf =
		    dynamic_cast<BoundedSurface*>(obj.get());
		  int dim = bd_sf->dimension();
		  Point trans_vec = -center_pt;
		  const double deg_eps = 1e-03;
		  BoundedUtils::translateBoundedSurf(trans_vec, *bd_sf, deg_eps);

		  // We must retesselate.
		  Tesselator* tess = data_.tesselator(i).get();
		  tess->tesselate();		  
		}
	      else
		{
		  std::cout << "Supporting translation of SplineSurface, LRSplineSurface & BoundedSurface only!"
			    << std::endl;
		}
#endif
	    }
	}
    }

    data_.enableUpdates();
    data_.updateObservers();
}


//===========================================================================
void gvApplicationVolAndLR::move_vertices_to_origin()
//===========================================================================
{
  data_.disableUpdates();
  // We first locate the boundingbox for the whole model.
  const BoundingBox& mod_bd_box = data_.boundingBox();
  Point low = mod_bd_box.low();
  Point high = mod_bd_box.high();
  Point transl_pt = -0.5*(low + high);
  vector<double> translate_vec(transl_pt.begin(), transl_pt.end());
  for (int i = 0; i < data_.numObjects(); ++i)
  {
      GeneralMesh* gen_mesh = getMesh(i);
      if (gen_mesh)
      {
	  std::cout << "Translating: " << translate_vec[0] << " " << translate_vec[1] << " " << translate_vec[2] << std::endl;
	  gen_mesh->translate(translate_vec);
      }
      else
      {
	  MESSAGE("Object does not support translation!");
      }
  }

  view_->setOriginFocalPoint(true);

  data_.enableUpdates();
  data_.updateObservers();

}


//===========================================================================
void gvApplicationVolAndLR::buildExtraGUI()
//===========================================================================
{
    //---------------------------------------------------------------------
    //------------ second menu item: View ------------------------------
    //---------------------------------------------------------------------

    view_menu_->addAction("Move vertices to origin", this, 
			    SLOT(move_vertices_to_origin()));

    //---------------------------------------------------------------------
    //------------ fifth menu item: Object --------------------------------
    //---------------------------------------------------------------------

    object_menu_->addAction("Translate to origin", this, 
			    SLOT(translate_to_origin()));

}

//===========================================================================
GeneralMesh* gvApplicationVolAndLR::getMesh(int object_id)
//===========================================================================
{
    Tesselator* tess = data_.tesselator(object_id).get();

    // All the meshes are derived from the class GeneralMesh.
    GeneralMesh* gen_mesh = NULL;

    // 1D: CurveTesselator: shared_ptr<LineStrip> mesh_;
    CurveTesselator* curve_tess = dynamic_cast<CurveTesselator*>(tess);
    if (curve_tess)
    {
	gen_mesh = curve_tess->getMesh().get();
    }

    // 2D: ParametricSurfaceTesselator: shared_ptr<GenericTriMesh> mesh_;
    ParametricSurfaceTesselator* param_surf_tess = dynamic_cast<ParametricSurfaceTesselator*>(tess);
    if (param_surf_tess)
    {
	gen_mesh = param_surf_tess->getMesh().get();
    }

    // 2D: RectangularSurfaceTesselator: shared_ptr<RegularMesh> mesh_;
    RectangularSurfaceTesselator* rect_surf_tess = dynamic_cast<RectangularSurfaceTesselator*>(tess);
    if (rect_surf_tess)
    {
	gen_mesh = rect_surf_tess->getMesh().get();
    }

    // 3D: RectangularVolumeTesselator (tesselating the hull => 2.5D): shared_ptr<RegularVolMesh> mesh_;
    RectangularVolumeTesselator* reg_vol_tess = dynamic_cast<RectangularVolumeTesselator*>(tess);
    if (reg_vol_tess)
    {
	gen_mesh = reg_vol_tess->getMesh().get();
    }

    return gen_mesh;
}
