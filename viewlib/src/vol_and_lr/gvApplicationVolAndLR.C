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
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include <assert.h>

#include <Qt/qmenubar.h>

using namespace Go;
using std::vector;


//===========================================================================
gvApplicationVolAndLR::gvApplicationVolAndLR(std::auto_ptr<DataHandler> dh,
					     QWidget * parent,
					     const char * name,
					     Qt::WFlags f)
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
void gvApplicationVolAndLR::translate_to_origin()
//===========================================================================
{
  data_.disableUpdates();
  bool was_translated = false;
  // We first locate the boundingbox for the whole model.
  const BoundingBox& mod_bd_box = data_.boundingBox();
  Point low = mod_bd_box.low();
  Point high = mod_bd_box.high();
  Point center_pt = 0.5*(low + high);
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
	      else
		{
		  std::cout << "Supporting translation of SplineSurface only!" << std::endl;
		}
	    }
	}
    }

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

    object_menu_->addAction("Translate to origin", this, 
			    SLOT(translate_to_origin()));

}
