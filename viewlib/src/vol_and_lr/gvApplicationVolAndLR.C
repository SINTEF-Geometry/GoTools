//===========================================================================
//                                                                           
// File: gvApplicationVolAndLR.C                                             
//                                                                           
// Created: Thu Jun 20 08:27:58 2013                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================


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
