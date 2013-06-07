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

#include "GoTools/igeslib/IGESconverter.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/Utils.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>


using std::ostream;
using std::istream;
using std::vector;
using std::string;
using namespace Go;

//*****************************************************************************
//
// NOTE: 
//
//*****************************************************************************

//-----------------------------------------------------------------------------
void IGESconverter::writedisp(ostream& os)
//-----------------------------------------------------------------------------
{
    MESSAGE_IF(!filled_with_data_, "Cannot write. I contain no data.");
    if (!filled_with_data_) return;

    os << "list dispgeom" << std::endl;
    int i;
    for (i=0; i< int(geom_.size()); i++)
      {
	
	if (geom_[i]->instanceType() == Class_SplineSurface)
	  writedispsurface(os, dynamic_cast<SplineSurface*>(geom_[i].get()));
	else if (geom_[i]->instanceType() == Class_SplineCurve)
	  writedispcurve(os, dynamic_cast<SplineCurve*>(geom_[i].get()));
	else if (geom_[i]->instanceType() == Class_BoundedSurface)
	  writedispboundedSurf(os,
			      dynamic_cast<BoundedSurface*>(geom_[i].get()));
    }
    os << "end" << std::endl;
}

//-----------------------------------------------------------------------------
void IGESconverter::readdisp(istream& is)
//-----------------------------------------------------------------------------
{
//     GO_ERROR("Not implemented yet.", InputError());

    MESSAGE("Implemented for a spline surface, without reference to format definition!");

    while (!is.eof()) {
	shared_ptr<GeomObject> geom_obj;
	string indicator, list_name, class_name;
	is >> indicator >> list_name >> class_name;

	if (class_name == "surf")
	    geom_obj = readdispsurface(is);
	else
	    THROW("Not implemented yet!");

	geom_.push_back(geom_obj);

	string surf_name;
	is >> surf_name;
	indicator.resize(0);
	is >> indicator;
	double x, y, z;
	while (indicator != "end") {
	    indicator.resize(0);
	    is >> x >> y >> z >> indicator;
	}

	Utils::eatwhite(is);
    }
    filled_with_data_ = true;
}

//-----------------------------------------------------------------------------
shared_ptr<SplineSurface> IGESconverter::readdispsurface(istream& is)
//-----------------------------------------------------------------------------
{
//     string text1;
//     string text2;
    int dim;
    int n1;
    int n2;
    int k1;
    int k2;
    int rat;

    //    is >> text1 >> text2
    is >> n1 >> n2 >> k1 >> k2 >> dim >> rat;

    //    text2.strip();
    //  ALWAYS_ERROR_IF(text2 != "surf",
    //	"Expecting input of surface, given " << text2 ".", InputError());

    bool rational = (rat == 2) ? true : false;

    int dim1 = dim + (rational==true);

    vector<double> knots1(k1+n1);
    vector<double> knots2(k2+n2);
    vector<double> coefs(n1*n2*dim1);

    int i, j;
    for (i = 0; i < k1+n1; ++i) // We read first knot vector.
	is >> knots1[i];

    for (i = 0; i < k2+n2; ++i) // We read second knot vector.
	is >> knots2[i];

    for (i = 0; i < n1*n2; ++i) // We read the coefs (and weights if rational).
	for (j = 0; j < dim1; ++j)
	    is >> coefs[i*dim+j];

    // We set values in our return surface.
    shared_ptr<SplineSurface> surf = shared_ptr<SplineSurface>
	(new SplineSurface(n1, n2, k1, k2, knots1.begin(), knots2.begin(),
			     coefs.begin(), dim, rational));

    return surf;
}

//-----------------------------------------------------------------------------
void IGESconverter::writedispsurface(ostream& os, SplineSurface* surf)
//-----------------------------------------------------------------------------
{
  int dim = surf->dimension();
  int in1 = surf->numCoefs_u();
  int in2 = surf->numCoefs_v();
  int ik1 = surf->order_u();
  int ik2 = surf->order_v();
  bool rational = surf->rational();
  int dim1 = dim + (rational==true);
  os << "surf " << in1 << " " << in2 << " " << ik1 << " ";
  os << ik2  << " " << dim << " ";
  os << (rational ? 2 : 1) << std::endl;

  int ki, kj;
  std::vector<double>::const_iterator it;
  for (it=surf->basis_u().begin(); it!=surf->basis_u().end(); it++)
    os << *it << std::endl;

  for (it=surf->basis_v().begin(); it!=surf->basis_v().end(); it++)
    os << *it << std::endl;

  std::vector<double>::const_iterator co = (rational) ? surf->rcoefs_begin()
    : surf->coefs_begin();
  for (ki=0; ki<in1*in2; ki++)
    {
      for (kj=0; kj<dim1; kj++)
	os << co[ki*dim1+kj] << " ";
      os << std::endl;
    }

  os << "ambient 0.1 0.1 0.1 " << std::endl;
  os << "diffuse 0.0 0.3 0.3 " << std::endl;
  os <<  std::endl;
}

//-----------------------------------------------------------------------------
void IGESconverter::writedispcurve(ostream& os, SplineCurve* crv)
//-----------------------------------------------------------------------------
{
  int dim = crv->dimension();
  int in = crv->numCoefs();
  int ik = crv->order();
  bool rational = crv->rational();
  int dim1 = dim + (rational==true);
  os << "curve " << in << " " << ik  << " " << dim << " ";
  os << (rational ? 2 : 1) << std::endl;

  int ki, kj;
  std::vector<double>::const_iterator it;
  for (it=crv->basis().begin(); it!=crv->basis().end(); it++)
    os << *it << std::endl;

  std::vector<double>::const_iterator co = (rational) ? crv->rcoefs_begin()
    : crv->coefs_begin();
  for (ki=0; ki<in; ki++)
    {
      for (kj=0; kj<dim1; kj++)
	os << co[ki*dim1+kj] << "  ";
      os << std::endl;
    }

  os << "diffuse 1.0 0.0 0.0 " << std::endl;
  os <<  std::endl;
}

//-----------------------------------------------------------------------------
void IGESconverter::writedispboundedSurf(ostream& os, BoundedSurface *bdsf)
//-----------------------------------------------------------------------------
{
  shared_ptr<ParamSurface> surf = bdsf->underlyingSurface();
  writedispsurface(os, dynamic_cast<SplineSurface*>(surf.get()));

  int ki;
  for (ki=0; ki<bdsf->numberOfLoops(); ki++)
    {
      shared_ptr<CurveLoop> loop = bdsf->loop(ki);
      std::vector< shared_ptr<ParamCurve> >::const_iterator it;
      for (it=loop->begin(); it!=loop->end(); it++)
	{
	  CurveOnSurface* surfcrv 
	    = dynamic_cast<CurveOnSurface*>((*it).get());
	  SplineCurve* crv 
	    = dynamic_cast<SplineCurve*>(surfcrv->spaceCurve().get());
	  writedispcurve(os, crv);
	}
    }
  os << std::endl << std::endl;
}
