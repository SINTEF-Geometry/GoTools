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


} // namespace Go
