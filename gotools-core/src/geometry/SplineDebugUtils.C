//===========================================================================
//                                                                           
// File: SplineDebugUtils.C                                                  
//                                                                           
// Created: Tue Sep 27 11:04:28 2005                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: SplineDebugUtils.C,v 1.1 2005-10-03 08:02:49 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/geometry/SplineDebugUtils.h"

#include <fstream>


using std::vector;
using std::endl;
using std::setprecision;

namespace Go
{


//===========================================================================
void writeSpaceParamCurve(const SplineCurve& pcurve, std::ostream& os, double z)
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
void writeTrimmedInfo(BoundedSurface& bd_sf,
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
void objToFile(GeomObject* geom_obj, char *to_file)
//===========================================================================
{
    if (geom_obj) {
	std::ofstream debug(to_file);
	geom_obj->writeStandardHeader(debug);
	geom_obj->write(debug);
    }
}


//===========================================================================
void objsToFile(vector<shared_ptr<GeomObject> >& geom_objs, char *to_file)
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
void writeSISLFormat(const SplineCurve& spline_cv, std::ostream& os)
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
