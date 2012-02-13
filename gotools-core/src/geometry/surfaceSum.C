//===========================================================================
//                                                                           
// File: surfaceSum.C                                                        
//                                                                           
// Created: Thu Jul 14 11:33:38 2005                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: surfaceSum.C,v 1.2 2005-09-12 12:40:41 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/geometry/GeometryTools.h"
#include <memory>

using std::vector;
using std::setprecision;
using std::endl;
using std::pair;
using std::make_pair;

namespace Go
{

  shared_ptr<SplineSurface>
  GeometryTools::surfaceSum(const SplineSurface& sf1, double fac1,
	     const SplineSurface& sf2, double fac2, double num_tol)

    //********************************************************************
    // Addition of two signed SplineSurfaces, i.e. this function can
    // also be used for subtraction. The surfaces is assumed to live on
    // the same parameter domain, but may have different knot vectors.
    //********************************************************************
  {
    // Check input
    ALWAYS_ERROR_IF(fabs(sf1.startparam_u() - sf2.startparam_u()) > num_tol ||
		    fabs(sf1.endparam_u() - sf2.endparam_u()) > num_tol ||
		    fabs(sf1.startparam_v() - sf2.startparam_v()) > num_tol ||
		    fabs(sf1.endparam_v() - sf2.endparam_v()) > num_tol,
		    "Inconsistent parameter domain.");

    // For the time being
    if (sf1.rational() || sf2.rational()) {
	THROW("Sum of rational surfaces is not implemented");
    }

    // Make copy of surfaces
    vector<shared_ptr<SplineSurface> > surfaces;
    surfaces.reserve(2);
    shared_ptr<SplineSurface> sf;
// #ifdef _MSC_VER
//     sf = shared_ptr<SplineSurface>(dynamic_cast<SplineSurface*>(sf1.clone()));
// #else
    sf = shared_ptr<SplineSurface>(sf1.clone());
// #endif
    surfaces.push_back(sf);
// #ifdef _MSC_VER
//     sf = shared_ptr<SplineSurface>(dynamic_cast<SplineSurface*>(sf2.clone()));
// #else
    sf = shared_ptr<SplineSurface>(sf2.clone());
// #endif
    surfaces.push_back(sf);

    // Make sure that the surfaces live on the same knot vector
    GeometryTools::unifySurfaceSplineSpace(surfaces, num_tol);

    // Add signed coefficients
    vector<double> coefs;
    int nmb_coefs_u = surfaces[0]->numCoefs_u();
    int nmb_coefs_v = surfaces[0]->numCoefs_v();
    int dim = surfaces[0]->dimension();
    coefs.resize(dim*nmb_coefs_u*nmb_coefs_v);
    int ki;
    std::vector<double>::iterator s1 = surfaces[0]->coefs_begin();
    std::vector<double>::iterator s2 = surfaces[1]->coefs_begin();
    for (ki=0; ki<dim*nmb_coefs_u*nmb_coefs_v; ki++)
      coefs[ki] = fac1*s1[ki] + fac2*s2[ki];

    // Create output curve
    shared_ptr<SplineSurface>
	surfacesum(new SplineSurface(nmb_coefs_u, nmb_coefs_v,
				     surfaces[0]->order_u(),
				     surfaces[0]->order_v(),
				     surfaces[0]->basis_u().begin(),
				     surfaces[0]->basis_v().begin(),
				     &coefs[0], dim, false));

    return surfacesum;
  }
}
