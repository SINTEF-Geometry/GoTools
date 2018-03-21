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

#include "GoTools/geometry/SISLconversion.h"
#include "sisl.h"

using namespace std;

namespace Go
{

SISLCurve* Curve2SISL( const SplineCurve& cv, bool copy)
{
    std::vector<double>::const_iterator coef;
    int kind;
    if (cv.rational()) {
	coef = cv.rcoefs_begin();
	kind = 2;
    } else {
	coef = cv.coefs_begin();
	kind = 1;
    }
    return newCurve(cv.numCoefs(), cv.order(),
		    const_cast<double*>(&(*(cv.basis().begin()))),
		    const_cast<double*>(&(*coef)),
		    kind, cv.dimension(), copy);
}



SISLCurve* Curve2SISL_rat( const SplineCurve& cv)
{
    std::vector<double>::const_iterator coef;
    int kind = 2;
    if (cv.rational()) 
      {
	coef = cv.rcoefs_begin();
	return newCurve(cv.numCoefs(), cv.order(),
			const_cast<double*>(&(*(cv.basis().begin()))),
			const_cast<double*>(&(*coef)),
			kind, cv.dimension(), 1);
      }
    else 
      {
	coef = cv.coefs_begin();
	int in = cv.numCoefs();
	int dim = cv.dimension();
	vector<double> sc(in*(dim+1));
	int ki;
	for (ki=0; ki<in; ++ki, coef+=dim)
	  {
	    std::copy(coef, coef+dim, sc.begin()+ki*(dim+1));
	    sc[ki*(dim+1)+dim] = 1.0;
	  }
	return newCurve(cv.numCoefs(), cv.order(),
			const_cast<double*>(&(*(cv.basis().begin()))),
			const_cast<double*>(&(sc[0])),
			kind, cv.dimension(), 1);
      }
}


SplineCurve* SISLCurve2Go( const SISLCurve* const cv)
{
    double* coefsstart;
    bool rational;
    if (cv->ikind == 2) {
	coefsstart = cv->rcoef;
	rational = true;
    } else if (cv->ikind == 1) {
	coefsstart = cv->ecoef;
	rational = false;
    } else {
	THROW("ikind must be 1 or 2.");
    }
    return new SplineCurve(cv->in, cv->ik, cv->et,
			     coefsstart, cv->idim, rational);
}

SISLSurf* GoSurf2SISL( const SplineSurface& sf, bool copy)
{
    std::vector<double>::const_iterator coefsstart;
    int ikind;
    if (sf.rational()) {
	coefsstart = sf.rcoefs_begin();
	ikind = 2;
    } else {
	coefsstart = sf.coefs_begin();
	ikind = 1;
    }
    return newSurf(sf.numCoefs_u(), sf.numCoefs_v(),
		   sf.order_u(), sf.order_v(),
		   const_cast<double*>(&(*(sf.basis_u().begin()))),
		   const_cast<double*>(&(*(sf.basis_v().begin()))),
		   const_cast<double*>(&(*coefsstart)),
		   ikind, sf.dimension(), copy);
}

SplineSurface* SISLSurf2Go( SISLSurf* sf)
{
    double* coefsstart;
    bool rational;
    if (sf->ikind == 2) {
	coefsstart = sf->rcoef;
	rational = true;
    } else if (sf->ikind == 1) {
	coefsstart = sf->ecoef;
	rational = false;
    } else {
	THROW("ikind must be 1 or 2.");
    }
    return new SplineSurface(sf->in1, sf->in2, sf->ik1, sf->ik2,
			       sf->et1, sf->et2,
			       coefsstart, sf->idim, rational);
}




} // namespace Go
