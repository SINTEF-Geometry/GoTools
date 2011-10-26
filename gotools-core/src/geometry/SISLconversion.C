//===========================================================================
//                                                                           
// File: SISLconversion.C                                                  
//                                                                           
// Created: Mon Nov 27 16:58:39 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: SISLconversion.C,v 1.2 2005-11-09 10:22:39 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/geometry/SISL_code.h"

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
