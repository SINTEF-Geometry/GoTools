//===========================================================================
//                                                                           
// File: AlgObj2DInt.C                                                       
//                                                                           
// Created: Fri Jan 28 13:29:29 2005                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: AlgObj2DInt.C,v 1.7 2006-03-08 09:31:20 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/intersections/AlgObj2DInt.h"
#include "GoTools/utils/errormacros.h"


using std::vector;


namespace Go {


//===========================================================================
AlgObj2DInt::AlgObj2DInt(int degree)
    : degree_(degree), power_basis_(false)
//===========================================================================
{
}


//===========================================================================
AlgObj2DInt::AlgObj2DInt(const vector<Alg2DElem>& terms)
    : terms_(terms), power_basis_(true)
//===========================================================================
{
    // We must run through terms_ to find the degree_.
    int degree = 0;
    for (size_t ki = 0; ki < terms_.size(); ++ki) {
	int sum_deg = terms_[ki].degrees_[0] +
	    terms_[ki].degrees_[1];
	if (sum_deg > degree) {
	    degree = sum_deg;
	}
    }

    degree_ = degree;
}


//===========================================================================
AlgObj2DInt::~AlgObj2DInt()
//===========================================================================
{
}


//===========================================================================
Alg2DElem AlgObj2DInt::term(int index)
//===========================================================================
{
    ASSERT(index < int(terms_.size()));
    return terms_[index];
}


} // namespace Go
