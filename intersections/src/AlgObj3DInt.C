//===========================================================================
//                                                                           
// File: AlgObj3DInt.C                                                       
//                                                                           
// Created: Fri Jan 28 13:29:16 2005                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: AlgObj3DInt.C,v 1.8 2006-05-04 12:19:00 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/intersections/AlgObj3DInt.h"
#include "GoTools/utils/errormacros.h"


using std::vector;


namespace Go
{


//===========================================================================
AlgObj3DInt::AlgObj3DInt(int degree)
    : degree_(degree), power_basis_(false)
//===========================================================================
{
}


//===========================================================================
AlgObj3DInt::AlgObj3DInt(const vector<Alg3DElem>& terms)
    : terms_(terms), power_basis_(true)
//===========================================================================
{
    // We must run through terms_ to find the degree_.
    int degree = 0;
    for (size_t ki = 0; ki < terms_.size(); ++ki) {
	int sum_deg = terms_[ki].degrees_[0] +
	    terms_[ki].degrees_[1] + terms_[ki].degrees_[2];
	if (sum_deg > degree) {
	    degree = sum_deg;
	}
    }

    degree_ = degree;

}


//===========================================================================
AlgObj3DInt::AlgObj3DInt(const BernsteinTetrahedralPoly& implicit,
			 const BaryCoordSystem3D& bc)
    : implicit_(implicit), bc_(bc)
//===========================================================================
{
    degree_ = implicit.degree();
}


//===========================================================================
AlgObj3DInt::~AlgObj3DInt()
//===========================================================================
{
}


//===========================================================================
Alg3DElem AlgObj3DInt::term(int index)
//===========================================================================
{
    ASSERT(index < int(terms_.size()));
    return terms_[index];
}


} // end namespace Go
