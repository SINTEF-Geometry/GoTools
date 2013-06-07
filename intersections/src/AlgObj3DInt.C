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
