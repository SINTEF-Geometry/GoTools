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

#ifndef _LUDECOMP_NEW_H
#define _LUDECOMP_NEW_H

#include <vector>
#include <stdexcept>

namespace Go 
{

/** LU decomposition algorithm, based on Crout's algorithm
 * \param mat The matrix to be decomposed.  The actually used class must support the 
 *            operation [][] and return 'double'.
 * \param num_rows Number of rows (which is also equal to the number of columns).
 * \param perm Should point to an n-sized array where the permutation is written.
 * \param parity Value upon function completion is 'true' if the number of row 
 *               interchanges was pair.  'false' otherwise.
 */
//===========================================================================
template<typename SquareMatrix> 
void LUDecomp(SquareMatrix& mat, int num_rows, int* perm, bool& parity);
//===========================================================================

/** Solve the system Ax = b for x, using LU decomposition of the matrix A.
 * Upon successful completion of the function, the matrix A will be LU decomposed,
 * and the solution x will be computed and stored at the memory area pointed to by the
 * last argument.
 * \param A   The system matrix A.  The actually used class must support the operation
 *            [][] and return 'double'.
 * \param num_unknowns Number of unknowns in the equation system.
 * \param vec At function invocation, 'vec' should point to an array of T, containing 
 *            the vector 'b'.  On successful completion of the function, this array 
 *            will contain the solution for x.
 */
//===========================================================================
template<typename SquareMatrix, typename T>
void LUsolveSystem(SquareMatrix& A, int num_unknowns, T* vec);
//===========================================================================

/** Using forward substitution to calculate x on the system Lx = b, where L is a
 * lower triangular matrix with unitary diagonal.
 * \param L The lower triangular matrix.  The actually used class must support the
 *          operation [][] and return 'double'.
 * \param x At function invocation, x should point to an array of T, containing the
 *          vector 'b'.  On successful completion of the function, this array will 
 *          contain the solution for x.
 * \param num_unknowns The system size (number of unknowns).
 */
//===========================================================================
template<typename SquareMatrix, typename T>
void forwardSubstitution(const SquareMatrix& L, T* x, int num_unknowns);
//===========================================================================

/** Using backward substitution to calculate x on the system Ux = b, where U is an
 * upper triangular matrix with unitary diagonal.
 * \param U The upper triangular matrix.  The actually used class must support the
 *          operation [][] and return 'double'.
 * \param x At function invocation, x should point to an array of T, containing the 
 *          vector 'b'.  On successful completion of the function, this array will 
 *          contain the solution for x.
 * \param num_unknowns The system size (number of unknowns).
 */
//===========================================================================
template<typename SquareMatrix, typename T>
void backwardSubstitution(const SquareMatrix& U, T* x, int num_unknowns);
//===========================================================================

/** Using forward substitution to calculate x on the system Lx = b, where L is a
 * lower triangular matrix with unitary diagonal.
 * \param L The lower triangular matrix.  The actually used class must support the
 *          operation [][] and return 'double'.
 * \param x At function invocation, x should point to a vector of doubles, containing the
 *          vector 'b'.  On successful completion of the function, this vector will 
 *          contain the solution for x.
 * \param num_unknowns The system size (number of unknowns).
 */
//===========================================================================
template<typename SquareMatrix>
void forwardSubstitution(const SquareMatrix& L, std::vector<double>* x, int num_unknowns);
//===========================================================================

/** Using backward substitution to calculate x on the system Ux = b, where U is an
 * upper triangular matrix with unitary diagonal.
 * \param U The upper triangular matrix.  The actually used class must support the
 *          operation [][] and return 'double'.
 * \param x At function invocation, x should point to a vector of doubles, containing the 
 *          vector 'b'.  On successful completion of the function, this vector will 
 *          contain the solution for x.
 * \param num_unknowns The system size (number of unknowns).
 */
//===========================================================================
template<typename SquareMatrix>
void backwardSubstitution(const SquareMatrix& U, std::vector<double>* x, int num_unknowns);
//===========================================================================

}; // end namespace Go

#include "GoTools/utils/LUDecomp_implementation.h"

#endif // _LUDECOMP_NEW_H

