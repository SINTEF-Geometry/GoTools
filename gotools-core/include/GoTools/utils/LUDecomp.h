//===========================================================================
//                                                                           
// File: LUDecomp.h                                                      
//                                                                           
// Created: Fri Mar 11 14:16:50 2005                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: LUDecomp.h,v 1.9 2007-02-27 09:20:06 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

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

