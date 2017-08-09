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

#ifndef _SOLVECG_H_
#define _SOLVECG_H_


//   -----------------------------------------------------------------------
//      Interface file for class SolveCG
//   -----------------------------------------------------------------------
//
//       Solve the equation system Ax=b where A is a symmetric
//       positive definite matrix using Conjugate Gradient Method.
//
//   -----------------------------------------------------------------------
//    Written by: Vibeke Skytt                            09.99
//    Based on  : PrCG.h written by Mike Floater
//   -----------------------------------------------------------------------

#include <vector>

namespace Go
{

/// Solve the equation system Ax=b where A is a symmetric
/// positive definite matrix using the Conjugate Gradient Method.
class SolveCG
{
public:

    /// Default constructor.
    SolveCG();

    /// Destructor.
    virtual ~SolveCG();

    /// Attach the left side of the equation system to the current
    /// object and represent the matrix as a sparse matrix. No test
    /// is applied on whether the matrix really is symmetric and positive
    /// definite.
    /// \param gmat the system matrix for the linear equations.
    ///             Size is nn*nn, stored columnwise.
    /// \param nn the number of unknowns in the system.
    void attachMatrix(double *gmat, int nn);

    /// Prepare for preconditioning.
    /// \param relaxfac relaxation parameter. Range: [0,0, 1.0].
    virtual void precondRILU(double relaxfac);

    /// Solve the equation system by conjugate gradient method.
    /// \param ex the solution vector.  The input should be the initial
    ///           guess.  Size is equal to nn.
    /// \param eb the right side of the equation. Size is equal to nn.
    /// \param nn the number of unknowns int the system.
    /// \return 0: success, 1: iterationcount exceeded, < 0: error.
    int solve(double *ex, double *eb, int nn);

    /// Set numerical tolerance used by the solver.
    /// \param tolerance numerical tolerance.
    void setTolerance(double tolerance = 1.0e-6)
    {tolerance_ = tolerance;}

    /// Set the maximal number of iterations to be used by the solver.
    /// \param max_iterations the maximal number of iterations.
    void setMaxIterations(int max_iterations)
    {max_iterations_ = max_iterations;}


protected:

    std::vector<double> A_;   // Sparse matrix containing the left side
                              // of the equation system.
    int nn_;           // Size of equation system, i.e. number of unknowns.
    int np_;           // Number of non-zero entries in the equation system.
    std::vector<int> irow_;  // The indexes in A_ and jcol_ of the
                          // first non-zeros of the nn_ rows.
    std::vector<int> jcol_;  // The np_ indexes j of the non-zero elements

    double  tolerance_; // The numerical tolerance deciding if we have reached a solution.
    int     max_iterations_; // The maximal number of iterations to be used by solver.

    // Parameters used in RILU preconditioning.

    std::vector<double> M_;  // Preconditioning matrix.
    double omega_;        // Relaxation parameter, in the range [0.0, 1.0].

    std::vector<int> diagonal_;  // Index of diagonal elements in the jcol
    int diagset_; // Whether the index of the diagonal elements has been set.

    /// Compute the matrix product sy = A_ * sx.
    /// \param sx the vector to be multiplied by the matrix.
    /// \param sy the resulting vector.
    template <typename RandomIterator1, typename RandomIterator2>
    void matrixProduct(RandomIterator1 sx, RandomIterator2 sy)
    {
	int kj, ki;
	for(kj=0; kj<nn_; kj++) {
	    sy[kj] = 0.0;
	    for(ki=irow_[kj]; ki<irow_[kj+1]; ki++) {
		sy[kj] += A_[ki] * sx[jcol_[ki]];
	    }
	}
    }

    /// Given an index in the full equation system, get the index in A_.
    int getIndex(int ki, int kj);

    /// Apply preconditioning matrix, i.e. solve the equation system
    /// M_*s = r, where M_ stores an LU-factorized matrix.
    /// \param r the input (right side) vector.
    /// \param s the output (unknown) vector.
    void forwBack(double *r, double *s);

    // Compute sy = A_^T * sx.
    void transposedMatrixProduct(double *sx, double *sy);

    /// Solve the equation system by conjugate gradient method
    /// using a given RILU (Relaxed Incomplete LU) preconditioner
    /// \param ex the solution vector.  The input should be the initial
    ///           guess.  Size is equal to nn.
    /// \param eb the right side of the equation. Size is equal to nn.
    /// \param nn the number of unknowns int the system.
    /// \return 0: success, 1: iterationcount exceeded, < 0: error.
    int solveRILU(double *ex, double *eb, int nn);

    // Solve the equation system by conjugate gradient method
    /// \param ex the solution vector.  The input should be the initial
    ///           guess.  Size is equal to nn.
    /// \param eb the right side of the equation. Size is equal to nn.
    /// \param nn the number of unknowns int the system.
    /// \return 0: success, 1: iterationcount exceeded, < 0: error.
    int solveStd(double *ex, double *eb, int nn);

    /// Print to file ("fM.m") the preconditioning matrix.
    void printPrecond();	// Print LU factorised preconditioning matrix

};

} // namespacew Go

#endif
