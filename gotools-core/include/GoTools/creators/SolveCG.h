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
    /// \param relaxfac relaxation parameter
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
    double omega_;        // Relaxation parameter.

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
    /// using a given RILU preconditioner.
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
