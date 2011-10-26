//===========================================================================
//                                                                           
// File: SplineApproximator.C                                              
//                                                                           
// Created: Wed Oct 18 14:15:21 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: SplineApproximator.C,v 1.11 2005-06-30 14:03:37 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/geometry/SplineApproximator.h"
#include "GoTools/utils/LUDecomp.h"
#include <vector>
#include <math.h>

using namespace std;

namespace Go
{

//===========================================================================
SplineApproximator::~SplineApproximator()
//===========================================================================
{
}


//===========================================================================
const BsplineBasis& SplineApproximator::basis()
//===========================================================================
{
    return basis_;
}


//===========================================================================
void SplineApproximator::interpolate(int num_points,
				       int dimension,
				       const double* param_start,
				       const double* data_start,
				       std::vector<double>& coefs)
//===========================================================================
{
    bool have_basis = false;
    if (basis_.numCoefs() >= 4 && basis_.order() == 4) {
	// We have been given an OK spline space
	have_basis = true;
	num_coefs_ = basis_.numCoefs();
    }

    // Check that we have reasonable conditions set and a good param sequence
    ALWAYS_ERROR_IF(num_coefs_ < 4,
		    "Please set the number of coefficients to 4 or more.");
    ALWAYS_ERROR_IF((num_points < 4) || (num_coefs_ > num_points),
		    "Insufficient number of points.");

    int i;
    for (i = 1; i < num_points; ++i) {
	ALWAYS_ERROR_IF(param_start[i] <= param_start[i-1],
			"Parameter sequence must be strictly increasing.");
    }

    // If we do not already have a basis, we must make one
    // Here, a uniform parametrization is hardcoded.
    int order = 4;
    if (!have_basis) {
	std::vector<double> knots;
	knots.reserve(num_coefs_ + order);
	double first = param_start[0];
	double last = param_start[num_points - 1];
	knots.insert(knots.end(), order, first);
	int numintknots = num_coefs_ - order;
	if (numintknots > 0) {
	    double inc = (last - first)/static_cast<double>(numintknots+1);
	    for (i = 0; i < numintknots; ++i) {
		knots.push_back(first + inc * (i + 1));
	    }
	}
	knots.insert(knots.end(), order, last);
	basis_ = BsplineBasis(num_coefs_, order, &knots[0]);
    }
    // make the approximation matrices
    // we are searching for a c such that ||b - Ac|| is minimised

    // make A
    vector<vector<double> > A(num_points, vector<double>(num_coefs_, 0));
    vector<double> weights(num_points, 1.0);
    int j;
    for (j = 0; j < num_points; ++j) {
	double tmp[4];
	basis_.computeBasisValues(param_start[j], tmp, 0);
	int column = basis_.lastKnotInterval() - order + 1;
	A[j][column] = weights[j]*tmp[0];
	A[j][column + 1] = weights[j]*tmp[1];
	A[j][column + 2] = weights[j]*tmp[2];
	A[j][column + 3] = weights[j]*tmp[3];
    }

    // make A into AtA (this could be optimised later)
    int row, col, k;
    vector<vector<double> > AtA(num_coefs_, vector<double>(num_coefs_, 0));
    for (row = 0; row < num_coefs_; ++row) {
	for (col = 0; col < num_coefs_; ++col) {
	    for (k = 0; k < num_points; ++k) {
		AtA[row][col] += A[k][row] * A[k][col];
	    }
	}
    }

    // make vector b (each component contains 'dimension' values)
    vector<vector<double> > b(num_points, vector<double>(dimension, 0));
    for (j = 0; j < num_points; ++j) {
	for (int dd = 0; dd < dimension; ++dd) {
	    b[j][dd] = data_start[j * dimension + dd];
	}
    }

    // make solution vector c
    vector< vector<double> > c(num_coefs_, vector<double>(dimension, 0));
    // filling solution vector with At*b, prior to solving the system AtA*c =At*b
    for (row = 0; row < num_coefs_; ++row) {
	for (col = 0; col < dimension; ++col) {
	    for (k = 0; k < num_points; ++k) {
		c[row][col] += A[k][row] * b[k][col];
	    }
	}
    }

    // solve for c
    LUsolveSystem(AtA, num_coefs_, &c[0]);

    // copy the data to coefs
    coefs.resize(dimension * num_coefs_);
    for (i = 0; i < num_coefs_; ++i) {
	for (int dd = 0; dd < dimension; ++dd) {
	    coefs[i * dimension + dd] = c[i][dd];
	    if (fabs(coefs[i * dimension + dd]) < 1e-14) {
		coefs[i * dimension + dd] = 0.0;
	    }
	}
    }

}


} // namespace Go

// Below is code dating from the time when this library depended on NEWMAT.
// This is no longer the case, but the code is kept here for future reference.

// // ------------- newmat dependent ------------------------
// #if 0
//     // Make the approximation matrices.
//     // We are searching for a c such that ||b - A*c|| is minimized.
//     Matrix A(num_points, num_coefs_);
//     A = 0;
//     Matrix c(num_coefs_, dimension);
//     c = 0;
//     Matrix b(num_points, dimension);
//     b = 0;

//     // Make A
//     // We will (for now) make a weight vector here:
//     std::vector<double> weights(num_points, 1.0);
//     if (!free_edges_) {
// 	weights.front() = 1.0;
// 	weights.back() = 1.0;
//     }
//     int j;
//     for (j = 0; j < num_points; ++j) {
// 	double tmp[4];
// 	basis_.computeBasisValues(param_start[j], tmp, 0);
// 	int column = basis_.lastKnotInterval() - order + 1;
// 	A.element(j, column) = weights[j]*tmp[0];
// 	A.element(j, column + 1) = weights[j]*tmp[1];
// 	A.element(j, column + 2) = weights[j]*tmp[2];
// 	A.element(j, column + 3) = weights[j]*tmp[3];
//     }

//     // Make b
//     for (j = 0; j < num_points; ++j) {
// 	for (int dd = 0; dd < dimension; ++dd) {
// 	    b.element(j, dd) = data_start[j*dimension + dd];
// 	}
//     }

//     // Solve for c
//     Matrix AtA = A.t()*A;
//     CroutMatrix AtALUfact = AtA;
//     ALWAYS_ERROR_IF(AtALUfact.IsSingular(), 
// 		    "Matrix is singular! This should never happen!");
//     c = AtALUfact.i()*A.t()*b;

//     // Copy the data to coefs
//     coefs.resize(dimension*num_coefs_);
//     for (i = 0; i < num_coefs_; ++i) {
// 	for (int dd = 0; dd < dimension; ++dd) {
// 	    coefs[i*dimension + dd] = c.element(i, dd);
// 	    if (fabs(coefs[i*dimension + dd]) < 1e-14)
// 		coefs[i*dimension + dd] = 0.0;
// 	}
//     }

// // ---------------- newmat independent ----------------------
// #else    
