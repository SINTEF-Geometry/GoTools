//===========================================================================
//                                                                           
// File: SplineInterpolator.C                                              
//                                                                           
// Created: Mon Oct 16 15:09:01 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: SplineInterpolator.C,v 1.26 2005-09-12 12:40:41 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/geometry/SplineInterpolator.h"

#include <vector>
#include "GoTools/utils/LUDecomp.h"
//#include "newmat.h"

using namespace std;

namespace Go
{

//===========================================================================
SplineInterpolator::~SplineInterpolator()
//===========================================================================
{
}


//===========================================================================
const BsplineBasis& SplineInterpolator::basis()
//===========================================================================
{
    return basis_;
}


//===========================================================================
void SplineInterpolator::interpolate(int num_points,
				       int dimension,
				       const double* param_start,
				       const double* data_start,
				       std::vector<double>& coefs)
//===========================================================================
{
    // Check that we have reasonable conditions set and a good param sequence
    ALWAYS_ERROR_IF(ctype_ == None, "No end conditions set.");

    for (int i = 1; i < num_points; ++i) {
	ALWAYS_ERROR_IF(param_start[i] <= param_start[i-1],
		    "Parameter sequence must be strictly increasing.");
    }
    ALWAYS_ERROR_IF(ctype_ == Hermite && (dimension != start_tangent_->dimension() ||
				      dimension != end_tangent_->dimension()),
		    "In Hermite interpolation, the end tangents must have the "
		    "same dimension as the interpolation data.");

    
    // First we make a knot vector and define the spline space
    int additional_coefs = (ctype_ == Free) ? 0 :
      (((ctype_ == NaturalAtStart && end_tangent_.get() == 0)
	|| (ctype_ == NaturalAtEnd && start_tangent_.get() == 0)) ? 1 : 2);
    int num_coefs = num_points + additional_coefs;
    int order = std::min(4, num_coefs);
    ALWAYS_ERROR_IF(num_coefs < 2,"Insufficient number of points.");
    std::vector<double> knots;
    knots.reserve(num_coefs + order);
    knots.insert(knots.end(), order, param_start[0]);
    if (ctype_ == Free) {
	knots.insert(knots.end(),
		     param_start + 2,
		     param_start + num_points - 2);
    } 
    else if (ctype_ == NaturalAtStart) {
	if (end_tangent_.get() != 0)
	    knots.insert(knots.end(), param_start + 1,
			 param_start + num_points - 1);
	else
	    knots.insert(knots.end(), param_start + 1,
			 param_start + num_points - 2);
    }
    else if (ctype_ == NaturalAtEnd) {
	if (start_tangent_.get() != 0)
	    knots.insert(knots.end(), param_start + 1,
			 param_start + num_points - 1);
	else
	    knots.insert(knots.end(), param_start + 2,
			 param_start + num_points - 1);
    }
    else { // ctype_ == Natural or Hermite
	knots.insert(knots.end(),
		     param_start + 1,
		     param_start + num_points - 1);
    }
    knots.insert(knots.end(), order, param_start[num_points-1]);
    basis_ = BsplineBasis(num_coefs, order, &knots[0]);

    // Create the interpolation matrix.
    // The first and last row (equation) depends on the boundary
    // conditions (for Hermite and Natural conditions) or are
    // nonexisting (for Free conditions, the matrix has two fewer
    // rows/equations).
    // The interpolation matrix is a tridiagonal matrix in the Free
    // and Hermite cases, and a tridiagonal matrix plus two elements
    // in the Natural case.

// This code dates from the time when this code was dependent on the NEWMAT
// library.  This is no longer the case, and the code has been substituted with
// the code below.  However, it is kept here for reference reasons.
//---------------------NEWMAT dependent-----------------------------
// #if 0
//     cerr << "NEWMAT DEPENDENT" << endl;
//     BandMatrix A(num_coefs, 2, 2);
//     A = 0;
//     ColumnVector c(num_coefs);
//     ColumnVector b(num_coefs);

//     // Boundary conditions
//     if (ctype_ == Hermite) {
// 	// Derivative conditions
// 	double tmp[8];
// 	basis_.computeBasisValues(param_start[0], tmp, 1);
// 	A.element(0, 0) = tmp[1]; // derivative of first B-spline
// 	A.element(0, 1) = tmp[3]; // derivative of second B-spline
// 	basis_.computeBasisValues(param_start[num_points-1], tmp, 1);
// 	A.element(num_coefs - 1, num_coefs - 2) = tmp[5];
// 	A.element(num_coefs - 1, num_coefs - 1) = tmp[7];
// 	// Boundary element conditions
// 	A.element(1, 0) = 1.0;
// 	A.element(num_coefs - 2, num_coefs - 1) = 1.0;
//     } else if (ctype_ == Natural) {
// 	// Derivative conditions
// 	double tmp[12];
// 	basis_.computeBasisValues(param_start[0], tmp, 2);
// 	A.element(0, 0) = tmp[2]; // second derivative of first B-spline
// 	A.element(0, 1) = tmp[5];
// 	A.element(0, 2) = tmp[8];
// 	basis_.computeBasisValues(param_start[num_points-1], tmp, 2);
// 	A.element(num_coefs - 1, num_coefs - 3) = tmp[5];
// 	A.element(num_coefs - 1, num_coefs - 2) = tmp[8];
// 	A.element(num_coefs - 1, num_coefs - 1) = tmp[11];
// 	// Boundary element conditions
// 	A.element(1, 0) = 1.0;
// 	A.element(num_coefs - 2, num_coefs - 1) = 1.0;
//     } else if (ctype_ == NaturalAtStart) {
// 	// Derivative condition
// 	double tmp[12];
// 	basis_.computeBasisValues(param_start[0], tmp, 2);
// 	A.element(0, 0) = tmp[2]; // second derivative of first B-spline
// 	A.element(0, 1) = tmp[5];
// 	A.element(0, 2) = tmp[8];
// 	if (end_tangent_.get() != 0) {
// 	    double tmp[8];
// 	    basis_.computeBasisValues(param_start[num_points-1], tmp, 1);
// 	    A.element(num_coefs - 1, num_coefs - 2) = tmp[5];
// 	    A.element(num_coefs - 1, num_coefs - 1) = tmp[7];

// 	    A.element(num_coefs - 2, num_coefs - 1) = 1.0;
// 	} else {
// 	    A.element(num_coefs - 1, num_coefs - 1) = 1.0;
// 	}
// 	// Boundary element conditions
// 	A.element(1, 0) = 1.0;
//     } else if (ctype_ == NaturalAtEnd) {
// 	// Derivative condition
// 	double tmp[12];
// 	basis_.computeBasisValues(param_start[num_points-1], tmp, 2);
// 	A.element(num_coefs - 1, num_coefs - 3) = tmp[5];
// 	A.element(num_coefs - 1, num_coefs - 2) = tmp[8];
// 	A.element(num_coefs - 1, num_coefs - 1) = tmp[11];
// 	if (start_tangent_.get() != 0) {
// 	    basis_.computeBasisValues(param_start[0], tmp, 1);
// 	    A.element(0, 0) = tmp[1]; // derivative of first B-spline
// 	    A.element(0, 1) = tmp[3]; // derivative of second B-spline

// 	    A.element(1, 0) = 1.0;
// 	} else {
// 	    A.element(0, 0) = 1.0;
// 	}
// 	// Boundary element conditions
// 	A.element(num_coefs - 2, num_coefs - 1) = 1.0;
//     } else if (ctype_ == Free) {
// 	// Boundary element conditions
// 	A.element(0, 0) = 1.0;
// 	A.element(num_coefs - 1, num_coefs - 1) = 1.0;
//     } else {
// 	THROW("Unknown boundary condition type: " << ctype_);
//     }

//     // Interior conditions
//     int rowoffset = ((ctype_ == Free) ||
// 		     ((ctype_ == NaturalAtEnd) && (start_tangent_.get() == 0)) ? 1 : 2);
//     int j;
//     for (j = 0; j < num_points-2; ++j) {
// 	double tmp[4];
// 	basis_.computeBasisValues(param_start[j+1], tmp, 0);
// 	int column = 1 + (basis_.lastKnotInterval() - 4);
// 	A.element(j + rowoffset, column) = tmp[0];
// 	A.element(j + rowoffset, column + 1) = tmp[1];
// 	A.element(j + rowoffset, column + 2) = tmp[2];
// 	A.element(j + rowoffset, column + 3) = tmp[3];
//     }

//     // Now we are ready to repeatedly solve Ac = b for # = dimension different
//     // right-hand-sides b.
//     BandLUMatrix ALUfact = A; // Computes LU factorization
//     ALWAYS_ERROR_IF(ALUfact.IsSingular(), 
// 		    "Matrix is singular! This should never happen!");

//     coefs.resize(dimension*num_coefs);

//     for (int dd = 0; dd < dimension; ++dd) {
// 	// Make the b vector
// 	int offset = (ctype_ == Free ||
// 		      ((ctype_ == NaturalAtEnd) && start_tangent_.get() == 0) ? 0 : 1);
// 	if (ctype_ == Hermite) {
// 	    b.element(0) = (*start_tangent_)[dd];
// 	    b.element(num_coefs - 1) = (*end_tangent_)[dd];
// 	} else if (ctype_ == Natural) {
// 	    b.element(0) = 0;
// 	    b.element(num_coefs - 1) = 0;
// 	} else if (ctype_ == NaturalAtStart) {
// 	    b.element(0) = 0;
// 	    if (end_tangent_.get() != 0)
// 		b.element(num_coefs - 1) = (*end_tangent_)[dd];
// 	} else if (ctype_ == NaturalAtEnd) {
// 	    b.element(num_coefs - 1) = 0;
// 	    if (start_tangent_.get() != 0)
// 		b.element(0) = (*start_tangent_)[dd];
// 	}
// 	for (j = 0; j < num_points; ++j) {
// 	    b.element(j+offset) = data_start[j*dimension + dd];
// 	}
// 	// Solve
// 	c = ALUfact.i() * b;
// 	// Copy results
// 	for (j = 0; j < num_coefs; ++j) {
// 	    coefs[j*dimension + dd] = c.element(j);
// 	}
//     }
// -------------------NEWMAT INDEPENDENT------------------------------
//#else
    vector<vector<double> > A(num_coefs, vector<double>(num_coefs, 0));
    vector<vector<double> > b(num_coefs, vector<double>(dimension));
    
    double tmp[12];
    // boundary conditions
    switch (ctype_) {
	case Hermite:
	    basis_.computeBasisValues(param_start[0], tmp, 1);
	    A[0][0] = tmp[1]; // derivative of first B-spline
	    A[0][1] = tmp[3]; // derivative of second B-spline
	    basis_.computeBasisValues(param_start[num_points-1], tmp, 1);
	    A[num_coefs - 1][num_coefs - 2] = tmp[5];
	    A[num_coefs - 1][num_coefs - 1] = tmp[7];
	    // Boundary element conditions
	    A[1][0] = 1.0;
	    A[num_coefs - 2][num_coefs - 1] = 1.0;
	    break;
	case Natural:
	    // Derivative conditions
	    basis_.computeBasisValues(param_start[0], tmp, 2);
	    A[0][0] = tmp[2]; // second derivative of first B-spline
	    A[0][1] = tmp[5];
	    A[0][2] = tmp[8];
	    basis_.computeBasisValues(param_start[num_points-1], tmp, 2);
	    A[num_coefs - 1][num_coefs - 3] = tmp[5];
	    A[num_coefs - 1][num_coefs - 2] = tmp[8];
	    A[num_coefs - 1][num_coefs - 1] = tmp[11];
	    // Boundary element conditions
	    A[1][0] = 1.0;
	    A[num_coefs - 2][num_coefs - 1] = 1.0;
	    break;
	case NaturalAtStart:
	    basis_.computeBasisValues(param_start[0], tmp, 2);
	    A[0][0] = tmp[2]; // second derivative of first B-spline
	    A[0][1] = tmp[5];
	    A[0][2] = tmp[8];
	    if (end_tangent_.get() != 0) {
		double tmp[8];
		basis_.computeBasisValues(param_start[num_points-1], tmp, 1);
		A[num_coefs - 1][num_coefs - 2] = tmp[5];
		A[num_coefs - 1][num_coefs - 1] = tmp[7];
		A[num_coefs - 2][num_coefs - 1] = 1.0;
	    } else {
		A[num_coefs - 1][num_coefs - 1] = 1.0;
	    }
	    // Boundary element conditions
	    A[1][0] = 1.0;
	    break;
	case NaturalAtEnd:
	    basis_.computeBasisValues(param_start[num_points-1], tmp, 2);
	    A[num_coefs - 1][num_coefs - 3] = tmp[5];
	    A[num_coefs - 1][num_coefs - 2] = tmp[8];
	    A[num_coefs - 1][num_coefs - 1] = tmp[11];
	    if (start_tangent_.get() != 0) {
		basis_.computeBasisValues(param_start[0], tmp, 1);
		A[0][0] = tmp[1]; // derivative of first B-spline
		A[0][1] = tmp[3]; // derivative of second B-spline
		A[1][0] = 1.0;
	    } else {
		A[0][0] = 1.0;
	    }
	    // Boundary element conditions
	    A[num_coefs - 2][num_coefs - 1] = 1.0;
	    break;
	case Free:
	    // Boundary element conditions
	    A[0][0] = 1.0;
	    A[num_coefs - 1][num_coefs - 1] = 1.0;
	    break;
	default:
	    THROW("Unknown boundary condition type." << ctype_);
    }
    
    // interior conditions
    int rowoffset = ((ctype_ == Free) ||
		     ((ctype_ == NaturalAtEnd) && (start_tangent_.get() == 0)) ? 1 : 2);
    int j;
    for (j = 0; j < num_points-2; ++j) {
	basis_.computeBasisValues(param_start[j+1], tmp, 0);
	int column = 1 + (basis_.lastKnotInterval() - 4);
	A[j + rowoffset][column] = tmp[0];
	A[j + rowoffset][column + 1] = tmp[1];
	A[j + rowoffset][column + 2] = tmp[2];
	A[j + rowoffset][column + 3] = tmp[3];
    }

    // make the b vectors boundary condition
    int offset = (ctype_ == Free ||
		  ((ctype_ == NaturalAtEnd) && start_tangent_.get() == 0) ? 0 : 1);
    switch(ctype_) {
	case Hermite:
	    copy(start_tangent_->begin(), start_tangent_->end(), b[0].begin());
	    copy(end_tangent_->begin(), end_tangent_->end(), b[num_coefs-1].begin());
	    break;
	case Natural:
	    b[0] = b[num_coefs-1] = vector<double>(dimension, 0);
	    break;
	case NaturalAtStart:
	    b[0] = vector<double>(dimension, 0);
	    if (end_tangent_.get() != 0)
		copy(end_tangent_->begin(), end_tangent_->end(), b[num_coefs-1].begin());
	    break;
	case NaturalAtEnd:
	    b[num_coefs-1] = vector<double>(dimension, 0);
	    if (start_tangent_.get() != 0)
		copy(start_tangent_->begin(), start_tangent_->end(), b[0].begin());
	    break;
	default:
	    // do nothing
	    break;
    }
    // fill in interior of the b vectors
    for (j = 0; j < num_points; ++j) {
	copy(&(data_start[j * dimension]),
	     &(data_start[(j+1) * dimension]),
	     b[j + offset].begin());
    }

    // computing the unknown vector A c = b.  b is overwritten by this unknown vector
    LUsolveSystem(A, num_coefs, &b[0]);

    // writing result to coefficients
    coefs.resize(dimension*num_coefs);
    for (int j = 0; j < num_coefs; ++j) {
	copy(b[j].begin(), b[j].end(), &coefs[j*dimension]);
    }

    //#endif
}


//===========================================================================
void SplineInterpolator::interpolate(const std::vector<double>& params,
				     const std::vector<double>& points,
				     const std::vector<int>& tangent_index,
				     const std::vector<double>& tangent_points,
				     int order,
				     std::vector<double>& coefs)
//===========================================================================
{
    makeBasis(params, tangent_index, order);
    interpolate(params, points, tangent_index, tangent_points, coefs);
}

// @@@ So far we're assuming that the user input fulfills demands wrt.
// degrees of freedom.
//===========================================================================
void SplineInterpolator::interpolate(const std::vector<double>& params,
				     const std::vector<double>& points,
				     const std::vector<int>& tangent_index,
				     const std::vector<double>& tangent_points,
				     std::vector<double>& coefs)
//===========================================================================
{
    ALWAYS_ERROR_IF(basis_set_ == false,
		"When using routine the basis_ must first be set/made.");

    int num_points = (int)params.size();
    int dimension = (int)points.size() / num_points;
    int num_coefs = basis_.numCoefs();
    int order = basis_.order();
    int tsize = (int)tangent_index.size();
    DEBUG_ERROR_IF(num_coefs != num_points + tsize,
	     "Inconsistency between size of basis and interpolation conditions.");
    DEBUG_ERROR_IF(num_coefs < order,
	     "Insufficient number of points.");

    coefs.resize(dimension*num_coefs);

    int i, j;
    // In the future there may be reason to want higher derivative information
    // in the points. Should present no problem; to be implemented when needed.



// This code dates from the time when this code was dependent on the NEWMAT
// library.  This is no longer the case, and the code has been substituted with
// the code below.  However, it is kept here for reference reasons.
//
// ---------------------newmat dependent ----------------------------
// #if 0
//     // @@@ sbr There should be performed an analysis regarding size of band matrix.
//     //    BandMatrix A(num_coefs, (order+1)/2, (order+1)/2);
//     // Ought to run through parameters (and tangents), checking against basis;
//     // then set up band matrix. This may fail. But then system should be singular.
//     // The failsafe solution would be to implement A as a full matrix...
//     BandMatrix A(num_coefs, order, order);
//     A = 0;
//     ColumnVector c(num_coefs);
//     ColumnVector b(num_coefs);

//     // Setting up the interpolation matrix A.
//     int ti = 0; // Index to first unused element of tangent_points.
//     for (i = 0; i < num_points; ++i) {
// 	bool der = ((tsize > ti) && (tangent_index[ti] == i)) ?
// 	    true : false; // true = using derivative info.
// 	std::vector<double> tmp(2*order);
// 	basis_.computeBasisValues(params[i], &tmp[0], 1);
// 	int ki = basis_.knotInterval(params[i]); // knot-interval of param.
// 	//	int column = 1 + (basis.lastKnotInterval() - 4);
// 	for (j = 0; j < order; ++j)
// 	    if ((ki-order+1+j>=0) && (ki-order+1+j<=num_coefs)) {
// 		A.element(i+ti, ki-order+1+j) = tmp[2*j];
// 		if (der)
// 		    A.element(i+ti+1, ki-order+1+j) = tmp[2*j+1];
// 	    }
// 	if (der)
// 	    ++ti;
//     }

//     // Now we are ready to repeatedly solve Ac = b for # = dimension different
//     // right-hand-sides b.
//     BandLUMatrix ALUfact = A; // Computes LU factorization
//     ALWAYS_ERROR_IF(ALUfact.IsSingular(), 
// 		    "Matrix is singular! This should never happen!");

//     coefs.resize(dimension*num_coefs);

//     for (int dd = 0; dd < dimension; ++dd) {
// 	// Make the b vector
// 	ti = 0;
// 	for (i = 0; i < num_points; ++i) {
// 	    bool der = ((tsize > ti) && (tangent_index[ti] == i)) ?
// 			true : false;
// 	    b.element(i+ti) = points[i*dimension + dd];
// 	    if (der) {
// 		b.element(i+ti+1) = tangent_points[ti*dimension + dd];
// 		++ti;
// 	    }
// 	}
// 	// Solve
// 	c = ALUfact.i() * b;
// 	// Copy results
// 	for (i = 0; i < num_coefs; ++i) {
// 	    coefs[i*dimension + dd] = c.element(i);
// 	}
//     }
//     //#else
//--------------------------- newmat independent ---------------------
    vector<vector<double> > A(num_coefs, vector<double>(num_coefs, 0));
    vector<vector<double> > b(num_coefs, vector<double>(dimension, 0));
    
    // setting up interpolation matrix A
    int ti = 0; // index to first unused element of tangent_points
    std::vector<double> tmp(2*order);
    for (i = 0; i < num_points; ++i) {
	bool der = ((tsize > ti) && (tangent_index[ti] == i)) ?
	    true : false; // true = using derivative info.
	double par = params[i];
	int ki = basis_.knotIntervalFuzzy(par); // knot-interval of param.
	basis_.computeBasisValues(params[i], &tmp[0], 1);
	//int ki = basis_.knotInterval(params[i]); // knot-interval of param.
	//	int column = 1 + (basis.lastKnotInterval() - 4);
	for (j = 0; j < order; ++j)
	    if ((ki-order+1+j>=0) && (ki-order+1+j<=num_coefs)) {
		A[i+ti][ki-order+1+j] = tmp[2*j];
		if (der)
		    A[i+ti+1][ki-order+1+j] = tmp[2*j+1];
	    }
	if (der)
	    ++ti;
    }

    // generating right-hand side
    ti = 0;
    for (i = 0; i < num_points; ++i) {
	bool der = ((tsize > ti) && (tangent_index[ti] == i)) ?
	    true : false;
	vector<double>::const_iterator pointit 
	    = points.begin() + i * dimension;
	copy(pointit, pointit + dimension, b[i+ti].begin());
	if (der) {
	    vector<double>::const_iterator tanptsit
		= tangent_points.begin() + ti * dimension;
	    copy(tanptsit, tanptsit + dimension, b[i+ti+1].begin());
	    ++ti;
	}
    }

    // Now we are ready to solve Ac = b.  b will be overwritten by solution
    LUsolveSystem(A, num_coefs, &b[0]);
    // copy results
    for (i = 0; i < num_coefs; ++i) {
	copy(b[i].begin(), b[i].end(), &coefs[i * dimension]);
    }
    //#endif
}


//===========================================================================
void SplineInterpolator::makeBasis(const std::vector<double>& params,
				   const std::vector<int>& tangent_index,
				   int order)
//===========================================================================
{
    int tsize = (int)tangent_index.size();
    int nmb_points = (int)params.size();
    DEBUG_ERROR_IF(nmb_points + tsize < order,
	     "Mismatch between interpolation conditions and order.");

    int i;
    int ti = 0; // Variable used to iterate through the tangent_index.
    int di;
    int k2 = order / 2; // Half the order.
    int kstop; // Control variable of the loop.
    int kpar;

    // Values in params are not repeated; derivatives are given by tangent_index.
    int nmb_int_cond = nmb_points + tsize;
    double startparam = params[0];
    double endparam = params[nmb_points - 1];
    // We remember start and end multiplicities.
    int start_mult = (tsize != 0 && tangent_index[0] == 0) ? 2 : 1;
    int end_mult = (tsize != 0 &&
		    tangent_index[tsize - 1] == nmb_points - 1) ?
	2 : 1;

    std::vector<double> knots(order + nmb_int_cond);
    int ki = 0; // Variabel used when setting values in the knot vector;
    // We make the knot vector k-regular in the startparam.
    for (ki = 0; ki < order; ++ki)
	knots[ki] = startparam;

    double dummy1, dummy2;
    if (order % 2 == 0) {
	// Even order: place the internal knots at the parameter values.
	dummy1 = 0.5 * (params[1] + startparam);
	dummy2 = 0.5 * (params[nmb_points - 2] + endparam);
	if (dummy1 == dummy2)
	    {
		dummy1 = (params[1] + startparam + startparam)/3.0;
		dummy2 = (params[nmb_points - 2] + endparam + endparam)/3.0;
	    }
	if (start_mult > 1)
	    ++ti;
	for (i = 0; i < start_mult - k2; ++i, ++ki)
	    knots[ki] = dummy1;

// #ifdef _MSC_VER
// 	int kss = (k2-start_mult > 0) ? k2-start_mult : 0;
// 	int kse = (k2-end_mult > 0) ? k2-end_mult : 0;
// #else
	int kss = std::max(0, k2 - start_mult);
	int kse = std::max(0, k2 - end_mult);
// #endif
	for (kpar = 1 + kss,
	       kstop = nmb_points - 1 - kse;
	     kpar < kstop; kpar++, ++ki) {
	  knots[ki] = params[kpar];
	  di = ki;
	    if (tsize != 0)
		while (tangent_index[ti] == kpar) {
		    knots[++ki] = knots[di];
		    ++ti;
		}
	}

	for (i = 0; i < end_mult - k2; ++i, ++ki)
	    knots[ki] = dummy2;
    } else {
	double delta;
	// The order is odd: place internal knots between parameter values.
	if (start_mult > 1)
	    ++ti;
	if (start_mult - k2 > 0)
	    {
		delta = (params[1] - startparam) / (start_mult - k2 + 1);
		dummy1 = startparam + delta;
		for (i = 0; i < start_mult - k2;
		     ++i, dummy1 += delta, ++ki)
		    knots[ki] = dummy1;
	    }

// #ifdef _MSC_VER
// 	int kss = (k2-start_mult > 0) ? k2-start_mult : 0;
// 	int kse = (k2-end_mult > 0) ? k2-end_mult : 0;
// #else
	int kss = std::max(0, k2 - start_mult);
	int kse = std::max(0, k2 - end_mult);
// #endif
	for (kpar = start_mult + kss,
		 kstop = nmb_points - end_mult - kse;
	     kpar < kstop; kpar++, ++ki) {
	    knots[ki] = 0.5 *(params[kpar] + params[kpar + 1]);
	    di = ki;
	    if (tsize != 0)
		while (tangent_index[ti] == kpar) {
		    knots[++ki] = knots[di];
		    ++ti;
		}
	}

	if (end_mult - k2 > 0)
	    {
// 		delta = (endparam - params[nmb_points -end_mult - 1]) /
		delta = (endparam - params[nmb_points - 2]) /
		    (end_mult - k2 + 1);
// 		dummy2 = params[nmb_points -end_mult - 1] +delta;
		dummy2 = params[nmb_points - 2] +delta;
		for (i = 0; i < end_mult - k2; ++i, dummy2 += delta, ++ki)
		    knots[ki] = dummy2;
	    }
    }

    for (ki = 0; ki < order; ki++)
	knots[nmb_int_cond+ki] = endparam;

    basis_ = BsplineBasis(order, knots.begin(), knots.end());
    basis_set_ = true;

}


} // namespace Go
