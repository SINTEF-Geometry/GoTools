//===========================================================================
//                                                                           
// File: ConstraintDefinitions.h                                             
//                                                                           
// Created: Thu Sep 22 10:08:12 2005                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: ConstraintDefinitions.h,v 1.1 2005-09-22 08:13:09 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _CONSTRAINTDEFINITIONS_H
#define _CONSTRAINTDEFINITIONS_H

#include <vector>


namespace Go
{

    /// Struct defining linear side constraints between control points in
    /// a surface.
    /// The elements of factor_ correspond to a linear combination of
    /// control points on the left side of the equation, whilst
    /// constant_term_ denotes the right side of the equation.
    typedef struct sideConstraint
    {
	int dim_;  // Dimension of coefficient (max 3).
	// For each coefficient involved in the constraint, the 
	// index of the coefficient is given and the factor corresponding
	// to the coefficient in the equation.
	std::vector<std::pair<int, double> > factor_;
	double constant_term_[3];  // The constant term in the current equation.
    } sideConstraint;


    /// Struct defining linear side constraints between control points in
    /// a set of surfaces.
    /// The elements of factor_ correspond to a linear combination of
    /// control points on the left side of the equation, whilst
    /// constant_term_ denotes the right side of the equation.
    typedef struct sideConstraintSet
    {
	int dim_;  // Dimension of coefficient (max 3).
	// For each coefficient involved in the constraint, the index
	// of the surface, the index of the coefficient and the 
	// factor corresponding to the coefficient in the equation is given.
	std::vector<std::pair<std::pair<int,int>, double> > factor_; 
	// The constant term in the current equation, on the right side of the
	//  equation. For dim < 3 not all elements are used.
	double constant_term_[3];
	/// Default constructor.
	sideConstraintSet()
	{
	    dim_ = 3;
	    constant_term_[0] = constant_term_[1] = constant_term_[2] = 0.0;
	}

	/// Constructor.
	/// \param dim dimension of geometric space.
	sideConstraintSet(int dim)
	{
	    dim_ = dim;
	    constant_term_[0] = constant_term_[1] = constant_term_[2] = 0.0;
	}
    } sideConstraintSet;

}

#endif // _CONSTRAINTDEFINITIONS_H

