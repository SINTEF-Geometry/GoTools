//===========================================================================
//
// File : BlockSolution.h
//
// Created: November, 2009
//
// Author: Vibeke Skytt, SINTEF, and Anh-Vu Vuong, TU Munich
//
// Revision: 
//
// Description:
//
//===========================================================================



#ifndef __BLOCKSOLUTION_H
#define __BLOCKSOLUTION_H


#include <vector>
#include "GoTools/utils/Point.h"
#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/topology/tpTopologyTable.h"



namespace Go
{

  class SfSolution;
  class VolSolution;

  // This class represents one solution in one block in a block-structured
  // isogeometric  model

  class BlockSolution
  {
  public:

    // Destructor
    virtual ~BlockSolution();

    // Cast as SfSolution, default is NULL
    virtual SfSolution* asSfSolution();

    // Cast as VolSolution, default is NULL
    virtual VolSolution* asVolSolution();

    // Get number of point type bounday conditions
    virtual int getNmbOfPointBdConditions() const = 0;
  

    // Given this block and its neighbour, check if the spline spaces matches
    virtual bool matchingSplineSpace(BlockSolution* other) const = 0;

    // Given this block and a neighbour, give a vector of corresponding
    // coefficient enumerations in the solutionspace
    // If the spline spaces don't match, no enumeration is returned
    virtual void 
      getMatchingCoefficients(BlockSolution* other,
			      std::vector<std::pair<int,int> >& enumeration,
			      int match_pos = 0) const = 0;

    // Given a block and one of its boundaries, give a vector of coefficient enumerations
    // in the solutionspace along the boundary
    virtual void
      getBoundaryCoefficients(int boundary,
			      std::vector<int>& enumeration) const = 0;

    // Given this block and its neighbour, make the spline spaces match
    virtual void makeMatchingSplineSpace(BlockSolution* other) = 0;

    // Functions used to modify data related to solutions
    // Increase degree in a given parameter direction. If new_degree is not larger 
    // than the current degree, no change is performed
    virtual void increaseDegree(int new_degree, int pardir) = 0;

    // Refine spline space of solution in a given parameter direction
    // Specify the knot interval, the knot is inserted in the middle
    virtual void insertKnots(const std::vector<int>& knot_intervals, int pardir) = 0;

    // Specify new knots. This function is not recommended in a multi-block
    // model
    virtual void insertKnots(const std::vector<double>& knots, int pardir) = 0;

    // Release scratch related to pre evaluated basis functions and surface.
    virtual void erasePreEvaluatedBasisFunctions() = 0;

    // Pre evaluate basis functions and derivative of basis functions to provide
    // input for the numerical integration
    // The gauss parameters with respect to the surface parameterization in both
    // parameter directions is input to the function
    // NB! Quadrature points should never lie on multiple knots if the multiplicity
    // is so high that it creates a C0 surface
    // NB! Refinement of the spline space or degree elvation will imply that the
    // pre evaluated values are removed, and this function must be called again
    virtual void performPreEvaluation(std::vector<std::vector<double> >& Gauss_par) = 0;

    // Return the value of the Jacobian determinant of the parameterization in a specified Gauss point.
    // Requires pre evaluation to be performed
    virtual double getJacobian(std::vector<int>& index_of_Gauss_point) const = 0;

    // Fetch position and derivatives of the geometry surface in a specified Gauss point
    // Requires pre evaluation to be performed
    virtual void valuesInGaussPoint(const std::vector<int>& index_of_Gauss_point,
				    std::vector<Point>& derivs) const = 0;  // Position, 

    // Attach coefficient information to specified solution
    virtual void setSolutionCoefficients(const std::vector<double>& coefs) = 0;

    // Get the total number of coefficients of the solution
    virtual int nmbCoefs() const = 0;

    // Get number of coefficients in one parameter direction of the solution
    // The first parameter direction has pardir = 0, etc
    virtual int nmbCoefs(int pardir) const = 0;

    // Get polynomial degree in one parameter direction of the solution
    // The first parameter direction has pardir = 0, etc
    virtual int degree(int pardir) const = 0;

    // Get knot vector of spline space in one paramenter direction of the solution
    // The first parameter direction has pardir = 0, etc
    virtual std::vector<double> knots(int pardir) const = 0;

    // Get vector of distinct knots in spline space in one paramenter direction of the solution
    // The first parameter direction has pardir = 0, etc
    virtual std::vector<double> distinctKnots(int pardir) const = 0;

    // Get B-spline basis in one paramenter direction of the solution
    // The first parameter direction has pardir = 0, etc
    virtual BsplineBasis basis(int pardir) const = 0;

    // Get dimension of solution space
    virtual int dimension() const = 0;

    // Update the conditions
    virtual void updateConditions() = 0;

    // Get tolerances
    virtual tpTolerances getTolerances() const = 0;

  };    // end class BlockSolution

} // end namespace Go


#endif    // #ifndef __BLOCKSOLUTION_H
