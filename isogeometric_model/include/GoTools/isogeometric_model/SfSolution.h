//===========================================================================
//
// File : SfSolution.h
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



#ifndef __SFSOLUTION_H
#define __SFSOLUTION_H


#include <vector>
#include <memory>
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/utils/Point.h"
#include "GoTools/isogeometric_model/BlockSolution.h"
#include "GoTools/isogeometric_model/SfBoundaryCondition.h"
#include "GoTools/isogeometric_model/SfPointBdCond.h"


namespace Go
{

  class IsogeometricSfBlock;

  struct preEvaluationSf
  {
    // Storage for input to and results of the grid evaluation for basis functions
    std::vector<double> gauss_par1_;  // Gauss points in 1. parameter direction
    std::vector<double> gauss_par2_;  // Gauss points in 2. parameter direction
    std::vector<double> basisvals_u_; // Non-zero basis functions and derivatives thereof, 1. par.dir
    std::vector<double> basisvals_v_; // Non-zero basis functions and derivatives thereof, 2. par.dir
    std::vector<int>    left_u_;      // Index of first non-zero basis function in 1. par. dir.
    std::vector<int>    left_v_;      // Index of first non-zero basis function in 2. par. dir.
  
    // Storage for grid evaluation of the geometry surface
    std::vector<double> points_;   // Position of the surface in the Gauss points
    std::vector<double> deriv_u_;  // 1. derivative of the surface in 1. par. dir. in the Gauss points
    std::vector<double> deriv_v_;  // 1. derivative of the surface in 2. par. dir. in the Gauss points
  };

  // This class represents one solution in one block in a block-structured
  // isogeometric surface model

  class SfSolution : public BlockSolution
  {
  public:

    // Constructor
    SfSolution(IsogeometricSfBlock* parent, shared_ptr<SplineSurface> sol_surf);

    // Destructor
    virtual ~SfSolution();

    // Cast as SfSolution
    virtual SfSolution* asSfSolution();

    // Add the boundary condition relevant for this block in the model
    // This function is called from IsogeometricSfModel
    void addBoundaryCondition(int edge_nmb, BdConditionType type, BdCondFunctor *fbd,
			      std::pair<double, double> end_par);

    // Add a constant value boundary condition relevant for this block in the model
    // This function is called from IsogeometricSfModel
    void addBoundaryCondition(int edge_nmb, BdConditionType type, const Point& const_val,
			      std::pair<double, double> end_par);

    // Add a point condition of Dirichlet type
    // This function is called from IsogeometricSfModel
    void addDirichletPointBdCond(double param[],
				 Point& condition_value);

    // Number of boundary conditions corresponding to this solution
    int getNmbOfBoundaryConditions() const;

    // Get a specified boundary conditions
    shared_ptr<SfBoundaryCondition> getBoundaryCondition(int index) const;

    // Get all boundary conditions related to a specified edge. 
    void 
      getEdgeBoundaryConditions(int edge_number, 
				std::vector<shared_ptr<SfBoundaryCondition> >& bd_cond) const;

    // Get all boundary conditions
    void 
      getEdgeBoundaryConditions(std::vector<shared_ptr<SfBoundaryCondition> >& bd_cond) const;

    // Get number of point type bounday conditions
    virtual int getNmbOfPointBdConditions() const;
  
    // Get a specified point type boundary condition
    shared_ptr<SfPointBdCond> getPointBdCondition(int index) const;

    // Get all point boundary conditions related to a specified edge. 
    void 
      getEdgePointBdConditions(int edge_number, 
			       std::vector<shared_ptr<SfPointBdCond> >& bd_cond) const;

    // Get all point boundary conditions 
    void 
      getPointBdCond(std::vector<shared_ptr<SfPointBdCond> >& bd_cond) const;

    // Given this block and its neighbour, check if the spline spaces matches
    virtual bool matchingSplineSpace(BlockSolution* other) const;

    // Given this block and a neighbour, give a vector of corresponding
    // coefficient enumerations in the solutionspace
    // match_pos tells which edge match to use (if the blocks are 'multineighbours')
    // If the spline spaces don't match, no enumeration is returned
    // If the blocks have less than match_pos+1 neighbour relations, no enumeration is returned
    virtual void 
      getMatchingCoefficients(BlockSolution* other,
			      std::vector<std::pair<int,int> >& enumeration,
			      int match_pos = 0) const;

    // Given a block and one of its boundaries, give a vector of coefficient enumerations
    // in the solutionspace along the boundary
    virtual void
      getBoundaryCoefficients(int boundary,
			      std::vector<int>& enumeration) const;

    // Given this block and its neighbour, make the spline spaces match
    virtual void makeMatchingSplineSpace(BlockSolution* other);

    // Functions used to modify data related to solutions
    // Increase degree in a given parameter direction. If new_degree is not larger 
    // than the current degree, no change is performed
    virtual void increaseDegree(int new_degree, int pardir);

    // Refine spline space of solution in a given parameter direction
    // Specify the knot interval, the knot is inserted in the middle
    virtual void insertKnots(const std::vector<int>& knot_intervals, int pardir);

    // Specify new knots. This function is not recommended in a multi-block
    // model
    virtual void insertKnots(const std::vector<double>& knots, int pardir);

    // Release scratch related to pre evaluated basis functions and surface.
    virtual void erasePreEvaluatedBasisFunctions();

    // Pre evaluate basis functions and derivative of basis functions to provide
    // input for the numerical integration
    // The gauss parameters with respect to the surface parameterization in both
    // parameter directions is input to the function
    // NB! Quadrature points should never lie on multiple knots if the multiplicity
    // is so high that it creates a C0 surface
    // NB! Refinement of the spline space or degree elvation will imply that the
    // pre evaluated values are removed, and this function must be called again
    virtual void performPreEvaluation(std::vector<std::vector<double> >& Gauss_par);

    // Get value and 1. derivative of all non-zero rational basis funtions
    // in the given Gauss point
    // Requires pre evaluation to be performed
    void getBasisFunctions(int index_of_Gauss_point1,
			   int index_of_Gauss_point2,
			   std::vector<double>& basisValues,
			   std::vector<double>& basisDerivs_u,
			   std::vector<double>& basisDerivs_v) const;

    // Not recommended, but provided if you really want it
    // Get value and 1. derivative of all non-zero rational basis funtions
    // in the given parameter value
    void getBasisFunctions(double param1,
			   double param2,
			   std::vector<double>& basisValues,
			   std::vector<double>& basisDerivs_u,
			   std::vector<double>& basisDerivs_v) const;

    // Return the value of the Jacobian determinant of the parameterization in a specified Gauss point.
    // Requires pre evaluation to be performed
    virtual double getJacobian(std::vector<int>& index_of_Gauss_point) const;

    // Fetch position and derivatives of the geometry surface in a specified Gauss point
    // Requires pre evaluation to be performed
    virtual void valuesInGaussPoint(const std::vector<int>& index_of_Gauss_point,
				    std::vector<Point>& derivs) const;  // Position, 

    // Attach coefficient information to specified solution
    virtual void setSolutionCoefficients(const std::vector<double>& coefs);

    // Return the surface representing a specified solution
    shared_ptr<SplineSurface> getSolutionSurface() const;

    // Return the geometry surface of the block
    shared_ptr<SplineSurface> getGeometrySurface() const;

    // Ensure minimum degree of solution space
    // The solution space will always have at least the degree of the
    // corresponding geometry surface
    void setMinimumDegree(int degree);

    // Refine the solution space in specified parameter direction
    // by a minimal refinement such that the geometry spline space
    // is a subspace of the solution spline space
    void refineToGeometry(int pardir);

    // Get the total number of coefficients in the block
    virtual int nmbCoefs() const;

    // Get number of coefficients in one parameter direction of the solution
    // The first parameter direction has pardir = 0, etc
    virtual int nmbCoefs(int pardir) const;

    // Get polynomial degree in one parameter direction of the solution
    // The first parameter direction has pardir = 0, etc
    virtual int degree(int pardir) const;

    // Get knot vector of spline space in one paramenter direction of the solution
    // The first parameter direction has pardir = 0, etc
    virtual std::vector<double> knots(int pardir) const;

    // Get vector of distinct knots in spline space in one paramenter direction of the solution
    // The first parameter direction has pardir = 0, etc
    virtual std::vector<double> distinctKnots(int pardir) const;

    // Get B-spline basis in one paramenter direction in the block
    // The first parameter direction has pardir = 0, etc
    virtual BsplineBasis basis(int pardir) const;

    // Get dimension of solution space
    virtual int dimension() const;

    // Update the conditions
    virtual void updateConditions();

    // Get parameter value of Gauss point.
    // pardir is either 0 for u-direction or 1 for v-direction
    double getGaussParameter(int index_of_Gauss_point, int pardir) const;

    // Get tolerances
    virtual tpTolerances getTolerances() const;

  private:
    // The solution stored as a surface
    // The coefficients are unknowns until the numerical analysis is performed
    shared_ptr<SplineSurface>  solution_;

    // Boundary conditions associated to this solution
    std::vector<shared_ptr<SfBoundaryCondition> > boundary_conditions_;

    // Point type boundary conditions associated to this solution
    std::vector<shared_ptr<SfPointBdCond> > point_bd_cond_;

    // Storage for input to and results of the grid evaluation for basis functions
    // and geometry
    shared_ptr<preEvaluationSf> evaluated_grid_;

    // Pointer to the block to which this boundary condition belongs
    IsogeometricSfBlock* parent_;

    // Get information about alle neighbouring edges to another
    // block solution, and if the spline spaces match
    void neighbourInfo(BlockSolution* other, std::vector<int>& edges, std::vector<int>& edges_other,
		       std::vector<bool>& equal_oriented, std::vector<bool>& space_matches) const;

  };   // end class SfSolution

} // end namespace Go


#endif    // #ifndef __SFSOLUTION_H
