//===========================================================================
//
// File : IsogeometricBlock.h
//
// Created: 30th of October, 2009
//
// Author: Vibeke Skytt, SINTEF  and Anh-Vu Vuong, TU Munich
//
// Revision: 
//
// Description:
//
//===========================================================================


#ifndef __ISOGEOMETRICBLOCK_H
#define __ISOGEOMETRICBLOCK_H


#include <vector>
#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/isogeometric_model/IsogeometricModel.h"




namespace Go
{

  class IsogeometricSfBlock;
  class IsogeometricVolBlock;

  // Abstract top class represening one block in a block-structured isogeometric model
  class IsogeometricBlock
  {
  public:
    // Constructor
    IsogeometricBlock(IsogeometricModel* model)
      : model_(model)
    {
    }

    // Destructor
    virtual ~IsogeometricBlock();

    // Return model object
    virtual IsogeometricModel* model()
    { return model_; }

    // Cast as IsogeometricSfBlock, default is NULL
    virtual IsogeometricSfBlock* asIsogeometricSfBlock();

    // Cast as IsogeometricVolBlock, default is NULL
    virtual IsogeometricVolBlock* asIsogeometricVolBlock();

    // Count the number of neighbouring volume blocks to this block
    virtual int nmbOfNeighbours() const = 0;

    // Return the neighbour along a specified boundary. If this boundary corresponds 
    // to an outer boundary, a zero point is returned
    virtual IsogeometricBlock* getNeighbour(int bd_nmb) const = 0;

    // Functions used to access data
    // Given this block and another one, check if they are neighbours
    virtual bool isNeighbour(IsogeometricBlock* other) const = 0;

    // Get the total number of coefficients in the block
    virtual int nmbCoefs() const = 0;

    // Get number of coefficients in one parameter direction in the block
    // The first parameter direction has pardir = 0, etc
    virtual int nmbCoefs(int pardir) const;

    // Get polynomial degree in one parameter direction in the block
    // The first parameter direction has pardir = 0, etc
    virtual int degree(int pardir) const;

    // Get knot vector of spline space in one paramenter direction in the block
    // The first parameter direction has pardir = 0, etc
    virtual std::vector<double> knots(int pardir) const;

    // Get vector of distinct knots in spline space in one paramenter direction in the block
    // The first parameter direction has pardir = 0, etc
    virtual std::vector<double> distinctKnots(int pardir) const;

    // Get B-spline basis in one paramenter direction in the block
    // The first parameter direction has pardir = 0, etc
    virtual BsplineBasis basis(int pardir) const = 0;

    // Refine the geometry volume
    // The solution spline space is refine accordingly such that the 
    // geometry space is always a sub space of the solution space
    // Insert a number of specified new knots in a specified parameter direction
    virtual void refineGeometry(std::vector<double> newknots, int pardir) = 0;

    // Refine the geometry model in a specified direction in such a way that the
    // spline space specified as input will be a sub space of the spline space of
    // geometry object. The minimum refinement to achieve this is chosen.
    virtual void refineGeometry(const BsplineBasis& other_basis, int pardir) = 0;

    // Increase degree of the geometry volume in a given parameter direction. 
    // If new_degree is not larger than the current degree, no change is performed
    virtual void increaseGeometryDegree(int new_degree, int pardir) = 0;


    // Release scratch related to pre evaluated basis functions and surface.
    virtual void erasePreEvaluatedBasisFunctions() = 0;

			
    // Fetch boundary conditions
    // Get the number of boundary conditions attached to this block
    virtual int getNmbOfBoundaryConditions() const = 0;

    // Get number of point type bounday conditions
    virtual int getNmbOfPointBdConditions() const = 0;

    // Get tolerances
    virtual tpTolerances getTolerances() const;

  private:

    IsogeometricModel* model_;    // The model object

  };    // end class IsogeometricBlock

} // end namespace Go


#endif    // #ifndef __ISOGEOMETRICBLOCK_H
