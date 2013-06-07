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
