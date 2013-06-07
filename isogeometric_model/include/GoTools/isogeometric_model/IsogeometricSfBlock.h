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

#ifndef __ISOGEOMETRICSFBLOCK_H
#define __ISOGEOMETRICSFBLOCK_H


#include <vector>
#include <memory>
#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/utils/Point.h"
#include "GoTools/isogeometric_model/IsogeometricBlock.h"
#include "GoTools/isogeometric_model/SfSolution.h"
#include "GoTools/isogeometric_model/SfBoundaryCondition.h"
#include "GoTools/isogeometric_model/SfPointBdCond.h"



namespace Go
{

  // This class represents one block in a block-structured isogeometric surface model
  // The block contains one spline surface describing the geometry and possibly
  // a number of associated solution spaces with boundary conditions. This 
  // information is collected in objects of type SfSolution

  class IsogeometricSfBlock : public IsogeometricBlock
  {
  public:
    // Functions used to create the block
    // Constructor
    IsogeometricSfBlock(IsogeometricModel* model,
			shared_ptr<SplineSurface> geom_sf,
			std::vector<int> solution_space_dimension,
			int index);

    // Destructor
    virtual ~IsogeometricSfBlock();

    // Cast as IsogeometricSfBlock
    virtual IsogeometricSfBlock* asIsogeometricSfBlock();

    // Multiblock. Add neighbourhood information
    // This function is called from IsogeometricSfModel and used in building the
    // complete block structured model.
    void addNeighbour(shared_ptr<IsogeometricSfBlock> neighbour,
		      int edge_nmb_this, int edge_nmb_other, bool equ_orient);

    // Count the number of neighbouring surface blocks to this block
    virtual int nmbOfNeighbours() const;

    // Return the neighbour along a specified edge. If this edge corresponds to an outer
    // boundary, a zero point is returned
    // edge_number: 0=umin, 1=umax, 2=vmin, 3=vmax
    virtual IsogeometricSfBlock* getNeighbour(int edge_nmb) const;

    // Functions used to access data
    // Given this block and another one, check if they are neighbours
    virtual bool isNeighbour(IsogeometricBlock* other) const;


    // Get the total number of coefficients in the block
    virtual int nmbCoefs() const;

    // Get B-spline basis in one paramenter direction in the block
    // The first parameter direction has pardir = 0, etc
    virtual BsplineBasis basis(int pardir) const;

    // Return the specified boundary curve, 
    // edge_number: 0=umin, 1=umax, 2=vmin, 3=vmax
    shared_ptr<SplineCurve> getGeomBoundaryCurve(int edge_number) const;

  
    // Given a point lying on a specified boundary, return the corresponding parameter
    // value of the corresponding boundary curve
    // edge_number: 0=umin, 1=umax, 2=vmin, 3=vmax
    double getParamOnBdCurve(int edge_number, const Point& position) const;

    // Fetch information about degenerate surface boundaries with respect 
    // to a given tolerance epsge
    // The edge numbers of any degenerate edges are collected
    bool geomIsDegenerate(std::vector<int>& degen_bd, double epsge);

    // Get the local enumeration of all coefficients belonging to degenerate
    // edges. The vector indexes of the degen_bd and enumeration correspond
    void getDegenEnumeration(std::vector<int>& degen_bd, 
			     std::vector<std::vector<int> >& enumeration,
			     double epsge);

    // Check if the surface block is periodic with respect to the tolerance
    // epsge
    // per[pardir] = -1 : Not closed or periodic
    //                0 : Closed or periodic with continuity of position
    //                1 : Periodic with continuity of position and derivative
    //               >1 : Higher order continuity
    bool geomIsPeriodic(int per[], double epsge);

    // Fetch enumeration along a periodic boundary in a specified parameter
    // direction. If the surface is not closed or periodic, the function 
    // returns false
    // Only the enumeration of the boundary coefficients are returned even
    // if the continuity is higher than C0 across the seam.
    bool getPeriodicEnumeration(int pardir,
				std::vector<std::pair<int, int> >& enumeration);

    // Refine the geometry surface
    // The solution spline space is refine accordingly such that the 
    // geometry space is always a sub space of the solution space
    // Insert a number of specified new knots in a specified parameter direction
    virtual void refineGeometry(std::vector<double> newknots, int pardir);

    // Refine the geometry model in a specified direction in such a way that the
    // spline space specified as input will be a sub space of the spline space of
    // geometry object. The minimum refinement to achieve this is chosen.
    virtual void refineGeometry(const BsplineBasis& other_basis, int pardir);

    // Increase degree of the geometry surface in a given parameter direction. 
    // If new_degree is not larger than the current degree, no change is performed
    virtual void increaseGeometryDegree(int new_degree, int pardir);

    // Update the current geometry surface with respect to a given boundary curve
    // If the spline space of the surface is not able to fit the new boundary
    // exactly, approximation is performed
    void updateGeometry(shared_ptr<SplineCurve> new_boundary, int edge_number);

    // Release scratch related to pre evaluated basis functions and surface.
    virtual void erasePreEvaluatedBasisFunctions();

    // 1. derivative in u-direction, 1. derivatuve in v_direction

    // Fetch boundary conditions
    // Get the number of boundary conditions attached to this block
    virtual int getNmbOfBoundaryConditions() const;

    // Get all boundary conditions related to a specified edge. Boundary conditions
    // related to all solution spaces are returned
    void 
      getEdgeBoundaryConditions(int edge_number, 
				std::vector<shared_ptr<SfBoundaryCondition> >& bd_cond) const;

    // Get number of point type bounday conditions
    virtual int getNmbOfPointBdConditions() const;

    // Get all point boundary conditions related to a specified face. Boundary 
    // conditions related to all solution spaces are returned
    void 
      getEdgePointBdConditions(int edge_number, 
			       std::vector<shared_ptr<SfPointBdCond> >& bd_cond) const;

    // Get specified solution space
    shared_ptr<SfSolution> getSolutionSpace(int solution_index);

    // Get geometry surface
    shared_ptr<SplineSurface> surface() const;

    // Ensure minimum degree of solution space
    // The solution space will always have at least the degree of the
    // corresponding geometry surface
    virtual void setMinimumDegree(int degree, int solutionspace_idx);

    // Update spline spaces of the solution to ensure consistency
    // Returns true if any update occured, false if not
    // Solution space index is a global value valid for all blocks in a model.
    virtual bool updateSolutionSplineSpace(int solutionspace_idx);

    // Get number of solution spaces
    int nmbSolutionSpaces() const;

    // The edge position of a boundary curve
    // Returns -1 if not possible to determine edge position. Otherwise returns
    //   4 if parameter direction is v-direction
    // + 2 if edge is constant parameter u=u_max or v=v_max
    // + 1 if orientation is reversed with respect to orientation on surface
    int getEdgeOrientation(shared_ptr<ParamCurve> crv, double tol);

    // Store list of neighbouring information between this block and another
    // edges will hold the edge position on this block for each match (0, 1, 2, 3
    // for edge u_min, u_max, v_min, v_max respectiveliy).
    // edges_other will hold the corresponding edge positions for the other block
    // equal_oriented will hold whether the surfaces are equally orineted at each match
    void getNeighbourInfo(IsogeometricSfBlock* other,
			  std::vector<int>& edges,
			  std::vector<int>& edges_other,
			  std::vector<bool>& equal_oriented);

  private:

    // The surface describing the geometry
    shared_ptr<SplineSurface> surface_;

    // The position index in the Model object
    int index_;

    // Solution spaces
    std::vector<shared_ptr<SfSolution> > solution_;

    // Adjacent blocks
    // neighbours_[0] = adjacent block along edge u = u_min
    // neighbours_[1] = adjacent block along edge u = u_max
    // neighbours_[2] = adjacent block along edge v = v_min
    // neighbours_[3] = adjacent block along edge v = v_max
    shared_ptr<IsogeometricSfBlock> neighbours_[4];

    // Adjacency position at other block:
    // neighb_edge_[i] = 0 : neighbours_[i] has this block as neighbour along edge u = u_min
    // neighb_edge_[i] = 1 : neighbours_[i] has this block as neighbour along edge u = u_max
    // neighb_edge_[i] = 2 : neighbours_[i] has this block as neighbour along edge v = v_min
    // neighb_edge_[i] = 3 : neighbours_[i] has this block as neighbour along edge v = v_max
    int neighb_edge_[4];

    // Indicate whether this block and neighbour are equally oriented along boundary
    bool equal_orientation_[4];

  };    // end class IsogeometricSfBlock

} // end namespace Go


#endif    // #ifndef __ISOGEOMETRICSFBLOCK_H
