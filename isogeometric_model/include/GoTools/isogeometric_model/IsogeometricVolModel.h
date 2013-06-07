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

#ifndef __ISOGEOMETRICVOLMODEL_H
#define __ISOGEOMETRICVOLMODEL_H


#include <vector>
#include <memory>
#include "GoTools/utils/Point.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/isogeometric_model/IsogeometricModel.h"
#include "GoTools/isogeometric_model/BdConditionType.h"
#include "GoTools/isogeometric_model/IsogeometricVolBlock.h"



namespace Go
{

  // This class is intended for storage and service functionality
  // for one isogeometric model of type structured multi-block volume model.
  // It contains geometry, boundary conditions, and the corresponding
  // solutions space(s)

  class IsogeometricVolModel : public IsogeometricModel
  {
  public:
    // Constructor
    // The geometry is described as a volume model, i.e. a volume set
    // and adjacency information. All volumes should be non-trimmed
    // If this is not the case, the constructor will throw an exception
    // solution_space_dimension specifies the dimension of each corresponding
    // solution volume
    IsogeometricVolModel(shared_ptr<VolumeModel> volmodel,
			 std::vector<int> solution_space_dimension);

    // Destructor
    virtual ~IsogeometricVolModel();

    // One boundary condition is added to the model
    // The area on the geometry model corresponding to this boundary
    // condition is specified as a polygon where each corner is given by
    // its geometric position and a pointer to the face on the volume 
    // boundary where this polygon corner lies.
    // The corners either belongs to the set of surface corners in the boundary
    // surfaces of the volume model, or they are points necessary to specify a
    // polygon that does not entirely follow the edges of the boundary surfaces
    // The functor (if given) can be evaluated and gives the value of the boundary
    // condition given the corresponding point on the geometry boundary curve
    // The functor is given primarily if the boundary condition is of type
    // Dirichlet and different from zero or constant Dirichlet. The argument to the
    // functor is of type Point and it returns a type Point. The dimension of the
    // argument is the same as the dimension of the geometry space, the output
    // dimension is that of the current solution space
    // In the case of constant Dirichlet, the constant value should be given
    // in the parameter constant_value
    // A boundary condition corresponds to one solution. The index of this
    // solution must be given.
    // The polygon must define a nonempty area.
    bool addBoundaryCond(std::vector<std::pair<ParamSurface*, Point> > polygon,
			 BdConditionType type, BdCondFunctor *fbd,
			 int solutionspace_idx, double *constant_value = 0);

    // A pointwise boundary condition of type Dirichlet is given. The position of the
    // boundary condition at the geometric model is given along with the boundary surface
    // where the condition is positioned. The value of the condition is given in the parameter
    // condition_value
    void addDirichletPointBdCond(ParamSurface* bd_surf, Point& pos,
				 Point& condition_value,
				 int solutionspace_idx);

    // Get the number of boundaries of the given volume model
    virtual int getNmbOfBoundaries() const;

    // Returns the outer boundary of the volume model
    // A surface model is a surface set including adjacency information
    std::vector<shared_ptr<ParamSurface> > getOuterBoundary() const;

    // Returns a specified boundary of the volume model
    std::vector<shared_ptr<ParamSurface> > getBoundary(int idx) const;

    // Ensure minimum degree of solution space
    virtual void setMinimumDegree(int degree, int solutionspace_idx);

    // Update spline spaces of the solution to ensure consistence
    virtual void updateSolutionSplineSpace();

    // Update spline spaces of the solution to ensure consistence
    virtual void updateSolutionSplineSpace(int solutionspace_idx);

    // Fetch all the single block defining this multi-block model
    void getIsogeometricBlocks(std::vector<shared_ptr<IsogeometricVolBlock> >& volblock);

  private:
    // The blocks which this block structured model consist of
    std::vector<shared_ptr<IsogeometricVolBlock> > vol_blocks_;

    /// The boundary surface models. The outer boundary surface model
    /// is the first element.
    std::vector<std::vector<shared_ptr<ParamSurface> > > boundary_surfaces_;

    // The block index of each boundary face segment
    std::vector<std::vector<int> > boundary_surface_block_;

    // The position of each boundary face
    // The value is a number from 0-11, given as p+o where
    // - p is 0, 2, 4, 6, 8 or 10 for pos (umin, umax, vmin, vmax, wmin, wmax)
    // - r is 0 for same orientation as on twin surface, 1 for reversed orientation
    std::vector<std::vector<int> > boundary_surface_pos_;

    // Update the geometry entity in the volume blocks to have corresponding
    // coefficients between neighbours
    void makeGeometrySplineSpaceConsistent();

    // Rebuild the array of boundary curves
    void buildBoundaryFaces(shared_ptr<VolumeModel> volmodel);

    // Get number of solution spaces
    int nmbSolutionSpaces() const;

  };   // end class IsogeometricVolModel

} // end namespace Go


#endif    // #ifndef __ISOGEOMETRICVOLMODEL_H
