//===========================================================================
//
// File : IsogeometricVolModel.h
//
// Created: 28th of October, 2009
//
// Author: Vibeke Skytt, SINTEF, and Anh-Vu Vuong, TU Munich
//
// Revision: 
//
// Description:
//
//===========================================================================



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
