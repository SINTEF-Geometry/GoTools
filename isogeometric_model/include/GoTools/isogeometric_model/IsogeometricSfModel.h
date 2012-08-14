//===========================================================================
//
// File : IsogeometricSfModel.h
//
// Created: 27th of October, 2009
//
// Author: Vibeke Skytt, SINTEF, and Anh-Vu Vuong, TU Munich
//
// Revision: 
//
// Description:
//
//===========================================================================



#ifndef __ISOGEOMETRICSFMODEL_H
#define __ISOGEOMETRICSFMODEL_H



#include <vector>
#include <memory>
#include "GoTools/utils/Point.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/isogeometric_model/IsogeometricModel.h"
#include "GoTools/isogeometric_model/BdConditionType.h"
#include "GoTools/isogeometric_model/IsogeometricSfBlock.h"


namespace Go
{

  // This class is intended for storage and service functionality
  // for one isogeometric model of type structured multi-block surface model.
  // It contains geometry, boundary conditions, and the corresponding
  // solutions space(s)

  class IsogeometricSfModel : public IsogeometricModel
  {
  public:
    // Constructor
    // The geometry is described as a surface model, i.e. a surface set
    // and adjacency information. All surfaces should be non-trimmed
    // If this is not the case, the constructor will throw an exception
    // solution_space_dimension specifies the dimension of each corresponding
    // solution surface
    IsogeometricSfModel(shared_ptr<SurfaceModel> sfmodel,
			std::vector<int> solution_space_dimension);

    // Destructor
    virtual ~IsogeometricSfModel();

    // The boundary pointer specifies whether the current boundary condition
    // is related to the outer boundary of the geometry model or some inner
    // boundary
    // The pos vector specifies the sequence of the boundary corresponding
    // to the current boundary conditions. All positions at joints between 
    // segments in the boundary belonging to the current boundary
    // condition must be included in the sequence. The sequence must be ordered.
    // The points are geometric points, i.e. not parametric points.
    // The segments are [pos_0, pos_1), [pos_1, pos_2), ...
    // The first and the last point in the sequence belong to the 
    // segment joints only if endpoint of the boundary condition corresponds
    // to one such joint
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
    bool addBoundaryCond(int boundary, std::vector<Point> pos,
			 BdConditionType type, BdCondFunctor *fbd,
			 int solutionspace_idx, double *constant_value = 0);

    // A pointwise boundary condition of type Dirichlet is given. The position of the
    // boundary condition at the geometric model is given along with the boundary 
    // where the condition is positioned. The value of the condition is given in 
    // the parameter condition_value
    void addDirichletPointBdCond(int boundary, Point& pos,
				 Point& condition_value,
				 int solutionspace_idx);

    // Get the number of boundaries of the given surface model
    virtual int getNmbOfBoundaries() const;

    // Return the outer boundary of the surface model
    CurveLoop getOuterBoundary() const;

    // Return the specified boundary curve
    // The outer boundary will always have index 0, inner holes in the surface
    // model have indices from 1
    CurveLoop getBoundary(int idx) const;

    // Ensure minimum degree of solution space
    // The solution space will always have at least the degree of the
    // corresponding geometry surface
    virtual void setMinimumDegree(int degree, int solutionspace_idx);

    // Update spline spaces of the solution to ensure consistency
    virtual void updateSolutionSplineSpace();

    // Update spline spaces of the solution to ensure consistency
    virtual void updateSolutionSplineSpace(int solutionspace_idx);

    // Fetch all the single block defining this multi-block model
    void getIsogeometricBlocks(std::vector<shared_ptr<IsogeometricSfBlock> >& sfblock);

  private:
    // The blocks which this block structured model consist of
    std::vector<shared_ptr<IsogeometricSfBlock> > sf_blocks_;

    /// The boundary curves. The outer curve loop is the first element.
    std::vector<CurveLoop> boundary_curves_;

    // The block index of each boundary curve segment
    std::vector<std::vector<int> > boundary_curve_block_;

    // The edge position of each boundary curve
    // The value is a number from 0-7, given as e+r where
    // - e is 0, 2, 4 or 6 for edge given by resp umin, umax, vmin, vmax
    // - r is 0 for same orientation as on surfaces, 1 for reversed orientation
    std::vector<std::vector<int> > boundary_curve_edge_;

    // Update the geometry entity in the surface blocks to have corresponding
    // coefficients between neighbours
    void makeGeometrySplineSpaceConsistent();

    // Rebuild the array of boundary curves
    void buildBoundaryCurves(shared_ptr<SurfaceModel> sfmodel);

    // Get number of solution spaces
    int nmbSolutionSpaces() const;

  };   // class IsogeometricSfModel

} // end namespace Go


#endif    // #ifndef __ISOGEOMETRICSFMODEL_H
