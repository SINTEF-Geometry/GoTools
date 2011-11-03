//===========================================================================
//
// File : IsogeometricVolBlock.h
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



#ifndef __ISOGEOMETRICVOLBLOCK_H
#define __ISOGEOMETRICVOLBLOCK_H


#include <vector>
#include <memory>
#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/utils/Point.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/isogeometric_model/IsogeometricBlock.h"
#include "GoTools/isogeometric_model/VolSolution.h"
#include "GoTools/isogeometric_model/VolBoundaryCondition.h"
#include "GoTools/isogeometric_model/VolPointBdCond.h"



namespace Go
{

  // This class represents one block in a block-structured isogeometric volume model
  class IsogeometricVolBlock : public IsogeometricBlock
  {
  public:
    // Functions used to create the block
    // Constructor
    IsogeometricVolBlock(IsogeometricModel* model,
			 std::shared_ptr<SplineVolume> geom_vol,
			 std::vector<int> solution_space_dimension,
			 int index);

    // Destructor
    virtual ~IsogeometricVolBlock();

    // Cast as IsogeometricVolBlock
    virtual IsogeometricVolBlock* asIsogeometricVolBlock();

    // Multiblock. Add neighbourhood information
    // This function is called from IsogeometricVolModel and used in building the
    // complete block structured model
    // same_dir_order: True if corr surfaces also have corr u-dir.
    void addNeighbour(std::shared_ptr<IsogeometricVolBlock> neighbour,
		      int face_nmb_this,
		      int face_nmb_other,
		      int orientation,
		      bool same_dir_order);

    // Count the number of neighbouring volume blocks to this block
    virtual int nmbOfNeighbours() const;

    // Return the neighbour along a specified boundary. If this bounary corresponds to an 
    // outer boundary, a zero point is returned
    // bd_nmb: 0=umin, 1=umax, 2=vmin, 3=vmax, 4=wmin, 5=wmax
    IsogeometricVolBlock* getNeighbour(int bd_nmb) const;

    // Functions used to access data
    // Given this block and another one, check if they are neighbours
    virtual bool isNeighbour(IsogeometricBlock* other) const;

    // Get the total number of coefficients in the block
    virtual int nmbCoefs() const;

    // Get B-spline basis in one paramenter direction in the block
    // The first parameter direction has pardir = 0, etc
    virtual BsplineBasis basis(int pardir) const;

    // Return the specified boundary surface, 
    // face_number: 0=umin, 1=umax, 2=vmin, 3=vmax, 4=wmin, 5=wmax
    std::shared_ptr<SplineSurface> getGeomBoundarySurface(int face_number) const;

    // Given a point lying on a specified boundary, return the corresponding parameter
    // value of the corresponding boundary surface
    // face_number: 0=umin, 1=umax, 2=vmin, 3=vmax, 4=wmin, 5=wmax
    std::vector<double> getParamOnBdSurf(int face_number, const Point& position) const;

    // Fetch information about degenerate volume boundaries with respect 
    // to a given tolerance epsge
    // The boundary face numbers of any degenerate edges are collected in
    // the first entity in the pair, the second entity specifies the type
    // of degeneracy
    // 1 - one surface boundary degenerate to a point,
    // 2 - surface degenerate to line, 3 - surface degenerate to point    
    bool geomIsDegenerate(std::vector<std::pair<int,int> >& degen_bd, 
			  double epsge);

    // Get the local enumeration of all coefficients belonging to degenerate
    // boundary faces. The vector indexes of the degen_bd and 
    // enumeration correspond. degen_bd is defined as before
    void getDegenEnumeration(std::vector<std::pair<int,int> >& degen_bd, 
			     std::vector<std::vector<int> >& enumeration,
			     double epsge);

    // Check if the volume block is periodic with respect to the tolerance
    // epsge
    // per[pardir] = -1 : Not closed or periodic
    //                0 : Closed or periodic with continuity of position
    //                1 : Periodic with continuity of position and derivative
    //                >1 : Higher order continuity
    bool geomIsPeriodic(int per[], double epsge);

    // Fetch enumeration along a periodic boundary in a specified parameter
    // direction. If the volume is not closed or periodic, the function 
    // returns false
    // Only the enumeration of the boundary coefficients are returned even
    // if the continuity is higher than C0 across the seam.
    bool getPeriodicEnumeration(int pardir,
				std::vector<std::pair<int, int> >& enumeration);

    // Refine the geometry volume
    // The solution spline space is refine accordingly such that the 
    // geometry space is always a sub space of the solution space
    // Insert a number of specified new knots in a specified parameter direction
    virtual void refineGeometry(std::vector<double> newknots, int pardir);

    // Refine the geometry model in a specified direction in such a way that the
    // spline space specified as input will be a sub space of the spline space of
    // geometry object. The minimum refinement to achieve this is chosen.
    virtual void refineGeometry(const BsplineBasis& other_basis, int pardir);

    // Increase degree of the geometry volume in a given parameter direction. 
    // If new_degree is not larger than the current degree, no change is performed
    virtual void increaseGeometryDegree(int new_degree, int pardir);

    // Update the current geometry volume with respect to a given boundary curve
    // If the spline space of the volume is not able to fit the new boundary
    // exactly, approximation is performed
    void updateGeometry(std::shared_ptr<SplineSurface> new_boundary, int face_number);

    // Release scratch related to pre evaluated basis functions and surface.
    virtual void erasePreEvaluatedBasisFunctions();

    // Fetch boundary conditions
    // Get the number of boundary conditions attached to this block
    virtual int getNmbOfBoundaryConditions() const;

    // Get a specified boundary condition
    std::shared_ptr<VolBoundaryCondition> getBoundaryCondition(int index) const;

    // Get all boundary conditions related to a specified face. Boundary conditions
    // related to all solution spaces are returned
    void getFaceBoundaryConditions(int face_number, 
				   std::vector<std::shared_ptr<VolBoundaryCondition> >& bd_cond) const;

    // Get number of point type bounday conditions
    virtual int getNmbOfPointBdConditions() const;
  
    // Get a specified point type boundary condition
    std::shared_ptr<VolPointBdCond> getPointBdCondition(int index) const;

    // Get all point boundary conditions related to a specified face. Boundary 
    // conditions related to all solution spaces are returned
    void getFacePointBdConditions(int face_number, 
				  std::vector<std::shared_ptr<VolPointBdCond> >& bd_cond) const;

    // Get specified solution space
    std::shared_ptr<VolSolution> getSolutionSpace(int solution_index);

    // Get geometry surface
    std::shared_ptr<SplineVolume> volume() const;

    // Ensure minimum degree of solution space
    // The solution space will always have at least the degree of the
    // corresponding geometry volume
    virtual void setMinimumDegree(int degree, int solutionspace_idx);

    // Update spline spaces of the solution to ensure consistency
    // Returns true if any update occured, false if not
    // Solution space index is a global value valid for all blocks in a model.
    virtual bool updateSolutionSplineSpace(int solutionspace_idx);

    // Get number of solution spaces
    int nmbSolutionSpaces() const;

    // The face position of a boundary surface
    // Returns -1 if not possible to determine face position.
    //   0, 4 & 8 for parameter direction u, v & w (respectively)
    // + 2 if end param (umax, vmax or wmax)
    // + 1 if orientation is reversed with respect to orientation on volume
    int getFaceOrientation(std::shared_ptr<ParamSurface> srf, double tol);

  private:

    // The volume describing the geometry
    std::shared_ptr<SplineVolume> volume_;

    // The position index in the Model object
    int index_;

    // Solution spaces
    std::vector<std::shared_ptr<VolSolution> > solution_;

    // Adjacency
    std::shared_ptr<IsogeometricVolBlock> neighbours_[6];

    // Adjacency position at other block:
    // neighb_face_[i] = 0 : neighbours_[i] has this block as neighbour along face u = u_min
    // neighb_face_[i] = 1 : neighbours_[i] has this block as neighbour along face u = u_max
    // neighb_face_[i] = 2 : neighbours_[i] has this block as neighbour along face v = v_min
    // neighb_face_[i] = 3 : neighbours_[i] has this block as neighbour along face v = v_max
    // neighb_face_[i] = 4 : neighbours_[i] has this block as neighbour along face w = w_min
    // neighb_face_[i] = 5 : neighbours_[i] has this block as neighbour along face w = w_max
    int neighb_face_[6];

    // Information about how the (neighbour) volumes are oriented
    // @@sbr201111 I guess is "turned" in the u-direction should be "reversed"?
    // orientation_[i] = 0 : The volumes are oriented in the same way
    // orientation_[i] = 1 : The neighbouring volume is turned in the u-direction
    //                       compared with the current volume
    // @@sbr201111 The u-direction, is that the direction of neighbouring volume, or the
    //             direction corresponding to u-direction of current volume?
    // orientation_[i] = 2 : The neighbouring volume is turned in the v-direction
    //                       compared with the current volume
    // orientation_[i] = 3 : The neighbouring volume is turned in the w-direction
    //                       compared with the current volume
    // orientation_[i] = 4 : The neighbouring volume is turned in two first parameter 
    //                       directions compared with the current volume
    // orientation_[i] = 5 : The neighbouring volume is turned in first and last parameter 
    //                       directions compared with the current volume
    // orientation_[i] = 6 : The neighbouring volume is turned in two last parameter 
    //                       directions compared with the current volume
    // orientation_[i] = 7 : The neighbouring volume is turned in all three parameter 
    //                       directions compared with the current volume
    int orientation_[6];

    /// True if u-directions concide for both boundary surfaces (regardless of orientation)
    bool same_dir_order_[6];

  };   // end class IsogeometricVolBlock

} // end namespace Go


#endif    // #ifndef __ISOGEOMETRICVOLBLOCK_H
