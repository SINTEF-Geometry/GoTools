//===========================================================================
//                                                                           
// File: Body.h                                                    
//                                                                           
// Created: September 2008
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description: Boundary representated solid model
//              The model will be face based.
//                                                                           
//===========================================================================

// VSK. September 2008
// Note that body requires a set of test functionality
// Check that each shell is closed
// Check that the shells don't intersect
// Check that the inner void shells is inside the outer shell
// These tests are not implemented yet
// If the contained shells have different topology tolerances, the body gets the
// largest tolerances. Should topological results as neighbourhood, gap, kink etc. 
// be updated with respect to the largest tolerances?
// The functionality is currently very lean, and will be extended when needed

#ifndef _BODY_H
#define _BODY_H

#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/utils/CoordinateSystem.h"
#include "GoTools/utils/BoundingBox.h"
#include <vector>

namespace Go
{

  /// \brief A boundary represented solid.

    class Body
	{
	public:

	    /// Constructor with coordinate system
	    Body(const CoordinateSystem<3> xyz);

	    /// Default constructor 
	    Body();

	    /// Constructor with multiple shells
	    Body(const CoordinateSystem<3> xyz, std::vector<std::shared_ptr<SurfaceModel> >& shells);

	    /// Constructor with multiple shells
	    Body(std::vector<std::shared_ptr<SurfaceModel> >& shells);

	    /// Constructor with one outer shell
	    Body(const CoordinateSystem<3> xyz, std::shared_ptr<SurfaceModel> shell);

	    /// Constructor with one outer shell
	    Body(std::shared_ptr<SurfaceModel> shell);

	    /// Destructor
	    ~Body();

	    /// Add another shell to the solid description. Used in
	    /// connection with building the solid
	    void addshell(std::shared_ptr<SurfaceModel> shell);

	    /// Fetch the outer shell belonging to this solid
	    std::shared_ptr<SurfaceModel> getOuterShell() const
		{
		  std::shared_ptr<SurfaceModel> dummy;
		  if (shells_.size() == 0)
		    return dummy;
		  else
		    return shells_[0];
		}

	    /// Fetch all shells belonging to this solid. The first
	    /// shell corresponds to the outer boundary of the solid,
	    /// any subsequent shells are void
	    std::vector<std::shared_ptr<SurfaceModel> > getAllShells() const
	      {
		return shells_;
	      }

	    /// The number of shells belonging to this solid
	    int nmbOfShells() const
	    {
	      return (int)(shells_.size());
	    }

	    /// Fetch a specified shell. Index zero corresponds to
	    /// the outer shell
	    std::shared_ptr<SurfaceModel> getShell(int idx) const
	    {
	       std::shared_ptr<SurfaceModel> dummy;
	      if (idx < 0 || int(shells_.size()) <= idx)
		return dummy;
	      else
		return shells_[idx];
	    }

	    /// Fetch all vertices
	    std::vector<std::shared_ptr<Vertex> > vertices() const;

	    /// The bounding box corresponding to this solid
	    virtual BoundingBox boundingBox() const;

	    /// Check if two bodies are neighbours
	    virtual bool areNeighbours(Body *other, std::shared_ptr<ftSurface>& bd_face1,
				       std::shared_ptr<ftSurface>& bd_face2, int adj_idx = 0) const;

	    /// Fetch all bodies adjacent to this one. This function is
	    /// relevant in an assembly system or in an isogeometric setting
	    void getAdjacentBodies(std::vector<Body*>& neighbours);

	    /// Check if a point lies inside this body
	    bool isInside(const Point& pnt);

	    /// Find the shell containing a given face (if any)
	    std::shared_ptr<SurfaceModel> getShell(ftSurface* face) const;

	protected:
	    /// Coordinate system of the body. Default is the xyz system.
	    /// The coordinate system is currently not used
	    CoordinateSystem<3> coordinate_;
	    
	    /// Outer and possibly inner void shells
	    std::vector<std::shared_ptr<SurfaceModel> > shells_;

	    /// Tolerances used in topology analysis
	    // These tolerances needs to be stored with the class as a topology
	    // structure may become obsolete if the tolerances change
	    // The tolerances are fetched from the related CompositeModels.
	    tpTolerances toptol_;
   
	    // Add back pointers from boundary face entities
	    void addBodyPointers();
	};

} // namespace Go

#endif // _BODY_H
