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

	    /* /// Constructor with coordinate system. NB! The coordinate system is not in use */
	    /* Body(const CoordinateSystem<3> xyz); */

	    /// Default constructor 
	    Body();

	    /* /// Constructor with multiple shells */
	    /* Body(const CoordinateSystem<3> xyz, std::vector<shared_ptr<SurfaceModel> >& shells); */

	    /// Constructor with multiple shells
	    Body(std::vector<shared_ptr<SurfaceModel> >& shells, 
		 int material_id=-1);

	    /* /// Constructor with one outer shell */
	    /* Body(const CoordinateSystem<3> xyz, shared_ptr<SurfaceModel> shell); */

	    /// Constructor with one outer shell
	    Body(shared_ptr<SurfaceModel> shell, int material_id=-1);

	    /// Destructor
	    ~Body();

	    /// Add another shell to the solid description. Used in
	    /// connection with building the solid
	    void addshell(shared_ptr<SurfaceModel> shell);

	    /// Fetch the outer shell belonging to this solid
	    shared_ptr<SurfaceModel> getOuterShell() const
		{
		  shared_ptr<SurfaceModel> dummy;
		  if (shells_.size() == 0)
		    return dummy;
		  else
		    return shells_[0];
		}

	    /// Fetch all shells belonging to this solid. The first
	    /// shell corresponds to the outer boundary of the solid,
	    /// any subsequent shells are void
	    std::vector<shared_ptr<SurfaceModel> > getAllShells() const
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
	    shared_ptr<SurfaceModel> getShell(int idx) const
	    {
	       shared_ptr<SurfaceModel> dummy;
	      if (idx < 0 || int(shells_.size()) <= idx)
		return dummy;
	      else
		return shells_[idx];
	    }

	    /// Material information
	    void setMaterial(int material_id)
	    {
	      material_id_ = material_id;
	    }

	    bool hasMaterialInfo() const
	    {
	      return (material_id_ >= 0);
	    }
	    
	    int getMaterial() const
	    {
	      return material_id_;
	    }

	    /// Total number of faces
	    int nmbOfFaces() const;

	    /// Fetch all vertices
	    std::vector<shared_ptr<Vertex> > vertices() const;

	    /// The bounding box corresponding to this solid
	    virtual BoundingBox boundingBox() const;

	    /// Check if two bodies are neighbours
	    virtual bool areNeighbours(Body *other, shared_ptr<ftSurface>& bd_face1,
				       shared_ptr<ftSurface>& bd_face2, int adj_idx = 0) const;

	    /// Fetch all bodies adjacent to this one. This function is
	    /// relevant in an assembly system or in an isogeometric setting
	    void getAdjacentBodies(std::vector<Body*>& neighbours);

	    /// Remove adjacency information to other bodies
	    void eraseBodyAdjacency();

	    /// Check if a point lies inside this body
	    bool isInside(const Point& pnt) const;

	    /// Find the shell containing a given face (if any)
	    shared_ptr<SurfaceModel> getShell(ftSurface* face) const;

	    /// Return topology tolerances
	    tpTolerances getTolerances()
	    {
	      return toptol_;
	    }

	protected:
	    /* /// Coordinate system of the body. Default is the xyz system. */
	    /* /// The coordinate system is currently not used */
	    /* CoordinateSystem<3> coordinate_; */
	    
	    /// Outer and possibly inner void shells
	    std::vector<shared_ptr<SurfaceModel> > shells_;

	    /// Possibility to store material information (-1 = not set)
	    int material_id_;

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
