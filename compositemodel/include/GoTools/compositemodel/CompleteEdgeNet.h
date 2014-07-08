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

#ifndef COMPLETEEDGENET_H
#define COMPLETEEDGENET_H

#include "GoTools/compositemodel/SurfaceModel.h"

namespace Go
{
  /// \brief Complete the edge net of a SurfaceModel. Fetches the wire
  /// frame model corresponding to a surface model, and extends it such
  /// that the extended wire frame will become the wire frame model corresponding
  /// to a volume model have the given surface model as its outer boundary.
  /// The extension to the wireframe is represented as pairs of vertices where
  /// these vertices lie at the endpoints of the missing edges.
  /// NB! This solution is currently not expected to handle all configurations.

  class CompleteEdgeNet
  {
  public:
    /// Constructor. The method is applied on a SurfaceModel
    /// \param sfmodel Pointer to a SurfaceModel
    CompleteEdgeNet(shared_ptr<SurfaceModel> sfmodel, 
		    bool perform_step2, bool smooth_connections);

    /// Destructor
    ~CompleteEdgeNet();

    /// Apply algorithm for completing the edge net and create
    /// a starting ground for a block structured model of solids
    /// \return Whether the edge net was completed or not. 
    bool perform(std::vector<std::pair<Point,Point> >& corr_vx_pts);

    /// Fetch modified model (after regularization)
    /// \return Pointer to the modified model
    shared_ptr<SurfaceModel> getRegularizedModel()
      {
	return model_;
      }

    /// Fetch new edges represented by their end vertices
    /// \return  Vector of pointers to end vertices of the edges
    std::vector<std::pair<shared_ptr<Vertex>, shared_ptr<Vertex> > >
      getMissingEdges()
      {
	return missing_edges_;
      }

  private:
    shared_ptr<SurfaceModel> model_;
    std::vector<std::pair<shared_ptr<Vertex>, shared_ptr<Vertex> > > missing_edges_;
    bool perform_step2_;
    bool smooth_connections_;

    /// Given a regular solid, add the edges required to make a block 
    /// structured model
    void addMissingEdges();

    void identifyVertexConnection(std::vector<shared_ptr<Vertex> > vxs,
				  size_t ki1, size_t ki2, 
				  int& ix1, int& ix2);

    void traverseEdges(std::vector<ftEdge*>& edges,
		       std::vector<ftEdge*>& curr_path,
		       std::vector<int>& curr_idx,
		       ftEdge *curr_edge,
		       shared_ptr<Vertex> vx,
		       bool search_end); 
 
    ftEdge* fetchNextEdge(ftEdge *curr_edge,
			  shared_ptr<Vertex> vx,
			  int& next_idx);

    bool regularizeEdgeLoop(std::vector<ftEdge*>& edges);

    void splitLoop(std::vector<ftEdge*>& edges,
		   std::vector<shared_ptr<Vertex> >& vxs,
		   bool to_add_edges,
		   std::vector<std::vector<ftEdge*> >& split_loops,
		   std::vector<std::vector<shared_ptr<Vertex> > >& split_vxs,
		   std::vector<bool>& add_edges_split);

    bool regularizeCurrLoop(std::vector<ftEdge*>& edges,
			    std::vector<shared_ptr<Vertex> >& vxs,
			    bool to_add_edges);

    double getVertexAngle(ftEdge *edge1, ftEdge *edge2);

    void addRemainingEdges();

    bool vertexInfo(shared_ptr<Vertex> vx, double& angle, Point& centre);

    void writePath(std::vector<ftEdge*>& edges, shared_ptr<Vertex> vx);

    std::vector<ftEdge*> getStartEdges();

    void addIdentifiedEdges(std::vector<std::pair<Point,Point> >& corr_vx_pts);

    bool betterConnectionInFace(shared_ptr<Vertex> source, Body *bd,
				ftEdge* edg1, ftEdge *edg2, shared_ptr<Vertex> dest);
  };

} // namespace Go

#endif // COMPLETEEDGENET_H
