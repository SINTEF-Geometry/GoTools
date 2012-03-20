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
    CompleteEdgeNet(shared_ptr<SurfaceModel> sfmodel);

    /// Destructor
    ~CompleteEdgeNet();

    /// Apply algorithm for completing the edge net and create
    /// a starting ground for a block structured model of solids
    /// \return Whether the edge net was completed or not. 
    bool perform();

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

    /// Given a regular solid, add the edges required to make a block 
    /// structured model
    void addMissingEdges();

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

    bool vertexInfo(shared_ptr<Vertex> vx, double& angle);

    void writePath(std::vector<ftEdge*>& edges, shared_ptr<Vertex> vx);

    std::vector<ftEdge*> getStartEdges();
  };

} // namespace Go

#endif // COMPLETEEDGENET_H
