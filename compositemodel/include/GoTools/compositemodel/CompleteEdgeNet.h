#ifndef COMPLETEEDGENET_H
#define COMPLETEEDGENET_H

#include "GoTools/compositemodel/SurfaceModel.h"

namespace Go
{
    /// Complete the edge net of a SurfaceModel.
  class CompleteEdgeNet
  {
  public:
    /// Constructor. The method is applied on a SurfaceModel
    /// \param sfmodel Pointer to a SurfaceModel
    CompleteEdgeNet(std::shared_ptr<SurfaceModel> sfmodel);

    /// Destructor
    ~CompleteEdgeNet();

    /// Apply algorithm for completing the edge net and create
    /// a starting ground for a block structured model of solids
    /// \return Whether the edge net was completed or not. 
    bool perform();

    /// Fetch modified model (after regularization)
    /// \return Pointer to the modified model
    std::shared_ptr<SurfaceModel> getRegularizedModel()
      {
	return model_;
      }

    /// Fetch new edges represented by their end vertices
    /// \return  Vector of pointers to end vertices of the edges
    std::vector<std::pair<std::shared_ptr<Vertex>, std::shared_ptr<Vertex> > >
      getMissingEdges()
      {
	return missing_edges_;
      }

  private:
    std::shared_ptr<SurfaceModel> model_;
    std::vector<std::pair<std::shared_ptr<Vertex>, std::shared_ptr<Vertex> > > missing_edges_;

    /// Given a regular solid, add the edges required to make a block 
    /// structured model
    void addMissingEdges();

    void traverseEdges(std::vector<ftEdge*>& edges,
		       std::vector<ftEdge*>& curr_path,
		       ftEdge *curr_edge,
		       std::shared_ptr<Vertex> vx,
		       bool search_end); 
 
    ftEdge* fetchNextEdge(ftEdge *curr_edge,
			  std::shared_ptr<Vertex> vx,
			  std::vector<ftEdge*> edges);

    void regularizeEdgeLoop(std::vector<ftEdge*> edges);

    double getVertexAngle(ftEdge *edge1, ftEdge *edge2);

    void writePath(std::vector<ftEdge*> edges, std::shared_ptr<Vertex> vx);
  };

} // namespace Go

#endif // COMPLETEEDGENET_H
