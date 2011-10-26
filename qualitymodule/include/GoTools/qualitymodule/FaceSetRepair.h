//===========================================================================
//                                                                           
// File: FaceSetRepair
//                                                                           
// Created: November 2009
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================
#ifndef _FACESETREPAIR_H
#define _FACESETREPAIR_H

#include "GoTools/qualitymodule/ModelRepair.h"
#include "GoTools/qualitymodule/QualityResults.h"
#include "GoTools/qualitymodule/testSuite.h"
#include <vector>

namespace Go
{
  class SurfaceModel;
  class FaceSetQuality;

  class FaceSetRepair : public ModelRepair
  {
  public:
    /// Constructor
    FaceSetRepair(std::shared_ptr<SurfaceModel> sfmodel);

    /// Constructor, given quality check class
    FaceSetRepair(std::shared_ptr<FaceSetQuality> quality);

    /// Destructor
    ~FaceSetRepair();


    /// Remove gaps in the model
    virtual void mendGaps();

    /// Remove identical and embedded faces
    virtual void identicalAndEmbeddedFaces();

    /// Handle identical vertices
    virtual void identicalVertices();

    /// Handle inconsistent face normals
    virtual void consistentFaceNormal();

    virtual void optimizeVertexPosition();

    virtual void mendEdgeDistance();
    
      std::shared_ptr<SurfaceModel> getAssociatedSfModel()
      {
	  return sfmodel_;
      }


  private:
    std::shared_ptr<SurfaceModel>  sfmodel_;
    std::shared_ptr<FaceSetQuality> quality_;

    bool vertex_update_;
    bool edges_update_;

    void gapTrimming(std::vector<std::pair<ftEdge*, ftEdge*> >& pos_discont,
		     double epsge, bool update_iso);

    void improveVertexPos(std::shared_ptr<Vertex> vx, double epsge);

    void getBoundaryEdges(ftEdge *e1, std::vector<ftEdge*>& edges, double epsge);

    void getBoundaryEdges(ftEdge *e1, ftEdge *e2, 
			  std::vector<std::pair<ftEdge*,ftEdge*> >& edges, 
			  double epsge);

    void updateOneSpline(std::vector<std::pair<ftEdge*, ftEdge*> >& pos_discont,
			 double epsge);
    void 
      updateMeetingSplines(std::vector<std::pair<ftEdge*, ftEdge*> >& pos_discont,
			   double epsge); 

    void 
      averageMeetingSplines(std::vector<std::pair<ftEdge*, ftEdge*> >& pos_discont,
			   double epsge); 

    void 
      modifyOneSpline(std::vector<std::pair<ftEdge*, ftEdge*> >& pos_discont,
		      double epsge);

    void 
      modifyTwoSplines(std::vector<std::pair<ftEdge*, ftEdge*> >& pos_discont,
		       double epsge);

    void 
      getEdges(ftEdge *edge1, std::vector<ftEdge*>& edges, double epsge);

    void 
      getEdges(ftEdge *edge1, ftEdge *edge2,
	       std::vector<std::pair<ftEdge*, ftEdge*> >& edges, double epsge);

};
} // namespace Go

#endif // _FACESETREPAIR_H
