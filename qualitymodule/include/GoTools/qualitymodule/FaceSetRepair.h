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
    FaceSetRepair(shared_ptr<SurfaceModel> sfmodel);

    /// Constructor, given quality check class
    FaceSetRepair(shared_ptr<FaceSetQuality> quality);

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
    
      shared_ptr<SurfaceModel> getAssociatedSfModel()
      {
	  return sfmodel_;
      }


  private:
    shared_ptr<SurfaceModel>  sfmodel_;
    shared_ptr<FaceSetQuality> quality_;

    bool vertex_update_;
    bool edges_update_;

    void gapTrimming(std::vector<std::pair<ftEdge*, ftEdge*> >& pos_discont,
		     double epsge, bool update_iso);

    void improveVertexPos(shared_ptr<Vertex> vx, double epsge);

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
