//===========================================================================
//                                                                           
// File: QualityResults
//                                                                           
// Created: April 2008
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _QUALITYRESULTS_H
#define _QUALITYRESULTS_H

#include "GoTools/qualitymodule/testSuite.h"
#include "GoTools/compositemodel/ftEdge.h"
#include "GoTools/compositemodel/Vertex.h"
#include <vector>

namespace Go
{
  class ftSurface;
  class ftCurve;
  class ftPoint;
  class PointOnCurve;
  class PointOnEdge;

  class QualityResults
  {
    friend class ModelQuality;
    friend class FaceSetQuality;
    friend class ModelRepair;
    friend class FaceSetRepair;

  public:

    // Destructor
    ~QualityResults();

  private:
    bool test_performed_[TEST_SUITE_SIZE];  
    double tolerance_used_[TEST_SUITE_SIZE];

    // Constructor
    QualityResults();

    void reset(testSuite whichtest);
    void performtest(testSuite whichtest, double tol);
    bool testPerformed(testSuite whichtest, double& tol);
	    
    // Result of test for degenerate surface boundaries
    std::vector<shared_ptr<ftSurface> > deg_sfs_;
      void addDegSf(shared_ptr<ftSurface> degface)
    {
      deg_sfs_.push_back(degface);
    }

    std::vector<shared_ptr<ftSurface> > getDegSfs()
	{
	    return deg_sfs_;
	}

    // Test on paralell or anti paralell derivative in surface corner
    std::vector<shared_ptr<ftPoint> > deg_sf_corners_;
    void addDegenerateSfCorner(shared_ptr<ftPoint> deg_corner)
    {
      deg_sf_corners_.push_back(deg_corner);
    }

    std::vector<shared_ptr<ftPoint> > getDegCorners()
	{
	    return deg_sf_corners_;
	}
    
    // Test on identical vertices
    std::vector<std::pair<shared_ptr<Vertex>,
      shared_ptr<Vertex> > > identical_vertices_;
    void addIdenticalVertices(std::pair<shared_ptr<Vertex>,
			   shared_ptr<Vertex> > vertex_pair)
    {
      identical_vertices_.push_back(vertex_pair);
    }

    std::vector<std::pair<shared_ptr<Vertex>, shared_ptr<Vertex> > > 
	getIdenticalVertices()
	{
	    return identical_vertices_;
	}

    // Test on identical edges
    std::vector<std::pair<shared_ptr<ftEdge>,
      shared_ptr<ftEdge> > > identical_edges_;
    void addIdenticalEdges(std::pair<shared_ptr<ftEdge>,
			   shared_ptr<ftEdge> > edge_pair)
    {
      identical_edges_.push_back(edge_pair);
    }

     std::vector<std::pair<shared_ptr<ftEdge>, shared_ptr<ftEdge> > > 
	 getIdenticalEdges()
	 {
	     return identical_edges_;
	 }

     // Test on embedded edges
    // The second edge in the pair is embedded in the first one
    std::vector<std::pair<shared_ptr<ftEdge>,
      shared_ptr<ftEdge> > > embedded_edges_;
    void addEmbeddedEdges(std::pair<shared_ptr<ftEdge>,
			  shared_ptr<ftEdge> > edge_pair)
    {
      embedded_edges_.push_back(edge_pair);
    }

     std::vector<std::pair<shared_ptr<ftEdge>, shared_ptr<ftEdge> > > 
	 getEmbeddedEdges()
	 {
	     return embedded_edges_;
	 }

     // Test on identical faces
    std::vector<std::pair<shared_ptr<ftSurface>,
      shared_ptr<ftSurface> > > identical_faces_;
    void addIdenticalFaces(std::pair<shared_ptr<ftSurface>,
			   shared_ptr<ftSurface> > face_pair)
    {
      identical_faces_.push_back(face_pair);
    }

    std::vector<std::pair<shared_ptr<ftSurface>,shared_ptr<ftSurface> > > 
	getIdenticalFaces()
	{
	    return identical_faces_;    
	}

// Test on embedded faces
    // The second face in the pair is embedded in the first one
    std::vector<std::pair<shared_ptr<ftSurface>,
      shared_ptr<ftSurface> > > embedded_faces_;
    void addEmbeddedFaces(std::pair<shared_ptr<ftSurface>,
			  shared_ptr<ftSurface> > face_pair)
    {
      embedded_faces_.push_back(face_pair);
    }

    std::vector<std::pair<shared_ptr<ftSurface>,shared_ptr<ftSurface> > > 
	getEmbeddedFaces()
	{
	    return embedded_faces_;    
	}

    // Test on mini surface
    std::vector<shared_ptr<ParamSurface> > mini_surface_;
    void addMiniSurface(shared_ptr<ParamSurface> surface)
	{
	    mini_surface_.push_back(surface);
	}

   std::vector<shared_ptr<ParamSurface> > 
       getMiniSurfaces()
       {
	   return mini_surface_;
       }
 
    // Test on mini edges
    std::vector<shared_ptr<ftEdge> > mini_edges_;
    void addMiniEdge(shared_ptr<ftEdge> mini_edge)
    {
	mini_edges_.push_back(mini_edge);
    }

    std::vector<shared_ptr<ftEdge> > 
	getMiniEdges()
	{
	    return mini_edges_;
	}

    // Test on mini face
    std::vector<shared_ptr<ftSurface> > mini_face_;
    void addMiniFace(shared_ptr<ftSurface> face)
	{
	    mini_face_.push_back(face);
	}

    std::vector<shared_ptr<ftSurface> > 
	getMiniFaces()
	{
	    return mini_face_;
	}

    // Test on vanishing surface normal
    std::vector<shared_ptr<ftPoint> > singular_points_;
    std::vector<shared_ptr<ftCurve> > singular_curves_;
    void addSingPnt(shared_ptr<ftPoint> sing_pnt)
    {
      singular_points_.push_back(sing_pnt);
    }
    void addSingCurve(shared_ptr<ftCurve> sing_curve)
    {
      singular_curves_.push_back(sing_curve);
    }

    std::vector<shared_ptr<ftPoint> > getSingPnts()
	{
	    return singular_points_;
	}

    std::vector<shared_ptr<ftCurve> > getSingCrvs()
	{
	    return singular_curves_;
	}

    // Test on vanishing curve tangent
    std::vector<shared_ptr<PointOnCurve> > sing_points_crv_;
    std::vector<std::pair<shared_ptr<PointOnCurve>, 
	shared_ptr<PointOnCurve> > > sing_curves_crv_;
    void addSingPntCrv(shared_ptr<PointOnCurve> sing_pnt)
    {
      sing_points_crv_.push_back(sing_pnt);
    }
    void addSingCurveCrv(std::pair<shared_ptr<PointOnCurve>, 
			 shared_ptr<PointOnCurve> > sing_curve)
    {
      sing_curves_crv_.push_back(sing_curve);
    }
    

    std::vector<shared_ptr<PointOnCurve> > getSingCrvPnts()
	{
	    return sing_points_crv_;
	}

    std::vector<std::pair<shared_ptr<PointOnCurve>, shared_ptr<PointOnCurve> > >
	getSingCrvCrvs()
	{
	    return sing_curves_crv_;
	}

    // Test on sliver faces
    std::vector<shared_ptr<ftSurface> > sliver_sfs_;
      void addSliverSf(shared_ptr<ftSurface> sliver_face)
    {
      sliver_sfs_.push_back(sliver_face);
    }

    std::vector<shared_ptr<ftSurface> > getSliverSfs()
	{
	    return sliver_sfs_;
	}

    //  Test on narrow region in face
    std::vector<std::pair<shared_ptr<PointOnEdge>,
	shared_ptr<PointOnEdge> > > narrow_region_;
    void addNarrowRegion(std::pair<shared_ptr<PointOnEdge>, 
			 shared_ptr<PointOnEdge> > narrow_region)
    {
	narrow_region_.push_back(narrow_region);
    }

    std::vector<std::pair<shared_ptr<PointOnEdge>, shared_ptr<PointOnEdge> > > 
	getNarrowRegion()
	{
	    return narrow_region_;
	}

    // Test on edge-vertex distance
    std::vector<std::pair<ftEdge*,
      shared_ptr<Vertex> > > edge_vertices_;
    void addDistantEdgeVertex(std::pair<ftEdge*, shared_ptr<Vertex> > edge_vertex)
    {
      edge_vertices_.push_back(edge_vertex);
    }

    std::vector<std::pair<ftEdge*,shared_ptr<Vertex> > > getDistantEdgeVertex()
	{
	    return edge_vertices_;
	}

    // Test on face-vertex distance
    std::vector<std::pair<ftSurface*,
      shared_ptr<Vertex> > > face_vertices_;
    void addDistantFaceVertex(std::pair<ftSurface*, shared_ptr<Vertex> > face_vertex)
    {
      face_vertices_.push_back(face_vertex);
    }

    std::vector<std::pair<ftSurface*,shared_ptr<Vertex> > > getDistantFaceVertex()
	{
	    return face_vertices_;
	}

    // Test on face-edge distance
    std::vector<std::pair<ftSurface*, ftEdge*> > face_edges_;
    void addDistantFaceEdge(std::pair<ftSurface*, ftEdge*> face_edge)
    {
      face_edges_.push_back(face_edge);
    }

    std::vector<std::pair<ftSurface*,ftEdge*> > getDistantFaceEdge()
	{
	    return face_edges_;
	}

    // Test on position discontinuity between edges
    std::vector<std::pair<ftEdge*, ftEdge*> > pos_discont_edges_;
    void addEdgePositionDiscont(std::pair<ftEdge*, ftEdge*>& pos_discont_edge)
    {
      pos_discont_edges_.push_back(pos_discont_edge);
    }
    std::vector<std::pair<ftEdge*, ftEdge*> > getEdgePosDiscont()
	{
	    return pos_discont_edges_;
	}

    // Test on tangent discontinuity between edges
    std::vector<std::pair<ftEdge*, ftEdge*> > tangent_discont_edges_;
    void addEdgeTangentDiscont(std::pair<ftEdge*, ftEdge*>& tangent_discont_edge)
    {
      tangent_discont_edges_.push_back(tangent_discont_edge);
    }
    std::vector<std::pair<ftEdge*, ftEdge*> > getEdgeTangDiscont()
	{
	    return tangent_discont_edges_;
	}

    // Test on position discontinuity between faces
    std::vector<std::pair<ftEdge*, ftEdge*> > pos_discont_faces_;
    void addFacePositionDiscont(std::pair<ftEdge*, ftEdge*>& pos_discont_face)
    {
      pos_discont_faces_.push_back(pos_discont_face);
    }
    std::vector<std::pair<ftEdge*, ftEdge*> > getFacePosDiscont()
	{
	    return pos_discont_faces_;
	}


    // Test on tangent discontinuity between faces
    std::vector<std::pair<ftEdge*, ftEdge*> > tangent_discont_faces_;
    void addFaceTangentDiscont(std::pair<ftEdge*, ftEdge*>& tangent_discont_face)
    {
      tangent_discont_faces_.push_back(tangent_discont_face);
    }
    std::vector<std::pair<ftEdge*, ftEdge*> > getFaceTangDiscont()
	{
	    return tangent_discont_faces_;
	}

    // Test on consistency of orientation within a loop
    std::vector<shared_ptr<ftEdge> > edge_in_loop_;
    void addInconsistentEdgeInLoop(shared_ptr<ftEdge> inconsistent_edge)
	{
	    edge_in_loop_.push_back(inconsistent_edge);
	}
    std::vector<shared_ptr<ftEdge> > getInconsistLoopEdge()
	{
	    return edge_in_loop_;
	}

    // Test on orientation of a boundary loop
    std::vector<shared_ptr<Loop> > loop_orientation_;
    void addInconsistentBdLoop(shared_ptr<Loop> inconsistent_loop)
	{
	    loop_orientation_.push_back(inconsistent_loop);
	}
    std::vector<shared_ptr<Loop> > getInconsistLoop()
	{
	    return loop_orientation_;
	}

    // Test on orientation of faces
    std::vector<shared_ptr<ftSurface> > face_orientation_;
    void addInconsistentFaceInSet(shared_ptr<ftSurface> inconsistent_face)
	{
	    face_orientation_.push_back(inconsistent_face);
	}
    std::vector<shared_ptr<ftSurface> > getInconsistFace()
	{
	    return face_orientation_;
	}

    // Test in C1- and G1-discontinuity internally in a curve
    std::vector<shared_ptr<ParamCurve> > c1_discont_cvs_;
    std::vector<shared_ptr<ParamCurve> > g1_discont_cvs_;
    void addC1DiscontCv(shared_ptr<ParamCurve> discont)
	{
	    c1_discont_cvs_.push_back(discont);
	}
    void addG1DiscontCv(shared_ptr<ParamCurve> discont)
	{
	    g1_discont_cvs_.push_back(discont);
	}
    std::vector<shared_ptr<ParamCurve> > getC1DiscontCvs()
	{
	    return c1_discont_cvs_;
	}
    std::vector<shared_ptr<ParamCurve> > getG1DiscontCvs()
	{
	    return g1_discont_cvs_;
	}

    // Test on G1-discontinuity internally in a face
    std::vector<shared_ptr<ftSurface> > g1_discont_sfs_;
    void addG1DiscontSf(shared_ptr<ftSurface> discont)
	{
	    g1_discont_sfs_.push_back(discont);
	}

    // Test on C1-discontinuity internally in a face
    std::vector<shared_ptr<ftSurface> > c1_discont_sfs_;
    void addC1DiscontSf(shared_ptr<ftSurface> discont)
	{
	    c1_discont_sfs_.push_back(discont);
	}
    std::vector<shared_ptr<ftSurface> > getG1DiscontSfs()
	{
	    return g1_discont_sfs_;
	}
    std::vector<shared_ptr<ftSurface> > getC1DiscontSfs()
	{
	    return c1_discont_sfs_;
	}

    // Minimum curvature radius, curve
    std::vector<std::pair<shared_ptr<PointOnCurve>, double> >  cv_curvature_;
    std::pair<shared_ptr<PointOnCurve>, double> minimum_cv_curvature_radius_;
    void smallCvCurvRad(std::pair<shared_ptr<PointOnCurve>, double> curv_rad)
    {
	cv_curvature_.push_back(curv_rad);
    }
    std::vector<std::pair<shared_ptr<PointOnCurve>, double> >  getSmallCvCurvatureR()
	{
	    return cv_curvature_;
	}
    void setMinimumCvCurvatureRadius(std::pair<shared_ptr<PointOnCurve>, double> curv_rad)
	{
	    minimum_cv_curvature_radius_ = curv_rad;
	}
    std::pair<shared_ptr<PointOnCurve>, double> getMinCvCurvatureR()
	{
	    return minimum_cv_curvature_radius_;
	}

    // Minimum curvature radius, surface
    std::vector<std::pair<shared_ptr<ftPoint>, double> >  sf_curvature_;
    std::pair<shared_ptr<ftPoint>, double> minimum_sf_curvature_radius_;
    void smallSfCurvRad(std::pair<shared_ptr<ftPoint>, double> curv_rad)
    {
	sf_curvature_.push_back(curv_rad);
    }
    std::vector<std::pair<shared_ptr<ftPoint>, double> >  getSmallSfCurvatureR()
	{
	    return sf_curvature_;
	}
    void setMinimumCurvatureRadius(std::pair<shared_ptr<ftPoint>, double> curv_rad)
	{
	    minimum_sf_curvature_radius_ = curv_rad;
	}
    std::pair<shared_ptr<ftPoint>, double> getMinSfCurvatureR()
	{
	    return minimum_sf_curvature_radius_;
	}

    // Acute angle between edges
    std::vector<std::pair<ftEdge*, ftEdge*> > edge_acute_angle_;
    void addEdgeAcuteAngle(std::pair<ftEdge*, ftEdge*> acute_angle)
	{
	    edge_acute_angle_.push_back(acute_angle);
	}
    std::vector<std::pair<ftEdge*, ftEdge*> > getEdgeAcuteAngle()
	{
	    return edge_acute_angle_;
	}

    // Acute angle between faces
    std::vector<std::pair<ftSurface*, ftSurface*> > face_acute_angle_;
    void addFaceAcuteAngle(std::pair<ftSurface*, ftSurface*> acute_angle)
	{
	    face_acute_angle_.push_back(acute_angle);
	}
    std::vector<std::pair<ftSurface*, ftSurface*> > getFaceAcuteAngle()
	{
	    return face_acute_angle_;
	}

    // Intersection between boundary loops
    std::vector<std::pair<shared_ptr<PointOnEdge>, 
	shared_ptr<PointOnEdge> > > loop_intersection_;
    void addLoopIntersection(std::pair<shared_ptr<PointOnEdge>, 
			     shared_ptr<PointOnEdge> > loop_int)
    {
	loop_intersection_.push_back(loop_int);
    }
    std::vector<std::pair<shared_ptr<PointOnEdge>,shared_ptr<PointOnEdge> > > 
	getIntersectingBdLoops()
	{
	    return loop_intersection_;
	}
	
    // Self intersection of boundary loop
    std::vector<std::pair<shared_ptr<PointOnEdge>, 
	shared_ptr<PointOnEdge> > > loop_self_intersection_;
    void addLoopSelfIntersection(std::pair<shared_ptr<PointOnEdge>, 
			     shared_ptr<PointOnEdge> > loop_int)
    {
	loop_self_intersection_.push_back(loop_int);
    }
    std::vector<std::pair<shared_ptr<PointOnEdge>, shared_ptr<PointOnEdge> > > 
	getSelfIntersectingBdLoops()
	{
	    return loop_self_intersection_;
	}
	
    // Curve with indistinct knots
    std::vector<shared_ptr<ParamCurve> > cv_indistinct_knots_;
    void addCvIndistinctKnot(shared_ptr<ParamCurve> crv)
	{
	    cv_indistinct_knots_.push_back(crv);
	}
    std::vector<shared_ptr<ParamCurve> > getCvIndistinctKnots()
	{
	    return cv_indistinct_knots_;
	}

    // Surface with indistinct knots
    std::vector<shared_ptr<ParamSurface> > sf_indistinct_knots_;
    void addSfIndistinctKnot(shared_ptr<ParamSurface> srf)
	{
	    sf_indistinct_knots_.push_back(srf);
	}
    std::vector<shared_ptr<ParamSurface> > getSfIndistinctKnots()
	{
	    return sf_indistinct_knots_;
	}

  };

} // namespace Go

#endif // _QUALITYRESULTS_H
