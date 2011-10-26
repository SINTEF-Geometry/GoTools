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

#include "GoTools/qualitymodule/QualityResults.h"


namespace Go
{

  //===========================================================================
  QualityResults::QualityResults()
  //===========================================================================  
  {
      for (int ki = 0; ki < TEST_SUITE_SIZE; ++ki)
      {
	  test_performed_[ki] = false;
	  tolerance_used_[ki] = -1.0;
      }
  }  


  //===========================================================================
  QualityResults::~QualityResults()
  //===========================================================================  
  {
  }

  //===========================================================================
  void QualityResults::reset(testSuite whichtest)
  //===========================================================================  
  {
    test_performed_[(int)whichtest] = false;
    tolerance_used_[(int)whichtest] = -1.0;

    switch (whichtest)
      {
      case IDENTICAL_VERTICES:
	  identical_vertices_.clear();
	break;
      case IDENTICAL_EDGES:
	  identical_edges_.clear();
	break;
      case EMBEDDED_EDGES:
	  embedded_edges_.clear();
	break;
      case IDENTICAL_FACES:
	identical_faces_.clear();
	break;
      case EMBEDDED_FACES:
	embedded_faces_.clear();
	break;
      case MINI_CURVE:
	break;
      case MINI_SURFACE:
	  mini_surface_.clear();
	break;
      case MINI_EDGE:
	  mini_edges_.clear();
	break;
      case MINI_FACE:
	  mini_face_.clear();
	break;
      case SLIVER_FACE:
	sliver_sfs_.clear();
	break;
      case NARROW_REGION:
	  narrow_region_.clear();
	break;
      case DEGEN_SRF_BD:
	deg_sfs_.clear();
	break;
      case DEGEN_SRF_CORNER:
	  deg_sf_corners_.clear();
	break;
      case VANISHING_TANGENT:
	  sing_points_crv_.clear();
	  sing_curves_crv_.clear();
	break;
      case VANISHING_NORMAL:
	{
	  singular_points_.clear();
	  singular_curves_.clear();
	  break;
	}
      case EDGE_VERTEX_DISTANCE:
	{
	  edge_vertices_.clear();
	  break;
	}
      case FACE_VERTEX_DISTANCE:
	{
	  face_vertices_.clear();
	  break;
	}
      case FACE_EDGE_DISTANCE:
	{
	  face_edges_.clear();
	  break;
	}
      case EDGE_POSITION_DISCONT:
	{
	  pos_discont_edges_.clear();
	  break;
	}
      case EDGE_TANGENTIAL_DISCONT:
	{
	  tangent_discont_edges_.clear();
	  break;
	}
      case FACE_POSITION_DISCONT:
	{
	  pos_discont_faces_.clear();
	  break;
	}
      case FACE_TANGENTIAL_DISCONT:
	{
	  tangent_discont_faces_.clear();
	  break;
	}
	  case LOOP_CONSISTENCY:
	  {
	      edge_in_loop_.clear();
	      break;
	  }
	  case LOOP_ORIENTATION:
	  {
	      loop_orientation_.clear();
	      break;
	  }
	  case FACE_ORIENTATION:
	  {
	      face_orientation_.clear();
	      break;
	  }
	  case CV_G1DISCONT:
	      {
		  g1_discont_cvs_.clear();
		  break;
	      }
	  case CV_C1DISCONT:
	      {
		  c1_discont_cvs_.clear();
		  break;
	      }
	  case SF_G1DISCONT:
	      {
		  g1_discont_sfs_.clear();
		  break;
	      }
	  case SF_C1DISCONT:
	      {
		  c1_discont_sfs_.clear();
		  break;
	      }
      case CV_CURVATURE_RADIUS:
	{
	    cv_curvature_.clear();
	  break;
	}
      case SF_CURVATURE_RADIUS:
	{
	  sf_curvature_.clear();
	  break;
	}
	  case EDGE_ACUTE_ANGLE:
	  {
	      edge_acute_angle_.clear();
	      break;
	  }
	  case FACE_ACUTE_ANGLE:
	  {
	      face_acute_angle_.clear();
	      break;
	  }
	  case LOOP_INTERSECTION:
	  {
	      loop_intersection_.clear();
	      break;
	  }
	  case LOOP_SELF_INTERSECTION:
	  {
	      loop_self_intersection_.clear();
	      break;
	  }
	  case INDISTINCT_KNOTS:
	  {
	      cv_indistinct_knots_.clear();
	      sf_indistinct_knots_.clear();
	      break;
	  }

      }
  }

  //===========================================================================
  void QualityResults::performtest(testSuite whichtest, double tol)
  //===========================================================================  
  {
      test_performed_[(int)whichtest] = true;
      tolerance_used_[(int)whichtest] = tol;
  }


  //===========================================================================
  bool QualityResults::testPerformed(testSuite whichtest, double& tol)
  //===========================================================================  
  {
      tol = tolerance_used_[(int)whichtest];
      return test_performed_[(int)whichtest];
  }


} // namespace Go



