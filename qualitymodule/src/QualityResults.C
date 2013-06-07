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



