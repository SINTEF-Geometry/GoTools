//===========================================================================
//                                                                           
// File: testSuite.h
//                                                                           
// Created: April 2008
//                                                                           
// Author: Vibeke Skytt 
//                                                                           
// Revision: 
//                                                                           
// Description: 
//
// Implementation files: 
//                                                                           
//===========================================================================

#ifndef _TESTSUITE_H
#define _TESTSUITE_H

//===========================================================================
//===========================================================================


namespace Go
{
    class SurfaceModel;

    enum testSuite
	{
	    IDENTICAL_VERTICES = 0,
	    IDENTICAL_EDGES = 1,
	    EMBEDDED_EDGES = 2,
	    IDENTICAL_FACES = 3,
	    EMBEDDED_FACES = 4,
	    MINI_CURVE = 5,
	    MINI_SURFACE = 6,
	    MINI_EDGE = 7,
	    MINI_FACE = 8,
	    SLIVER_FACE = 9,
	    NARROW_REGION = 10,
	    DEGEN_SRF_BD = 11,
	    DEGEN_SRF_CORNER = 12,
	    VANISHING_TANGENT = 13,
	    VANISHING_NORMAL = 14,
	    EDGE_VERTEX_DISTANCE = 15,
	    FACE_VERTEX_DISTANCE = 16,
	    FACE_EDGE_DISTANCE = 17,
	    EDGE_POSITION_DISCONT = 18,
	    EDGE_TANGENTIAL_DISCONT = 19,
	    FACE_POSITION_DISCONT = 20,
	    FACE_TANGENTIAL_DISCONT = 21,
	    LOOP_CONSISTENCY = 22,
	    LOOP_ORIENTATION = 23,
	    FACE_ORIENTATION = 24,
	    CV_G1DISCONT = 25,
	    CV_C1DISCONT = 26,
	    SF_G1DISCONT = 27,
	    SF_C1DISCONT = 28,
	    CV_CURVATURE_RADIUS = 29,
	    SF_CURVATURE_RADIUS = 30,
	    EDGE_ACUTE_ANGLE = 31,
	    FACE_ACUTE_ANGLE = 32,
	    LOOP_INTERSECTION = 33,
	    LOOP_SELF_INTERSECTION = 34,
	    INDISTINCT_KNOTS = 35
	};
#define TEST_SUITE_SIZE 36

} // namespace Go

#endif //  _TESTSUITE_H

       
