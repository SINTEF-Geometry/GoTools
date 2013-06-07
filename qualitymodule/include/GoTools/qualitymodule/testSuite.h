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

       
