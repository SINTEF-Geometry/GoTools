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

#ifndef _GVUTILITIES_H
#define _GVUTILITIES_H

/// Draw a cylinder, or possibly a cone.
void draw_cylinder(double x0, double y0, double z0,
                   double x1, double y1, double z1,
                   double radius, double radius2, int n);
/// Draw a set of axes.
void draw_gl_axes(int n, double r, double radius, double rim, double l);
void draw_gl_axes(double relscale);

typedef double FLOAT_TYPE;
typedef FLOAT_TYPE MATRIX3[9];
typedef FLOAT_TYPE MATRIX4[16];
FLOAT_TYPE m3_det( MATRIX3 mat );
void m3_identity( MATRIX3 mat );
void m3_inverse( MATRIX3 mr, MATRIX3 ma );
void m4_submat( MATRIX4 mr, MATRIX3 mb, int i, int j );
FLOAT_TYPE m4_det( MATRIX4 mr );
int m4_inverse( MATRIX4 mr, MATRIX4 ma );

#endif // _GVUTILITIES_H


