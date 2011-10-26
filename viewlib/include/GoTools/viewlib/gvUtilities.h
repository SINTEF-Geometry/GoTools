//===========================================================================
//                                                                           
// File: gvUtilities.h                                                       
//                                                                           
// Created: Mon Jun 25 15:11:32 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: gvUtilities.h,v 1.1 2007-04-17 12:25:45 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

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


