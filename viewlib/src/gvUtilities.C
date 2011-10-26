//===========================================================================
//                                                                           
// File: gvUtilities.C                                                       
//                                                                           
// Created: Mon Jun 25 15:12:11 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: gvUtilities.C,v 1.1 2007-04-17 12:25:57 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/viewlib/gvUtilities.h"
#include "GoTools/utils/Values.h"
#ifdef _MSC_VER
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif
#include <math.h>



//----------------------------------------------------------------------
//
// Drawing a cylinder, or possibly a cone.
//
//----------------------------------------------------------------------

void draw_cylinder(double x0, double y0, double z0,
                   double x1, double y1, double z1,
                   double radius, double radius2, int n)
{
  int i;
  double y, z, r;
  
  glPushMatrix();
  glTranslatef((GLfloat)x0, (GLfloat)y0, (GLfloat)z0);
  /*
    Now, we have to move P=(x1-x0, y1-y0, z1-z0) to (r, 0, 0) where
    r=length(x1-x0, y1-y0, z1-z0). First, rotate P to Q in the xy-plane:
  */
  /* Hvorfor blir det ikke riktig med neg. rotasjon her? */
  glRotatef(GLfloat(atan2(x1-x0, z1-z0)/M_PI*180.0-90.0),
            0.0, 1.0, 0.0);
/*  printf("%f %f %f\t",
         x1-x0,
         z1-z0, 
         atan2(x1-x0, z1-z0)/M_PI*180.0);*/
  /*
    Now, we have to move Q=(sqrt((x1-x0)^2+(z1-z0)^2), y1-y0, 0) to
    (r, 0, 0).
  */
  r=sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0));
  /* Hvorfor blir det ikke riktig med neg. rotasjon her? */
  glRotatef(GLfloat(atan2(y1-y0,
        sqrt((x1-x0)*(x1-x0)+(z1-z0)*(z1-z0)))/M_PI*180.0),
            0.0, 0.0, 1.0);
/*  printf("%f %f %f\n",
         y1-y0,
         sqrt((x1-x0)*(x1-x0)+(z1-z0)*(z1-z0)),
         atan2(y1-y0,
sqrt((x1-x0)*(x1-x0)+(z1-z0)*(z1-z0)))/M_PI*180.0);*/
  /*
    The next step is to draw the cylinder. Or cone. Normals not perfect
for
    the cone. Impossible at pointy end with tri-fan?
  */
  if (radius==radius2)
    {
      glBegin(GL_QUAD_STRIP);
      for (i=0; i<=n; i++)
        {
          y=radius*cos((2.0*M_PI*i)/n);
          z=radius*sin((2.0*M_PI*i)/n);
          glNormal3f(0.0, GLfloat(y/radius), GLfloat(z/radius));
          glVertex3f(0.0, GLfloat(y), GLfloat(z));
          glNormal3f(0.0, GLfloat(y/radius), GLfloat(z/radius));
          glVertex3f(  GLfloat(r), GLfloat(y), GLfloat(z));
        }
      glEnd();
    }
  else
    {
      glBegin(GL_TRIANGLE_FAN);
      glNormal3f(1.0, 0.0, 0.0);
      glVertex3f(float(r), 0.0, 0.0);
      for (i=0; i<=n; i++)
        {
          y=radius*cos((2.0*M_PI*i)/n);
          z=radius*sin((2.0*M_PI*i)/n);
          glNormal3f(0.0, float(-y/radius), float(-z/radius));
          glVertex3f(0.0, float(y), float(z));
        }
      glEnd();
    }
  /*
    Finally, we restore the transformation stack.
  */
  glPopMatrix();
}




//----------------------------------------------------------------------
//
// Drawing a set of axes.
//
//----------------------------------------------------------------------

void draw_gl_axes(int n, double r, double radius, double rim, double l)
{
  // int n=10;
  // double r=0.7;
  // double radius=0.01, rim=0.04, l=0.1;
  
    
    draw_cylinder(0.0, 0.0, -r, 0.0, 0.0, r, radius, radius, n);
    draw_cylinder(0.0, 0.0, r-0.1, 0.0, 0.0, r+0.3, rim, 0.0, n);
    draw_cylinder(-0.5*l, l, r+0.3, 0.5*l, l, r+0.3,
		  radius, radius, n);
    draw_cylinder(-0.5*l, l*2.0, r+0.3, 0.5*l, l*2.0, r+0.3,
		  radius, radius, n);
    draw_cylinder(-0.5*l, l, r+0.3, 0.5*l, l*2.0, r+0.3,
		  radius, radius, n);
  
    draw_cylinder(-r, 0.0, 0.0, r, 0.0, 0.0, radius, radius, n);
    draw_cylinder(r-0.1, 0.0, 0.0, r+0.3, 0.0, 0.0, rim, 0.0, n);
    draw_cylinder(r+0.3, l, 0.5*l, r+0.3, 2*l, -0.5*l,
		  radius, radius, n);
    draw_cylinder(r+0.3, l, -0.5*l, r+0.3, 2*l, 0.5*l,
		  radius, radius, n);
  
    draw_cylinder(0.0, -r, 0.0, 0.0, r, 0.0, radius, radius, n);
    draw_cylinder(0.0, r-0.1, 0.0, 0.0, r+0.3, 0.0, rim, 0.0, n);
    draw_cylinder(0.0, r+0.3, 2.0*l, 0.0, r+0.3, 1.5*l,
		  radius, radius, n);
    draw_cylinder(0.0, r+0.3, 1.5*l, -0.5*l, r+0.3, l,
		  radius, radius, n);
    draw_cylinder(0.0, r+0.3, 1.5*l, 0.5*l, r+0.3, l,
		  radius, radius, n);
  
}

void draw_gl_axes(double relscale)
{
    draw_gl_axes(10,
		 relscale*1.0,
		 relscale*0.015,
		 relscale*0.02,
		 relscale*0.07);
}


// From http://skal.planet-d.net/demo/matrixfaq.htm

FLOAT_TYPE m3_det( MATRIX3 mat )
{
    FLOAT_TYPE det;

    det = mat[0] * ( mat[4]*mat[8] - mat[7]*mat[5] )
	- mat[1] * ( mat[3]*mat[8] - mat[6]*mat[5] )
	+ mat[2] * ( mat[3]*mat[7] - mat[6]*mat[4] );

    return( det );
}


void m3_identity( MATRIX3 mat )
{
    mat[0] = mat[4] = mat[8] = 1.0;
    mat[1] = mat[2] = mat[3] = mat[5] = mat[6] = mat[7] = 0.0;
}

void m3_inverse( MATRIX3 mr, MATRIX3 ma )
{
    FLOAT_TYPE det = m3_det( ma );

    if ( fabs( det ) < 0.0005 )
	{
	    m3_identity( ma );
	    return;
	}

    mr[0] =    ma[4]*ma[8] - ma[5]*ma[7]   / det;
    mr[1] = -( ma[1]*ma[8] - ma[7]*ma[2] ) / det;
    mr[2] =    ma[1]*ma[5] - ma[4]*ma[2]   / det;

    mr[3] = -( ma[3]*ma[8] - ma[5]*ma[6] ) / det;
    mr[4] =    ma[0]*ma[8] - ma[6]*ma[2]   / det;
    mr[5] = -( ma[0]*ma[5] - ma[3]*ma[2] ) / det;

    mr[6] =    ma[3]*ma[7] - ma[6]*ma[4]   / det;
    mr[7] = -( ma[0]*ma[7] - ma[6]*ma[1] ) / det;
    mr[8] =    ma[0]*ma[4] - ma[1]*ma[3]   / det;
}

void m4_submat( MATRIX4 mr, MATRIX3 mb, int i, int j )
{
    // idst and jdst is arbitrarily initialized to avoid complaints
    // from the compiler. But strictly speaking we should check if i
    // and j are acceptable... @@@jbt.

    for (int ti = 0; ti < 4; ti++ ) {
	int idst = i;
	if ( ti < i )
	    idst = ti;
	else
	    if ( ti > i )
		idst = ti-1;
	
	for (int tj = 0; tj < 4; tj++ ) {
	    int jdst = j;
	    if ( tj < j )
		jdst = tj;
	    else
		if ( tj > j )
		    jdst = tj-1;

	    if ( ti != i && tj != j )
		mb[idst*3 + jdst] = mr[ti*4 + tj ];
	}
    }
}

FLOAT_TYPE m4_det( MATRIX4 mr )
{
    FLOAT_TYPE  det, result = 0, i = 1;
    MATRIX3 msub3;
    int     n;

    for ( n = 0; n < 4; n++, i *= -1 )
        {
	    m4_submat( mr, msub3, 0, n );

	    det     = m3_det( msub3 );
	    result += mr[n] * det * i;
        }

    return( result );
}

int m4_inverse( MATRIX4 mr, MATRIX4 ma )
{
    FLOAT_TYPE  mdet = m4_det( ma );
    MATRIX3 mtemp;
    int     i, j, sign;

    if ( fabs( mdet ) < 0.0005 )
        return( 0 );

    for ( i = 0; i < 4; i++ )
        for ( j = 0; j < 4; j++ )
	    {
		sign = 1 - ( (i +j) % 2 ) * 2;

		m4_submat( ma, mtemp, i, j );

		mr[i+j*4] = ( m3_det( mtemp ) * sign ) / mdet;
	    }

    return( 1 );
}

