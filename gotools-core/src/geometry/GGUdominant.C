//==========================================================================
//                                                                          
// File: GGUdominant.C                                                       
//                                                                          
// Created: Wed Mar  5 14:36:55 2003                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: GGUdominant.C,v 1.4 2005-02-22 13:16:22 afr Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include "GoTools/geometry/GeometryTools.h"


using namespace std;


namespace Go {


//==========================================================================
void findDominant(const SplineSurface& surface,
		  Vector3D& dominant_u, Vector3D& dominant_v)
//==========================================================================
{
    int nu = surface.numCoefs_u();
    int nv = surface.numCoefs_v();
    vector<double>::const_iterator start = surface.coefs_begin();
    Vector3D temp;
    // Dominant in u-direction
    dominant_u = Vector3D(0.0, 0.0, 0.0);
    for (int j = 0; j < nv; ++j) {
	for (int dd = 0; dd < 3; ++dd) {
	    temp[dd] = *(start + 3*(nu*j + (nu-1)) + dd)
		- *(start + 3*(nu*j) + dd);
	}
	dominant_u += temp;
    }
    // Dominant in v-direction
    dominant_v = Vector3D(0.0, 0.0, 0.0);
    for (int i = 0; i < nu; ++i) {
	for (int dd = 0; dd < 3; ++dd) {
	    temp[dd] = *(start + 3*(nu*(nv-1) + i) + dd)
		- *(start + 3*i + dd);
	}
	dominant_v += temp;
    }

    return;
}


//==========================================================================
bool negativeProj(const SplineSurface& surface,
		  const Array<Vector3D, 2>& refvector,
		  const double eps)
//==========================================================================
{
    int num_u = surface.numCoefs_u();
    int num_v = surface.numCoefs_v();
    Vector3D temp;
    int i = 0, j = 0;
    while (i < num_u-1) {
	j = 0;
	while (j < num_v) {
	    temp[0] = *(surface.coefs_begin() + 3*(num_u*j + i+1))
		- *(surface.coefs_begin() + 3*(num_u*j + i));
	    temp[1] = *(surface.coefs_begin() + 3*(num_u*j + i+1) + 1)
		- *(surface.coefs_begin() + 3*(num_u*j + i) + 1);
	    temp[2] = *(surface.coefs_begin() + 3*(num_u*j + i+1) + 2)
		- *(surface.coefs_begin() + 3*(num_u*j + i) + 2);
	    // Positive tolerance means that there must be a small
	    // _nonzero_ negative projection before it is reported as
	    // negative!
	    if (temp * refvector[0] < -eps)
		return true;
	    ++j;
	}
	++i;
    }
    i = 0;
    while (i < num_u) {
	j = 0;
	while (j < num_v-1) {
	    temp[0] = *(surface.coefs_begin() + 3*(num_u*(j+1) + i))
		- *(surface.coefs_begin() + 3*(num_u*j + i));
	    temp[1] = *(surface.coefs_begin() + 3*(num_u*(j+1) + i) + 1)
		- *(surface.coefs_begin() + 3*(num_u*j + i) +1);
	    temp[2] = *(surface.coefs_begin() + 3*(num_u*(j+1) + i) + 2)
		- *(surface.coefs_begin() + 3*(num_u*j + i) + 2);
	    // Positive tolerance means that there must be a small
	    // _nonzero_ negative projection before it is reported as
	    // negative!
	    if (temp * refvector[1] < -eps)
		return true;
	    ++j;
	}
	++i;
    }

    return false;
}



//==========================================================================


} // namespace Go
