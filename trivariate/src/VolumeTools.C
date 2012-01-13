//===========================================================================
//                                                                           
// File: GeometryTools.C                                                     
//                                                                           
// Created: Fri Feb 18 12:53:43 2005                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: VolumeTools.C,v 1.2 2008-11-27 12:59:18 kfp Exp $
//                                                                           
//===========================================================================


#include "GoTools/trivariate/VolumeTools.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/trivariate/SurfaceOnVolume.h"
#include "GoTools/creators/HermiteAppC.h"
#include "GoTools/trivariate/VolumeParameterCurve.h"
#include "GoTools/trivariate/VolumeSpaceCurve.h"
#include <fstream>

//#define DEBUG

namespace Go
{

using std::vector;


//---------------------------------------------------------------------------
int analyzePeriodicity(const SplineVolume& sf, int direction, double knot_tol)
//---------------------------------------------------------------------------
{
    if (direction < 0 || direction > 1) {
	THROW("Error in direction parameter. Must be 0 or 1.");
    }
    // Make a hypercurve from this volume.
    // Direction parameter is one more for representVolume...() :-P
    shared_ptr<SplineCurve> cv
	= representVolumeAsCurve(sf, direction + 1);
    return analyzePeriodicity(*cv, knot_tol);
}


// Describe a volume as a curve in a given direction
//===========================================================================
shared_ptr<SplineCurve>
representVolumeAsCurve(const SplineVolume& volume, int cv_dir)
//===========================================================================
{
  int dim = volume.dimension();
  int kdim = dim + (volume.rational() ? 1 : 0);
  const std::vector<double>::const_iterator co = (volume.rational()) ?
    volume.rcoefs_begin() : volume.coefs_begin();
  std::vector<double> huge_curve_coefs;
  std::vector<double>::const_iterator coefstart; 
  int new_dim;
  if (cv_dir == 0)
    {
      int nu = volume.numCoefs(0);
      int nv = volume.numCoefs(1);
      int nw = volume.numCoefs(2);
      huge_curve_coefs.resize(nu*nv*nw*kdim);
      int pos = 0;
      for (int k = 0; k < nw; ++k)
	for (int j = 0; j < nv; ++j)
	  for (int i = 0; i < nu; ++i)
	    for (int d = 0; d < kdim; ++d)
	      huge_curve_coefs[((i*nw+k)*nv+j)*kdim+d] = co[pos++];
      coefstart = huge_curve_coefs.begin();
      new_dim = volume.numCoefs(1) * volume.numCoefs(2) * kdim;
    }
  else if (cv_dir == 1)
    {
      int nu = volume.numCoefs(0);
      int nv = volume.numCoefs(1);
      int nw = volume.numCoefs(2);
      huge_curve_coefs.resize(nu*nv*nw*kdim);
      int pos = 0;
      for (int k = 0; k < nw; ++k)
	for (int j = 0; j < nv; ++j)
	  for (int i = 0; i < nu; ++i)
	    for (int d = 0; d < kdim; ++d)
	      huge_curve_coefs[((j*nw+k)*nu+i)*kdim+d] = co[pos++];
      coefstart = huge_curve_coefs.begin();
      new_dim = volume.numCoefs(0) * volume.numCoefs(2) * kdim;
    }
  else
    {
      coefstart = co;
      new_dim = volume.numCoefs(0) * volume.numCoefs(1) * kdim;
    }

  const BsplineBasis& bas = volume.basis(cv_dir);
  int num = bas.numCoefs();
  std::vector<double>::const_iterator knotstart = bas.begin();

  shared_ptr<SplineCurve> curve(new SplineCurve(num, bas.order(),
						bas.begin(), coefstart,
						new_dim,
						false));
  return curve;
}



// Describe a curve as a volume in a given direction
//===========================================================================
shared_ptr<SplineVolume>
representCurveAsVolume(const SplineCurve& curve,
			int cv_dir,
			const BsplineBasis& other_bas1,
			const BsplineBasis& other_bas2,
			bool rational)
//===========================================================================
{
  if (curve.rational())
    THROW("It does not make sense to have a rational hypercurve.");

  int kdim = curve.dimension() / (other_bas1.numCoefs()*other_bas2.numCoefs());

  std::vector<double>::const_iterator co
      = curve.coefs_begin();
  std::vector<double> vol_coefs;
  std::vector<double>::const_iterator coefstart;
  const BsplineBasis& bas = curve.basis();

  if (cv_dir == 0)
    {
      int nu = bas.numCoefs();
      int nv = other_bas1.numCoefs();
      int nw = other_bas2.numCoefs();
      vol_coefs.resize(nu*nv*nw*kdim);
      int pos = 0;
      for (int k = 0; k < nw; ++k)
	for (int j = 0; j < nv; ++j)
	  for (int i = 0; i < nu; ++i)
	    for (int d = 0; d < kdim; ++d)
	      vol_coefs[pos++] = co[((i*nw+k)*nv+j)*kdim+d];
      coefstart = vol_coefs.begin();
    }
  else if (cv_dir == 1)
    {
      int nu = other_bas1.numCoefs();
      int nv = bas.numCoefs();
      int nw = other_bas2.numCoefs();
      vol_coefs.resize(nu*nv*nw*kdim);
      int pos = 0;
      for (int k = 0; k < nw; ++k)
	for (int j = 0; j < nv; ++j)
	  for (int i = 0; i < nu; ++i)
	    for (int d = 0; d < kdim; ++d)
	      vol_coefs[pos++] = co[((j*nw+k)*nu+i)*kdim+d];
      coefstart = vol_coefs.begin();
    }
  else
    coefstart = co;

  int dim = rational ? kdim - 1 : kdim;

  if (cv_dir == 0)
    { 
      shared_ptr<SplineVolume> volume
	(new SplineVolume(bas, other_bas1, other_bas2,
			  coefstart, dim, rational));
      return volume;
    }
  else if (cv_dir == 1)
    {
      shared_ptr<SplineVolume> volume
	(new SplineVolume(other_bas1, bas, other_bas2,
			  coefstart, dim, rational));
      return volume;
    }
  else
    {
      shared_ptr<SplineVolume> volume
	(new SplineVolume(other_bas1, other_bas2, bas,
			  coefstart, dim, rational));
      return volume;
    }

}


// Describe a volume as a surface in given directions
//===========================================================================
shared_ptr<SplineSurface>
representVolumeAsSurface(const SplineVolume& volume,
			 int sf_dir1,
			 int sf_dir2)
//===========================================================================
{
  int dim = volume.dimension();
  int kdim = dim + (volume.rational() ? 1 : 0);
  const std::vector<double>::const_iterator co = (volume.rational()) ?
    volume.rcoefs_begin() : volume.coefs_begin();
  std::vector<double> huge_surf_coefs;
  std::vector<double>::const_iterator coefstart; 

  int ijump, jjump, kjump;
  int nu = volume.numCoefs(0);
  int nv = volume.numCoefs(1);
  int nw = volume.numCoefs(2);
  int new_dim;

  if (sf_dir1 == sf_dir2)
    THROW("Must have two different parameter directions to describe a volume as a surface.");
  if (sf_dir1 < 0 || sf_dir1 > 2)
    THROW("sf_dir1 must be 0, 1 or 2");
  if (sf_dir2 < 0 || sf_dir2 > 2)
    THROW("sf_dir2 must be 0, 1 or 2");

  if (sf_dir1 == 0)
    {
      if (sf_dir2 == 1)
	{
	  kjump = 1;
	  ijump = nw;
	  jjump = nw*nu;
	  new_dim = nw * kdim;
	}
      else   // sf_dir2 == 2
	{
	  jjump = 1;
	  ijump = nv;
	  kjump = nv*nu;
	  new_dim = nv * kdim;
	}
    }
  else if (sf_dir1 == 1)
    {
      if (sf_dir2 == 0)
	{
	  kjump = 1;
	  jjump = nw;
	  ijump = nw*nv;
	  new_dim = nw * kdim;
	}
      else   // sf_dir2 == 2
	{
	  ijump = 1;
	  jjump = nu;
	  kjump = nu*nv;
	  new_dim = nu * kdim;
	}
    }
  else  // sf_dir1 == 2
    {
      if (sf_dir2 == 0)
	{
	  jjump = 1;
	  kjump = nv;
	  ijump = nv*nw;
	  new_dim = nv * kdim;
	}
      else   // sf_dir2 == 1
	{
	  ijump = 1;
	  kjump = nu;
	  jjump = nu*nw;
	  new_dim = nu * kdim;
	}
    }

  huge_surf_coefs.resize(nu*nv*nw*kdim);
  int pos = 0;
  for (int k = 0; k < nw; ++k)
    for (int j = 0; j < nv; ++j)
      for (int i = 0; i < nu; ++i)
	for (int d = 0; d < kdim; ++d)
	  huge_surf_coefs[(i*ijump + j*jjump + k*kjump)*kdim+d] = co[pos++];
  coefstart = huge_surf_coefs.begin();

  const BsplineBasis& bas1 = volume.basis(sf_dir1);
  const BsplineBasis& bas2 = volume.basis(sf_dir2);

  shared_ptr<SplineSurface> surface(new SplineSurface(bas1.numCoefs(), bas2.numCoefs(),
						      bas1.order(), bas2.order(),
						      bas1.begin(), bas2.begin(),
						      coefstart,
						      new_dim,
						      false));
  return surface;
}


// Describe a surface as a volume in given directions
//===========================================================================
shared_ptr<SplineVolume>
representSurfaceAsVolume(const SplineSurface& surface,
			 int sf_dir1,
			 int sf_dir2,
			 const BsplineBasis& other_bas,
			 bool rational)
//===========================================================================
{
  if (surface.rational())
    THROW("It does not make sense to have a rational hypersurface.");

  int kdim = surface.dimension() / other_bas.numCoefs();

  std::vector<double>::const_iterator co
      = surface.coefs_begin();
  std::vector<double> vol_coefs;
  std::vector<double>::const_iterator coefstart;
  const BsplineBasis& bas1 = surface.basis_u();
  const BsplineBasis& bas2 = surface.basis_v();

  int ijump, jjump, kjump;
  int nu, nv, nw;

  if (sf_dir1 == sf_dir2)
    THROW("Must have two different parameter directions to describe a surface as a volume.");
  if (sf_dir1 < 0 || sf_dir1 > 2)
    THROW("sf_dir1 must be 0, 1 or 2");
  if (sf_dir2 < 0 || sf_dir2 > 2)
    THROW("sf_dir2 must be 0, 1 or 2");

  if (sf_dir1 == 0)
    {
      if (sf_dir2 == 1)
	{
	  nu = surface.numCoefs_u();
	  nv = surface.numCoefs_v();
	  nw = other_bas.numCoefs();
	  kjump = 1;
	  ijump = nw;
	  jjump = nw*nu;
	}
      else   // sf_dir2 == 2
	{
	  nu = surface.numCoefs_u();
	  nv = other_bas.numCoefs();
	  nw = surface.numCoefs_v();
	  jjump = 1;
	  ijump = nv;
	  kjump = nv*nu;
	}
    }
  else if (sf_dir1 == 1)
    {
      if (sf_dir2 == 0)
	{
	  nu = surface.numCoefs_v();
	  nv = surface.numCoefs_u();
	  nw = other_bas.numCoefs();
	  kjump = 1;
	  jjump = nw;
	  ijump = nw*nv;
	}
      else   // sf_dir2 == 2
	{
	  nu = other_bas.numCoefs();
	  nv = surface.numCoefs_u();
	  nw = surface.numCoefs_v();
	  ijump = 1;
	  jjump = nu;
	  kjump = nu*nv;
	}
    }
  else  // sf_dir1 == 2
    {
      if (sf_dir2 == 0)
	{
	  nu = surface.numCoefs_v();
	  nv = other_bas.numCoefs();
	  nw = surface.numCoefs_u();
	  jjump = 1;
	  kjump = nv;
	  ijump = nv*nw;
	}
      else   // sf_dir2 == 1
	{
	  nu = other_bas.numCoefs();
	  nv = surface.numCoefs_v();
	  nw = surface.numCoefs_u();
	  ijump = 1;
	  kjump = nu;
	  jjump = nu*nw;
	}
    }

  vol_coefs.resize(nu*nv*nw*kdim);
  int pos = 0;
  for (int k = 0; k < nw; ++k)
    for (int j = 0; j < nv; ++j)
      for (int i = 0; i < nu; ++i)
	for (int d = 0; d < kdim; ++d)
	  vol_coefs[pos++] = co[(i*ijump + j*jjump + k*kjump)*kdim+d];
  coefstart = vol_coefs.begin();

  int dim = rational ? kdim - 1 : kdim;

  if (sf_dir1 == 0)
    {
      if (sf_dir2 == 1)
	{
	  shared_ptr<SplineVolume> volume
	    (new SplineVolume(bas1, bas2, other_bas,
			      coefstart, dim, rational));
	  return volume;
	}
      else   // sf_dir2 == 2
	{
	  shared_ptr<SplineVolume> volume
	    (new SplineVolume(bas1, other_bas, bas2,
			      coefstart, dim, rational));
	  return volume;
	}
    }
  else if (sf_dir1 == 1)
    {
      if (sf_dir2 == 0)
	{
	  shared_ptr<SplineVolume> volume
	    (new SplineVolume(bas2, bas1, other_bas,
			      coefstart, dim, rational));
	  return volume;
	}
      else   // sf_dir2 == 2
	{
	  shared_ptr<SplineVolume> volume
	    (new SplineVolume(other_bas, bas1, bas2,
			      coefstart, dim, rational));
	  return volume;
	}
    }
  else  // sf_dir1 == 2
    {
      if (sf_dir2 == 0)
	{
	  shared_ptr<SplineVolume> volume
	    (new SplineVolume(bas2, other_bas, bas1,
			      coefstart, dim, rational));
	  return volume;
	}
      else   // sf_dir2 == 1
	{
	  shared_ptr<SplineVolume> volume
	    (new SplineVolume(other_bas, bas2, bas1,
			      coefstart, dim, rational));
	  return volume;
	}
    }
}

  //===========================================================================
  bool getVolAdjacencyInfo(shared_ptr<ParamVolume> vol1,
			   shared_ptr<SurfaceOnVolume> vol_sf1,
			   shared_ptr<ParamVolume> vol2,
			   shared_ptr<SurfaceOnVolume> vol_sf2,
			   double tol,
			   int& bd1, int& bd2, int& orientation,
			   bool& same_seq)
  //===========================================================================
  {
    // bd1, bd2:
    // 0 = umin, 1 = umax, 2 = vmin,  3 = vmax, 4 = wmin, 5 = wmax
    int orientation1, orientation2;
    bool swap1, swap2;
    bd1 = vol_sf1->whichBoundary(tol, orientation1, swap1);
    bd2 = vol_sf2->whichBoundary(tol, orientation2, swap2);
    if (bd1 < 0 || bd2 < 0)
      return false;  // Adjacency not along boundary

    shared_ptr<SplineVolume> svol1 =
      dynamic_pointer_cast<SplineVolume, ParamVolume>(vol_sf1->getVolume());
    shared_ptr<SplineVolume> svol2 =
      dynamic_pointer_cast<SplineVolume, ParamVolume>(vol_sf2->getVolume());
    if (!svol1.get() || !svol2.get())
      return false;  // Cannot handle

    // For the first boundary surface, define parameter values in the
    // volume a little bit inside the corners (to avoid getting
    // undefined situations due to degeneracies)
    const Array<double,6> pardomain1 = svol1->parameterSpan();
    const Array<double,6> pardomain2 = svol2->parameterSpan();
    double parvals1[12], parvals2[12];
    int idx1, idx2, idx3;
    if (bd1 == 0 || bd1 == 1)
      {
	idx1 = 0;
	idx2 = 1;
	idx3 = 2;
      }
    else if (bd1 == 2 || bd1 == 3)
      {
	idx1 = 1;
	idx2 = 0;
	idx3 = 2;
      }
    else
      {
	idx1 = 2;
	idx2 = 0;
	idx3 = 1;
      }
    double tdel2 = 0.1*(pardomain1[2*idx2+1] - pardomain1[2*idx2]);
    double tdel3 = 0.1*(pardomain1[2*idx3+1] - pardomain1[2*idx3]);
    parvals1[idx1] = parvals1[idx1+3] = parvals1[idx1+6] = parvals1[idx1+9] = 
      vol_sf1->getConstVal();
    parvals1[idx2] = parvals1[idx2+9] = pardomain1[2*idx2] + tdel2;
    parvals1[idx2+3] = parvals1[idx2+6] = pardomain1[2*idx2+1] - tdel2;
    parvals1[idx3] = parvals1[idx3+3] = pardomain1[2*idx3] + tdel3;
    parvals1[idx3+6] = parvals1[idx3+9] = pardomain1[2*idx3+1] - tdel3;

    // For the second boundary surface, define parameter values in the
    // volume a little bit inside the corners     
    int idx4, idx5, idx6;
    if (bd2 == 0 || bd2 == 1)
      {
	idx4 = 0;
	idx5 = 1;
	idx6 = 2;
      }
    else if (bd2 == 2 || bd2 == 3)
      {
	idx4 = 1;
	idx5 = 0;
	idx6 = 2;
      }
    else
      {
	idx4 = 2;
	idx5 = 0;
	idx6 = 1;
      }
    double tdel5 = 0.1*(pardomain2[2*idx5+1] - pardomain2[2*idx5]);
    double tdel6 = 0.1*(pardomain2[2*idx6+1] - pardomain2[2*idx6]);
     parvals2[idx4] = parvals2[idx4+3] = parvals2[idx4+6] = parvals2[idx4+9] = 
      vol_sf2->getConstVal();
    parvals2[idx5] = parvals2[idx5+9] = pardomain2[2*idx5] + tdel5;
    parvals2[idx5+3] = parvals2[idx5+6] = pardomain2[2*idx5+1] - tdel5;
    parvals2[idx6] = parvals2[idx6+3] = pardomain2[2*idx6] + tdel6;
    parvals2[idx6+6] = parvals2[idx6+9] = pardomain2[2*idx6+1] - tdel6;

    // Evaluate the volume corresponding to each boundary surface in the 
    // modified corners and inbetween these corners
    Point vp1[8], vp2[8];
    int ki, kj;
    double par[3];
    for (ki=0; ki<4; ++ki)
      {
	kj = (ki+1)%4;
	svol1->point(vp1[2*ki], parvals1[3*ki], parvals1[3*ki+1], 
		     parvals1[3*ki+2]);
	par[idx1] = parvals1[3*ki+idx1];
	if (ki == 0 || ki == 2)
	  {
	    par[idx2] = 0.5*(parvals1[3*ki+idx2] + parvals1[3*kj+idx2]);
	    par[idx3] = parvals1[3*ki+idx3];
	  }
	else
	  {
	    par[idx2] = parvals1[3*ki+idx2];
	    par[idx3] = 0.5*(parvals1[3*ki+idx3] + parvals1[3*kj+idx3]);
	  }
	svol1->point(vp1[2*ki+1], par[0], par[1], par[2]);
      } 

    for (ki=0; ki<4; ++ki)
      {
	kj = (ki+1)%4;
	svol2->point(vp2[2*ki], parvals2[3*ki], parvals2[3*ki+1], 
		     parvals2[3*ki+2]);
	par[idx4] = parvals2[3*ki+idx4];
	if (ki == 0 || ki == 2)
	  {
	    par[idx5] = 0.5*(parvals2[3*ki+idx5] + parvals2[3*kj+idx5]);
	    par[idx6] = parvals2[3*ki+idx6];
	  }
	else
	  {
	    par[idx5] = parvals2[3*ki+idx5];
	    par[idx6] = 0.5*(parvals2[3*ki+idx6] + parvals2[3*kj+idx6]);
	  }
	svol2->point(vp2[2*ki+1], par[0], par[1], par[2]); 
      }
    
    // Compute sum of distances between corresponding points for all 
    // possible configurations of swapping and reversing of parameter
    // directions between the two boundary surfaces
    double dp[8];
    // Same
    dp[0] = vp1[0].dist(vp2[0]) + vp1[1].dist(vp2[1]) + vp1[2].dist(vp2[2]) +
      vp1[3].dist(vp2[3]) + vp1[4].dist(vp2[4]) + vp1[5].dist(vp2[5]) +
      vp1[6].dist(vp2[6]) + vp1[7].dist(vp2[7]);

    // First parameter direction opposite
    dp[1] = vp1[0].dist(vp2[2]) + vp1[1].dist(vp2[1]) + vp1[2].dist(vp2[0]) +
      vp1[3].dist(vp2[7]) + vp1[4].dist(vp2[6]) + vp1[5].dist(vp2[5]) +
      vp1[6].dist(vp2[4]) + vp1[7].dist(vp2[3]);

    // Second parameter direction opposite
    dp[2] = vp1[0].dist(vp2[6]) + vp1[1].dist(vp2[5]) + vp1[2].dist(vp2[4]) +
      vp1[3].dist(vp2[3]) + vp1[4].dist(vp2[2]) + vp1[5].dist(vp2[1]) +
      vp1[6].dist(vp2[0]) + vp1[7].dist(vp2[7]);

    // Both parameter directions opposite
    dp[3] = vp1[0].dist(vp2[4]) + vp1[1].dist(vp2[5]) + vp1[2].dist(vp2[6]) +
      vp1[3].dist(vp2[7]) + vp1[4].dist(vp2[0]) + vp1[5].dist(vp2[1]) +
      vp1[6].dist(vp2[2]) + vp1[7].dist(vp2[3]);

    // Swap parameter directions
    dp[4] = vp1[0].dist(vp2[0]) + vp1[1].dist(vp2[7]) + vp1[2].dist(vp2[6]) +
      vp1[3].dist(vp2[5]) + vp1[4].dist(vp2[4]) + vp1[5].dist(vp2[3]) +
      vp1[6].dist(vp2[2]) + vp1[7].dist(vp2[1]);

    // Swap and opposite first parameter direction
    dp[5] = vp1[0].dist(vp2[6]) + vp1[1].dist(vp2[7]) + vp1[2].dist(vp2[0]) +
      vp1[3].dist(vp2[1]) + vp1[4].dist(vp2[2]) + vp1[5].dist(vp2[3]) +
      vp1[6].dist(vp2[4]) + vp1[7].dist(vp2[5]);

    // Swap and opposite second parameter direction
    dp[6] = vp1[0].dist(vp2[2]) + vp1[1].dist(vp2[3]) + vp1[2].dist(vp2[4]) +
      vp1[3].dist(vp2[5]) + vp1[4].dist(vp2[6]) + vp1[5].dist(vp2[7]) +
      vp1[6].dist(vp2[0]) + vp1[7].dist(vp2[1]);

    // Swap and both parameter directions opposite
    dp[7] = vp1[0].dist(vp2[4]) + vp1[1].dist(vp2[3]) + vp1[2].dist(vp2[2]) +
      vp1[3].dist(vp2[1]) + vp1[4].dist(vp2[0]) + vp1[5].dist(vp2[7]) +
      vp1[6].dist(vp2[6]) + vp1[7].dist(vp2[5]);

    // Find minimum distance. Due to possible squewed parameterizations, we
    // cannot expect it to be zero
    int idm = 0;
    double dmin = dp[0];
    for (ki=1; ki<8; ++ki)
      if (dp[ki] < dmin)
	{
	  idm = ki;
	  dmin = dp[ki];
	}

    // Extract configuration information
    same_seq = (idm < 4) ? true : false;
    orientation = 0;
    if (idm == 1 || idm == 5)
      orientation = 1;
    if (idm == 2 || idm == 6)
      orientation = 2;
    if (idm == 3 || idm == 7)
      orientation = 3;

    return true;
  }


 //===========================================================================
  bool getCorrCoefVolEnum(shared_ptr<SplineVolume> vol1,
			  shared_ptr<SplineVolume> vol2,
			  int bd1, int bd2, int orientation,
			  bool same_seq, 
			  vector<std::pair<int, int> >& enumeration)
 //===========================================================================
  {
    int kn1[3], kn2[3];
    int ki;
    for (ki=0; ki<3; ++ki)
      {
	kn1[ki] = vol1->numCoefs(ki);
	kn2[ki] = vol2->numCoefs(ki);
      }

    int ix1 = bd1/2;
    int ix2 = bd2/2;
    int ki1 = (ix1 == 0) ? 1 : 0;
    int ki2 = (ix1 == 2) ? 1 : 2;
    int kj1 = (ix2 == 0) ? 1 : 0;
    int kj2 = (ix2 == 2) ? 1 : 2;

    if (kn1[ki1] != ((same_seq) ? kn2[kj1] : kn2[kj2]))
      return false;
    if (kn1[ki2] != ((same_seq) ? kn2[kj2] : kn2[kj1]))
      return false;

    enumeration.resize(kn1[ki1]*kn1[ki2]);
    int start1 = (bd1 == 0 || bd1 == 2 || bd1 == 4) ? 0 :
      ((bd1 == 1) ? kn1[0]-1 : ((bd1 == 3) ? kn1[0]*(kn1[1]-1) : 
				kn1[0]*kn1[1]*(kn1[2]-1)));
    int del1_1, del1_2;
    if (ix1 == 0)
      {
	del1_1 = kn1[0];
	del1_2 = kn1[0]*kn1[1];
      }
    else if (ix1 == 1)
      {
	del1_1 = 1;
	del1_2 = kn1[0]*kn1[1];
      }
    else
      {
	del1_1 = 1;
	del1_2 = kn1[0];
      }

    int start2 = (bd2 == 0 || bd2 == 2 || bd2 == 4) ? 0 :
      ((bd2 == 1) ? kn2[0]-1 : ((bd2 == 3) ? kn2[0]*(kn2[1]-1) : 
				kn2[0]*kn2[1]*(kn2[2]-1)));
    int del2_1, del2_2;
    if (ix2 == 0)
      {
	del2_1 = kn2[0];
	del2_2 = kn2[0]*kn2[1];
      }
    else if (ix2 == 1)
      {
	del2_1 = 1;
	del2_2 = kn2[0]*kn2[1];
      }
    else
      {
	del2_1 = 1;
	del2_2 = kn2[0];
      }

    if (!same_seq)
      {
	std::swap(del2_1, del2_2);
	std::swap(kn2[kj1], kn2[kj2]);
      }

    if (orientation == 1 || orientation == 3)
      {
	start2 += del2_1*(kn2[kj1]-1);
	del2_1 *= -1;
      }
    if (orientation == 2 || orientation == 3)
      {
	start2 += del2_2*(kn2[kj2]-1);
	del2_2 *= -1;
      }

    int kr, kh, idx1, idx2;
    for (kr=0; kr<kn1[ki1]*kn1[ki2]; start1+=del1_2, start2+=del2_2)
      for (kh=0, idx1=start1, idx2=start2; kh<kn1[ki1];
	   idx1+=del1_1, idx2+=del2_1, ++kh, ++kr)
	enumeration[kr] = std::make_pair(idx1, idx2);
	
    return true;
  }

 //===========================================================================
bool getVolCoefEnumeration(shared_ptr<SplineVolume> vol, int bd,
			     vector<int>& enumeration) 
//===========================================================================
{
  if (bd < 0 || bd > 5)
    return false;

  int kn[3];
  int ki;
  for (ki=0; ki<3; ++ki)
    kn[ki] = vol->numCoefs(ki);

  int ix = bd/2;
  int ki1 = (ix == 0) ? 1 : 0;
  int ki2 = (ix == 2) ? 1 : 2;

  enumeration.resize(kn[ki1]*kn[ki2]);
  int start = (bd == 0 || bd == 2 || bd == 4) ? 0 :
    ((bd == 1) ? kn[0]-1 : ((bd == 3) ? kn[0]*(kn[1]-1) : 
			    kn[0]*kn[1]*(kn[2]-1)));
  int del1, del2;
  if (ix == 0)
    {
      del1 = kn[0];
      del2 = kn[0]*kn[1];
    }
  else if (ix == 1)
    {
      del1 = 1;
      del2 = kn[0]*kn[1];
    }
  else
    {
      del1 = 1;
      del2 = kn[0];
    }

  int kr, kh, idx;
  for (kr=0; kr<kn[ki1]*kn[ki2]; start+=del2)
    for (kh=0, idx=start; kh<kn[ki1]; idx+=del1, ++kh, ++kr)
      enumeration[kr] = idx;
	
    return true;
}

 //===========================================================================
bool getVolBdCoefEnumeration(shared_ptr<SplineVolume> vol, int bd,
			     int bd_cv, vector<int>& enumeration) 
//===========================================================================
{
  if (bd < 0 || bd > 5)
    return false;
  if (bd_cv < 0 && bd_cv > 3)
    return false;

  int kn[3];
  int ki;
  for (ki=0; ki<3; ++ki)
    kn[ki] = vol->numCoefs(ki);

  int ix = bd/2;
  int ki1 = (ix == 0) ? 1 : 0;
  int ki2 = (ix == 2) ? 1 : 2;
  int kj = (bd_cv <= 1) ? ki2 : ki1;
  enumeration.resize(kn[kj]);

  // Set pointer to boundary surface
  int start = (bd == 0 || bd == 2 || bd == 4) ? 0 :
    ((bd == 1) ? kn[0]-1 : ((bd == 3) ? kn[0]*(kn[1]-1) : 
			    kn[0]*kn[1]*(kn[2]-1)));

  // Adjust pointer to curve on surface and set
  int del1, del2, del;
  if (ix == 0)
    {
      del1 = kn[0];
      del2 = kn[0]*kn[1];
      if (bd_cv == 1)
	start += kn[0]*(kn[1]-1);
      else if (bd_cv == 3)
	start += kn[0]*kn[1]*(kn[2]-1);
    }
  else if (ix == 1)
    {
      del1 = 1;
      del2 = kn[0]*kn[1];
      if (bd_cv == 1)
	start += kn[0]-1;
      else if (bd_cv == 3)
	start += kn[0]*kn[1]*(kn[2]-1);
     }
  else
    {
      del1 = 1;
      del2 = kn[0];
      if (bd_cv == 1)
	start += kn[0]-1;
      else if (bd_cv == 3)
	start += kn[0]*(kn[1]-1);
    }
  del = (bd_cv <= 1) ? del2 : del1;

   int kh, idx;
  for (kh=0, idx=start; kh<kn[kj]; idx+=del, ++kh)
    enumeration[kh] = idx;
	
    return true;
}

 //===========================================================================
   vector<shared_ptr<ParamSurface> > getBoundarySurfaces(shared_ptr<ParamVolume> vol)
 //===========================================================================
   {
    vector<shared_ptr<ParamSurface> > bd_sfs = vol->getAllBoundarySurfaces();

    if (bd_sfs.size() != 6)
      return bd_sfs;   

    const Array<double,6> params = vol->parameterSpan();

    vector<shared_ptr<ParamSurface> > sfs;
    for (size_t ki=0; ki<bd_sfs.size(); ++ki)
      {
	int dir = (int)ki/2 + 1;
	shared_ptr<ParamSurface> curr_bd = 
	  shared_ptr<ParamSurface>(new SurfaceOnVolume(vol, bd_sfs[ki],
						       dir, params[(int)ki], (int)ki,
						       false));
	sfs.push_back(curr_bd);
      }
    return sfs;
  }

 //===========================================================================
   vector<shared_ptr<ParamSurface> > 
   getOrientedBoundarySurfaces(shared_ptr<ParamVolume> vol)
 //===========================================================================
   {
    vector<shared_ptr<ParamSurface> > bd_sfs = vol->getAllBoundarySurfaces();
    // bd_sfs[1]->swapParameterDirection();
    // bd_sfs[2]->swapParameterDirection();
    // bd_sfs[5]->swapParameterDirection();
    bd_sfs[0]->swapParameterDirection();
    bd_sfs[3]->swapParameterDirection();
    bd_sfs[4]->swapParameterDirection();

#ifdef DEBUG
    shared_ptr<SplineVolume> vol2 = 
      dynamic_pointer_cast<SplineVolume,ParamVolume>(vol);
    if (vol2)
      {
	bool lefthanded = vol2->isLeftHanded();
	std::cout << "Volume lefthanded: " << lefthanded << std::endl;
      }
    std::ofstream of("volume_boundaries.g2");
    for (int kr=0; kr<6; ++kr)
      {
	bd_sfs[kr]->writeStandardHeader(of);
	bd_sfs[kr]->write(of);
      }
#endif

    if (bd_sfs.size() != 6)
      return bd_sfs;   

    const Array<double,6> params = vol->parameterSpan();

    vector<shared_ptr<ParamSurface> > sfs;
    for (size_t ki=0; ki<bd_sfs.size(); ++ki)
      {
	int dir = (int)ki/2 + 1;
	//bool swap = (ki == 1 || ki == 2 || ki == 5);
	bool swap = (ki == 0 || ki == 3 || ki == 4);
	shared_ptr<ParamSurface> curr_bd = 
	  shared_ptr<ParamSurface>(new SurfaceOnVolume(vol, bd_sfs[ki],
						       dir, params[(int)ki], (int)ki,
						       swap));
	sfs.push_back(curr_bd);
      }
    return sfs;
  }

 //===========================================================================
   shared_ptr<SurfaceOnVolume> 
   getBoundarySurface(shared_ptr<SplineVolume> vol, int idx)
 //===========================================================================
   {
     shared_ptr<SplineSurface> bd_sf = vol->getBoundarySurface(idx);

    const Array<double,6> params = vol->parameterSpan();

    int dir = idx/2 + 1;
    shared_ptr<SurfaceOnVolume> volsf = 
      shared_ptr<SurfaceOnVolume>(new SurfaceOnVolume(vol, bd_sf,
						      dir, params[idx], idx,
						      false));
    return volsf;
  }

 //===========================================================================
   shared_ptr<SurfaceOnVolume> 
   getOrientedBoundarySurface(shared_ptr<SplineVolume> vol, int idx)
 //===========================================================================
   {
     shared_ptr<SplineSurface> bd_sf = vol->getBoundarySurface(idx);
     bool swap = false;
     //if (idx == 1 || idx == 2 || idx == 5)
     if (idx == 0 || idx == 3 || idx == 4)
       {
	 bd_sf->swapParameterDirection();
	 swap = true;
       }

    const Array<double,6> params = vol->parameterSpan();

    int dir = idx/2 + 1;
    shared_ptr<SurfaceOnVolume> volsf = 
      shared_ptr<SurfaceOnVolume>(new SurfaceOnVolume(vol, bd_sf,
						      dir, params[idx], idx,
						      swap));
    return volsf;
  }
  //===========================================================================
 void volCommonSplineSpace(shared_ptr<SplineVolume> vol1, int bd1,
			   shared_ptr<SplineVolume> vol2, int bd2,
			   int orientation, bool same_seq)
 //===========================================================================
 {
   // Reorient vol2 to get the same orientation as vol1. First represent
   // the volumes as surfaces
   int dir1 = (bd1 > 1) ? 0 : 1;
   int dir2 = (bd1 >= 4) ? 1 : 2;
   int dir3 = (bd2 > 1) ? 0 : 1;
   int dir4 = (bd2 >= 4) ? 1 : 2;
   BsplineBasis basis1 = (bd1 <= 1) ? vol1->basis(0) : 
     ((bd1 <= 3) ? vol1->basis(1) : vol1->basis(2));
   BsplineBasis basis2 = (bd2 <= 1) ? vol2->basis(0) : 
     ((bd2 <= 3) ? vol2->basis(1) : vol2->basis(2));
   bool rat1 = vol1->rational();
   bool rat2 = vol2->rational();

   shared_ptr<SplineSurface> srf1 = representVolumeAsSurface(*vol1,
							     dir1, dir2);
   shared_ptr<SplineSurface> srf2 = representVolumeAsSurface(*vol2,
							     dir3, dir4);

   // Ensure same sequence of parameter directions
   if (!same_seq)
     srf2->swapParameterDirection();

   // Ensure same orientation of parameter directions
   if (orientation == 1 || orientation == 3)
     srf2->reverseParameterDirection(true);

   if (orientation == 2 || orientation == 3)
     srf2->reverseParameterDirection(false);
   
   // Set same parameter domain
    double umin1 = srf1->startparam_u();
    double vmin1 = srf1->startparam_v();
    double umax1 = srf1->endparam_u();
    double vmax1 = srf1->endparam_v();
    double umin2 = srf2->startparam_u();
    double vmin2 = srf2->startparam_v();
    double umax2 = srf2->endparam_u();
    double vmax2 = srf2->endparam_v();
    double ptol = 1.0e-12;
    if (fabs(umin1-umin2) > ptol || fabs(umax1-umax2) > ptol)
    {
	double umin = 0.5*(umin1 + umin2);
	double umax = 0.5*(umax1 + umax2);
	double vmin = 0.5*(vmin1 + vmin2);
	double vmax = 0.5*(vmax1 + vmax2);
    
	srf1->setParameterDomain(umin, umax, vmin, vmax);
	srf2->setParameterDomain(umin, umax, vmin, vmax);
    }
   
    // Ensure the same spline spaces
    vector<shared_ptr<SplineSurface> > sfs(2);
    sfs[0] = srf1;
    sfs[1] = srf2;
    unifySurfaceSplineSpace(sfs, ptol, 0);
    srf1 = sfs[0];
    srf2 = sfs[1];

    // Reverse all translations
    srf1->setParameterDomain(umin1, umax1, vmin1, vmax1);
    srf2->setParameterDomain(umin2, umax2, vmin2, vmax2);

    if (orientation == 1 || orientation == 3)
      srf2->reverseParameterDirection(true);

    if (orientation == 2 || orientation == 3)
      srf2->reverseParameterDirection(false);
    
    if (!same_seq)
      srf2->swapParameterDirection();

    // Represent surfaces as volumes
    shared_ptr<SplineVolume> vol3 = 
      representSurfaceAsVolume(*srf1, dir1, dir2, basis1, rat1);
    shared_ptr<SplineVolume> vol4 = 
      representSurfaceAsVolume(*srf2, dir3, dir4, basis2, rat2);
    vol1->swap(*vol3);
    vol2->swap(*vol4);
 }

//===========================================================================
shared_ptr<SplineCurve> liftVolParamCurve(shared_ptr<ParamCurve> pcurve, 
					  shared_ptr<ParamVolume> vol,
					  double tol)
//===========================================================================
{
  // Represent the parameter curve as an evaluator based curve
  shared_ptr<VolumeSpaceCurve> vol_space =
    shared_ptr<VolumeSpaceCurve>(new VolumeSpaceCurve(vol, pcurve));

  // Approximate
  HermiteAppC approx(vol_space.get(), tol, tol);
  approx.refineApproximation();
  return approx.getCurve();
}

//===========================================================================
shared_ptr<SplineCurve> projectVolParamCurve(shared_ptr<ParamCurve> spacecurve, 
					     shared_ptr<ParamVolume> vol,
					     double tol)
//===========================================================================
{
  // Represent the parameter curve as an evaluator based curve
  shared_ptr<VolumeParameterCurve> vol_par =
    shared_ptr<VolumeParameterCurve>(new VolumeParameterCurve(vol, 
							      spacecurve));

  // Approximate
  HermiteAppC approx(vol_par.get(), tol, tol);
  approx.refineApproximation();
  shared_ptr<SplineCurve> crv = approx.getCurve();
  if (!crv.get())
    {
      // Approximation failed. Check for seam curves
      // Check if the volume is closed in any direction
      shared_ptr<SplineVolume> vol2 = 
	dynamic_pointer_cast<SplineVolume, ParamVolume>(vol);
      if (!vol2.get())
	return crv;   // Should not happen

      int closed[3];
      int ki;
      for (ki=0; ki<3; ++ki)
	closed[ki] = vol2->volumePeriodicity(ki, tol);

      // Make candidate start and end points
      vector<Point> startpt;
      vector<Point> endpt;
      
      Point pos1 = vol_par->eval(vol_par->start());
      Point pos2 = vol_par->eval(vol_par->end());

      for (ki=0; ki<3; ++ki)
	{
	  if (closed[ki] < 0)
	    continue;

	  if (pos1[ki] - vol2->startparam(ki) < tol)
	    {
	      Point pos = pos1;
	      pos[ki] = vol2->endparam(ki);
	      startpt.push_back(pos);
	    }

	  if (pos2[ki] - vol2->startparam(ki) < tol)
	    {
	      Point pos = pos2;
	      pos[ki] = vol2->endparam(ki);
	      endpt.push_back(pos);
	    }
	  
 	  if (vol2->endparam(ki) - pos1[ki] < tol)
	    {
	      Point pos = pos1;
	      pos[ki] = vol2->startparam(ki);
	      startpt.push_back(pos);
	    }

	  if (vol2->endparam(ki) - pos2[ki] < tol)
	    {
	      Point pos = pos2;
	      pos[ki] = vol2->startparam(ki);
	      endpt.push_back(pos);
	    }
	}
      startpt.push_back(pos1);
      endpt.push_back(pos2);
	  
      for (size_t kj=0; kj<startpt.size(); ++kj)
	{
	  shared_ptr<Point> pt1(new Point(startpt[kj].begin(), startpt[kj].end()));
	  for (size_t kr=0; kr<endpt.size(); ++kr)
	    {
	      shared_ptr<Point> pt2(new Point(endpt[kr].begin(), endpt[kr].end()));
	      shared_ptr<VolumeParameterCurve> vol_par2 =
		shared_ptr<VolumeParameterCurve>(new VolumeParameterCurve(vol, 
									  spacecurve,
									  pt1, pt2));
	      // Approximate
	      HermiteAppC approx2(vol_par2.get(), tol, tol);
	      approx2.refineApproximation();
	      shared_ptr<SplineCurve> crv2 = approx2.getCurve();
	      if (crv2.get())
		return crv2;
	    }
	}
	      
    }
  return crv;
}

} // namespace Go
