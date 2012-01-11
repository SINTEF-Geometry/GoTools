//===========================================================================
//                                                                           
// File: ModifySurf.C                                                     
//                                                                           
// Created: Mars 2010
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/creators/ModifySurf.h"
#include "GoTools/geometry/RectDomain.h"
#include "GoTools/creators/ApproxSurf.h"
#include "GoTools/creators/SmoothSurfSet.h"
#include "GoTools/geometry/Utils.h"

using namespace Go;
using std::vector;
using std::max;
using std::min;
using std::make_pair;

//===========================================================================
void ModifySurf::replaceBoundary(shared_ptr<SplineSurface> surf,
				 shared_ptr<SplineCurve> curve,
				 int bd_idx, double tol)
//===========================================================================
{
  // Sample points in the initial surface
  RectDomain dom = surf->parameterDomain();
  int dim = surf->dimension();

  vector<double> sample_pts0;
  vector<double> sample_pts;
  vector<double> sample_pars;
  vector<double> par_u, par_v;
  int nmb_sample = 20;
  surf->gridEvaluator(nmb_sample, nmb_sample, sample_pts0, par_u, par_v);

  // Replace boundary curve
  surf->replaceBoundaryCurve(bd_idx, curve);

  // Collect approxation points
  sample_pars.reserve(2*nmb_sample*nmb_sample);
  sample_pts.reserve(sample_pts0.size());
  int ki, kj, kr;
  for (kj=2, kr=kj*nmb_sample; kj<nmb_sample-2; ++kj, kr=kj*nmb_sample)
    for (ki=2, kr+=2; ki<nmb_sample-2; ++ki, ++kr)
      {
	Point curr(sample_pts0.begin()+kr*dim, sample_pts0.begin()+(kr+1)*dim);
	double u1, v1, d1;
	Point pnt2;
	double seed[2];
	seed[0] = par_u[ki];
	seed[1] = par_v[kj];
	surf->closestPoint(curr, u1, v1, pnt2, d1, tol, NULL, seed);

	sample_pars.push_back(u1);
	sample_pars.push_back(v1);
	sample_pts.insert(sample_pts.end(), sample_pts0.begin()+kr*dim,
			  sample_pts0.begin()+(kr+1)*dim);
      }

  static bool mod_inner = true;
  if (mod_inner)
    {
  // Let the modified surface approximate the sample points. Parameter 
  // iteration is performed
  ApproxSurf approx(surf, sample_pts, sample_pars, surf->dimension(), tol);

  // Keep the boundary
  approx.setFixBoundary(true);

  // Do not refine surfaces
  approx.setDoRefine(false);

  double maxdist, avdist;
  int nmb_out;
  shared_ptr<SplineSurface> mod_srf = approx.getApproxSurf(maxdist, avdist, nmb_out);
				
  if (mod_srf.get())
    surf->swap(*mod_srf);
    }
}


//===========================================================================
bool 
ModifySurf::enforceCoefCoLinearity(shared_ptr<SplineSurface> sf1, int bd1, 
				   shared_ptr<SplineSurface> sf2, int bd2, 
				   double tol, 
				   std::vector<std::vector<int> >& enumeration)
//===========================================================================
{
  // Collect surfaces
  int dim = sf1->dimension();
  vector<double>::iterator c1 = sf1->coefs_begin();
  vector<double>::iterator c2 = sf2->coefs_begin();
  vector<shared_ptr<SplineSurface> > sfs(2);
  sfs[0] = sf1;
  sfs[1] = sf2;

  // Define continuity constraints
  vector<sideConstraintSet> constraints(2*enumeration.size());
  for (size_t ki=0; ki<enumeration.size(); ++ki)
    {
      // Equality of coefficients at the boundary
      sideConstraintSet curr(dim);
      curr.factor_.push_back(make_pair(make_pair(0, enumeration[ki][1]), 1.0));
      curr.factor_.push_back(make_pair(make_pair(1, enumeration[ki][2]), -1.0));
      constraints[2*ki] = curr;

      // Co-linearity of coefficients at the boundary and the adjacent
      // rows
      // Compute distance between coeffients on either side of the
      // boundary
      double d1 = sqrt(distance_squared(&c1[enumeration[ki][0]*dim],
					&c1[(enumeration[ki][0]+1)*dim],
					&c1[enumeration[ki][1]*dim]));
      double d2 = sqrt(distance_squared(&c2[enumeration[ki][2]*dim],
					&c2[(enumeration[ki][2]+1)*dim],
					&c2[enumeration[ki][3]*dim]));
      sideConstraintSet curr2(dim);
      curr2.factor_.push_back(make_pair(make_pair(0, enumeration[ki][0]), -d2));
      curr2.factor_.push_back(make_pair(make_pair(0, enumeration[ki][1]), d2));
      curr2.factor_.push_back(make_pair(make_pair(1, enumeration[ki][2]), d1));
      curr2.factor_.push_back(make_pair(make_pair(1, enumeration[ki][3]), -d1));
      constraints[2*ki+1] = curr2;
    }

  // Define which coefficients that may change
  vector<vector<int> > coef_known(2);
  int kn1 = sf1->numCoefs_u();
  int kn2 = sf1->numCoefs_v();
  int kn3 = sf2->numCoefs_u();
  int kn4 = sf2->numCoefs_v();

  // 1. surface
  vector<int> ck1(kn1*kn2, 1.0); // Initially all coefficients are fixed

  // Release coefficients close to the boundary
  int min_nmb = 2;
  int max_nmb = (bd1 == 0 || bd1 == 1) ? kn1 - 2 : kn2 - 2;
  int nmb = 4;
  nmb = std::max(min_nmb, std::min(max_nmb, nmb));
  if (nmb > max_nmb)
    return false;  // Not enough degrees of freedom

  int ki, kj;
  int start_u = (bd1 == 1) ? kn1-nmb-1 : 0;
  int start_v = (bd1 == 3) ? kn2-nmb-1 : 0;
  int end_u = (bd1 == 0) ? nmb : kn1;
  int end_v = (bd1 == 2) ? nmb : kn2;
  if (bd1 == 0 || bd1 == 1)
    {
      start_v += 1;
      end_v -= 1;
    }
  else
    {
      start_u += 1;
      end_u -= 1;
    }
  for (kj=start_v; kj<end_v; ++kj)
    for (ki=start_u; ki<end_u; ++ki)
      ck1[kj*kn1+ki] = 0;

  coef_known[0] = ck1;

  // 2. surface
  vector<int> ck2(kn3*kn4, 1.0); // Initially all coefficients are fixed

  // Release coefficients close to the boundary
  max_nmb = (bd2 == 0 || bd2 == 1)  ? kn3 - 2 : kn4 - 2;
  nmb = 4;
  nmb = std::max(min_nmb, std::min(max_nmb, nmb));
  if (nmb > max_nmb)
    return false;  // Not enough degrees of freedom

  start_u = (bd2 == 1) ? kn3-nmb-1 : 0;
  start_v = (bd2 == 3) ? kn4-nmb-1 : 0;
  end_u = (bd2 == 0) ? nmb : kn3;
  end_v = (bd2 == 2) ? nmb : kn4;
  if (bd2 == 0 || bd2 == 1)
    {
      start_v += 1;
      end_v -= 1;
    }
  else
    {
      start_u += 1;
      end_u -= 1;
    }
  for (kj=start_v; kj<end_v; ++kj)
    for (ki=start_u; ki<end_u; ++ki)
      ck2[kj*kn3+ki] = 0;

  coef_known[1] = ck2;

  // Perform smoothing across the boundary
  // First set weights
  double w1 = 0.0;
  double w2 = 0.0025;
  double w3 = 0.0025;
  double w_orig = 0.05;

  SmoothSurfSet smooth(true);
  smooth.attach(sfs, coef_known, (int)constraints.size());

  smooth.setOptimize(w1, w2, w3);

  smooth.approxOrig(w_orig);

  smooth.setSideConstraints(constraints);

  vector<shared_ptr<SplineSurface> > sfs2;
  int stat = smooth.equationSolve(sfs2);
  if (stat < 0)
    return false;

  // Modify input surfaces
  sf1->swap(*sfs2[0]);
  sf2->swap(*sfs2[1]);
  return true;
}


