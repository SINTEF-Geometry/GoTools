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


using namespace Go;
using std::vector;
using std::max;
using std::min;

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
