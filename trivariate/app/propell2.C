#include <vector>
#include <fstream>
#include <string.h>
#include "GoTools/utils/Point.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/creators/ApproxCurve.h"
#include "GoTools/creators/ApproxCrvToSeqs.h"
#include "GoTools/trivariate/LoftVolumeCreator.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "sisl.h"

using std::ifstream;
using std::ofstream;
using std::endl;
using std::vector;
using namespace Go;

void parameterizePointSequences(vector<vector<Point> >& pnt_seq,
				vector<int>& mid_idx,
				vector<vector<double> >& param)
{
  param.resize(pnt_seq.size());

  // First perform independent chord lenght parameterization for
  // each point sequence
  size_t ki;
  for (ki=0; ki<pnt_seq.size(); ++ki)
    {
      vector<double> parvals(pnt_seq[ki].size());
      parvals[0] = 0.0;
      for (size_t kj=1; kj<pnt_seq[ki].size(); ++kj)
	{
	  double dist = pnt_seq[ki][kj-1].dist(pnt_seq[ki][kj]);
	  parvals[kj] = parvals[kj-1] + dist;
	}
      param[ki] = parvals;
    }

  // Compute average values for the start, mid and end parameter
  double start = 0.0;
  double mid = 0.0;
  double end = 0.0;
  for (ki=0; ki<param.size(); ++ki)
    {
      start += param[ki][0];
      mid += param[ki][mid_idx[ki]];
      end += param[ki][param[ki].size()-1];
    }
  start /= (double)(param.size());
  mid /= (double)(param.size());
  end /= (double)(param.size());

  // Reparameterize to get common start, mid and end parameter
  for (ki=0; ki<param.size(); ++ki)
    {
      int kj;
      double ta = param[ki][0];
      double tb = param[ki][mid_idx[ki]];
      double tc = param[ki][param[ki].size()-1];
      for (kj=0; kj<mid_idx[ki]; ++kj)
	param[ki][kj] = start + (param[ki][kj] - ta)*(mid - start)/(tb - ta);
      for (kj=mid_idx[ki]; kj<(int)param[ki].size(); ++kj)
	param[ki][kj] = mid + (param[ki][kj] - tb)*(end - mid)/(tc - tb);
    }
}

int main(int argc, char* argv[] )
{

  if (argc != 7)
    {
      std::cout << "Usage: infile tolerance c1approxfac split_idx approx_blade outfile" << endl;
      //      exit(-1);
    }

  // Open inputgeometry file
  ifstream is(argv[1]);

  double eps = atof(argv[2]);
  double c1 = atof(argv[3]);
  int split_idx = atoi(argv[4]);
  int approx_blade = atoi(argv[5]);

  // Open outfile
  ofstream os(argv[6]);

  vector<double> height;
  vector<double> twist;
  vector<double> chord;
  vector<double> cent;
  vector<double> orig;
  vector<int> newset;

  vector<vector<Point> > pnt_seq;
  vector<int> mid_idx;
  // Read 2D data
  char prev_file[50];
  while (!is.eof())
    {
      double rad, tw, chd, cen, org;
      is >> rad >> tw >> chd >> cen >> org;
      height.push_back(rad);
      twist.push_back(tw);
      chord.push_back(chd);
      cent.push_back(cen);
      orig.push_back(org);

      char file[50];
      is >> file;
#ifndef _MSC_VER
      int compare = strcasecmp(prev_file, file);
#else
      int compare = _stricmp(prev_file, file);
#endif // _MSC_VER
      std::cout << file << ", compare: " << compare << std::endl;
      std::ifstream is2(file);

      if (compare != 0)
	{
	  vector<Point> seq;
	  double xmin = 1.0e8;
	  double xmax = -1.0e8;
	  int max_idx = 0;
	  int ki=0;
	  while (!is2.eof())
	    {
	      double x, y;
	      is2 >> x >> y;
	      xmin = std::min(xmin, x);
	      if (x > xmax)
		{
		  xmax = x;
		  max_idx = ki;
		}
	      Point pnt(x, y);
	      seq.push_back(pnt);
	      ki++;
	    }
	  vector<Point> seq2;
	  seq2.insert(seq2.end(), seq.begin()+max_idx, seq.end()-1);
	  mid_idx.push_back((int)seq2.size()-1);
	  seq2.insert(seq2.end(), seq.begin()+1, seq.begin()+max_idx+1);
	  pnt_seq.push_back(seq2);
	  newset.push_back(1);
	}
      else
	newset.push_back(0);

      memcpy(prev_file, file, (size_t)(50*sizeof(char)));
    }
  
  // Parameterize point sequences
  vector<vector<double> > param;
  parameterizePointSequences(pnt_seq, mid_idx, param);
  double split_par = param[0][mid_idx[0]];

  // Make 2D curves in standard position
  // First make initial knot vector
  int in = 5; //8;
  int ik = 4;
  vector<double> knots(in+ik);
  for (int ki=0; ki<ik; ++ki)
    knots[ki] = param[0][0];
  //for (int ki=0; ki<ik; ++ki)
    knots[ik] = param[0][mid_idx[0]];
  for (int ki=0; ki<ik; ++ki)
    knots[in+ki] = param[0][param[0].size()-1];

  // bool reset = true;
  vector<vector<double> > all_pnts;
  for (size_t kj=0; kj<pnt_seq.size(); ++kj)
    //for (int kj=pnt_seq.size()-1; kj>=0; --kj)
    {
      vector<double> pnt;
      for (size_t kr=0; kr<pnt_seq[kj].size(); ++kr)
	pnt.insert(pnt.end(), pnt_seq[kj][kr].begin(), pnt_seq[kj][kr].end());
      all_pnts.push_back(pnt);
    }
  
  ApproxCrvToSeqs approx(all_pnts, param, 2, eps, in, ik, knots);

  approx.unsetSmooth();
  approx.setC1Approx(c1);

  double maxdist, avdist;
  int iter = 8;
  vector<shared_ptr<SplineCurve> > tmp_crvs = 
    approx.getApproxCurves(maxdist, avdist, iter);
  std::cout << "Curve approx, max dist: " << maxdist;
  std::cout << ", average: " << avdist << std::endl;

  // Transform curves to correct position
  vector<shared_ptr<SplineCurve> > section_crvs(height.size());
  
  shared_ptr<SplineCurve> curr;
  size_t kj;
  int idx;
  for (kj=0, idx=-1; kj<height.size(); ++kj)
    {
      if (newset[kj])
	curr = tmp_crvs[++idx];

      // Make transformed coefficients. The curve is non-rational
      Point med = 0.5*(pnt_seq[idx][0] - pnt_seq[idx][mid_idx[idx]]);
      int in = curr->numCoefs();
      int dim = curr->dimension();
      vector<double>::iterator ctr = curr->coefs_begin();
      vector<double> coefs(in*(dim+1));
      for (int kh=0; kh<in; kh++, ctr+=dim)
	{
	  Point pnt0 = Point(*ctr, *(ctr+1));
	  pnt0 -= med;
	  pnt0 *= chord[kj];
	  pnt0 += med;
	  coefs[kh*(dim+1)] = pnt0[0];
	  coefs[kh*(dim+1)+1] = pnt0[1];
	  coefs[kh*(dim+1)+2] = height[kj];
	}
      section_crvs[kj] = 
	shared_ptr<SplineCurve>(new SplineCurve(in, curr->order(), curr->knotsBegin(),
						coefs.begin(), dim+1));
      section_crvs[kj]->writeStandardHeader(os);
      section_crvs[kj]->write(os);
    }
      
      

  double tol = 1.0e-4;
  GeometryTools::unifyCurveSplineSpace(section_crvs, tol);

  // Make lofted surfaces
  vector<SISLCurve*> sisl_cvs(section_crvs.size());
  vector<int> ntype(section_crvs.size(), 1);
  for (kj=0; kj<section_crvs.size(); ++kj)
    sisl_cvs[kj] = Curve2SISL(*section_crvs[kj], false);
  int status = 0;
  double *gpar1 = NULL;
  double *gpar2 = NULL;
  double *gpar3 = NULL;
  SISLSurf *sisl_surf1 = NULL;
  SISLSurf *sisl_surf2 = NULL;
  SISLSurf *sisl_surf3 = NULL;
  shared_ptr<SplineSurface> surf1;
  shared_ptr<SplineSurface> surf2;
  shared_ptr<SplineSurface> surf3;
  double start_par;

  int nmb1 = (split_idx > 1 && split_idx < (int)sisl_cvs.size()-2) ? split_idx : 
      (int)sisl_cvs.size();
  s1538(nmb1, &sisl_cvs[0], &ntype[0], 0.0, 1, 4, 
	0, &sisl_surf1, &gpar1, &status);
//   s1538(sisl_cvs.size(), &sisl_cvs[0], &ntype[0], 0.0, 1, 3, 
// 	0, &sisl_surf, &gpar, &status);

    if (status >= 0)
    {
      surf1 = shared_ptr<SplineSurface>(SISLSurf2Go(sisl_surf1));
      surf1->writeStandardHeader(os);
      surf1->write(os);
    }

    if (nmb1 < (int)sisl_cvs.size())
    {
      if (approx_blade)
	{
	  vector<double> cv_par;
	  cv_par.insert(cv_par.end(), height.begin()+nmb1, 
			height.begin()+section_crvs.size());
	  double tol = 1.0e-4;
	  vector<shared_ptr<SplineCurve> > blade_crvs;
	  blade_crvs.insert(blade_crvs.end(), section_crvs.begin()+nmb1,
			    section_crvs.end());
	  GeometryTools::unifyCurveSplineSpace(blade_crvs, tol);
	  vector<double> cv_coefs;
	  for (size_t kj=0; kj<blade_crvs.size(); ++kj)
	    cv_coefs.insert(cv_coefs.end(), blade_crvs[kj]->coefs_begin(),
			    blade_crvs[kj]->coefs_end());

	  int cvdim = blade_crvs[0]->numCoefs()*section_crvs[0]->dimension();
	  ApproxCurve approxsf(cv_coefs, cv_par, cvdim, 10.0*eps, 4, 4);
	  double maxdist, avdist;
	  int iter = 3;
	  shared_ptr<SplineCurve> crvsf = 
	    approxsf.getApproxCurve(maxdist, avdist, iter);
	  std::cout << "Surface curve, max dist: " << maxdist;
	  std::cout << ", average: " << avdist << std::endl;

	  surf2 =
	      GeometryTools::representCurveAsSurface(*crvsf, 2, blade_crvs[0]->basis(), false);
	  start_par = height[nmb1];
 	}
      else
	{
	    s1538((int)sisl_cvs.size()-nmb1, &sisl_cvs[nmb1], &ntype[nmb1], 0.0, 1, 4, 
		0, &sisl_surf2, &gpar2, &status);
	  if (status >= 0)
	    {
	      surf2 = shared_ptr<SplineSurface>(SISLSurf2Go(sisl_surf2));
	      start_par = gpar2[0];
	    }
	}
      surf2->writeStandardHeader(os);
      surf2->write(os);
    }
      
  for (kj=0; kj<sisl_cvs.size(); ++kj)
    freeCurve(sisl_cvs[kj]);

  if (surf2.get())
    {
      vector<SISLCurve*> bd_cvs(4);
      SplineCurve *bd1, *cross1, *bd2, *cross2;
      surf1->constParamCurve(gpar1[nmb1-1], true, bd1, cross1);
      surf2->constParamCurve(start_par, true, bd2, cross2);

      bd_cvs[0] = Curve2SISL(*bd1, false);
      bd_cvs[1] = Curve2SISL(*cross1, false);
      bd_cvs[2] = Curve2SISL(*bd2, false);
      bd_cvs[3] = Curve2SISL(*cross2, false);
      vector<int> ntype2(4);
      ntype2[0] = ntype2[2] = 1;
      ntype2[1] = ntype2[3] = 4;
      // double epar[2];
      s1538(4, &bd_cvs[0], &ntype2[0], 0.0, 1, 4, 
	    0, &sisl_surf3, &gpar3, &status);
      if (status >= 0)
	{
	  surf3 = shared_ptr<SplineSurface>(SISLSurf2Go(sisl_surf3));
	  surf3->writeStandardHeader(os);
	  surf3->write(os);
	}

      delete bd1;
      delete bd2;
      delete cross1;
      delete cross2;
      for (int ki=0; ki<4; ++ki)
	freeCurve(bd_cvs[ki]);
    }

  if (sisl_surf1) freeSurf(sisl_surf1);
  if (sisl_surf2) freeSurf(sisl_surf2);
  if (sisl_surf3) freeSurf(sisl_surf3);
  if (gpar1) free(gpar1);
  if (gpar2) free(gpar2);
  if (gpar2) free(gpar3);

  if (surf2.get())
    {
      shared_ptr<SplineSurface> surf4 = shared_ptr<SplineSurface>(surf1->clone());

      double dist1, dist2;
      surf1->appendSurface((ParamSurface*)surf3.get(), 2, 1, dist1);
      std::cout << "Frist append, attempt 1: " << dist1 << std::endl;
      surf1->appendSurface((ParamSurface*)surf2.get(), 2, 1, dist1);
      std::cout << "Second append, attempt 1: " << dist1 << std::endl;
      surf1->writeStandardHeader(os);
      surf1->write(os);

      surf4->appendSurface((ParamSurface*)surf3.get(), 2, 2, dist2);
      std::cout << "Frist append, attempt 2: " << dist2 << std::endl;
      surf4->appendSurface((ParamSurface*)surf2.get(), 2, 2, dist2);
      std::cout << "Second append, attempt 2: " << dist2 << std::endl;
      surf4->writeStandardHeader(os);
      surf4->write(os);
    }

  // Make boundary surfaces for the volume
  shared_ptr<SplineSurface> sub1 = 
    shared_ptr<SplineSurface>(surf1->subSurface(surf1->startparam_u(),
						surf1->startparam_v(),
						split_par,
						surf1->endparam_v()));
  shared_ptr<SplineSurface> sub2 = 
    shared_ptr<SplineSurface>(surf1->subSurface(split_par,
						surf1->startparam_v(),
						surf1->endparam_u(),
						surf1->endparam_v()));
  sub2->reverseParameterDirection(true);

  // Create volume by loft
  vector<shared_ptr<SplineSurface> > bd_sfs(2);
  bd_sfs[0] = sub1;
  bd_sfs[1] = sub2;
  shared_ptr<SplineVolume> vol = 
    shared_ptr<SplineVolume>(LoftVolumeCreator::loftVolume(bd_sfs.begin(),
							   2));

//   vol->writeStandardHeader(os);
//   vol->write(os);

  std::ofstream of2("propeller_vol.g2");
  vol->writeStandardHeader(of2);
  vol->write(of2);
}

      
