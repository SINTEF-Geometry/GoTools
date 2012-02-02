#include <fstream>
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/trivariate/LoftVolumeCreator.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/trivariate/VolumeTools.h"
#include <memory>

using namespace Go;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::vector;

int main(int argc, char* argv[] )
{
  if (argc != 11)
      cout << "Usage: " << "<infile> <hull thickness> <nmb stiffeners> <nmb stiffeners2> <height fraction> <length> <thickness fraction> <split par1> <split par2> <outfile> " << endl;
  // Open input surface file
  ifstream is(argv[1]);
  ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

  ObjectHeader head;
  is >> head;

  // Read volume from file
  SplineVolume vol;
  is >> vol;

  double hull_thick = atof(argv[2]);

  // Number of stiffeners
  int nmb_stiff = atoi(argv[3]);
  int nmb_stiff2 = atoi(argv[4]);

  // Hight of stiffener
  double hight_frac = atof(argv[5]);
  if (hight_frac > 1.0)
    hight_frac = 1.0;
  
  // Thickness of stiffener
  double len = atof(argv[6]);
  double thickness = atof(argv[7]);
  double thick = thickness/len;

  vector<double> split_pars(2);
  split_pars[0] = atof(argv[8]);
  split_pars[1] = atof(argv[9]);

  // Open outfile
  ofstream os(argv[10]);
  ALWAYS_ERROR_IF(os.bad(), "Bad or no output filename");

  // Parameter domain
  double u1 = vol.startparam(0);
  double u2 = vol.endparam(0);
  // double v1 = vol.startparam(1);
  double v2 = vol.endparam(1);
  double w1 = vol.startparam(2);
  double w2 = vol.endparam(2);

  int ki, kr, kc;
  double eps = 1.0e-10;  // More or less dummy
  double t1, t2;
  double span1 = vol.knotSpan(0, vol.order(0)-1);
  t1 = u1 + (1.0 - hight_frac)*span1;
  t2 = u2 - (1.0 - hight_frac)*span1;

  // Make deck
  double param[12];
  param[6] = param[0] = t1; 
  param[10] = param[7] = param[4] = param[1] = v2;
  param[5] = param[2] = w1;
  param[9] = param[3] = t2; 
  param[11] = param[8] = w2;

  double coef[24];
  Point tmp_pt;
  for (ki=0; ki<4; ki++)
    {
      vol.point(tmp_pt, param[ki*3], param[ki*3+1], param[ki*3+2]);
      for (kr=0; kr<3; kr++)
	coef[3*ki+kr] = tmp_pt[kr];
    }
  for (ki=0; ki<4; ki++)
    {
      coef[12+3*ki] = coef[3*ki];
      coef[12+3*ki+1] = coef[3*ki+1];
      coef[12+3*ki+2] = coef[3*ki+2] + thickness;
    }
  double et[4], et2[4];
  et[0] = et[1] = 0.0;
  et[2] = et[3] = 1.0;
  et2[0] = et2[1] = w1;
  et2[2] = et2[3] = w2;
  shared_ptr<SplineVolume> deck = 
    shared_ptr<SplineVolume>(new SplineVolume(2, 2, 2, 2, 2, 2,
					      et, et2, et, coef, 3));
  deck->raiseOrder(2, 0, 0);

  // Find parameter values at which to split volume in length direction
  if (nmb_stiff2 < 2)
    {
      std::cout << "Not enough stiffeners" << std::endl;
      exit(-1);
    }

  // Make parameters for split hull volumes
  vector<double> vol_split(6+2*nmb_stiff2);
  vol_split[1] = t1;
  vol_split[2] = split_pars[0];
  vol_split[vol_split.size()-3] = split_pars[1];
  vol_split[vol_split.size()-2] = t2;

  shared_ptr<SplineSurface> inner_hull = 
    shared_ptr<SplineSurface>(vol.constParamSurface(v2, 1));
  shared_ptr<SplineCurve> inner_cv =
    shared_ptr<SplineCurve>(inner_hull->constParamCurve(w1, true));
  Point pt1(coef+12, coef+15);
  Point pt2(coef+15, coef+18);
  double len_across = pt1.dist(pt2);

  Point clo_pt1, clo_pt2;
  double clo_par1, clo_par2;
  double clo_d1, clo_d2;
  inner_cv->closestPoint(pt1, inner_cv->startparam(), inner_cv->endparam(),
			 clo_par1, clo_pt1, clo_d1, &t1);
  inner_cv->closestPoint(pt2, inner_cv->startparam(), inner_cv->endparam(),
			 clo_par2, clo_pt2, clo_d2, &t2);
  vol_split[0] = clo_par1;
  vol_split[vol_split.size()-1] = clo_par2;

  // Create stiffeners in the length direction
  double thick2 = thickness/len_across;
  double del_stiff = 1.0/(double)(nmb_stiff2+1);
  Point pt3, pt4;
  vector<shared_ptr<SplineVolume> > stiffeners2;
  vector<double> deck_split;
  for (ki=0; ki<nmb_stiff2; ki++)
    {
      double ta = (ki+1)*del_stiff - 0.5*thick2;
      double tb = (ki+1)*del_stiff + 0.5*thick2;
      deck_split.push_back(ta);
      deck_split.push_back(tb);

      deck->point(pt1, ta, w1, deck->startparam(2));
      deck->point(pt2, tb, w1, deck->startparam(2));
      deck->point(pt3, ta, w2, deck->startparam(2));
      deck->point(pt4, tb, w2, deck->startparam(2));
      for (kr=0; kr<3; kr++)
	{
	  coef[kr] = coef[6+kr] = pt1[kr];
	  coef[3+kr] = coef[9+kr] = pt2[kr];
	  coef[12+kr] = coef[18+kr] = pt3[kr];
	  coef[15+kr] = coef[21+kr] = pt4[kr];
	}
      coef[2] = coef[5] = coef[14] = coef[17] = hull_thick;
		  
      pt1 = Point(coef, coef+3);
      pt2 = Point(coef+3, coef+6);
      inner_cv->closestPoint(pt1, inner_cv->startparam(), inner_cv->endparam(),
			     clo_par1, clo_pt1, clo_d1);
      inner_cv->closestPoint(pt2, inner_cv->startparam(), inner_cv->endparam(),
			     clo_par2, clo_pt2, clo_d2);

      vol_split[3+2*ki] = clo_par1;
      vol_split[3+2*ki+1] = clo_par2;

      double et3[4];
      et3[0] = et3[1] = ta;
      et3[2] = et3[3] = tb;
      shared_ptr<SplineVolume> stiff2 = 
	shared_ptr<SplineVolume>(new SplineVolume(2, 2, 2, 2, 2, 2,
						  et3, et, et2, coef, 3));
      stiff2->raiseOrder(2, 0, 0);

      stiffeners2.push_back(stiff2);
    }

  // Split hull volumes
  vector<shared_ptr<SplineVolume> > sub_vol = vol.split(vol_split, 0);

  // Split deck
  vector<shared_ptr<SplineVolume> > sub_deck = deck->split(deck_split, 0);

  ofstream os2("bd_tmp.g2");

//   for (ki=0; ki<(int)sub_vol.size(); ki++)
//     {
//       sub_vol[ki]->writeStandardHeader(os);
//       sub_vol[ki]->write(os);
//       vector<shared_ptr<SplineSurface> > bd_sfs = 
// 	sub_vol[ki]->getBoundarySurfaces();
//       for (size_t kj=0; kj<bd_sfs.size(); ++kj)
// 	{
// 	  bd_sfs[kj]->writeStandardHeader(os2);
// 	  bd_sfs[kj]->write(os2);
// 	}
//     }

//   for (ki=0; ki<(int)sub_deck.size(); ki++)
//     {
//       sub_deck[ki]->writeStandardHeader(os);
//       sub_deck[ki]->write(os);
//       vector<shared_ptr<SplineSurface> > bd_sfs = 
// 	sub_deck[ki]->getBoundarySurfaces();
//       for (size_t kj=0; kj<bd_sfs.size(); ++kj)
// 	{
// 	  bd_sfs[kj]->writeStandardHeader(os2);
// 	  bd_sfs[kj]->write(os2);
// 	}
//     }


//    for (ki=0; ki<(int)stiffeners2.size(); ki++)
//     {
//       stiffeners2[ki]->writeStandardHeader(os);
//       stiffeners2[ki]->write(os);
//       vector<shared_ptr<SplineSurface> > bd_sfs = 
// 	stiffeners2[ki]->getBoundarySurfaces();
//       for (size_t kj=0; kj<bd_sfs.size(); ++kj)
// 	{
// 	  bd_sfs[kj]->writeStandardHeader(os2);
// 	  bd_sfs[kj]->write(os2);
// 	}
//     }


   // Find parameter values at which to split volume in across the hull
  if (nmb_stiff < 2)
    {
      std::cout << "Not enough stiffeners" << std::endl;
      exit(-1);
    }

  // Collection of all sub volumes
  vector<shared_ptr<SplineVolume> > volumes;
  vector<shared_ptr<SplineVolume> > tmp_volumes;

  // Reorganize current volumes in vertical and horizontal direction
  // Top volumes
  vector<shared_ptr<SplineVolume> > above_deck(4);
  above_deck[0] = sub_vol[0];
  above_deck[1] = sub_vol[1];
  above_deck[2] = sub_vol[sub_vol.size()-2];
  above_deck[3] = sub_vol[sub_vol.size()-1];
  for (ki=0; ki<4; ++ki)
    above_deck[ki]->raiseOrder(0, 2, 0);

  // Vertical
  vector<shared_ptr<SplineVolume> > vertical(2+stiffeners2.size());
  vertical[0] = sub_vol[2];
  vertical[0]->reverseParameterDirection(0);
  vertical[0]->swapParameterDirection(0, 1);
  for (ki=0; ki<(int)stiffeners2.size(); ki++)
    vertical[1+ki] = stiffeners2[ki];
  vertical[1+ki] = sub_vol[sub_vol.size()-3];
  vertical[1+ki]->swapParameterDirection(0, 1);
  vertical[1+ki]->reverseParameterDirection(0);

  // Make sure that the vertical volumes have the same spline space
  vector<shared_ptr<SplineCurve> > vertical_cv;
  for (kr=0; kr<(int)vertical.size(); kr++)
    {
      vertical[kr]->setParameterDomain(0.0, 1.0, 0.0, 1.0, 0.0, 1.0); // To avoid trouble
      vertical_cv.push_back(VolumeTools::representVolumeAsCurve(*vertical[kr].get(), 1));
    }
  unifyCurveSplineSpace(vertical_cv, eps);
  for (kr=0; kr<(int)vertical.size(); kr++)
    vertical[kr] = VolumeTools::representCurveAsVolume(*vertical_cv[kr].get(), 1,
					  vertical[kr]->basis(0),
					  vertical[kr]->basis(2),
					  vertical[kr]->rational());

  // Horizontal, from left to right
  vector<shared_ptr<SplineVolume> > horizontal(2*(stiffeners2.size()+1));
  for (ki=0, kr=0; ki<(int)sub_deck.size(); ki+=2)
    {
      horizontal[kr++] = sub_vol[3+ki];
      horizontal[kr] = sub_deck[ki];
      horizontal[kr++]->swapParameterDirection(1, 2);
    }

  // Make sure that two and two horizonal volumes have the same spline space
  for (ki=0; ki<(int)horizontal.size(); ki+=2)
    {
      vector<shared_ptr<SplineCurve> > horizontal_cv;
      for (kr=0; kr<2; kr++)
	{
	  horizontal[ki+kr]->setParameterDomain(0.0, 1.0, 0.0, 1.0, 0.0, 1.0); // To avoid trouble
	  horizontal_cv.push_back(VolumeTools::representVolumeAsCurve(*horizontal[ki+kr].get(), 
							 0));
	}
      unifyCurveSplineSpace(horizontal_cv, eps);
      for (kr=0; kr<2; kr++)
	horizontal[ki+kr] = VolumeTools::representCurveAsVolume(*horizontal_cv[kr].get(), 
						   0,
						   horizontal[ki+kr]->basis(1),
						   horizontal[ki+kr]->basis(2),
						   horizontal[ki+kr]->rational());
    }     

  // Left over pieces
  vector<shared_ptr<SplineVolume> > pieces1;
  vector<shared_ptr<SplineVolume> > pieces2;
  for (ki=1; ki<(int)sub_deck.size(); ki+=2)
    {
      pieces1.push_back(sub_deck[ki]);
      pieces1[pieces1.size()-1]->swapParameterDirection(1, 2);
    }
  for (ki=4; ki<(int)sub_vol.size()-4; ki+=2)
    pieces2.push_back(sub_vol[ki]);

  vector<double> parvals(2*nmb_stiff-2, 0.0);
  double wdel = (w2 - thick - w1 - thick)/(double)(nmb_stiff - 1);
  double wcurr;

  parvals[0] = wcurr = w1 + thick;
  for (ki=1; ki<nmb_stiff-1; ki++)
    {
      wcurr += wdel;
      parvals[2*ki-1] = wcurr - 0.5*thick;
      parvals[2*ki] = wcurr + 0.5*thick;
    }
  parvals[2*ki-1] = w2 - thick;

  // Vertical
  // Split first wall
  vector<shared_ptr<SplineVolume> > sub_vert1 = 
    vertical[0]->split(parvals,2);
  volumes.insert(volumes.end(),sub_vert1.begin(), sub_vert1.end());


  for (ki=1; ki<(int)vertical.size(); ki++)
    {
      // Split current wall
      vector<shared_ptr<SplineVolume> > sub_vert2 = 
	vertical[ki]->split(parvals,2);
      volumes.insert(volumes.end(),sub_vert2.begin(), sub_vert2.end());
      if (ki == 1)
	tmp_volumes.push_back(sub_vert2[0]);

      // Split bottom
      vector<shared_ptr<SplineVolume> > bottom = 
	horizontal[2*(ki-1)]->split(parvals,2);
      volumes.insert(volumes.end(), bottom.begin(), bottom.end());
      if (ki <= 2)
	tmp_volumes.push_back(bottom[0]);

      // Split deck
      vector<shared_ptr<SplineVolume> > top = 
	horizontal[2*(ki-1)+1]->split(parvals,2);
      volumes.insert(volumes.end(), top.begin(), top.end());

      // Fetch boundary curves for coons patch
      for (kc=0; kc < (int)sub_vert1.size(); kc+=2)
	{
	  vector<shared_ptr<SplineCurve> > sub_cvs1;
	  vector<shared_ptr<SplineCurve> > sub_cvs2;

	  // Left surface
	  shared_ptr<SplineSurface> sf = 
	    shared_ptr<SplineSurface>(sub_vert1[kc]->getBoundarySurface(1));
	  shared_ptr<SplineCurve> cv1 =
	    shared_ptr<SplineCurve>(sf->constParamCurve(sf->startparam_v(),
							true));
	  shared_ptr<SplineCurve> cv2 =
	    shared_ptr<SplineCurve>(sf->constParamCurve(sf->endparam_v(),
							true));
	  cv1->reverseParameterDirection();
	  cv2->reverseParameterDirection();
	  sub_cvs1.push_back(cv1);
	  sub_cvs2.push_back(cv2);

	  ofstream os3("bdsf_tmp.g2");
	  sf->writeStandardHeader(os3);
	  sf->write(os3);

	  // Bottom surface
	  sf = shared_ptr<SplineSurface>(bottom[kc]->getBoundarySurface(3));
	  cv1 = shared_ptr<SplineCurve>(sf->constParamCurve(sf->startparam_v(),
							    true));
	  cv2 = shared_ptr<SplineCurve>(sf->constParamCurve(sf->endparam_v(),
							    true));
	  sub_cvs1.push_back(cv1);
	  sub_cvs2.push_back(cv2);

	  sf->writeStandardHeader(os3);
	  sf->write(os3);

	  // Right surface
	  sf = shared_ptr<SplineSurface>(sub_vert2[kc]->getBoundarySurface(0));
	  cv1 = shared_ptr<SplineCurve>(sf->constParamCurve(sf->startparam_v(),
							    true));
	  cv2 = shared_ptr<SplineCurve>(sf->constParamCurve(sf->endparam_v(),
							    true));
	  sub_cvs1.push_back(cv1);
	  sub_cvs2.push_back(cv2);

	  sf->writeStandardHeader(os3);
	  sf->write(os3);

	  // Top surface
	  sf = shared_ptr<SplineSurface>(top[kc]->getBoundarySurface(2));
	  cv1 = shared_ptr<SplineCurve>(sf->constParamCurve(sf->startparam_v(),
							    true));
	  cv2 = shared_ptr<SplineCurve>(sf->constParamCurve(sf->endparam_v(),
							    true));
	  cv1->reverseParameterDirection();
	  cv2->reverseParameterDirection();
	  sub_cvs1.push_back(cv1);
	  sub_cvs2.push_back(cv2);
	  
	  sf->writeStandardHeader(os3);
	  sf->write(os3);

	  // Make loop
	  vector<shared_ptr<ParamCurve> > tmp_loop1;
	  vector<shared_ptr<ParamCurve> > tmp_loop2;
	  for (size_t kj=0; kj<sub_cvs1.size(); ++kj)
	    tmp_loop1.push_back(sub_cvs1[kj]);
	  for (size_t kj=0; kj<sub_cvs2.size(); ++kj)
	    tmp_loop2.push_back(sub_cvs2[kj]);

	  CurveLoop loop1(tmp_loop1, eps);
	  CurveLoop loop2(tmp_loop2, eps);

	  // Make Coons patches
	  SplineSurface *coons1 = CoonsPatchGen::createCoonsPatch(loop1);
	  SplineSurface *coons2 = CoonsPatchGen::createCoonsPatch(loop2);

	  // Make stiffener volume
	  vector<shared_ptr<SplineSurface> > stiff_sfs(2);
	  stiff_sfs[0] = shared_ptr<SplineSurface>(coons1);
	  stiff_sfs[1] = shared_ptr<SplineSurface>(coons2);

	  shared_ptr<SplineVolume> stiff = 
	    shared_ptr<SplineVolume>(LoftVolumeCreator::loftVolume(stiff_sfs.begin(),
								   2));
	  volumes.push_back(stiff);
	}
      sub_vert1 = sub_vert2;
    }

  // Split part above deck
  vector<shared_ptr<SplineVolume> > tmp_sub;
  for (ki=0; ki<4; ki++)
    {
      tmp_sub = above_deck[ki]->split(parvals,2);
      volumes.insert(volumes.end(), tmp_sub.begin(), tmp_sub.end());
    }

  // Split remaining pieces
  for (ki=0; ki<(int)pieces1.size(); ki++)
    {
      tmp_sub = pieces1[ki]->split(parvals,2);
      volumes.insert(volumes.end(), tmp_sub.begin(), tmp_sub.end());
    }
  for (ki=0; ki<(int)pieces2.size(); ki++)
    {
      tmp_sub = pieces2[ki]->split(parvals,2);
      volumes.insert(volumes.end(), tmp_sub.begin(), tmp_sub.end());
      if (ki == 0)
	tmp_volumes.push_back(tmp_sub[0]);
    }

  // Output
  for (ki=0; ki<(int)volumes.size(); ki++)
    {
      volumes[ki]->writeStandardHeader(os);
      volumes[ki]->write(os);
      vector<shared_ptr<SplineSurface> > bd_sfs = 
	volumes[ki]->getBoundarySurfaces();
      for (size_t kj=0; kj<bd_sfs.size(); ++kj)
	{
	  bd_sfs[kj]->writeStandardHeader(os2);
	  bd_sfs[kj]->write(os2);
	}
    }

  ofstream os3("tmp_vol.g2");
  ofstream os4("tmp_bdsf.g2");
  for (ki=0; ki<(int)tmp_volumes.size(); ki++)
    {
      tmp_volumes[ki]->writeStandardHeader(os3);
      tmp_volumes[ki]->write(os3);
      vector<shared_ptr<SplineSurface> > bd_sfs = 
	tmp_volumes[ki]->getBoundarySurfaces();
      for (size_t kj=0; kj<bd_sfs.size(); ++kj)
	{
	  bd_sfs[kj]->writeStandardHeader(os4);
	  bd_sfs[kj]->write(os4);
	}
    }


}

