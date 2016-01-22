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

#include "GoTools/creators/HermiteAppS.h"
#include "GoTools/creators/IntCrvEvaluator.h"
#include "GoTools/geometry/GapRemoval.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/SplineUtils.h"
#include "GoTools/creators/AdaptCurve.h"
#include "GoTools/creators/SpaceIntCrv.h"
#include "GoTools/creators/TrimCurve.h"
#include "GoTools/creators/ApproxSurf.h"
#include "GoTools/creators/SmoothCurveSet.h"
#include "GoTools/creators/SmoothSurf.h"
#include <fstream>

using std::vector;
using std::setprecision;
using std::endl;
using std::pair;
using std::make_pair;

namespace Go
{

//===========================================================================
void
GapRemoval::removeGapSpline(shared_ptr<SplineSurface>& srf1, 
			    shared_ptr<CurveOnSurface>& bd_cv1,
			    double start1, double end1,
			    shared_ptr<SplineSurface>& srf2, 
			    shared_ptr<CurveOnSurface>& bd_cv2,
			    double start2, double end2, 
			    Point vertex1, Point vertex2,
			    double epsge, bool *same_orientation)
//===========================================================================
{
    // T-connection configurations between adjacent surfaces are currently
    // not handled in all cases. Check that the edge represent an entire boundary for
    // both adjacent surfaces
  Point f1_p1 = bd_cv1->faceParameter(start1);
  Point f1_p2 = bd_cv1->faceParameter(end1);
  Point f2_p1 = bd_cv2->faceParameter(start2);
  Point f2_p2 = bd_cv2->faceParameter(end2);

  RectDomain dom1 = srf1->parameterDomain();
  RectDomain dom2 = srf2->parameterDomain();

  double ptol = 1.0e-10;
  int bd1, bd2;  // Specifies the surface boundaries corresponding to 
  // the current edges
  // 0 = umin, 1 = umax, 2 = vmin,  3 = vmax
  bool same_orient1, same_orient2;
  bd1 = bd_cv1->whichBoundary(epsge, same_orient1);
  bd2 = bd_cv2->whichBoundary(epsge, same_orient2);

  if (bd1 < 0 || bd2 < 0)
    return;  // Unexpected situation

  // @@@ VSK, We should have a special treatment of the surface corners to
  // avoid creating new gaps towards other surfaces meeting in the corner
  // That is not implemented yet.
  shared_ptr<SplineSurface> s1 = srf1;
  shared_ptr<SplineSurface> s2 = srf2;

  bool atcorner1 = (fabs(bd_cv1->startparam() - start1) < ptol &&
    fabs(bd_cv1->endparam() - end1) < ptol);
  bool atcorner2 = (fabs(bd_cv2->startparam() - start2) < ptol &&
		  fabs(bd_cv2->endparam() - end2) < ptol);
  bool keep_first = false, keep_second = false;
  if (!atcorner1 && !atcorner2)
    return;  // Specific T-situation. Currently not handled
  else if (!atcorner1)
    {
      // Pick the relevant part of surface one, modify surface two
      double u1 = (bd1 > 1) ? start1 : dom1.umin();
      double v1 = (bd1 <= 1) ? start1 : dom1.vmin();
      double u2 = (bd1 > 1) ? end1 : dom1.umax();
      double v2 = (bd1 <= 1) ? end1 : dom1.vmax();
      shared_ptr<SplineSurface> srf3(s1->subSurface(u1, v1, u2, v2, ptol));
      s1 = srf3;
      dom1 = s1->containingDomain();
      keep_first = true;
    }
  else if (!atcorner2)
    {
      // Pick the relevant part of surface two, modify surface one
      double u1 = (bd2 > 1) ? start2 : dom2.umin();
      double v1 = (bd2 <= 1) ? start2 : dom2.vmin();
      double u2 = (bd2 > 1) ? end2 : dom2.umax();
      double v2 = (bd2 <= 1) ? end2 : dom2.vmax();
      shared_ptr<SplineSurface> srf3(s2->subSurface(u1, v1, u2, v2, ptol));
      s2 = srf3;
      dom2 = s2->containingDomain();
      keep_second = true;
    }

  bool opposite = false;
  double t1 = (bd1 == 0 ||  bd1 == 1) ? f1_p1[1] : f1_p1[0];
  double t2 = (bd1 == 0 ||  bd1 == 1) ? f1_p2[1] : f1_p2[0];
  double t3 = (bd2 == 0 ||  bd2 == 1) ? f2_p1[1] : f2_p1[0];
  double t4 = (bd2 == 0 ||  bd2 == 1) ? f2_p2[1] : f2_p2[0];
  if ((t2 - t1)*(t4 -t3) < 0.0)
    opposite = true;
  if ((same_orient1 && !same_orient2) || (!same_orient1 && same_orient2))
    opposite = !opposite;
  if (same_orientation != NULL && !(*same_orientation))
    opposite = !opposite;

  try {
    GeometryTools::averageBoundaryCoefs(s1, bd1, keep_first, s2, bd2, keep_second, true, 
			 vertex1, true, vertex2, opposite);
  }
  catch(...)
    {
      return;
    }

  // Update boundary curves
  bool updated1, updated2;
  updated1 = bd_cv1->updateIsoCurves();
  updated2 = bd_cv2->updateIsoCurves();

  double mdist1, mdist2;
  int nmb_sample = 200;
  checkBoundaryDist(bd_cv1, bd_cv2, start1, end1, start2, end2,
		    nmb_sample, mdist1, mdist2);
  if (getenv("DEBUG") && (*getenv("DEBUG")) == '1')
    {
  std::cout << "removeGapSpline, distances: " << mdist1 << ", ";
  std::cout << mdist2 << std::endl;
    }
}

//===========================================================================
double
GapRemoval::removeGapTrim(shared_ptr<CurveOnSurface>& bd_cv1,
			  double start1, double end1,
			  shared_ptr<CurveOnSurface>& bd_cv2,
			  double start2, double end2, Point vertex1, Point vertex2,
			  double epsge)
//===========================================================================
{
  // It is two possibilities for a gap between two trimmed surfaces,
  // the surfaces don't lie close enough or the trimming curves don't lie
  // close enough. This function attempts to close the gap by 
  // regenerating the trimming curves
  double max_dist = 0.0;

  // Check if any of the surfaces are boundary trimmed
  bool same_orient1, same_orient2;
  int bd1 = bd_cv1->whichBoundary(epsge, same_orient1);
  int bd2 = bd_cv2->whichBoundary(epsge, same_orient2);
  if (bd1 != -1 && bd2 != -1)
    return 1.0e8;   // Not an appropriate method. Return a large numbeer

  // Check if the boundary curves have a corresponding orientation
  Point p11 = bd_cv1->ParamCurve::point(start1);
  Point p12 = bd_cv1->ParamCurve::point(end1);
  Point p21 = bd_cv2->ParamCurve::point(start2);
  Point p22 = bd_cv2->ParamCurve::point(end2);
  bool same_orientation = 
    (p11.dist(p21) + p12.dist(p22) < p11.dist(p22) + p12.dist(p21));
  int keep_crv = (bd1 >= 0) ? 1 : ((bd2 >= 0) ? 2 : 0);
  shared_ptr<IntCrvEvaluator> int_crv = 
    shared_ptr<IntCrvEvaluator>(new IntCrvEvaluator(bd_cv1, start1, end1,
						    bd_cv2, start2, end2,
						    same_orientation,
						    keep_crv));
  
  vector<int> dims(3);
  dims[0] = bd_cv1->dimension();
  dims[1] = dims[2] = 2;  // Curves in the parameter domain

  // Uses the parameterization of the 1. trimming curve
  vector<double> initpars;
  shared_ptr<SplineCurve> gcrv = 
    dynamic_pointer_cast<SplineCurve, ParamCurve>(bd_cv1->spaceCurve());
  if (gcrv.get())
    gcrv->basis().knotsSimple(initpars);
  size_t ki;
  for (ki=0; ki<initpars.size();)
    {
      if (initpars[ki] < start1 || initpars[ki] > end1)
	initpars.erase(initpars.begin() + ki);
      else
	ki++;
    }

  // Define approximator
  shared_ptr<HermiteAppS> approximator;
  if (initpars.size() > 2)
    approximator = 
      shared_ptr<HermiteAppS>(new HermiteAppS(int_crv.get(), &initpars[0], 
					      (int)initpars.size(),
					      epsge, epsge, dims));
  else
    approximator = 
      shared_ptr<HermiteAppS>(new HermiteAppS(int_crv.get(), epsge, epsge, dims));

  // // Check initial points
  // double init_dist = int_crv->getMaxErr();

  try {
    // Approximate
    approximator->refineApproximation();
    max_dist = int_crv->getMaxErr();
  }
  catch (...)
    {
      max_dist = int_crv->getMaxErr();
//       std::ofstream out("distant_sfs.g2");
//       srf1->writeStandardHeader(out);
//       srf1->write(out);
//       srf2->writeStandardHeader(out);
//       srf2->write(out);

    }

  if (max_dist < epsge)
    {
      // Fetch curves
      vector<shared_ptr<SplineCurve> > crvs = approximator->getCurves();

      // Check if the new curves are large enough
      Point p13 = bd_cv1->ParamCurve::point(bd_cv1->startparam());
      Point p14 = bd_cv1->ParamCurve::point(bd_cv1->endparam());
      Point p23 = bd_cv2->ParamCurve::point(bd_cv2->startparam());
      Point p24 = bd_cv2->ParamCurve::point(bd_cv2->endparam());
      
      shared_ptr<SplineCurve> space_crv1 = crvs[0];
      shared_ptr<SplineCurve> space_crv2 = 
	shared_ptr<SplineCurve>(crvs[0]->clone());
      shared_ptr<SplineCurve> param_crv1 = crvs[1];
      shared_ptr<SplineCurve> param_crv2 = crvs[2];

      // Ensure that the parameter interval of the new curves related
      // to the second boundary curve corresponds to the boundary curve
      space_crv2->setParameterInterval(start2, end2);
      param_crv2->setParameterInterval(start2, end2);

      if (p11.dist(p13) > epsge)
	{
	  // Boundary curve 1 must be extended at the start. Fetch the
	  // appropriate piece from the initial curve
	  shared_ptr<CurveOnSurface> sub1 =
	    shared_ptr<CurveOnSurface>(bd_cv1->subCurve(bd_cv1->startparam(),
							start1));
	  shared_ptr<SplineCurve> space = 
	    shared_ptr<SplineCurve>(sub1->spaceCurve()->geometryCurve());
	  shared_ptr<SplineCurve> param = 
	    shared_ptr<SplineCurve>(sub1->parameterCurve()->geometryCurve());
	  double dist1, dist2;
	  space->appendCurve(space_crv1.get(), 1, dist1, false);
	  param->appendCurve(param_crv1.get(), 1, dist2, false);
	  space_crv1 = space;
	  param_crv1 = param;
	}

      if (p12.dist(p14) > epsge)
	{
	  // Boundary curve 1 must be extended at the end. Fetch the
	  // appropriate piece from the initial curve
	  shared_ptr<CurveOnSurface> sub1 =
	    shared_ptr<CurveOnSurface>(bd_cv1->subCurve(end1,
							bd_cv1->endparam()));
	  shared_ptr<SplineCurve> space = 
	    shared_ptr<SplineCurve>(sub1->spaceCurve()->geometryCurve());
	  shared_ptr<SplineCurve> param = 
	    shared_ptr<SplineCurve>(sub1->parameterCurve()->geometryCurve());
	  double dist1, dist2;
	  space_crv1->appendCurve(space.get(), 1, dist1, false);
	  param_crv1->appendCurve(param.get(), 1, dist2, false);
	}

      if (p21.dist(p23) > epsge)
	{
	  // Boundary curve 2 must be extended at the start. Fetch the
	  // appropriate piece from the initial curve
	  shared_ptr<CurveOnSurface> sub1 =
	    shared_ptr<CurveOnSurface>(bd_cv2->subCurve(bd_cv2->startparam(),
							start2));
	  shared_ptr<SplineCurve> space = 
	    shared_ptr<SplineCurve>(sub1->spaceCurve()->geometryCurve());
	  shared_ptr<SplineCurve> param = 
	    shared_ptr<SplineCurve>(sub1->parameterCurve()->geometryCurve());
	  double dist1, dist2;
	  space->appendCurve(space_crv2.get(), 1, dist1, false);
	  param->appendCurve(param_crv2.get(), 1, dist2, false);
	  space_crv2 = space;
	  param_crv2 = param;
	}

      if (p22.dist(p24) > epsge)
	{
	  // Boundary curve 2 must be extended at the end. Fetch the
	  // appropriate piece from the initial curve
	  shared_ptr<CurveOnSurface> sub1 =
	    shared_ptr<CurveOnSurface>(bd_cv2->subCurve(end2,
							bd_cv2->endparam()));
	  shared_ptr<SplineCurve> space = 
	    shared_ptr<SplineCurve>(sub1->spaceCurve()->geometryCurve());
	  shared_ptr<SplineCurve> param = 
	    shared_ptr<SplineCurve>(sub1->parameterCurve()->geometryCurve());
	  double dist1, dist2;
	  space_crv2->appendCurve(space.get(), 1, dist1, false);
	  param_crv2->appendCurve(param.get(), 1, dist2, false);
	}

       // Update trimming curves
      if (keep_crv != 1)
	{
	  bd_cv1->setCurves(space_crv1, param_crv1);
	  bd_cv1->setCurveTypeInfo(2);
	}
      if (keep_crv != 2)
	{
	  bd_cv2->setCurves(space_crv2, param_crv2);
	  bd_cv2->setCurveTypeInfo(2);
	}
    }

  double mdist1, mdist2;
  int nmb_sample = 200;
  checkBoundaryDist(bd_cv1, bd_cv2, start1, end1, start2, end2,
		    nmb_sample, mdist1, mdist2);
  if (getenv("DEBUG") && (*getenv("DEBUG")) == '1')
    {
  std::cout << "removeGapTrim, distances: " << mdist1 << ", ";
  std::cout << mdist2 << ", computed maxdist: " << max_dist << std::endl;
    }

  return max_dist;
}

//===========================================================================
bool
GapRemoval::removeGapSplineTrim(shared_ptr<SplineSurface>& srf1, 
				vector<shared_ptr<CurveOnSurface> >& bd_cv1,
				vector<double> start1, vector<double> end1,
				vector<shared_ptr<CurveOnSurface> >& bd_cv2,
				vector<double> start2, vector<double> end2, 
				Point vertex1, Point vertex2, double epsge)
//===========================================================================
{
  // Adapt a spline surface to a number of surfaces along associated
  // trimming curves

  if (bd_cv1.size() == 0)
    return false;  

  // Check input
  ASSERT(bd_cv1.size() == bd_cv2.size());
  size_t ki;
  for (ki=0; ki<bd_cv1.size(); ++ki)
    if (bd_cv1[ki]->underlyingSurface().get() != srf1.get())
      THROW("Iconsistent data structure");

  // Update trimming curves corresponding to the other surfaces
  for (ki=0; ki<bd_cv2.size(); ++ki)
    if (bd_cv2[ki].get())
      bd_cv2[ki]->updateCurves(0.5*epsge);

  // Fetch boundary curve of spline surface
  // double ptol = 1.0e-10;
  int bd1;
  bool same_orient;
  bd1 = bd_cv1[0]->whichBoundary(epsge, same_orient);
  int idx[] = {3,1,0,2};
  int ccw_nmb = idx[bd1];
  shared_ptr<SplineCurve> crv = 
    shared_ptr<SplineCurve>(srf1->edgeCurve(ccw_nmb));

  // Check continuity at the sub boundaries
  // First get curve parameter corrsponding to the boundary curve
  size_t nmb_bd = bd_cv1.size();
  Point face_par1 = bd_cv1[0]->faceParameter(start1[0]);
  Point face_par2 = bd_cv1[nmb_bd-1]->faceParameter(end1[nmb_bd-1]);
  int dir = (bd1 == 0 || bd1 == 1) ? 1 : 0;
  double par1 = face_par1[dir];
  double par2 = face_par2[dir];
  (void)crv->basis().knotIntervalFuzzy(par1);
  (void)crv->basis().knotIntervalFuzzy(par2);
  int cont1 = crv->order() - crv->basis().knotMultiplicity(par1) - 1;
  int cont2 = crv->order() - crv->basis().knotMultiplicity(par2) - 1;

  // Compute boundary conditions to the updated boundary curve piece
  vector<Point> der1(cont1+1);
  vector<Point> der2(cont2+1);
  if (cont1 >= 0)
    crv->point(der1, par1, cont1);
  else if (par1 == start1[0])
    der1.push_back(vertex1);
  else if (par1 == end1[nmb_bd-1])
    der1.push_back(vertex2);
  if (cont2 >= 0)
    crv->point(der2, par2, cont2);
  else if (par2 == start2[0])
    der2.push_back(vertex1);
  else if (par2 == end2[nmb_bd-1])
    der2.push_back(vertex2);

  // Make a curve piece approximating the trimming curve
  bool is_turned = false;
  if (par1 > par2)
    {
      std::swap(par1,par2);
      std::swap(cont1, cont2);
      is_turned = true;
    }

  shared_ptr<SplineCurve> sub_crv = 
    shared_ptr<SplineCurve>(crv->subCurve(par1, par2));

  // Make evaluator based curve
  vector<bool> opposite(nmb_bd, false);
  for (size_t idx=0; idx<nmb_bd; ++idx)
    {
      if (bd_cv2[idx].get() == 0)
	  continue;

      // bool same1;
      // int bd_1 = bd_cv1[idx]->whichBoundary(epsge, same1);

      Point pt1 = bd_cv1[idx]->ParamCurve::point(start1[idx]);
      Point pt2 = bd_cv1[idx]->ParamCurve::point(end1[idx]);
      Point pt3 = bd_cv2[idx]->ParamCurve::point(start2[idx]);
      Point pt4 = bd_cv2[idx]->ParamCurve::point(end2[idx]);
      if (pt1.dist(pt3) + pt2.dist(pt4) > pt1.dist(pt4) + pt2.dist(pt3))
	opposite[idx] = true;
    }

  shared_ptr<ParamCurve> boundary = sub_crv;
  shared_ptr<SpaceIntCrv> int_crv = 
    shared_ptr<SpaceIntCrv>(new SpaceIntCrv(boundary, (bd1 <= 1) ? 1 : 0,
					    bd_cv1, start1, end1, bd_cv2, 
					    start2, end2, opposite, same_orient));

  // Define approximator
  shared_ptr<AdaptCurve> adapt = 
    shared_ptr<AdaptCurve>(new AdaptCurve(int_crv.get(), 0.5*epsge, sub_crv));

  // Set boundary conditions
  adapt->setEndPoints(der1, der2);

  // Approximate
  adapt->approximate();

  // Fetch result
  double maxdist, avdist;
  sub_crv = adapt->getAdaptCurve(maxdist, avdist);
  if (maxdist > epsge)
    {
      MESSAGE("Approximation accuracy not reached");
      return false;
    }

  if (getenv("DEBUG") && (*getenv("DEBUG")) == '1')
    {
  std::ofstream out1("cv_mod.g2");
  sub_crv->writeStandardHeader(out1);
  sub_crv->write(out1);
    }

  // Create a new curve where the current boundary curve piece is replaced
  // with the updated curve
  shared_ptr<SplineCurve> new_bd = replaceCurvePiece(crv, sub_crv, par1, cont1, 
						     par2, cont2);

  if (getenv("DEBUG") && (*getenv("DEBUG")) == '1')
    {
  std::ofstream out2("cv_mod2.g2");
  new_bd->writeStandardHeader(out2);
  new_bd->write(out2);
    }

  // Replace surface coefficients
  bool replaced;
  replaced = srf1->replaceBoundaryCurve(bd1, new_bd);

  if (getenv("DEBUG") && (*getenv("DEBUG")) == '1')
    {
      std::ofstream out3("sf_mod.g2");
      srf1->writeStandardHeader(out3);
      srf1->write(out3);
    }

  // Update input boundary curve
  bool updated;
  double parval;
  if (bd1 == 0)
    parval = srf1->startparam_u();
  else if (bd1 == 1)
    parval = srf1->endparam_u();
  else if (bd1 == 2)
    parval = srf1->startparam_v();
  else 
    parval = srf1->endparam_v();

  std::ofstream out4("bd_mod2.g2");

  double mdist1, mdist2;
  int nmb_sample = 200;
  for (ki=0; ki<nmb_bd; ++ki)
    {
      updated = bd_cv1[ki]->updateIsoCurves((bd1<=1) ? 1 : 2, parval, bd1);

  if (getenv("DEBUG") && (*getenv("DEBUG")) == '1')
    {
      bd_cv1[ki]->spaceCurve()->writeStandardHeader(out4);
      bd_cv1[ki]->spaceCurve()->write(out4);
    }

      if (bd_cv2[ki].get())
	{
	  checkBoundaryDist(bd_cv1[ki], bd_cv2[ki], start1[ki], 
			    end1[ki], start2[ki], end2[ki],
			    nmb_sample, mdist1, mdist2);
  if (getenv("DEBUG") && (*getenv("DEBUG")) == '1')
    {
	  std::cout << "removeGapSplineTrim, distances: " << mdist1 << ", ";
	  std::cout << mdist2 << ", from approx: " << maxdist;
	  std::cout << ", " << avdist << std::endl << std::endl;
    }
	}
    }

  return true;
}

//===========================================================================
  void GapRemoval::modifySplineSf(shared_ptr<ParamSurface>& psrf1, 
				  vector<shared_ptr<CurveOnSurface> >& bd_cv1,
				  vector<double> start1, vector<double> end1,
				  shared_ptr<ParamSurface>& srf2, 
				  vector<shared_ptr<CurveOnSurface> >& bd_cv2,
				  vector<double> start2, vector<double> end2, 
				  double epsge)
//===========================================================================
  {
    // Estimate curve lengths to get a reasonable distribution of
    // sample points
    size_t idx;
    size_t nmb_bd = bd_cv1.size();
    vector<double> len(nmb_bd);
    double len_sum = 0.0;
    for (idx=0; idx<nmb_bd; ++idx)
      {
	len[idx] = bd_cv1[idx]->estimatedCurveLength();
	len_sum += len[idx];
      }
    
    vector<double> points;
    vector<double> parvals;
    int nmb_sample = 200;
    for (idx=0; idx<nmb_bd; ++idx)
      {
	// Define evaluator based curve corresponding to the unspecified surface
	shared_ptr<TrimCurve> trim_crv = 
	  shared_ptr<TrimCurve>(new TrimCurve(bd_cv2[idx].get(),start2[idx], 
					      end2[idx]));

	// Check if the boundary curves have a corresponding orientation
	Point p11 = bd_cv1[idx]->ParamCurve::point(bd_cv1[idx]->startparam());
	Point p12 = bd_cv1[idx]->ParamCurve::point(bd_cv1[idx]->endparam());
	Point p21 = bd_cv2[idx]->ParamCurve::point(bd_cv1[idx]->startparam());
	Point p22 = bd_cv2[idx]->ParamCurve::point(bd_cv1[idx]->endparam());
	bool same_orientation = 
	  (p11.dist(p21) + p12.dist(p22) < p11.dist(p22) + p12.dist(p21));

	// Evaluate the intersection curve in a dense set of points
	double t1 = trim_crv->start();
	double t2 = trim_crv->end();
	int curr_sample = (int)(nmb_sample*len[idx]/len_sum);
	double tint = (t2 - t1)/(int)(curr_sample-1);
	double tpar;
	int ki;
	double fac = (end1[idx] - start1[idx])/(t2 - t1);
	for (tpar=t1, ki=0; ki<curr_sample; ++ki, tpar+=tint)
	  {
	    vector<Point> pts = trim_crv->eval(tpar);
	    points.insert(points.end(), pts[0].begin(), pts[0].end());

	    // Get corrsponding paramter values in the spline surface
	    double guess = (same_orientation) ? start1[idx] + (tpar - t1)*fac :
	      end1[idx] - (tpar - t1)*fac;
	    double tpar2, clo_d;
	    Point clo_pt;
	    bd_cv1[idx]->closestPoint(pts[0], start1[idx], end1[idx], tpar2, 
				      clo_pt, clo_d, &guess);
	    Point param = bd_cv1[idx]->faceParameter(tpar2);
	    parvals.insert(parvals.end(), param.begin(), param.end());
	  }
      }

  // Modify spline surface 
    // Fetch spline surface
    vector<shared_ptr<CurveOnSurface> > all_bd1;
    shared_ptr<SplineSurface> srf1 = getSplineAndBd(psrf1, all_bd1);

    // Add sample points along the other boundaries
    vector<double> sample_pts1;
    vector<double> sample_pars1;
    vector<int> bd_idx1;
    getBoundarySamples(all_bd1, bd_cv1, sample_pts1, sample_pars1,
		       bd_idx1, epsge);

  // Make approximator
    // Combine point conditions
    vector<double> all_pts1;
    vector<double> all_pars1;
    all_pts1.insert(all_pts1.end(), points.begin(), points.end());
    all_pts1.insert(all_pts1.end(), sample_pts1.begin(), sample_pts1.end());
    all_pars1.insert(all_pars1.end(), parvals.begin(), parvals.end());
    all_pars1.insert(all_pars1.end(), sample_pars1.begin(), 
		     sample_pars1.end()); 
   int dim = srf1->dimension();
  ApproxSurf approx(srf1, all_pts1, all_pars1, dim, epsge,
		    0, true, true, (int)sample_pars1.size()/2);

  // Fix boundary curves not crossed by the trimming curve
  vector<int> edge_fix1(4, 0);
  int ccw_idx[] = {3,1,0,2};
  for (size_t kj=0; kj<bd_idx1.size(); kj++)
    edge_fix1[ccw_idx[bd_idx1[kj]]] = 2;

  approx.edgeFix(&edge_fix1[0]);

  double mdist1, mdist2;
  int nsample = 200;

  double maxdist, avdist;
  int nmb_out;
  shared_ptr<SplineSurface> srf1_2 = approx.getApproxSurf(maxdist, avdist,
							   nmb_out);
  if (getenv("DEBUG") && (*getenv("DEBUG")) == '1')
    {
  std::ofstream out1("spline_mod2.g2");
  srf1->writeStandardHeader(out1);
  srf1->write(out1);
  srf1_2->writeStandardHeader(out1);
  srf1_2->write(out1);
  out1 << "400 1 0 4 255 0 0 255" << std::endl;
  out1 << points.size()/3 << std::endl;
  size_t kr;
  for (kr=0; kr<points.size(); kr+=3)
    out1 << points[kr] << "  " << points[kr+1] << "  " << points[kr+2] << std::endl;
  out1 << "400 1 0 4 0 255 0 255" << std::endl;
  out1 << sample_pts1.size()/3 << std::endl;
  for (kr=0; kr<sample_pts1.size(); kr+=3)
    out1 << sample_pts1[kr] << "  " << sample_pts1[kr+1] << "  " << sample_pts1[kr+2] << std::endl;  
    }

  srf1->swap(*srf1_2.get());
  for (idx=0; idx<nmb_bd; ++idx)
    {
      bd_cv1[idx]->updateCurves(epsge);

      checkBoundaryDist(bd_cv1[idx], bd_cv2[idx], start1[idx], 
			end1[idx], start2[idx], end2[idx],
			nsample, mdist1, mdist2);
  if (getenv("DEBUG") && (*getenv("DEBUG")) == '1')
    {
      std::cout << "modifySplineSf, distances: " << mdist1 << ", ";
      std::cout << mdist2 << std::endl;
    }
    }

  }

//===========================================================================
  void GapRemoval::modifySplines(shared_ptr<ParamSurface>& psrf1, 
				 vector<shared_ptr<CurveOnSurface> >& bd_cv1,
				 vector<double>& start1, vector<double>& end1,
				 shared_ptr<ParamSurface>& psrf2, 
				 vector<shared_ptr<CurveOnSurface> >& bd_cv2,
				 vector<double>& start2, vector<double>& end2, 
				 vector<Point>& vertex, double epsge)
//===========================================================================
  {
    // Estimate curve lengths to get a reasonable distribution of
    // sample points
    std::ofstream out0("in_cvs.g2");
    size_t idx, kj;
    size_t nmb_bd = bd_cv1.size();
    vector<double> len(nmb_bd);
    double len_sum = 0.0;
    for (idx=0; idx<nmb_bd; ++idx)
      {
	len[idx] = bd_cv1[idx]->estimatedCurveLength();
	len_sum += len[idx];

  if (getenv("DEBUG") && (*getenv("DEBUG")) == '1')
    {
	bd_cv1[idx]->spaceCurve()->writeStandardHeader(out0);
	bd_cv1[idx]->spaceCurve()->write(out0);
	bd_cv2[idx]->spaceCurve()->writeStandardHeader(out0);
	bd_cv2[idx]->spaceCurve()->write(out0);
    }
      }
    
    vector<double> points;
    vector<double> parval1;
    vector<double> parval2;
    int nmb_sample = 500; //200;
    for (idx=0; idx<nmb_bd; ++idx)
      {
	// Check if the boundary curves have a corresponding orientation
	Point p11 = bd_cv1[idx]->ParamCurve::point(bd_cv1[idx]->startparam());
	Point p12 = bd_cv1[idx]->ParamCurve::point(bd_cv1[idx]->endparam());
	Point p21 = bd_cv2[idx]->ParamCurve::point(bd_cv2[idx]->startparam());
	Point p22 = bd_cv2[idx]->ParamCurve::point(bd_cv2[idx]->endparam());
	bool same_orientation = 
	  (p11.dist(p21) + p12.dist(p22) < p11.dist(p22) + p12.dist(p21));

	// Evaluator based intersection curve
	shared_ptr<IntCrvEvaluator> int_crv = 
	  shared_ptr<IntCrvEvaluator>(new IntCrvEvaluator(bd_cv1[idx], start1[idx], 
							  end1[idx], bd_cv2[idx], 
							  start2[idx], end2[idx],
							  same_orientation));

	// Evaluate the intersection curve in a dense set of points
	double t1 = int_crv->start();
	double t2 = int_crv->end();
	int curr_sample = (int)(nmb_sample*len[idx]/len_sum);
	double tint = (t2 - t1)/(int)(curr_sample-1);
	double tpar;
	int ki;
	for (tpar=t1, ki=0; ki<curr_sample; ++ki, tpar+=tint)
	  {
	    vector<Point> pts = int_crv->eval(tpar);
	    points.insert(points.end(), pts[0].begin(), pts[0].end());
	    parval1.insert(parval1.end(), pts[1].begin(), pts[1].end());
	    parval2.insert(parval2.end(), pts[2].begin(), pts[2].end());
	  }
      }

  // Modify surface 1
    // Fetch spline surface
    vector<shared_ptr<CurveOnSurface> > all_bd1;
    shared_ptr<SplineSurface> srf1 = getSplineAndBd(psrf1, all_bd1);

    // Add sample points along the other boundaries
    vector<double> sample_pts1;
    vector<double> sample_pars1;
    vector<int> bd_idx1;
    getBoundarySamples(all_bd1, bd_cv1, sample_pts1, sample_pars1,
		       bd_idx1, epsge);

  // Make approximator
    // Combine point conditions
    vector<double> all_pts1;
    vector<double> all_pars1;
    all_pts1.insert(all_pts1.end(), points.begin(), points.end());
    all_pts1.insert(all_pts1.end(), sample_pts1.begin(), sample_pts1.end());
    all_pars1.insert(all_pars1.end(), parval1.begin(), parval1.end());
    all_pars1.insert(all_pars1.end(), sample_pars1.begin(), 
		     sample_pars1.end()); 
  int dim = srf1->dimension();
  ApproxSurf approx1(srf1, all_pts1, all_pars1, dim, 0.3*epsge,
		     0, true, true, (int)sample_pars1.size()/2, false);

  // Fix boundary curves not crossed by the trimming curve
  vector<int> edge_fix1(4, 0);
  int ccw_idx[] = {3,1,0,2};
  for (kj=0; kj<bd_idx1.size(); kj++)
    edge_fix1[ccw_idx[bd_idx1[kj]]] = 2;

  approx1.edgeFix(&edge_fix1[0]);

  double maxdist1, avdist1;
  int nmb_out1;
  int max_iter = 7;
  int nmb_keep = 3;
  shared_ptr<SplineSurface> srf1_2 = approx1.getApproxSurf(maxdist1, avdist1,
							   nmb_out1, max_iter,
							   nmb_keep);

  std::ofstream out1("spline_mod1.g2");
  if (getenv("DEBUG") && (*getenv("DEBUG")) == '1')
    {
  srf1->writeStandardHeader(out1);
  srf1->write(out1);
  srf1_2->writeStandardHeader(out1);
  srf1_2->write(out1);
  out1 << "400 1 0 4 255 0 0 255" << std::endl;
  out1 << points.size()/3 << std::endl;
  size_t kr;
  for (kr=0; kr<points.size(); kr+=3)
    out1 << points[kr] << "  " << points[kr+1] << "  " << points[kr+2] << std::endl;
  out1 << "400 1 0 4 0 255 0 255" << std::endl;
  out1 << sample_pts1.size()/3 << std::endl;
  for (kr=0; kr<sample_pts1.size(); kr+=3)
    out1 << sample_pts1[kr] << "  " << sample_pts1[kr+1] << "  " << sample_pts1[kr+2] << std::endl;  
    }

  srf1->swap(*srf1_2.get());
//   for (idx=0; idx<nmb_bd; ++idx)
//     {
//       bd_cv1[idx]->updateCurves(0.3*epsge);
//     }

  // Modify surface 2
    // Fetch spline surface
  vector<shared_ptr<CurveOnSurface> > all_bd2;
  shared_ptr<SplineSurface> srf2 = getSplineAndBd(psrf2, all_bd2);

  // Add sample points along the other boundaries
  vector<double> sample_pts2;
  vector<double> sample_pars2;
  vector<int> bd_idx2;
  getBoundarySamples(all_bd2, bd_cv2, sample_pts2, sample_pars2,
		     bd_idx2, epsge);

  // Make approximator
  vector<double> all_pts2;
  vector<double> all_pars2;
  all_pts2.insert(all_pts2.end(), points.begin(), points.end());
  all_pts2.insert(all_pts2.end(), sample_pts2.begin(), sample_pts2.end());
  all_pars2.insert(all_pars2.end(), parval2.begin(), parval2.end());
  all_pars2.insert(all_pars2.end(), sample_pars2.begin(), 
		     sample_pars2.end()); 
  ApproxSurf approx2(srf2, all_pts2, all_pars2, dim, 0.3*epsge,
		     0, true, true, (int)sample_pars2.size()/2, false);

  // Fix boundary curves not crossed by the trimming curve
  vector<int> edge_fix2(4, 0);
  for (kj=0; kj<bd_idx2.size(); kj++)
    edge_fix2[ccw_idx[bd_idx2[kj]]] = 2;

  approx2.edgeFix(&edge_fix2[0]);

  double maxdist2, avdist2;
  int nmb_out2;
  shared_ptr<SplineSurface> srf2_2 = approx2.getApproxSurf(maxdist2, avdist2,
							   nmb_out2, max_iter,
							   nmb_keep);
  if (getenv("DEBUG") && (*getenv("DEBUG")) == '1')
    {
  srf2->writeStandardHeader(out1);
  srf2->write(out1);
  srf2_2->writeStandardHeader(out1);
  srf2_2->write(out1);
  out1 << "400 1 0 4 0 255 0 255" << std::endl;
  out1 << sample_pts2.size()/3 << std::endl;
  for (size_t kr=0; kr<sample_pts2.size(); kr+=3)
    out1 << sample_pts2[kr] << "  " << sample_pts2[kr+1] << "  " << sample_pts2[kr+2] << std::endl;
    }
   srf2->swap(*srf2_2.get());

  double mdist1, mdist2;
  int nsample = 200;

  for (idx=0; idx<nmb_bd; ++idx)
    {
//       bd_cv2[idx]->updateCurves(0.3*epsge);
	double app_err;
      app_err = removeGapTrim(bd_cv1[idx], start1[idx], end1[idx],
			      bd_cv2[idx], start2[idx], end2[idx],
			      vertex[idx], vertex[idx+1], 0.5*epsge);
//       else
// 	{
// 	  bd_cv1[idx]->updateCurves(0.3*epsge);
// 	  bd_cv2[idx]->updateCurves(0.3*epsge);
// 	}

      checkBoundaryDist(bd_cv1[idx], bd_cv2[idx], start1[idx], 
			end1[idx], start2[idx], end2[idx],
			nsample, mdist1, mdist2);
  if (getenv("DEBUG") && (*getenv("DEBUG")) == '1')
    {
      std::cout << "modifySplines, distances: " << mdist1 << ", ";
      std::cout << mdist2 << std::endl;
      std::cout << "First max: " << maxdist1 << ", med: " << avdist1;
      std::cout << ", out: " << nmb_out1 << std::endl;
      std::cout << "Second max: " << maxdist2 << ", med: " << avdist2;
      std::cout << ", out: " << nmb_out2 << std::endl;
    }
    }

  }

//===========================================================================
void
GapRemoval::removeGapSpline2(vector<shared_ptr<CurveOnSurface> >& bd_cv1,
			     vector<double>& start1, vector<double>& end1,
			     vector<shared_ptr<CurveOnSurface> >& bd_cv2,
			     vector<double>& start2, vector<double>& end2, 
			     vector<Point>& vertex, double epsge)
//===========================================================================
{
  // Make sure that the surfaces have the same degree along
  // the common boundary curves
  // First collect all surfaces and corresponding boundary information
  size_t nmb_bd = bd_cv1.size();
  ASSERT(bd_cv2.size() == nmb_bd);
  vector<pair<shared_ptr<SplineSurface>, pair<int,bool> > > sf_bd1(nmb_bd);
  vector<pair<shared_ptr<SplineSurface>, pair<int,bool> > > sf_bd2(nmb_bd);
  int order = 0;
  size_t idx;
  double ptol = 1.0e-10;
  for (idx=0; idx<nmb_bd; ++idx)
    {
      shared_ptr<SplineSurface> s1 = 
	dynamic_pointer_cast<SplineSurface,ParamSurface>(bd_cv1[idx]->underlyingSurface());
      shared_ptr<SplineSurface> s2 = 
	dynamic_pointer_cast<SplineSurface,ParamSurface>(bd_cv2[idx]->underlyingSurface());
      if (!s1.get() || !s2.get())
	continue;

      int bd1, bd2;  // Specifies the surface boundaries corresponding to 
      // the current edges
      // 0 = umin, 1 = umax, 2 = vmin,  3 = vmax
      bool same_orient1, same_orient2;
      bd1 = bd_cv1[idx]->whichBoundary(epsge, same_orient1);
      bd2 = bd_cv2[idx]->whichBoundary(epsge, same_orient2);
      sf_bd1[idx] = make_pair(s1, make_pair(bd1, same_orient1));
      sf_bd2[idx] = make_pair(s2, make_pair(bd2, same_orient2));
      int order1 = (bd1 == 0 || bd1 == 1) ? s1->order_v() : s1->order_u();
      int order2 = (bd2 == 0 || bd2 == 1) ? s2->order_v() : s2->order_u();
      order = std::max(order, std::max(order1, order2));
    }

  // Raise degree when necessary
  for (idx=0; idx<nmb_bd; idx++)
    {
      int bd1 = sf_bd1[idx].second.first;
      int bd2 = sf_bd2[idx].second.first;
      shared_ptr<SplineSurface> s1 = sf_bd1[idx].first;
      shared_ptr<SplineSurface> s2 = sf_bd2[idx].first;
      int order1 = (bd1 == 0 || bd1 == 1) ? s1->order_v() : s1->order_u();
      int order2 = (bd2 == 0 || bd2 == 1) ? s2->order_v() : s2->order_u();
      if (order1 < order)
	s1->raiseOrder((bd1 == 2 || bd1 == 3) ? order-order1 : 0,
		       (bd1 == 0 || bd1 == 1) ? order-order1 : 0);
      if (order2 < order)
	s2->raiseOrder((bd2 == 2 || bd2 == 3) ? order-order2 : 0,
		       (bd2 == 0 || bd2 == 1) ? order-order2 : 0);
    }
      

  // Fetch boundary curves of spline surfaces
  // Count also the number of different curves and make index arrays
  // related to those
  int indx[] = {3,1,0,2};
  vector<shared_ptr<SplineCurve> > crvs(2*nmb_bd);
  vector<int> cp(2*nmb_bd);
  int ncv1=0, ncv2=0;
  for (idx=0; idx<nmb_bd; idx++)
    {
      int bd1 = sf_bd1[idx].second.first;
      int bd2 = sf_bd2[idx].second.first;
      shared_ptr<SplineSurface> s1 = sf_bd1[idx].first;
      shared_ptr<SplineSurface> s2 = sf_bd2[idx].first;
      if (idx == 0 || sf_bd1[idx-1].first != s1 ||
	  sf_bd1[idx-1].second.first != bd1)
	{
	  crvs[2*idx] = shared_ptr<SplineCurve>(s1->edgeCurve(indx[bd1]));
	  cp[2*idx] = ncv1;
	  ncv1++;
	}
      else
	{
	  crvs[2*idx] = crvs[2*(idx-1)];
	  cp[2*idx] = cp[2*(idx-1)];
	}
      if (idx == 0 || sf_bd2[idx-1].first != s2 || 
	  sf_bd2[idx-1].second.first != bd2)
	{
	  crvs[2*idx+1] = shared_ptr<SplineCurve>(s2->edgeCurve(indx[bd2]));
	  cp[2*idx+1] = ncv2;
	  ncv2++;
	}
      else
	{
	  crvs[2*idx+1] = crvs[2*idx-1];
	  cp[2*idx+1] = cp[2*idx-1];
	}
    }

  // Remove outer pieces in the curve chains where the curve has no partner
  shared_ptr<SplineCurve> c1_1 = crvs[0];
  shared_ptr<SplineCurve> c1_2 = crvs[2*(nmb_bd-1)];
  shared_ptr<SplineCurve> c2_1 = crvs[1];
  shared_ptr<SplineCurve> c2_2 = crvs[2*nmb_bd-1];
  int bd1 = sf_bd1[0].second.first;
  // int bd2 = sf_bd2[0].second.first;
  int dir = (bd1 <= 1) ? 1 : 0;
  Point p1 = bd_cv1[0]->faceParameter(start1[0]);
  Point p2 = bd_cv2[0]->faceParameter(start2[0]);
  Point p3 = bd_cv1[nmb_bd-1]->faceParameter(end1[nmb_bd-1]);
  Point p4 = bd_cv2[nmb_bd-1]->faceParameter(end2[nmb_bd-1]);
  bool keep_der1=false, keep_der2=false, keep_der3=false, keep_der4=false;
  if ((p3[dir] - p1[dir])*(p4[dir] - p2[dir]) <= 0.0)
    {
      p2 = bd_cv2[0]->faceParameter(end2[0]);
      p4 = bd_cv2[nmb_bd-1]->faceParameter(start2[nmb_bd-1]);
    }
  double subpar1, subpar2, subpar3, subpar4;
  if (std::min(p1[dir],p3[dir]) > c1_1->startparam()+DEFAULT_PARAMETER_EPSILON)
    {
      subpar1 = std::min(p1[dir],p3[dir]);
      shared_ptr<SplineCurve> sub_crv = 
	shared_ptr<SplineCurve>(c1_1->subCurve(subpar1, c1_1->endparam()));
      for (size_t kr=0; kr<nmb_bd; ++kr)
	{
	  if (crvs[2*kr].get() != c1_1.get())
	    break;
	  crvs[2*kr] = sub_crv;
	}
      keep_der1 = true;
    }
  dir = (sf_bd2[0].second.first == 0 || sf_bd2[0].second.first == 1) ? 1 : 0;
  if (std::min(p2[dir],p4[dir]) > c2_1->startparam()+DEFAULT_PARAMETER_EPSILON)
    {
      subpar2 = std::min(p2[dir],p4[dir]);
      shared_ptr<SplineCurve> sub_crv = 
	shared_ptr<SplineCurve>(c2_1->subCurve(subpar2, c2_1->endparam()));
      for (size_t kr=0; kr<nmb_bd; ++kr)
	{
	  if (crvs[2*kr+1].get() != c2_1.get())
	    break;
	  crvs[2*kr+1] = sub_crv;
	}
      keep_der3 = true;
    }
  dir = (sf_bd1[nmb_bd-1].second.first == 0 || 
	 sf_bd1[nmb_bd-1].second.first == 1) ? 1 : 0;
  if (std::max(p1[dir],p3[dir]) < c1_2->endparam()-DEFAULT_PARAMETER_EPSILON)
    {
      subpar3 = std::max(p1[dir],p3[dir]);
      shared_ptr<SplineCurve> sub_crv = 
	shared_ptr<SplineCurve>(c1_2->subCurve(c1_2->startparam(), subpar3));
      for (int kr=nmb_bd-1; kr>=0; --kr)
	{
	  if (crvs[2*kr].get() != c1_2.get())
	    break;
	  crvs[2*kr] = sub_crv;
	}
      keep_der2 = true;
    }
  dir = (sf_bd2[nmb_bd-1].second.first == 0 || 
	 sf_bd2[nmb_bd-1].second.first == 1) ? 1 : 0;
  if (std::max(p2[dir],p4[dir]) < c2_2->endparam()-DEFAULT_PARAMETER_EPSILON)
    {
      subpar4 = std::max(p2[dir],p4[dir]);
      shared_ptr<SplineCurve> sub_crv = 
	shared_ptr<SplineCurve>(c2_2->subCurve(c2_2->startparam(), subpar4));
      for (int kr=nmb_bd-1; kr>=0; --kr)
	{
	  if (crvs[2*kr+1].get() != c2_2.get())
	    break;
	  crvs[2*kr+1] = sub_crv;
	}
      keep_der4 = true;
    }

  vector<int> seam(ncv1+ncv2, 0);
  vector<vector<int> > coef_known(ncv1+ncv2);
  vector<shared_ptr<sideConstraintSet> > constraints;
  vector<vector<double> > inter_point(ncv1+ncv2);
  vector<vector<double> > inter_par(ncv1+ncv2);
  vector<vector<int> > inter_der(ncv1+ncv2);
  vector<bool> turn(idx);
  bool interpolate = false;

  for (idx=0; idx<nmb_bd; idx++)
    {
      // Ensure same orientation of curves to simplify definition of constraints
      Point f1_p1 = bd_cv1[idx]->faceParameter(start1[idx]);
      Point f1_p2 = bd_cv1[idx]->faceParameter(end1[idx]);
      Point f2_p1 = bd_cv2[idx]->faceParameter(start2[idx]);
      Point f2_p2 = bd_cv2[idx]->faceParameter(end2[idx]);

      bool opposite = false;
      int bd1 = sf_bd1[idx].second.first;
      int bd2 = sf_bd2[idx].second.first;
      double t1 = (bd1 == 0 ||  bd1 == 1) ? f1_p1[1] : f1_p1[0];
      double t2 = (bd1 == 0 ||  bd1 == 1) ? f1_p2[1] : f1_p2[0];
      double t3 = (bd2 == 0 ||  bd2 == 1) ? f2_p1[1] : f2_p1[0];
      double t4 = (bd2 == 0 ||  bd2 == 1) ? f2_p2[1] : f2_p2[0];
      if ((t2 - t1)*(t4 -t3) < 0.0)
	opposite = true;

      bool same_orient1 = sf_bd1[idx].second.second;
      bool same_orient2 = sf_bd2[idx].second.second;
      if ((same_orient1 && !same_orient2) || (!same_orient1 && same_orient2))
	opposite = !opposite;
      turn[idx] = opposite;
      if (!same_orient1)
	std::swap(t1, t2);
      if (!same_orient2)
	std::swap(t3, t4);
      if (opposite && (idx>0 && crvs[2*idx+1].get() != crvs[2*idx-1].get()))
	  crvs[2*idx+1]->reverseParameterDirection();
	
 
      if (fabs(t1 - crvs[2*idx]->startparam()) < ptol)
	t1 = crvs[2*idx]->startparam();
      if (fabs(t2 - crvs[2*idx]->endparam()) < ptol)
	t2 = crvs[2*idx]->endparam();
      if (fabs(t3 - crvs[2*idx+1]->startparam()) < ptol)
	t3 = crvs[2*idx+1]->startparam();
      if (fabs(t4 - crvs[2*idx+1]->endparam()) < ptol)
	t4 = crvs[2*idx+1]->endparam();

     // Set coefficient or interpolation condition for for the vertices
      bool atcorner1[2];
      bool atcorner2[2];
      atcorner1[0] = (fabs(crvs[2*idx]->startparam() - t1) < ptol);
      atcorner2[0] = (fabs(crvs[2*idx]->endparam() - t2) < ptol);
      atcorner1[1] = (fabs(crvs[2*idx+1]->startparam() - t3) < ptol);
      atcorner2[1] = (fabs(crvs[2*idx+1]->endparam() - t4) < ptol);
      Point vx1 = (same_orient1) ? vertex[idx] : vertex[idx+1];
      Point vx2 = (same_orient1) ? vertex[idx+1] : vertex[idx];
      for (int ki=0; ki<2; ++ki)
	{
	  if (atcorner1[ki])
	    {
	      crvs[2*idx+ki]->makeKnotStartRegular();
	      crvs[2*idx+ki]->replaceEndPoint(vx1, true);
	    }
	  else
	    {
	      int ix = ki*(ncv1+cp[2*idx+1]) + (1-ki)*cp[2*idx];
	      inter_point[ix].insert(inter_point[ix].end(), 
					   vx1.begin(), vx1.end());
	      inter_par[ix].push_back((ki==0) ? t1 : t3);
	      inter_der[ix].push_back(0);
	      interpolate = true;
	    }

	  if (atcorner2[ki])
	    {
	      crvs[2*idx+ki]->makeKnotEndRegular();
	      crvs[2*idx+ki]->replaceEndPoint(vx2, false);
	    }
	  else
	    {
	      int ix = ki*(ncv1+cp[2*idx+1]) + (1-ki)*cp[2*idx];
	      inter_point[ix].insert(inter_point[ix].end(), 
					   vx2.begin(), vx2.end());
	      inter_par[ix].push_back((ki==0) ? t2 : t4);
	      inter_der[ix].push_back(0);
	      interpolate = true;
	    }
	}


      // Ensure same parameter domain to simplify definition of constraints
      double fac = (t2 - t1)/(t4 - t3);
      double init1 = crvs[2*idx+1]->startparam();
      double init2 = crvs[2*idx+1]->endparam();
      double par1 = t1 - (t3 - init1)*fac;
      double par2 = t2 + (init2 - t4)*fac;
      crvs[2*idx+1]->setParameterInterval(par1, par2);

       // Define constraints between coefficients to remove the gap
      vector<shared_ptr<sideConstraintSet> > curr_constraints =
	getCoefConstraints(crvs[2*idx], cp[2*idx], crvs[2*idx+1], 
			   ncv1+cp[2*idx+1], ptol);
      constraints.insert(constraints.end(), curr_constraints.begin(),
			 curr_constraints.end());

      // Redo modifications parameter interval
      crvs[2*idx+1]->setParameterInterval(init1, init2);
      for (int ki=0; ki<2; ++ki)
	{
	  // Set information about free and fixed coefficients
	  int ix = ki*(ncv1+cp[2*idx+1]) + (1-ki)*cp[2*idx];
	  coef_known[ix].resize(crvs[2*idx+ki]->numCoefs());
	  std::fill(coef_known[ix].begin(), coef_known[ix].end(), 0);
	  coef_known[ix][0] = 
	    coef_known[ix][coef_known[ix].size()-1] = 1;
	}
    }
  if (keep_der1)
    coef_known[0][1] = 1;
  if (keep_der2)
    coef_known[0][coef_known[0].size()-2] = 1;
  if (keep_der3)
    coef_known[nmb_bd-1][1] = 1;
  if (keep_der4)
    coef_known[nmb_bd-1][coef_known[0].size()-2] = 1;

  //smooth.attach(crvs, seam, coef_known, constraints.size());
  // Modify curves. First make a unique set of curves
  vector<shared_ptr<SplineCurve> > crvs2(ncv1+ncv2);
  size_t kr, kh;
  for (idx=0, kr=0, kh=0; idx<nmb_bd; ++idx)
    {
      if (idx == 0 || crvs[2*idx].get() != crvs[2*(idx-1)].get())
	crvs2[kr++] = crvs[2*idx];
      if (idx == 0 || crvs[2*idx+1].get() != crvs[2*idx-1].get())
	{
	  crvs2[ncv1+kh] = crvs[2*idx+1];
	  kh++;
	}
    }
  
  if (getenv("DEBUG") && (*getenv("DEBUG")) == '1')
    {
  std::ofstream out("cv_mod3.g2");
  for (size_t kr=0; kr<crvs2.size(); ++kr)
    {
      crvs2[kr]->writeStandardHeader(out);
      crvs2[kr]->write(out);
    }
  out << "400 1 0 4 255 0 0 255" << std::endl;
  out << vertex.size() << std::endl;
  for (size_t kr=0; kr<vertex.size(); ++kr)
    out << vertex[kr] << " " << std::endl;
    }

  // Perform modification
  SmoothCurveSet smooth;
  smooth.attach(crvs2, seam, coef_known);

  double wgt_orig = 0.8;
  double wgt2 = 0.5*(1.0 - wgt_orig);
  double wgt3 = 0.5*(1.0 - wgt_orig);
  smooth.setOptimize(0.0, wgt2, wgt3);

  smooth.setApproxOrig(wgt_orig);

  int status;
  if (interpolate)
    smooth.setInterpolationConditions(inter_point, inter_par, inter_der,
				      false, wgt_orig, &status);

  // Must from somewhere call updateSideConstraints first (or inside setSideConstraints)
  bool exact = true; // false;
  if (exact)
    smooth.setSideConstraints(constraints, exact);
  else
    {
      double constraint_wgt = 1.0;
      smooth.setApproxSideConstraints(constraints, constraint_wgt);
    }

  vector<shared_ptr<SplineCurve> > updated_crvs;
  status = smooth.equationSolve(updated_crvs);
  if (status != 0)
    {
  if (getenv("DEBUG") && (*getenv("DEBUG")) == '1')
    {
    std::cout << "Something wrong with curve smoothing" << std::endl;
    }
    }

  if (getenv("DEBUG") && (*getenv("DEBUG")) == '1')
    {
  std::ofstream out2("cv_mod3_2.g2");
  for (size_t kr=0; kr<updated_crvs.size(); ++kr)
    {
      out2 << "100 1 0 4 255 0 0 255" << std::endl;
      updated_crvs[kr]->write(out2);
    }
    }


  if ((int)updated_crvs.size() == ncv1+ncv2)
    {
      for (idx=0; idx<nmb_bd; ++idx)
	{
	  if (turn[idx])
	    updated_crvs[ncv1+cp[2*idx+1]]->reverseParameterDirection();
	}

      // Check if the first or last curves must be extended at their endpoints
      if (crvs2[0].get() != c1_1.get())
	updated_crvs[0] = replaceCurvePiece(c1_1, updated_crvs[0], subpar1, 1,
					    c1_1->endparam(), 1);

      if (crvs2[ncv1].get() != c2_1.get())
	updated_crvs[ncv1] = replaceCurvePiece(c2_1, updated_crvs[ncv1], 
					       subpar2, 1,
					       c2_1->endparam(), 1);

      if (crvs2[ncv1-1].get() != c1_2.get())
	updated_crvs[ncv1-1] = replaceCurvePiece(c1_2, updated_crvs[ncv1-1], 
						 c1_2->startparam(), 1, 
						 subpar3, 1);


      if (crvs2[ncv1+ncv2-1].get() != c2_2.get())
	updated_crvs[ncv1+ncv2-1] = replaceCurvePiece(c2_2, 
						      updated_crvs[ncv1+ncv2-1], 
						      c2_2->startparam(), 1,
						      subpar4, 1);

      
      // Update surfaces
      std::ofstream out1("spline_mod3.g2");
  if (getenv("DEBUG") && (*getenv("DEBUG")) == '1')
    {
      for (size_t kr=0; kr<updated_crvs.size(); ++kr)
	{
	  out1 << "100 1 0 4 0 255 0 255" << std::endl;
	  updated_crvs[kr]->write(out1);
	}
    }
      for (idx=0; idx<nmb_bd; ++idx)
	{
	  int bd1 = sf_bd1[idx].second.first;
	  int bd2 = sf_bd2[idx].second.first;
	  shared_ptr<SplineSurface> s1 = sf_bd1[idx].first;
	  shared_ptr<SplineSurface> s2 = sf_bd2[idx].first;
      
  if (getenv("DEBUG") && (*getenv("DEBUG")) == '1')
    {
	  s1->writeStandardHeader(out1);
	  s1->write(out1);
	  s2->writeStandardHeader(out1);
	  s2->write(out1);
    }
 	  if (idx == 0 || cp[2*idx] != cp[2*(idx-1)])
	    s1->replaceBoundaryCurve(bd1, updated_crvs[cp[2*idx]]);
	  if (idx == 0 || cp[2*idx+1] != cp[2*idx-1])
	    s2->replaceBoundaryCurve(bd2, updated_crvs[ncv1+cp[2*idx+1]]);

  if (getenv("DEBUG") && (*getenv("DEBUG")) == '1')
    {
	  s1->writeStandardHeader(out1);
	  s1->write(out1);
	  s2->writeStandardHeader(out1);
	  s2->write(out1);
    }
	}
    }


  // Update boundary curves
  double mdist1, mdist2;
  int nmb_sample = 200;

  for (idx=0; idx<nmb_bd; ++idx)
    {
	bool updated1, updated2;
      updated1 = bd_cv1[idx]->updateIsoCurves();
      updated2 = bd_cv2[idx]->updateIsoCurves();

      checkBoundaryDist(bd_cv1[idx], bd_cv2[idx], start1[idx], 
			end1[idx], start2[idx], end2[idx],
			nmb_sample, mdist1, mdist2);
   if (getenv("DEBUG") && (*getenv("DEBUG")) == '1')
    {
     std::cout << "modifySplines, distances: " << mdist1 << ", ";
      std::cout << mdist2 << std::endl;
    }
    }
}


//===========================================================================
  vector<shared_ptr<sideConstraintSet> > 
  GapRemoval::getCoefConstraints(shared_ptr<SplineCurve>& crv1, int idx1,
				 shared_ptr<SplineCurve>& crv2, int idx2,
				 double tol)
//===========================================================================
  {
    vector<shared_ptr<sideConstraintSet> > constraints;

    if (!(crv1->order() == crv2->order()))
      return constraints;

    // Compute union knot vector
    vector<double> union_knots;
    vector<BsplineBasis> bbasis(2);
    bbasis[0] = crv1->basis();
    bbasis[1] = crv2->basis();

    GeometryTools::makeUnionKnots(bbasis, tol, union_knots);
    
    int kn = (int)union_knots.size() - crv1->order();  // Number of coefficients in
    // refined curve
    constraints.resize(kn); // Number of constraints
    
    // Set constraints related to crv1
    int ki, kj, kr, kmy;
    vector<double>::iterator st1 = crv1->knotsBegin();
    int kpl, kfi, kla;
    vector<double> galfa(crv1->order());
    for (ki=0, kmy=0; ki<kn; ++ki)
      {
	constraints[ki] = 
	  shared_ptr<sideConstraintSet>(new sideConstraintSet(crv1->dimension()));
	while (kmy < crv1->numCoefs()+crv1->order()-1 &&
	       st1[kmy+1] <= union_knots[ki])
	  kmy++;
	if (kmy == crv1->numCoefs()+crv1->order())
	    break;

	if (union_knots[ki] < st1[kmy])
	  continue;

	SplineUtils::osloalg(ki, kmy, crv1->order(), crv1->numCoefs(), 
		&kpl, &kfi, &kla,
		&union_knots[0], &st1[0], &galfa[0]);

	for (kj=kfi, kr=kfi+kpl; kj<=kla; ++kj, ++kr)
	  constraints[ki]->factor_.push_back(std::make_pair(std::make_pair(idx1, kj),
							    galfa[kr]));
      }
		
    // Set constraints related to crv2
    vector<double>::iterator st2 = crv2->knotsBegin();
    for (ki=0, kmy=0; ki<kn; ++ki)
      {
	while (kmy < crv2->numCoefs()+crv2->order()-1 &&
		st2[kmy+1] <= union_knots[ki])
	  kmy++;

	if (kmy == crv2->numCoefs()+crv2->order())
	  break;

	if (union_knots[ki] < st2[kmy])
	  continue;

	SplineUtils::osloalg(ki, kmy, crv2->order(), crv2->numCoefs(), 
		&kpl, &kfi, &kla,
		&union_knots[0], &st2[0], &galfa[0]);

	for (kj=kfi, kr=kfi+kpl; kj<=kla; ++kj, ++kr)
	  constraints[ki]->factor_.push_back(std::make_pair(std::make_pair(idx2, kj),
							    -galfa[kr]));
      }

    // Remove constraints containing only one curve
    for (ki=0; ki < (int)constraints.size();)
      {
	if (constraints[ki]->factor_.size() == 0)
	  {
	    constraints.erase(constraints.begin()+ki);
	  }
	else
	  {
	    int idx1 = constraints[ki]->factor_[0].first.first;
	    size_t kj;
	    for (kj=1; kj<constraints[ki]->factor_.size(); ++kj)
	      {
		int idx2 = constraints[ki]->factor_[kj].first.first;
		if (idx1 != idx2)
		  break;
	      }
	    if (kj == constraints[ki]->factor_.size())
	      {
		constraints.erase(constraints.begin()+ki);
	      }
	    else
	      ki++;
	  }
      }
	
	  
	
    return constraints;
  }

//===========================================================================
shared_ptr<SplineCurve>
  GapRemoval::replaceCurvePiece(shared_ptr<SplineCurve> crv,
				shared_ptr<SplineCurve> sub_crv,
				double par1, int cont1, 
				double par2, int cont2)
//===========================================================================
{
  crv->basis().check();
  sub_crv->basis().check();

  // Enforce only a reasonable continuity at the joints
  cont1 = std::min(cont1, 2);
  cont2 = std::min(cont2, 2);

  // Create a new curve where the current boundary curve piece is replaced
  // with the updated curve
  // First make C0 curve
  shared_ptr<SplineCurve> new_crv;
  (void)crv->basis().knotIntervalFuzzy(par1);
  (void)crv->basis().knotIntervalFuzzy(par2);
  int ins1 = crv->order() - crv->basis().knotMultiplicity(par1) - 1;
  int ins2 = crv->order() - crv->basis().knotMultiplicity(par2) - 1;
  double dist1, dist2;
  if (ins1 > 0)
    {
      new_crv = shared_ptr<SplineCurve>(crv->subCurve(crv->startparam(), par1));
      new_crv->appendCurve(sub_crv.get(), cont1, dist1, false);
    }
  else
    new_crv = sub_crv;
  if (ins2 > 0)
    {
      shared_ptr<SplineCurve> tmp_crv = 
	shared_ptr<SplineCurve>(crv->subCurve(par2, crv->endparam()));
      new_crv->appendCurve(tmp_crv.get(), cont2, dist2, false);
    }

  new_crv->basis().check();

  return new_crv;
}

//===========================================================================
shared_ptr<SplineSurface> 
GapRemoval::getSplineAndBd(shared_ptr<ParamSurface> psurf,
			   vector<shared_ptr<CurveOnSurface> >& bd_crvs)
//===========================================================================
  {
    shared_ptr<BoundedSurface> bd_sf = 
      dynamic_pointer_cast<BoundedSurface, ParamSurface>(psurf);
    shared_ptr<SplineSurface> srf =
      dynamic_pointer_cast<SplineSurface, ParamSurface>(psurf);

    bool is_spline = false;
    if (bd_sf.get())
      {
	is_spline = bd_sf->hasUnderlyingSpline(srf);
	vector<CurveLoop> loops = bd_sf->allBoundaryLoops();
	for (size_t ki=0; ki<loops.size(); ++ki)
	  {
	    int nmb_cvs = loops[ki].size();
	    for (int kj=0; kj<nmb_cvs; ++kj)
	      {
		shared_ptr<CurveOnSurface> cv =
		  dynamic_pointer_cast<CurveOnSurface, ParamCurve>(loops[ki][kj]);
		if (cv.get())
		  bd_crvs.push_back(cv);
	      }
	  }
      }
    return srf;
  }

//===========================================================================
void 
GapRemoval::getBoundarySamples(vector<shared_ptr<CurveOnSurface> >& all_bd,
			       vector<shared_ptr<CurveOnSurface> >& bd_cvs,
			       vector<double>& pts, vector<double>& pars,
			       vector<int>& bd_idx, double epsge)
//===========================================================================
{
  size_t idx, ki;
  size_t nmb_bd = all_bd.size();
  vector<double> len(nmb_bd);
  double len_sum = 0.0;
  for (idx=0; idx<nmb_bd; ++idx)
    {
      len[idx] = all_bd[idx]->estimatedCurveLength();
      len_sum += len[idx];
    }
  int nmb_sample = 1500; //800;

  // Note boundaries that cannot be fixed
  vector<int> bd_use;
  for (idx=0; idx<bd_cvs.size(); ++idx)
    {
      bool same;
      int boundary = bd_cvs[idx]->whichBoundary(epsge, same);
      if (boundary >= 0)
	bd_use.push_back(boundary);
    }
	
  for (idx=0; idx<nmb_bd; ++idx)
    {
      for (ki=0; ki<bd_cvs.size(); ++ki)
	if (bd_cvs[ki] == all_bd[idx])
	  break;

      if (ki < bd_cvs.size())
	continue;  // Curve sampled elsewhere

      // Check if the curve follows a boundary. In that case
      // report the boundary instead of sampling the curve
      bool same;
      int boundary = all_bd[idx]->whichBoundary(epsge, same);
      
      // Check if the boundary is already in use
      for (ki=0; ki<bd_use.size(); ++ki)
	if (boundary == bd_use[ki])
	  {
	    boundary = -1;
	    break;
	  }

//       if (boundary >= 0)
//       {
// 	bd_idx.push_back(boundary);
// 	continue;
//       }

      int curr_sample = nmb_sample*((int)len[idx])/((int)len_sum);
      double t1 = all_bd[idx]->startparam();
      double t2 = all_bd[idx]->endparam();
      double tint = (t2 - t1)/(int)(curr_sample-1);
      double tpar;
      int kj;

      // Define evaluator based curve corresponding to the unspecified surface
      shared_ptr<TrimCurve> trim_crv = 
	shared_ptr<TrimCurve>(new TrimCurve(all_bd[idx].get(), t1, t2));
      for (tpar=t1, kj=0; kj<curr_sample; ++kj, tpar+=tint)
	{
	  vector<Point> pnt = trim_crv->eval(tpar);
	  pts.insert(pts.end(), pnt[0].begin(), pnt[0].end());
	  pars.insert(pars.end(), pnt[1].begin(), pnt[1].end());
	}
    }
	
}

//===========================================================================
bool 
GapRemoval::modifyAtVertex(shared_ptr<SplineSurface> srf,
			   Point face_param, Point vertex,
			   double epsge)
//===========================================================================
{
  // Improve vertex position
  double u_par, v_par, dist;
  Point pos;
  srf->closestPoint(vertex, u_par, v_par, pos, dist, epsge, NULL, 
		    face_param.begin());

  if (dist < 0.5*epsge)
    return false;  // No point in updating

  // Store initial surface
  shared_ptr<SplineSurface> tmp_srf = shared_ptr<SplineSurface>(srf->clone());

  // Check configuration
  RectDomain dom = srf->parameterDomain();
  Array<double,2> par(u_par, v_par);
  BoundingBox box = srf->boundingBox();
  double eps2d = 2.0*epsge*(dom.lowerLeft().dist(dom.upperRight()))/(box.low().dist(box.high()));
  int bd_idx;

  static double smooth_wgt = 0.01;
  srf->getBoundaryIdx(pos, epsge, bd_idx);
  int in1 = srf->numCoefs_u();
  int in2 = srf->numCoefs_v();
  if (bd_idx != -1 && dom.isOnCorner(par, eps2d))
    {
      // Change corner coefficient to correspond with the vertex
      int idx;
      if (u_par - srf->startparam_u() < srf->endparam_u() - u_par)
	{
	  if (v_par - srf->startparam_v() < srf->endparam_v() - v_par)
	    idx = 0;
	  else 
	    idx = in1*(in2-1);
	}
      else
	{
	  if (v_par - srf->startparam_v() < srf->endparam_v() - v_par)
	    idx = in1-1;
	  else 
	    idx = in1*in2-1;
	}
      if (srf->rational())
	{
	  vector<double>::iterator rc = srf->rcoefs_begin();
	  vector<double>::iterator cc = srf->coefs_begin();
	  int dim = srf->dimension();
	  for (int ki=0; ki<dim; ++ki)
	    {
	      cc[idx*dim+ki] = vertex[ki];
	      rc[idx*(dim+1)+ki] = vertex[ki]*rc[idx*(dim+1)+dim];
	    }
	}
      else
	{
	  vector<double>::iterator cc = srf->coefs_begin();
	  int dim = srf->dimension();
	  for (int ki=0; ki<dim; ++ki)
	      cc[idx*dim+ki] = vertex[ki];
	}

    }
  else if (bd_idx != -1)
    {
      int bd = (bd_idx <= 1) ? bd_idx+2 : bd_idx-2;
      int idx[] = {3,1,0,2};
      int ccw_nmb = idx[bd];
      shared_ptr<SplineCurve> crv = 
	shared_ptr<SplineCurve>(srf->edgeCurve(ccw_nmb));
      double parval = (bd <= 1) ? v_par : u_par;

      int c1, c2;
      crv->basis().coefsAffectingParam(parval, c1, c2);
      int ki1 = std::max(1,c1);
      int ki2 = std::min(crv->numCoefs()-2,c2);
      if (ki2 >= ki1)
	{
	  SmoothCurve smooth;
	  vector<int> coef_known(crv->numCoefs(), 1);

	  for (int ki=ki1; ki<=ki2; ++ki)
	    coef_known[ki] = 0;

	  smooth.attach(crv, &coef_known[0]);
	  smooth.setOptim(0.0, 0.5*smooth_wgt, 0.5*smooth_wgt);
	  vector<double> pnt(vertex.begin(), vertex.end());
	  vector<double> par_pnt(1, parval);
	  vector<double> wgt_pnt(1, 1.0);
	  smooth.setLeastSquares(pnt, par_pnt, wgt_pnt, 1.0-smooth_wgt);

	  shared_ptr<SplineCurve> mod_crv;
	  smooth.equationSolve(mod_crv);

	  bool replaced;
	  replaced = srf->replaceBoundaryCurve(bd, mod_crv, false);
	}
    }
  else
    {
      int c1, c2, c3, c4;
      srf->basis_u().coefsAffectingParam(u_par, c1, c2);
      srf->basis_v().coefsAffectingParam(v_par, c3, c4);
      if (std::min(in2-2,c4) >= std::max(1,c3) &&
	  std::min(in1-2,c2) >= std::max(1,c1))
	{
	  SmoothSurf smooth;
	  vector<int> coef_known(in1*in2, 1);

	  for (int kj=std::max(1,c3); kj<=std::min(in2-2,c4); ++kj)
	    for (int ki=std::max(1,c1); ki<=std::min(in1-2,c2); ++ki)
	      coef_known[kj*in1+ki] = 0;

	  static double orig_wgt = 0.05;
	  int seam[2];
	  seam[0] = seam[1] = 0;
	  smooth.attach(srf, seam, &coef_known[0]);
	  smooth.setOptimize(0.0, 0.5*smooth_wgt, 0.5*smooth_wgt);
	  vector<double> pnt(vertex.begin(), vertex.end());
	  vector<double> par_pnt(2);
	  par_pnt[0] = u_par;
	  par_pnt[1] = v_par;
	  vector<double> wgt_pnt(1, 1.0);
	  smooth.setLeastSquares(pnt, par_pnt, wgt_pnt, 1.0-smooth_wgt-orig_wgt);
	  smooth.approxOrig(orig_wgt);

	  shared_ptr<SplineSurface> mod_sf;
	  smooth.equationSolve(mod_sf);

	  srf->swap(*mod_sf);
	}
    }

  // Test result
  Point pos2 = srf->ParamSurface::point(u_par, v_par);
  double dist2 = pos2.dist(pos);
  double dist3 = pos2.dist(vertex);
  double dist4 = vertex.dist(pos);
  if (getenv("DEBUG") && (*getenv("DEBUG")) == '1')
    {
  std::cout << "Vertex: " << vertex << ", new pos: " << pos2;
  std::cout << ", dist: " << dist3 << ", to prev: " << dist2;
  std::cout << ", initial: " << dist4 << std::endl;
    }

  if (dist3 > dist || dist3 > epsge)
    {
      // Return to initial surface
      srf->swap(*tmp_srf);
      return false;
    }

  return true;
}


//===========================================================================
void 
GapRemoval::checkBoundaryDist(shared_ptr<CurveOnSurface> bd1,
			      shared_ptr<CurveOnSurface> bd2,
			      double start1, double end1,
			      double start2, double end2,
			      int nmb_sample, double& mdist1,
			      double& mdist2)
//===========================================================================
{
  int ki;
  double tint = (end1 - start1)/(double)(nmb_sample-1);
  double par1, par2;
  Point pt1, pt2;
  double dist;
  double avdist1=0.0, avdist2=0.0;
  double dd1 = 0.0, dd2 = 0.0;
  
  shared_ptr<TrimCurve> trim_crv1 = 
    shared_ptr<TrimCurve>(new TrimCurve(bd1.get(), start1, end1));
  shared_ptr<TrimCurve> trim_crv2 = 
    shared_ptr<TrimCurve>(new TrimCurve(bd2.get(), start2, end2));

  // Just to test
  vector<shared_ptr<CurveOnSurface> > bdvec1(1, bd1);
  vector<shared_ptr<CurveOnSurface> > bdvec2(1, bd2);
  vector<double> startvec1(1, start1);
  vector<double> startvec2(1, start2);
  vector<double> endvec1(1, end1);
  vector<double> endvec2(1, end2);
  vector<bool> opp(1, false);
  bool same1;
  int bd_1 = bd1->whichBoundary(1.0e-4, same1);
  shared_ptr<SpaceIntCrv> mid = 
    shared_ptr<SpaceIntCrv>(new SpaceIntCrv(bd1->spaceCurve(), (bd_1 <= 1) ? 1 : 0,
					    bdvec1, startvec1, endvec1, bdvec2,
					    startvec2, endvec2, opp, true));
  mdist1 = mdist2 = 0.0;
  for (ki=0, par1=start1; ki<nmb_sample; ++ki, par1+=tint)
    {
      pt1 = bd1->ParamCurve::point(par1);
      bd2->closestPoint(pt1, start2, end2, par2, pt2, dist);

      mdist1 = std::max(mdist1, dist);
      avdist1 += dist;

      vector<Point> der1 = trim_crv1->eval(par1);
      vector<Point> der2 = trim_crv2->eval(par2);
      dist = der1[0].dist(der2[0]);
      mdist2 = std::max(mdist2, dist);
      avdist2 += dist;

      Point mid_pt = mid->eval(par1);

      dd1 = std::max(dd1, pt1.dist(der1[0]));
      dd2 = std::max(dd2, pt2.dist(der2[0]));
    }
  avdist1 /= (double)nmb_sample;
  avdist2 /= (double)nmb_sample;

  if (getenv("DEBUG") && (*getenv("DEBUG")) == '1')
    {
  std::cout << "Average distances: " << avdist1 << ", " << avdist2;
  std::cout << ", between same curve: " << dd1 << ", " << dd2 << std::endl;
    }
}

  } // end namespace Go
