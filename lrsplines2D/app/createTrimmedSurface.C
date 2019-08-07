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

#include "GoTools/lrsplines2D/TrimUtils.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/utils/Point.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/RectDomain.h"
#include "GoTools/geometry/LoopUtils.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/lrsplines2D/TrimCrvUtils.h"
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

using namespace Go;

using std::cout;
using std::endl;
using std::vector;
using std::string;


int main(int argc, char* argv[])
{
  if (argc != 5)
    {
      cout << "Usage: " << argv[0] << " sf.g2 point_cloud.g2 tightness bd_sf.g2" << endl;
      return -1;
    }

    
  std::ifstream filein_sf(argv[1]);
  std::ifstream filein_points(argv[2]);
  int tightness = atoi(argv[3]);
  std::ofstream fileout_bd_sf(argv[4]);

  GoTools go_tools;
  go_tools.init();
  Registrator<LRSplineSurface> r293;

  ObjectHeader header;
  header.read(filein_sf);
  shared_ptr<LRSplineSurface> sf;
  shared_ptr<GeomObject> geom_obj(Factory::createObject(header.classType()));
  geom_obj->read(filein_sf);
  sf = dynamic_pointer_cast<LRSplineSurface, GeomObject>(geom_obj);
  if (!sf.get())
    {
      std::cerr << "Input file contains no LR B-spline surface" << std::endl;
      exit(-1);
    }

  header.read(filein_points);
  PointCloud3D points;
  points.read(filein_points);

  // Check for correspondance between the point set and the parameter
  // domain of the surface

  BoundingBox box = points.boundingBox();
  Point low = box.low();
  Point high = box.high();

  RectDomain dom = sf->containingDomain();
  double umin = dom.umin();
  double umax = dom.umax();
  double vmin = dom.vmin();
  double vmax = dom.vmax(); 
  
  if (umin > low[0] || umax < high[0] ||  vmin > low[1] || vmax < high[1])
    {
      std::cout << " Point cloud not corresponding to surface domain. Exiting " << std::endl;
      exit(-1);
    }

  // Translate to origo
  Vector3D vec(-0.5*(umin+umax), -0.5*(vmin+vmax), 0.0);
  points.translate(vec);

  std::ofstream of1("translated_cloud.g2");
  points.writeStandardHeader(of1);
  points.write(of1);

  // Update parameter domain
  Point tmp_vec(vec.begin(), vec.end());
  TrimCrvUtils::translateSurfaceDomain(sf.get(), tmp_vec);

  if (sf->dimension() == 3)
    {
      sf->translate(tmp_vec);
    }

  // Set parameters for computations of trimming sequence
  int max_rec;
  int nmb_div;
  if (tightness <= 2)
    {
      max_rec = 1;
      nmb_div = (tightness == 2) ? 20 : 15;
    }
  else if (tightness <= 5)
    {
      max_rec = 2;
      nmb_div = (tightness == 3) ? 8 : ((tightness == 4) ? 12 : 15);
    }
  else
    {
      max_rec = 3;
      nmb_div = (tightness == 6) ? 10 : ((tightness == 7) ? 12 : 15);
    }

  // Compute trimming seqence
  vector<double> points2(points.rawData(), 
			 points.rawData()+points.dimension()*points.numPoints());
  vector<vector<double> > seqs;
  TrimUtils trimutil(&points2[0], points.numPoints(), 1);
  trimutil.computeTrimSeqs(max_rec, nmb_div, seqs);
  
  double udel, vdel;
  trimutil.getDomainLengths(udel, vdel);
  std::cout << "Minimum sub domains, diag = "<< trimutil.getDomainDiag();
  std::cout << ", udel = " << udel << ", vdel = " << vdel << std::endl;

  std::ofstream of0("translated_seq.g2");
  for (size_t kr=0; kr<seqs.size(); ++kr)
    {
      if (seqs[kr].size() <= 2)
	continue;
      of0 << "410 1 0 0" << std::endl;
      of0 << seqs[kr].size()/2-4 << std::endl;
      size_t ki, kj;
      for (ki=0; ki<seqs[kr].size()-4; ki+=2)
	{
	  for (kj=0; kj<2; ++kj)
	    of0 << seqs[kr][ki+kj] << "  ";
	  of0 << 0 << " ";
	  for (; kj<4; ++kj)
	    of0 << seqs[kr][ki+kj] << "  ";
	  of0 << 0 << std::endl;
	}
      //of0 << seqs[kr][ki] << " " << seqs[kr][ki+1] << " " << 0 << " ";
      //of0 << seqs[kr][0] << " " << seqs[kr][1] << " " << 0 << std::endl;
    }

  // Compute trimming loop
  // First split the outer loop in kinks
  double eps = std::max(udel, vdel);
  double kink_tol = 5e-01; // 0.1 deg => 5.7 degrees.
  vector<vector<double> > split_seqs = 
    TrimCrvUtils::splitTrimPoints(seqs[0], eps, kink_tol);
  
  // Create trimming curves
  const int par_dim = 2;
  const int max_iter = 5;
  vector<shared_ptr<SplineCurve> > par_cvs;
  for (size_t ki = 0; ki < split_seqs.size(); ++ki)
    {
      shared_ptr<SplineCurve> spline_cv_appr_2d
	(TrimCrvUtils::approximateTrimPts(split_seqs[ki], par_dim, eps, 
					  max_iter));
      par_cvs.push_back(spline_cv_appr_2d);
    }

  std::ofstream of01("trimming_cvs.g2");
  for (size_t kr=0; kr<par_cvs.size(); ++kr)
    {
      par_cvs[kr]->writeStandardHeader(of01);
      par_cvs[kr]->write(of01);
    }

  // The curve should be CCW.
  const double int_tol = 1e-06;
  bool loop_is_ccw = LoopUtils::loopIsCCW(par_cvs, eps, int_tol);
  if (!loop_is_ccw)
    {
        MESSAGE("We should change direction of the loop cv!");
	for (size_t ki = 0; ki < par_cvs.size(); ++ki)
	  {
	    par_cvs[ki]->reverseParameterDirection();
	  }
	reverse(par_cvs.begin(), par_cvs.end());
    }

  TrimCrvUtils::moveCurveCoefsInsideSurfDomain(sf.get(), par_cvs);

    bool use_linear_segments = false;
    vector<shared_ptr<ParamCurve> > par_loop;
    if (use_linear_segments)
    {
	for (size_t ki = 0; ki < par_cvs.size(); ++ki)
	{
	    vector<shared_ptr<Line> > line_segments = TrimCrvUtils::approximateCurve(*par_cvs[ki], eps);
	    par_loop.insert(par_loop.end(), line_segments.begin(), line_segments.end());
	}
    }
    else
    {
	par_loop.insert(par_loop.end(), par_cvs.begin(), par_cvs.end());
    }

    vector<shared_ptr<CurveOnSurface> > loop;
    for (size_t ki = 0; ki < par_loop.size(); ++ki)
    {
	shared_ptr<CurveOnSurface> cv_on_sf(new CurveOnSurface(sf, par_loop[ki], true));
	loop.push_back(cv_on_sf);
    }
    const bool fix_trim_cvs = false;
    const double epsgeo_bd_sf = 1e-03;
    BoundedSurface bd_sf(sf, loop, epsgeo_bd_sf, fix_trim_cvs);
    int valid_state = 0;
    bool is_valid = bd_sf.isValid(valid_state);
    if (!is_valid)
    {
	MESSAGE("Created invalid BoundedSurface, valid_state = " << valid_state);
    }
    else
    {
	MESSAGE("Surface is valid!");
     }

    std::ofstream of("translated_trimmed.g2");
    bd_sf.writeStandardHeader(of);
    bd_sf.write(of);

    if (bd_sf.dimension() == 1)
      {
	std::ofstream of2("translated_trimmed_3D.g2");
	shared_ptr<BoundedSurface> tmp_bd(bd_sf.clone());
	shared_ptr<ParamSurface> tmp_sf = tmp_bd->underlyingSurface();
	if (tmp_sf.get())
	  {
	    shared_ptr<LRSplineSurface> tmp_lr = 
	      dynamic_pointer_cast<LRSplineSurface, ParamSurface>(tmp_sf);
	    tmp_lr->to3D();
	    tmp_bd->writeStandardHeader(of2);
	    tmp_bd->write(of2);
	    // std::ofstream of22("translatedtp_trimmed_3D.g2");
	    // bool ok = tmp_bd->makeUnderlyingSpline();
	    // tmp_bd->writeStandardHeader(of22);
	    // tmp_bd->write(of22);
	  }
      }

    // Translate back
    RectDomain dom2 = bd_sf.containingDomain();
    double umin2 = dom2.umin();
    double umax2 = dom2.umax();
    double vmin2 = dom2.vmin();
    double vmax2 = dom2.vmax(); 
    bd_sf.setParameterDomain(umin2-vec[0], umax2-vec[0],
			     vmin2-vec[1], vmax2-vec[1]);
    bd_sf.writeStandardHeader(fileout_bd_sf);
    bd_sf.write(fileout_bd_sf);

  
}
