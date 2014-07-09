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

#include "GoTools/compositemodel/IntResultsCompCv.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/tesselator/LineStrip.h"
#include "GoTools/tesselator/CurveTesselator.h"

using std::vector;

namespace Go
{
  //===========================================================================
  // Constructor
  IntResultsCompCv::IntResultsCompCv(CompositeCurve* compcv, const ftLine& line)
  //===========================================================================
  : IntResultsModel(CompositeCurve_Line)
  {
    addLineInfo(line);
  }

  //===========================================================================
  // Constructor
  IntResultsCompCv::IntResultsCompCv(CompositeCurve* compcv, const ftPlane& plane)
  //===========================================================================
  : IntResultsModel(CompositeCurve_Plane)
  {
    addPlaneInfo(plane);
  }

  //===========================================================================
  // Destructor
  IntResultsCompCv::~IntResultsCompCv()
  //===========================================================================
  {
  }

  //===========================================================================
  void IntResultsCompCv::addIntPt(shared_ptr<ParamCurve> cv, double* parval)
  //===========================================================================
  {
    if (numpar_ == 1)
      int_pts_1cv_.push_back(PointOnCurve(cv, *parval));
  }


  //===========================================================================
  void IntResultsCompCv::addIntCv(shared_ptr<ParamCurve> cv, double* startpar,
			     double *endpar)
  //===========================================================================
  {
    if (numpar_ == 1)
      int_seg_1cv_.push_back(std::make_pair(PointOnCurve(cv, *startpar),
				       PointOnCurve(cv, *endpar)));
  }

  //===========================================================================
  void 
  IntResultsCompCv::getIntersectionPoints(std::vector<PointOnCurve>& int_points) const
  //===========================================================================
  {
    int_points = int_pts_1cv_;
  }

  //===========================================================================
  void 
  IntResultsCompCv::getIntersectionCurves(std::vector<std::pair<PointOnCurve, PointOnCurve> >& int_crvs) const
  //===========================================================================
  {
    int_crvs = int_seg_1cv_;
  }

  //===========================================================================
  void IntResultsCompCv::tesselate(std::vector<shared_ptr<LineStrip> >& meshes,
				   PointCloud3D& points) const
  //===========================================================================
  {
    int res = 100;
    tesselate(res, meshes, points);
  }

  //===========================================================================
  void IntResultsCompCv::tesselate(int resolution,
				   std::vector<shared_ptr<LineStrip> >& meshes,
				   PointCloud3D& points) const
  //===========================================================================
  {
    if (hasIntCurves())
      {
	for (size_t ki=0; ki<int_seg_1cv_.size(); ++ki)
	  {
	    double t1 = int_seg_1cv_[ki].first.getPar();
	    double t2 = int_seg_1cv_[ki].second.getPar();
	    shared_ptr<ParamCurve> cv = int_seg_1cv_[ki].first.getCurve();
	    shared_ptr<ParamCurve> sub_cv = 
	      shared_ptr<ParamCurve>(cv->subCurve(std::min(t1,t2),
						  std::max(t1,t2)));
	    
	    CurveTesselator tesselator(*sub_cv.get());
	    tesselator.changeRes(resolution);
	    shared_ptr<LineStrip> mesh = tesselator.getMesh();
	    meshes.push_back(mesh);
	  }
      }

    vector<double> coords;
    for (size_t ki=0; ki<int_pts_1cv_.size(); ++ki)
      {
	Point pt = int_pts_1cv_[ki].getPos();
	coords.insert(coords.end(), pt.begin(), pt.end());
      }

    points = PointCloud3D(&coords[0], (int)int_pts_1cv_.size());
  }

  //===========================================================================
  void IntResultsCompCv::tesselate(double density,
				   std::vector<shared_ptr<LineStrip> >& meshes,
				   PointCloud3D& points) const
  //===========================================================================
  {
    int min_nmb = 5;
    int max_nmb = (int)(1000000.0/(int)int_seg_1cv_.size());
    if (hasIntCurves())
      {
	for (size_t ki=0; ki<int_seg_1cv_.size(); ++ki)
	  {
	    double t1 = int_seg_1cv_[ki].first.getPar();
	    double t2 = int_seg_1cv_[ki].second.getPar();
	    shared_ptr<ParamCurve> cv = int_seg_1cv_[ki].first.getCurve();
	    shared_ptr<ParamCurve> sub_cv = 
	      shared_ptr<ParamCurve>(cv->subCurve(std::min(t1,t2),
						  std::max(t1,t2)));
	    double len = sub_cv->estimatedCurveLength();
	    int res = (int)(len/density);
	    res = std::max(min_nmb, std::min(res, max_nmb));
	    
	    CurveTesselator tesselator(*sub_cv.get());
	    tesselator.changeRes(res);
	    shared_ptr<LineStrip> mesh = tesselator.getMesh();
	    meshes.push_back(mesh);
	  }
      }

    vector<double> coords;
    for (size_t ki=0; ki<int_pts_1cv_.size(); ++ki)
      {
	Point pt = int_pts_1cv_[ki].getPos();
	coords.insert(coords.end(), pt.begin(), pt.end());
      }

    points = PointCloud3D(&coords[0], (int)int_pts_1cv_.size());
  }

} // namespace Go
