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

#ifndef _GAPREMOVAL_H
#define _GAPREMOVAL_H

#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/creators/ConstraintDefinitions.h"

namespace Go {

  /// Functionality for removal of gaps between two adjacent surfaces of various types

namespace GapRemoval
{

  /// Both surfaces are non-trimmed spline surfaces and the
  /// common boundary are described by CurveOnSurface entities
  /// Remove gap by averaging coefficients at the common boundary.
  /// May increase the associated spline spaces.
  void
  removeGapSpline(shared_ptr<SplineSurface>& srf1, 
		  shared_ptr<CurveOnSurface>& bd_cv1,
		  double start1, double end1,
		  shared_ptr<SplineSurface>& srf2, 
		  shared_ptr<CurveOnSurface>& bd_cv2,
		  double start2, double end2, Point vertex1, Point vertex2,
		  double epsge, bool *same_orientation = NULL);

  /// Both surfaces are trimmed of some type and the
  /// common boundary are described by CurveOnSurface entities
  double
  removeGapTrim(shared_ptr<CurveOnSurface>& bd_cv1,
		double start1, double end1,
		shared_ptr<CurveOnSurface>& bd_cv2,
		double start2, double end2, Point vertex1, Point vertex2,
		double epsge);

  /// Both surfaces are trimmed and have an undarlying spline surface. The
  /// common boundary are described by CurveOnSurface entities limited
  /// by bounding parameters. Sub curves of the same curve on surface entity 
  /// may occur several times, but are split to have correspondance between
  /// the two surfaces along the commonboundary.
  bool
  removeGapSplineTrim(shared_ptr<SplineSurface>& srf1, 
		      std::vector<shared_ptr<CurveOnSurface> >& bd_cv1,
		      std::vector<double> start1, 
		      std::vector<double> end1,
		      std::vector<shared_ptr<CurveOnSurface> >& bd_cv2,
		      std::vector<double> start2, 
		      std::vector<double> end2, Point vertex1, 
		      Point vertex2, double epsge);

  /// Modify the underlying spline surface of srf1 to improve the adaption
  /// to given trimming curves. 
  void modifySplineSf(shared_ptr<ParamSurface>& srf1, 
		      std::vector<shared_ptr<CurveOnSurface> >& bd_cv1,
		      std::vector<double> start1, std::vector<double> end1,
		      shared_ptr<ParamSurface>& srf2, 
		      std::vector<shared_ptr<CurveOnSurface> >& bd_cv2,
		      std::vector<double> start2, std::vector<double> end2, 
		      double epsge);

  /// Modify the underlying spline surface of srf1 to improve the adaption
  /// to given trimming curves.
  void modifySplines(shared_ptr<ParamSurface>& srf1, 
		     std::vector<shared_ptr<CurveOnSurface> >& bd_cv1,
		     std::vector<double>& start1, std::vector<double>& end1,
		     shared_ptr<ParamSurface>& srf2, 
		     std::vector<shared_ptr<CurveOnSurface> >& bd_cv2,
		     std::vector<double>& start2, std::vector<double>& end2, 
		     std::vector<Point>& vertex, double epsge);
  /// Both surfaces are non-trimmed spline surfaces and the
  /// common boundary are described by CurveOnSurface entities
  /// Remove gap by defining conditions between  coefficients at the commont
  /// boundary between the two surface and modify both surfaces accordingly
  void
    removeGapSpline2(std::vector<shared_ptr<CurveOnSurface> >& bd_cv1,
		     std::vector<double>& start1, std::vector<double>& end1,
		     std::vector<shared_ptr<CurveOnSurface> >& bd_cv2,
		     std::vector<double>& start2, std::vector<double>& end2, 
		     std::vector<Point>& vertex, double epsge);

  /// Define linear constraints between coefficients on two curves entities
  /// representing the same trace to make the trace of the two curves identical.
  std::vector<shared_ptr<sideConstraintSet> > 
    getCoefConstraints(shared_ptr<SplineCurve>& crv1, int idx1,
		       shared_ptr<SplineCurve>& crv2, int idx2, 
		       double tol);

  /// Replace a sub curve (between par1 and par2) in a given spline curve by 
  /// a modified curve piece. The continuity at the joints are given by cont1 
  /// and cont2
  shared_ptr<SplineCurve>
    replaceCurvePiece(shared_ptr<SplineCurve> crv,
		      shared_ptr<SplineCurve> sub_crv,
		      double par1, int cont1, 
		      double par2, int cont2);

  /// Fetch the underlying spline surface and trimming curves from a
  /// bounded surface
  shared_ptr<SplineSurface> 
    getSplineAndBd(shared_ptr<ParamSurface> psurf,
		   std::vector<shared_ptr<CurveOnSurface> >& bd_crvs);

  /// Given the trimming curves of a bounded surface. Sort the curves followin
  /// boundaries of the underlying surface into the vector bd_cvs. The remaining
  /// curves are sampled and the corresponding points and parameter values are
  /// returned in the vectors pts and pars.
  void 
    getBoundarySamples(std::vector<shared_ptr<CurveOnSurface> >& all_bd,
		       std::vector<shared_ptr<CurveOnSurface> >& bd_cvs,
		       std::vector<double>& pts, std::vector<double>& pars,
		       std::vector<int>& bd_idx, double epsge);

  /// Modify a surface to adapt to a given point, vertex, in the parameter value
  /// face_param.
  bool modifyAtVertex(shared_ptr<SplineSurface> srf,
		      Point face_param, Point vertex,
		      double epsge);

  /// Compute the distance between trimming curves on adjacent surfaces
  /// along a common boundary.
  /// \param bd1 trimming curve on first surface
  /// \param bd2 trimming curve on second surface
  /// \param start1 start parameter to the relevant piece of bd1
  /// \param end1 end parameter to the relevant piece of bd1
  /// \param start2 start parameter to the relevant piece of bd2
  /// \param end2 end parameter to the relevant piece of bd2
  /// \param mdist1 maximum distance between given trimming curves
  /// \param mdist2 maximum distance between the trimming curves projected
  /// onto the associated surface
void 
  checkBoundaryDist(shared_ptr<CurveOnSurface> bd1,
		    shared_ptr<CurveOnSurface> bd2,
		    double start1, double end1,
		    double start2, double end2,
		    int nmb_sample, double& mdist1,
		    double& mdist2);

} // of namespace CreatorsUtils.

}; // end namespace Go



#endif // _GAPREMOVAL_H
