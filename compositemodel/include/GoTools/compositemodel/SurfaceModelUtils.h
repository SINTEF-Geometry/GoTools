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

#ifndef _SURFACEMODELUTILS_H
#define _SURFACEMODELUTILS_H

#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/compositemodel/Vertex.h"
#include "GoTools/topology/tpTolerances.h"
#include "GoTools/tesselator/GeneralMesh.h"
#include "GoTools/geometry/Line.h"

namespace Go
{
  /// Utility functionality for surface models/surface collections. 
  /// Intended for internal use in operations on surface models

  class ftSurface;
  class SurfaceModel;
  class Body;
  class BoundedSurface;
  class CurveOnSurface;
  class ftPointSet;

  namespace SurfaceModelUtils
  {
    /// Check if the surface may be closed. In that case split it
    /// into non-closed pieces
    std::vector<shared_ptr<ParamSurface> > 
      checkClosedFaces(shared_ptr<ParamSurface> surface, double tol);

    /// Extract faces that share the same underlying surface
    /// Either the underlying surface is exactly the same, or the faces
    /// relate to the same elementary surface
    void sameUnderlyingSurf(std::vector<shared_ptr<ftSurface> >& sf_set,
			    double tol, double angtol,
			    std::vector<std::vector<shared_ptr<ftSurface> > >& faces,
			    std::vector<shared_ptr<ParamSurface> >& under_sfs);

    /// Extend an underlying surface to be able to serve as a support for all
    /// surfaces in sf_set.
    /// NB! Only supporting elementary surfaces
    shared_ptr<ParamSurface>
      extendedUnderlyingSurface(std::vector<shared_ptr<ftSurface> >& sf_set,
				double tol, double angtol);

    /// Merge surfaces into larger one when possible
    void simplifySurfaceModel(shared_ptr<SurfaceModel>& model, int degree);

  /// Simplify current shell by merging surfaces with a smooth connection
  /// Will produce trimmed surfaces
    void simplifySurfaceModel2(shared_ptr<SurfaceModel>& model, int degree,
			       bool remove_joints=true);

    /// Check if the surfaces corresponding to face1 and face2 can be merged 
    /// into one larger surface
    int mergeSituation(ftSurface* face1, ftSurface* face2,
		       shared_ptr<Vertex> vx1, shared_ptr<Vertex> vx2,
		       int& dir1, double& val1, bool& atstart1, 
		       int& dir2, double& val2, bool& atstart2, 
		       std::pair<Point, Point>& co_par1, 
		       std::pair<Point, Point>& co_par2, double eps);
    
    /// Fetch candidate faces for merge
    std::vector<ftSurface*> getMergeCandFaces(shared_ptr<ftSurface> curr,
					      std::vector<std::pair<shared_ptr<Vertex>,
					      shared_ptr<Vertex> > >& common_vxs,
					      double angtol);

    /// Estimate the size of a possible merged surface from the surfaces corresponding
    /// to face1 and face2
    void estMergedSfSize(ftSurface* face1, ftSurface* face2,
			 shared_ptr<Vertex> vx1,shared_ptr<Vertex> vx2,
			 double& len_frac, double& other_frac, double& sf_reg,
			 double neighbour, double bend);

    /// Split surface model surfaces according to intersections between the models
    /// and sort according to inside/outside of both models
    void sortTrimmedSurfaces(std::vector<std::vector<shared_ptr<CurveOnSurface> > >& cvs1,
			     std::vector<shared_ptr<ParamSurface> >& sfs1,
			     std::vector<int>& at_bd1,
			     Body *model1,
			     std::vector<std::vector<shared_ptr<CurveOnSurface> > >& cvs2,
			     std::vector<shared_ptr<ParamSurface> >& sfs2,
			     std::vector<int>& at_bd2,
			     Body *model2, double eps, double angtol,
			     std::vector<std::vector<std::pair<shared_ptr<ParamSurface>, int> > >& groups,
			     SurfaceModel *shell1=NULL, SurfaceModel *shell2=0);

    /// Intersect surface with a line
    void
      intersectLine(shared_ptr<ParamSurface>& surface,
		    Point pnt, Point dir, double tol,
		    std::vector<std::pair<Point,Point> >& result,
		    std::vector<std::pair<shared_ptr<ParamCurve>, 
		    shared_ptr<ParamCurve> > >& line_seg);

    /// Compute extremal point in a given direction
    bool extremalPoint(shared_ptr<ParamSurface>& surface, Point dir, 
		       tpTolerances& toptol, Point& ext_pnt, 
		       double ext_par[]);

    void setResolutionFromDensity(shared_ptr<ParamSurface> surf,
				  double density, double tol2d,
				  int min_nmb, int max_nmb,
				  int& u_res, int& v_res);

    void tesselateOneSrf(shared_ptr<ParamSurface> surf,
			 shared_ptr<GeneralMesh>& mesh,
			 double tol2d, int n=20, int m=20);

    void triangulateFaces(std::vector<shared_ptr<ftSurface> >& faces,
			  shared_ptr<ftPointSet>& triang, double tol);

    void triangulateModel(shared_ptr<SurfaceModel>& model, double density,
			  shared_ptr<ftPointSet>& triang);

    void samplePointsModel(shared_ptr<SurfaceModel>& model,
			   double density,
			   shared_ptr<ftPointSet>& triang);
    
    void getFaceInnerSamplePoints(shared_ptr<ftSurface>& face, 
				  double density,
				  shared_ptr<ftPointSet>& points);

    void getBoundarySamplePoints(shared_ptr<SurfaceModel>& model, 
				 double density,
				 shared_ptr<ftPointSet>& points);
    
    void performModelTriangulate(shared_ptr<SurfaceModel>& model, 
				 shared_ptr<ftPointSet>& triang,
				 double density);
    void 
      reduceUnderlyingSurface(shared_ptr<BoundedSurface>& bd_sf,
			      std::vector<shared_ptr<CurveOnSurface> >& cvs);

  }
}
#endif
