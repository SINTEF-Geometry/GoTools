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

#ifndef __SURFACEMODELUTILS_H
#define __SURFACEMODELUTILS_H

#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/compositemodel/Vertex.h"

namespace Go
{
  /// Utility functionality for surface models/surface collections. 
  /// Intended for internal use in operations on surface models

  class ftSurface;
  class SurfaceModel;
  class Body;
  class BoundedSurface;
  class CurveOnSurface;

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

    /// Check if the surfaces corresponding to face1 and fac2 can be merged into
    /// one larger surface
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
			     std::vector<bool>& at_bd1,
			     Body *model1,
			     std::vector<std::vector<shared_ptr<CurveOnSurface> > >& cvs2,
			     std::vector<shared_ptr<ParamSurface> >& sfs2,
			     std::vector<bool>& at_bd2,
			     Body *model2, double eps, double angtol,
			     std::vector<std::vector<shared_ptr<ParamSurface> > >& groups);
  }
}
#endif
