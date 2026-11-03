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

#ifndef DEFINEREFS3D_H
#define DEFINEREFS3D_H

#include "GoTools/lrsplines3D/LRSplineVolume.h"

//==============================================================================
namespace Go
//==============================================================================
{
  /// Given a specifiction for how to refine with respect to an LRBSpline3D
  /// or Element3D, prepare input for LRSplineVolume::refine

  class LRBSpline3D;

  enum RefineStrategy3D
    {
      notSet = 0,
      fullSpan = 1,
      minSpanSize = 2,
      minSpanError = 3,
      minSpanCombined = 4,
      structuredMesh = 5,
      restrictedMesh = 6,
    };
  
  namespace DefineRefs3D
  {
    /// For a set of B-splines, return specifiction of new mesh rectangles to insert
    /// in order to perform structured mesh refinement
    /// \param adjust if true: elongate existing mesh lines if possible
    ///               if false: always insert new mesh lines in the middle of
    ///               existing knots
    /// \param reduced if true: skip some mesh lines if a complete structured mesh
    ///                         would lead to close mesh lines
    void refineStructuredMesh(const LRSplineVolume& vol,
			      std::vector<LRBSpline3D*>& source,
			      std::vector<LRSplineVolume::Refinement3D>& refs_x,
			      std::vector<LRSplineVolume::Refinement3D>& refs_y,
			      std::vector<LRSplineVolume::Refinement3D>& refs_z,
			      bool adjust = true, bool reduced = false);

    /// Add a new mesh rectangle to the collection of such. Combine with
    /// previously defined rectangle if possible
    /// Note that a combination is not unique
    void appendRef(std::vector<LRSplineVolume::Refinement3D>& refs,
		   LRSplineVolume::Refinement3D& curr_ref, double tol);

    void refineStructuredMesh(const LRSplineVolume& vol,
			      LRBSpline3D* bspline, Direction3D fixdir,
			      int mult,
			      std::vector<LRSplineVolume::Refinement3D>& refs);

   void refineFullSpan(const LRSplineVolume& vol,
		       double upar, double vpar, double wpar,
		       Direction3D fixdir, int mult,
		       LRSplineVolume::Refinement3D& refs);

    void refineFromElement(const LRSplineVolume& vol,
			   Element3D *elem, Direction3D fixdir,
			   RefineStrategy3D strategy, int mult,
			   LRSplineVolume::Refinement3D& refs);

    
  } // end namespace DefineRefs3D
} // end namespace Go


#endif
