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

#ifndef _BLENDINGCURVE_
#define _BLENDINGCURVE_

#include "GoTools/creators/EvalCurve.h"
#include "GoTools/geometry/SplineCurve.h"
#include <vector>
#include <memory>


namespace Go
{

/// This curve represent an offset curve from a given space curve,
/// along a direction obtained by blending two 'cross-tangent curves',
/// and with an offset distance which is a linear function interpolating
/// the cross-tangent length at the start and end of the curve.

class CrossTangentOffset : public EvalCurve
{
 public:

  /// Constructor, taking a curve from which we take the offset, and four other
  /// curves used to calculate the offset direction and magnitude.  Two of these
  /// curves are the cross-tangent curves and the two other are blending functions
  /// (dimension 1).  To find the offset direction at a given point, the two 
  /// cross-tangent curves are evaluated at the specified parameter, multiplied by
  /// their respective blending functions and added together.  The offset length 
  /// is computed by linearly interpolating the length of this blended cross-tangent
  /// at the start and end parameter of the curve.
  /// \param poscurve the curve from which we take the offset
  /// \param tangcv1 the first cross-tangent curve
  /// \param tangcv2 the second cross-tangent curve
  /// \param blend1 the blending function for the first cross-tangent curve
  /// \param blend2 the blending function for the second cross-tangent curve
  CrossTangentOffset(shared_ptr<SplineCurve>& poscurve,
		       shared_ptr<SplineCurve>& tangcv1,
		       shared_ptr<SplineCurve>& tangcv2,
		       shared_ptr<SplineCurve>& blend1,
		       shared_ptr<SplineCurve>& blend2);


  /// virtual destructor enables safe inheritance
  virtual ~CrossTangentOffset();

  // Inherited from EvalCurve
  virtual Point eval( double t) const ;

  // Inherited from EvalCurve
  virtual void eval( double t, int n, Point der[]) const ; // n = order of diff

  // Inherited from EvalCurve
  virtual double start() const ;

  // Inherited from EvalCurve
  virtual double end() const ;

  // Inherited from EvalCurve
  virtual int dim() const ;

  /// Inherited from EvalCurve::approximationOK().  For this class, both tolerances
  /// are used.
  /// \param par the parameter at which to check the curve
  /// \param approxpos the position we want to check whether or not the curve
  ///                  approximates for parameter 'par'.
  /// \param tol1 spatial approximation tolerance.  If the evaluated position is
  ///             outside this tolerance, 'false' is returned.  If it inside this
  ///             tolerance "by far", 'true' is returned.  Otherwise, 'tol2' is taken
  ///             into account.
  /// \param tol2 This tolerance is taken into account when the evaluated point is 
  ///             within 'tol1' from 'approxpos', but not convincingly so.  In that case
  ///             this tolerance is used to check whether the evaluated cross tangent lies
  ///             in the plane spanned by the tangent curves at this point.  (It is thus
  ///             an \em angle tolerance).
  /// \return 'true' if the curve approximates the point at the parameter, 'false'
  ///         otherwise.
  virtual bool approximationOK(double par, Point approxpos,
			       double tol1, double tol2) const ;

 private:
  const shared_ptr<SplineCurve> poscurve_;
  std::vector<shared_ptr<SplineCurve> > tangcurves_;  // will always have 2 elements
  std::vector<shared_ptr<SplineCurve> > blends_; // will always have two elemnets
  shared_ptr<SplineCurve> length_;

  void evalcrtan(double t, int n, Point der[]) const ;
  Point evalcrtan(double t) const ;
};

}
#endif
