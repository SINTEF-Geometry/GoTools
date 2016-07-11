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

#ifndef _PARAMVOLUME_H
#define _PARAMVOLUME_H

#include "GoTools/utils/config.h"
#include "GoTools/utils/DirectionCone.h"
#include "GoTools/geometry/GeomObject.h"
#include <memory>
#include <vector>


namespace Go
{

  class ParamSurface;
  class SplineVolume;

  /** Base class for parametric Volumes in Go
   *
   */
class ParamVolume : public GeomObject
{
public:
    /// Virtual destructor, enables safe inheritance.
    virtual ~ParamVolume();

    /// make a clone of this volume and return a pointer to it (user is 
    /// responsible for clearing up memory afterwards).
    /// \return pointer to cloned object
    virtual ParamVolume* clone() const = 0;

    /// Creates a DirectionCone covering all tangents to 
    /// this volume along a given parameter direction.
    /// \param pardir if 1, the DirectionCone will be defined on basis 
    ///        of the volume's tangents along the first parameter direction. If 2,
    ///        the second parameter direction will be used. If 3, the third
    ///        parameter direction will be used.
    /// \return a DirectionCone (not necessarily the smallest) containing all tangents
    ///         to this volume along the specified parameter direction.
    virtual DirectionCone tangentCone(int pardir) const = 0;

    /// Return the parameter domain of the volume.  This is an array containing the start
    /// and end parameter values.  The spline's parameter i has its start value at the
    /// array position (2i) and its end parameter value at the array position (2i+1)
    /// \return An array describing the parametric domain of the volume 
    virtual const Array<double,6> parameterSpan() const = 0;

    /// Evaluates the volume's position for a given parameter triple.
    /// \param pt the result of the evaluation is written here 
    /// \param upar the first parameter
    /// \param vpar the second parameter
    /// \param wpar the third parameter
    virtual void point(Point& pt, double upar, double vpar, double wpar) const = 0;

    /// Evaluates the volume's position and a certain number of derivatives
    /// for a given parameter triple.
    /// \param pts the vector containing the evaluated values. Its size must be set
    ///            by the user prior to calling this function. Upon completion of the
    ///            function, its first entry is the volume's position at the given
    ///            parameter pair. Then, if 'derivs' > 0, the two next entries will
    ///            be the volume tangents along the first, second and third parameter
    ///            direction. The next three entries are the second- and cross
    ///            derivatives, in the order (du2, dudv, dudw, dv2, dvdw, dw2), and
    ///            similar for even higher derivatives.
    /// \param upar the first parameter 
    /// \param vpar the second parameter
    /// \param wpar the third parameter
    /// \param derivs number of requested derivatives
    /// \param u_from_right specify whether derivatives along the first parameter are
    ///                     to be calculated from the right ('true', default) or from
    ///                     the left ('false')
    /// \param v_from_right specify whether derivatives along the second parameter are
    ///                     to be calculated from the right ('true', default) or from
    ///                     the left ('false')
    /// \param w_from_right specify whether derivatives along the third parameter are
    ///                     to be calculated from the right ('true', default) or from
    ///                     the left ('false')
    /// \param resolution tolerance used when determining whether parameters are located 
    ///                   at special values of the parameter domain (in particualar; knot
    ///                   values in case of spline objects.
    virtual void point(std::vector<Point>& pts, 
		       double upar, double vpar, double wpar,
		       int derivs,
		       bool u_from_right = true,
		       bool v_from_right = true,
		       bool w_from_right = true,
		       double resolution = 1.0e-12) const = 0;

    /// Determine the parameter value of the start of the 'next
    /// segment' from a parameter value, along a given parameter
    /// direction.  A 'segment' is here defined as a parameter
    /// interval in which there will be no discontinuities in
    /// derivatives or other artifacts.  For spline objects, a segment
    /// will typically be the interval between two consecutive,
    /// non-coincident knots.
    /// \param dir the parameter direction in which we search for the
    ///            next segment (0, 1 or 2)
    /// \param par the parameter value starting from which we search
    ///            for the start value of the next segment
    /// \param forward define whether we shall move forward ('true')
    ///                or backwards when searching along this parameter
    /// \param tol tolerance used for determining whether the 'par' is
    ///            already located \em on the next segment value
    /// \return the value of the start value of the next segment (or
    ///         the end of the previous segment, if we are moving
    ///         backwards...)
    virtual double nextSegmentVal(int dir, double par, bool forward, double tol) const = 0;

    /// Iterates to the closest point to pt on the volume. 
    /// \param pt the point to find the closest point to
    /// \param clo_u u parameter of the closest point
    /// \param clo_v v parameter of the closest point
    /// \param clo_w w parameter of the closest point
    /// \param clo_pt the geometric position of the closest point
    /// \param clo_dist the distance between pt and clo_pt
    /// \param epsilon parameter tolerance (will in any case not be higher than
    ///                sqrt(machine_precision) x magnitude of solution
    /// \param seed pointer to parameter values where iteration starts.
    virtual void closestPoint(const Point& pt,
			      double&        clo_u,
			      double&        clo_v,
			      double&        clo_w,
			      Point&         clo_pt,
			      double&        clo_dist,
			      double         epsilon,
			      double   *seed = 0) const = 0;

    /// Reverses the direction of the basis in input direction.
    /// \param pardir which parameter direction to reverse
    virtual void reverseParameterDirection(int pardir) = 0;

    /// Swaps two parameter directions
    /// \param pardir1 One of the parameter directions (0, 1 or 2) to be swapped
    /// \param pardir2 The other parameter direction (0, 1 or 2) to be swapped
    virtual void swapParameterDirection(int pardir1, int pardir2) = 0;

    /// Fetch all boundary surfaces corresponding to the volume.
    virtual 
      std::vector<shared_ptr<ParamSurface> > getAllBoundarySurfaces() const = 0;

    /// Generate and return a ParamSurface that represents a constant parameter 
    /// surface on the volume
    /// \param parameter value of the fixed parameter
    /// \param pardir 0 if the surface is constant in the u-parameter,
    ///               1 if the surface is constant in the v-parameter,
    ///               2 if the surface is constant in the w-parameter.
    /// \return pointer to a newly constructed SplineSurface representing the 
    ///         specified constant parameter surface.  It is the user's reponsibility
    ///         to delete it when it is no longer needed.
    virtual ParamSurface* constParamSurface(double parameter,
					    int pardir) const
    {
      // For the time being ...
      return 0;
    }

    /// Check if the volume is of type spline
    virtual bool isSpline() const
    {
      return false;  // Default behaviour, overridden in the spline case
    }

    /// Return the spline volume represented by this volume, if any
    virtual SplineVolume* asSplineVolume() 
    {
      return 0;  // Default behaviour
    }

      /// Translate
    virtual void translate(const Point& vec) = 0;
};


} // namespace Go


#endif // _PARAMVOLUME_H

