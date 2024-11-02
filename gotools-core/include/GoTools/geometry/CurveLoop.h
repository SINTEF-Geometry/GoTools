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

#ifndef _CURVELOOP_H
#define _CURVELOOP_H

// #include <memory>
#include <memory>
#include <vector>
#include "GoTools/utils/config.h"
#include "GoTools/geometry/ParamCurve.h"


namespace Go
{
//===========================================================================
/// Computes the largest gap in the loop specified by the vector of curves
template <class PtrToCurveType>
inline double computeLoopGap(const std::vector< PtrToCurveType >& curves)
//===========================================================================
{

    // Here, we should check that the given curves indeed are forming a
    // loop, so every endpoint is within space_epsilon of the start of the
    // next curve.
    // Also, we make sure that all curves have the same dimension
    // and that the curves are of the same type
    if ((curves.size() == 0) || (curves[0].get() == NULL))
        return -1.0;

    int dim = curves[0]->dimension();
//     ClassType type = curves[0]->instanceType();
    int n = (int)curves.size();
    int i;
    for (i = 1; i < n; ++i) {
	if (curves[i].get() == NULL)
	    return -1.0;
	if (curves[i]->dimension() != dim) {
	    THROW("Curves do not have the same dimension.");
	}
// 	if (curves[i]->instanceType() != type) {
// 	    THROW("Not all curves are of the same type.");
// 	}
    }
    Point startp(dim);
    Point endp(dim);
    double maxdist = -1.0;
    double dist;
    for (i = 1; i < n; ++i) {
	curves[i-1]->point(endp, curves[i-1]->endparam());
	curves[i]->point(startp, curves[i]->startparam());
	dist = endp.dist(startp);
	if (dist > maxdist) maxdist = dist;
    }
    curves[n-1]->point(endp, curves[n-1]->endparam());
    curves[0]->point(startp, curves[0]->startparam());
    dist = endp.dist(startp);
    if (dist > maxdist) maxdist = dist;
    return maxdist;
}


    /** CurveLoop represents a closed loop defined by the composition of 
     *  a set of curves.  The start point of curve (i+1) should coincide
     *  (within a certain tolerance) with the end point of curve (i), and the
     *  end point of the last curve should coincide with the start point of
     *  the first curve.  CurveLoops are useful in, for example, describing
     *  the boundary of a surface.
     */

 class CurveOnSurface;


class GO_API CurveLoop
{
public:
    /// Create an empty loop
    CurveLoop();

    /// Create a loop based on a given set of curves. 
    /// \param curves a vector with shared pointers to the curves defining the
    ///        CurveLoop.  The start point of each curve should match the 
    ///        end point of the previous curve within the given tolerance. 
    ///        The same is the same for the start point on the first curve and
    ///        the end point on the last curve).  Moreover, it is expected that
    ///        all curves are of the same type, and that they all lie in the 
    ///        same space (ie. have the same dimension).
    /// \param space_epsilon the given tolerance for defining coincidence between
    ///                      start/end points on curves.
    CurveLoop(const std::vector< shared_ptr<ParamCurve> >& curves,
	      double space_epsilon, bool allow_fix=true);

    CurveLoop(const std::vector< shared_ptr<CurveOnSurface> >& curves,
	      double space_epsilon, bool allow_fix=true);

    /// Virtual destructor allows safe inheritance
    virtual ~CurveLoop();

    /// Quick swap of one CurveLoop with another.
    void swap(CurveLoop& other);
    
    /// Reset the CurveLoop based on a given vector of curves. (But keep the 
    /// previously set tolerance value.
    /// \param curves a vector with shared pointers to the curves defining the
    ///        CurveLoop.  The start point of each curve should match the 
    ///        end point of the previous curve within the given tolerance. 
    ///        The same is the same for the start point on the first curve and
    ///        the end point on the last curve).  Moreover, it is expected that
    ///        all curves are of the same type, and that they all lie in the 
    ///        same space (ie. have the same dimension).
    void setCurves(const std::vector< shared_ptr<ParamCurve> >& curves,
		   bool allow_fix=true);

    /// Reverse the direction of all curves and their mutual order.
    void turnOrientation();

    /// Set the tolerance (used to determine whether the start/end points on curves
    /// are coincident) to a given value
    /// \param space_epsilon set the tolerance to this value
    void setSpaceEpsilon(const double space_epsilon);

    /// Get the tolerance value (used to determine whether the start/end points on 
    /// curves are coincident).
    double getSpaceEpsilon() const;

    /// Query the number of curves constituting the CurveLoop
    /// \return the number of curves constituting the CurveLoop
    int size() const { return (int)curves_.size(); }

    /// Dimension
    int dimension()
    {
      return curves_.size() == 0 ? 0 : curves_[0]->dimension();
    }

    /// Get a const iterator to the first curve in the CurveLoop
    /// \return a const iterator to the first curve in the CurveLoop
    std::vector< shared_ptr<ParamCurve> >::const_iterator begin() const
    { return curves_.begin(); }
    
    /// Get a const iterator to one-past-the-last curve in the CurveLoop
    /// \return a const iterator to the one-past-the-last curve in the CurveLoop
    std::vector< shared_ptr<ParamCurve> >::const_iterator end() const
    { return curves_.end(); }

    /// Get a iterator to the first curve in the CurveLoop
    /// \return a iterator to the first curve in the CurveLoop
    std::vector< shared_ptr<ParamCurve> >::iterator begin()
    { return curves_.begin(); }

    /// Get a iterator to one-past-the-last curve in the CurveLoop
    /// \return a iterator to the one-past-the-last curve in the CurveLoop
    std::vector< shared_ptr<ParamCurve> >::iterator end()
    { return curves_.end(); }

    /// Get a shared pointer to the i'th curve in the CurveLoop
    /// \param index the index of the requested curve
    /// \return a shared pointer to the requested curve.
    shared_ptr<ParamCurve> operator[] (int index) const;

    /// Find the closest point on the curve loop to a point specified 
    /// by the user.
    /// \param pt The point given by the user.  We want to determine the closest
    ///           point to this on the CurveLoop.
    /// \param clo_ind Upon return: the index of the curve segment on which the 
    ///                closest point was found.
    /// \param clo_par Upon return: the parameter of the detected closest point on 
    ///                the curve containing it.
    /// \param clo_pt  Upon return: the geometric position of the detected closest 
    ///                point
    /// \param clo_dist Upon return: the distance to the detected closest point.
    void closestPoint(const Point& pt, int& clo_ind, double& clo_par, 
		      Point& clo_pt, double& clo_dist) const;

    /// View the loop as a curve in the parameter space and compute the
    /// closest point in the loop.
    /// This function is only interesting if the curves constituting the CurveLoop
    /// are of type "CurveOnSurface".  In that case, the curves have a 2D 
    /// representation in the parametric domain of the surface, as well as a 3D 
    /// representation in geometric space.  This function will search for the closest
    /// point in \em parametric space to a \em parametric (2D) point specified 
    /// by the user.
    /// \param pt The point given by the user.  It should be a 2D point refering to 
    ///           the parametric domain of a surface.  We want to determine the closest
    ///           point to this on the CurveLoop.
    /// \param clo_ind Upon return: the index of the curve segment on which the closest
    ///                point was found.
    /// \param clo_par Upon return: the parameter of the detected closest point on the
    ///                curve segment containing it (a scalar value)
    /// \param clo_pt Upon return: the curve segment's closest point represented as a 
    ///               2D parameter pair in the domain of the underlying surface.
    /// \param clo_dist Upon return: the distance between the point given by the
    ///                 user and the found closest point, measured in the parametric 
    ///                 domain of the underlying surface.
    void closestParPoint(const Point& pt, int& clo_ind, double& clo_par, 
			 Point& clo_pt, double& clo_dist) const;

    /// Return the curves in the loop as a vector 
    std::vector<shared_ptr<ParamCurve> > getCurves()
    { return curves_; }

    /// Return joint points between curves
    std::vector<Point> getCorners() const;

    /// Collect curves with smooth connections
    void getSmoothCurves(std::vector<std::vector<shared_ptr<ParamCurve> > >& curves,
			 double angtol);

    /// If the curves do not form a simple loop within the input
    /// tolerance, the object is invalid.
    bool isValid() const;

    /// If loop is invalid, we try to fix loop. If we fail we alter
    /// nothing. Return value is true if loop is valid. max_gap is
    /// set to the largest gap in the loop.
    bool fixInvalidLoop(double& max_gap);

//     double maxGap(int nmb_seg_samples);

    /// Simplify this curve loop if possible by reducing the number of curves
    bool simplify(double tol, double ang_tol, double& max_dist);

    /// Split specified curve at a given parameter and update corresponding
    /// loop. Return new curves
    std::vector<shared_ptr<ParamCurve> > split(int idx, double par);

    /// Maximum distance between subsequent curves
    double getMaxCurveDist() const
    { return computeLoopGap(curves_); }

    /// Remove curve and tighten gap if possible
    /// Return value: 0 - curve not found
    ///               1 - curve removed, adjacent curves not changed
    ///               2 - curve removed, adjacent curves updated
    int removeCrvAndFix(shared_ptr<CurveOnSurface> cv);

    /// If both parameter and space curve are given for a segment, and
    /// they do not match, one of them recreated
    /// Missing curves are created
    void fixMismatchCurves(double eps);

    void analyze(); // Setting the valid_state_ for the loop.

private:
    std::vector< shared_ptr<ParamCurve> > curves_;
    double space_epsilon_;
    // The curves should form a loop.
    int valid_state_; //  0 = not validated.
                      //  1 = valid
                      // -1 = not valid (the segments do not form a loop)
                      // -2 = not valid (the CurveOnSurface par segments do not form a loop)


};




} // namespace Go

#endif // _CURVELOOP_H

