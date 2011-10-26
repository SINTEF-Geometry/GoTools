//===========================================================================
//                                                                           
// File: SplineInterpolator.h                                              
//                                                                           
// Created: Mon Oct 16 14:59:57 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: SplineInterpolator.h,v 1.16 2009-05-13 07:30:49 vsk Exp $
//                                                                           
// Description: Does spline interpolation.
//                                                                           
//===========================================================================

#ifndef _SPLINEINTERPOLATOR_H
#define _SPLINEINTERPOLATOR_H


#include "GoTools/geometry/Interpolator.h"
#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/utils/Point.h"
#include <memory>


namespace Go
{

    /** An Interpolator that generates a spline curve
     *  interpolating the given dataset.
     */
class GO_API SplineInterpolator : public Interpolator
{
public:
    /// Constructor takes no arguments
    SplineInterpolator()
	: ctype_(None), basis_set_(false)
	{
	}

    /// Virtual destructor ensures safe inheritance
    virtual ~SplineInterpolator();

    // inherited from Interpolator
    virtual const BsplineBasis& basis();

    /// This interpolating function also takes user-specified tangent information 
    /// into account, but does \em not rely on the previously specified end 
    /// conditions.  It constructs its own basis based on the given 
    /// parametrization and order, so neither does it care about any basis
    /// previously set by the user.  It carries out order-1 spline interpolation.
    /// \param params vector containing the parameters for the data points.  Its
    ///               size is equal to the total number of datapoints.
    /// \param points vector containing the coordinates of the data points.  Its
    ///               size is equal to the total number of datapoints multiplied
    ///               with the spatial dimension.  (NB: This is how the algorithm
    ///               deduces the spatial dimension to use: divide the size of
    ///               'points' with the size of 'params').
    /// \param tangent_index a vector containing the indexes of those data points
    ///                      that has tangents associated with them.
    /// \param tangent_points a vector containing the tangents associated with 
    ///                       those data points that are referred to in 'param_index'.
    ///                       Its size is equal to the total number of datapoints that
    ///                       has tangents, multiplied by the dimension of the space.
    /// \param order the order of the spline basis to generate.
    /// \param coefs Upon function completion, this vector will hold the coordinates
    ///              of the control points of the generated spline curve.  (Its basis
    ///              can be obtained by calling the \ref basis() function).
    void interpolate(const std::vector<double>& params,
		     const std::vector<double>& points,
		     const std::vector<int>& tangent_index,
		     const std::vector<double>& tangent_points,
		     int order,
		     std::vector<double>& coefs);

    /// Does the same as the  \ref interpolate() function above, but expects the basis
    /// to have been specified in advance (by \ref makeBasis() or \ref setBasis()).
    /// If the basis is not consistent with the given data set, an exception will be
    /// thrown. (The number of point must at least be equal to the basis' order, and 
    /// the number of basis functions must be equal to the number of points plus the 
    /// number of tangents).
    /// For a description of the parameter list, see \ref interpolate().
    void interpolate(const std::vector<double>& params,
		     const std::vector<double>& points,
		     const std::vector<int>& tangent_index,
		     const std::vector<double>& tangent_points,
		     std::vector<double>& coefs);

    /// The interpolating function, as inherited by \ref Interpolator. 
    /// Does cubic spline interpolation of points (does not care
    /// about tangents).  It constructs its own basis based on the
    /// given parametrization, so it does not care about any basis
    /// previously set by the user.  Previously specified end conditions will be
    /// taken into account.
    /// For parameter list, see \ref Interpolator.
    virtual void interpolate(int num_points, int dimension,
			     const double* param_start,
			     const double* data_start,
			     std::vector<double>& coefs);

    /// CondType Enumerator specifying possible boundary conditions.
    /// \verbatim
    /// None           boundary conditions have not been specified yet. This will lead 
    ///                to an error if the virtual interpolate() function is called. (The
    ///                other two interpolate() functions do not care about end conditions.
    /// Hermite        The tangents are imposed at the start and end of the curve.
    /// Natural        No imposed tangents, but curvature is imposed to be zero at the
    ///                start and end of the curve.
    /// Free           No imposed conditions, neither at start nor end of curve.
    /// NaturalAtStart Curvature imposed to be zero at start of curve.  End of curve
    ///                will either have no conditions at all, or have an imposed
    ///                tangent (depending on whether the  latter has been specified or not.
    /// NaturalAtEnd   Curvature imposed to be zero at end of curve.  Start of curve
    ///                will either have no conditions, or have an imposed tangent
    ///                (depending on whether the latter has been specified or not).
    /// \endverbatim
    enum CondType { None, Hermite, Natural, Free, NaturalAtStart, NaturalAtEnd };

    /// Set the endpoint conditions to 'Hermite'; impose tangents at start and end of curve.
    /// \b Note: endpoint conditions are considered by the \em virtual interpolate()
    ///          function.  The other interpolate functions disregard this setting.
    /// \param start_tangent the imposed tangent at the start of the curve
    /// \param end_tangent the imposed tangent at the end of the curve.
    void setHermiteConditions(const Point& start_tangent,
			      const Point& end_tangent) {
	ctype_ = Hermite;
	start_tangent_ = std::shared_ptr<Point>(new Point(start_tangent));
	end_tangent_ =  std::shared_ptr<Point>(new Point(end_tangent));
    }

    /// Set the endpoint conditions to 'Natural'; zero curvature at start and end of curve.
    /// \b Note: endpoint conditions are considered by the \em virtual interpolate()
    ///          function.  The other interpolate functions disregard this setting.
    void setNaturalConditions()	{ ctype_ = Natural; }

    /// Set the endpoint condition at start of curve to 'Natural'.  (Condition at end of 
    /// curve will be 'Hermite' if a tangent has been previously specified, or 'Free'
    /// otherwise).
    /// \b Note: endpoint conditions are considered by the \em virtual interpolate()
    ///          function.  The other interpolate functions disregard this setting.
    void setNaturalStartCondition() { ctype_ = NaturalAtStart;}

    /// Set the endpoint condition at end of curve to 'Natural'.  (Condition at start
    /// of curve will be 'Hermite' if a tangent has been previously specified, or 'Free'
    /// otherwise).
    /// \b Note: endpoint conditions are considered by the \em virtual interpolate()
    ///          function.  The other interpolate functions disregard this setting.
    void setNaturalEndCondition() { ctype_ = NaturalAtEnd; }

    /// Set the endpoint conditions to 'Free' (meaning no conditions at all).
    /// \b Note: endpoint conditions are considered by the \em virtual interpolate()
    ///          function.  The other interpolate functions disregard this setting.
    void setFreeConditions() { ctype_ = Free; }

    /// Specify tangents for endpoints of curve.  One or both of the shared pointers
    /// may be zero, which means that the corresponding endpoint condition will be set
    /// to 'Natural'; otherwise it will be set to 'Hermite'.
    /// \b Note: endpoint conditions are considered by the \em virtual interpolate()
    ///          function.  The other interpolate functions disregard this setting.
    void setEndTangents(std::shared_ptr<Point>& start_tangent,
			std::shared_ptr<Point>& end_tangent)
    {
	start_tangent_ = start_tangent;
	end_tangent_ = end_tangent;
	if (start_tangent.get() != 0 && end_tangent.get() != 0)
	    ctype_ = Hermite;
	else if (start_tangent.get() != 0)
	    ctype_ = NaturalAtEnd;
	else if (end_tangent.get() != 0)
	    ctype_ = NaturalAtStart;
	else
	    ctype_ = Natural;
    }

    /// Query the type of endpoint conditions currenltly specified.
    /// \return the currently specified CondType.
    CondType getCondType() {return ctype_; }

    /// Given a set of interpolation conditions (currently position and tangent
    /// information) and order, specify a fitting BsplineBasis (stored internally).
    /// \param params vector containing the parametrization of the expected datapoints
    ///               (one parameter per point).
    /// \param tangent_index vector containing indexes of those datapoints where tangent
    ///                      conditions will be imposed as well.
    /// \param order Order of the requested BsplineBasis.
    void makeBasis(const std::vector<double>& params,
		   const std::vector<int>& tangent_index,
		   int order);

    /// WARNING! basis must be suited to interpolation conditions!
    /// If there is no need to end up with a specific basis, makeBasis() is safer.

    /// Specify a user-defined BpslineBasis to be used for interpolation. 
    /// \param basis the specified basis.  NB: it must be suited to the interpolation
    ///              conditions that are going to be used!  If there is no need to
    ///              end up with a specific basis, makeBasis() is safer.
    void setBasis(const BsplineBasis& basis) {
	basis_ = basis;
	basis_set_ = true;
    }

private:
    CondType ctype_;
    std::shared_ptr<Point> start_tangent_;
    std::shared_ptr<Point> end_tangent_;
    BsplineBasis basis_;
    bool basis_set_;
};



} // namespace Go


#endif // _SPLINEINTERPOLATOR_H

