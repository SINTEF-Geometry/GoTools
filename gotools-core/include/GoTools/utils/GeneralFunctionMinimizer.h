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

#ifndef _GENERAL_FUNCTIONMINIMIZER_H
#define _GENERAL_FUNCTIONMINIMIZER_H

#include "GoTools/utils/Point.h"
#include <vector>
#include <limits>

namespace Go {

template<class Functor> class FunctionMinimizer;


/// This is the algorithm for minimising a function taking multiple parameters, using the
/// conjugated gradient method.  It is used in conjunction with the \ref FunctionMinimizer 
/// class.  Documentation for both can be found here: \ref FunctionMinimizer
//===========================================================================
template<class Functor>
void minimise_conjugated_gradient(FunctionMinimizer<Functor>& dfmin);
//===========================================================================

/// This is the FunctionMinimizer class that can be used ex. with
///  \ref minimise_conjugated_gradient.  Together, they provide functionality for
/// minimizing a function with an arbitrary number of parameters, using linear minimization
/// along a set of directions which is chosen by the conjugated gradient algorithm, although
/// the FunctionMinimizer can also be used independently for minimizing a function along any
/// user-defined direction.  The function to minimize must be presented in the form of a
/// Functor, which has the following member functions:
///
/// \verbatim
/// double operator()(const double* arg) const; // to evaluate the function at the parameters
///                                             // pointed to by 'arg'.
/// void grad(const double* arg, double* grad) const; // to evaluate the function
///                                                   // gradient at the parameters
///                                                   // pointed to by 'arg'.  The resulting
///                                                   // vector will be written to 'grad'.
/// double minPar(int n) const; // return the lower bound of the n'th function parameter.
/// double maxPar(int n) const; // return the upper bound of the n'th function parameter.
/// \endverbatim
/// 
/// Minimizing a function whose domain is R^n is a two-fold problem.  One is of linear 
/// minimization along a specified direction.  Another is to choose a set of directions to 
/// minimize along. The first problem is taken care of by the class FunctionMinimizer,
/// which wraps around the user-defined function and contains information about the current
/// point of evaluation and minimizing along a given direction.  The second problem can 
/// be taken care of by the supplied function 'minimise_conjugated_gradient', or the user
/// can use FunctionMinimizer in combination with his or her own direction-choosing algorithm.
///
/// The user first should wrap his or hers Functor into a FunctionMinimizer, where the
/// starting seed and tolerances should also be specified.  The FunctionMinimizer can be
/// used separately to minimize along directions that the user specify (using the 'minimize'
/// member function), or alternatively can be given to the supplied nonmember function
/// 'minimise_conjugated_gradient', which will try to find a local minimum starting from
/// the supplied seed.  
/// 
/// Example:
/// Let us assume that the user has a function taking four arguments.  He defines the 
/// functor MyFunction containing the required member functions, and instanciates it:
/// \verbatim 
/// MyFunction f(...constructor arguments...)
/// \endverbatim
/// He now wraps it in a FunctionMinimizer, giving the seed where he wants to start the
/// search:
/// \verbatim
/// FunctionMinimizer<MyFunction> fmini(4, f, seed_p, 1.0e-8);
/// \endverbatim
/// This function minimizer can now be used with the 'minimise_conjugated_gradient'-algorithm
/// in order to search for a (local) minimum:
/// \verbatim
/// minimise_conjugated_gradient(fmini, 5);
/// \endverbatim
/// The obtained 4-tuple of parameters corresponding to the found minimum can be obtained by:
/// \verbatim
/// const double* par_ptr = fmini.getPar();
/// \endverbatim
template<class Functor>
class FunctionMinimizer
{
 public:
    /// Constructor. 'minimization_tol' should in general be kept no lower than the 
    /// square root of machine precision. 'seed' is a pointer to an array containing the
    /// start values for the search in the 'num_param' parameters.
    FunctionMinimizer(int num_param, const Functor& fun, 
		      const double* const seed, 
		      double tol = std::sqrt(std::numeric_limits<double>::epsilon()));

    ~FunctionMinimizer() 
    {}
    
    /// Move the current parameter point in order to minimize the
    /// function along an unoriented line. The line is given by the
    /// current parameter point and the parameter 'dir'.
    /// \param dir a direction parallel to the search line
    /// \retval hit_domain_edge indicates whether the found minimum is
    /// on the boundary of the domain.
    double minimize(const Point& dir, bool& hit_domain_edge, bool rerun = false);


    /// Move the current parameter point in order to minimize the
    /// function along a ray. The ray is given by the
    /// current parameter point and the parameter 'dir'.
    /// In contrast to minimize(), linminBrent() requires that the
    /// minimum is bracketed. The method uses an algorithm inspired
    /// by Brent's method (with certain modifications). This method is
    /// rapid when the function is well-behaved. Otherwise it will use
    /// a slow but robust 'golden mean'-search. The current parameter
    /// point is moved to the minimum found along the given
    /// ray. Before calling this function, the minimum must be
    /// bracketed, and the brackets must be given to the function
    /// through the 'brackets' variable, and the corresponding
    /// function values through the 'fval_brak' variable.  The
    /// function value at the new minimum is returned.
    double linminBrent(const Point& dir,
		       const double* bracket,
		       const double* fval_brak);


    /// evaluating distance function at the current point of evaluation
    inline double fval() const; 

    /// evaluating the distance function for an arbitrary parameter value (not the 
    /// parameter pair).
    inline double fval(const Point& param) const;

    /// evaluate the gradient (2 components) of the distance function at the currently set 
    /// parameter pair of evaluation.  
    inline void grad(Point& result) const;

    /// evaluate the gradient (2 components) of the distance function for an arbitrary
    /// parameter value (not the current parameter pair).  
    inline void grad(const Point& param, Point& result) const;

    /// movie the current parameter pair a certain distance ('multiplier') in a certain
    /// direction ('dir', which has 2 components)
    inline void moveUV(const Point& dir, double multiplier);

    /// Returns 'true' if the current parameter for the first/second curve is equal to the 
    /// start/end-parameter for the curve.
    bool atMin(int param_ix) const {return at_min_[param_ix];}
    bool atMax(int param_ix) const {return at_max_[param_ix];}

    /// Get the current value for parameter nb. 'param_ix'.
    double getPar(int param_ix) const {return par_[param_ix];}

    /// Get a pointer to the array containing the current values of the parameters.
    const double* getPar() const {return par_.begin();}

    /// Get the number of function parameters.
    int numPars() const {return par_.size();}
			      
 private:
    const Functor fun_;

    Point par_;
    std::vector<bool> at_min_;
    std::vector<bool> at_max_;

    mutable bool cached_value_updated_;
    mutable bool cached_grad_updated_;

    mutable double cached_value_;
    mutable Point cached_grad_;

    const double minimization_tol_;
    std::vector<double> param_tol_;

    static const double root_machine_precision_;
    static const double rel_tol_;
    static const double perturbation_;
    static const double golden_ratio_;
    static const double default_partition_;
    static const double tiny_;
    static const int max_iter_;


    // update at_umin_, at_umax_, at_vmin_ and at_vmax_.  Should be called whenever the
    // current parameters ('cur_uv_') are changed.
    inline void checkBorder();

    // clear cached values.  Should be called whenever the current parameters ('cur_uv_') 
    // are changed.
    inline void resetCache() const;

    // determine how far we are allowed to step from the current parameter point in a certain
    // direction ('dir') before going out of the valid parameter domain.
    double determineMaxSteplength(const Point& dir) const;

    // Bracket the function minimum from the current parameter point along a certain 
    // direction ('dir'). The maximum allowed steplength is given in 'max_step'.  
    // Upon completion, two things might happen:
    // Possibility 1 : the function was able to bracket a minimum.  The position of the 
    //                 bracketing points (from current parameters, along 'dir') is 
    //                 returned in 'bracket' (3 values), and the corresponding distance
    //                 function values are returned in fval_brak.  The function returns 'true'.
    // Possibility 2 : the function hit the 'max_step' before it was able to bracket any
    //                 minimum.  The function returns 'false'.
    bool bracketInterval(const Point& dir,  
			 const double max_step,
			 //const double working_tolerance,
			 double* bracket,    // points to 3-array
			 double* fval_brak); // points to 3-array

    // given a direction 'dir', this function generates a direction 'result' that is a scaled
    // version of 'dir', with a possible sign change to make sure that the direction is 
    // "downhill" (negative scalar product with function gradient).  Returns 'true' on success,
    // and 'false' if the direction was practically perpendicular to the function gradient.
    bool orientDirection(const Point& dir, Point& result);

    // Returns -1 if the scalar product of p1 and p2 is negative, +1 if the scalar product
    // is positive and 0 if it close to 0.  In this case, 'close to 0' means that it is 
    // impossible to judge the derivatives sign when taking into account the numerical
    // noise associated with the magnitude of the points' components.
    int scalarProductSign(const Point& p1, const Point& p2);

    // Estimate the minimum of a parabola where two points and one tangent is known.  The 
    // two points are the current parameter point and the point obtained by marching from this
    // point a distance of 'dir' multiplied by 'p2'.  The function value in this second point
    // must also be given by 'fp2'.  The tangent known is the one for the current parameter 
    // point.  The abscissa for the minimum point of the parabola is returned.
    double parabolicEstimate(double p2, double fp2, const Point& dir);

    // Estimate the minimum of a parabola where three points are known.  If the minimum is found
    // to lie within the distance 'max_from_b' from b, the point is returned as 'u', and the
    // function returns 'true'.  Otherwise, it returns 'false', and 'u' is not calculated.
    static bool parabolicFit(double a, double fa, 
			     double b, double fb, 
			     double c, double fc, 
			     double max_from_b, double&u);

    // determine the smallest detectable change around the value 'c'
    static double numerical_tolerance(double c);

    // Determine the numerical tolerance around a given point in a given direction
    static double numerical_tolerance(const Point& x, const Point& dir);

    // Determine the numerical tolerance around a given point in a given direction, also
    // considering the change in function value (fx is the function at x, dfx is the derivative)
    static double numerical_tolerance(const Point& x, const Point& dir,
				      const double fx, const double dfx);


    // A helper function for linminBrent.  A new point has been evaluated (u, fu), and we must
    // update the other points used by the algorithm.
    static void adjustBrackets(double& a, double& fa,
			       double& b, double& fb,
			       double& c, double& fc,
			       double& b2, double& fb2,
			       double& b3, double& fb3,
			       double u, double fu);

    static void shift3(double& a, double& b, double& c, double d) 
	{ a = b; b = c; c = d;}
    static void shift2(double& a, double& b, double c) 
	{ a = b; b = c;}
};

}; // end namespace Go

#include "GoTools/utils/GeneralFunctionMinimizer_implementation.h"

#endif // _GENERAL_FUNCTIONMINIMIZER_H

