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

#ifndef _SMOOTHCURVE_H_
#define _SMOOTHCURVE_H_


//   -----------------------------------------------------------------------
//      Interface file for class SmoothCurve
//   -----------------------------------------------------------------------
//
//       Used to modify a tensor product B-spline surface with respect
//       to conditions on smoothness, editing constraints and boundary
//       conditions.
//
//       Implementation of the member functions are given in the
//       following files:
//
//          1. SmoothCurve.C
//             a. SmoothCurve()
//             b. ~SmoothCurve()
//             c. attach()
//             d. setOptim()
//             e. setLeastSquares()
//             f. equationSolve()
//
//   -----------------------------------------------------------------------
//    Written by: Vibeke Skytt                            28-04.1998.
//   -----------------------------------------------------------------------

#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/creators/ConstraintDefinitions.h"
//#include "newmat.h"
#include <memory>

#include <vector>


namespace Go
{

/// This class is used to generate a spline curve from a given spline
/// space that approximates a set of weighed data points.  In addition,
/// user can set other conditions on things like curve smoothness or 
/// periodicity.
/// To use this object, you must first instanciate it, then run
/// the 'attach()' function, the 'setOptim()' function, the 'setLeastSquares()'
/// function and finally the 'equationSolve()' function.  Optionally, the functions
/// 'setPeriodicity()' and 'setSideConstraints()' can be run before 'equationSolve()',
/// if the corresponding conditions should be speficied.
class SmoothCurve 
{

 public:

    /// Constructor.  Initialises the object with a given spatial dimension
    /// (default is 3).
    /// \param dim dimension of the geometry space.
    SmoothCurve(int dim = 3);

    /// Destructor.
    ~SmoothCurve();

    /// Spefify the spline space we want to use, 

    /// Initializes data given by an intermediate curve.
    /// \param incurve curve defining the initial spline curve.
    /// \param coef_known array defining the free coefficients.
    ///                   Size equal to the number of coefficients in incurve.
    ///                   Value 0 = not known, value 1 = known.
    /// \param numSideConstraints the number of linear side constrains for
    ///                           the approximation.
    /// \return status value: 0 = OK, negative not OK.
    int attach(const shared_ptr<SplineCurve>& incurve,
	       int coef_known[],
	       int numSideConstraints = 0);

    /// Compute the smoothing part of the equation system.
    /// The sum of the weights should lie in the unit interval.
    /// \param weight1 weight for smoothing with respect to the 1st derivative.
    /// \param weight2 weight for smoothing with respect to the 2nd derivative.
    /// \param weight3 weight for smoothing with respect to the 3rd derivative.
    void setOptim(const double weight1, const double weight2,
		  const double weight3);

    /// Compute matrices for least squares approximation.
    /// \param pnts the input points (for 3D: (x0, y0, z0, x1, y1, z1, ...))
    ///             defining the approximation part of the smoothing problem.
    /// \param param_pnts the parameters for the input pnts.
    /// \param pnt_weights weight attached to each point to be approximated.
    ///                    Should lie in the unit interval.
    ///                    Typically they are all 1.0.
    /// \param weight multiplier of all weights in the pnt_weights vector
    void setLeastSquares(std::vector<double>& pnts,
			 std::vector<double>& param_pnts,
			 std::vector<double>&  pnt_weights,
			 double weight);

    /// Set periodicity constraints in one par. dir.
    /// \param cont the wanted continuity across the seam.
    /// \param weight weight given to approximative continuity conditions.
    ///               Should lie inside the unit interval.
    void setPeriodicity(int cont, double weight);

    /// Set linear side constraints (linear equation involving the free coefficients)
    /// to the minimization problem. The problem is solved using the method of
    /// Lagrange multipliers.
    /// \param constraints the linear side constraints.
    void setSideConstraints(std::vector<sideConstraint>& constraints);

    /// Solve equation system, and produce output curve.
    /// \param curve the curve resulting from the smoothing process.
    void equationSolve(shared_ptr<SplineCurve>& curve);


 private:
    int idim_;             // Dimension of geometry space.
    int kdim_;             // Dimension of homogeneous space.
    int ider_;             // Maximum derivative involved in the computations. 
//     int cont_bound[2];    // Fixed continuity at each endpoint, i.e. the number
//     // of coefficients not to be changed.
    int cont_seam_;        // Number of rows affected by continuity at the seam
    int kcond_; // Number of unknowns in equation system (# of coefs + kconstraint).
    int kconstraint_; // # of constraints.

    int *coefknown_;       // Array indicating the status of coefficients, i.e.
    // free, fixed, not involved, equal to a given
    // coefficients.
    std::vector<int> pivot_;    // Array giving the position of the free coefficients
    // in the equation system.   

    // The input curve
    shared_ptr<SplineCurve> qcurve_;

    // Parameters defining the spline space.
    int kk_;         // Order of curve.
    int kn_;         // Number of coefficients of curve.                
    std::vector<double>::const_iterator st_;     // Pointer to knot  of curve.

    // Rational case
    bool rational_;    // Is curve rational?
    shared_ptr<SplineCurve> bspline_curve_;  // 1-dimensional rational curve with
                                                    // 0's as control values, and weights from 
                                                    // input curve. Used to calculate B-splines
                                                    // for the rational case

    // Parameters used to define the specific input spline curve.
    std::vector<double>::iterator scoef_; // Pointer to coefficients.          

    // Parameters used in integration
    double ***integral_;   // Array used to store integrals of inner product

    // Storage of the equation system.
    //Matrix gmat_;       // Matrix at left side of equation system.
    //ColumnVector gright_;    // Right side of equation system. 
    std::vector<std::vector<double> > gmat_; // Matrix at left side of equation system
    std::vector<double> gright_; // Right side of equation system

    /// Allocate storage for values in integral.
    // \return return value 0 = OK, negative return value = failure.
    int prepareIntegral();

    /// Release the space occupied by the integral_.
    void releaseScratch();

    /// Ensure that the expected continuity at the seam is satisfied.
    void adjustAtSeam();

    /// Set constraints on approximative C0-continuity at the seam.
    /// \param weight the weight should lie inside the unit interval.
    void setC0AtSeam(double weight);

    /// Set constraints on approximative C1-continuity at the seam.
    /// \param weight the weight should lie inside the unit interval.
    void setC1AtSeam(double weight);
};

}// end namespace Go

#endif
