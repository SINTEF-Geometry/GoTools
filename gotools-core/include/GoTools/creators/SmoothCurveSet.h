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

#ifndef _SMOOTHCURVESET_H
#define _SMOOTHCURVESET_H

#include "GoTools/creators/ConstraintDefinitions.h"
#include "GoTools/geometry/SplineCurve.h"

namespace Go
{

  /// Side constraint on modification of curve set. Defines a specific
  /// relations between two curves in given parameter values
struct cvSetConstraint
{
  // Let pos1 be the evaluation of the first cv wrt to the cv1_der_'th
  // derivative, pos2 evaluation of the 2nd wrt to the cv2_der_'th derivative.
  // Let furthermore sign = (opp_) ? -1 : 1.
  // Then we should have: pos1 = sign*pos2.

  cvSetConstraint(int cv1_id, double cv1_par, int cv1_der,
		  int cv2_id, double cv2_par, int cv2_der, bool opp)
    : cv1_id_(cv1_id), cv1_par_(cv1_par), cv1_der_(cv1_der),
      cv2_id_(cv2_id), cv2_par_(cv2_par), cv2_der_(cv2_der), opp_(opp)
    {;}

  /// Index of first curve.
  int cv1_id_; 
  /// Parameter in first curve.
  double cv1_par_;
  /// The derivative of first curve involved in the constraint. 
  int cv1_der_; 
  /// Index of second curve.
  int cv2_id_;
  /// Parameter in second curve.
  double cv2_par_;
  /// The derivative of second curve involved in the constraint. 
  int cv2_der_;
  /// If true the evaluation of the second curve should be negated.
  bool opp_; 

};

/// Smoothing, point interpolation and point approximation applied to a set of curves
/// while maintaining a set of continuity conditions between the curves, exact or
/// approximative.
class SmoothCurveSet
{
 private:
     /// Struct for storing integral information of a curve
   typedef struct integralInfo
    {
	// Parameters used in integration
	std::vector<double> vec_;
	double ***integral_;  // Array used to store integrals of inner product
	// of derivatives of B-splines 
	// The 1st index runs over the derivatives, the 2nd & 3rd in 
        // the B-spline parameters
	bool integralset_; // Whether integral1 & integral2 have been computed.
	int der_; // The number of derivatives to compute.
      
	/// Constructor.
	integralInfo() 
	{ integral_ = 0; integralset_ = false; der_ = -1; }

	/// Destructor.
	~integralInfo()
	{ erase(); }

	/// Resize/set the struct variables based in input agrguments.
	/// \param ider the new number of derivatives to compute.
	/// \param in1 number of coefficients in the u direction.
	/// \param in2 number of coefficients in the v direction.
      void resize(int ider, int in)
	{
	    int ki, kj;
	    vec_.resize((ider+1)*in*in);
	    std::fill(vec_.begin(), vec_.end(), 0.0);

	    integral_ = new double**[ider+1];

	    for (ki=0; ki<=ider; ki++)
		{
		    integral_[ki] = new double*[in];

		    for (kj=0; kj<in; kj++)
			integral_[ki][kj] = &vec_[(ki*in+kj)*in];

		}

	    integralset_ = 0;
	    der_ = ider;
	}

	/// Free the memory of the arrays in the struct.
	void erase()
	{
	    int ki;
	    for (ki=0; ki<=der_; ki++)
		{
		    delete [] integral_[ki];
		}
	    delete [] integral_;
	    integralset_ = false;
	}
    } integralInfo;

public:
   /// Default constructor to the class
  SmoothCurveSet();        
  // SmoothCurveSet. Initializes class variable.

  /// Destructor.
  ~SmoothCurveSet();        

  /// Initializes data given by an intermediate set of curves.
  /// For each curve there exists a vector coef_known (of size equal to the number of
  /// coefficients in the curve)
  /// Input is array of iterators to first element.
  int attach(std::vector<shared_ptr<SplineCurve> >& incvs,
	     std::vector<int>& seem,
	     std::vector<std::vector<int> >& coef_known,
	     int numSideConstraints = 0);

  // @@@ VSK. Can it be relevant to use different weights for 
  // different curves, and how to define the weights in that case?
  /// Compute the smoothing part of the equation system.
  int setOptimize(double weight1, double weight2, double weight3);

  /// Compute matrices for least squares approximation.
  /// Each curve is assigned a number of data points with corresponding
  /// parameter values and weights
  int setLeastSquares(const std::vector<std::vector<double> >& pnts,
		      const std::vector<std::vector<double> >& param_pnts,
		      const std::vector<std::vector<double> >& pnt_weights,
		      double weight);

//   int setApproxSideConstraints(sideConstraintSetPntrArray&
// 			       constraints,
// 			       double weight);

  /// Add term for approximation of the original curves
  void setApproxOrig(double weight);

/*   // Compute matrices for approximation of normal directions. */
/*   // The number of std::vectors corresponds to number of sfs in set. */
/*   int setNormalCond(const std::vector<std::vector<double> >& pnts, */
/* 		    const std::vector<std::vector<double> >& param_pnts, */
/* 		    const std::vector<std::vector<double> >& pnt_weights, */
/* 		    double weight); */

  /// We add the interpolation conditions as linear side constraints.
  /// Assuming the degrees of freedom are sufficient (i.e. that the input
  /// curve provided by the user has enough knots).
  // Well, if the user wants to approximate the interpolation pts
  // there is a setLeastSquares routine which does just that (and it even
  // allows separate weights).
  void setInterpolationConditions(const std::vector<std::vector<double> >& pnts,
				  const std::vector<std::vector<double> >& param_pnts,
				  const std::vector<std::vector<int> >& der,
				  bool appr_constraints, double appr_wgt,
				  int* jstat);

  /// Set linear side constraints between the coefs in (possibly different)
  /// input cvs.
  int 
    setCvSetConstraints(const std::vector<shared_ptr<cvSetConstraint> >& cv_set_constraints,
			bool appr_constraints, double appr_wgt);

  /// We may have side constraints which are not suitable for exact equality as
  /// spline solution space may not be large enough. We therefore allow using
  /// least squares to minimize the error.
  /// This applies in particular to constraint involving higher order
  /// derivatives.
  /// Assuming input is preprocessed (all coefs in constraints are free).
  int setApproxSideConstraints(std::vector<shared_ptr<sideConstraintSet> >& constraints,
			       double weight);

  /// Solve equation system, and produce output curves.
  int equationSolve(std::vector<shared_ptr<SplineCurve> >& curves);

  /// The contribution to the equation system from the approximation of 
  /// normal directions.
  int setOrthCond(const std::vector<std::vector<double> >& pnts,
		  const std::vector<std::vector<double> >& param_pnts,
		  double weight);

  /// Add side constraints to the functional (Lagrange multiplier).
  /// Assuming input is preprocessed (all coefs in constraints are free).
  /// If replace_constraints==true the old constraints are removed prior
  /// to adding new constraints.
  void setSideConstraints(std::vector<shared_ptr<sideConstraintSet> >& constraints,
			  bool replace_constraints);


private:


  std::vector<shared_ptr<integralInfo> > cv_integral_; // size nmb_cvs

  int idim_;                // Dimension of geometry space.
  int kdim_;                // Normal conditions.
  int ider_;                // Maximum derivative involved in the computations. 
  std::vector<int> cont_seem_;  // Number of rows affected by continuity
                                    // at the seem.

  // The input curves
  std::vector<shared_ptr<SplineCurve> > cvs_;

  const int copyCoef_;


  // Parameters used to define the specific input curve.
  std::vector<std::vector<double> > coef_array_; // Array with curve coefficients.

  std::vector<std::vector<int> > coefknown_;
  std::vector<std::vector<int> > pivot_;

  // No coefs are assumed to be known.
  int kncond_;
  int knconstraint_; // Number of side constraints.

  int kpointer_; // Used to differ corresponding coefs + whether coef is known.

  // Storage of the equation system.
  std::vector<double> gmat_;         // Matrix at left side of equation system.  
  std::vector<double> gright_;       // Right side of equation system. 

  // Set pointers between identical coefficients at a periodic seem
  // (i.e. c0 cont).
  // If possible, update fixed coefficients at the seem.
  void preparePeriodicity(int cvidx, int seem);

  // Set periodicity constraints for cvs with seem[cvidx] > 1.
  // Expects that the gmat and gright have been initialized.
  int setPeriodicity();

//   // Set constraints on approximative C1-continuity at a seem
//   void setC1AtSeem(int cvidx, double weight);
  
//   // Set constraints on approximative C2-continuity at a seem
//   void setC2AtSeem(int cvidx, double weight);

  // Update constraints by adding known coefs to the right hand side
  // of constraint expression.
  int updateSideConstraints(std::vector<shared_ptr<sideConstraintSet> >& constraints,
			    const std::vector<std::vector<int> >& coef_known);



  // Extract the linear side contraints expression in der in cv in tpar.
  std::vector<std::pair<std::pair<int,int>, double > >
    getSideConstraint(int cv_id,
		      double tpar,
		      int der,
		      int sign,
		      int* jstat);

  int get_min_deriv(shared_ptr<SplineCurve> cv, double support_mult);


    // We update weights according to spline space of curves.
    // Size of weights should be 4 (i.e. smoothing & appr terms).
    int setWeights(double weights[], double new_weights[]);

    int set_weights(shared_ptr<SplineCurve> cv, double support_mult,
		    double weights[], double new_weights[]);
    void spline_space_cont(shared_ptr<SplineCurve> cv, int& nmbc);

};

} // end namespace Go

#endif // _SMOOTHCURVESET_H

