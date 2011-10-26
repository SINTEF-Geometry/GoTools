//===========================================================================
//                                                                           
// File: SmoothCurveSet.h                                                  
//                                                                           
// Created: 
//                                                                           
// Author: Sverre Briseid, Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _SMOOTHCURVESET_H
#define _SMOOTHCURVESET_H

#include "GoTools/creators/ConstraintDefinitions.h"
#include "GoTools/geometry/SplineCurve.h"

namespace Go
{

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

  int cv1_id_; // Index of first curve.
  double cv1_par_; // Parameter in first curve.
  int cv1_der_; // The derivative of first curve involved in the constraint.
  int cv2_id_;
  double cv2_par_;
  int cv2_der_;
  bool opp_; // If true the evaluation of the second curve should be negated.

};


class SmoothCurveSet
{
 private:
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

  SmoothCurveSet();         // Default constructor to the class
  // SmoothCurveSet. Initializes class variable.

  ~SmoothCurveSet();        // Destructor.

  // Initializes data given by an intermediate surface.
  // For each sf there exists a vector coef_known (of size kn1*kn2)
  // Input is array of iterators to first element.
  int attach(std::vector<std::shared_ptr<SplineCurve> >& incvs,
	     std::vector<int>& seem,
	     std::vector<std::vector<int> >& coef_known,
	     int numSideConstraints = 0);

  // @@@ VSK. Can it be relevant to use different weights for 
  // different curves, and how to define the weights in that case?
  // Compute the smoothing part of the equation system.
  int setOptimize(double weight1, double weight2, double weight3);

  // Compute matrices for least squares approximation.
  int setLeastSquares(const std::vector<std::vector<double> >& pnts,
		      const std::vector<std::vector<double> >& param_pnts,
		      const std::vector<std::vector<double> >& pnt_weights,
		      double weight);

//   int setApproxSideConstraints(sideConstraintSetPntrArray&
// 			       constraints,
// 			       double weight);

  void setApproxOrig(double weight);

/*   // Compute matrices for approximation of normal directions. */
/*   // The number of std::vectors corresponds to number of sfs in set. */
/*   int setNormalCond(const std::vector<std::vector<double> >& pnts, */
/* 		    const std::vector<std::vector<double> >& param_pnts, */
/* 		    const std::vector<std::vector<double> >& pnt_weights, */
/* 		    double weight); */

  // We add the interpolation conditions as linear side constraints.
  // Assuming the degrees of freedom are sufficient (i.e. that the input
  // curve provided by the user has enough knots).
  // Well, if the user wants to approximate the interpolation pts
  // there is a setLeastSquares routine which does just that (and it even
  // allows separate weights).
  void setInterpolationConditions(const std::vector<std::vector<double> >& pnts,
				  const std::vector<std::vector<double> >& param_pnts,
				  const std::vector<std::vector<int> >& der,
				  bool appr_constraints, double appr_wgt,
				  int* jstat);

  // Set linear side constraints between the coefs in (possibly different)
  // input cvs.
  int 
    setCvSetConstraints(const std::vector<std::shared_ptr<cvSetConstraint> >& cv_set_constraints,
			bool appr_constraints, double appr_wgt);

  // We may have side constraints which are not suitable for exact equality as
  // spline solution space may not be large enough. We therefore allow using
  // least squares to minimize the error.
  // This applies in particular to constraint involving higher order
  // derivatives.
  // Assuming input is preprocessed (all coefs in constraints are free).
  int setApproxSideConstraints(std::vector<std::shared_ptr<sideConstraintSet> >& constraints,
			       double weight);

  // Solve equation system, and produce output curves.
  int equationSolve(std::vector<std::shared_ptr<SplineCurve> >& curves);

  int setOrthCond(const std::vector<std::vector<double> >& pnts,
		  const std::vector<std::vector<double> >& param_pnts,
		  double weight);

  // Add side constraints to the functional (Lagrange multiplier).
  // Assuming input is preprocessed (all coefs in constraints are free).
  // If replace_constraints==true the old constraints are removed prior
  // to adding new constraints.
  void setSideConstraints(std::vector<std::shared_ptr<sideConstraintSet> >& constraints,
			  bool replace_constraints);


private:


  std::vector<std::shared_ptr<integralInfo> > cv_integral_; // size nmb_cvs

  int idim_;                // Dimension of geometry space.
  int kdim_;                // Normal conditions.
  int ider_;                // Maximum derivative involved in the computations. 
  std::vector<int> cont_seem_;  // Number of rows affected by continuity
                                    // at the seem.

  // The input curves
  std::vector<std::shared_ptr<SplineCurve> > cvs_;

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
  int updateSideConstraints(std::vector<std::shared_ptr<sideConstraintSet> >& constraints,
			    const std::vector<std::vector<int> >& coef_known);



  // Extract the linear side contraints expression in der in cv in tpar.
  std::vector<std::pair<std::pair<int,int>, double > >
    getSideConstraint(int cv_id,
		      double tpar,
		      int der,
		      int sign,
		      int* jstat);

  int get_min_deriv(std::shared_ptr<SplineCurve> cv, double support_mult);


    // We update weights according to spline space of curves.
    // Size of weights should be 4 (i.e. smoothing & appr terms).
    int setWeights(double weights[], double new_weights[]);

    int set_weights(std::shared_ptr<SplineCurve> cv, double support_mult,
		    double weights[], double new_weights[]);
    void spline_space_cont(std::shared_ptr<SplineCurve> cv, int& nmbc);

};

} // end namespace Go

#endif // _SMOOTHCURVESET_H

