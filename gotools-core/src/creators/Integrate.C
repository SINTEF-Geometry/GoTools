#include "GoTools/creators/Integrate.h"
#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/geometry/SplineCurve.h"

using namespace Go;
using std::shared_ptr;
using std::dynamic_pointer_cast;
using std::vector;
using std::max;
using std::min;
using std::swap;

namespace {

const int indices[] = {1, 2, 3, 4, 5, 8};

const double w_2_0 = 0.5555555556;
const double w_2_1 = 0.8888888889;
const double w_3_0 = 0.3478548451;
const double w_3_1 = 0.6521451549;
const double w_4_0 = 0.2369268851;
const double w_4_1 = 0.4786286705;
const double w_4_2 = 0.5688888889;
const double w_5_0 = 0.1012285363;
const double w_5_1 = 0.2223810345;
const double w_5_2 = 0.3137066459;
const double w_5_3 = 0.3626837834;
    
const double weight[][8] =
    { {   1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0},
      {   1.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0},
      { w_2_0,  w_2_1,  w_2_0,    0.0,    0.0,    0.0,    0.0,    0.0},
      { w_3_0,  w_3_1,  w_3_1,  w_3_0,    0.0,    0.0,    0.0,    0.0},
      { w_4_0,  w_4_1,  w_4_2,  w_4_1,  w_4_0,    0.0,    0.0,    0.0},
      { w_5_0,  w_5_1,  w_5_2,  w_5_3,  w_5_3,  w_5_2,  w_5_1,  w_5_0} };

const double s_0_0 = 0.5;
const double s_1_0 = -0.5773502692;
const double s_2_0 = -0.7745966692;
const double s_3_0 = -0.8611363116;
const double s_3_1 = -0.3399810436;
const double s_4_0 = -0.9061798459;
const double s_4_1 = -0.5384693101;
const double s_5_0 = -0.9602898565;
const double s_5_1 = -0.7966664774;
const double s_5_2 = -0.5255324099;
const double s_5_3 = -0.1834346425;

const double sample[][8] =
    { { s_0_0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0},
      { s_1_0, -s_1_0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0},
      { s_2_0,    0.0, -s_2_0,    0.0,    0.0,    0.0,    0.0,    0.0},
      { s_3_0,  s_3_1, -s_3_1, -s_3_0,    0.0,    0.0,    0.0,    0.0},
      { s_4_0,  s_4_1,    0.0, -s_4_1, -s_4_0,    0.0,    0.0,    0.0},
      { s_5_0,  s_5_1,  s_5_2,  s_5_3, -s_5_3, -s_5_2, -s_5_1, -s_5_0} };


}; // anonymous namespace containing the numerical variables needed

namespace Go {

//===========================================================================
void GaussQuadValues(const BsplineBasis& basis,  // B-spline basis
		     vector<double>& parameters,  // Parameter values for all points
		     vector<double>& par_weights)  // Weight for each parameter
//===========================================================================

   //--------------------------------------------------------------------------
   //     Purpose : Store parameters and weights used for numerical integration
   //               by Gauss quadrature. All parameter values for all Bezier segments
   //               are stored in parameters, while the weights within an interval is
   //               stored only once in weights. Thus, the lenght of weights gives
   //               the number of samples for each Bezier segment.
   //--------------------------------------------------------------------------
{
   int ord = basis.order();
   int ncoef = basis.numCoefs();
   int kind = (ord-1 < 5) ? ord-1 : 5;

   int w_size = indices[kind];
   parameters.resize(w_size * (ncoef-ord+1));
   par_weights.resize(w_size);

   for (int i = 0; i < w_size; ++i)
     par_weights[i] = 0.5 * weight[kind][i];

   vector<double>::const_iterator it = basis.begin();
   for (int i = 0, par_pos = 0; i < ncoef-ord+1; ++i)
     {
       double int_start = it[i+ord-1];
       double int_end = it[i+ord];
       for (int j = 0; j < w_size; ++j, ++par_pos)
	 parameters[par_pos] = 0.5 * (sample[kind][j]*(int_end-int_start) + int_end + int_start);
     }
}


//==============================================================================
void  GaussQuadInner(const BsplineBasis& basis,  // B-spline basis
		     int ider,    // Number of derivatives to compute.
		     double lim1, // Start of parameter interval.
		     double lim2, // End of parameter interval.
		     double*** integral) // Computed integrals.
//==============================================================================

   //--------------------------------------------------------------------------
   //     Purpose : Compute all definite integrals of inner products of
   //		    derivatives of B-splines up to a given order where the
   //               differentiation is of the same order for both B-splines.
   //               The interval of integration are equal to the parameter
   //               intervals of the surface in the current par. dir.
   //
   //     Calls   : BsplineBasis::computeBasisValues  -
   //                                      Compute derivatives of B-splines.
   //
   //     Written by : Vibeke Skytt,  SINTEF SI,  09.93. (04.02)
   //--------------------------------------------------------------------------
{
   int ki, kp, k1, k2, kl, kr;
   int kder;
   int kind;
   int kleft;
   double ta, tb;
   double tval;
   double tpar;
   int ik = basis.order();
   int in = basis.numCoefs();
   vector<double>::const_iterator et = basis.begin();
   vector<double> sbder(ik*(ider+1), 0.0);

   /* Traverse all knot intervals inside the limits of integration, computing
      the integral of the inner product of two B-splines defined in this
      interval.    */

   kind = (ik-1 < 5) ? ik-1 : 5;
   for (kl=ik-1; kl<in; kl++)
   {
      ta = et[kl];
      tb = et[kl+1];
      if (tb <= lim1 || ta >= lim2) continue;

      for (kr=0; kr<indices[kind]; kr++)
      {
	 /* Compute parameter value in which to evaluate B-splines.  */

	 tpar = 0.5*(sample[kind][kr]*(tb-ta) + tb + ta);

	 /* Evaluate B-splines and derivatives of B-splines. */

         basis.computeBasisValues(tpar, &sbder[0], ider);
         kleft = basis.lastKnotInterval();

	 for (ki=kleft-ik+1, k1=0; ki<=kleft; ki++, k1++)
	    for (kp=ki, k2=k1; kp<=kleft; kp++, k2++)
	       for (kder=0; kder<=ider; kder++)
	       {
		  tval = (double)0.5*(tb-ta)*weight[kind][kr]*
		     sbder[k1*(ider+1)+kder]*sbder[k2*(ider+1)+kder];

		  integral[kder][ki][kp] += tval;
		  if (kp > ki) integral[kder][kp][ki] += tval;
	       }
      }
   }
}




//==============================================================================
void GaussQuadInnerFlat(const BsplineBasis& basis, // B-spline basis.
			int derivs, // Number of derivatives to compute.
			int start_der, // First derivative to compute.
			int gap,   // Difference between derivation order.
			double lim1, // Start of parameter interval.
			double lim2, // End of parameter interval.
			vector<double>& integral) // Computed integrals.
//==============================================================================

//--------------------------------------------------------------------------
//     Purpose : Compute all definite integrals of inner products of
//		    derivatives of B-splines up to a given order where the
//               gap between the order of differentiation on the first and
//               second B-spline is constant.
//               The computation starts at a specified derivation depths,
//               assuming the lower derivation order integrals are
//               already found.
//               The interval of integration are equal to the parameter
//               intervals of the surface in the current par. dir.
//
//     Calls   : BsplineBasis::computeBasisValues  -
//                                      Compute derivatives of B-splines.
//
//     Written by : Kjell Fredrik Pettersen, SINTEF IKT, 2009-12-16, based on
//                  GaussQuadInner()
//--------------------------------------------------------------------------
{
  int order = basis.order();
  int n_coefs = basis.numCoefs();
  vector<double>::const_iterator et = basis.begin();
  vector<double> sbder(order*(derivs+gap+1), 0.0);

  /* Traverse all knot intervals inside the limits of integration, computing
     the integral of the inner product of two B-splines defined in this
     interval.    */

  int kind = (order-1 < 5) ? order-1 : 5;
  for (int k = order - 1; k < n_coefs; ++k)   // For every Bezier segment
    {
      double ta = et[k];
      double tb = et[k + 1];
      if (tb <= lim1 || ta >= lim2) continue;

      for (int l = 0; l < indices[kind]; ++l)  // For every sample value in segment
	{
	  /* Compute parameter value in which to evaluate B-splines.  */

	  double tpar = 0.5 * (sample[kind][l]*(tb-ta) + tb + ta);

	  /* Evaluate B-splines and derivatives of B-splines. */

	  basis.computeBasisValues(tpar, &sbder[0], derivs+gap);
	  int left = basis.lastKnotInterval();

	  for (int i = left - order + 1, k1 = 0; i <= left; ++i, ++k1)  // For every first B-spline
	    for (int j = i, k2 = k1; j <= left; ++j, ++k2)     // For every second B-spline
	      for (int der = start_der; der <= derivs; ++der)  // For every derivation order
		{
		  double tval = 0.5 * (tb-ta) * weight[kind][l]
		    * sbder[k1*(derivs+gap+1)+der+gap] * sbder[k2*(derivs+gap+1)+der];
		  double tval_oposite = 0.5 * (tb-ta) * weight[kind][l]
		    * sbder[k1*(derivs+gap+1)+der] * sbder[k2*(derivs+gap+1)+der+gap];

		  integral[j-i+order-1 + (2*order-1)*(i + n_coefs*der)] += tval;
		  if (j > i)
		    integral[i-j+order-1 + (2*order-1)*(j + n_coefs*der)] += tval_oposite;
		}
	}
    }
}



//==============================================================================
void GaussQuadInner2(const BsplineBasis& basis,  // B-spline basis
		     int ider,    // Number of derivatives to compute.
		     double lim1, // Start of parameter interval.
		     double lim2, // End of parameter interval.
		     double** integral) // Computed integrals.
//==============================================================================

   //--------------------------------------------------------------------------
   //     Purpose : Compute all definite integrals of inner products of
   //		    derivatives of B-splines of a given order where the
   //               differentiation is of the same order for both B-splines.
   //               The interval of integration are equal to the parameter
   //               intervals of the surface in the current par. dir.
   //
   //     Calls   : BsplineBasis::computeBasisValues  -
   //                                      Compute derivatives of B-splines.
   //
   //     Written by : Vibeke Skytt,  SINTEF SI,  08.96. (04.02)
   //--------------------------------------------------------------------------
{
   int ki, kp, k1, k2, kl, kr;
   int kind;
   double ta, tb;
   double tval;
   double tpar;
   int kleft;
   int ik = basis.order();
   int in = basis.numCoefs();
   vector<double>::const_iterator et = basis.begin();
   vector<double> sbder(ik*(ider+1), 0.0);

   /* Traverse all knot intervals inside the limits of integration, computing
      the integral of the inner product of two B-splines defined in this
      interval.    */

   kind = (ik-1 < 5) ? ik-1 : 5;
   for (kl=ik-1; kl<in; kl++)
   {
      ta = et[kl];
      tb = et[kl+1];
      if (tb <= lim1 || ta >= lim2) continue;

      for (kr=0; kr<indices[kind]; kr++)
      {
	 /* Compute parameter value in which to evaluate B-splines.  */

	 tpar = 0.5*(sample[kind][kr]*(tb-ta) + tb + ta);

	 /* Evaluate B-splines and derivatives of B-splines. */

         basis.computeBasisValues(tpar, &sbder[0], ider);
         kleft = basis.lastKnotInterval();

	 for (ki=kleft-ik+1, k1=0; ki<=kleft; ki++, k1++)
	    for (kp=ki, k2=k1; kp<=kleft; kp++, k2++)
            {
              tval = (double)0.5*(tb-ta)*weight[kind][kr]*
                     sbder[k1*(ider+1)+ider]*sbder[k2*(ider+1)+ider];

              integral[ki][kp] += tval;
              if (kp > ki) integral[kp][ki] += tval;
            }
      }
   }
}


//==============================================================================
void  GaussQuadInnerRational(const BsplineBasis& basis,  // B-spline basis
			     int ider,    // Number of derivatives to compute.
			     double lim1, // Start of parameter interval.
			     double lim2, // End of parameter interval.
			     shared_ptr<SplineCurve> bspline_curve, // 1-dim rational 0-function defining weights
			     double*** integral) // Computed integrals.
//==============================================================================

   //--------------------------------------------------------------------------
   //     Purpose : Compute all definite integrals of inner products of
   //		    derivatives of rational B-splines up to a given order where the
   //               differentiation is of the same order for both B-splines.
   //               The interval of integration are equal to the parameter
   //               intervals of the surface in the current par. dir.
   //
   //     Calls   : BsplineBasis::computeBasisValues  -
   //                                      Compute derivatives of B-splines.
   //
   //     Written by : Kjell Fredrik Pettersen, 2009-11-26, based on
   //                  GaussQuadInnerRational by Vibeke Skytt,  SINTEF SI,  09.93. (04.02)
   //--------------------------------------------------------------------------
{
   int ki, kp, k1, k2, kl, kr;
   int kder;
   int kind;
   int kleft;
   double ta, tb;
   double tval;
   double tpar;
   int ik = basis.order();
   int in = basis.numCoefs();
   vector<double>::const_iterator et = basis.begin();
   vector<double> sbder(ik*(ider+1), 0.0);

   vector<double>::iterator bspl_it = bspline_curve->rcoefs_begin();
   vector<Point> pts(ider+1);
   for (int i = 0; i <= ider; ++i)
     pts[i] = Point(1);
   
   /* Traverse all knot intervals inside the limits of integration, computing
      the integral of the inner product of two B-splines defined in this
      interval.    */

   kind = (ik-1 < 5) ? ik-1 : 5;
   for (kl=ik-1; kl<in; kl++)
   {
      ta = et[kl];
      tb = et[kl+1];
      if (tb <= lim1 || ta >= lim2) continue;

      for (kr=0; kr<indices[kind]; kr++)
      {
	 /* Compute parameter value in which to evaluate B-splines.  */

	 tpar = 0.5*(sample[kind][kr]*(tb-ta) + tb + ta);

	 /* Evaluate B-splines and derivatives of B-splines. */

         kleft = basis.knotInterval(tpar);
	 for (ki=kleft-ik+1, k1=0; ki<=kleft; ki++, k1++)
	   {
	     bspl_it[ki*2] = 1.0;
	     bspline_curve->point(pts, tpar, ider);
	     for (int i = 0; i <= ider; ++i)
	       sbder[k1*(ider+1)+i] = pts[i][0];
	     bspl_it[ki*2] = 0.0;
	   }

	 for (ki=kleft-ik+1, k1=0; ki<=kleft; ki++, k1++)
	    for (kp=ki, k2=k1; kp<=kleft; kp++, k2++)
	       for (kder=0; kder<=ider; kder++)
	       {
		  tval = (double)0.5*(tb-ta)*weight[kind][kr]*
		     sbder[k1*(ider+1)+kder]*sbder[k2*(ider+1)+kder];

		  integral[kder][ki][kp] += tval;
		  if (kp > ki) integral[kder][kp][ki] += tval;
	       }
      }
   }
}



}; // end namespace Go
