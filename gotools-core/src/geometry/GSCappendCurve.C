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
//#define DEBUG

#include "GoTools/geometry/SplineCurve.h"
#include <algorithm>
#include <math.h>
#include <fstream>

#include <iterator> // For back_inserter.  This one should be required by VC++ and GCC as well...

#ifdef __BORLANDC__
#include <iterator> // For back_inserter.  This one should be required by VC++ and GCC as well...
#endif

using std::back_inserter;

namespace Go{


//===========================================================================
void SplineCurve::appendCurve(ParamCurve* other_curve, int continuity,
			      double& dist, bool repar, double tol)
//===========================================================================
{
    SplineCurve* other_cv = dynamic_cast<SplineCurve*>(other_curve);
    ALWAYS_ERROR_IF(other_cv == 0,
		"Given an empty curve or not a SplineCurve.");
    ALWAYS_ERROR_IF(dim_ != other_cv->dimension(),
		    "The curves must lie in the same space.");

    ALWAYS_ERROR_IF(continuity < -1 || continuity + 1 > order(),
		    "Specified continuity not attainable.");

#ifdef DEBUG
    // DEBUG OUTPUT
    std::ofstream of0("basis0.g2");
    writeStandardHeader(of0);
    write(of0);
    other_cv->writeStandardHeader(of0);
    other_cv->write(of0);
#endif

    // Making sure the curves have the same order. Raise if necessary.
    int diff_order = order() - other_cv->order();
    if (diff_order > 0)
      other_cv->raiseOrder(diff_order);
    else if (diff_order < 0)
      raiseOrder(abs(diff_order));

    // Make sure that we have k-regularity at meeting ends.
    makeKnotEndRegular();
    other_cv->makeKnotStartRegular();

    // Ensure that either none of the curves or both are rational 
    if (rational_ && !other_cv->rational())
      other_cv->representAsRational();
    if (!rational_ && other_cv->rational())
      representAsRational();

 #ifdef DEBUG
   // DEBUG OUTPUT
    std::ofstream of1("basis1.g2");
    writeStandardHeader(of1);
    write(of1);
    other_cv->writeStandardHeader(of1);
    other_cv->write(of1);
#endif

    if (rational_)
      {
	// Set end weight to 1
	setBdWeight(1.0, false);
	other_cv->setBdWeight(1.0, true);
      }

    // Reparametrization (translatation and mult.) of other_cv->basis().knots_ .
    if (repar && continuity > 0) {

      if (rational_)
	{
	  // The weights corresponding to the two coefficients close to the
	  // joints should be equal
	  equalBdWeights(false);
	  other_cv->equalBdWeights(true);
	}
    }

 #ifdef DEBUG
   // DEBUG OUTPUT
    std::ofstream of1_2("basis1_2.g2");
    writeStandardHeader(of1_2);
    write(of1_2);
    other_cv->writeStandardHeader(of1_2);
    other_cv->write(of1_2);
#endif

 #ifdef DEBUG
    // DEBUG OUTPUT
    std::ofstream of2("basis2.dat");
    basis_.write(of2);
    of2 << std::endl;
    other_cv->basis_.write(of2);
    of2 << std::endl;
#endif

    if (continuity > -1 && rational_)
      {
	// Make sure that the rational curves have equal weights in the end
	int k2 = (numCoefs() - 1) * (dim_+1);
	double frac = rcoefs_[k2+dim_]/other_cv->rcoefs_[dim_];
	int kn = other_cv->numCoefs();
	for (int j=0; j<kn*(dim_+1); ++j)
	  other_cv->rcoefs_[j] *= frac;
      }

#ifdef DEBUG
    // DEBUG OUTPUT
    std::ofstream of3("basis3.dat");
    basis_.write(of3);
    of3 << std::endl;
    other_cv->basis_.write(of3);
    of3 << std::endl;
#endif

    if (repar && continuity > 0) {

	double sum1 = 0;
	double sum2 = 0;
	if (rational_)
	  {
	    int k2 = (numCoefs() - 1) * (dim_+1);
	    int k1 = (numCoefs() - 2) * (dim_+1);
	    for (int j = 0; j < dim_; ++j)
	      {
		double t0 = (rcoefs_[k2 + j] - rcoefs_[k1 + j])*rcoefs_[k2 + dim_] -
		  rcoefs_[k2 + j]*(rcoefs_[k2 + dim_] - rcoefs_[k1 + dim_]);
		sum1 += t0*t0;
	      }
	    sum1 = sqrt(sum1)/(rcoefs_[k2+dim_]*rcoefs_[k2+dim_]);

	    k2 = dim_+1;
	    k1 = 0;
	    for (int j = 0; j < dim_; ++j)
	      {
		double t0 = (other_cv->rcoefs_[k2 + j] - 
			     other_cv->rcoefs_[k1 + j])*other_cv->rcoefs_[k1 + dim_] -
		  other_cv->rcoefs_[k1 + j]*(other_cv->rcoefs_[k2 + dim_] - 
					     other_cv->rcoefs_[k1 + dim_]);
		sum2 += t0*t0;
	      }
	    sum2 = sqrt(sum2)/(other_cv->rcoefs_[k1+dim_]*other_cv->rcoefs_[k1+dim_]);
	  }
	else
	  {
	    for (int j = 0; j < dim_; ++j) {
	      sum1 += (coefs_[(numCoefs() - 1) * dim_ + j] -
		       coefs_[(numCoefs() - 2) * dim_ + j]) *
		(coefs_[(numCoefs() - 1) * dim_ + j] -
		 coefs_[(numCoefs() - 2) * dim_ + j]);
	      sum2 += (other_cv->coefs_[dim_ + j] - other_cv->coefs_[j]) *
		(other_cv->coefs_[dim_ + j] - other_cv->coefs_[j]);
	    }
	    sum1 = sqrt(sum1);
	    sum2 = sqrt(sum2);
	  }

#ifdef DEBUG
    // DEBUG OUTPUT
    std::ofstream of4("basis4.dat");
    basis_.write(of4);
    of4 << std::endl;
    other_cv->basis_.write(of4);
    of4 << std::endl;
#endif

	if (sum1 > 1.0e-14) { // @@sbr We should have a universal noise-tolerance.
	  double del1 = basis_.begin()[numCoefs()] - basis_.begin()[numCoefs() - 1];
	  double del2 = other_cv->basis_.begin()[order()] -
	    other_cv->basis_.begin()[order() - 1];
	  double k = sum2*del1/(sum1*del2);
	    other_cv->basis_.rescale(endparam(), endparam() +
				     k * (other_cv->basis_.begin()
					  [other_cv->numCoefs() + order() - 1] -
					  other_cv->basis_.startparam()));
	} else {
	    MESSAGE("Curve seems to be degenerated in end pt!");
	}
    } else {
        other_cv->basis_.rescale(endparam(),
				 endparam() +
				 (other_cv->basis_.begin()
				  [other_cv->numCoefs() + order() - 1] -
				  other_cv->basis_.startparam()));
    }

    // Join the curve-segments (i.e. set endpoints equal), given that...
    if (continuity != -1) {
	for (int j = 0; j < dim_; ++j) {
	    other_cv->coefs_[j] = 
		coefs_[(numCoefs() - 1)*dim_ + j] =
		(coefs_[(numCoefs() - 1)*dim_ + j] + other_cv->coefs_[j])/2;
	    }
    }

    double tpar = basis_.endparam();
    int ti = numCoefs() + order() - 1; // Index of last occurence of tpar.

#ifdef DEBUG
    // DEBUG OUTPUT
    std::ofstream of5("basis5.dat");
    basis_.write(of5);
    of5 << std::endl;
    other_cv->basis_.write(of5);
    of5 << std::endl;
#endif

    // Add other_cv's coefs.
    if (rational_)
      {
	// Ensure identity of the weight in the joint
	other_cv->setBdWeight(rcoefs_[rcoefs_.size()-1], true);
      
	rcoefs_.insert(rcoefs_end(), other_cv->rcoefs_begin(), 
		       other_cv->rcoefs_end());
      }
    else
      coefs_.insert(coefs_end(), other_cv->coefs_begin(), other_cv->coefs_end());

#ifdef DEBUG
    // DEBUG OUTPUT
    std::ofstream of("basis.dat");
    basis_.write(of);
    of << std::endl;
    other_cv->basis_.write(of);
    of << std::endl;
#endif

    // Make an updated basis_ .
    std::vector<double> new_knotvector;
    std::copy(basis_.begin(), basis_.end(), std::back_inserter(new_knotvector));
    std::copy(other_cv->basis_.begin() + order(), other_cv->basis_.end(),
	 std::back_inserter(new_knotvector));
    basis_ = BsplineBasis(order(), new_knotvector.begin(), new_knotvector.end());

    if (rational_)
      updateCoefsFromRcoefs();

    SplineCurve orig_curve = *this; // Save curve for later estimates.

    // Obtain wanted smoothness.
    int i;
    try {
    for (i = 0; i < continuity + 1; ++i)
	removeKnot(tpar);
    }
    catch (...)
    {
	// Leave the knots
    }

    // Estimate distance between curve and smoothed curve: 
    // Raise (copy of) smoothed curve to original spline space 
    // and calculate max distance between corresponding spline-coefs.
    SplineCurve raised_smooth_curve = *this;
    std::vector<double> knots;
    for (i = 0; i < continuity + 1; ++i) knots.push_back(tpar);
    raised_smooth_curve.insertKnot(knots);
    double sum, root_sum;
    dist = 0;
    for (i = std::max(0, ti - (continuity + 1) - order()); 
	 i < ti - order() + continuity + 1; ++i) {
	sum = 0;
	for (int j = 0; j < dim_; ++j)
	    sum += (orig_curve.coefs_[i * dim_ + j] - 
		    raised_smooth_curve.coefs_[i * dim_ + j])
		    * (orig_curve.coefs_[i * dim_ + j] - 
		       raised_smooth_curve.coefs_[i * dim_ + j]);
	// to avoid use of the max function, which is likely to cause
	// trouble with the Microsoft Visual C++ Compiler, the following 
	// two lines are added, and the third one is commented out.
	root_sum = sqrt(sum);
	dist = dist > root_sum ? dist : root_sum;
	//dist = std::max(dist, sqrt(sum));
    }
}


//===========================================================================
void SplineCurve::appendCurve(ParamCurve* cv, bool repar)
//===========================================================================
{
    // For the time being assuming C1 as default.
    int cont = 1;
    double dist_dummy = 0;
    appendCurve(cv, cont, dist_dummy, repar);
}


//===========================================================================
void SplineCurve::makeKnotStartRegular()
//===========================================================================
{
    // Testing whether knotstart is already d+1-regular.
    if (basis_.begin()[0] < basis_.begin()[order() - 1]) {
	
	double tpar = basis_.startparam();
	int ti = order() - 1; // Index of last occurence of tpar (in other_curve).
	int mt = 1; // Multiplicity of tpar.
	
	while ((basis_.begin()[ti - mt] == tpar) && (mt < order())) ++mt;
	std::vector<double> new_knots;
	for (int i = 0; i < order() - mt; ++i) new_knots.push_back(tpar);
	insertKnot(new_knots);
	coefs_.erase(coefs_begin(), coefs_begin() + (order() - mt) * dim_);
	basis_ = BsplineBasis(order(), basis_.begin() + order() - mt,
				basis_.end());
    }
}

//===========================================================================
void SplineCurve::makeKnotEndRegular()
//===========================================================================
{
    // Testing whether knotstart is already d+1-regular.
    if (basis_.begin()[numCoefs()] < basis_.begin()[numCoefs() + order() - 1]) {

	double tpar = basis_.endparam();
	int ti = numCoefs(); // Index of first occurence of tpar.
	int mt = 1; // Multiplicity of tpar.

	while ((basis_.begin()[ti + mt] == tpar) && (mt < order())) ++mt;
	std::vector<double> new_knots;
	for (int i = 0; i < order() - mt; ++i) new_knots.push_back(tpar);
	insertKnot(new_knots);
	coefs_.erase(coefs_begin() + (numCoefs() - order() + mt) * dim_,
		     coefs_begin() + numCoefs() * dim_);
	basis_ = BsplineBasis(order(), basis_.begin(),
				basis_.begin() + numCoefs() + mt);
    }
}



} // namespace Go;
