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

#include "GoTools/geometry/BsplineBasis.h"
#include <algorithm>
#include <iomanip>
#include <assert.h>
#include <math.h>


using namespace Go;
using std::vector;
using std::streamsize;


//-----------------------------------------------------------------------------
BsplineBasis::~BsplineBasis()
//-----------------------------------------------------------------------------
{
}

//-----------------------------------------------------------------------------
void BsplineBasis::swap(BsplineBasis& other)
//-----------------------------------------------------------------------------
{
    std::swap(num_coefs_, other.num_coefs_);
    std::swap(order_, other.order_);
    knots_.swap(other.knots_);
    std::swap(last_knot_interval_, other.last_knot_interval_);
}

//-----------------------------------------------------------------------------
void BsplineBasis::read(std::istream& is)
//-----------------------------------------------------------------------------
{
    bool is_good = is.good();
    if (!is_good) {
	THROW("Invalid object header!");
    }
    is >> num_coefs_ >> order_;
    is_good = is.good();
    if (!is_good) {
	THROW("Invalid object header!");
    }
    knots_.resize(num_coefs_+order_);
    for (int i=0; i<num_coefs_+order_; ++i)
	is >> knots_[i];
    last_knot_interval_ = order_-1;
    CHECK(this);

    is_good = is.good();
    if (!is_good) {
	THROW("Invalid object header!");
    }
}


//-----------------------------------------------------------------------------
void BsplineBasis::write(std::ostream& os) const
//-----------------------------------------------------------------------------
{
    streamsize prev = os.precision(15);
    os << num_coefs_ << ' ';
    os << order_ << '\n';
    os << knots_[0];
    for (int i = 1; i<num_coefs_+order_; ++i)
	os << ' ' << knots_[i];
    os << '\n';
    os.precision(prev);   // Reset precision to it's previous value
}


//-----------------------------------------------------------------------------
void BsplineBasis::read_bin(std::istream& is)
//-----------------------------------------------------------------------------
{
    assert(sizeof(int) == 4);
    assert(sizeof(double) == 8);

    // reading number of coefficients
    char* dp = reinterpret_cast<char*>(&num_coefs_);
    is.read(dp, sizeof(int));
    
    // reading order
    dp = reinterpret_cast<char*>(&order_);
    is.read(dp, sizeof(int));

    // reading knotvector
    knots_.resize(num_coefs_ + order_);
    dp = reinterpret_cast<char*>(&knots_[0]);
    is.read(dp, sizeof(double) * (num_coefs_ + order_));
    last_knot_interval_ = order_-1;
    CHECK(this);

}

//-----------------------------------------------------------------------------
void BsplineBasis::write_bin(std::ostream& os) const
//-----------------------------------------------------------------------------
{
    assert(sizeof(int) == 4);
    assert(sizeof(double) == 8);
    
    // writing number of coefficients
    const char* dp = reinterpret_cast<const char*>(&num_coefs_);
    os.write(dp, sizeof(int));
    // writing order
    dp = reinterpret_cast<const char*>(&order_);
    os.write(dp, sizeof(int));

    // writing knotvector
    dp = reinterpret_cast<const char*>(&knots_[0]);
    os.write(dp, sizeof(double) * (num_coefs_ + order_) );
}

//-----------------------------------------------------------------------------
void BsplineBasis::reverseParameterDirection()
//-----------------------------------------------------------------------------
{
    // We want the new knot vector to be a mirror image of the old one,
    // translated so that is starts and ends in the same values as before.

    int num = num_coefs_ + order_;
    double start = knots_[0];
    double end = knots_[num - 1];
    // We do the mirroring in place by doing it from
    // both ends at the same time.
    // If we have an n-regular knot vector, we detect it,
    // and skip the outer knots.
    int startind = 1; // "1-regular" knot vector. Trivial case.
    // int startind = 0; // "1-regular" knot vector. Trivial case.
    while ((knots_[startind] == knots_[startind-1])
    	   && (knots_[num - startind - 1] == knots_[num - startind]))
    	++startind;
    double tmp_last, tmp1, tmp2;
    int i, j, k;
    for (i = startind; i <= (num-1)/2; ++i) {
	tmp_last = knots_[num - 1 - i];
	tmp1 = (start+end)-knots_[i];
	tmp2 = (start+end)-tmp_last;
      // tmp1 = knots_[i] + knots_[num-1-i];
      // tmp2 = knots_[i];
      // tmp2 = tmp1 - tmp2;
      // tmp1 = tmp1 - knots_[num-1-i];
	if (fabs(tmp2-tmp1) < 1.0e-12)
	    tmp1 = tmp2;
	knots_[num - 1 - i] = tmp1;
	knots_[i] = tmp2;
	// knots_[i] = tmp1;
	// knots_[num - 1 - i] = tmp2;
    }

    // Ensure exact knot multiplicity
    double tol = 1.0e-15;
    for (i=0; i<num; i=j)
      {
	for (j=i+1; j<num; ++j)
	  if (fabs(knots_[j] - knots_[i]) > tol)
	    break;

	double med = 0.0;
	for (k=i; k<j; ++k)
	  med += knots_[k];
	med /= (double)(j-i);
	
	for (k=i; k<j; ++k)
	  knots_[k] = med;
      }
	
    CHECK(this);
}


//-----------------------------------------------------------------------------
void BsplineBasis::rescale(double new_start, double new_end)
//-----------------------------------------------------------------------------
{
    int n = num_coefs_+order_;
    double old_start = knots_[order_ - 1];
    double old_end = knots_[num_coefs_];
    for (int i=0; i<n; ++i) {
	double factor = (knots_[i] - old_start)/(old_end - old_start);
	knots_[i] = new_start*(1.0-factor) + new_end*factor;
    }
    //CHECK(this); 
}



//-----------------------------------------------------------------------------
void BsplineBasis::insertKnot(double apar)
//-----------------------------------------------------------------------------
{
    ++num_coefs_;
    std::vector<double>::iterator where
	= std::lower_bound(knots_.begin(), knots_.end(), apar);
    knots_.insert(where, apar);
}


//-----------------------------------------------------------------------------
void BsplineBasis::increaseOrder(int order)
//-----------------------------------------------------------------------------
{
  if (order_ >= order)
    return;  // Polynomial order already high enough

  int nn = order - order_;  // Number of new knots to add for each position
  
  vector<double> knotval;
  knotsSimple(knotval);
  vector<double> newknots(knotval.size()*nn);
  int ki;
  size_t kj, kr;
  for (kr=0, kj=0; kr<knotval.size(); ++kr)
    for (ki=0; ki<nn; ++ki)
      newknots[kj++] = knotval[kr];

  vector<double> knots2(knots_.size() + newknots.size());;
  std::merge(knots_.begin(), knots_.end(),
	     newknots.begin(), newknots.end(),
	     knots2.begin());
  order_ = order;
  knots_ = knots2;
  num_coefs_ = (int)knots_.size() - order_;
  last_knot_interval_ = order_-1;
}

//-----------------------------------------------------------------------------
void BsplineBasis::insertKnot(const std::vector<double>& new_knots)
//-----------------------------------------------------------------------------
{
    for(size_t i = 0; i < new_knots.size(); ++i) {
	insertKnot(new_knots[i]);
    }
}

//-----------------------------------------------------------------------------
void BsplineBasis::removeKnot(double old_knot)
//-----------------------------------------------------------------------------
{
    --num_coefs_;
    std::vector<double>::iterator where =
	std::lower_bound(knots_.begin(), knots_.end(), old_knot);
    knots_.erase(where);
}

// Validity checking
//-----------------------------------------------------------------------------
void BsplineBasis::check() const
//-----------------------------------------------------------------------------
{
    ALWAYS_ERROR_IF(order_ <= 0, "Order must be positive.");
    ALWAYS_ERROR_IF(num_coefs_ < order_,
		    "Number of vertices must be at least equal to order.");

    int n = num_coefs_ + order_;
    int ks = (int)knots_.size();
    ALWAYS_ERROR_IF(ks != n,
		    "The length of the knot vector must equal order + number of vertices.");

    // Check that the knot vector is increasing
    int i;
    for (i=0; i<n-1; i++)
	ALWAYS_ERROR_IF(knots_[i+1]<knots_[i],
			"Knot vector must be nondecreasing.");

    // Check that knot multiplicity is <= order
    int m = 1;
    for (i=0; i<n-1; i++) {
	if (knots_[i+1]==knots_[i])
	    m++;
	else
	    m = 1;

	if (m > order_) {
	  THROW("Knot multiplicity (" << knots_[i] << ") exceeds order.");
	}
    }
}

//-----------------------------------------------------------------------------
bool BsplineBasis::isOK() const
//-----------------------------------------------------------------------------
{
    if (order_ <= 0)
	return false;
    if (num_coefs_ < order_)
	return false;

    int n = num_coefs_ + order_;
    int ks = (int)knots_.size();
    if (ks != n)
	return false;

    // Check that the knot vector is increasing
    int i;
    for (i=0; i<n-1; i++)
	if (knots_[i+1]<knots_[i])
	    return false;

    // Check that knot multiplicity is <= order
    int m = 1;
    for (i=0; i<n-1; i++) {
	if (knots_[i+1]==knots_[i])
	    m++;
	else
	    m = 1;

	if (m > order_) 
	    return false;
    }

    return true;  // No errors found
}

//-----------------------------------------------------------------------------
bool BsplineBasis::indistinctKnots(double tol, std::vector<double>& first_knotval) const
//-----------------------------------------------------------------------------
{
    // First check if the knot vector has any formal errors
    bool is_ok = isOK();
    if (!is_ok)
	return true;  // Error in knot vector

    // Travers knot vector and check for distinct knot with distance less than the
    // given tolerance
    double tk = knots_[0];
    bool indistinct = false;
    for (size_t ki=1; ki<knots_.size(); ++ki)
    {
	if (knots_[ki] > tk)
	{
	    if (knots_[ki] - tk < tol)
	    {
		first_knotval.push_back(tk);
		indistinct = true;
	    }
	    tk = knots_[ki];
	}
    }
    return indistinct;
}

//-----------------------------------------------------------------------------
bool BsplineBasis::sameSplineSpace(const BsplineBasis& other, 
				   double tol) const
//-----------------------------------------------------------------------------
{
  if (num_coefs_ != other.num_coefs_)
    return false;

  if (order_ != other.order_)
    return false;

  for (size_t ki=0; ki<knots_.size(); ++ki)
    if (fabs(knots_[ki] - other.knots_[ki]) > tol)
      return false;

  return true;
}

//-----------------------------------------------------------------------------
vector<double> BsplineBasis::missingKnots(const BsplineBasis& other, 
					  double tol) const
//-----------------------------------------------------------------------------
{
  vector<double> diff;
  size_t ki, kj;
  for (ki=0; ki<knots_.size(); )
    for (kj=0; kj<other.knots_.size(); )
      {
	if (knots_[ki] < other.knots_[kj]-tol)
	  ki++;
	else if (knots_[ki] > other.knots_[kj]+tol)
	  {
	    diff.push_back(other.knots_[kj]);
	    kj++;
	  }
	else
	  {
	    ki++;
	    kj++;
	  }
      }
  return diff;
}

//-----------------------------------------------------------------------------
BsplineBasis BsplineBasis::subBasis(double tmin, double tmax,
				    double knot_diff_tol) const
//-----------------------------------------------------------------------------
{
    int start_ind = knotIntervalFuzzy(tmin, knot_diff_tol);
    int end_ind = knotIntervalFuzzy(tmax, knot_diff_tol);

    std::vector<double> new_knots(knots_.begin() + start_ind + 1, knots_.begin() + end_ind);
    int start_mult = (knots_[start_ind] == tmin) ? knotMultiplicity(tmin) - 1 : 0;
    new_knots.insert(new_knots.begin(), order_ - start_mult, tmin);
    new_knots.insert(new_knots.end(), order_, tmax);

    return BsplineBasis(order_, new_knots.begin(), new_knots.end());
}


//-----------------------------------------------------------------------------
BsplineBasis BsplineBasis::extendedBasis(int order) const
//-----------------------------------------------------------------------------
{
  // Check input
  ASSERT(order >= order_);

  std::vector<double> new_knots(knots_.begin(), knots_.end());  
  double ta = new_knots[order_-1];
  double tb = new_knots[num_coefs_];
  new_knots.insert(new_knots.begin()+order_-1, order-order_, ta);

  int kn = num_coefs_ + order - order_;
  new_knots.insert(new_knots.begin()+kn, order-order_, tb);

  return  BsplineBasis(order, new_knots.begin(), new_knots.end());
}

//-----------------------------------------------------------------------------
int BsplineBasis::endMultiplicity(bool atstart) const
//-----------------------------------------------------------------------------
{
  int mult = 1;
  int ki;
  if (atstart)
    {
      for (ki=order_-1; ki>0 && knots_[ki] == knots_[ki-1]; ki--, mult++);
    }
  else
    {
      for (ki=num_coefs_; ki<num_coefs_+order_-1 && 
	     knots_[ki] == knots_[ki+1]; ki++, mult++);
    }
// #ifdef _MSC_VER
//   return (mult < order_) ? mult : order_;
// #else
  return std::min(mult,order_);
// #endif
}


//-----------------------------------------------------------------------------
int BsplineBasis::knotMultiplicity(const double parval) const
//-----------------------------------------------------------------------------
{
  if (parval == knots_[num_coefs_])
    return endMultiplicity(false);

  int index = knotInterval(parval);

  if (knots_[index] != parval)
    return 0;

  int mult = 1;
  while ((index - mult > -1) && knots_[index] == knots_[index-mult])
    ++mult;

  return mult;
}


//===========================================================================
void BsplineBasis::cNDiscontinuities(std::vector<double>& cNDisconts, int depth) const
//===========================================================================
{
  double startpar = knots_[0];
  double par = startpar;
  int mult = 1;

  for (int i = 1; i < num_coefs_+order_; ++i)
    {
      if (par == knots_[i])
	++mult;
      else
	{
	  if (par == startpar || mult >= order_ - depth)
	    cNDisconts.push_back(par);
	  par = knots_[i];
	  mult = 1;
	}
    }

  cNDisconts.push_back(par);
}

//===========================================================================
int BsplineBasis::getMinContinuity() const
//===========================================================================
{
  int min_cont = order_;
  int ki, kj;
  for (ki=order_; ki<num_coefs_; ki+=kj)
    {
      for (kj=ki+1; kj<num_coefs_; ++kj)
	if (knots_[kj] > knots_[ki])
	  break;
      min_cont = std::min(min_cont, order_-(kj-ki)-1);
    }
  return min_cont;
}

//===========================================================================
void BsplineBasis::knotMultiplicities(vector<int>& multiplicities) const
//===========================================================================
{
    vector<double> knot_values;
    knotsSimple(knot_values);
    multiplicities.clear();
    for (size_t i = 0; i < knot_values.size(); ++i) {
        multiplicities.push_back(knotMultiplicity(knot_values[i]));
    }
    return;
}
  

//===========================================================================
void BsplineBasis::knotsSimple(std::vector<double>& result) const
//===========================================================================
{
  double startpar = knots_[0];
  double par = startpar;

  for (int i = 1; i < num_coefs_+order_; ++i)
    {
      if (par != knots_[i])
	{
	  result.push_back(par);
	  par = knots_[i];
	}
    }

  result.push_back(par);
}
  

//===========================================================================
int BsplineBasis::numElem() const
//===========================================================================
{
  int num = 0;
  double startpar = knots_[0];
  double par = startpar;

  for (int i = 1; i < num_coefs_+order_; ++i)
    {
      if (par != knots_[i])
	{
	  num++;
	  par = knots_[i];
	}
    }

  return num;
}
  


//-----------------------------------------------------------------------------
void BsplineBasis::coefsAffectingParam(double tpar, int& first_coef, int& last_coef) const
//-----------------------------------------------------------------------------

{
  int kleft = knotInterval(tpar);
  // first_coef should be the index of first coef affecting tpar.
  if (tpar == endparam())
    {
      first_coef = num_coefs_ - 1;
      last_coef = first_coef;
    }
  else
    {
      first_coef = kleft - order_ + 1;
      // last_coef should be the index of last coef affecting tpar.
      if (tpar == startparam())
	{
	  last_coef = first_coef;
	}
      else
	{
	  int knot_mult = 0;
	  while ((knots_[kleft-knot_mult] == tpar) && knot_mult < order_)
	    ++knot_mult;
	  last_coef = first_coef + order_ - 1 - knot_mult;
	}
    }
};
