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

#ifndef _BSPLINEBASIS_H
#define _BSPLINEBASIS_H

#include <assert.h>
#include <vector>
#include "GoTools/utils/Values.h"
#include "GoTools/geometry/Streamable.h"
#include "GoTools/utils/config.h"

#ifdef NDEBUG
#define CHECK(object)
#else
#define CHECK(object) (object)->check()
#endif

namespace Go
{

    /** Class representing a B-spline basis of a spline space.
     *  This basis is defined by its order, its dimension (number of 
     *  basis functions) and its knotvector.
     */

class GO_API BsplineBasis : public Streamable
{
public:
    /// Default constructor, making an empty basis
    BsplineBasis()
	: num_coefs_(0), order_(0), last_knot_interval_(0)
    {
    }

    /// Constructor
    /// \param number Number of basis functions in space
    /// \param order order of the basis functions (polynomial degree + 1)
    /// \param knotstart iterator pointing to start of defining knotvector
    /// The length of the knotvector is deduced automatically, and should be
    /// equal to 'order' + 'number'.
    template <typename RandomIterator>
    BsplineBasis(int number,
		 int order,
		 RandomIterator knotstart)
	: num_coefs_(number), order_(order),
	  knots_(knotstart, knotstart + number + order),
	  last_knot_interval_(order - 1)
    {
	check();
    }

    /// Constructor
    /// \param order order of the basis functions (polynomial degree + 1)
    /// \param knotstart iterator pointing to start of defining knotvector
    /// \param knotend iterator pointing to one-past-end of defining knotvector.
    /// The number of basis functions will be deduced automatically by subtracting
    /// the 'order' from the length of the knotvector.
    template <typename RandomIterator>
    BsplineBasis(int order,
		   RandomIterator knotstart,
		   RandomIterator knotend)
	: num_coefs_((int)(knotend - knotstart - order)), order_(order),
	  knots_(knotstart, knotend),
	  last_knot_interval_(order - 1)
    {
	check();
    }

    /// Virtual destructor, enables safe inheritance. 
    virtual ~BsplineBasis();


    /// Change the defining data of an already-existing BsplineBasis.
    /// \param number Number of basis functions in space
    /// \param order of the basis functions (polynomial degree + 1)
    /// \param knotstart iterator pointing to start of defining knotvector.
    template <typename RandomIterator>
    void setData(int number, int order, RandomIterator knotstart)
    {
	num_coefs_ = number;
	order_ = order;
	knots_.clear();
	knots_.insert(knots_.end(), knotstart, knotstart + number + order);
	last_knot_interval_ = order - 1;
    }

    /// Quick swap of two BsplineBasises.
    /// \param other the BsplineBasis to swap with.
    void swap(BsplineBasis& other);

    /// Read the BsplineBasis from a stream (ASCII mode).
    /// \param is the stream from which the BsplineBasis is read.
    virtual void read(std::istream& is);

    /// Write the BsplineBasis to a stream (ASCII mode).
    /// \param os the stream to which the BsplineBasis is written
    virtual void write(std::ostream& os) const;

    /// Read the BsplineBasis from a stream (binary mode).
    /// \param is the stream from which the BsplineBasis is read.
    virtual void read_bin(std::istream& is);

    /// Write the BsplineBasis to a stream (binary mode).
    /// \param os the stream to which the BsplineBasis is written
    virtual void write_bin(std::ostream& os) const;

    /// Find the interval in which the parameter value 't' lies.
    /// \param t the parameter value to test.
    int knotInterval(double t) const;

    /// Create a vector containing the basis values in a given parameter.
    /// \param t the parameter at which to evaluate the basis functions
    /// \param derivs the number of function derivatives to calculate for each nonzero
    ///               basis function at 't'.
    /// \return a vector containing order() * ('derivs'+1) values, which are
    /// the function values and derivatives up to order 'derivs' of the nonzero basis 
    /// functions at parameter value 't'.  The order in which these values are stored
    /// is as follows: first all function values for the 'order' nonzero basis functions,
    /// then their first derivatives, then their second, etc..
    std::vector<double> computeBasisValues(double t, int derivs=0) const;

    /// Space for the basis values is allocated outside and passed in.
    /// Needed space is equal to order()*('derivs'+1) doubles. 
    /// \param t the parameter at which to evaluate the basis function
    /// \param basisvals_start pointer to the memory area where the result will be
    ///                        written.  This must be allocated by the user.  
    ///                        The order in which the result values are written is similar
    ///                        to computeBasisValues(double, int)
    /// \param derivs number of derivatives that should be evaluated for each nonzero
    ///               basis function (derivs = 0 => only function values will be computed).
    /// \param resolution accuracy for determining whether a given parametric value lies
    ///                   exactly 'on' a knot.
    void computeBasisValues(double t,
			    double* basisvals_start,
			    int derivs = 0,
			    double resolution=1.0e-12) const; 

    /// Compute basis values for many points simultaneously.
    /// \param parvals_start pointer to the start of list of parameters where you 
    ///                      want to evaluate the basis functions
    /// \param parvals_end pointer to one-past-end of list of parameters where you 
    ///                    want to evaluate the basis functions
    /// \param basisvals_start pointer to range where the result of the evaluations
    ///                        will be written.  The necessary memory is equal to
    ///                        'derivs' multiplied by order() multiplied by length of 
    ///                        parameter list, and should be allocated by user before calling
    ///                        the function.  The order in which these values are stored is:
    ///                        first all values concerning the first parameter in the input
    ///                        range, then all values concerning the second parameter in the 
    ///                        input range, etc.  For each such parameter, the concerned
    ///                        values are stored in the order described in 
    ///                        computeBasisValues(double, int)
    /// \param knotinter_start pointer to range where the indexes to the parameter intervals
    ///                        in which the evaluated parameter values are located.  ie. 
    ///                        the index at knotinter_start[i] refer to the parameter interval
    ///                        where the i'th parameter in the range
    ///                        [parvals_start, parvals_end] is located.  The total size of this
    ///                        range is therefore equal to (parvals_end - parvals_start), and
    ///                        should be allocated by the user.
    /// \param  derivs         number of derivatives that should be evaluated for each nonzero
    ///                        basis function (derivs = 0 => only function values will be 
    ///                        computed).
    void computeBasisValues(const double* parvals_start,
			    const double* parvals_end,
			    double* basisvals_start,
			    int* knotinter_start,
			    int derivs = 0) const;

    /// This function is similar to computeBasisValues(double, int), except that the 
    /// values are calculated from the left, as opposed to the default right-evaluation.
    /// \see computeBasisValues()
    std::vector<double> computeBasisValuesLeft(double tval, int derivs ) const;

    /// This function is similar to computeBasisValues(double, double*, int), except that the
    /// values are calculated from the  left, as opposed to the default right-evaluation.  
    /// \param resolution tolerance deciding whether a parameter value is situated 
    /// "directly on" a knot.  (In such situations, it might matter whether we evaluate
    /// from left or from left).
    /// \see computeBasisValues()
    /// \relates computeBasisValues()
    void computeBasisValuesLeft(double tval, 
				double* basisvals_start,
				int derivs,
				double resolution=1.0e-12) const;

    /// This function is similar to computeBasisValues(const double*, const double*, 
    /// double*, int*, int), except that the values are calculated from the left, as opposed
    /// to the default right-evaluation.
    void computeBasisValuesLeft(const double* parvals_start,
				const double* parvals_end,
				double* basisvals_start,
				int* knotinter_start,
				int derivs) const;

    /// We want the new knot vector to be a mirror image of the old one,
    /// translated so that is starts and ends in the same values as before.
    void reverseParameterDirection();

    /// Rescales the knotvector so that parameter values lay in the specified interval.
    /// \param new_start start of specified interval.
    /// \param new_end end of specified interval.
    void rescale (double new_start, double new_end);

    /// Insert a knot into the knotvector.
    /// \param apar parameter value of the new knot.
    void insertKnot(double apar);

    /// Adapt the knot vector to an increased polynomial order
    /// The given order is the expected new polynomial order of the basis
    void increaseOrder(int order);

    /// Find the interval in which a given parameter lies.  
    /// \param t the parameter we want to locate within an interval
    /// \param tol if 't' is within DEFAULT_PARAMETER_EPSILON of an existing knot,
    ///            't' will be changed to that knot.  To some callers of this function,
    ///            that may be the primary wanted effect of this function.
    int knotIntervalFuzzy(double& t, double tol = DEFAULT_PARAMETER_EPSILON) const;

    /// Insert several knots into the knotvector
    /// \param new_knots a STL vector containing the new knots to insert into the vector
    void insertKnot(const std::vector<double>& new_knots);

    /// Remove a knot from the knotvector.
    /// \param old_knot  Parameter value indicating the knot to be removed.  The 
    ///                  corresponding knot is the largest one which is less than or equal
    ///                  to 'old_knot'.
    void removeKnot(double old_knot);

    /// Return the index of the left knot of the last knot interval accessed for some
    /// reason (evaluation, etc.) in this basis.
    inline int lastKnotInterval() const;
    
    /// Get the greville parameter starting from a certain knot
    /// \param index the index of the knot from which the greville parameter should be 
    ///              computed.
    /// \return the corresponding greville parameter
    inline double grevilleParameter(int index) const;

    /// If there is a knot in the knotvector with the value 'parval', its multiplicity
    /// is returned.  Otherwise, 0 i returned.
    /// \param parval the parameter value for which we want to know the knot multiplicity
    /// \return the multiplicity in case the knot exists, 0 otherwise.
    int knotMultiplicity(const double parval) const;

    /// Get the vector of knot multiplicities. To get the vector of corresponding
    /// knot values, use knotsSimple().
    /// \param multiplicities the resulting vector of knot multiplicities
    void knotMultiplicities(std::vector<int>& muliplicities) const;

    /// return the multiplicity of the first or the last knot
    /// \param atstart if this is 'true' then the multiplicity of the first knot will be
    /// returned.  Otherwise, the multiplicity of the last knot will be returned. 
    /// \return the corresponding knot multiplicity 
    int endMultiplicity(bool atstart) const;

    /// Store all knots without multiplicity int one array
    /// \param result upon function return, this vector holds all the knots
    void knotsSimple(std::vector<double>& result) const;

    /// Query the number of elements, number of non-degenerate knot intervals,
    /// in the knot vector
    int numElem() const;

    /// Find the first knot, the last knot, and all knots where not all B-splines have
    /// C^N-continuity, i.e. knots with multiplicity at least the order of the B-splines
    /// minus N
    /// \param cNdisconts upon function return, this vector holds all the knots found
    /// \param depth the value of N, the degre of smoothness
    void cNDiscontinuities(std::vector<double>& cNDisconts, int depth) const;

    /// The minimum continuity at inner knots
    int getMinContinuity() const;

    /// query whether the basis has a k-regular knotvector (start and end knot has 
    /// multiplicity equal to order().
    /// \return whether this basis has a k-regular knotvector or not.
    bool isKreg() const
    {
	return (endMultiplicity(true) == order_ &&
		endMultiplicity(false) == order_);
    }

    /// Get the range of basis functions that are involved when evaluating at a 
    /// given parameter. 
    /// \param tpar the parameter where we want to determine the influential range
    /// \param first_coef on function completion, will point to start of relevant range
    /// \param last_coef on function completion, will point to end of relevant range
    void coefsAffectingParam(double tpar, int& first_coef, int& last_coef) const;
    
    /// Query the number of basis functions
    /// \return the number of basis functions
    int numCoefs() const
    { return num_coefs_; }
    
    /// Query the order of the basis functions
    /// \return the order of the basis functions
    int order() const
    { return order_; }

    /// Query the start value of the parameter range
    /// \return the start value of the parameter range
    double startparam() const
    { return knots_[order_-1]; }

    /// Query the end value of the parameter range
    /// \return the end value of the parameter range
    double endparam() const
    { return knots_[num_coefs_]; }

    /// Get a const-iterator to the beginning of the knotvector
    /// \return an iterator to the beginning of the knotvector
    std::vector<double>::const_iterator begin() const
    { return knots_.begin(); }

    /// Get a const-iterator to the end of the knotvector
    /// \return an iterator to the end of the knotvector
    std::vector<double>::const_iterator end() const
	{ return knots_.end(); }

    /// Get an iterator to the beginning of the knotvector
    /// \return an iterator to the beginning of the knotvector
    std::vector<double>::iterator begin()
	{ return knots_.begin(); }
    /// Get an iterator to the end of the knotvector
    /// \return an iterator to the end of the knotvector
    std::vector<double>::iterator end()
	{ return knots_.end(); }

    std::vector<double> getKnots()
      {
	return knots_;
      }

    /// Checks the validity of the object. Throws at error.
    void check() const;

    /// Checks the validity of the object. Does not throw.
    bool isOK() const;

    /// Look for indistinct knots
    bool indistinctKnots(double tol, std::vector<double>& first_knotval) const;

    /// Check equality of two spline spaces
    bool sameSplineSpace(const BsplineBasis& other, 
			 double tol = DEFAULT_PARAMETER_EPSILON) const;

    /// Return knots existing in the other spline space, but not in this one
    std::vector<double> 
      missingKnots(const BsplineBasis& other, 
		   double tol = DEFAULT_PARAMETER_EPSILON) const;

    /// Return the basis for the subspace limited by the interval [tmin,tmax].
    BsplineBasis subBasis(double tmin, double tmax, double knot_diff_tol = 1e-05) const;

    /// Make a basis function with the same inner knots, but with increased degree
    /// NB! The continuity of the space is increased
    /// If the basis is periodic, the new knots will be inserted at the start and end
    /// of the parameter interval, i.e. the continuity at the seam will be reduced
    /// Order is equal to polynomial degree + 1
    /// The function throws if the specified order is less than the existing one
    BsplineBasis extendedBasis(int order) const;

private:
    // Data members
    int num_coefs_;
    int order_;
    std::vector<double> knots_;
    mutable int last_knot_interval_;


};


inline int BsplineBasis::lastKnotInterval() const
{
    return last_knot_interval_;
}


inline double BsplineBasis::grevilleParameter(int index) const
{
    assert(order_ > 1); // The greville abscissae is the average over 'order_ - 1' knots.
    double greville = 0.0;
    for (int i = 1; i < order(); ++i)
        greville += knots_[index+i];

    return greville/(order() - 1.0);
}


} // namespace Go


#endif // This is what is 'ended': #ifndef _BSPLINEBASIS_H
