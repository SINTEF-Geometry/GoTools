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

#ifndef _SPLINEVOLUME_H
#define _SPLINEVOLUME_H


#include "GoTools/trivariate/ParamVolume.h"
#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/RectDomain.h"
#include "GoTools/utils/ScratchVect.h"


namespace Go
{

class Interpolator;
class SplineCurve;
class DirectionCone;

/// Structure for storage of results of grid evaluation of the basis function of a spline volume.
/// Positional evaluation information in one parameter value
struct BasisPts
{
  /// Parameter tripple in which the basis functions are evaluated
    double param[3];   
  /// Index of the knot interval where the parameter value is situated for all
  /// parameter directions. The indices of the non-zero basis functions are
  /// left_idx[i]-order[i]+1, ..., left_idx[i] for i=0,1,2
    int left_idx[3];   
  /// The value of all basis functions, size equal to (degree_u+1)*(degree_v+1)*(degree_w+1)
    std::vector< double > basisValues; 

  /// Constructor
    void preparePts(double u, double v, double w, int idx_u, int idx_v, int idx_w,
		    int size)
	{
	    param[0] = u;
	    param[1] = v;
	    param[2] = w;
	    left_idx[0] = idx_u;
	    left_idx[1] = idx_v;
	    left_idx[2] = idx_w;
	    basisValues.resize(size);
	}
};

/// Structure for storage of results of grid evaluation of the basis function of s spline volume.
/// Evaluation of position and first derivatives in one parameter value
struct BasisDerivs
{
  /// Parameter tripple in which the basis functions are evaluated
    double param[3];
  /// Index of the knot interval where the parameter value is situated for all
  /// parameter directions. The indices of the non-zero basis functions are
  /// left_idx[i]-order[i]+1, ..., left_idx[i] for i=0,1,2
  int left_idx[3];   
    /// The value of all basis functions, size equal to (degree_u+1)*(degree_v+1)*(degree_w+1)
    std::vector< double > basisValues; 
    /// the derivative of all basis functions in u direction, same size as previous
    std::vector< double > basisDerivs_u; 
    /// the derivative of all basis functions in v direction, same size as previous
    std::vector< double > basisDerivs_v;
    /// the derivative of all basis functions in w direction, same size as previous
    std::vector< double > basisDerivs_w;

  /// Constructor
    void prepareDerivs(double u, double v, double w, int idx_u, int idx_v, int idx_w,
		       int size)
	{
	    param[0] = u;
	    param[1] = v;
	    param[2] = w;
	    left_idx[0] = idx_u;
	    left_idx[1] = idx_v;
	    left_idx[2] = idx_w;
	    basisValues.resize(size);
	    basisDerivs_u.resize(size);
	    basisDerivs_v.resize(size);
	    basisDerivs_w.resize(size);
	}
};

/// Structure for storage of results of grid evaluation of the basis function of a 
/// spline volume.
/// Evaluation of position, first and second derivatives in one parameter value
struct BasisDerivs2
{   
  /// Parameter tripple in which the basis functions are evaluated
    double param[3];
  /// Index of the knot interval where the parameter value is situated for all
  /// parameter directions. The indices of the non-zero basis functions are
  /// left_idx[i]-order[i]+1, ..., left_idx[i] for i=0,1,2
     int left_idx[3];
    /// The value of all basis functions, size equal to (degree_u+1)*(degree_v+1)*(degree_w+1)
    std::vector< double > basisValues; 

    /// the derivative of all basis functions in u direction, same size as previous
    std::vector< double > basisDerivs_u;
    /// the derivative of all basis functions in v direction, same size as previous
    std::vector< double > basisDerivs_v;
    /// the derivative of all basis functions in w direction, same size as previous
    std::vector< double > basisDerivs_w;

    /// the second derivative of all basis functions twice in u direction, same size as previous
    std::vector< double > basisDerivs_uu;
    /// the second derivative of all basis functions in u and v direction, same size as previous
    std::vector< double > basisDerivs_uv;
    /// the second derivative of all basis functions in u and w direction, same size as previous
    std::vector< double > basisDerivs_uw;
    /// the second derivative of all basis functions twice in v direction, same size as previous
    std::vector< double > basisDerivs_vv;
    /// the second derivative of all basis functions in v and w direction, same size as previous
    std::vector< double > basisDerivs_vw;
    /// the second derivative of all basis functions twice in w direction, same size as previous
    std::vector< double > basisDerivs_ww;

  /// Constructor
    void prepareDerivs(double u, double v, double w, int idx_u, int idx_v, int idx_w,
		       int size)
	{
	    param[0] = u;
	    param[1] = v;
	    param[2] = w;
	    left_idx[0] = idx_u;
	    left_idx[1] = idx_v;
	    left_idx[2] = idx_w;
	    basisValues.resize(size);
	    basisDerivs_u.resize(size);
	    basisDerivs_v.resize(size);
	    basisDerivs_w.resize(size);
	    basisDerivs_uu.resize(size);
	    basisDerivs_uv.resize(size);
	    basisDerivs_uw.resize(size);
	    basisDerivs_vv.resize(size);
	    basisDerivs_vw.resize(size);
	    basisDerivs_ww.resize(size);
	}
};

#define SPLINE_VOLUME_PERIOD_INFO_SIZE 3

/// \brief SplineVolume provides methodes for storing,
/// reading and manipulating rational and non-rational
/// B-spline volumes.
///
/// Non-rational B-spline volumes represented on the form
/// \f[ \sum_{i=1}^{n_1}\sum_{j=1}^{n_2}\sum_{k=1}^{n_3} P_{i,j,k} B_{i,o_1}(u) B_{j,o_2}(v) B_{k,o_3}(w), \f]
/// where the B-spline coefficients are stored in a vector called coefs.
/// The coefs are stored as
/// \f[ P_{0,0,0}, P_{1,0,0}, \ldots , P_{n_1, n_2,n_3}, \f]
/// where \f$ P_{i,j,k} \f$ is represented as dim doubles.
///
/// NURBS volumes are represented on the form
/// \f[ \frac{\sum_{i=1}^{n_1}\sum_{j=1}^{n_2} \sum_{k=1}^{n_3} h_{i,j,k}P_{i,j,k} B_{i,o_1}(u) B_{j,o_2}(v) B_{k,o_3}(w)} {\sum_{i=1}^{n_1}\sum_{j=1}^{n_2} \sum_{k=1}^{n_3} h_{i,j,k} B_{i,o_1}(u) B_{j,o_2}(v) B_{k,o_3}(w)} \f]
/// where \f$ B_{i,o} \f$ is the i'th non-rational B-spline of order o. The Projected
/// coefficients are stored in the coefs vector, i.e.
/// \f[ h_{0,0,0} \cdot P_{0,0,0}, h_{1,0,0} \cdot P_{1,0,0}, \ldots , h_{n_1, n_2, n_3} \cdot P_{n_1, n_2, n_3}. \f]
/// In addition the coefficients for the volume in projective space are kept, i.e
/// \f[ P_{0,0,0}, h_{0,0,0}, P_{1,0,0}, h_{1,0,0}, \ldots , P_{n_1, n_2, n_3}, h_{n_1, n_2, n_3}. \f]
/// With this representation volume is in the convex hull of the coefficients
/// stored in coefs both for rational and non rational volumes. 

class SplineVolume : public ParamVolume
{
public:
    /// Creates an uninitialized SplineVolume, which can only be assigned to 
    /// or read(...) into.
    SplineVolume()
	: dim_(-1), rational_(false)
    {
	for (int ki=0; ki<SPLINE_VOLUME_PERIOD_INFO_SIZE; ++ki)
	    periodicity_info_[ki] = std::make_pair(-1, 0.0);
    }

    /// Create a SplineVolume by explicitly providing all spline-related 
    /// information.
    /// \param number1 number of control points along the u-parameter
    /// \param number2 number of control points along the v-parameter
    /// \param number3 number of control points along the w-parameter
    /// \param order1 BsplineBasis order along the u-parameter
    /// \param order2 BsplineBasis order along the v-parameter
    /// \param order3 BsplineBasis order along the w-parameter
    /// \param knot1start pointer to the array describing the knotvector 
    ///                   for the u-parameter
    /// \param knot2start pointer to the array describing the knotvector
    ///                   for the v-parameter
    /// \param knot3start pointer to the array describing the knotvector
    ///                   for the w-parameter
    /// \param coefsstart pointer to the array where the control points
    ///                   are consecutively stored.  The storage order is
    ///                   such that control points along the u-parameter have
    ///                   the shortest stride (stored right after each other),
    ///                   then comes the v-parameter and finally the w-parameter.
    ///                   If the volume is rational, pay attention to the
    ///                   comments below.
    /// \param dim        dimension of the space in which the volume lies 
    ///                   (usually 3).  
    /// \param rational Specify whether the volume is rational or not.
    ///                 If the volume is rational, coefficients must be in
    ///                 the following format: 
    ///                 wP1 wP2 .... wPdim w.   Ie. a (dim+1)-dimensional form.
    ///                 (This is the same form that is used within SISL).
    template <typename RandomIterator1,
	      typename RandomIterator2,
              typename RandomIterator3,
	      typename RandomIterator4>
    SplineVolume(int number1,
		    int number2,
		    int number3,
		    int order1,
		    int order2,
		    int order3,
		    RandomIterator1 knot1start,
		    RandomIterator2 knot2start,
		    RandomIterator3 knot3start,
		    RandomIterator4 coefsstart,
		    int dim,
		    bool rational = false)
	: dim_(dim), rational_(rational),
          basis_u_(number1, order1, knot1start),
          basis_v_(number2, order2, knot2start),
          basis_w_(number3, order3, knot3start)
    {
	if (rational) {
	    int n = (dim+1)*number1*number2*number3;
	    rcoefs_.resize(n);
	    std::copy(coefsstart, coefsstart + n, rcoefs_.begin());
	    coefs_.resize(dim*number1*number2*number3);
	    updateCoefsFromRcoefs();
	} else {
	    int n = dim*number1*number2*number3;
	    coefs_.resize(n);
	    std::copy(coefsstart, coefsstart + n, coefs_.begin());
	}

	for (int ki=0; ki<SPLINE_VOLUME_PERIOD_INFO_SIZE; ++ki)
	    periodicity_info_[ki] = std::make_pair(-1, 0.0);
    }

    /// Create a SplineVolume by explicitly providing all spline-related 
    /// information.
    /// \param basis_u the BsplineBasis to be used along the u-direction
    /// \param basis_v the BsplineBasis to be used along the v-direction
    /// \param basis_w the BsplineBasis to be used along the v-direction
    /// \param coefsstart pointer to the array where the control points
    ///                   are consecutively stored. The storage order is
    ///                   such that control points along the u-parameter have
    ///                   the shortest stride (stored right after each other),
    ///                   then comes the v-parameter and finally the w-parameter.
    ///                   If the volume is rational, pay attention to the
    ///                   comments below.
    /// \param dim dimension of the space in which the volume lies (usually 3).
    /// \param rational Specify whether the volume is rational or not.
    ///                 If the volume is rational, coefficients must be in the
    ///                 following format:
    ///                 wP1 wP2 .... wPdim w.   Ie. a (dim+1)-dimensional form.
    ///                 (This is the same form that is used within SISL).
    template <typename RandomIterator>
    SplineVolume(const BsplineBasis& basis_u,
		    const BsplineBasis& basis_v,
		    const BsplineBasis& basis_w,
		    RandomIterator coefsstart,
		    int dim,
		    bool rational = false)
	: dim_(dim), rational_(rational),
	  basis_u_(basis_u),
          basis_v_(basis_v),
	  basis_w_(basis_w)
    {
	int number1 = basis_u.numCoefs();
	int number2 = basis_v.numCoefs();
	int number3 = basis_w.numCoefs();
	if (rational) {
	  int n = (dim+1)*number1*number2*number3;
	    rcoefs_.resize(n);
	    std::copy(coefsstart, coefsstart + n, rcoefs_.begin());
	    coefs_.resize(dim*number1*number2*number3);
	    updateCoefsFromRcoefs();
	} else {
	    int n = dim*number1*number2*number3;
	    coefs_.resize(n);
	    std::copy(coefsstart, coefsstart + n, coefs_.begin());
	}
	for (int ki=0; ki<SPLINE_VOLUME_PERIOD_INFO_SIZE; ++ki)
	    periodicity_info_[ki] = std::make_pair(-1, 0.0);
    }

    /// Virtual destructor, enables safe inheritance.
    virtual ~SplineVolume();

    // inherited from Streamable
    virtual void read (std::istream& is);

    // inherited from Streamable
    virtual void write (std::ostream& os) const;

    // inherited from GeomObject
    virtual BoundingBox boundingBox() const;

    // inherited from GeomObject
    /// Return the dimension of the space in which the object lies (3 is expected). 
    virtual int dimension() const;
    
    /// quick swap of two SplineVolume objects with each other
    void swap(SplineVolume& other);

    // inherited from GeomObject
    virtual ClassType instanceType() const;

    // inherited from GeomObject
    static ClassType classType()
    { return Class_SplineVolume; }

    virtual SplineVolume* clone() const
    { return new SplineVolume(*this); }

    // inherited from ParamVolume
    // (u_min, u_max, vmin, ...)
    virtual const Array<double,6> parameterSpan() const;

    // inherited from ParamVolume
    virtual void point(Point& pt, double upar, double vpar, double wpar) const;

    // Output: Partial derivatives up to order derivs (pts[0]=S(u,v),
    // pts[1]=dS/du=S_u, pts[2]=S_v, pts[3]=S_w,
    // pts[4]=S_uu, pts[5]=S_uv, pts[6]=S_uw, pts[7]=S_vv, pts[8]=S_vw, pts[9]=S_ww,
    // then S_uuu, S_uuv, S_uuw, S_uvv, S_uvw, S_uww, S_vvv, S_vvw, S_vww, S_www, ...)
    // inherited from ParamVolume
    virtual void point(std::vector<Point>& pts, 
		       double upar, double vpar, double wpar,
		       int derivs,
		       bool u_from_right = true,
		       bool v_from_right = true,
		       bool w_from_right = true,
		       double resolution = 1.0e-12) const;

    /// Get the start value for the specified parameter direction.
    /// \param i the parameter direction
    /// \return the start value for the parameter direction given by the parameter pardir
    virtual double startparam(int i) const;

    /// Get the end value for the specified parameter direction.
    /// \param i the parameter direction
    /// \return the end value for the parameter direction given by the parameter pardir
    virtual double endparam(int i) const;

    // inherited from ParamVolume
    virtual DirectionCone tangentCone(int pardir) const;

    /// Get the derivative volume. Expresses the i,j,k-th derivative of
    /// a spline volume as a spline volume.
    /// \param ider1 what derivative to use along the u parameter
    /// \param ider2 what derivative to use along the v parameter
    /// \param ider3 what derivative to use along the w parameter
    /// \return pointer to a newly constructed SplineVolume which is 
    ///         the derivative volume.  User assumes ownership of this
    ///         object, and is responsible for destroying it.
    SplineVolume* derivVolume(int ider1, int ider2, int ider3) const;

    /// Get a SplineVolume which represents a part of 'this' SplineVolume
    /// \param from_upar start value for u-parameter in the sub-volume to be 
    ///                  generated 
    /// \param from_vpar start value for v-parameter in the sub-volume to be
    ///                  generated
    /// \param from_wpar start value for w-parameter in the sub-volume to be
    ///                  generated
    /// \param to_upar end value for u-parameter in the sub-volume to be generated
    /// \param to_vpar end value for v-parameter in the sub-volume to be generated
    /// \param to_wpar end value for w-parameter in the sub-volume to be generated
    /// \param fuzzy tolerance used to determine whether given parameter values
    ///              are located directly \em knot values.
    SplineVolume* subVolume(double from_upar,
			    double from_vpar,
			    double from_wpar,
			    double to_upar,
			    double to_vpar,
			    double to_wpar,
			    double fuzzy =
			    DEFAULT_PARAMETER_EPSILON) const;

    // Split volume in specified parameter values in a given direction
    std::vector<shared_ptr<SplineVolume> > 
      split(std::vector<double>& param,
	    int pardir,
	    double fuzzy = DEFAULT_PARAMETER_EPSILON) const; 

    // inherited from ParamVolume
    virtual void closestPoint(const Point& pt,
			      double&        clo_u,
			      double&        clo_v, 
			      double&        clo_w, 
			      Point&         clo_pt,
			      double&        clo_dist,
			      double         epsilon,
			      double   *seed = 0) const;

    /// Returns the corner closest to a given point together with
    /// the associated enumeration of the corner coefficient.
    /// In degenerate cases, the enumeration will reflect an arbitrary 
    /// corner
    int closestCorner(const Point& pt, double& upar, double& vpar,
		      double& wpar, Point& corner, double& dist) const;

    /// Appends a volume to the end of 'this' volume along a specified parameter
    /// direction.
    /// \param vol the volume to append to 'this' volume
    /// \param join_dir the parameter direction along which to extend 'this' volume
    ///                 by appending the 'vol' volume. If it is equal to 0, we will
    ///                 extend the volume along the u-parameter; if it equals 1, we
    ///                 will extend the volume along the v-parameter; if it equals
    ///                 2, we will extend the volume along the w-parameter.
    /// \param cont continuity at jont, legal values are 0 and 1 (i,.e. C^0 or
    ///             C^1-continuity).
    /// \param repar The parametrization of the concerned parameter in the other 
    ///              volume will \em always be shifted so that it starts where the
    ///              parametrization of 'this' volume ends.  However, if 'repar' is 
    ///              set to 'true', it will also be \em scaled as a function of position
    ///              of control points close to the transition.
    void appendVolume(SplineVolume* vol, int join_dir,
		       int cont, bool repar=true);

    // inherited from ParamVolume
    virtual void reverseParameterDirection(int pardir);

    // inherited from ParamVolume
    virtual void swapParameterDirection(int pardir1, int pardir2);

    /// get one of the BsplineBasises of the volume
    /// \param pardir specify whether to return the BsplineBasis for the first
    ///          parameter (0), for the second parameter (1) or the third
    ///          (2) parameter.
    /// \return const reference to the requested BsplineBasis.
    const BsplineBasis& basis(int pardir) const
    {
      if (pardir == 0) return basis_u_;
      if (pardir == 1) return basis_v_;
      return basis_w_;
    }

    /// Query the number of control points along the specified parameter direction
    /// \return the number of control points along the specified parameter direction
    int numCoefs(int pardir) const
    {
      if (pardir == 0) return basis_u_.numCoefs();
      if (pardir == 1) return basis_v_.numCoefs();
      return basis_w_.numCoefs();
    }

    /// Query the order of the BsplineBasis for the specified parameter
    /// \return  the order of the BsplineBasis for the specified parameter
    int order(int pardir) const
    {
      if (pardir == 0) return basis_u_.order();
      if (pardir == 1) return basis_v_.order();
      return basis_w_.order();
    }

    /// Query the number of elements in the SplineVolume
    int numElem() const
    {
      return basis_u_.numElem()*basis_v_.numElem()*basis_v_.numElem();
    }
    
    /// Query the number of elements in one parameter direction of 
    // the SplineVolume 
    /// pardir = 0: u-direction, pardir = 1: vdirection, pardir = 2: wdirection
    int numElem(int pardir) const
    {
      if (pardir == 0) return basis_u_.numElem();
      if (pardir == 1) return basis_v_.numElem();
      return basis_w_.numElem();
    }

    /// Query whether the volume is rational
    /// \return 'true' if the volume is rational, 'false' otherwise
    bool rational() const
    { return rational_; }
    
    /// Get an iterator to the start of the internal array of non-rational control 
    /// points.
    /// \return an (nonconst) iterator to the start of the internal array of non-
    ///         rational control points
    std::vector<double>::iterator coefs_begin()
    { return coefs_.begin(); }

    /// Get an iterator to the one-past-end position of the internal array of non-
    /// rational control points
    /// \return an (nonconst) iterator to the one-past-end position of the internal
    ///         array of non-rational control points
    std::vector<double>::iterator coefs_end()
    { return coefs_.end(); }

    /// Get a const iterator to the start of the internal array of non-rational
    /// control points.
    /// \return a const iterator to the start of the internal array of non-rational
    ///         control points.
    std::vector<double>::const_iterator coefs_begin() const
    { return coefs_.begin(); }

    /// Get a const iterator to the one-past-end position of the internal array of
    /// non-rational control points.
    /// \return a const iterator to the one-past-end position of the internal array of
    ///         non-rational control points.
    std::vector<double>::const_iterator coefs_end() const
    { return coefs_.end(); }

    /// Get an iterator to the start of the internal array of \em rational control 
    /// points.
    /// \return an (nonconst) iterator ro the start of the internal array of rational
    ///         control points.
    std::vector<double>::iterator rcoefs_begin()
    { return rcoefs_.begin(); }

    /// Get an iterator to the one-past-end position of the internal array of 
    /// \em rational control points.
    /// \return an (nonconst) iterator to the start of the internal array of rational
    ///         control points.
    std::vector<double>::iterator rcoefs_end()
    { return rcoefs_.end(); }

    /// Get a const iterator to the start of the internal array of \em rational
    /// control points.
    /// \return a const iterator to the start of the internal array of rational
    ///         control points
    std::vector<double>::const_iterator rcoefs_begin() const
    { return rcoefs_.begin(); }

    /// Get a const iterator to the one-past-end position of the internal array
    /// of \em rational control points.
    /// \return a const iterator to the one-past-end position of the internal array 
    ///         of rational control points.
    std::vector<double>::const_iterator rcoefs_end() const
    { return rcoefs_.end(); }

    /// Get an iterator to the start of the internal array of active control 
    /// points.
    /// \return an (nonconst) iterator to the start of the internal array of 
    ///         rational or non-rational control points
    std::vector<double>::iterator ctrl_begin()
    { return rational_ ? rcoefs_.begin() : coefs_.begin(); }

    /// Get an iterator to the one-past-end position of the internal array of 
    /// active control points
    /// \return an (nonconst) iterator to the one-past-end position of the internal
    ///         array of rational or non-rational control points
    std::vector<double>::iterator ctrl_end()
    { return rational_ ? rcoefs_.end() : coefs_.end(); }

    /// Get a const iterator to the start of the internal array of active
    /// control points.
    /// \return a const iterator to the start of the internal array of rational or 
    ///         non-rational control points.
    std::vector<double>::const_iterator ctrl_begin() const
    { return rational_ ? rcoefs_.begin() : coefs_.begin(); }

    /// Get a const iterator to the one-past-end position of the internal array of
    /// active control points.
    /// \return a const iterator to the one-past-end position of the internal array of
    ///         rational or non-rational control points.
    std::vector<double>::const_iterator ctrl_end() const
    { return rational_ ? rcoefs_.end() : coefs_.end(); }

    /// Replace one specified coefficient (local enumeration)
    void replaceCoefficient(int ix, Point coef);

     /// Return all weights corresponding to this volume, a non-rational volume
    /// has all weights equal to one
    void getWeights(std::vector<double>& weights) const;

    /// Set the parameter domain to a given box
    /// \param u1 new min. value of first parameter span
    /// \param u2 new max. value of first parameter span
    /// \param v1 new min. value of second parameter span
    /// \param v2 new max. value of second parameter span
    /// \param w1 new min. value of third parameter span
    /// \param w2 new max. value of third parameter span
    void setParameterDomain(double u1, double u2, double v1, double v2, double w1, double w2);

    /// Insert a new knot in the knotvector of the given parameter direction
    /// \param pardir parameter direction in which to insert the knots
    ///               (u=0, v=1, w=2)
    /// \param apar the parameter value at which a new knot will be inserted
    void insertKnot(int pardir, double apar);
    
    /// Insert new knots in the knotvector of the given parameter direction
    /// \param pardir parameter direction in which to insert the knots
    ///               (u=0, v=1, w=2)
    /// \param new_knots a vector containing the parameter values of the
    ///                  new knots to insert.
    void insertKnot(int pardir, const std::vector<double>& new_knots);

    /// Remove a knot from the knotvector of the given parameter direction. 
    /// \param pardir the parameter direction (0, 1 or 2)
    /// \param tpar the parameter value of the knot to be removed
    void removeKnot(int pardir, double tpar);

    /// Inserts knots in the specified knot vector, such that all knots
    /// have multiplicity order
    /// \param pardir the parameter direction (0, 1, or 2)
    void makeBernsteinKnots(int pardir);

    /// Returns the number of knot intervals in the specified knot vector.
    /// \param pardir parameter direction
    /// \return the number of knot intervals in the knotvector for the
    ///         specified parameter direction
    int numberOfPatches(int pardir) const;

    /// Returns the size of the knot interval (knot[iknot],knot[iknot+1])
    /// in the specified direction. An index outside the legal range will 
    /// result in a zero knot span
    double knotSpan(int pardir, int iknot) const;

    /// Raise the order of the spline volume as indicated by parameters.
    /// \param raise_u the order of the BsplineBasis associated with the first
    ///                parameter will be raised this many times.
    /// \param raise_v the order of the BsplineBasis associated with the second
    ///                parameter will be raised this many times.
    /// \param raise_w the order of the BsplineBasis associated with the third
    ///                parameter will be raised this many times.
    void raiseOrder(int raise_u, int raise_v, int raise_w);

    /// Generate and return a SplineSurface that represents a constant parameter 
    /// surface on the volume
    /// \param parameter value of the fixed parameter
    /// \param pardir 0 if the surface is constant in the u-parameter,
    ///               1 if the surface is constant in the v-parameter,
    ///               2 if the surface is constant in the w-parameter.
    /// \return pointer to a newly constructed SplineSurface representing the 
    ///         specified constant parameter surface.  It is the user's reponsibility
    ///         to delete it when it is no longer needed.
    virtual SplineSurface* constParamSurface(double parameter,
					     int pardir) const;

    /// Fetch all boundary surfaces corresponding to the volume.
    virtual std::vector<shared_ptr<ParamSurface> > 
	getAllBoundarySurfaces() const;

    /// Fetch all 6 boundary surfaces corresponding to the volume.
    /// Sequence: u_min, u_max, v_min, v_max, w_min, w_max
    std::vector<shared_ptr<SplineSurface> > 
	getBoundarySurfaces(bool do_clear = false) const;

    /// Fetch one boundary surface
    shared_ptr<SplineSurface> getBoundarySurface(int idx,
						 bool do_clear = false) const
      {
	(void)getBoundarySurfaces(do_clear);
	return bd_sfs_[idx];
      }

    /// Fetch parameter boundaries of a specified element
    /// elem_par - parameter values of element boundaries, sequence umin,
    /// umax, vmin, vmax, wmin, wmax
    void getElementBdPar(int elem_ix, double elem_par[]) const;

    /// Fetch all boundary surfaces of a specified element
    /// elem_par - parameter values of element boundaries, sequence umin,
    /// umax, vmin, vmax, wmin, wmax
    std::vector<shared_ptr<SplineSurface> > getElementBdSfs(int elem_ix,
							    double elem_par[]) const;

    /// Evaluate points and derivatives on an entire grid, taking computational
    /// advantage over calculating all these values simultaneously rather than
    /// one-by-one.
    /// \param num_u number of values to evaluate along first parameter direction
    /// \param num_v number of values to evaluate along second parameter direction
    /// \param num_w number of values to evaluate along third parameter directon
    /// \param points upon function return, this vector holds all the evaluated
    ///               points
    /// \param der_u upon function return, this vector holds all first derivatives
    ///              in the first parameter direction
    /// \param der_v upon function return, this vector holds all first derivatives
    ///              in the second parameter direction
    /// \param der_w upon function return, this vector holds all first derivatives
    ///              in the third parameter direction
    /// \param param_u upon function return, this vector holds all the numerical
    ///                values for the first parameter where evaluation has taken place
    /// \param param_v upon function return, this vector holds all the numerical
    ///                values for the second parameter where evaluation has taken place.
    /// \param param_w upon function return, this vector holds all the numerical
    ///                values for the third parameter where evaluation has taken place.
    /// \param evaluate_from_right specifies directional derivatives, true=right, false=left
    ///                            same behaviour for all parameter directions
    void gridEvaluator (int num_u, int num_v, int num_w,
			std::vector< double > &points,
			std::vector< double > &der_u,
			std::vector< double > &der_v,
			std::vector< double > &der_w,
			std::vector< double > &param_u,
			std::vector< double > &param_v,
			std::vector< double > &param_w,
			bool evaluate_from_right = true) const;

    /// Evaluate points and derivatives on an entire grid, taking computational
    /// advantage over calculating all these values simultaneously rather than
    /// one-by-one.
    /// \param param_u this vector holds all the numerical values for the first parameter 
    ///                where evaluation will tak place
    /// \param param_v this vector holds all the numerical values for the second parameter 
    ///                where evaluation will tak place
    /// \param param_w this vector holds all the numerical values for the third parameter 
    ///                where evaluation will tak place
    /// \param points upon function return, this vector holds all the evaluated
    ///               points
    /// \param der_u upon function return, this vector holds all first derivatives
    ///              in the first parameter direction
    /// \param der_v upon function return, this vector holds all first derivatives
    ///              in the second parameter direction
    /// \param der_w upon function return, this vector holds all first derivatives
    ///              in the third parameter direction
    /// \param evaluate_from_right specifies directional derivatives, true=right, false=left
    ///                            same behaviour for all parameter directions
    void gridEvaluator (const std::vector< double > &param_u,
			const std::vector< double > &param_v,
			const std::vector< double > &param_w,
			std::vector< double > &points,
			std::vector< double > &der_u,
			std::vector< double > &der_v,
			std::vector< double > &der_w,
			bool evaluate_from_right = true) const;

    /// Evaluate points on an entire grid, taking computational advantage over
    /// calculating all these values simultaneously rather than one-by-one.
    /// \param num_u number of values to evaluate along first parameter direction
    /// \param num_v number of values to evaluate along second parameter direction
    /// \param num_w number of values to evaluate along third parameter directon
    /// \param points upon function return, this vector holds all the evaluated
    ///               points
    /// \param param_u upon function return, this vector holds all the numerical
    ///                values for the first parameter where evaluation has taken place
    /// \param param_v upon function return, this vector holds all the numerical
    ///                values for the second parameter where evaluation has taken place.
    /// \param param_w upon function return, this vector holds all the numerical
    ///                values for the third parameter where evaluation has taken place.
    void gridEvaluator (int num_u, int num_v, int num_w,
			std::vector< double > &points,
			std::vector< double > &param_u,
			std::vector< double > &param_v,
			std::vector< double > &param_w) const;

    /// Evaluate points on an entire grid, taking computational advantage over 
    /// calculating all these values simultaneously rather than one-by-one.
    /// \param param_u this vector holds all the numerical values for the first parameter 
    ///                where evaluation will tak place
    /// \param param_v this vector holds all the numerical values for the second parameter 
    ///                where evaluation will tak place
    /// \param param_w this vector holds all the numerical values for the third parameter 
    ///                where evaluation will tak place
    /// \param points upon function return, this vector holds all the evaluated
    ///               points
    void gridEvaluator (const std::vector< double > &param_u,
			const std::vector< double > &param_v,
			const std::vector< double > &param_w,
			std::vector< double > &points) const;

    /// Evaluate positions and first derivatives of all basis values in a given parameter tripple
    /// For non-rationals this is an interface to BsplineBasis::computeBasisValues 
    /// where the basis values in each parameter direction are multiplied to 
    /// compute B_i(u)*B_j(v)*B_k(w), for rationals the routine evaluates the rational
    /// basis functions, i.e. the basis functions are divided by the denominator of the
    /// volume in the given parameter direction
    /// \param param the parameter tripple in which to compute
    /// \param basisValues the value of all basis functions, size equal to 
    ///                    (degree_u+1)*(degree_v+1)*(degree_w+1)
    /// \param basisDerivs_u the derivative of all basis functions, same size as previous
    /// \param basisDerivs_v the derivative of all basis functions, same size as previous
    /// \param basisDerivs_w the derivative of all basis functions, same size as previous
    /// \param evaluate_from_right specifies directional derivatives, true=right, false=left
    void computeBasis(double param[], 
		      std::vector< double > &basisValues,
		      std::vector< double > &basisDerivs_u,
		      std::vector< double > &basisDerivs_v,
		      std::vector< double > &basisDerivs_w,
		      bool evaluate_from_right = true) const;

    /// Evaluate positions of all basis values in a specified
    /// grid.  For non-rationals this is an interface to BsplineBasis::computeBasisValues 
    /// where the basis values in each parameter direction are multiplied to 
    /// compute B_i(u)*B_j(v)*B_k(w), for rationals the routine evaluates the rational 
    /// basis functions, i.e. the basis functions are divided by the denominator of the volume in 
    /// the given parameter direction
    /// \param param_u this vector holds all the numerical values for the first parameter 
    ///                where evaluation will tak place
    /// \param param_v this vector holds all the numerical values for the second parameter 
    ///                where evaluation will tak place
    /// \param param_w this vector holds all the numerical values for the third parameter 
    ///                where evaluation will tak place
    /// \param basisValues the value of all basis functions. The first matrix dimension 
    ///                will correspond to all evaluation parameters, i.e. the size is 
    ///                param_u.size()*param_v*size()*param_w.size(). The second 
    ///                dimension correspond to all coefficients, i.e. 
    ///                nmb_coef_u*nmb_coef_v*nmb_coef_w. Most entries in the matrix
    ///                will be zero
    typedef std::vector<double>  Dvector; // for convenience to shorten
    typedef std::vector<Dvector> Dmatrix; // function argument lists...

    // Size of returned vectors: kk1*kk2*kk3, where kk1 is order_u etc.
    void computeBasis(const std::vector<double>::const_iterator& bas_vals_u,
		      const std::vector<double>::const_iterator& bas_vals_v,
		      const std::vector<double>::const_iterator& bas_vals_w,
		      int left_u,
		      int left_v,
		      int left_w,
		      std::vector<double>& basisValues,
		      std::vector<double>& basisDerivs_u,
		      std::vector<double>& basisDerivs_v,
		      std::vector<double>& basisDerivs_w) const;

    void computeBasisGrid(const Dvector& param_u,
			  const Dvector& param_v,
			  const Dvector& param_w,
			  Dmatrix& basisValues) const; 

    void computeBasis(double param_u,
		      double param_v,
		      double param_w,
		      BasisPts& result) const;

    void computeBasis(double param_u,
		      double param_v,
		      double param_w,
		      BasisDerivs& result,
		      bool evaluate_from_right = true) const;

    void computeBasis(double param_u,
		      double param_v,
		      double param_w,
		      BasisDerivs2& result,
		      bool evaluate_from_right = true) const;

    void computeBasisGrid(const Dvector& param_u,
			  const Dvector& param_v,
			  const Dvector& param_w,
			  std::vector<BasisPts>& result) const; 

    /// Evaluate positions and first derivatives of all basis values in a specified
    /// grid.  For non-rationals this is an interface to BsplineBasis::computeBasisValues 
    /// where the basis values in each parameter direction are multiplied to 
    /// compute B_i(u)*B_j(v)*B_k(w), for rationals the routine evaluates the rational 
    /// basis functions, i.e. the basis functions are divided by the denominator of the volume in 
    /// the given parameter direction
    /// \param param_u this vector holds all the numerical values for the first parameter 
    ///                where evaluation will tak place
    /// \param param_v this vector holds all the numerical values for the second parameter 
    ///                where evaluation will tak place
    /// \param param_w this vector holds all the numerical values for the third parameter 
    ///                where evaluation will tak place
    /// \param basisValues the value of all basis functions. The first matrix dimension 
    ///                will correspond to all evaluation parameters, i.e. the size is 
    ///                param_u.size()*param_v*size()*param_w.size(). The second 
    ///                dimension correspond to all coefficients, i.e. 
    ///                nmb_coef_u*nmb_coef_v*nmb_coef_w. Most entries in the matrix
    ///                will be zero
    /// \param basisDerivs_u the derivative of all basis functions, size as above
    /// \param basisDerivs_v the derivative of all basis functions
    /// \param basisDerivs_w the derivative of all basis functions
    /// \param evaluate_from_right specifies directional derivatives, true=right, false=left

    void computeBasisGrid(const Dvector& param_u,
			  const Dvector& param_v,
			  const Dvector& param_w,
			  Dmatrix& basisValues,
			  Dmatrix& basisDerivs_u,
			  Dmatrix& basisDerivs_v,
			  Dmatrix& basisDerivs_w,
			  bool evaluate_from_right = true) const; 


    void computeBasisGrid(const Dvector& param_u,
			  const Dvector& param_v,
			  const Dvector& param_w,
			  std::vector<BasisDerivs>& result,
			  bool evaluate_from_right = true) const; 

    void computeBasisGrid(const Dvector& param_u,
			  const Dvector& param_v,
			  const Dvector& param_w,
			  std::vector<BasisDerivs2>& result,
			  bool evaluate_from_right = true) const;


    // inherited from ParamVolume
    virtual double nextSegmentVal(int dir, double par, bool forward, double tol) const;

    /// Check if a volume is open, closed or parametric in a given parameter direction
    /// \param pardir The parameter direction in which to check periodicity
    /// \param epsilon Tolerance used in computation
    /// \ return 0 - Open, i.e. multiple knots in end parameters and not closed
    ///          1 - Closed, multiple knots in end parameters 
    ///        > 1 - Periodic knot vector with position and (result-1) derivative equal across the seem
    int volumePeriodicity(int pardir, double epsilon) const;

    /// Check degeneracy
    /// \param tol Tolerance used in test
    /// \param is_degenerate 0 - Not degenerate, 1 - surface degenerate in specified boundary,
    ///                      2 - surface degenerate to line, 3 - surface degenerate to point    
    void checkDegeneracy(double tol, int is_degenerate[]) const;

    /// Specify degeneracy
    /// \param which_sf 0=u_min, 1=u_max, 2=v_min, 3=v_max, 4=w_min, 5=w_max
    /// \param type 0 - Not degenerate, 1 - surface degenerate in specified boundary,
    ///             2 - surface degenerate to line, 3 - surface degenerate to point    
    /// \param b true if surface is degenerate along bottom curve (v=min for surface in parameterized 
    //                                                             u and v)
    /// \param r true if surface is degenerate along right curve
    /// \param t true if surface is degenerate along top curve
    /// \param l true if surface is degenerate along left curve
    /// \param tol Tolerance used in test
    bool isDegenerate(int which_sf, int& type, bool& b, bool& r,
		      bool& t, bool& l, double tol) const;

    /// Check if the volume is of type spline
    virtual bool isSpline() const
    {
      return true;  // This is a spline
    }

    /// Return the spline surface represented by this surface
    virtual SplineVolume* asSplineVolume() 
    {
        // We return a copy of this object, to avoid differing memory handling depending on volume type.
        return clone();
    }

    // inherited from ParamVolume
    virtual void translate(const Point& vec);

    /// Scale the coefficients of the volume with a given factor
    void scale(double fac);

    /// Adds the given deformation vector to the coefficients.
    /// Similar as translate, but with separate vector for each point
    void deform(const std::vector<double>& vec, int vdim = 0);

    /// Add coefficients from another volume. Weights are not summed for rational cases
    /// Nothing is done and exception is raised if
    /// - Spline spaces are different in any paramter direction (order or knot vectors are not identical)
    /// - The geometry spaces of the volumes have different dimension
    /// - One volume is rational while the other is not
    /// - The weights are not equal within a given tolerance (only if volumes are rational)
    /// \param other The other spline volume with coefficients to be added into this volume
    /// \param tol tolerance used to test if weights are considered equal in rational case
    void add(const SplineVolume* other, double tol = 1.0e-10);

    /// Test if the volume is lefthanded. Returns false if volume does
    /// not lie in 3-dimensional space. The test is done by taking the
    /// XXX of the derivatives at the centre point. Assumes volume is
    /// not self intersecting, and has linearly independent derivatives
    /// at the centre.
    bool isLeftHanded();

private:

    /// Degeneracy information regarding one boundary surface of the current spline volume.
    struct degenerate_info
    {
      /// Whether or not the degeneracy information is computed
      bool is_set_;
      /// 0 - Not degenerate, 1 - surface degenerate in specified boundary, 
      /// 2 - surface degenerate to line, 3 - surface degenerate to point
      int  type_;  
      /// Indicates if this boundary surface has degenerate 
      /// boundaries, b_ = bottom (vmin), r_ = right (umax), t_ = top (vmax), 
      // l_ = left (umin)
      bool b_, r_, t_, l_; 
      /// Tolerance used in computations
      double tol_;  
	
      /// Default constructor. The information is not computed
      degenerate_info()
      { is_set_ = false; }
    };

    /// Canonical data
    /// Dimension of the geometry space
    int dim_;
    /// Whether or not the volume is rational
    bool rational_; 
    /// Spline space in first parameter direction
    BsplineBasis basis_u_; 
    /// Spline space in second parameter direction
    BsplineBasis basis_v_; 
    /// Spline space in third parameter direction
    BsplineBasis basis_w_; 
    /// Like ecoef in SISL
    std::vector<double> coefs_; 
    /// Like rcoef in SISL, only used if rational
    std::vector<double> rcoefs_;

    /// Generated data
    #define SPLINE_VOLUME_BD_SFS_SIZE 6
    mutable shared_ptr<SplineSurface> bd_sfs_[SPLINE_VOLUME_BD_SFS_SIZE];  // Sequence: u_min, u_max, v_min, v_max, 
                                                                                  // w_min, w_max
    mutable std::pair<int,double> periodicity_info_[SPLINE_VOLUME_PERIOD_INFO_SIZE ];   // Initalize first to -1 means not checked,
                                                                                        // second is tolerance used in computation
    #define SPLINE_VOLUME_DEGEN_SIZE 6
    mutable degenerate_info degen_[SPLINE_VOLUME_DEGEN_SIZE];    // For each surface, sequence as before


    /// Helper functions
    void updateCoefsFromRcoefs();
    std::vector<double>& activeCoefs() { return rational_ ? rcoefs_ : coefs_; }

    void getPlaneNormals(Point pnt, Point vec, Point& norm1, Point& norm2) const;

    void pointsGrid(const std::vector< double > &param_u,
		    const std::vector< double > &param_v,
		    const std::vector< double > &param_w,
		    int derivs,
		    std::vector< double > &points,
		    bool evaluate_from_right = true) const;

 public:
    void pointsGrid(int numu, int numv, int numw,
		    std::vector<double>::iterator basisvals_u,
		    std::vector<double>::iterator basisvals_v,
		    std::vector<double>::iterator basisvals_w,
		    std::vector<int>::iterator knotinter_u,
		    std::vector<int>::iterator knotinter_v,
		    std::vector<int>::iterator knotinter_w,
		    int derivs,
		    std::vector< double > &points) const;

 private:
    void accumulateBasis(double* basisvals_u, int uorder,
			 double* basisvals_v, int vorder,
			 double* basisvals_w, int worder,
			 double* weights, 
			 double* basisValues) const;
			 
    void accumulateBasis(const double* basisvals_u, int uorder,
			 const double* basisvals_v, int vorder,
			 const double* basisvals_w, int worder,
			 const double* weights, 
			 double* basisValues,
			 double* basisDerivs_u,
			 double* basisDerivs_v,
			 double* basisDerivs_w) const;
			 
    void accumulateBasis(double* basisvals_u, int uorder,
			 double* basisvals_v, int vorder,
			 double* basisvals_w, int worder,
			 double* weights, 
			 double* basisValues,
			 double* basisDerivs_u,
			 double* basisDerivs_v,
			 double* basisDerivs_w,
			 double* basisDerivs_uu,
			 double* basisDerivs_uv,
			 double* basisDerivs_uw,
			 double* basisDerivs_vv,
			 double* basisDerivs_vw,
			 double* basisDerivs_ww) const;

    void raiseOrder_wdir(int raise);

    void getSeed(const Point& pt, double par[]) const;
			 

};


} // namespace Go




#endif // _SPLINEVOLUME_H

