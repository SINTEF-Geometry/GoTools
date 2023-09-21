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

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineUtils.h"
#include <array>

using namespace std;

namespace Go
{
  namespace
  {
    /// Functor that scales the input argument.
    class ScaleBy// : public std::unary_function<double, double>
    {
      double m_scale;

    public:
      ScaleBy(const double& scale) : m_scale(scale) {}

      double operator()(const double& value) { return m_scale * value; }
    };
  } // anonymous namespace

//===========================================================================
// this is an attempt to optimize the point evaluation routine.  
// If it is not satisfactory for some reason, the original routine exists commented-out
// below
void SplineSurface::point(Point& result, double upar, double vpar) const
//===========================================================================
{
    result.resize(dim_);
    const int uorder = order_u();
    const int vorder = order_v();
    const int unum = numCoefs_u();
    int kdim = rational_ ? dim_ + 1 : dim_;

#ifdef _OPENMP
    ScratchVect<double, 10> Bu(uorder);
    ScratchVect<double, 10> Bv(vorder);
    ScratchVect<double, 4> tempPt(kdim);
    ScratchVect<double, 4> tempResult(kdim);
#else
    static ScratchVect<double, 10> Bu(uorder);
    static ScratchVect<double, 10> Bv(vorder);
    static ScratchVect<double, 4> tempPt(kdim);
    static ScratchVect<double, 4> tempResult(kdim);
#endif

    Bu.resize(uorder);
    Bv.resize(vorder);
    tempPt.resize(kdim);
    tempResult.resize(kdim);

    // compute tbe basis values and get some data about the spline spaces
    basis_u_.computeBasisValues(upar, Bu.begin());
    basis_v_.computeBasisValues(vpar, Bv.begin());
    const int uleft = basis_u_.lastKnotInterval();
    const int vleft = basis_v_.lastKnotInterval();
    
    // compute the tensor product value
    const int start_ix =  (uleft - uorder + 1 + unum * (vleft - vorder + 1)) * kdim;

    register double* ptemp;
    register const double* co_ptr = rational_ ? &rcoefs_[start_ix] : &coefs_[start_ix];
    fill(tempResult.begin(), tempResult.end(), double(0));

    for (register double* bval_v_ptr = Bv.begin(); bval_v_ptr != Bv.end(); ++bval_v_ptr) {
	register const double bval_v = *bval_v_ptr;
	fill(tempPt.begin(), tempPt.end(), 0);
	for (register double* bval_u_ptr = Bu.begin(); bval_u_ptr != Bu.end(); ++bval_u_ptr) {
	    register const double bval_u = *bval_u_ptr;
	    for (ptemp = tempPt.begin(); ptemp != tempPt.end(); ++ptemp) {
		*ptemp += bval_u * (*co_ptr++);
	    }
	}
	ptemp = tempPt.begin();
	for (register double* p = tempResult.begin(); p != tempResult.end(); ++p) {
	    *p += (*ptemp++) * bval_v;
	}
	co_ptr += kdim * (unum - uorder);
    }

    copy(tempResult.begin(), tempResult.begin() + dim_, result.begin());
    if (rational_) {
	const double w_inv = double(1) / tempResult[kdim - 1];
	transform(result.begin(), result.end(), result.begin(), ScaleBy(w_inv));
    }
}

// this was the original point evaluation routine.

// //===========================================================================
// void SplineSurface::point(Point& result, double upar, double vpar) const
// //===========================================================================
//  {
//     if (result.dimension() != dim_)
// 	result.resize(dim_);

//     // Take care of the rational case
//     const std::vector<double>& co = rational_ ? rcoefs_ : coefs_;
//     int kdim = dim_ + (rational_ ? 1 : 0);

//     // Make temporary storage for the basis values and a temporary
//     // computation cache.
//     Go::ScratchVect<double, 10> b0(basis_u_.order());
//     Go::ScratchVect<double, 10> b1(basis_v_.order());
//     Go::ScratchVect<double, 4> temp(kdim);
//     Go::ScratchVect<double, 4> restemp(kdim);
//     std::fill(restemp.begin(), restemp.end(), 0.0);

//     // Compute the basis values and get some data about the spline spaces
//     basis_u_.computeBasisValues(upar, &b0[0]);
//     int uleft = basis_u_.lastKnotInterval();
//     int uorder = basis_u_.order();
//     int unum = basis_u_.numCoefs();
//     basis_v_.computeBasisValues(vpar, &b1[0]);
//     int vleft = basis_v_.lastKnotInterval();
//     int vorder = basis_v_.order();

//     // Compute the tensor product value
//     int coefind = uleft-uorder+1 + unum*(vleft-vorder+1);
//     for (int jj = 0; jj < vorder; ++jj) {
// 	std::fill(temp.begin(), temp.end(), 0.0);
// 	for (int ii = 0; ii < uorder; ++ii) {
// 	    for (int dd = 0; dd < kdim; ++dd) {
// 		temp[dd] += b0[ii] * co[coefind*kdim + dd];
// 	    }
// 	    coefind += 1;
// 	}
// 	for (int dd = 0; dd < kdim; ++dd)
// 	    restemp[dd] += temp[dd]*b1[jj];
// 	coefind += unum - uorder;
//     }

//     // Copy from restemp to result
//     if (rational_) {
// 	for (int dd = 0; dd < dim_; ++dd) {
// 	    result[dd] = restemp[dd]/restemp[kdim-1];
// 	}
//     } else {
// 	for (int dd = 0; dd < dim_; ++dd) {
// 	    result[dd] = restemp[dd];
// 	}
//     }
// }


//===========================================================================
void
SplineSurface::point(std::vector<Point>& result, double upar, double vpar,
		     int derivs, bool u_from_right, bool v_from_right,
		     double resolution) const
//===========================================================================
{
    DEBUG_ERROR_IF(derivs < 0, "Negative number of derivatives makes no sense.");
    int totpts = (derivs + 1)*(derivs + 2)/2;
    DEBUG_ERROR_IF((int)result.size() < totpts, "The vector of points must have sufficient size.");

    for (int i = 0; i < totpts; ++i) {
	if (result[i].dimension() != dim_) {
	    result[i].resize(dim_);
	}
    }

    if (derivs == 0) {
	point(result[0], upar, vpar);
	return;
    }

    // Take care of the rational case
    const std::vector<double>& co = rational_ ? rcoefs_ : coefs_;
    int kdim = dim_ + (rational_ ? 1 : 0);

    // Make temporary storage for the basis values and a temporary
    // computation cache.
    Go::ScratchVect<double, 30> b0(basis_u_.order() * (derivs+1));
    Go::ScratchVect<double, 30> b1(basis_v_.order() * (derivs+1));
    Go::ScratchVect<double, 30> temp(kdim * totpts);
    Go::ScratchVect<double, 30> restemp(kdim * totpts);
    std::fill(restemp.begin(), restemp.end(), 0.0);
    // Compute the basis values and get some data about the spline spaces
    if (u_from_right) {
	basis_u_.computeBasisValues(upar, &b0[0], derivs, resolution);
    } else {
	basis_u_.computeBasisValuesLeft(upar, &b0[0], derivs, resolution);
    }
    int uleft = basis_u_.lastKnotInterval();
    int uorder = basis_u_.order();
    int unum = basis_u_.numCoefs();
    if (v_from_right) {
	basis_v_.computeBasisValues(vpar, &b1[0], derivs, resolution);
    } else {
	basis_v_.computeBasisValuesLeft(vpar, &b1[0], derivs, resolution);
    }
    int vleft = basis_v_.lastKnotInterval();
    int vorder = basis_v_.order();
    // Compute the tensor product value
    int coefind = uleft-uorder+1 + unum*(vleft-vorder+1);
    int derivs_plus1=derivs+1;
    for (int jj = 0; jj < vorder; ++jj) {
      int jjd=jj*(derivs_plus1);
	std::fill(temp.begin(), temp.end(), 0.0);
		
	for (int ii = 0; ii < uorder; ++ii) {
	  int iid=ii*(derivs_plus1);
	  const double *co_p=&co[coefind*kdim];
	  for (int dd = 0; dd < kdim; ++dd,++co_p) {
	      int temp_ind=dd;
		for (int vder = 0; vder < derivs_plus1; ++vder) {
		    for (int uder = 0; uder < vder+1; ++uder) {
			temp[temp_ind]
			  += b0[iid+vder - uder]*(*co_p);
			temp_ind+=kdim;
		    }
		}
	    }
	    coefind += 1;
	}

	for (int dd = 0; dd < kdim; ++dd) {
	    int dercount = 0;
	    for (int vder = 0; vder < derivs_plus1; ++vder) {
		for (int uder = 0; uder < vder + 1; ++uder) {
		    restemp[dercount*kdim + dd] 
			+= temp[dercount*kdim + dd]*b1[uder + jjd];
		    ++dercount;
		}
	    }
	}

	coefind += unum - uorder;
    }
    // Copy from restemp to result
    if (rational_) {
	std::vector<double> restemp2(totpts*dim_);
	SplineUtils::surface_ratder(&restemp[0], dim_, derivs, &restemp2[0]);
	for (int i = 0; i < totpts; ++i) {
	    for (int dd = 0; dd < dim_; ++dd) {
		result[i][dd] = restemp2[i*dim_ + dd];
	    }
	}
    } else {
      double* restemp_it=restemp.begin();
	for (int i = 0; i < totpts; ++i) {
	    for (int dd = 0; dd < dim_; ++dd) {
		result[i][dd] = *restemp_it;
		++restemp_it;
	    }
	}
    }
        
}

#define NOT_FINISHED_YET
#ifdef NOT_FINISHED_YET
//===========================================================================
void SplineSurface::pointsGrid(int m1, int m2, int derivs,
				 const double* basisvals1,
				 const double* basisvals2,
				 const int* knotint1,
				 const int* knotint2,
				 double* result,
				 double* normals) const
//===========================================================================
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Evaluate this surface over an m1 * m2
*              grid, assuming that the B-splines have been
*              pre-evaluated, by s1504, over that grid.
*              The knots et1 and et2 and grid points (x[i],y[j]) are not needed.
*              Compute ider derivatives.
*
* INPUT      : ps1    - Pointer to the surface to evaluate.
*              ider   - Number of derivatives to calculate.
*                       < 0 : No derivative calculated.
*                       = 0 : Position calculated.
*                       = 1 : Position and first derivative calculated.
*                       etc.
*              m1     - Number of grid points in first direction.
*              m2     - Number of grid points in first direction.
*              ileft1 - Array of indexes to the intervals in the knotvector
*                       in the first parameter direction where each subsequence
*                       of k1 non-zero B-splines are located.
*              ileft2 - Array of indexes to the intervals in the knotvector
*                       in the second parameter direction where each subsequence
*                       of k2 non-zero B-splines are located.
*              ebder1 - Triple array of dimension [(ider+1)*k1*m1] containing
*                       values of the k1 nonzero B-splines and their
*                       derivatives at the points x[0],...,x[m1-1]
*                       (i.e. pre-evaluated B-splines).
*                       These numbers are stored in the following order:
*                       First the (ider+1) derivatives of the first nonzero
*                       B-spline at x[0]. Then the (ider+1) derivatives of
*                       the second nonzero B-spline at x[0], etc.
*                       Later we repeat for x[1],... etc.
*
*              ebder2 - Triple array of dimension [(ider+1)*k2*m2] containing
*                       values of the k2 nonzero B-splines and their
*                       derivatives at the points y[0],...,y[m2-1]
*
* OUTPUT     : eder   - Array where the derivatives of the surface
*                       are placed, dimension
*                         idim * ((ider+1)(ider+2) / 2) * m1 * m2.
*                       The sequence is position,
*                       first derivative in first parameter direction,
*                       first derivative in second parameter direction,
*                       (2,0) derivative, (1,1) derivative, (0,2)
*                       derivative, etc. at point (x[0],y[0]),
*                       followed by the same information at (x[1],y[0]),
*                       etc.
*              norm   - Normals of surface. Is calculated if ider >= 1.
*                       Dimension is idim*m1*m2.
*                       The normals are normalized. @afr changed
*              jstat  - status messages
*                          = 2      : Surface is degenerate
*                                     at some point, normal
*                                     has zero length.
*                          = 1      : Surface is close to
*                                     degenerate at some point
*                                     Angle between tangents,
*                                     less than angular tolerance.
*                          = 0      : ok
*                          < 0      : error
*
* METHOD     : The code and method is similar to that of s1421 except that
*              the B-splines and their derivatives have already been
*              calculated (and are given in ebder1 and ebder2) and
*              we evaluate over a whole grid rather than at one point.
*              The method is to find the control points and control derivatives
*              of each isocurve in x (fixed y value). We then
*              evaluate at each x value along the isocurve.
*
* CALLS      : s6err     - errormacros.handling routine
*              s6strider - Make derivative of rational expression
*
* WRITTEN BY : Michael Floater, SINTEF, May 1998.
*********************************************************************
*/
{
    const double* ebder1 = basisvals1;
    const double* ebder2 = basisvals2;
    double* eder = result;
    double* norm = normals;
    const int* ileft1 = knotint1;
    const int* ileft2 = knotint2;
    int ider = derivs;
    int i1,i2;          /* Loop variables. */
    int i1pos,i2pos;    /* Offset indexes. */
    int kn1,kn2;        /* The number of B-splines accociated with the knot
			   vectors st1 and st2.                            */
    int kk1,kk2;        /* The polynomial order of the surface in the two
			   directions.                                     */
    int kdim;           /* The dimension of the space in which the surface
			   lies. Equivalently, the number of components
			   of each B-spline coefficient.                   */
    int kleft1,kleft2;  /* Local versions of knot intervals.            */
    int ki,kx,kjh;      /* Control variables in for loops and for stepping
			   through arrays.                                 */
    int kih2;           /* Index for stepping through ebder2. */
    int ky,kl,kl1,kl2;  /* Control variables in for loops and for stepping
			   through arrays.                                 */
    const double* scoef;         /* The B-spline coefficients of the surface.
			   This is an array of dimension [kn2*kn1*kdim].   */
    double tt;          /* Dummy variable used for holding an array element
			   in a for loop.                                  */
    int size;           /* Space occupied by points and derivs at one eval. */
    int sizeh;          /* Space occupied by homogeneous points and derivs . */
    int size1,size2;    /* Useful variables. */
    /*int ederpos; */      /* Index of position in eder. */
    /*int normpos; */      /* Index of position in norm. */
    
    int knumb2;         /* Necessary size of ew */
    
    int tot,temp;       /* Temporary variables. */



    kn1 = basis_u_.numCoefs();
    kn2 = basis_v_.numCoefs();
    kk1 = basis_u_.order();
    kk2 = basis_v_.order();
    
    kdim = dim_;

    if (rational_) {
	scoef = &rcoefs_[0];
	kdim +=1;
    } else {
	scoef = &coefs_[0];
    }

//  @ afr: The following should have been declared outside this func.
//    Dvec eder(dimension() * ((ider+1)*(ider+2) / 2) * m1 * m2);
//    Dvec norm(dimension() * m1 * m2);

    sizeh = kdim*(ider+1)*(ider+2)/2;
    double* sder = new double[sizeh];
    double* enorm = new double[dim_];

    size = dim_*(ider+1)*(ider+2)/2;
    //    int offset = m1*m2*kdim;
    size1 = (ider+1)*kk1;
    size2 = (ider+1)*kk2;

    // Allocate space for B-spline values and derivatives and one work array.

    knumb2 = kn1*(ider+1)*kdim;
    double* ew = new double[knumb2];
    //    double* ew_iterator = ew; // Iterator used for clearing ew.

    double*  norm_iterator = norm;
    double*  eder_iterator = eder;

    // Run through grid points in the y direction.
    for(i2=0, i2pos=0; i2<m2; i2++, i2pos += size2) {

	kleft2 = ileft2[i2];

	/* Compute the control points (and control derivatives
	   of the v = x[i2] isocurve. */

	/* Set all the elements of ew to 0. */
	double* ew_iterator;
	for (ew_iterator=ew; ew_iterator<ew+knumb2; ew_iterator++)
	    *ew_iterator = 0;

	/* ki steps through the appropriate kk2 rows of B-spline coefficients
	   while kih2 steps through the B-spline value and derivatives for the
	   B-spline given by ki.            */

	kih2 = i2pos;

	for (ki=kleft2-kk2+1; ki<=kleft2; ki++) {
	    /* kx counts through the ider+1 derivatives to be computed.
	       kjh steps through ew once for each ki to accumulate the
	       contribution from the different B-splines.
	       kl1 points to the first component of the first B-spline
	       coefficient in row no. ki of the B-spline coefficient matrix.
	       */

	    kjh = 0; kl1 = kdim*ki*kn1;
	    for (kx=0; kx<=ider; kx++) {
		/* The value of the B-spline derivative is stored in tt while
		   kl2 steps through the kdim components of all the B-spline
		   coefficients that multiply nonzero B-splines along st1.
		   */

		tt = ebder2[kih2++]; kl2 = kl1;
		for (kl=0; kl<kdim*kn1; kl++,kjh++,kl2++) {
		    ew[kjh] += scoef[kl2]*tt;
		}
	    }
	}

	/* Run through grid points in the x direction
	   evaluating along the iso-curve defined by y = y[i2]. */
	for(i1=0, i1pos=0; i1<m1; i1++, i1pos += size1) {
	    kleft1 = ileft1[i1];

	    /* Set all the elements of sder to 0. */
	    for(ki=0; ki<sizeh; ki++) sder[ki] = 0;

	    for(ky=0; ky<=ider; ky++) {
		kl1 = kdim * (ky * kn1 + kleft1 - kk1 + 1);
		for(kx=0; kx<=ider-ky; kx++) {
		    tot = kx + ky;
		    temp = ((tot * (tot+1)) >> 1) + ky;
		    kjh = temp * kdim;

		    for(ki=0; ki<kk1; ki++) {
			tt = ebder1[i1pos + (ider+1) * ki + kx];
			for(kl=0; kl<kdim; kl++) {
			    sder[kjh+kl] += ew[kl1 + kdim * ki + kl] * tt;
			}
		    }
		}
	    }

	    /* If rational surface calculate the derivatives based on
	       derivatives in homogenous coordinates */
      
//  	    if (kind_ == 2 || kind_ == 4) {
//  //	        @ this must be implemented!  
//  //	        s6strider(sder,ps1->idim,ider,eder+ederpos,&kstat);
//  	    } else {
//  		for (ki=0; ki<sizeh; ki++) { 
//  		    *(eder_iterator+ki) = sder[ki];
//  		}
//  	    }
	    if (rational_) {
		SplineUtils::surface_ratder(sder,dim_, ider, eder_iterator);
//	        s6strider(sder,dim_,eder_iterator,&kstat);
	    } else {
		for (ki=0; ki<sizeh; ki++) { 
		    *(eder_iterator+ki) = sder[ki];
		}
	    }
      
	    /* Calculate normal if idim==3 and ider>0. */

	    if (ider>0 && kdim ==3) {
//		enorm = GoCrossProduct(eder_iterator + 3,
//				       eder_iterator + 6);
//		GoNormalize(enorm);
//		copy(enorm,enorm+3,norm_iterator);
		double* i1 = eder + 3;
		double* i2 = eder + 6;
		norm_iterator[0] = i1[1]*i2[2] - i1[2]*i2[1];
		norm_iterator[1] = i1[2]*i2[0] - i1[0]*i2[2];
		norm_iterator[2] = i1[0]*i2[1] - i1[1]*i2[0];
		double ssum = norm_iterator[0]*norm_iterator[0]
		    + norm_iterator[1]*norm_iterator[1]
		    + norm_iterator[2]*norm_iterator[2];		
		double invl = 1/sqrt(ssum);
		norm_iterator[0] = norm_iterator[0]*invl;
		norm_iterator[1] = norm_iterator[1]*invl;
		norm_iterator[2] = norm_iterator[2]*invl;
	    }

	    eder_iterator += size;
	    norm_iterator += kdim;
	}

    }


  /* Free memory. */

    delete [] ew;
    delete [] sder;
    delete [] enorm;
    


  /* Not enough memory. */

  /* kdim less than 1. */

  /* Polynomial order less than 1. */

  /* Fewer B-splines than the order. */

  /* Illegal derivative requested. */

  /* Error in lower level routine.  */

  return;
}
#endif


  void SplineSurface::evalGrid(int num_u, int num_v, 
		double umin, double umax, 
		double vmin, double vmax,
		std::vector<double>& points,
		double nodata_val) const
  {
    vector<double> param_u;
    vector<double> param_v;
    gridEvaluator(num_u, num_v, points, param_u, param_v,
		  umin, umax, vmin, vmax);
  }

void SplineSurface::gridEvaluator(int num_u, int num_v,
				  std::vector<double>& points,
				  std::vector<double>& normals,
				  std::vector<double>& param_u,
				  std::vector<double>& param_v,
				  bool normalize) const
{
    ASSERT(dimension() == 3);
    ASSERT(num_u > 1 && num_v > 1);
    int kdim = rational_ ? dim_+1 : dim_;
    int ncoef_u = numCoefs_u();
    // int ncoef_v = numCoefs_v();
    int ord_u = order_u();
    int ord_v = order_v();

    const double start_u = startparam_u();
    const double start_v = startparam_v();
    const double diff_u = endparam_u() - start_u;
    const double diff_v = endparam_v() - start_v;
    const double du=diff_u/(num_u-1);
    const double dv=diff_v/(num_v-1);

    param_u.resize(num_u);
    param_v.resize(num_v);

    for(int step_u = 0; step_u < num_u; ++step_u)
	param_u[step_u] = start_u+double(step_u)*du;
    for(int step_v = 0; step_v < num_v; ++step_v)
	param_v[step_v] = start_v+double(step_v)*dv;

    vector<BasisDerivsSf> basis_func;

    computeBasisGrid(param_u, param_v, basis_func);
    vector<double> der_u, der_v;
    points.resize(0);
    normals.resize(0);

    vector<double> cf = rational_ ? rcoefs_ : coefs_;

    double inv_weight = 1.0;
    vector<BasisDerivsSf>::const_iterator it = basis_func.begin();
    for (int j = 0; j < num_v; ++j)
      {
	int coef_pos_v = ((*it).left_idx[1]-ord_v+1) * ncoef_u * kdim;
	for (int i = 0; i < num_u; ++i, ++it)
	  {
	    int coef_pos = coef_pos_v + ((*it).left_idx[0]-ord_u+1) * kdim;
	    Point pt_coord(0.0, 0.0, 0.0);
	    Point du_coord(0.0, 0.0, 0.0);
	    Point dv_coord(0.0, 0.0, 0.0);
	    int bas_pos = 0;
	    // double denom = 0.0;
	    for (int bas_v = 0; bas_v < ord_v; ++bas_v)
	      {
		for (int bas_u = 0; bas_u < ord_u; ++bas_u, ++bas_pos)
		  {
		    if (rational_)
		      inv_weight = 1.0 / cf[coef_pos + dim_];
		    for (int k = 0; k < dim_; ++k, ++coef_pos)
		      {
			double current_coef = cf[coef_pos];
			pt_coord[k] += current_coef * (*it).basisValues[bas_pos] * inv_weight;
			du_coord[k] += current_coef * (*it).basisDerivs_u[bas_pos] * inv_weight;
			dv_coord[k] += current_coef * (*it).basisDerivs_v[bas_pos] * inv_weight;
		      }
		    if (rational_)
		      ++coef_pos;
		  }
		coef_pos += (ncoef_u - ord_u)*kdim;
	      }
	    Point norm = du_coord % dv_coord;
	    if (normalize)
	      norm.normalize();
	    for (int k = 0; k < 3; ++k)
	      {
		points.push_back(pt_coord[k]);
		normals.push_back(norm[k]);
	      }
	  }
      }
}

// Same as above, but no normals.
void SplineSurface::gridEvaluator(int num_u, int num_v,
				  std::vector<double>& points,
				  std::vector<double>& param_u,
				  std::vector<double>& param_v,
				  double start_u, double end_u,
				  double start_v, double end_v) const
{
  //    ASSERT(dimension() == 3);
    ASSERT(num_u > 1 && num_v > 1);
    const double diff_u = end_u - start_u;
    const double diff_v = end_v - start_v;
    const double du=diff_u/(num_u-1);
    const double dv=diff_v/(num_v-1);

    param_u.resize(num_u);
    param_v.resize(num_v);

    for(int step_u = 0; step_u < num_u; ++step_u)
	param_u[step_u] = start_u+double(step_u)*du;
    for(int step_v = 0; step_v < num_v; ++step_v)
	param_v[step_v] = start_v+double(step_v)*dv;

    gridEvaluator(points,param_u,param_v);
}

// Added as separate method by kmo for usage in ICADA project.
void SplineSurface::gridEvaluator(std::vector<double>& points,
				  const std::vector<double>& param_u,
				  const std::vector<double>& param_v) const
{
    int num_u = (int)param_u.size();
    int num_v = (int)param_v.size();
    vector<double> basisvals_u(num_u * basis_u_.order());
    vector<double> basisvals_v(num_v * basis_v_.order());
    vector<int>    knotinter_u(num_u * basis_u_.order());
    vector<int>    knotinter_v(num_v * basis_v_.order());

    basis_u_.computeBasisValues(&param_u[0], &param_u[0]+param_u.size(),
				&basisvals_u[0], &knotinter_u[0], 0);
    basis_v_.computeBasisValues(&param_v[0], &param_v[0]+param_v.size(),
				&basisvals_v[0], &knotinter_v[0], 0);

    points.resize(num_u * num_v * dim_);
    pointsGrid(num_u, num_v, 0,
	       &basisvals_u[0], &basisvals_v[0],
	       &knotinter_u[0], &knotinter_v[0],
	       &points[0], 0);
}


//===========================================================================
void SplineSurface::gridEvaluator(const vector<double>& params_u,
				  const vector<double>& params_v,
				  vector<double>& points,
				  vector<double>& derivs_u,
				  vector<double>& derivs_v,
				  bool evaluate_from_right) const
//===========================================================================
{
  int kdim = rational_ ? dim_+1 : dim_;
  int num_u = (int)params_u.size();
  int num_v = (int)params_v.size();
  int ncoef_u = numCoefs_u();
  int ord_u = order_u();
  int ord_v = order_v();

  vector<BasisDerivsSf> basis_func;
  computeBasisGrid(params_u, params_v, basis_func, evaluate_from_right);

  points.resize(0);
  derivs_u.resize(0);
  derivs_v.resize(0);

  vector<double>::const_iterator cf_it = ctrl_begin();
  double inv_weight = 1.0;
  vector<BasisDerivsSf>::const_iterator it = basis_func.begin();

  for (int j = 0; j < num_v; ++j)
    {
      int coef_pos_v = ((*it).left_idx[1]-ord_v+1) * ncoef_u * kdim;
      for (int i = 0; i < num_u; ++i, ++it)
	{
	  int coef_pos = coef_pos_v + ((*it).left_idx[0]-ord_u+1) * kdim;
	  vector<double> pt_coord(dim_, 0.0);
	  vector<double> du_coord(dim_, 0.0);
	  vector<double> dv_coord(dim_, 0.0);
	  int bas_pos = 0;

	  for (int bas_v = 0; bas_v < ord_v; ++bas_v)
	    {
	      for (int bas_u = 0; bas_u < ord_u; ++bas_u, ++bas_pos)
		{
		  if (rational_)
		    inv_weight = 1.0 / cf_it[coef_pos + dim_];
		  for (int k = 0; k < dim_; ++k, ++coef_pos)
		    {
		      double current_coef = cf_it[coef_pos];
		      pt_coord[k] += current_coef * (*it).basisValues[bas_pos] * inv_weight;
		      du_coord[k] += current_coef * (*it).basisDerivs_u[bas_pos] * inv_weight;
		      dv_coord[k] += current_coef * (*it).basisDerivs_v[bas_pos] * inv_weight;
		    }
		  if (rational_)
		    ++coef_pos;
		}
	      coef_pos += (ncoef_u - ord_u)*kdim;
	    }
	  for (int k = 0; k < dim_; ++k)
	    {
	      points.push_back(pt_coord[k]);
	      derivs_u.push_back(du_coord[k]);
	      derivs_v.push_back(dv_coord[k]);
	    }
	}
    }
}


// This thingie is intended to copy and adjust the existing
//  pointsGrid() method, so it would return correctly formatted results.
void SplineSurface::pointsGridNoDerivs(int m1, int m2,
				       const double* basisvals1,
				       const double* basisvals2,
				       const int* knotint1,
				       const int* knotint2,
				       double* result,
				       double* normals,
				       bool normalize) const
{
    int i1,i2;          /* Loop variables. */
    int i1pos,i2pos;    /* Offset indexes. */
    int kn1,kn2;        /* The number of B-splines accociated with the knot
			   vectors st1 and st2.                            */
    int kk1,kk2;        /* The polynomial order of the surface in the two
			   directions.                                     */
    int kleft1,kleft2;  /* Local versions of knot intervals.            */
    int ki,kjh;         /* Control variables in for loops and for stepping
			   through arrays.                                 */
    int kih2;           /* Index for stepping through basisvals2. */
    int kl,kl1,kl2;     /* Control variables in for loops and for stepping
			   through arrays.                                 */
    const double* scoef;         /* The B-spline coefficients of the surface.
			   This is an array of dimension [kn2*kn1*kdim].   */
    double tt;          /* Dummy variable used for holding an array element
			   in a for loop.                                  */
    int size1,size2;    /* Useful variables. */
    
    int knumb2;         /* Necessary size of ew */
    



    kn1 = basis_u_.numCoefs();
    kn2 = basis_v_.numCoefs();
    kk1 = basis_u_.order();
    kk2 = basis_v_.order();
    
    // Check that dim_=3 && rational_=false;

    scoef = &coefs_[0];


    ScratchVect<double, 9> sder(9);
    ScratchVect<double, 3> enorm(3);

    //    int offset = m1*m2*9;
    size1 = 2*kk1;
    size2 = 2*kk2;

    // Allocate space for B-spline values and derivatives and one work array.

    knumb2 = kn1*6;
    ScratchVect<double, 300> ew(knumb2);

    double*  normals_iterator = normals;
    double*  result_iterator = result;

    // Run through grid points in the y direction.
    for(i2=0, i2pos=0; i2<m2; i2++, i2pos += size2) {

	kleft2 = knotint2[i2];

	/* Compute the control points (and control derivatives
	   of the v = x[i2] isocurve. */

	/* Set all the elements of ew to 0. */
	for (int i = 0; i<knumb2; i++)
	    ew[i] = 0;

	/* ki steps through the appropriate kk2 rows of B-spline coefficients
	   while kih2 steps through the B-spline value and derivatives for the
	   B-spline given by ki.            */

	kih2 = i2pos;

	for (ki=kleft2-kk2+1; ki<=kleft2; ki++) {
	    /* kx counts through the ider+1 derivatives to be computed.
	       kjh steps through ew once for each ki to accumulate the
	       contribution from the different B-splines.
	       kl1 points to the first component of the first B-spline
	       coefficient in row no. ki of the B-spline coefficient matrix.
	       */

	    kjh = 0; kl1 = 3*ki*kn1;
		/* The value of the B-spline derivative is stored in tt while
		   kl2 steps through the kdim components of all the B-spline
		   coefficients that multiply nonzero B-splines along st1.
		   */
	    tt = basisvals2[kih2++]; kl2 = kl1;
	    for (kl=0; kl<3*kn1; kl++,kjh++,kl2++) {
		ew[kjh] += scoef[kl2]*tt;
	    }
	    tt = basisvals2[kih2++]; kl2 = kl1;
	    for (kl=0; kl<3*kn1; kl++,kjh++,kl2++) {
		ew[kjh] += scoef[kl2]*tt;
	    }
	}

	/* Run through grid points in the x direction
	   evaluating along the iso-curve defined by y = y[i2]. */
	for(i1=0, i1pos=0; i1<m1; i1++, i1pos += size1) {
	    kleft1 = knotint1[i1];

	    /* Set all the elements of sder to 0. */
	    sder[0] = sder[1] = sder[2] = sder[3] = sder[4] = 0;
	    sder[5] = sder[6] = sder[7] = sder[8] = 0;

	    kl1 = 3 * (kleft1 - kk1 + 1);
	    for(ki=0; ki<kk1; ki++) {
		tt = basisvals1[i1pos + ki*2];
		sder[0] += ew[kl1 + ki*3 + 0] * tt;
		sder[1] += ew[kl1 + ki*3 + 1] * tt;
		sder[2] += ew[kl1 + ki*3 + 2] * tt;
	    }
	    for(ki=0; ki<kk1; ki++) {
		tt = basisvals1[i1pos + ki*2 + 1];
		sder[3] += ew[kl1 + ki*3 + 0] * tt;
		sder[4] += ew[kl1 + ki*3 + 1] * tt;
		sder[5] += ew[kl1 + ki*3 + 2] * tt;
	    }

	    kl1 = 3 * (kn1 + kleft1 - kk1 + 1);
	    for(ki=0; ki<kk1; ki++) {
		tt = basisvals1[i1pos + ki*2];
		sder[6] += ew[kl1 + ki*3 + 0] * tt;
		sder[7] += ew[kl1 + ki*3 + 1] * tt;
		sder[8] += ew[kl1 + ki*3 + 2] * tt;
	    }
	    kjh = 6;

	    *(result_iterator++) = sder[0];
	    *(result_iterator++) = sder[1];
	    *(result_iterator++) = sder[2];
      
	    double* i1 = &sder[3];
	    double* i2 = &sder[6];
	    normals_iterator[0] = i1[1]*i2[2] - i1[2]*i2[1];
	    normals_iterator[1] = i1[2]*i2[0] - i1[0]*i2[2];
	    normals_iterator[2] = i1[0]*i2[1] - i1[1]*i2[0];
	    if (normalize)
	      {
		double ssum = normals_iterator[0]*normals_iterator[0]
		  + normals_iterator[1]*normals_iterator[1]
		  + normals_iterator[2]*normals_iterator[2];		
		double invl = 1/sqrt(ssum);
		normals_iterator[0] = normals_iterator[0]*invl;
		normals_iterator[1] = normals_iterator[1]*invl;
		normals_iterator[2] = normals_iterator[2]*invl;
	      }
	    normals_iterator += 3;
	}

    }
    
  return;
}


//===========================================================================
void SplineSurface::computeBasis(double param[], 
				 std::vector< double > &basisValues,
				 std::vector< double > &basisDerivs_u,
				 std::vector< double > &basisDerivs_v,
				 bool evaluate_from_right) const
//===========================================================================
{
    vector<double> basisvals_u(2 * basis_u_.order());
    vector<double> basisvals_v(2 * basis_v_.order());

    // Compute basis values
    if (evaluate_from_right)
      {
	basis_u_.computeBasisValues(param[0], &basisvals_u[0], 1);
	basis_v_.computeBasisValues(param[1], &basisvals_v[0], 1);
      }
    else 
      {
	basis_u_.computeBasisValuesLeft(param[0], &basisvals_u[0], 1);
	basis_v_.computeBasisValuesLeft(param[1], &basisvals_v[0], 1);
      }

    computeBasis(basisvals_u.begin(), basisvals_v.begin(),
		 basis_u_.lastKnotInterval(),
		 basis_v_.lastKnotInterval(),
		 basisValues,
		 basisDerivs_u,
		 basisDerivs_v);
}

//===========================================================================
void SplineSurface::computeBasis(const vector<double>::const_iterator& bas_vals_u,
				 const vector<double>::const_iterator& bas_vals_v,
				 int left_u,
				 int left_v,
				 vector<double>& basisValues,
				 vector<double>& basisDerivs_u,
				 vector<double>& basisDerivs_v) const
//===========================================================================
{
    int kk1 = basis_u_.order();
    int kk2 = basis_v_.order();
    int nn1 = basis_u_.numCoefs();

    // Accumulate
    int ki, kj, kr;
    basisValues.resize(kk1*kk2);
    basisDerivs_u.resize(kk1*kk2);
    basisDerivs_v.resize(kk1*kk2);
    vector<double> weights(kk1*kk2);
    if (rational_)
    {
	int kdim = dim_ + 1;
	int uleft = left_u - kk1 + 1;
	int vleft = left_v - kk2 + 1;
	for (kj=vleft, kr=0; kj<vleft+kk2; ++kj)
	    for (ki=uleft; ki<uleft+kk1; ++ki)
		weights[kr++] = rcoefs_[(kj*nn1+ki)*kdim+dim_];
    }

    accumulateBasis(bas_vals_u, bas_vals_v,
		    weights, basisValues,
		    basisDerivs_u, basisDerivs_v);

}

//===========================================================================
void SplineSurface::computeBasisGrid(const Dvector& param_u,
				     const Dvector& param_v,
				     Dmatrix& basisValues) const 
//===========================================================================
{
  int numu = (int)param_u.size();
  int numv = (int)param_v.size();
  int uorder = basis_u_.order();
  int vorder = basis_v_.order();
  int ucoefs = basis_u_.numCoefs();
  int vcoefs = basis_v_.numCoefs();

  vector<double> basisvals_u(numu * uorder);
  vector<double> basisvals_v(numv * vorder);
  vector<int>    left_u(numu);
  vector<int>    left_v(numv);

  // Compute basis values
  basis_u_.computeBasisValues(&param_u[0], &param_u[0]+param_u.size(),
			      &basisvals_u[0], &left_u[0]);
  basis_v_.computeBasisValues(&param_v[0], &param_v[0]+param_v.size(),
			      &basisvals_v[0], &left_v[0]);

  // Initiate to zero
  int ki;
  int num_coefs = ucoefs*vcoefs;
  int num_par = numu*numv;
  basisValues.resize(num_par);
  for (ki=0; ki<num_par; ++ki)
    basisValues[ki].assign(num_coefs, 0.0);

  // Fetch all weights
  vector<double> weights, currw;
  if (rational_)
  {
      currw.resize(uorder*vorder);
      weights.resize(ucoefs*vcoefs);
      getWeights(weights);
  }

  // For all points
  int kj, kh, kv;
  int idx1, idx2, idx4, idx5;
  int numorder = uorder*vorder;
  vector<double> tmpVal(numorder);
  for (kh=0, kj=0, idx2=0; kj<numv; ++kj, idx2+=vorder)
  {
      for (ki=0, idx1=0; ki<numu; ++ki, ++kh, idx1+=uorder)
      {
	  int uleft = left_u[ki] - uorder + 1;
	  int vleft = left_v[kj] - vorder + 1;
	  if (rational_)
	  {
	      // Collect relevant weights
	      vector<double>::iterator wgt = weights.begin() + vleft*ucoefs;
	      vector<double>::iterator currwgt = currw.begin();
	      for (kv=0; kv<vorder; ++kv, wgt+=ucoefs, currwgt+=uorder)
	      {
		  std::copy(wgt+uleft, wgt+uleft+uorder, currwgt);
	      }
	  }
      
	  accumulateBasis(basisvals_u.begin() + idx1, basisvals_v.begin() + idx2,
			  currw, tmpVal);

	  // Copy results into output array
	  for (kv=0, idx4=0, idx5=vleft*ucoefs+uleft; kv<vorder; ++kv, idx4+=uorder, idx5+=ucoefs)
	  {
	      std::copy(tmpVal.begin()+idx4, tmpVal.begin()+idx4+uorder, 
			basisValues[kh].begin()+idx5);
	  }
      }
  }
  
}

//===========================================================================
void SplineSurface::computeBasisGrid(const Dvector& param_u,
				     const Dvector& param_v,
				     Dmatrix& basisValues,
				     Dmatrix& basisDerivs_u,
				     Dmatrix& basisDerivs_v,
				     bool evaluate_from_right) const 
//===========================================================================
{
  int derivs = 1;  // Compute position  and 1. derivative
  int numu = (int)param_u.size();
  int numv = (int)param_v.size();
  int uorder = basis_u_.order();
  int vorder = basis_v_.order();
  int ucoefs = basis_u_.numCoefs();
  int vcoefs = basis_v_.numCoefs();

  vector<double> basisvals_u(numu * uorder * (derivs + 1));
  vector<double> basisvals_v(numv * vorder * (derivs + 1));
  vector<int>    left_u(numu);
  vector<int>    left_v(numv);

  // Compute basis values
  if (evaluate_from_right)
    {
      basis_u_.computeBasisValues(&param_u[0], &param_u[0]+param_u.size(),
				  &basisvals_u[0], &left_u[0], derivs);
      basis_v_.computeBasisValues(&param_v[0], &param_v[0]+param_v.size(),
				  &basisvals_v[0], &left_v[0], derivs);
    }
  else 
    {
      basis_u_.computeBasisValuesLeft(&param_u[0], &param_u[0]+param_u.size(),
				      &basisvals_u[0], &left_u[0], derivs);
      basis_v_.computeBasisValuesLeft(&param_v[0], &param_v[0]+param_v.size(),
				      &basisvals_v[0], &left_v[0], derivs);
    }

  // Initiate to zero
  int ki;
  int num_coefs = ucoefs*vcoefs;
  int num_par = numu*numv;
  basisValues.resize(num_par);
  basisDerivs_u.resize(num_par);
  basisDerivs_v.resize(num_par);
  for (ki=0; ki<num_par; ++ki)
  {
      basisValues[ki].assign(num_coefs, 0.0);
      basisDerivs_u[ki].assign(num_coefs, 0.0);
      basisDerivs_v[ki].assign(num_coefs, 0.0);
  }

  // Fetch all weights
  vector<double> weights, currw;
  if (rational_)
  {
      currw.resize(uorder*vorder);
      weights.resize(ucoefs*vcoefs);
      getWeights(weights);
  }

  // For all points
  int kj, kh, kv;
  int idx1, idx2, idx4, idx5;
  int numorder = uorder*vorder;
  vector<double> tmpVal(numorder), tmpDer_u(numorder), tmpDer_v(numorder);
  for (kh=0, kj=0, idx2=0; kj<numv; ++kj, idx2+=2*vorder)
  {
      for (ki=0, idx1=0; ki<numu; ++ki, ++kh, idx1+=2*uorder)
      {
	  int uleft = left_u[ki] - uorder + 1;
	  int vleft = left_v[kj] - vorder + 1;
	  if (rational_)
	  {
	      // Collect relevant weights
	      vector<double>::iterator wgt = weights.begin() + vleft*ucoefs;
	      vector<double>::iterator currwgt = currw.begin();
	      for (kv=0; kv<vorder; ++kv, wgt+=ucoefs, currwgt+=uorder)
	      {
		  std::copy(wgt+uleft, wgt+uleft+uorder, currwgt);
	      }
	  }
      
	  accumulateBasis(basisvals_u.begin() + idx1, basisvals_v.begin() + idx2,
			  currw, tmpVal, tmpDer_u, tmpDer_v);

	  // Copy results into output array
	  for (kv=0, idx4=0, idx5=vleft*ucoefs+uleft; kv<vorder; ++kv, idx4+=uorder, idx5+=ucoefs)
	  {
	      std::copy(tmpVal.begin()+idx4, tmpVal.begin()+idx4+uorder, 
			basisValues[kh].begin()+idx5);
	      std::copy(tmpDer_u.begin()+idx4, tmpDer_u.begin()+idx4+uorder, 
			basisDerivs_u[kh].begin()+idx5);
	      std::copy(tmpDer_v.begin()+idx4, tmpDer_v.begin()+idx4+uorder, 
			basisDerivs_v[kh].begin()+idx5);
	  }
      }
  }
  
}

//===========================================================================
void SplineSurface::computeBasis(double param_u,
				 double param_v,
				 BasisPtsSf& result) const
//===========================================================================
{
    int uorder = basis_u_.order();
    int vorder = basis_v_.order();
    int nn1 = basis_u_.numCoefs();
    vector<double> basisvals_u(uorder);
    vector<double> basisvals_v(vorder);

    // Compute basis values
    basis_u_.computeBasisValues(param_u, &basisvals_u[0], 0);
    basis_v_.computeBasisValues(param_v, &basisvals_v[0], 0);

    int ulast = basis_u_.lastKnotInterval();
    int vlast = basis_v_.lastKnotInterval();
    result.preparePts(param_u, param_v, ulast, vlast,
		      uorder*vorder);

    vector<double> weights;
   if (rational_)
    {
      // Collect relevant weights
      int kr, ki, kj;
      int kdim = dim_ + 1;
      int uleft = ulast - uorder + 1;
      int vleft = vlast - vorder + 1;
      weights.resize(uorder*vorder);
      for (kj=vleft, kr=0; kj<vleft+vorder; ++kj)
	for (ki=uleft; ki<uleft+uorder; ++ki)
	  weights[kr++] = rcoefs_[(kj*nn1+ki)*kdim+dim_];
    }
      
   accumulateBasis(basisvals_u.begin(), basisvals_v.begin(),
		  weights, result.basisValues);
}

//===========================================================================
void SplineSurface::computeBasis(double param_u,
				 double param_v,
				 BasisDerivsSf& result,
				 bool evaluate_from_right) const
//===========================================================================
{
  int derivs = 1;  // Compute position  and 1. derivative
  int uorder = basis_u_.order();
  int vorder = basis_v_.order();
  int nn1 = basis_u_.numCoefs();
  vector<double> basisvals_u(uorder * (derivs + 1));
  vector<double> basisvals_v(vorder * (derivs + 1));

  // Compute basis values
  if (evaluate_from_right)
    {
      basis_u_.computeBasisValues(param_u, &basisvals_u[0], derivs);
      basis_v_.computeBasisValues(param_v, &basisvals_v[0], derivs);
    }
  else
    {
      basis_u_.computeBasisValuesLeft(param_u, &basisvals_u[0], derivs);
      basis_v_.computeBasisValuesLeft(param_v, &basisvals_v[0], derivs);
    }

  int ulast = basis_u_.lastKnotInterval();
  int vlast = basis_v_.lastKnotInterval();
  result.prepareDerivs(param_u, param_v, ulast, vlast,
		       uorder*vorder);

  vector<double> weights;
  if (rational_)
    {
      // Collect relevant weights
      int kr, ki, kj;
      int kdim = dim_ + 1;
      int uleft = ulast - uorder + 1;
      int vleft = vlast - vorder + 1;
      weights.resize(uorder*vorder);
      for (kj=vleft, kr=0; kj<vleft+vorder; ++kj)
	for (ki=uleft; ki<uleft+uorder; ++ki)
	  weights[kr++] = rcoefs_[(kj*nn1+ki)*kdim+dim_];
    }
 
  accumulateBasis(basisvals_u.begin(), basisvals_v.begin(),
		  weights, result.basisValues, 
		  result.basisDerivs_u, result.basisDerivs_v);
}

//===========================================================================
void SplineSurface::computeBasis(double param_u,
				 double param_v,
				 BasisDerivsSf2& result,
				 bool evaluate_from_right) const
//===========================================================================
{
  int derivs = 2;  // Compute position, 1. and 2. derivative
  int uorder = basis_u_.order();
  int vorder = basis_v_.order();
  int nn1 = basis_u_.numCoefs();
  vector<double> basisvals_u(uorder * (derivs + 1));
  vector<double> basisvals_v(vorder * (derivs + 1));

  // Compute basis values
  if (evaluate_from_right)
    {
      basis_u_.computeBasisValues(param_u, &basisvals_u[0], derivs);
      basis_v_.computeBasisValues(param_v, &basisvals_v[0], derivs);
    }
  else
    {
      basis_u_.computeBasisValuesLeft(param_u, &basisvals_u[0], derivs);
      basis_v_.computeBasisValuesLeft(param_v, &basisvals_v[0], derivs);
    }

  int ulast = basis_u_.lastKnotInterval();
  int vlast = basis_v_.lastKnotInterval();
  result.prepareDerivs(param_u, param_v, ulast, vlast,
		       uorder*vorder);

  vector<double> weights;
  if (rational_)
    {
      // Collect relevant weights
      int kr, ki, kj;
      int kdim = dim_ + 1;
      int uleft = ulast - uorder + 1;
      int vleft = vlast - vorder + 1;
      weights.resize(uorder*vorder);
      for (kj=vleft, kr=0; kj<vleft+vorder; ++kj)
	for (ki=uleft; ki<uleft+uorder; ++ki)
	  weights[kr++] = rcoefs_[(kj*nn1+ki)*kdim+dim_];
    }
 
  accumulateBasis(basisvals_u.begin(), basisvals_v.begin(),
		  weights, result.basisValues, 
		  result.basisDerivs_u, result.basisDerivs_v,
		  result.basisDerivs_uu, result.basisDerivs_uv, result.basisDerivs_vv);
}


//===========================================================================
void SplineSurface::computeBasis(double param_u,
                                 double param_v,
                                 BasisDerivsSf3& result,
                                 bool evaluate_from_right) const
//===========================================================================
{
  int derivs = 3;  // Compute position, 1. and 2. derivative
  int uorder = basis_u_.order();
  int vorder = basis_v_.order();
  int nn1 = basis_u_.numCoefs();
  vector<double> basisvals_u(uorder * (derivs + 1));
  vector<double> basisvals_v(vorder * (derivs + 1));

  // Compute basis values
  if (evaluate_from_right)
    {
      basis_u_.computeBasisValues(param_u, &basisvals_u[0], derivs);
      basis_v_.computeBasisValues(param_v, &basisvals_v[0], derivs);
    }
  else
    {
      basis_u_.computeBasisValuesLeft(param_u, &basisvals_u[0], derivs);
      basis_v_.computeBasisValuesLeft(param_v, &basisvals_v[0], derivs);
    }

  int ulast = basis_u_.lastKnotInterval();
  int vlast = basis_v_.lastKnotInterval();
  result.prepareDerivs(param_u, param_v, ulast, vlast,
		                   uorder*vorder);

  vector<double> weights;
  if (rational_)
  {
    // Collect relevant weights
    int kr, ki, kj;
    int kdim = dim_ + 1;
    int uleft = ulast - uorder + 1;
    int vleft = vlast - vorder + 1;
    weights.resize(uorder*vorder);
    for (kj=vleft, kr=0; kj<vleft+vorder; ++kj)
      for (ki=uleft; ki<uleft+uorder; ++ki)
        weights[kr++] = rcoefs_[(kj*nn1+ki)*kdim+dim_];
  }

  accumulateBasis(basisvals_u.begin(), basisvals_v.begin(), weights, result);
}


//===========================================================================
void SplineSurface::computeBasis(double param_u,
				 double param_v,
                                 int derivs,
				 BasisDerivsSfU& result,
				 bool evaluate_from_right) const
//===========================================================================
{
  int uorder = basis_u_.order();
  int vorder = basis_v_.order();
  int nn1 = basis_u_.numCoefs();
  vector<double> basisvals_u(uorder * (derivs + 1));
  vector<double> basisvals_v(vorder * (derivs + 1));

  // Compute basis values
  if (evaluate_from_right)
    {
      basis_u_.computeBasisValues(param_u, &basisvals_u[0], derivs);
      basis_v_.computeBasisValues(param_v, &basisvals_v[0], derivs);
    }
  else
    {
      basis_u_.computeBasisValuesLeft(param_u, &basisvals_u[0], derivs);
      basis_v_.computeBasisValuesLeft(param_v, &basisvals_v[0], derivs);
    }

  int ulast = basis_u_.lastKnotInterval();
  int vlast = basis_v_.lastKnotInterval();
  result.prepareDerivs(param_u, param_v, ulast, vlast,
		       derivs, uorder*vorder);

  vector<double> weights;
  if (rational_)
    {
      // Collect relevant weights
      int kr, ki, kj;
      int kdim = dim_ + 1;
      int uleft = ulast - uorder + 1;
      int vleft = vlast - vorder + 1;
      weights.resize(uorder*vorder);
      for (kj=vleft, kr=0; kj<vleft+vorder; ++kj)
	for (ki=uleft; ki<uleft+uorder; ++ki)
	  weights[kr++] = rcoefs_[(kj*nn1+ki)*kdim+dim_];
    }

  accumulateBasis(basisvals_u, basisvals_v, weights, result);
}

//===========================================================================
void SplineSurface::computeBasisGrid(const Dvector& param_u,
				     const Dvector& param_v,
				     vector<BasisPtsSf>& result) const
//===========================================================================
{
  int numu = (int)param_u.size();
  int numv = (int)param_v.size();
  int uorder = basis_u_.order();
  int vorder = basis_v_.order();
  int ucoefs = basis_u_.numCoefs();
  int vcoefs = basis_v_.numCoefs();

  vector<double> basisvals_u(numu * uorder);
  vector<double> basisvals_v(numv * vorder);
  vector<int>    left_u(numu);
  vector<int>    left_v(numv);

  // Compute basis values
  basis_u_.computeBasisValues(&param_u[0], &param_u[0]+param_u.size(),
			      &basisvals_u[0], &left_u[0]);
  basis_v_.computeBasisValues(&param_v[0], &param_v[0]+param_v.size(),
			      &basisvals_v[0], &left_v[0]);

  // Initiate output
  result.resize(numu*numv);

  // Fetch all weights
  vector<double> weights, currw;
  if (rational_)
  {
      currw.resize(uorder*vorder);
      weights.resize(ucoefs*vcoefs);
      getWeights(weights);
  }

  // For all points
  int ki, kj, kh, kv;
  int idx1, idx2;
  for (kh=0, kj=0, idx2=0; kj<numv; ++kj, idx2+=vorder)
  {
      for (ki=0, idx1=0; ki<numu; ++ki, ++kh, idx1+=uorder)
      {
	  result[kh].preparePts(param_u[ki], param_v[kj], 
				left_u[ki], left_v[kj], 
				uorder*vorder);


	  if (rational_)
	  {
	      // Collect relevant weights
	      int uleft = left_u[ki] - uorder + 1;
	      int vleft = left_v[kj] - vorder + 1;
	      vector<double>::iterator wgt = weights.begin() + vleft*ucoefs;
	      vector<double>::iterator currwgt = currw.begin();
	      for (kv=0; kv<vorder; ++kv, wgt+=ucoefs, currwgt+=uorder)
	      {
		  std::copy(wgt+uleft, wgt+uleft+uorder, currwgt);
	      }
	  }
      
	  accumulateBasis(basisvals_u.begin() + idx1, basisvals_v.begin() + idx2,
			  currw, result[kh].basisValues);
	  }
      }
}
  

//===========================================================================
void SplineSurface::computeBasisGrid(const Dvector& param_u,
				     const Dvector& param_v,
				     vector<BasisDerivsSf>& result,
				     bool evaluate_from_right) const
//===========================================================================
{
  int derivs = 1;  // Compute position  and 1. derivative
  int numu = (int)param_u.size();
  int numv = (int)param_v.size();
  int uorder = basis_u_.order();
  int vorder = basis_v_.order();
  int ucoefs = basis_u_.numCoefs();
  int vcoefs = basis_v_.numCoefs();

  vector<double> basisvals_u(numu * uorder * (derivs + 1));
  vector<double> basisvals_v(numv * vorder * (derivs + 1));
  vector<int>    left_u(numu);
  vector<int>    left_v(numv);

  // Compute basis values
  if (evaluate_from_right)
    {
      basis_u_.computeBasisValues(&param_u[0], &param_u[0]+param_u.size(),
				  &basisvals_u[0], &left_u[0], derivs);
      basis_v_.computeBasisValues(&param_v[0], &param_v[0]+param_v.size(),
				  &basisvals_v[0], &left_v[0], derivs);
    }
  else 
    {
      basis_u_.computeBasisValuesLeft(&param_u[0], &param_u[0]+param_u.size(),
				      &basisvals_u[0], &left_u[0], derivs);
      basis_v_.computeBasisValuesLeft(&param_v[0], &param_v[0]+param_v.size(),
				      &basisvals_v[0], &left_v[0], derivs);
    }

  // Initiate output
  result.resize(numu*numv);

  // Fetch all weights
  vector<double> weights, currw;
  if (rational_)
  {
      currw.resize(uorder*vorder);
      weights.resize(ucoefs*vcoefs);
      getWeights(weights);
  }

  // For all points
  int ki, kj, kh, kv;
  int idx1, idx2;
  for (kh=0, kj=0, idx2=0; kj<numv; ++kj, idx2+=2*vorder)
  {
      for (ki=0, idx1=0; ki<numu; ++ki, ++kh, idx1+=2*uorder)
      {
	  result[kh].prepareDerivs(param_u[ki], param_v[kj], 
				   left_u[ki], left_v[kj], 
				   uorder*vorder);


	  if (rational_)
	  {
	      // Collect relevant weights
	      int uleft = left_u[ki] - uorder + 1;
	      int vleft = left_v[kj] - vorder + 1;
	      vector<double>::iterator wgt = weights.begin() + vleft*ucoefs;
	      vector<double>::iterator currwgt = currw.begin();
	      for (kv=0; kv<vorder; ++kv, wgt+=ucoefs, currwgt+=uorder)
	      {
		  std::copy(wgt+uleft, wgt+uleft+uorder, currwgt);
	      }
	  }
      
	  accumulateBasis(basisvals_u.begin() + idx1, basisvals_v.begin() + idx2,
			  currw, result[kh].basisValues, 
			  result[kh].basisDerivs_u, result[kh].basisDerivs_v);
	  }
      }
}
  

//===========================================================================
void SplineSurface::computeBasisGrid(const Dvector& param_u,
				     const Dvector& param_v,
				     vector<BasisDerivsSf2>& result,
				     bool evaluate_from_right) const
//===========================================================================
{
  int derivs = 2;  // Compute position, 1. and 2. derivative
  int numu = (int)param_u.size();
  int numv = (int)param_v.size();
  int uorder = basis_u_.order();
  int vorder = basis_v_.order();
  int ucoefs = basis_u_.numCoefs();
  int vcoefs = basis_v_.numCoefs();

  vector<double> basisvals_u(numu * uorder * (derivs + 1));
  vector<double> basisvals_v(numv * vorder * (derivs + 1));
  vector<int>    left_u(numu);
  vector<int>    left_v(numv);

  // Compute basis values
  if (evaluate_from_right)
    {
      basis_u_.computeBasisValues(&param_u[0], &param_u[0]+param_u.size(),
				  &basisvals_u[0], &left_u[0], derivs);
      basis_v_.computeBasisValues(&param_v[0], &param_v[0]+param_v.size(),
				  &basisvals_v[0], &left_v[0], derivs);
    }
  else 
    {
      basis_u_.computeBasisValuesLeft(&param_u[0], &param_u[0]+param_u.size(),
				      &basisvals_u[0], &left_u[0], derivs);
      basis_v_.computeBasisValuesLeft(&param_v[0], &param_v[0]+param_v.size(),
				      &basisvals_v[0], &left_v[0], derivs);
    }

  // Initiate output
  result.resize(numu*numv);

  // Fetch all weights
  vector<double> weights, currw;
  if (rational_)
  {
      currw.resize(uorder*vorder);
      weights.resize(ucoefs*vcoefs);
      getWeights(weights);
  }

  // For all points
  int ki, kj, kh, kv;
  int idx1, idx2;
  for (kh=0, kj=0, idx2=0; kj<numv; ++kj, idx2+=3*vorder)
  {
      for (ki=0, idx1=0; ki<numu; ++ki, ++kh, idx1+=3*uorder)
      {
	  result[kh].prepareDerivs(param_u[ki], param_v[kj], 
				   left_u[ki], left_v[kj], 
				   uorder*vorder);

	  if (rational_)
	  {
	      // Collect relevant weights
	      int uleft = left_u[ki] - uorder + 1;
	      int vleft = left_v[kj] - vorder + 1;
	      vector<double>::iterator wgt = weights.begin() + vleft*ucoefs;
	      vector<double>::iterator currwgt = currw.begin();
	      for (kv=0; kv<vorder; ++kv, wgt+=ucoefs, currwgt+=uorder)
	      {
		  std::copy(wgt+uleft, wgt+uleft+uorder, currwgt);
	      }
	  }
      
	  accumulateBasis(basisvals_u.begin() + idx1, basisvals_v.begin() + idx2,
			  currw, result[kh].basisValues, 
			  result[kh].basisDerivs_u, result[kh].basisDerivs_v,
			  result[kh].basisDerivs_uu, result[kh].basisDerivs_uv, result[kh].basisDerivs_vv);
      }
  }
  
}


//===========================================================================
void SplineSurface::computeBasisGrid(const Dvector& param_u,
                                     const Dvector& param_v,
                                     vector<BasisDerivsSf3>& result,
                                     bool evaluate_from_right) const
//===========================================================================
{
  int derivs = 3;  // Compute position, 1., 2. and 3. derivative
  int numu = (int)param_u.size();
  int numv = (int)param_v.size();
  int uorder = basis_u_.order();
  int vorder = basis_v_.order();
  int ucoefs = basis_u_.numCoefs();
  int vcoefs = basis_v_.numCoefs();

  vector<double> basisvals_u(numu * uorder * (derivs + 1));
  vector<double> basisvals_v(numv * vorder * (derivs + 1));
  vector<int>    left_u(numu);
  vector<int>    left_v(numv);

  // Compute basis values
  if (evaluate_from_right)
  {
    basis_u_.computeBasisValues(&param_u[0], &param_u[0]+param_u.size(),
                                &basisvals_u[0], &left_u[0], derivs);
    basis_v_.computeBasisValues(&param_v[0], &param_v[0]+param_v.size(),
                                &basisvals_v[0], &left_v[0], derivs);
  }
  else
  {
    basis_u_.computeBasisValuesLeft(&param_u[0], &param_u[0]+param_u.size(),
                                    &basisvals_u[0], &left_u[0], derivs);
    basis_v_.computeBasisValuesLeft(&param_v[0], &param_v[0]+param_v.size(),
                                    &basisvals_v[0], &left_v[0], derivs);
  }

  // Initiate output
  result.resize(numu*numv);

  // Fetch all weights
  vector<double> weights, currw;
  if (rational_)
  {
      currw.resize(uorder*vorder);
      weights.resize(ucoefs*vcoefs);
      getWeights(weights);
  }

  // For all points
  int ki, kj, kh, kv;
  int idx1, idx2;
  for (kh=0, kj=0, idx2=0; kj<numv; ++kj, idx2+=4*vorder)
  {
    for (ki=0, idx1=0; ki<numu; ++ki, ++kh, idx1+=4*uorder)
    {
      result[kh].prepareDerivs(param_u[ki], param_v[kj],
                               left_u[ki], left_v[kj],
                               uorder*vorder);

      if (rational_)
      {
        // Collect relevant weights
        int uleft = left_u[ki] - uorder + 1;
        int vleft = left_v[kj] - vorder + 1;
        vector<double>::iterator wgt = weights.begin() + vleft*ucoefs;
        vector<double>::iterator currwgt = currw.begin();
        for (kv=0; kv<vorder; ++kv, wgt+=ucoefs, currwgt+=uorder)
        {
          std::copy(wgt+uleft, wgt+uleft+uorder, currwgt);
        }
      }

      accumulateBasis(basisvals_u.begin() + idx1, basisvals_v.begin() + idx2,
                      currw, result[kh]);
    }
  }
}



//===========================================================================
void SplineSurface::accumulateBasis(const vector<double>::const_iterator& basisvals_u,
				    const vector<double>::const_iterator& basisvals_v,
				    const vector<double>& weights,
				    vector<double>& basisValues) const
//===========================================================================
{
    int ki, kj, kr;
    int uorder = basis_u_.order();
    int vorder = basis_v_.order();
    if (rational_)
    {
	// Compute denominator and derivatives thereof
        double sum = 0.0;
	for (kj=0, kr=0; kj<vorder; ++kj)
	    for (ki=0; ki<uorder; ++ki, ++kr)
	        sum += basisvals_u[ki]*basisvals_v[kj]*weights[kr];

	// Compute rational expression
	// double sum2 = sum*sum;
	for (kj=0, kr=0; kj<vorder; ++kj)
	    for (ki=0; ki<uorder; ++ki, ++kr)
		basisValues[kr] = basisvals_u[ki]*basisvals_v[kj]*weights[kr]/sum;
    }
    else
    {
	// Multiply basis values in the three parameter directions
	for (kj=0, kr=0; kj<vorder; ++kj)
	    for (ki=0; ki<uorder; ++ki, ++kr)
		basisValues[kr] = basisvals_u[ki]*basisvals_v[kj];
    }

}

//===========================================================================
void SplineSurface::accumulateBasis(const vector<double>::const_iterator& basisvals_u,
				    const vector<double>::const_iterator& basisvals_v,
				    const vector<double>& weights,
				    vector<double>& basisValues,
				    vector<double>& basisDerivs_u,
				    vector<double>& basisDerivs_v) const
//===========================================================================
{
    int ki, kj, kr;
    int uorder = basis_u_.order();
    int vorder = basis_v_.order();
    if (rational_)
    {
	// Compute denominator and derivatives thereof
	double sum = 0.0, dusum = 0.0, dvsum = 0.0;
	for (kj=0, kr=0; kj<vorder; ++kj)
	    for (ki=0; ki<uorder; ++ki, ++kr)
	    {
		sum += basisvals_u[2*ki]*basisvals_v[2*kj]*weights[kr];
		dusum += basisvals_u[2*ki+1]*basisvals_v[2*kj]*weights[kr];
		dvsum += basisvals_u[2*ki]*basisvals_v[2*kj+1]*weights[kr];
	    }

	// Compute rational expression
	double sum2 = sum*sum;
	for (kj=0, kr=0; kj<vorder; ++kj)
	{
	    double tmp2 = (basisvals_v[2*kj+1]*sum - basisvals_v[2*kj]*dvsum)/sum2;
	    double fac = basisvals_v[2*kj]/sum2;
	    for (ki=0; ki<uorder; ++ki, ++kr)
	    {
		basisValues[kr] = basisvals_u[2*ki]*basisvals_v[2*kj]*weights[kr]/sum;
		double tmp1 = (basisvals_u[2*ki+1]*sum - basisvals_u[2*ki]*dusum)*weights[kr]*fac;
		    basisDerivs_u[kr] = tmp1;
		    basisDerivs_v[kr] = tmp2*weights[kr]*basisvals_u[2*ki];
	    }
	}
    }
    else
    {
	// Multiply basis values in the two parameter directions
	for (kj=0, kr=0; kj<vorder; ++kj)
	    for (ki=0; ki<uorder; ++ki, ++kr)
	    {
		basisValues[kr] = basisvals_u[2*ki]*basisvals_v[2*kj];
		basisDerivs_u[kr] = basisvals_u[2*ki+1]*basisvals_v[2*kj];
		basisDerivs_v[kr] = basisvals_u[2*ki]*basisvals_v[2*kj+1];
	    }
    }

}


//===========================================================================
void SplineSurface::accumulateBasis(const vector<double>::const_iterator& basisvals_u,
				    const vector<double>::const_iterator& basisvals_v,
				    const vector<double>& weights,
				    vector<double>& basisValues,
				    vector<double>& basisDerivs_u,
				    vector<double>& basisDerivs_v,
				    vector<double>& basisDerivs_uu,
				    vector<double>& basisDerivs_uv,
				    vector<double>& basisDerivs_vv) const
//===========================================================================
{
    int ki, kj, kr;
    int uorder = basis_u_.order();
    int vorder = basis_v_.order();
    if (rational_)
    {
	// Compute denominator and derivatives thereof
	double
	  sum = 0.0, dusum = 0.0, dvsum = 0.0,
	  duusum = 0.0, duvsum = 0.0, dvvsum = 0.0;
	for (kj=0, kr=0; kj<vorder; ++kj)
	    for (ki=0; ki<uorder; ++ki, ++kr)
	    {
	        sum += basisvals_u[3*ki]*basisvals_v[3*kj]*weights[kr];
		dusum += basisvals_u[3*ki+1]*basisvals_v[3*kj]*weights[kr];
		dvsum += basisvals_u[3*ki]*basisvals_v[3*kj+1]*weights[kr];
		duusum += basisvals_u[3*ki+2]*basisvals_v[3*kj]*weights[kr];
		duvsum += basisvals_u[3*ki+1]*basisvals_v[3*kj+1]*weights[kr];
		dvvsum += basisvals_u[3*ki]*basisvals_v[3*kj+2]*weights[kr];
	    }

	// Compute rational expression
	double sum2 = sum*sum;
	double sum3 = sum*sum2;
	for (kj=0, kr=0; kj<vorder; ++kj)
	{
	    double tmp_v_v = (basisvals_v[3*kj+1]*sum - basisvals_v[3*kj]*dvsum)/sum2;
	    double tmp_v_vv = (basisvals_v[3*kj+2]*sum2 + 2.0*basisvals_v[3*kj]*dvsum*dvsum
			       - basisvals_v[3*kj]*dvvsum*sum - 2.0*basisvals_v[3*kj+1]*dvsum*sum) / sum3;
	    double fac_v2 = basisvals_v[3*kj] / sum2;
	    double fac_v3 = basisvals_v[3*kj] / sum3;
	    for (ki=0; ki<uorder; ++ki, ++kr)
	    {

		basisValues[kr] = basisvals_u[3*ki]*basisvals_v[3*kj]*weights[kr]/sum;

		basisDerivs_u[kr] = (basisvals_u[3*ki+1]*sum - basisvals_u[3*ki]*dusum) * weights[kr] * fac_v2;
		basisDerivs_v[kr] = tmp_v_v * weights[kr] * basisvals_u[3*ki];

		basisDerivs_uu[kr] = (basisvals_u[3*ki+2]*sum2 + 2.0*basisvals_u[3*ki]*dusum*dusum
				      - basisvals_u[3*ki]*duusum*sum - 2.0*basisvals_u[3*ki+1]*dusum*sum) * weights[kr] * fac_v3;
		basisDerivs_uv[kr] = (basisvals_u[3*ki+1]*basisvals_v[3*kj+1]*sum2 + 2.0*basisvals_u[3*ki]*basisvals_v[3*kj]*dusum*dvsum
				      - basisvals_u[3*ki+1]*basisvals_v[3*kj]*dvsum*sum - basisvals_u[3*ki]*basisvals_v[3*kj+1]*dusum*sum
				      - basisvals_u[3*ki]*basisvals_v[3*kj]*duvsum*sum) * weights[kr] / sum3;
		basisDerivs_vv[kr] = tmp_v_vv * weights[kr] * basisvals_u[3*ki];
	    }
	}
    }
    else
    {
	// Multiply basis values in the two parameter directions
	for (kj=0, kr=0; kj<vorder; ++kj)
	    for (ki=0; ki<uorder; ++ki, ++kr)
	    {
		basisValues[kr] = basisvals_u[3*ki]*basisvals_v[3*kj];
		basisDerivs_u[kr] = basisvals_u[3*ki+1]*basisvals_v[3*kj];
		basisDerivs_v[kr] = basisvals_u[3*ki]*basisvals_v[3*kj+1];
		basisDerivs_uu[kr] = basisvals_u[3*ki+2]*basisvals_v[3*kj];
		basisDerivs_uv[kr] = basisvals_u[3*ki+1]*basisvals_v[3*kj+1];
		basisDerivs_vv[kr] = basisvals_u[3*ki]*basisvals_v[3*kj+2];
	    }
    }

}


//===========================================================================
void SplineSurface::accumulateBasis(const vector<double>::const_iterator& basisvals_u,
                                    const vector<double>::const_iterator& basisvals_v,
                                    const vector<double>& weights,
                                    BasisDerivsSf3& result) const
//===========================================================================
{
    size_t uorder = basis_u_.order();
    size_t vorder = basis_v_.order();
    if (rational_)
    {
      const int val    = 0;
      const int derx   = 1;
      const int derxx  = 2;
      const int derxxx = 3;
      const int derxy  = 4;
      const int derxxy = 5;
      const int derxyy = 6;
      const int dery   = 7;
      const int deryy  = 8;
      const int deryyy = 9;

      std::array<double, 10> W{};

      size_t kr = 0;
      for (size_t kj=0; kj<vorder; ++kj)
        for (size_t ki=0; ki<uorder; ++ki, ++kr) {
          W[val]    += basisvals_u[4*ki]   * basisvals_v[4*kj]   * weights[kr];
          W[derx]   += basisvals_u[4*ki+1] * basisvals_v[4*kj]   * weights[kr];
          W[derxx]  += basisvals_u[4*ki+2] * basisvals_v[4*kj]   * weights[kr];
          W[derxxx] += basisvals_u[4*ki+3] * basisvals_v[4*kj]   * weights[kr];
          W[derxy]  += basisvals_u[4*ki+1] * basisvals_v[4*kj+1] * weights[kr];
          W[derxxy] += basisvals_u[4*ki+2] * basisvals_v[4*kj+1] * weights[kr];
          W[derxyy] += basisvals_u[4*ki+1] * basisvals_v[4*kj+2] * weights[kr];
          W[dery]   += basisvals_u[4*ki]   * basisvals_v[4*kj+1] * weights[kr];
          W[deryy]  += basisvals_u[4*ki]   * basisvals_v[4*kj+2] * weights[kr];
          W[deryyy] += basisvals_u[4*ki]   * basisvals_v[4*kj+3] * weights[kr];
        }

      kr = 0;
      for (size_t kj=0; kj<vorder; ++kj)
        for (size_t ki=0; ki<uorder; ++ki, ++kr) {
          const double& Ni = basisvals_u[4*ki];
          const double& Nid = basisvals_u[4*ki+1];
          const double& Nidd = basisvals_u[4*ki+2];
          const double& Niddd = basisvals_u[4*ki+3];
          const double& Nj = basisvals_v[4*kj];
          const double& Njd = basisvals_v[4*kj+1];
          const double& Njdd = basisvals_v[4*kj+2];
          const double& Njddd = basisvals_v[4*kj+3];

          result.basisValues[kr] = Ni*Nj*weights[kr]/W[val];

          double nom2 = weights[kr]/pow(W[val],2.0);
          double H1 = (Nid*Nj*W[val] - Ni*Nj*W[derx]);
          double H2 = (Ni*Njd*W[val] - Ni*Nj*W[dery]);
          result.basisDerivs_u[kr] = H1 * nom2;
          result.basisDerivs_v[kr] = H2 * nom2;

          double dH1dx = Nidd*Nj*W[val] - Ni*Nj*W[derxx];
          double dH1dy = Nid*Njd*W[val] + Nid*Nj*W[dery] - Ni*Njd*W[derx] - Ni*Nj*W[derxy];
          double dH2dx = Nid*Njd*W[val] + Ni*Njd*W[derx] - Nid*Nj*W[dery] - Ni*Nj*W[derxy];
          double dH2dy = Ni*Njdd*W[val] - Ni*Nj*W[deryy];
          double nom3 = weights[kr]/pow(W[val],3.0);
          double G1 = (dH1dx*W[val] - 2*H1*W[derx]);
          double G2 = (dH2dy*W[val] - 2*H2*W[dery]);
          result.basisDerivs_uu[kr] = G1 * nom3;
          result.basisDerivs_vv[kr] = G2 * nom3;
          result.basisDerivs_uv[kr] = (dH1dy*W[val] - 2*H1*W[dery]) * nom3;

          double d2H1dx2  = Niddd*Nj*W[val] + Nidd*Nj*W[derx] - Nid*Nj*W[derxx] - Ni*Nj*W[derxxx];
          double d2H1dxdy = Nidd*Njd*W[val] + Nidd*Nj*W[dery] - Ni*Njd*W[derxx] - Ni*Nj*W[derxxy];
          double d2H2dy2  = Ni*Njddd*W[val] + Ni*Njdd*W[dery] - Ni*Njd*W[deryy] - Ni*Nj*W[deryyy];
          double d2H2dxdy = Nid*Njdd*W[val] + Ni*Njdd*W[derx] - Nid*Nj*W[deryy] - Ni*Nj*W[derxyy];

          double dG1dx = d2H1dx2*W[val] + dH1dx*W[derx] - 2*dH1dx*W[derx] - 2*H1*W[derxx];
          double dG1dy = d2H1dxdy*W[val] + dH1dx*W[dery] - 2*dH1dy*W[derx] - 2*H1*W[derxy];
          double dG2dx = d2H2dxdy*W[val] + dH2dy*W[derx] - 2*dH2dx*W[dery] - 2*H2*W[derxy];
          double dG2dy = d2H2dy2*W[val] + dH2dy*W[dery] - 2*dH2dy*W[dery] - 2*H2*W[deryy];

          double nom4 = weights[kr] / pow(W[val],4);
          result.basisDerivs_uuu[kr] = (dG1dx*W[val] - 3*G1*W[derx]) * nom4;
          result.basisDerivs_vvv[kr] = (dG2dy*W[val] - 3*G2*W[dery]) * nom4;
          result.basisDerivs_uuv[kr] = (dG1dy*W[val] - 3*G1*W[dery]) * nom4;
          result.basisDerivs_uvv[kr] = (dG2dx*W[val] - 3*G2*W[derx]) * nom4;
        }
    }
    else
    {
      // Multiply basis values in the two parameter directions
      size_t kr = 0;
      for (size_t kj=0; kj<vorder; ++kj)
        for (size_t ki=0; ki<uorder; ++ki, ++kr)
        {
          result.basisValues[kr] = basisvals_u[4*ki]*basisvals_v[4*kj];
          result.basisDerivs_u[kr] = basisvals_u[4*ki+1]*basisvals_v[4*kj];
          result.basisDerivs_v[kr] = basisvals_u[4*ki]*basisvals_v[4*kj+1];
          result.basisDerivs_uu[kr] = basisvals_u[4*ki+2]*basisvals_v[4*kj];
          result.basisDerivs_uv[kr] = basisvals_u[4*ki+1]*basisvals_v[4*kj+1];
          result.basisDerivs_vv[kr] = basisvals_u[4*ki]*basisvals_v[4*kj+2];
          result.basisDerivs_uuu[kr] = basisvals_u[4*ki+3]*basisvals_v[4*kj];
          result.basisDerivs_uuv[kr] = basisvals_u[4*ki+2]*basisvals_v[4*kj+1];
          result.basisDerivs_uvv[kr] = basisvals_u[4*ki+1]*basisvals_v[4*kj+2];
          result.basisDerivs_vvv[kr] = basisvals_u[4*ki]*basisvals_v[4*kj+3];
        }
    }

}


void SplineSurface::accumulateBasis(const vector<double>& basisvals_u,
				    const vector<double>& basisvals_v,
				    const vector<double>& weights,
                                    BasisDerivsSfU& output) const
{
  int ki, kj, kr;
  int uorder = basis_u_.order();
  int vorder = basis_v_.order();
  if (rational_)
  {
    std::cerr << "Not implemented for rationals!" << std::endl;
    ASSERT(0);
  }
  else
  {
    // Multiply basis values in the two parameter directions
    int nperfunc = output.derivs+1;
    for (kj=0, kr=0; kj<vorder; ++kj)
      for (ki=0; ki<uorder; ++ki, ++kr)
      {
        output.values[kr][0] = basisvals_u[ki*nperfunc]*basisvals_v[kj*nperfunc];
        for (size_t i=0;i<output.derivs;++i) {
          output.values[kr][2*i+1] = basisvals_u[ki*nperfunc+i+1]*basisvals_v[kj*nperfunc];
          output.values[kr][2*i+2] = basisvals_u[ki*nperfunc]*basisvals_v[kj*nperfunc+i+1];
        }
      }
  }
}


} // namespace Go



