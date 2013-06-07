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
#include "GoTools/geometry/SplineCurve.h"
#include <memory>
#include <algorithm>

using std::vector;

namespace Go {

// Anonymous namespace
namespace {


//==========================================================================
inline void s1927(const vector<double>& w1, int nur,
		  int ik, const vector<int>& ed,
		  const vector<double>& w2, int nrc,
		  const vector<double>& w3, int nlr,
		  vector<double>& ex, const vector<double>& ey)
//==========================================================================
{
    int nn = nur + nlr;
    int nlc = nn - nrc;

    // Allocate output array ex
    ex.resize(nn);
    fill(ex.begin(), ex.end(), 0.0);

    // Solve L*z = ey
    int ii, jj;
    for (ii = 0; ii < nur; ii++) {
	int di = ed[ii];
	double wii = w1[(di - 1) * nur + ii];
	double sum = ey[ii];
	if (di > 1) {
	    int dim = di - 1;
	    int midi = ii - di + 1;
	    for (jj = 0; jj < dim; jj++)
		sum -= w1[jj * nur + ii] * (ex[jj + midi]);
	}
	ex[ii] = sum / wii;
    }

    // Solve filled part of L*z = ey
    for (ii = nur; ii < nn; ii++) {
	int mur = ii - nur;
	double wii = w3[ii * nlr + mur];
	double sum = ey[ii];
	if (ii >= 1) {
	    for (jj = 0; jj < ii; jj++)
		sum -= w3[jj * nlr + mur] * ex[jj];
	}
	ex[ii] = sum / wii;
    }

    // Solve U*ex = z   ; Jump if filled part of U is exhausted
    for (ii = nn - 2; ii >= nur; ii--) {
	double sum = ex[ii];
	int mur = ii - nur;
	for (jj = ii + 1; jj < nn; jj++)
	    sum -= w3[jj * nlr + mur] * ex[jj];
	ex[ii] = sum;
    }

    // Test if w2 contains diagonal elements
    if (nlc < nn) {
	for (ii = nur - 1; ii >= 0; ii--) {
	    double sum = ex[ii];
	    for (jj = nlc; jj < nn; jj++)
		sum -= w2[(jj - nlc) * nur + ii] * ex[jj];
	    ex[ii] = sum;
	}
    }
    for (ii = nur - 1; ii >= 0; ii--) {
	int di = ed[ii];
	if (di < ik) {
	    double sum = ex[ii];
	    int midi = ii - di + 1;
	    for (int jj = di; jj < ik; jj++)
		sum -= w1[jj * nur + ii] * (ex[jj + midi]);
	    ex[ii] = sum;
	}
    }

    return;
}


//==========================================================================
inline void s1926(vector<double>& w1, int nur,
		  int ik, const vector<int>& ed,
		  vector<double>& w2, int nrc,
		  vector<double>& w3, int nlr)
//==========================================================================
{
    int nn = nur + nlr;
    int nlc = nn - nrc;

    // Elimination scheme, jump if band part of W is completed
    int ii, jj, ll;
    for (ii = 0; ii < nur; ii++) {
	int di = ed[ii];
	double wii = w1[(di - 1) * nur + ii];
	// Jump if W(ii,jj) is trivially zero, jj = ii+1,ii+2,...,nlr
	if (di < ik) {
	    for (jj = di; jj < ik; jj++)
		w1[jj * nur + ii] /= wii;
	    // Perform elimination row by row
	    int midi = ii - di;
	    for (ll = ii + 1;; ll++) {
		// Jump if ii-th element of rows of band-part has been
		// eliminated
		if (ll >= nur)
		    break;
		int midl = ll - ed[ll];
		// Jump if W(ii,jj) is trivially zero, jj = ll,ll+1,...,nur
		if (midl >= ii)
		    break;
		int korr = midl - midi;
		double wli = w1[(di - korr - 1) * nur + ll];
		for (jj = di; jj < ik; jj++)
		    w1[(jj - korr) * nur + ll] += -w1[jj * nur + ii] * wli;
	    }
	    //  Eliminate ii-th column of w3 using ii-th row from w1
	    if (nlr > 0) {
		for (ll = 0; ll < nlr; ll++) {
		    double wli = w3[ii * nlr + ll];
		    for (jj = di; jj < ik; jj++)
			w3[(jj + midi + 1) * nlr + ll]
			    -= w1[jj * nur + ii] * wli;
		}
	    }
	}
    }

    // Apply the above elimination scheme on w2
    if (nrc > 0) {
	// Jump if band part of W is completed or if system error
	// occures, i.e. if w2 contains some diagonal elements of W)
	for (ii = 0; ii < nur; ii++) {
	    int di = ed[ii];
	    double wii = w1[(di - 1) * nur + ii];
	    for (jj = 0; jj < nrc; jj++)
		w2[jj * nur + ii] /= wii;
	    // Perform elimination row by row
	    int midi = ii - di;
	    for (ll = ii + 1;; ll++) {
		//  Jump if ii-th element of rows of band-part has
		//  been eliminated
		if (ll >= nur)
		    break;
		int midl = ll - ed[ll];
		// Jump if W(ii,jj) is trivially zero, jj = ll,ll+1,...,nur
		if (midl >= ii)
		    break;
		int korr = midl - midi;
		double wli = w1[(di - korr - 1) * nur + ll];
		for (jj = 0; jj < nrc; jj++)
		    w2[jj * nur + ll] -= w2[jj * nur + ii] * wli;
	    }
	    //  Eliminate ii-th column of w3 using ii-th row from w2
	    for (ll = 0; ll < nlr; ll++) {
		double wli = w3[ii * nlr + ll];
		for (jj = nlc; jj < nn; jj++)
		    w3[jj * nlr + ll] -= w2[(jj - nlc) * nur + ii] * wli;
	    }
	}
    }

    // Eliminate w3-part of W
    for (ii = nur; ii < nn; ii++) {
	// 1 <= ii <= nn
	int mur = ii - nur;
	double wii = w3[ii * nlr + mur];
	//  1 <= ii < nn
	for (jj = ii + 1; jj < nn; jj++)
	    w3[jj * nlr + mur] /= wii;
	for (ll = mur + 1; ll < nlr; ll++) {
	    double wli = w3[ii * nlr + ll];
	    for (jj = ii + 1; jj < nn; jj++)
		w3[jj * nlr + ll] -= w3[jj * nlr + mur] * wli;
	}
    }

    return;
}


//==========================================================================
inline void s1897(const vector<double>& et, int ik, double ax, int left,
		  int deriv, vector<double>& ebiatx)
//==========================================================================
{
    vector<double> edltr(ik);
    vector<double> edltl(ik);

    ebiatx[0] = 1.0;

    int j;
    for (j = 1; j <= deriv; j++) {
	edltr[j - 1] = et[left + j] - ax;
	edltl[j - 1] = ax - et[left + 1 - j];
	double fak = double(j);

	double saved = 0.0;
	for (int count = 1; count <= j; count++) {
	    double dummy = edltr[count - 1] + edltl[j - count];

	    double term = fak * ebiatx[count - 1] / dummy;
	    ebiatx[count - 1] = saved - term;
	    saved = term;
	}
	ebiatx[j] = saved;
    }

    for (j = deriv + 1; j < ik; j++) {
	edltr[j - 1] = et[left + j] - ax;
	edltl[j - 1] = ax - et[left + 1 - j];
	double fak = j / double(j - deriv);

	double saved = 0.0;
	for (int count = 1; count <= j; count++) {
	    double dummy = edltr[count - 1] + edltl[j - count];

	    double term = fak * ebiatx[count - 1] / dummy;
	    ebiatx[count - 1] = saved + edltr[count - 1] * term;
	    saved = edltl[j - count] * term;
	}
	ebiatx[j] = saved;
    }

    return;
}


//==========================================================================
inline void s1925(const vector<double>& etau, const vector<double>& epoint,
		  int inbpnt,
		  const vector<int>& eder, const vector<double>& et,
		  vector<double>& ebcoef, int in, int ik,
		  int iright, int dim,
		  vector<double>& ew1, int nur, vector<int>& ed,
		  vector<double>& ew2, int inrc, vector<double>& ew3,
		  int inlr)
//==========================================================================
{
    int nn = nur + inlr;
    int nlc = inbpnt - inrc;
    int kk = (ik <= nlc ? ik : nlc);
    double taudel = 0.0;
    int leftdel = 0;
    int leftmin = ik - 1;
    int leftmax = in -1;

    // Band part of W

    int left = leftmin;
    int esize = std::max(ik, nn);
    vector<double> ebder(esize);
  
    int ii;
    for (ii = 0; ii < nur; ii++) {
	double taui = etau[ii];
	int ideri = eder[ii];

	// Locate left so that  et[left] <= taui < et[left+1]
	while (left < leftmax && et[left + 1] <= taui)
	    left++;

	// et(left-1) <= taui < et(left)
	ed[ii] = ii - (left - ik);

	int iadd = (0 >= ik - left - 1 ? 0 : ik - left - 1);
	int iadd_save = iadd;
	int kmiadd = ik - iadd;

	// Compute the value and ideri first derivative of the ik
	// (possibly) nonzero B-spline associated with the knot vector
	// et at a point

	if (iadd > 0) {
	    s1897(et, ik, taui + taudel, left + leftdel,
		  ideri, ebder);

	    ed[ii] -= iadd;
	    int ish = inrc - iadd;
	    for (int jj = 0; jj < iadd; jj++) {
		ew1[(jj + kmiadd) * nur + ii] = 0.0;
		ew2[(jj + ish) * nur + ii] = ebder[jj];
	    }
	} else {
	    s1897(et, ik, taui, left, ideri, ebder);
	}

	int isub = (0 >= ii - ed[ii] + kmiadd - nlc + 1
		    ? 0 : ii - ed[ii] + kmiadd - nlc + 1);
	int kmisub = ik - isub;
	if (isub > 0) {
	    int ish = isub - (kmiadd - kk);
	    ed[ii] += ish;
	    iadd -= ish;
	    int stop = (ik <= kmisub + inrc ? ik : kmisub + inrc);

	    for (int jj = kmisub; jj < stop; jj++) {
		ew1[(jj - kmisub) * nur + ii] = 0.0;
		ew2[(jj - kmisub) * nur + ii] = ebder[jj];
	    }
	}
	for (int jj = iadd_save; jj < kmisub; jj++)
	    ew1[(jj - iadd) * nur + ii] = ebder[jj];

    }

    // Band part of W is now completed

    if (ii < inbpnt) {
	// Will compute lower, filled part of W for closed curve
	// interpolation

	int store = inlr * inbpnt;
	for (int jj = 0; jj < store; jj++)
	    ew3[jj] = 0.0;

	// Repeat until filled part of W is completed
	for (; ii < inbpnt; ii++) {
	    double taui = etau[ii];

	    // Locate left so that  et[left] <= taui < et[left+1]
	    while (left < in -1 && et[left + 1] <= taui)
		left++;

	    // et(left-1) <= taui < et(left)
	    int ideri = eder[ii];

	    // Compute the value and the ideri first derivatives of
	    // the ik (possibly) nonzero B-spline associated with the
	    // knot vector et at the point (taui)
	    s1897(et, ik, taui, left, ideri, ebder);

	    int imnur = ii - nur;
	    int lfmkm = left - ik;
	    for (int jj = 0; jj < ik; jj++) {
		int isum = jj + lfmkm + 1;
		int kmod;
		if (isum >= 0)
		    kmod = isum % inbpnt;
		if (isum < 0)
		    kmod = (isum + 1) % inbpnt + in -1;
		ew3[kmod * inlr + imnur] = ebder[jj];
	    }
	}
    }

    // W is now contained in ew1, ew2 and ew3 as required by the
    // subroutine s1898

    s1926(ew1, nur, kk, ed, ew2, inrc, ew3, inlr);

    int store = iright * dim * inbpnt;
    int jj;
    for (jj = 0; jj < store; jj++)
	ebcoef[jj] = epoint[jj];

    // epoint is now properly contained in ebcoef.
    // Solve interpolation equations
    int kl, dim1;
    for (kl = 0; kl < iright; kl++) {
	for (dim1 = 0; dim1 < dim; dim1++) {
	    int store = inbpnt * dim * kl + dim1;
	    for (jj = 0; jj < nn; jj++, store += dim)
		ebder[jj] = ebcoef[store];

	    vector<double> mcoef;
	    s1927(ew1, nur, kk, ed, ew2, inrc, ew3,
		   inlr, mcoef, ebder);

	    store = inbpnt * dim * kl + dim1;
	    for (jj = 0; jj < nn; jj++, store += dim)
		ebcoef[store] = mcoef[jj];

	}
    }

    return;
}


//==========================================================================
inline void s1890(const vector<double>& oknots, int oik, int oin,
		  vector<double>& par, vector<int>& der)
//==========================================================================
{
    int count1, count2;		// Loop control variables
    int start, stop;
    int numb;			// Number of wrong parameters

    double pvl;			// Single parameter value
    double delta;		// Used for correcting wrong
				// parameter values 


    par.resize(oin);
    der.resize(oin);
    fill(der.begin(), der.end(), 0);

    // P R O D U C E  P A R A M E T E R   V A L U E S.
    // First we produce parameter values by a simple algorithm. The
    // parameter values calculated in a wrong way are then corrected.

    par[0] = oknots[oik - 1];
    par[oin - 1] = oknots[oin];
  
    for (count1 = 2; count1 < oin; count1++) {
	stop = count1 + oik;
	double sum = 0.0;
	for (count2 = count1; count2 <= stop; count2++)
	    sum += oknots[count2 - 1];
	par[count1 - 1] = sum / (oik + 1);
    }

    // Find second distinct knot value
    pvl = oknots[oik - 1];
    for (count1 = oik; oknots[count1] <= pvl; count1++) ;

    // Find number of parameter values with wrong value at start of curve
    pvl = (oknots[oik - 1] + oknots[count1]) / 2.0;
    for (numb = 0, start = 1; par[start] <= pvl; start++, numb++);

    if (numb > 0) {
	delta = (pvl - par[0]) / (numb + 1);
	// Fill inn missing parameter values
	pvl = par[0] + delta;
	for (count1 = 1; count1 <= numb; count1++) {
	    par[count1] = pvl;
	    pvl += delta;
	}
    }

    // Find last but one distinct knot value
    pvl = oknots[oin];
    for (count1 = oin - 1; oknots[count1] >= pvl; count1--);

    // Find end parameters in wrong interval
    pvl = (oknots[count1] + oknots[oin + 1]) / 2.0;
    for (numb = 0, stop = oin - 2; par[stop] >= pvl; stop--, numb++);
    if (numb > 0) {
	delta = (par[oin - 1] - pvl) / (numb + 1);
	pvl = par[oin - 1] - delta;
	for (count1 = 1; count1 <= numb; count1++) {
	    par[oin - 1 - count1] = pvl;
	    pvl -= delta;
	}
    }

    return;
}


//==========================================================================
inline void s1891(vector<double>& etau, vector<double>& epoint,
		  int idim, int inbpnt, int iright,
		  vector<int>& eder, vector<double>& et,
		  vector<double>& ebcoef, int& in, int ik,
		  int inlr, int inrc)
//==========================================================================
{
    // Indicate dimension of B-spline
    in = inbpnt; // + ik - 1;

    ebcoef.resize(in * idim * iright);
    fill(ebcoef.begin(), ebcoef.end(), 0.0);

    int nur = inbpnt - inlr;

    // Allocate arrays ew1, ew2, ew3, ed

    int inlx = (1 >= inlr ? 1 : inlr);
    int inrx = (1 >= inrc ? 1 : inrc);
    int limit1 = (ik * nur) + (inrx * nur) + (inlx * inbpnt);
  
    vector<double> ewarray(limit1 + 1, 0.0);
  
    vector<double> ew1 = ewarray;
    vector<double> ew2(ew1.begin() + (ik * nur), ew1.end());
    vector<double> ew3(ew2.begin() + (inrx * nur), ew2.end());

    vector<int> ed(nur, 0);
  
    s1925(etau, epoint, inbpnt, eder, et, ebcoef,
	  in, ik, iright, idim, ew1, nur, ed, ew2, inrc, ew3, inlr);

    return;
}


} // Anonymous namespace


//==========================================================================
SplineSurface* SplineSurface::derivSurface(int ider1, int ider2) const
//==========================================================================
{
    ALWAYS_ERROR_IF(ider1 < 0 || ider2 < 0,
		    "Trying to take a negative number of derivatives!");

    SplineSurface* result;

    int k1, k2, n1, n2;
    int i, j;
    int dim = dimension();

    if (rational()) {

      	// NURBS surface.
	int rdim = dim + 1;

	// Make the denominator surface.
	int ncoefs = numCoefs_u() * numCoefs_v();
	vector<double> ratcoef(ncoefs);
	for (j = 0; j < ncoefs ; ++j) {
	    ratcoef[j] = rcoefs_[(j+1) * rdim - 1];
	}
	SplineSurface rat(numCoefs_u(), numCoefs_v(),
			    order_u(), order_v(),
			    basis_u_.begin(), basis_v_.begin(),
			    ratcoef.begin(), 1);

	// Make resolution for testing of knot equality. 
	const double rel_par_res = 3.0e-8;
	double eps = fabs(basis_u_.begin()[numCoefs_u()]
			  - basis_u_.begin()[order_u() - 1]) * rel_par_res;

	// Make the new knot vector st1 in first direction. 
	k1 = (ider1 + ider2 + 1) * order_u() - (ider1 + ider2);
	vector<double> st1((2 + numCoefs_u() - order_u()) * k1);

	n1 = 0;
	int newcurr;
	for (newcurr = 0; newcurr < k1; ++newcurr) {
	    st1[newcurr] = basis_u_.begin()[order_u() - 1];
	    n1++;
	}

	int oldcurr = order_u();
	int oldprev = oldcurr;
	double limit = basis_u_.begin()[order_u()] - eps;
	int multadd = (ider1 + ider2) * order_u() - ider2;
	int mult = 0;
	while (basis_u_.begin()[oldcurr] < limit) {
	    mult = 0;
	    while ((basis_u_.begin()[oldcurr]
		    - basis_u_.begin()[oldprev]) < eps) {
		oldcurr++;
		mult++;
	    }
	    mult += multadd;
	    if (mult > k1) {
		mult = k1;
	    }
	    for (j = 0; j < mult; j++) {
		st1[newcurr + j] = basis_u_.begin()[oldprev];
		n1++;
	    }
	    newcurr += mult;
	    oldprev = oldcurr;
	}
	for (j = 0; j < k1; j++) {
	    st1[newcurr + j] = basis_u_.begin()[numCoefs_u()];
	}

	// Resize new knot vector st1
	st1.resize(n1 + k1);

	// Calculate parameter values and derivate indicators. 
	//	int kstat = 0;	

	vector<double> par1;
	vector<int> der1;
	s1890(st1, k1, n1, par1, der1);

	// Knot vector in second parameter direction 

	// Make resolution for testing of knot equality. 
	eps = fabs(basis_v_.begin()[numCoefs_v()]
		   - basis_v_.begin()[order_v() - 1]) * rel_par_res;


	// Make the new knot vector st2. 
	k2 = (ider1 + ider2 + 1) * order_v() - (ider1 + ider2);
	vector<double> st2((2 + numCoefs_v() - order_v()) * k2);

	n2 = 0;
	for (newcurr = 0; newcurr < k2; newcurr++) {
	    st2[newcurr] = basis_v_.begin()[order_v() - 1];
	    n2++;
	}

	multadd = (ider1 + ider2) * order_v() - ider1;
	oldcurr = order_v();
	oldprev = oldcurr;
	limit = basis_v_.begin()[numCoefs_v()] - eps;
	while (basis_v_.begin()[oldcurr] < limit) {
	    mult = 0;
	    while ((basis_v_.begin()[oldcurr]
		    - basis_v_.begin()[oldprev]) < eps) {
		oldcurr++;
		mult++;
	    }
	    mult += multadd;
	    if (mult > k2) {
		mult = k2;
	    }
	    for (j = 0; j < mult; j++) {
		st2[newcurr + j] = basis_v_.begin()[oldprev];
		n2++;
	    }
	    newcurr += mult;
	    oldprev = oldcurr;
	}

	for (j = 0; j < k2; j++)
	    st2[newcurr + j] = basis_v_.begin()[numCoefs_v()];

	// Resize new knot vector st2
	st2.resize(n2 + k2);

	// Calculate parameter values and derivate indicators. 
	vector<double> par2;
	vector<int> der2;
	s1890(st2, k2, n2, par2, der2);


	// ------------------------------- 
	// Calculate interpolation points. 
	// ------------------------------- 

	vector<double> tau(rdim * n1 * n2);

	int m = 0;
	Point denom_pt;
// 	int maxder = (ider1 >= ider2 ? ider1 : ider2);
	int sumder = ider1 + ider2;
// 	int numder = (maxder+1) * (maxder+2) / 2;
// 	vector<Point> pts(numder);
	for (j = 0; j < n1; j++) {
	    for (int l = 0; l < n2; l++) {
		rat.point(denom_pt, par1[j], par2[l]);
		double denom = pow(denom_pt[0], (ider1 + ider2 + 1));
		vector<Point> pts =
		    ParamSurface::point(par1[j], par2[l], sumder); //maxder);
		for (i = 0; i < dim; i++) {
		    tau[m++] =
			pts[(ider1+ider2)*(ider1+ider2+1)/2 + ider2][i]*denom;
// 		    tau[m++] = pts[(ider1+ider2+1)/2 + ider1][i] * denom;
		}
		tau[m++] = denom;
	    }
	}

	// Solve the interpolation equation in the second parameter
	// direction.
	vector<double> scoef;
	int n;
	int inlr = 0;
	int inrc = 0;
	s1891(par2, tau, rdim, n2, n1, der2, st2,
		  scoef, n, k2, inlr, inrc);

	// Transpose scoef and put in tau. 
	m = 0;
	for (j = 0; j < n2; j++) {
	    int indx1 = j * rdim;
	    for (int l = 0; l < n1; l++) {
		int indx2 = l * n2 * rdim;
		for (i = 0; i < rdim; i++) {
		    tau[m++] = scoef[indx2 + indx1 + i];
		}
	    }
	}

	// Solve the interpolation equation in the first parameter
	// direction.
	s1891(par1, tau, rdim, n1, n2, der1, st1,
		  scoef, n, k1, inlr, inrc);
	
	result = new SplineSurface(n1, n2, k1, k2,
				     st1.begin(), st2.begin(),
				     scoef.begin(), dim, true);

    } else {

	// Not NURBS. 

	k1 = order_u();
	n1 = numCoefs_u();

	// Create curve representing the surface a a curve in the
	// second parameter direction, copy input arrays 
	int kdim = numCoefs_u() * dim;
	SplineCurve qc1(numCoefs_v(), order_v(),
			  basis_v_.begin(), coefs_begin(), kdim);

	// Make the derivative in the second parameter direction 
	shared_ptr<SplineCurve> pqc2(qc1.derivCurve(ider2));
	SplineCurve qc2 = *pqc2;

	// Remember new knot vector in second parameter direction 
	k2 = qc2.order();
	n2 = qc2.numCoefs();
	vector<double> st2;
	st2.insert(st2.end(), qc2.basis().begin(), qc2.basis().end());

	// Allocate space for turned parameter directions 
	vector<double> coefs((n1*n2 + n2*ider1) * dim);

	// Turn parameter directions
	for (i = 0; i < n1; i++) {
	    for (j = 0; j < n2; j++) {
		for (int k = 0; k < dim; k++) {
		    coefs[(i*n2 + j)*dim + k]
			= qc2.coefs_begin()[(j*n1 + i)*dim + k];
		}
	    }
	}

	// Represent the surface as curve using the first knot vector 
	kdim = n2 * dim;
	SplineCurve qc3(numCoefs_u(), order_u(),
			  basis_u().begin(), coefs.begin(), kdim);

	// Make the derivative in the first parameter direction 
	shared_ptr<SplineCurve> pqc4(qc3.derivCurve(ider1));
	SplineCurve qc4 = *pqc4;

	// Remember new knot vector in first parameter direction 
	k1 = qc4.order();
	n1 = qc4.numCoefs();
	vector<double> st1;
	st1.insert(st1.end(), qc4.basis().begin(), qc4.basis().end());

	// Turn parameter directions of coefficients to match surface 
	for (i = 0; i < n2; i++) {
	    for (j = 0; j < n1; j++) {
		for (int k = 0; k < dim; k++) {
		    coefs[(i*n1 + j)*dim + k]
			= qc4.coefs_begin()[(j*n2 + i)*dim + k];
		}
	    }
	}
	// Create surface object containing the differentiated of the
	// surface
	result = new SplineSurface(n1, n2, k1, k2,
				     st1.begin(), st2.begin(),
				     coefs.begin(), dim);
    }
    return result;
}




} // namespace Go;
