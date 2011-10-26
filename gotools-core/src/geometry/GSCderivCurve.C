//===========================================================================
//                                                                           
// File: GSCderivCurve.C                                                     
//                                                                           
// Created: Tue Aug 14 11:09:03 2001                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@math.sintef.no>
//                                                                           
// Revision: $Id: GSCderivCurve.C,v 1.14 2005-10-20 11:46:28 sbr Exp $
//                                                                           
// Description: Functions to differentiate a spline curve.
//              Ported from sisl (s1720.c).
//===========================================================================


#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/utils/Point.h"


using std::vector;


namespace Go
{


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
    const int nn = nur + inlr;
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
	    for (jj = 0; jj < nn; jj++, store += dim) {
		ebder[jj] = ebcoef[store];
	    }

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
    // The returned der vector will consist of oin 0.0.
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


//===========================================================================
SplineCurve* SplineCurve::derivCurve(int ider) const
//===========================================================================
{
    ALWAYS_ERROR_IF(ider < 0,
		    "Trying to derivate a spline curve "
		    "a negative number of times!");

    MESSAGE_IF(ider >= order(),
	       "Trying to derivate a spline curve "
	       "more times than there are degrees! "
	       "Returning a curve of order 1.");

    // Maybe not the most interesting assignment, but no reason not to
    // handle it.
    if (ider == 0) return new SplineCurve(numCoefs(), order(),
					  basis().begin(), coefs_begin(),
					  dimension());

    // We run through knots to verify that the inner continuity of the curve is acceptable.
    int max_inner_knot_mult = 0;
    int knot_mult = 0;
    for (int ki = order(); ki < numCoefs(); ki += knot_mult) {
	knot_mult = basis_.knotMultiplicity(basis_.begin()[ki]);
	if (knot_mult > max_inner_knot_mult) {
	    max_inner_knot_mult = knot_mult;
	}
    }
    if (order() - ider - max_inner_knot_mult < 1) {
	// @@sbr The message should only be issued after analysis of the geometric continuity.
//  	MESSAGE("Inner continuity not high enough for the requested derivative curve!");
// 	THROW("Inner continuity not high enough for the requested derivative curve!");
    }

    // @@ Maybe there should be more consistency when it comes to almost 0.
    double dequal = 1e-13;
    double dzero = 0.0;
    int k;                   // Order of the output curve.
    int n;                   // Number of the vertices in output curves.
    int dim = dimension();   // Dimension of the space in which the
                             // curves lies.
    double tdel;             // Help variabel.
    std::vector<double>::iterator s1, s2, s3;      // Pointers used in loop.

    SplineCurve* deriv_curve;
    if (rational()) {
	// NURBS curve.
	int rdim = dim_ + (rational_ ? 1 : 0);
	int num_coefs = numCoefs();
	const double *st = &(basis_.begin()[0]);

	// Make the denominator curve.
	vector<double> rat_coefs(num_coefs, 0.0);
	for (int ki = 0; ki < num_coefs; ++ki) {
	    rat_coefs[ki] = rcoefs_[(ki + 1)*rdim - 1];
	}

	SplineCurve denom_cv = SplineCurve(num_coefs, order(), basis_.begin(),
					   rat_coefs.begin(), 1, false);

	// Make resolution for testing of knot equality.
	double rel_par_res = 3.0e-8;
	double eps = fabs(st[num_coefs] - st[order()-1])*rel_par_res;

	// Make the new knot vector.

	int nmb_new_coefs = 0;
	int kk = (ider + 1)*order() - ider;
	int multadd = ider*order();
	vector<double> new_knots((2 + num_coefs - order())*kk, 0.0);
	int new_curr;
	for (new_curr = 0; new_curr < kk; ++new_curr) {
	    new_knots[new_curr] = st[order()-1];
	    nmb_new_coefs++;
	}

	int old_curr = order();
	int old_prev = old_curr;
	while ((st[old_curr] + eps) < st[num_coefs]) {
	    int mult = 0;
	    while ((st[old_curr] - st[old_prev]) < eps) {
		old_curr++;
		mult++;
	    }
	    mult += multadd;
	    if (mult > kk) {
		mult = kk;
	    }
	    for (int kj = 0; kj < mult; ++kj) {
		new_knots[new_curr + kj] = st[old_prev];
		nmb_new_coefs++;
	    }
	    new_curr += mult;
	    old_prev = old_curr;
	}

	for (int kj = 0; kj < kk; ++kj) {
	    new_knots[new_curr+kj] = st[num_coefs];
	}

	// Calculate parameter values and derivate indicators.
	vector<double> par;
	vector<int> der;
	s1890(new_knots, kk, nmb_new_coefs, par, der);

	// 	s1890(st, kk, nmb_new_coefs, &par, &der, &kstat);
	// 	if (kstat < 0) goto error;

	// Calculate interpolation points.
	vector<double> tau(nmb_new_coefs*rdim, 0.0);
	for (int ki = 0; ki < nmb_new_coefs; ++ ki) {
	    Point denom_pt = denom_cv.ParamCurve::point(par[ki]);
	    double denom = pow(denom_pt[0], (ider + 1));

	    vector<Point> der_pt = ParamCurve::point(par[ki], ider);
	    for (int kj = 0; kj < dim; ++kj) {
		tau[ki*rdim + kj] = der_pt[ider][kj]*denom;
	    }
	    tau[ki*rdim + dim] = denom;
	}

	/* Make the new curve description. */
	vector<double> der_coefs;
	s1891(par, tau, rdim, nmb_new_coefs, 1, der, new_knots,
	      der_coefs, nmb_new_coefs, kk, 0, 0);

	deriv_curve = new SplineCurve(nmb_new_coefs, kk, new_knots.begin(),
				      der_coefs.begin(), dim, true);


#ifdef GEOMETRY_DEBUG
	// We try to see if the interpolation in the homogenuous space was a success.
	SplineCurve hom_curve(nmb_new_coefs, kk, new_knots.begin(),
			      der_coefs.begin(), rdim, false);
	Point tau_pt(4);
	for (size_t ki = 0; ki < par.size(); ++ki) {
	    Point der_pt = hom_curve.ParamCurve::point(par[ki]);
	    tau_pt.setValue(&tau[ki*rdim]);
	    double dist = der_pt.dist(tau_pt);
	    if (dist > 0.01) {
		std::cout << "dist: " << std::endl;
	    }
	}
#endif // GEOMETRY_DEBUG

    } else {
	// Not a rational case.
	if (ider >= order()) {
	    n = numCoefs() + order() -1;
	    k = 1;
	} else {
	    n = numCoefs() + ider;
	    k = order() - ider;
	}

	// Allocating the new arrays to the new curve.
	std::vector<double> new_knots(n + k, 0.0); // The first new knot-vector.
	std::vector<double> new_coefs(n*dim, 0.0);
	std::copy(basis().begin(), basis().end(), // + 
		  new_knots.begin());
	if (ider < order()) {
	    std::copy(coefs_begin(), coefs_end(),
		      new_coefs.begin());
	}

	// Here we are computing a new coefficient vector for each round.
	int j;
	if (ider < order())
	    for (j = 1; j <= ider; ++j) {
		k = order() - j;      // The new order of the curve.

		s1 = new_coefs.begin()
		    + (numCoefs() - 1 + j) * dim; // Last new vertice.

		// Here I just refere to the referenses.
		// The new vertices are computig from back to front,
		// and the "last" dim vertices is being computed outside
		// the main loop because we did not have new_coefs[-1],
		// instead we use zeroes.
		for (s3 = new_knots.begin() + numCoefs() - 1 + j;
		     s3 > new_knots.begin(); --s3, s1 -= 2 * dim) {
		    tdel = s3[k] - *s3;

		    //		if (DNEQUAL(tdel,DZERO))
		    if (fabs(tdel - dzero) > dequal)
			for (s2 = s1 + dim; s1 < s2; s1++)
			    *s1 = (*s1 - s1[-dim]) * k/tdel;
		    else
			for (s2 = s1 + dim; s1 < s2; s1++) *s1 = dzero; //DZERO;
		}

		tdel = s3[k] - *s3;

		//	    if (DNEQUAL(tdel,DZERO))
		if (fabs(tdel - dzero) > dequal)
		    for (s2=s1+dim; s1<s2; s1++) *s1 = *s1*k/tdel;
		else
		    for (s2=s1+dim; s1<s2; s1++) *s1 = dzero; //DZERO;
	    }

	// Allocating new curve-object.

	// Unnecessary interior knots and coefs must be removed prior to
	// curve creation.
	// See s1705.c.
	int nmb_new_coefs = 0;
	for (s1 = new_coefs.begin(), s2 = new_knots.begin(), s3 = s2 + n;
	     s2 < s3; s1 += dim, ++s2)
	    if (s2[k] > *s2) {
		// Here we copy necessary vertices to compress the vector. */
		for (j=0; j<dim; j++)
		    new_coefs[nmb_new_coefs*dim + j] = s1[j];

		// Here we copy necessary knots to compress the vector.
		new_knots[nmb_new_coefs] = *s2;

		// Updating number of copied knots.
		++nmb_new_coefs;
	    }

	// Copy the last k knots.
	for (j=0; j<k; ++j) 
	    new_knots[nmb_new_coefs+j] = s3[j];
    

#ifdef GEOMETRY_DEBUG
	if (nmb_new_coefs + k > (int)(new_knots.size())) {
	    std::cout << "nmb_new_coefs: " << nmb_new_coefs << ", k: " << k << 
		", new_knots.size(): " << new_knots.size() << std::endl;
	}
	if (nmb_new_coefs*dim > (int)(new_coefs.size())) {
	    std::cout << "nmb_new_coefs: " << nmb_new_coefs << 
		", new_coefs.size(): " << new_coefs.size() << std::endl;
	}
#endif // GEOMETRY_DEBUG

	std::vector<double> der_knots(new_knots.begin(), new_knots.begin() + nmb_new_coefs + k);
	std::vector<double> der_coefs(new_coefs.begin(), new_coefs.begin() + nmb_new_coefs*dim);

	deriv_curve = new SplineCurve(nmb_new_coefs, k, der_knots.begin(), der_coefs.begin(), dim);

    }

    return deriv_curve;
}


} // namespace Go;
