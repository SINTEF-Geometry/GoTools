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

#include "GoTools/implicitization/BernsteinMulti.h"
#include "GoTools/implicitization/BernsteinPoly.h"
#include "GoTools/implicitization/Binomial.h"
#include "GoTools/utils/errormacros.h"
#include <algorithm>


using namespace std;


namespace Go {


//===========================================================================
double BernsteinMulti::operator() (double u, double v) const
//===========================================================================
{
    vector<double> coefs = coefs_;

    // The de Casteljau algorithm
    double u1 = 1.0 - u;
    for (int nu = degu_; nu > 0; --nu) {
	int ind = 0;
	int jnd = 0;
	for (int iv = 0; iv <= degv_; ++iv) {
	    for (int iu = 0; iu < nu; ++iu) {
		coefs[ind] = u1 * coefs[jnd] + u * coefs[jnd+1];
		++ind;
		++jnd;
	    }
	    ++jnd;
	}
    }
    double v1 = 1.0 - v;
    for (int nv = degv_; nv > 0; --nv) {
	for (int iv = 0; iv < nv; ++iv) {
	    coefs[iv] = v1 * coefs[iv] + v * coefs[iv+1];
	}
    }

    return coefs[0];
}


//==========================================================================
bool BernsteinMulti::isZero(const double eps) const
//==========================================================================
{
    for (int i = 0; i < int(coefs_.size()); ++i) {
	if (fabs(coefs_[i]) > eps) {
	    return false;
	}
    }
    return true;
}


//==========================================================================
bool BernsteinMulti::isStrictlyPositive(const double eps) const
//==========================================================================
{
    for (int i = 0; i < int(coefs_.size()); ++i) {
	if (coefs_[i] <= eps) {
	    return false;
	}
    }
    return true;
}


//==========================================================================
bool BernsteinMulti::isStrictlyNegative(const double eps) const
//==========================================================================
{
    for (int i = 0; i < int(coefs_.size()); ++i) {
	if (coefs_[i] >= -eps) {
	    return false;
	}
    }
    return true;
}


//==========================================================================
bool BernsteinMulti::isNonNegative(const double eps) const
//==========================================================================
{
    for (int i = 0; i < int(coefs_.size()); ++i) {
	if (coefs_[i] < -eps) {
	    return false;
	}
    }
    return true;
}


//==========================================================================
double BernsteinMulti::norm() const
//==========================================================================
{
    // This defines the L_1 norm  <- not quite correct, norm() is a
    // bound for the L_1 norm, but not identical to it (if the
    // polynomial changes sign).
    double res = 0.0;
    for (int i = 0; i < int(coefs_.size()); ++i)
	res += fabs(coefs_[i]);

    return res /= (int)coefs_.size();
}


//==========================================================================
void BernsteinMulti::normalize()
//==========================================================================
{
    double n = norm();
    for (int i = 0; i < int(coefs_.size()); ++i)
	coefs_[i] /= n;
    return;
}


//==========================================================================
double BernsteinMulti::mean() const
//==========================================================================
{
    double res = 0.0;
    for (int i = 0; i < int(coefs_.size()); ++i)
	res += coefs_[i];

    return res /= (int)coefs_.size();
}


//==========================================================================
BernsteinMulti BernsteinMulti::deriv(int der1, int der2) const
//==========================================================================
{
    ALWAYS_ERROR_IF(der1 < 0 || der2 < 0,
		"Differentiating a BernsteinMulti " 
		    "a negative number of times.");


    if (der1 == 0 && der2 == 0)
	return *this;

    if (der1 > degu_ || der2 > degv_)
	return BernsteinMulti(0.0);

    int du = degu_;
    int dv = degv_;
    static vector<double> coefs;
    coefs = coefs_;

    // Differentiate in v-direction
    for (int n = 0; n < der2; ++n) {
	int ind = 0;
	for (int j = 0; j < dv; ++j) {
	    for (int i = 0; i <= du; ++i) {
		coefs[ind] = coefs[ind+du+1] - coefs[ind];
		coefs[ind] *= dv;
		++ind;
	    }
	}
	--dv;
    }
    // Differentiate in u-direction
    for (int m = 0; m < der1; ++m) {
	int ind = 0;
	int jnd = 0;
	for (int j = 0; j <= dv; ++j) {
	    for (int i = 0; i < du; ++i) {
		coefs[ind] = coefs[jnd+1] - coefs[jnd];
		coefs[ind] *= du;
		++ind;
		++jnd;
	    }
	    ++jnd;
	}
	--du;
    }
    coefs.erase(coefs.end()-((degv_+1)*der1+(degu_+1-der1)*der2),
		coefs.end());

    return BernsteinMulti(du, dv, coefs.begin(), coefs.end());
}


//===========================================================================
BernsteinMulti BernsteinMulti::detHess() const
//===========================================================================
{
    // Describes Gaussian curvature of the graph of the polynomial
    return deriv(2, 0) * deriv(0, 2) - deriv(1, 1) * deriv(1, 1);
}


//===========================================================================
BernsteinMulti BernsteinMulti::traceHess() const
//===========================================================================
{
    // Describes mean curvature (modulo a factor 2) of the graph of
    // the polynomial
    return deriv(2, 0) + deriv(0, 2);
}

//==========================================================================
double BernsteinMulti::blossom(const vector<double>& uvec,
 			       const vector<double>& vvec) const
//==========================================================================
{
    int du = degreeU();
    int dv = degreeV();

    ALWAYS_ERROR_IF(int(uvec.size()) != du || int(vvec.size()) != dv,
		    "Vectors of arguments must have size degreeU() "
		    "and degreeV().");


    vector<double> coefs = coefs_;

    // First take sums in the u-directions ...
    for (int nu = du; nu > 0; --nu) {
	double b = uvec[du-nu];
	double a = 1.0 - b;
	int ind = 0;
	int jnd = 0;
	for (int iv = 0; iv <= dv; ++iv) {
	    for (int iu = 0; iu < nu; ++iu) {
		coefs[ind] = a * coefs[jnd] + b * coefs[jnd+1];
		++ind;
		++jnd;
	    }
	    ++jnd;
	}
    }
    // ... then in the v-direction
    for (int nv = dv; nv > 0; --nv) {
	double b = vvec[dv-nv];
	double a = 1.0 - b;
	for (int iv = 0; iv < nv; ++iv) {
	    coefs[iv] = a * coefs[iv] + b * coefs[iv+1];
	}
    }

    return coefs[0];

}


//===========================================================================
BernsteinMulti BernsteinMulti::pickDomain(double u0, double u1,
					  double v0, double v1) const
//===========================================================================
{
    vector<double> coefs((degu_+1) * (degv_+1));

    // Use blossoming to calculate new control points
    vector<double> uvec(degu_);
    vector<double> vvec(degv_, v0);
    int ind = 0;
    for (int j = 0; j <= degv_; ++j) {
	fill(uvec.begin(), uvec.end(), u0);
	if (j != 0)
	    vvec[degv_-j] = v1;
	for (int i = 0; i <= degu_; ++i) {
	    if (i != 0)
		uvec[degu_-i] = u1;
	    coefs[ind] = blossom(uvec, vvec);
	    ++ind;
	}
    }

    return BernsteinMulti(degu_, degv_, coefs.begin(), coefs.end());
}


//===========================================================================
BernsteinPoly BernsteinMulti::pickLine(Array<double,2> a,
				       Array<double,2> b) const
//===========================================================================
{
    typedef vector<double>::iterator iter;

    int du = degreeU();
    int dv = degreeV();
    int D = du + dv;

    // First pick out the domain with corners at (au,av) and
    // (bu,bv). The final line will then be the line with enpoints on
    // these corners.
    BernsteinMulti tmp = pickDomain(a[0], b[0], a[1], b[1]);

    static Binomial binom;

    // Preprocessing coefficients by multiplying binomial coefs
    iter pt = tmp.coefs_.begin();
    iter bin_dv_j = binom[dv];
    int j;
    for (j = 0; j <= dv; ++j) {
	iter bin_du_i = binom[du];
	for (int i = 0; i <= du; ++i) {
	    *pt *= *bin_du_i * *bin_dv_j;
	    ++pt;
	    ++bin_du_i;
	}
	++bin_dv_j;
    }

    // Calculating new coefficients
    static vector<double> coefs(50);
    coefs.resize((D+1));
    iter ct = coefs.begin();
    fill(ct, coefs.end(), 0.0);
    iter rt = ct;

    pt = tmp.coefs_.begin();
    for (j = 0; j <= dv; ++j) {
	ct = rt;
	for (int i = 0; i <= du; ++i) {
	    *ct += *pt;
	    ++ct;
	    ++pt;
	}
	++rt;
    }

    // Postprocessing with more binomial coefs
    ct = coefs.begin();
    iter bin_D_i = binom[D];
    for (int i = 0; i <= D; ++i) {
	*ct /= *bin_D_i;
	++ct;
	++bin_D_i;
    }

    return BernsteinPoly(coefs.begin(), coefs.end());
}


//===========================================================================
void BernsteinMulti::degreeElevate(int du, int dv)
//===========================================================================
{
    *this *= BernsteinMulti(du, dv, vector<double>((du+1) * (dv+1), 1.0));

    return;
}


//===========================================================================
BernsteinMulti&
BernsteinMulti::operator*= (const BernsteinMulti& multi)
//===========================================================================
{
    // Orders of left multi
    int mu = degu_;
    int mv = degv_;
    // Orders of right multi
    int nu = multi.degu_;
    int nv = multi.degv_;
    // Total orders of product multi
    int Nu = mu + nu;
    int Nv = mv + nv;

    // Coefficient vectors for *this and multi
    static vector<double> p;
    p.resize((mu+1) * (mv+1));
    static vector<double> q;
    q.resize((nu+1) * (nv+1));

    typedef vector<double>::iterator iter;
    typedef vector<double>::const_iterator const_iter;

    static Binomial binom;
    binom(Nu > Nv ? Nu : Nv, 0);

    // Preprocessing the coefficients by multiplying in binomial
    // coefficients
    iter pt = p.begin();
    iter ct = coefs_.begin();
    Binomial::iter bin_mv_iv = binom[mv];
    int iv;
    for (iv = 0; iv <= mv; ++iv) {
	Binomial::iter bin_mu_iu = binom[mu];
	for (int iu = 0; iu <= mu; ++iu) {
	    *pt = *bin_mu_iu * *bin_mv_iv * *ct;
	    ++pt;
	    ++ct;
	    ++bin_mu_iu;
	}
	++bin_mv_iv;
    }
    iter qt = q.begin();
    const_iter mct = multi.coefs_.begin();
    Binomial::iter bin_nv_iv = binom[nv];
    for (iv = 0; iv <= nv; ++iv) {
	Binomial::iter bin_nu_iu = binom[nu];
	for (int iu = 0; iu <= nu; ++iu) {
	    *qt = *bin_nu_iu * *bin_nv_iv * *mct;
	    ++qt;
	    ++mct;
	    ++bin_nu_iu;
	}
	++bin_nv_iv;
    }

    // Sum over products
    coefs_.resize((Nu+1) * (Nv+1));
    ct = coefs_.begin();
    fill(ct, coefs_.end(), 0.0);
    iter rt = ct;

    pt = p.begin();
    for (iv = 0; iv <= mv; ++iv) {
	for (int iu = 0; iu <= mu; ++iu) {
	    ct = rt;
	    qt = q.begin();
	    for (int jv = 0; jv <= nv; ++jv) {
		for (int ju = 0; ju <= nu; ++ju) {
		    *ct += *pt * *qt;
		    ++ct;
		    ++qt;
		}
		if (jv < nv) {
		    ct += Nu - nu;
		}
	    }
	    ++pt;
	    ++rt;
	}
	if (iv < mv) {
	    rt += Nu - mu;
	}
    }

    // Postprocessing with more binomial coefs
    ct = coefs_.begin();
    Binomial::iter bin_Nv_iv = binom[Nv];
    for (iv = 0; iv <= Nv; ++iv) {
	Binomial::iter bin_Nu_iu = binom[Nu];
	for (int iu = 0; iu <= Nu; ++iu) {
	    *ct /= *bin_Nu_iu * *bin_Nv_iv;
	    ++ct;
	    ++bin_Nu_iu;
	}
	++bin_Nv_iv;
    }

    degu_ = Nu;
    degv_ = Nv;

    return *this;
}


//===========================================================================
BernsteinMulti& BernsteinMulti::operator*= (double c)
//===========================================================================
{
    for (int i = 0; i < int(coefs_.size()); ++i)
	coefs_[i] *= c;
    return *this;
}


//==========================================================================
BernsteinMulti&
BernsteinMulti::operator+= (const BernsteinMulti& multi)
//==========================================================================
{
    int mu = degu_;
    int mv = degv_;
    int nu = multi.degu_;
    int nv = multi.degv_;

    int maxu = max(mu, nu);
    int maxv = max(mv, nv);

    static BernsteinMulti tmp;
    tmp = multi;

    if (mu < maxu || mv < maxv)
	degreeElevate(maxu-mu, maxv-mv);
    if (nu < maxu || nv < maxv)
	tmp.degreeElevate(maxu-nu, maxv-nv);

    for (int i = 0; i < int(coefs_.size()); ++i)
	coefs_[i] += tmp.coefs_[i];

    return *this;
}


//===========================================================================
BernsteinMulti& BernsteinMulti::operator+= (double c)
//===========================================================================
{
    *this += BernsteinMulti(c);
    return *this;
}


//===========================================================================
BernsteinMulti&
BernsteinMulti::operator-= (const BernsteinMulti& multi)
//===========================================================================
{
    *this += -1.0 * multi;
    return *this;
}


//===========================================================================
BernsteinMulti& BernsteinMulti::operator-= (double c)
//===========================================================================
{
    *this += -c;
    return *this;
}


//===========================================================================
BernsteinMulti& BernsteinMulti::operator/= (double c)
//===========================================================================
{
    for (int i = 0; i < int(coefs_.size()); ++i)
	coefs_[i] /= c;
    return *this;
}


//===========================================================================
BernsteinPoly BernsteinMulti::bindU(double u) const
//===========================================================================
{
    double one_minus_u = 1.0 - u;
    int du = degreeU();
    int dv = degreeV();
    vector<double> tmp = coefs_;
    vector<double> newcoefs(dv + 1);
    for (int i = 0; i < dv+1; ++i) {
	for (int nu = du; nu > 0; --nu) {
	    for (int iu = 0; iu < nu; ++iu) {
		tmp[(du+1)*i + iu]
		    = one_minus_u * tmp[(du+1)*i + iu]
		    + u * tmp[(du+1)*i + iu + 1];
	    }
	}
	newcoefs[i] = tmp[(du+1)*i];
    }

    return BernsteinPoly(newcoefs);
}


//===========================================================================
BernsteinPoly BernsteinMulti::bindV(double v) const
//===========================================================================
{
    double one_minus_v = 1.0 - v;
    int du = degreeU();
    int dv = degreeV();
    vector<double> tmp = coefs_;
    vector<double> newcoefs(du + 1);
    for (int i = 0; i < du+1; ++i) {
	for (int nv = dv; nv > 0; --nv) {
	    for (int iv = 0; iv < nv; ++iv) {
		tmp[(du+1)*iv + i]
		    = one_minus_v * tmp[(du+1)*iv + i]
		    + v * tmp[(du+1)*(iv+1) + i];
	    }
	}
	newcoefs[i] = tmp[i];
    }

    return BernsteinPoly(newcoefs);
}


//===========================================================================
void BernsteinMulti::read(istream& is)
//===========================================================================
{
//     int degu, degv;
    is >> degu_ >> degv_;
    coefs_.resize((degu_+1) * (degv_+1));
    for (int i = 0; i < int(coefs_.size()); ++i) {
	is >> coefs_[i];
    }
}


//===========================================================================
void BernsteinMulti::write(ostream& os) const
//===========================================================================
{
    os << degu_ << '\n' << degv_ << '\n';
    os << coefs_[0];
    for (int i = 1; i < int(coefs_.size()); ++i) {
	os << ' ' << coefs_[i];
    }
}


//===========================================================================


} // namespace Go

