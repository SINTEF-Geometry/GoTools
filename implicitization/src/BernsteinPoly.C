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

#include "GoTools/implicitization/BernsteinPoly.h"
#include "GoTools/implicitization/Binomial.h"
#include "GoTools/utils/errormacros.h"
#include <cmath>
#include <algorithm>


using namespace std;


namespace Go {


//===========================================================================
double BernsteinPoly::operator() (double t) const
//===========================================================================
{
    int d = degree();
    vector<double> coefs = coefs_;

    // The de Casteljau algorithm
    double t1 = 1.0 - t;
    for (int n = d; n > 0; --n) {
	for (int i = 0; i < n; ++i) {
	    coefs[i] = t1 * coefs[i] + t * coefs[i+1];
	}
    }

    return coefs[0];
}


//===========================================================================
bool BernsteinPoly::isZero(const double eps) const
//===========================================================================
{
    for (int i = 0; i < int(coefs_.size()); ++i) {
	if (fabs(coefs_[i]) > eps) {
	    return false;
	}
    }
    return true;
}


//==========================================================================
bool BernsteinPoly::isStrictlyPositive(const double eps) const
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
double BernsteinPoly::norm() const
//==========================================================================
{
    double res = 0.0;
    for (int i = 0; i < int(coefs_.size()); ++i)
        res += fabs(coefs_[i]);

    return res /= (int)coefs_.size();
}


//===========================================================================
void BernsteinPoly::normalize()
//===========================================================================
{
    double n = norm();
    for (int i = 0; i < int(coefs_.size()); ++i)
        coefs_[i] /= n;
    return;
}


//===========================================================================
BernsteinPoly BernsteinPoly::deriv(int der) const
//===========================================================================
{
    ALWAYS_ERROR_IF(der < 0,
		"Differentiating a BernsteinPoly " 
		    "a negative number of times.");


    if (der == 0)
	return *this;

    if (der > degree())
	return BernsteinPoly(0.0);

    int d = degree();
    static vector<double> coefs;
    coefs = coefs_;
    for (int n = 0; n < der; ++n) {
	for (int i = 0; i < d; ++i) {
	    coefs[i] = coefs[i+1] - coefs[i];
	    coefs[i] *= d;
	}
	--d;
    }
    coefs.erase(coefs.end()-der, coefs.end());

    return BernsteinPoly(coefs.begin(), coefs.end());
}


//===========================================================================
double BernsteinPoly::integral(double a, double b) const
//===========================================================================
{
    vector<double> coefs;
    double len = 1.0;
    if (a != 0.0 || b != 1.0) {
	coefs = pickInterval(a, b).coefs_;
	len = fabs(b - a);
    } else {
	coefs = coefs_;
    }

    double res = 0.0;
    for (int i = 0; i < int(coefs.size()); ++i)
	res += coefs[i];

    return res * len / (int)coefs.size();
}


//===========================================================================
double BernsteinPoly::blossom(const vector<double>& tvec) const
//===========================================================================
{
    ALWAYS_ERROR_IF(int(tvec.size()) != degree(),
		    "Vector of arguments must have size degree().");


    int d = degree();
    vector<double> coefs = coefs_;
    for (int n = d; n > 0; --n) {
	double b = tvec[d-n];
	double a = 1.0 - b;
	for (int i = 0; i < n; ++i) {
	    coefs[i] = a * coefs[i] + b * coefs[i+1];
	}
    }

    return coefs[0];
}


//===========================================================================
BernsteinPoly BernsteinPoly::pickInterval(double a, double b) const
//===========================================================================
{
    int deg = degree();
    vector<double> coefs(deg + 1);
    vector<double> tvec(deg, a);
    coefs[0] = blossom(tvec);
    for (int i = 1; i <= deg; ++i) {
	tvec[deg-i] = b;
	coefs[i] = blossom(tvec);
    }

    return BernsteinPoly(coefs.begin(), coefs.end());
}


//===========================================================================
void BernsteinPoly::degreeElevate(int d)
//===========================================================================
{
    *this *= BernsteinPoly(vector<double>(d+1, 1.0));

    return;
}


//==========================================================================
BernsteinPoly& BernsteinPoly::operator*= (const BernsteinPoly& poly)
//==========================================================================
{
    // Degrees
    int m = degree();
    int n = poly.degree();
    int N = m + n;

    // Coefficients
    static vector<double> p;
    p.resize(m+1);
    static vector<double> q;
    q.resize(n+1);

    typedef vector<double>::iterator iter;
    typedef vector<double>::const_iterator const_iter;

    static Binomial binom;

    // Preprocessing coefficients by multiplying binomial coefs
    iter pt = p.begin();
    iter ct = coefs_.begin();
    iter bin_m_i = binom[m];
    int i;
    for (i = 0; i <= m; ++i) {
	*pt = *bin_m_i;
	*pt *= *ct;
	++pt;
	++ct;
	++bin_m_i;
    }
    iter qt = q.begin();
    const_iter ppt = poly.coefs_.begin();
    iter bin_n_i = binom[n];
    for (i = 0; i <= n; ++i) {
	*qt = *bin_n_i * *ppt;
	++qt;
	++ppt;
	++bin_n_i;
    }

    // Summing over products
    coefs_.assign(N + 1, 0.0);
    ct = coefs_.begin();
    iter rt = ct;

    pt = p.begin();
    for (i = 0; i <= m; ++i) {
	ct = rt;
	qt = q.begin();
	for (int j = 0; j <= n; ++j) {
	    *ct += *pt * *qt;
	    ++ct;
	    ++qt;
	}
	++pt;
	++rt;
    }

    // Postprocessing with more binomial coefs
    ct = coefs_.begin();
    iter bin_N_i = binom[N];
    for (i = 0; i <= N; ++i) {
	*ct /= *bin_N_i;
	++ct;
	++bin_N_i;
    }

    return *this;
}


//===========================================================================
BernsteinPoly& BernsteinPoly::operator*= (double c)
//===========================================================================
{
    for (int i = 0; i < int(coefs_.size()); ++i)
	coefs_[i] *= c;
    return *this;
}


//===========================================================================
BernsteinPoly& BernsteinPoly::operator+= (const BernsteinPoly& poly)
//===========================================================================
{
    int m = degree();
    int n = poly.degree();

    int maxdeg = max(m, n);

    static BernsteinPoly tmp;
    tmp = poly;

    if (m < maxdeg)
	degreeElevate(maxdeg-m);
    if (n < maxdeg)
	tmp.degreeElevate(maxdeg-n);

    for (int i = 0; i <= maxdeg; ++i)
	coefs_[i] += tmp.coefs_[i];

    return *this;
}


//===========================================================================
BernsteinPoly& BernsteinPoly::operator+= (double c)
//===========================================================================
{
    *this += BernsteinPoly(c);
    return *this;
}


//===========================================================================
BernsteinPoly& BernsteinPoly::operator-= (const BernsteinPoly& poly)
//===========================================================================
{
    *this += ((-1.0) * poly);
    return *this;
}


//===========================================================================
BernsteinPoly& BernsteinPoly::operator-= (double c)
//===========================================================================
{
    *this += -c;
    return *this;
}


//===========================================================================
BernsteinPoly& BernsteinPoly::operator/= (double c)
//===========================================================================
{
    for (int i = 0; i < int(coefs_.size()); ++i)
	coefs_[i] /= c;
    return *this;
}


//===========================================================================
void BernsteinPoly::read(istream& is)
//===========================================================================
{
    int deg;
    is >> deg;
    coefs_.resize(deg + 1);
    for (int i = 0; i < int(coefs_.size()); ++i) {
	is >> coefs_[i];
    }
}

//==========================================================================
void BernsteinPoly::write(ostream& os) const
//==========================================================================
{
    os << degree() << '\n';
    os << coefs_[0];
    for (int i = 1; i < int(coefs_.size()); ++i) {
	os << ' ' << coefs_[i];
    }
}


//==========================================================================


} // namespace Go
