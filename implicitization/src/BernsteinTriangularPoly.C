//==========================================================================
//                                                                          
// File: BernsteinTriangularPoly.C                                           
//                                                                          
// Created: Mon Apr 22 09:24:06 2002                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: BernsteinTriangularPoly.C,v 1.21 2007-05-09 17:57:04 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include "GoTools/implicitization/BernsteinTriangularPoly.h"
#include "GoTools/utils/binom.h"
#include "GoTools/utils/errormacros.h"
#include <algorithm>
#include <math.h>


using namespace std;


namespace Go {


//==========================================================================
double BernsteinTriangularPoly::norm() const
//==========================================================================
{
    double res = 0.0;
    for (int i = 0; i < int(coefs_.size()); ++i)
        res += fabs(coefs_[i]);

    return res /= (int)coefs_.size();
}


//===========================================================================
void BernsteinTriangularPoly::normalize()
//===========================================================================
{
    double n = norm();
    for (int i = 0; i < int(coefs_.size()); ++i)
        coefs_[i] /= n;
    return;
}


//===========================================================================
void BernsteinTriangularPoly::deriv(int der, const Vector3D& d,
				    BernsteinTriangularPoly& btp) const
//===========================================================================
{
    // It is the user's responsibility to ensure that d[0]+d[1]+d[2]=0

    ALWAYS_ERROR_IF(der < 0,
		"Differentiating a BernsteinTriangularPoly " 
		    "a negative number of times.");


    if (der == 0) {
	btp = *this;
	return;
    }

    VecType tmp = coefs_;
    for (int r = 1; r <= der; ++r) {
	int m = -1;
	for (int i = 0; i <= deg_-r; ++i) {
	    for (int l = 0; l <= i; ++l) {
		m++;
		tmp[m] = d[0] * tmp[m]
		    + d[1] * tmp[m + 1 + i]
		    + d[2] * tmp[m + 2 + i];
		tmp[m] *= deg_ - r + 1;
	    }
	}
    }

    btp.deg_ = deg_ - der;
    int n = (btp.deg_ + 1) * (btp.deg_ + 2) / 2;
    btp.coefs_.assign(tmp.begin(), tmp.begin()+n);

    return;
}


//===========================================================================
BernsteinTriangularPoly&
BernsteinTriangularPoly::operator*= (const BernsteinTriangularPoly& poly)
//===========================================================================
{
    // Aliases for the old coefficients
    const VecType& p = coefs_;
    const VecType& q = poly.coefs_;
    int psz = (int)coefs_.size();
    int qsz = (int)poly.coefs_.size();

    int d1 = degree();
    int d2 = poly.degree();
    int d = d1 + d2;

    // A vector for the new coefficients
    VecType g((d+1)*(d+2)/2);

    int n, n1, n2;
    for (int i = d; i >= 0; --i) {
	for (int j = d-i; j >= 0; --j) {
	    // n is the index of the linear array given (i,j,k) with
	    // i+j+k=d.
	    n = (d-i+1)*(d-i+2)/2 - j - 1;
	    g[n] = 0.0;
	    for (int i1 = d1; i1 >= 0; --i1) {
		for (int j1 = d1-i1; j1 >= 0; --j1) {
		    n1 = (d1-i1+1)*(d1-i1+2)/2 - j1 - 1;
		    n2 = (d2-(i-i1)+1)*(d2-(i-i1)+2)/2 - (j-j1) - 1;
		    if (n1 < 0 || n2 < 0 || n1 >= psz || n2 >= qsz)
			continue;
		    g[n] += p[n1] * q[n2]
			* trinomial(d1, i1, j1)
			* trinomial(d2, i-i1, j-j1)
			/ trinomial(d, i, j);
		}
	    }
	}
    }

    deg_ = d;
    swap(coefs_, g);
    return *this;
}


//===========================================================================
BernsteinTriangularPoly& BernsteinTriangularPoly::operator*= (double c)
//===========================================================================
{
    for (int i = 0; i < int(coefs_.size()); ++i)
	coefs_[i] *= c;
    return *this;
}


//===========================================================================
BernsteinTriangularPoly&
BernsteinTriangularPoly::operator+= (const BernsteinTriangularPoly& poly)
//===========================================================================
{
    int d = degree();
    int d1 = poly.degree();

    if (d >= d1) {
	double value = 0.0;
	// n and n1 are linear indices - see comment for operator*=
	int n, n1;
	for (int i = d; i >= 0; --i) {
	    for (int j = d-i; j >= 0; --j) {
		value = 0.0;
		for (int i1 = d1; i1 >= 0; --i1) {
		    for (int j1 = d1-i1; j1 >= 0; --j1) {
			n1 = (d1-i1+1)*(d1-i1+2)/2 - j1 - 1;
			value += poly.coefs_[n1]
			    * binom(i, i1) * binom(j, j1)
			    * binom(d-i-j, d1-i1-j1)
			    / binom(d, d1);
		    }
		}
		n = (d-i+1)*(d-i+2)/2 - j - 1;
		coefs_[n] += value;
	    }
	}
    } else {
	BernsteinTriangularPoly tmp = poly;
	tmp += *this;
	*this = tmp;
    }

    return *this;
}


//===========================================================================
BernsteinTriangularPoly& BernsteinTriangularPoly::operator+= (double c)
//===========================================================================
{
    *this += BernsteinTriangularPoly(c);
    return *this;
}


//===========================================================================
BernsteinTriangularPoly&
BernsteinTriangularPoly::operator-= (const BernsteinTriangularPoly& poly)
//===========================================================================
{
    *this += ((-1.0) * poly);
    return *this;
}


//===========================================================================
BernsteinTriangularPoly& BernsteinTriangularPoly::operator-= (double c)
//===========================================================================
{
    *this += -c;
    return *this;
}


//===========================================================================
BernsteinTriangularPoly& BernsteinTriangularPoly::operator/= (double c)
//===========================================================================
{
    for (int i = 0; i < int(coefs_.size()); ++i)
	coefs_[i] /= c;
    return *this;
}


//===========================================================================
void BernsteinTriangularPoly::read(istream& is)
//===========================================================================
{
    bool is_good = is.good();
    if (!is_good) {
	THROW("Invalid input file!");
    }

    is >> deg_;
    int sz = (deg_+1)*(deg_+2)/2;
    coefs_.resize(sz);
    for (int i = 0; i < sz; ++i) {
	is >> coefs_[i];
    }
}

//===========================================================================
void BernsteinTriangularPoly::write(ostream& os) const
//===========================================================================
{
    os << degree() << '\n';
    os << coefs_[0];
    for (int i = 1; i < int(coefs_.size()); ++i) {
	os << ' ' << coefs_[i];
    }
    os << endl;
}


//==========================================================================


} // namespace Go
