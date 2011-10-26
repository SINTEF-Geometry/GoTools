//==========================================================================
//                                                                          
// File: BernsteinTetrahedralPoly.C                                          
//                                                                          
// Created: Tue Jan 21 09:55:36 2003                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: BernsteinTetrahedralPoly.C,v 1.19 2007-05-09 17:57:04 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include "GoTools/implicitization/BernsteinTetrahedralPoly.h"
#include "GoTools/implicitization/BernsteinPoly.h"
#include "GoTools/utils/binom.h"
#include "GoTools/utils/errormacros.h"
#include <algorithm>
#include <math.h>


using namespace std;


namespace Go {


//==========================================================================
double BernsteinTetrahedralPoly::norm() const
//==========================================================================
{
    double res = 0.0;
    for (int i = 0; i < int(coefs_.size()); ++i)
        res += fabs(coefs_[i]);

    return res /= (int)coefs_.size();
}


//===========================================================================
void BernsteinTetrahedralPoly::normalize()
//===========================================================================
{
    double n = norm();
    for (int i = 0; i < int(coefs_.size()); ++i)
        coefs_[i] /= n;
    return;
}


//===========================================================================
void BernsteinTetrahedralPoly::deriv(int der, const Vector4D& d,
				     BernsteinTetrahedralPoly& btp) const
//===========================================================================
{
    // It is the user's responsibility to ensure that d[0]+d[1]+d[2]+d[3]=0

    ALWAYS_ERROR_IF(der < 0,
		"Differentiating a BernsteinTetrahedralPoly " 
		    "a negative number of times.");


    if (der == 0) {
	btp = *this;
	return;
    }

    vector<double> tmp = coefs_;
    for (int r = 1; r <= der; ++r) {
	int m = 0;
	for (int i = 0; i <= deg_-r; ++i) {
	    int k = (i+1) * (i+2) / 2;
	    for (int j = 0; j <= i; ++j) {
		for (int l = 0; l <= j; ++l) {
		    tmp[m] = d[0] * tmp[m]
			+ d[1] * tmp[m + k]
			+ d[2] * tmp[m + 1 + j + k]
			+ d[3] * tmp[m + 2 + j + k];
		    tmp[m] *= deg_ - r + 1;
		    m++;
		}
	    }
	}
    }

    btp.deg_ = deg_ - der;
    int n = (btp.deg_ + 1) * (btp.deg_ + 2) * (btp.deg_ + 3) / 6;
    btp.coefs_.assign(tmp.begin(), tmp.begin()+n);

    return;
}


//===========================================================================
BernsteinPoly
BernsteinTetrahedralPoly::pickLine(const Array<double, 4>& a,
				   const Array<double, 4>&b) const
//===========================================================================
{
    vector<double> coefs(deg_ + 1);
    vector<Array<double, 4> > uvec(deg_, a);
    coefs[0] = blossom(uvec);
    for (int i = 1; i <= deg_; ++i) {
	uvec[deg_-i] = b;
	coefs[i] = blossom(uvec);
    }

    return BernsteinPoly(coefs.begin(), coefs.end());


}


//===========================================================================
BernsteinTetrahedralPoly&
BernsteinTetrahedralPoly::operator*= (const BernsteinTetrahedralPoly& poly)
//===========================================================================
{
    // Aliases for the old coefficients
    const vector<double>& p = coefs_;
    const vector<double>& q = poly.coefs_;
    int psz = (int)coefs_.size();
    int qsz = (int)poly.coefs_.size();

    int d1 = degree();
    int d2 = poly.degree();
    int d = d1 + d2;

    // A vector for the new coefficients
    vector<double> g((d+1)*(d+2)*(d+3)/6);

    for (int i = d; i >= 0; --i) {
	for (int j = d-i; j >= 0; --j) {
	    for (int k = d-i-j; k >= 0; --k) {
		// n is the index og the linear array given (i,j,k,l)
		// with i+j+k+l=d.
		int n = (d-i)*(d-i+1)*(d-i+2)/6
		    + (d-i-j)*(d-i-j+1)/2
		    + d-i-j-k;
		g[n] = 0.0;
		for (int i1 = d1; i1 >= 0; --i1) {
		    for (int j1 = d1-i1; j1 >= 0; --j1) {
			for (int k1 = d1-i1-j1; k1 >= 0; --k1) {
			    int n1 = (d1-i1)*(d1-i1+1)*(d1-i1+2)/6
				+ (d1-i1-j1)*(d1-i1-j1+1)/2
				+ d1-i1-j1-k1;
			    int n2 = (d2-(i-i1))*(d2-(i-i1)+1)*(d2-(i-i1)+2)/6
				+ (d2-(i-i1)-(j-j1))*(d2-(i-i1)-(j-j1)+1)/2
				+ d2-(i-i1)-(j-j1)-(k-k1);
			    if (n1 < 0 || n2 < 0 || n1 >= psz || n2 >= qsz)
				continue;
			    g[n] += p[n1] * q[n2]
				* quadrinomial(d1, i1, j1, k1)
				* quadrinomial(d2, i-i1, j-j1, k-k1)
				/ quadrinomial(d, i, j, k);
			}
		    }
		}
	    }
	}
    }

    deg_ = d;
    swap(coefs_, g);
    return *this;
}


//===========================================================================
BernsteinTetrahedralPoly& BernsteinTetrahedralPoly::operator*= (double c)
//===========================================================================
{
    for (int i = 0; i < int(coefs_.size()); ++i)
	coefs_[i] *= c;
    return *this;
}


//===========================================================================
BernsteinTetrahedralPoly&
BernsteinTetrahedralPoly::operator+= (const BernsteinTetrahedralPoly& poly)
//===========================================================================
{
    THROW("operator+= not yet implemented for BernsteinTetrahedralPoly.");

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
	BernsteinTetrahedralPoly tmp = poly;
	tmp += *this;
	*this = tmp;
    }

    return *this;
}


//===========================================================================
BernsteinTetrahedralPoly& BernsteinTetrahedralPoly::operator+= (double c)
//===========================================================================
{
    *this += BernsteinTetrahedralPoly(c);
    return *this;
}


//===========================================================================
BernsteinTetrahedralPoly&
BernsteinTetrahedralPoly::operator-= (const BernsteinTetrahedralPoly& poly)
//===========================================================================
{
    *this += ((-1.0) * poly);
    return *this;
}


//===========================================================================
BernsteinTetrahedralPoly& BernsteinTetrahedralPoly::operator-= (double c)
//===========================================================================
{
    *this += -c;
    return *this;
}


//===========================================================================
BernsteinTetrahedralPoly& BernsteinTetrahedralPoly::operator/= (double c)
//===========================================================================
{
    for (int i = 0; i < int(coefs_.size()); ++i)
	coefs_[i] /= c;
    return *this;
}


//===========================================================================
void BernsteinTetrahedralPoly::read(istream& is)
//===========================================================================
{
    bool is_good = is.good();
    if (!is_good) {
	THROW("Invalid input file!");
    }

    is >> deg_;
    coefs_.resize((deg_ + 1)*(deg_ + 2)*(deg_ + 3)/6);
    for (int i = 0; i < int(coefs_.size()); ++i) {
	is >> coefs_[i];
    }
}

//===========================================================================
void BernsteinTetrahedralPoly::write(ostream& os) const
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
