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

#include "GoTools/parametrization/PrMatSparse.h"
#include "GoTools/utils/errormacros.h"
#include <cmath>

using namespace std;

//-----------------------------------------------------------------------------
PrMatSparse::PrMatSparse(int m, int n, int num_nonzero)
//-----------------------------------------------------------------------------
    : m_(m), n_(n), p_(num_nonzero),
      irow_(m+1), jcol_(num_nonzero), a_(num_nonzero)
{
  irow_[m_] = num_nonzero;
}

PrMatSparse::PrMatSparse(int m, int n, int num_nonzero,
			 const int* irow, const int* jcol, const double* data)
    : m_(m), n_(n), p_(num_nonzero),
      irow_(irow, irow + m + 1), jcol_(jcol, jcol + num_nonzero), 
      a_(data, data + num_nonzero)
{
    ALWAYS_ERROR_IF(irow_[m_] != num_nonzero, "Invalid data given for sparse matrix.");
}

//-----------------------------------------------------------------------------
PrMatSparse::~PrMatSparse()
//-----------------------------------------------------------------------------
{
//     cerr << "PrMatSparse::~PrMatSparse()" << endl;
}


//-----------------------------------------------------------------------------
void PrMatSparse::redim(int m, int n, int num_nonzero)
//-----------------------------------------------------------------------------
{
  m_ = m;
  n_ = n;
  p_ = num_nonzero;
  irow_.resize(m+1);
  jcol_.resize(num_nonzero);
  a_.resize(num_nonzero);
  irow_[m_] = num_nonzero;
}


//-----------------------------------------------------------------------------
void PrMatSparse::setToMatrix(const PrMatrix& m, double tol)
//-----------------------------------------------------------------------------
{
    int R = m.rows();
    int C = m.colmns();
    vector<int> irow, jcol;
    vector<double> data;
    for (int r = 0; r < R; ++r) {
	irow.push_back((int)data.size());
	for (int c = 0; c < C; ++c) {
	    double x = m(r,c);
	    if (fabs(x) <= tol) {
	    } else {
		data.push_back(x);
		jcol.push_back(c);
	    }
	}
    }
    irow.push_back((int)data.size());
    ASSERT(int(irow.size()) == R + 1);
    ASSERT(int(data.size()) == irow[R]);

    m_ = R;
    n_ = C;
    p_ = (int)data.size();
    irow_.swap(irow);
    jcol_.swap(jcol);
    a_.swap(data);
}


//-----------------------------------------------------------------------------
void PrMatSparse::prod(const PrVec& x, PrVec& y) const // Find y = Ax
//-----------------------------------------------------------------------------
{
  if(x.size() != n_ || y.size() != m_)
  {
    MESSAGE("Error in PrMatSparse::prod");
    MESSAGE("Matrix and vectors have incompatible sizes");
    return;
  }

  int i,k;
  for(i=0; i<m_; i++)
  {
    y(i) = 0.0;
    for(k=irow(i); k<irow(i+1); k++)
    {
      y(i) += (*this)(k) * x(jcol(k));
    }
  }
}

//-----------------------------------------------------------------------------
void PrMatSparse::print(std::ostream& os)
//-----------------------------------------------------------------------------
{
  int k;
  os << "irow = " << std::endl;
  for(k=0; k<int(irow_.size()); k++)
  {
    os << irow(k) << std::endl;
  }
  os << std::endl;

  os << "jcol = " << std::endl;
  for(k=0; k<p_; k++)
  {
    os << jcol(k) << std::endl;
  }
  os << std::endl;

  os << "a = " << std::endl;
  for(k=0; k<p_; k++)
  {
    os << (*this)(k) << std::endl;
  }
  os << std::endl;
}
//-----------------------------------------------------------------------------
void PrMatSparse::printFull(std::ostream& os)
//-----------------------------------------------------------------------------
{
    // Do one row at a time
    vector<double> rowcontent(n_);
    for (int row = 0; row < m_; ++row) {
	std::fill(rowcontent.begin(), rowcontent.end(), 0.0);
	for (int i = irow_[row]; i < irow_[row+1]; ++i) {
	    rowcontent[jcol_[i]] = a_[i];
	}
	for (int col = 0; col < n_; ++col) {
	    os << rowcontent[col] << "  ";
	}
	os << '\n';
    }
}



//-----------------------------------------------------------------------------
void PrMatSparse::read(std::istream& is)
//-----------------------------------------------------------------------------
{
    int ind = 0;
    double d;
    for (int i = 0; i < m_; ++i) {
	irow_[i] = ind;
	for (int j = 0; j < n_; ++j) {
	    is >> d;
	    if (d != 0.0) {
		ALWAYS_ERROR_IF(ind > p_, "Too many nonzero elements!");

		a_[ind] = d;
		jcol_[ind] = j;
		++ind;
	    }
	}
    }
}


//-----------------------------------------------------------------------------
int PrMatSparse::rows() const 
//-----------------------------------------------------------------------------
{
    return m_;
}


//-----------------------------------------------------------------------------
int PrMatSparse::colmns() const
//-----------------------------------------------------------------------------
{
    return n_;
}


//-----------------------------------------------------------------------------
double PrMatSparse::operator () (int i, int j) const
//-----------------------------------------------------------------------------
{

  for(int k=irow(i); k<irow(i+1); k++)
  {
    if(jcol(k) == j) return (*this)(k);
  }
  return 0.0;
}


//-----------------------------------------------------------------------------
void PrMatSparse::matProd(PrMatSparse& B, PrMatSparse& C) const
//-----------------------------------------------------------------------------
{
  // this routine multiplies two sparse matrices, "A=(this)" times "B" and
  // stores the result in "C"

  // first check: do the dimensions of "A" and "B" match?
  if(B.rows() != n_)
  {
    MESSAGE("Error in PrMatSparse::matProd");
    MESSAGE("Matrices have incompatible sizes");
    return;
  }

  int i;
  int numNonZeros = 0;
  vector< vector<int> > k_idx (m_);
  vector< vector<double> > C_ik (m_);

  // compute all C_ik's
  for (i=0; i<m_; i++) {                        // for all rows i of A

    for (int j=irow(i); j<irow(i+1); j++) {         // for all columns jj in
      int jj = jcol(j);                             // row i of A

      for (int k=B.irow(jj); k<B.irow(jj+1); k++) { // for all columns kk in
	int kk = B.jcol(k);                         // row jj of B
	
	// is kk already in k_idx[i] ? and if so, where ?
	int idx = -1;
	for (size_t l=0; l<k_idx[i].size(); l++) {
	  if (k_idx[i][l] == kk) {
	      idx = (int)l;
	    break;
	  }
	}

	// if kk is a "new" column index, append k_idx[i] with it
	if (idx < 0) {
	  k_idx[i].push_back(kk);
	  C_ik[i].push_back(0.0);
	  idx = (int)C_ik[i].size()-1;
	}

	// accumulate the dot products for all C(i,k)'s
	C_ik[i][idx] += a_[j] * B(k);
      }
    }
    numNonZeros += (int)k_idx[i].size();
  }

  // transfer the result to the PrMatSparse structure
  C = PrMatSparse(m_, B.colmns(), numNonZeros);
  int offset = 0;

  for (i=0; i<m_; i++) {
    C.irow(i) = offset;

    for (size_t k=0; k<k_idx[i].size(); k++) {
      C(offset) = C_ik[i][k];
      C.jcol(offset) = k_idx[i][k];
      offset++;
    }
  }
}

