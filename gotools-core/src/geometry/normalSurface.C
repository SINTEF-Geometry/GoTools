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
#include "GoTools/utils/binom.h"
#include "GoTools/geometry/GeometryTools.h"

using namespace std;

namespace Go {



// Anonymous namespace
namespace {

void outerprod1(double vec1[], double vec2[], double prod[])
{
  prod[0] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
  prod[1] = vec1[0]*vec2[2] - vec1[2]*vec2[0];
  prod[2] = vec1[0]*vec2[3] - vec1[3]*vec2[0];
  prod[3] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
  prod[4] = vec1[1]*vec2[3] - vec1[3]*vec2[1];
  prod[5] = vec1[2]*vec2[3] - vec1[3]*vec2[2];
}

void outerprod2(double vec1[], double vec2[], double prod[])
{
  prod[0] =  vec1[5]*vec2[1] - vec1[4]*vec2[2] + vec1[3]*vec2[3];
  prod[1] = -vec1[5]*vec2[0] + vec1[2]*vec2[2] - vec1[1]*vec2[3];
  prod[2] =  vec1[4]*vec2[0] - vec1[2]*vec2[1] + vec1[0]*vec2[3];
}

void outerprod3(double vec1[], double vec2[], double prod[])
{
  prod[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
  prod[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
  prod[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
}


} // Anonymous namespace

//==========================================================================
SplineSurface* SplineSurface::normalSurface() const
//==========================================================================
{

  const SplineSurface* tsurf;       // Intermediate surfaces      
  const SplineSurface* psurf;

  vector<SplineSurface> bezpatch;   // Bezier patches as SplineSurfaces
  vector<SplineSurface> bezpatchdu;
  vector<SplineSurface> bezpatchdv;

  const int dim = dimension();   // Dimension of geometry space.
  int hdim;                      // Dimension of homogenous space.
  const int n = order_u();       // Order of input surface in 1. par. dir
  const int m = order_v();       // Order of input surface in 2. par. dir

  double vec3[3], vec6[6];       // Help vectors in cross product computations.

  // Ensure k-regularity of surface.
  if (basis_u().isKreg() && basis_v().isKreg())
    tsurf=this;
  else
    tsurf=subSurface(startparam_u(),startparam_v(),endparam_u(),endparam_v());


  if (rational()) {
    // NURBS surface. Make surface in homogenuous space.
    hdim = dim+1;
    psurf = new SplineSurface(numCoefs_u(),numCoefs_v(),order_u(),order_v(),
			      tsurf->basis_u().begin(), tsurf->basis_v().begin(),
			      tsurf->rcoefs_begin(), dim+1, false);

    // Convert the original surface to Bezier patches.
    GeometryTools::splitSurfaceIntoPatches(*psurf, bezpatch);
  }
  else {
    hdim = dim;
    psurf = tsurf;
  }


  // Differentiate the surface in both parameter directions.
  SplineSurface* dp_du = psurf->derivSurface(1,0); 
  SplineSurface* dp_dv = psurf->derivSurface(0,1);

  // Convert the derivative surfaces to Bezier patches.
  GeometryTools::splitSurfaceIntoPatches(*dp_du, bezpatchdu);
  GeometryTools::splitSurfaceIntoPatches(*dp_dv, bezpatchdv);

  // Calculate number of Bezier patches.
  int numpat_u = psurf->numberOfPatches_u();
  int numpat_v = psurf->numberOfPatches_v();

  // Calculate the order of normal surface.
  int order1, order2;
  if (rational()) {
    order1 = 3*order_u() - 4;
    order2 = 3*order_v() - 4;
  }
  else {
    order1 = 2*order_u() - 2;
    order2 = 2*order_v() - 2;
  }


  // Calculate number of coefficients of normal surface.
  int numc_u = numpat_u*order1;
  int numc_v = numpat_v*order2;


  // Knot vectors of normal surface.
  vector<double> knots1(numc_u+order1);
  vector<double> knots2(numc_v+order2);

  // Coefficients (vertices) of normal surface.
  vector<double> coefs(dim*numc_u*numc_v);
  
  // scratch vectors for coefficients of scaled basis functions.
  vector<double> bp_d,bp_f,bp_n;
  if (rational()) {
    bp_d.resize((2*n-1)*m*12);
    bp_f.resize(order2*dim);
    bp_n.resize(dim*order1*(order2+1));
  }
  else
    bp_n.resize(dim*order1*order2);    

  double bint, bfac;

  // Loop through all the Bezier patches and compute the normal surface patch
  // for each Bezier patch.

  int pu, pv, p, vb, kj;
  int i, j, k, l, d;
  int ka, kb, kc, kd,ke,kn,kr;
  int k1,k2,k3,k4,k5,k6,k7;
  double tstart_u, tend_u, tstart_v, tend_v;
  
  for (pv=0; pv<numpat_v; ++pv) {
    vb = pv*numpat_u;
    kr = pv*order2*numc_u*dim;

    for (pu=0; pu<numpat_u; ++pu,kr+=order1*dim) {
      p = vb + pu;
      
      tstart_u = bezpatchdu[p].startparam_u();
      tend_u   = bezpatchdu[p].endparam_u();

      tstart_v = bezpatchdv[p].startparam_v();
      tend_v   = bezpatchdv[p].endparam_v();

      // Store knots
      if (pu == 0) {
	if (pv == 0) {
	  for (kj=0; kj<order2; ++kj)
	       knots2[kj] = tstart_v;
	}
	for (kj=0; kj<order2; ++kj)
	  knots2[(pv+1)*order2+kj] = tend_v;
      }

      if (pv == 0) {
	if (pu == 0) {
	  for (kj=0; kj<order1; ++kj)
	    knots1[kj] = tstart_u;
	}
	for (kj=0; kj<order1; ++kj)
	  knots1[(pu+1)*order1+kj] = tend_u;
      }


      vector<double>& bp_a = bezpatchdu[p].activeCoefs();
      vector<double>& bp_b = bezpatchdv[p].activeCoefs(); 
 
      if (!rational()) {

	// A nonrational surface. 
	// Compute representations of Pu and Pv in a scaled basis.

	for (j=0,ka=0; j<m; ++j) {
	  bint=binom(m-1,j);
	  for (i=0; i<n-1; ++i,ka+=dim) {
	    bfac=binom(n-2,i)*bint;
	    if (bfac != 1.0) {
	      for (d=0; d<dim; ++d)
		bp_a[ka+d] *= bfac;
	    }
	  }
	}

	for (l=0,kb=0; l<m-1; ++l) {
	  bint=binom(m-2,l);
	  for (k=0; k<n; ++k,kb+=dim) {
	    bfac=binom(n-1,k)*bint; 
	    if (bfac != 1.0) {
	      for (d=0; d<dim; ++d)
		bp_b[kb+d] *= bfac;
	    }
	  }
	}

	// Multiplication
	for (kn=0; kn<int(bp_n.size()); ++kn)
	  bp_n[kn] = 0.0;

	for (j=0,ka=0; j<m; ++j)
	  for (i=0; i<n-1; ++i,ka+=dim)
	    for (l=0,kb=0; l<m-1; ++l) {
	      kn=(j+l)*order1;
	      for (k=0; k<n; ++k,kb+=dim) {
		outerprod3(&bp_a[ka], &bp_b[kb], vec3);
		kn=((j+l)*order1+i+k)*dim;
		for (d=0; d<dim; d++)
		  bp_n[kn+d] += vec3[d];
	      }
	    } //l

	// Conversion of the result to the Bernstein basis.
	for (j=0,kn=0; j<order2; ++j) {
	  bint = binom(m+m-3,j);
	  for (i=0; i<order1; ++i,kn+=dim) {
	    bfac = 1.0/(binom(n+n-3,i)*bint);
	    if (bfac != 1.0) {
	      for (d=0; d<dim; d++)
		bp_n[kn+d] *= bfac;
	    }
	  }
	}
      }

      else {
	// A rational surface.
	// Compute representations of P, Pu and Pv in a scaled basis

	const int hh = (3*m)/2 - 2;
	const int kkk1=2*n-1;
	const int kkk2=2*m-1;
	vector<double>& bp_c = bezpatch[p].activeCoefs();

	// Conversion to the scaled bases.
	for (j=0,ka=0,kc=0; j<m; ++j,kc+=hdim) {
	  bint=binom(m-1,j);
	  for (i=0; i<n-1; ++i,ka+=hdim,kc+=hdim) {
	    bfac=binom(n-2,i)*bint;
	    if (bfac != 1.0)
	      for (d=0; d<hdim; ++d) {
		bp_a[ka+d] *= bfac;
		bp_c[kc+d] *= bfac;
	      }
	  }
	}

	for (j=0,kb=0; j<m-1; ++j) {
	  bint=binom(m-2,j);
	  for (i=0; i<n; ++i,kb+=hdim) {
	    bfac=binom(n-1,i)*bint;
	    if (bfac != 1.0)
	      for (d=0; d<hdim; ++d)
		bp_b[kb+d] *= bfac;
	  }
	}
	  
	//                    ^
	// Multiplication of  P and Pu (page 706)
	for (kd=0; kd<int(bp_d.size()); ++kd)
	  bp_d[kd] = 0.0;

	for (j=0,kc=0; j<m; ++j,kc+=hdim) {
	  for (i=0; i<n-1; ++i,kc+=hdim)
	    for (l=0,ka=0; l<m;++l) {
	      k1=(j+l)*(2*n-2);
	      for(k=0; k<n-1; ++k,ka+=hdim) {
		outerprod1(&bp_c[kc], &bp_a[ka], vec6);
		kd=(k1+i+k)*6;
		for (d=0; d<6; ++d)
		  bp_d[kd+d] += vec6[d];
	      }
	    }
	}

	//                    ^
	// Multiplication of (P ^ Pu) and Pv (page 707)
	for (ke=0; ke<int(bp_n.size()); ++ke)
	  bp_n[ke] = 0.0;

	for (j=0,k1=0; j<kkk2; ++j,k1+=kkk1-1) {
	  k7 = k1*6;
	  for (l=0,k3=0; l<m-1;++l,k3+=n) {
	    if (j+l != hh) {
	      k5 = (j+l)*order1;
	      for (i=0, k2=k7; i<kkk1-2; ++i, k2+=6)
		for (k=0, k4=k3*hdim, k6=(k5+i)*dim; k<n;
		     ++k, k4+=hdim, k6+=dim) {
		  outerprod2(&bp_d[k2], &bp_b[k4], vec3);
		  for (d=0; d<dim; ++d)
		    bp_n[k6+d] += vec3[d];
		}
	    }
	  }
	}

	// Degree reduction and convertion to Bezier basis.

	k4 = order1*dim;
	k5 = (order2-1)*dim;
	k6 = order1*order2*dim;
	for (i=0, k1=0; i<order1; ++i, k1+=dim) {
	  bint=binom(order1-1,i);

	  // Degree reduction
	  for (d=0; d<dim; ++d)
	    bp_f[d] = bp_n[k1+d];
	  
	  for (j=1, k2=dim, k7=k4; j<hh; ++j, k2+=dim, k7+=k4)
	    for (d=0; d<dim; d++)
	      bp_f[k2+d] = bp_n[k7+k1+d] - bp_f[k2-dim+d];
	  
	  for (d=0; d<dim; d++)
	    bp_f[k5+d] = bp_n[k6+k1+d];
	  
	  for (j=order2-2, k2=j*dim, k7=k6-k4; j>=hh; --j, k2-=dim, k7-=k4)
	    for (d=0; d<dim; ++d)
	      bp_f[k2+d] = bp_n[k7+k1+d] - bp_f[k2+dim+d];

	  // Convertion to Bezier basis.
	  for (j=0, k3=0, k7=0; j<order2; ++j, k3+=k4, k7+=dim) {
	    bfac = 1.0/(bint*binom(order2-1,j));
	    for (d=0; d<dim; ++d)
	      bp_n[k3+k1+d] = bp_f[k7+d]*bfac;
	  }
	}

      }

      // Add the contribution from the current patch of the
      // normal surface to the total normal surface.

      k1=numc_u*dim;
      for (j=0, k4=0,k2=0; j<order2; ++j,k2+=k1) {
	for (i=0,k3=0; i<order1; ++i,k4+=dim,k3+=dim)
	  for (d=0; d<dim; ++d)
	    coefs[kr+k2+k3+d] = bp_n[k4+d];
      }
      
    }
  }
	  

  // Create surface object containing the normal surface

  SplineSurface* result;
  result = new SplineSurface(numc_u, numc_v, order1, order2,
			     knots1.begin(), knots2.begin(),
			     coefs.begin(), dim, psurf->rational());
  
// Free memory.
  if (tsurf != this)
    delete tsurf;

  if (rational())
    delete psurf;

  delete dp_du;
  delete dp_dv;


  return result;
}

} // namespace Go
