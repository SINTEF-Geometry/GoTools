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
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/geometry/SplineDebugUtils.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRProjection.h"
#include "sislP.h"
#include <iostream>
#include <fstream>
#include <numeric>

using namespace Go;
using std::vector;


void Chebyshev(int degree1, int degree2, vector<double>& pol)
{
  int ki, kj, kr;
  vector<vector<double> > Tc(6);
  Tc[0].push_back(1.0);
  Tc[1].insert(Tc[1].end(), {0.0, 1.0});
  Tc[2].insert(Tc[2].end(), {-1.0, 0.0, 2.0});
  Tc[3].insert(Tc[3].end(), {0.0, -3.0, 0.0, 4.0});
  Tc[4].insert(Tc[4].end(), {1.0, 0.0, -8.0, 0.0, 8.0});
  Tc[5].insert(Tc[5].end(), {0.0, 5.0, 0.0, -20.0, 0.0, 16.0});
  int deg1 = std::min(degree1, 6);
  int deg2 = std::min(degree2, 6);
  for (kj=0, kr=0; kj<=deg2; ++kj)
    for (ki=0; ki<=deg1; ++ki)
      {
	pol[kr++] = Tc[deg1][ki]*Tc[deg2][kj];
      }

  int num = 50;
  std::ofstream ofx("Cpolx.g2");
  std::ofstream ofy("Cpoly.g2");
  ofx << "410 1 0 0" << std::endl;
  ofx << num-1 << std::endl;
  double x = -1;
  double x2 = -1;
  double del = 2.0/(double)(num-1);
  double val = Tc[deg1][0];
  for (ki=1; ki<=deg1; ++ki, x2*=x)
    val += Tc[deg1][ki]*x2;
  Point prev(x, val, 0.0);
  for (kj=1, x+=del; kj<num; ++kj, x+=del)
    {
      val = Tc[deg1][0];
      x2 = x;
      for (ki=1; ki<=deg1; ++ki, x2*=x)
	val += Tc[deg1][ki]*x2;
      Point curr(x, val, 0.0);
      ofx << prev << " " << curr << std::endl;
      prev = curr;
    }
  
  ofy << "410 1 0 0" << std::endl;
  ofy << num-1 << std::endl;
  double y = -1;
  double y2 = -1;
  val = Tc[deg2][0];
  for (ki=1; ki<=deg2; ++ki, y2*=y)
    val += Tc[deg2][ki]*y2;
   prev.setValue(y, val, 0.0);
  for (kj=1, y+=del; kj<num; ++kj, y+=del)
    {
      val = Tc[deg2][0];
      y2 = y;
      for (ki=1; ki<=deg2; ++ki, y2*=y)
	val += Tc[deg2][ki]*y2;
      Point curr(y, val, 0.0);
      ofy << prev << " " << curr << std::endl;
      prev = curr;
    }
}

void polynomialTerms(double u, double v, int num1, int num2, int max_num,
		     double *terms)
{
  vector<double> u_terms(num1+1), v_terms(num2+1);
  u_terms[0] = v_terms[0] = 1.0;

  for (int ki=1; ki<=num1; ++ki)
    u_terms[ki] = u_terms[ki-1]*u;
  for (int ki=1; ki<=num2; ++ki)
     v_terms[ki] = v_terms[ki-1]*v;

  for (int kj=0, kr=0; kj<=num2; ++kj)
    for (int ki=0; ki<=num1; ++ki)
      {
	if (ki+kj <= max_num)
	  terms[kr++] = u_terms[ki]*v_terms[kj];
      }
}

shared_ptr<LRSplineSurface> refineSurf(shared_ptr<LRSplineSurface> surf_in, std::istream& is)
{
  shared_ptr<LRSplineSurface> surf(surf_in->clone());
  
  double dom1[4];
  dom1[0] = surf->startparam_u();
  dom1[1] = surf->endparam_u();
  dom1[2] = surf->startparam_v();
  dom1[3] = surf->endparam_v();

  double dom2[4];
  is >> dom2[0] >> dom2[1] >> dom2[2] >> dom2[3];

  double fac1 = (dom1[1] - dom1[0])/(dom2[1] - dom2[0]);
  double fac2 = (dom1[3] - dom1[2])/(dom2[3] - dom2[2]);
  
  int nmb_refs;
  is >> nmb_refs;
  for (int ki=0; ki<nmb_refs; ++ki)
    {
      double parval, start, end;
      int dir;
      int mult;
      is >> parval >> start >> end >> dir >> mult;

      if (dir == 0)
	{
	  parval = dom1[0] + (parval - dom2[0])*fac1;
	  start = dom1[2] + (start - dom2[2])*fac2;
	  end = dom1[2] + (end - dom2[2])*fac2;
	}
      else
	{
	  parval = dom1[2] + (parval - dom2[2])*fac2;
	  start = dom1[0] + (start - dom2[0])*fac1;
	  end = dom1[0] + (end - dom2[0])*fac1;
	}
	
      surf->refine((dir==0) ? XFIXED : YFIXED, parval, start, end, mult, true);
    }

  return surf;
}
				       

int main(int argc, char *argv[])
{
if (argc != 6) {
    std::cout << "Parameters: lr func in (.g2), refinements in, refinements in 2, ref 3, num_samples"  << std::endl;
    exit(-1);
  }

  std::ifstream sf_in(argv[1]);
  std::ifstream ref_in(argv[2]);
  std::ifstream ref_in2(argv[3]);
  std::ifstream ref_in3(argv[4]);
  int nsample = atoi(argv[5]);
  
  ObjectHeader header;
  header.read(sf_in);
  shared_ptr<LRSplineSurface> lrsf(new LRSplineSurface());
  lrsf->read(sf_in);
  if (lrsf->dimension() != 1)
    {
      Point zero(1);
      zero[0] = 0.0;
      for (auto bspl=lrsf->basisFunctionsBegin();
	   bspl!=lrsf->basisFunctionsEnd(); ++bspl)
	{
	  Point cf = bspl->second->Coef();
	  double gamma = bspl->second->gamma();
	  bspl->second->setCoefAndGamma(zero, gamma);
	}
    }
  double range[4];
  range[0] = range[2] = -1;
  range[1] = range[3] = 1;
  lrsf->setParameterDomain(range[0], range[1], range[2], range[3]);
  shared_ptr<LRSplineSurface> lrsf2(lrsf->clone());
  shared_ptr<LRSplineSurface> lrsf3(lrsf->clone());

  double eps = 1.0e-9;
  int dim = 1;
  int degree1 = 4;
  int degree2 = 3;
  int tot_degree = degree1 + degree2;
  double ucoef[5] = {1.0, -7, 35.0/3.0, -7, 1.0};
  //double ucoef[4] = {-1.0, 5.0, -5.0, 1.0};
  double vcoef[4] = {-1.0, 5.0, -5.0, 1.0};
  double coef[20];
  for (int kb=0, kc=0; kb<=degree2; ++kb)
    for (int ka=0; ka<=degree1; ++ka, ++kc)
      coef[kc] = ucoef[ka]*vcoef[kb];

  double uknots[10] = {-1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  //double uknots[8] = {-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0};
  double vknots[8] = {-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0};
  double uknotstp[16] = {-1.0, -1.0, -1.0, -1.0, -1.0, -5.0/7.0, -3.0/7.0,
			 -1.0/7.0, 1.0/7.0, 3.0/7.0, 5.0/7.0,
			 1.0, 1.0, 1.0, 1.0, 1.0};
  double vknotstp[13] = {-1.0, -1.0, -1.0, -1.0, -2.0/3.0, -1.0/3.0, 0.0, 1.0/3.0,
			 2.0/3.0, 1.0, 1.0, 1.0, 1.0};
  double ucoeftp[11] = {1.0, -1.0/7.0, -197.0/147.0, -1019.0/1029.0,
			4099.0/7203.0, 9379.0/7203.0, 4099.0/7203.0, 
			-1019.0/1029.0, -197.0/147.0, -1.0/7.0, 1.0};
  double vcoeftp[9] = {-1.0, 0.0, 10.0/9.0, 1.0, 0.0, -1.0, -10.0/9.0, 0.0, 1.0};
  double coeftp[99];
  for (int kb=0, kc=0; kb<9; ++kb)
    for (int ka=0; ka<11; ++ka, ++kc)
      coeftp[kc] = ucoeftp[ka]*vcoeftp[kb];

  vector<double> nk1 = {-5.0/7.0, -3.0/7.0, -1.0/7.0, 1.0/7.0, 3.0/7.0, 5.0/7.0};
  vector<double> nk2 = {-2.0/3.0, -1.0/3.0, 0.0, 1.0/3.0, 2.0/3.0};
  
  shared_ptr<SplineCurve> ucv(new SplineCurve(degree1+1, degree1+1,
					      &uknots[0], &ucoef[0],
					      dim));
  shared_ptr<SplineCurve> vcv(new SplineCurve(degree2+1, degree2+1,
					      &vknots[0], &vcoef[0],
					      dim));
  shared_ptr<SplineSurface> bsurf(new SplineSurface(degree1+1, degree2+1,
						    degree1+1, degree2+1,
						    &uknots[0], &vknots[0],
						    &coef[0], dim));
  shared_ptr<SplineSurface> bsurftp(new SplineSurface(11, 9,
						    degree1+1, degree2+1,
						    &uknotstp[0], &vknotstp[0],
						    &coeftp[0], dim));
  bsurf->setParameterDomain(range[0], range[1], range[2], range[3]);
  shared_ptr<SplineCurve> ucvtp(new SplineCurve(11, degree1+1,
						&uknotstp[0], &ucoeftp[0],
						dim));
  shared_ptr<SplineCurve> vcvtp(new SplineCurve(9, degree2+1,
						&vknotstp[0], &vcoeftp[0],
						dim));
    std::ofstream ofv1("vcv.g2");
  SplineDebugUtils::writeSpace1DCurve(*vcvtp, ofv1);
    std::ofstream ofu1("ucv.g2");
  SplineDebugUtils::writeSpace1DCurve(*ucvtp, ofu1);

  ucv->insertKnot(nk1);
  vcv->insertKnot(nk2);
  
  shared_ptr<LRSplineSurface> lrsurf(new LRSplineSurface(bsurf.get(), 1.0e-9));
  shared_ptr<LRSplineSurface> lrsurf2 = refineSurf(lrsurf, ref_in);
  shared_ptr<LRSplineSurface> lrsurftp(new LRSplineSurface(bsurftp.get(), 1.0e-9));
  shared_ptr<LRSplineSurface> lrsurftp2 = refineSurf(lrsurftp, ref_in2);
  shared_ptr<LRSplineSurface> lrb2 = refineSurf(lrsurf, ref_in3);
  

  std::ofstream of("ChebBsurf.g2");
  lrsurf->writeStandardHeader(of);
  lrsurf->write(of);
  std::ofstream of1("ChebBsurf2.g2");
  lrsurf2->writeStandardHeader(of1);
  lrsurf2->write(of1);
  std::ofstream oftp("ChebBsurftp.g2");
  lrsurftp->writeStandardHeader(oftp);
  lrsurftp->write(oftp);
   std::ofstream oftp2("ChebBsurftp2.g2");
  lrsurftp2->writeStandardHeader(oftp2);
  lrsurftp2->write(oftp2);
   std::ofstream ofb2("ChebBlrtp.g2");
  lrb2->writeStandardHeader(ofb2);
  lrb2->write(ofb2);
  std::ofstream ofc1("Chebc1.g2");
  SplineDebugUtils::writeSpace1DCurve(*ucv, ofc1);
  std::ofstream ofc2("Chebc2.g2");
  SplineDebugUtils::writeSpace1DCurve(*vcv, ofc2);

  int num_terms = (degree1+1)*(degree2+1);
  vector<double> pol(num_terms, 0.0);
  Chebyshev(degree1, degree2, pol);
  for (int ka=0; ka<num_terms; ++ka)
    std::cout << pol[ka] << " ";
  std::cout << std::endl;
 
  int max_level = 10;
  for (int level=0; level<max_level; ++level)
    {
      int num_update = 0;
      for (auto bspl=lrsf->basisFunctionsBegin();
	   bspl!=lrsf->basisFunctionsEnd(); ++bspl)
	{
	  int blevel = bspl->second->getNestLevel();
	  if (blevel != level)
	    continue;

	  double umin = bspl->second->umin();
	  double umax = bspl->second->umax();
	  double vmin = bspl->second->vmin();
	  double vmax = bspl->second->vmax();
	  shared_ptr<SplineSurface> sub(bsurftp->subSurface(umin, vmin,
							  umax, vmax));

	  BSplineUniLR* uni_u = bspl->second->getUnivariate(XFIXED);
	  BSplineUniLR* uni_v = bspl->second->getUnivariate(YFIXED);
	  int startmult1 = bspl->second->endmult_u(true);
	  int endmult1 = bspl->second->endmult_u(false);
	  int startmult2 = bspl->second->endmult_v(true);
	  int endmult2 = bspl->second->endmult_v(false);
	  int order1 = uni_u->degree() + 1;
	  int order2 = uni_v->degree() + 1;
	  vector<double> kval1 = uni_u->getKnots();
	  int kn1 = order1-startmult1-endmult1+1;
	  vector<double> knots1;
	  vector<double> kval2 = uni_v->getKnots();
	  int kn2 = order2-startmult2-endmult2+1;
	  vector<double> knots2;
	  for (int ka=0; ka<kn1; ++ka)
	    {
	      int kb;
	      for (kb=0; kb<16; ++kb)
		if (fabs(kval1[ka+startmult1]-uknotstp[kb]) < eps)
		  break;
	      if (kb == 16)
		knots1.push_back(kval1[ka+startmult1]);
	    }
	  for (int ka=0; ka<kn2; ++ka)
	    {
	      int kb;
	      for (kb=0; kb<13; ++kb)
		if (fabs(kval2[ka+startmult2]-vknotstp[kb]) < eps)
		  break;
	      if (kb == 13)
		knots2.push_back(kval2[ka+startmult2]);
	    }

	  if (knots1.size() > 0)
	    sub->insertKnot_u(knots1);
	  if (knots2.size() > 0)
	    sub->insertKnot_v(knots2);
	  int k1 = order1 - startmult1;
	  int k2 = order2 - startmult2;
	  int kk = k2*sub->numCoefs_u() + k1;
	  Point coef = Point(sub->coefs_begin()+kk*dim,
			     sub->coefs_begin()+(kk+1)*dim);

	  std::cout << "knots 1: ";
	  for (size_t kh=0; kh<kval1.size(); ++kh)
	    std::cout << kval1[kh] << " ";
	  std::cout << std::endl;
	  std::cout << "knots 2: ";
	  for (size_t kh=0; kh<kval2.size(); ++kh)
	    std::cout << kval2[kh] << " ";
	  std::cout << std::endl;
	  std::cout << "nested depth: " << blevel << std::endl;
	  std::cout << "coefficient: " << coef << std::endl;

	  if (blevel > 0)
	    {
	      //std::cout << "Level: " << blevel << std::endl;
		
	      bool OK = bspl->second->adaptProjCoef(coef);
	      std::cout << "Adjusted coefficient " << coef << std::endl;
	      if (!OK)
		std::cout << "Method1. Adapt coefficient failed" << std::endl;
	    }
	  std::cout << std::endl;
	  lrsf->setCoef(coef, bspl->second.get());
	  num_update++;
	}
      if (num_update == 0)
	break;
    }

  std::ofstream of2("lrproj.g2");
  lrsf->writeStandardHeader(of2);
  lrsf->write(of2);

  double maxdist01 = 0.0, avdist01 = 0.0;
   double udel = 2.0/(double)(nsample-1);
   double vdel = 2.0/(double)(nsample-1);
   double u, v;
   int ki, kj;
   double fac = 1.0/(double)(nsample*nsample);
   for (kj=0, v=-1; kj<nsample; ++kj, v+=vdel)
     {
       v = std::min(v, 1.0);
       for (ki=0, u=-1; ki<nsample; ++ki, u+=udel)
	 {
	   u = std::min(u, 1.0);

	   Point pos1 = bsurftp->ParamSurface::point(u, v);
	   Point pos2 = lrsurftp2->ParamSurface::point(u, v);
	   double dd = pos1.dist(pos2);
	   maxdist01 = std::max(maxdist01, dd);
	   avdist01 += fac*dd;
	 }
     }
   std::cout << "Tp surf ref, max: " << maxdist01 << ", avdist01: " << avdist01 << std::endl;

  double maxdist02 = 0.0, avdist02 = 0.0;
   for (kj=0, v=-1; kj<nsample; ++kj, v+=vdel)
     {
       v = std::min(v, 1.0);
       for (ki=0, u=-1; ki<nsample; ++ki, u+=udel)
	 {
	   u = std::min(u, 1.0);

	   Point pos1 = bsurftp->ParamSurface::point(u, v);
	   Point pos2 = lrsurf2->ParamSurface::point(u, v);
	   double dd = pos1.dist(pos2);
	   maxdist02 = std::max(maxdist02, dd);
	   avdist02 += fac*dd;
	 }
     }
   std::cout << "LRref surf, max: " << maxdist02 << ", avdist02: " << avdist02 << std::endl;

   double maxdist = 0.0, avdist = 0.0;
   // double udel = 2.0/(double)(nsample-1);
   // double vdel = 2.0/(double)(nsample-1);
   // double u, v;
   // int ki, kj;
   //double fac = 1.0/(double)(nsample*nsample);
   for (kj=0, v=-1; kj<nsample; ++kj, v+=vdel)
     {
       v = std::min(v, 1.0);
       for (ki=0, u=-1; ki<nsample; ++ki, u+=udel)
	 {
	   u = std::min(u, 1.0);

	   Point pos1 = bsurftp->ParamSurface::point(u, v);
	   Point pos2 = lrsf->ParamSurface::point(u, v);
	   double dd = pos1.dist(pos2);
	   maxdist = std::max(maxdist, dd);
	   avdist += fac*dd;
	 }
     }

   std::cout << "Compare in sample points" << std::endl;
   std::cout << "Knot insertion. maxdist1: " << maxdist << ", avdist1: " << avdist << std::endl;

   int deg1 = lrsf2->degree(XFIXED);
   int deg2 = lrsf2->degree(YFIXED);
   vector<int> k_ix1(deg1), k_ix2(deg2);
   int db1 = 1, db2 = 1;
   for (int ki=1; ki<=deg1; ++ki)
     {
       db1 *= ki;
       k_ix1[ki-1] = ki;
     }
   for (int kj=1; kj<=deg2; ++kj)
     {
       db2 *= kj;
       k_ix2[kj-1] = kj;
     }
   
   vector<double> tmp(num_terms);
  for (int level=0; level<max_level; ++level)
    {
      int num_update = 0;
      for (auto bspl=lrsf2->basisFunctionsBegin(), bspl2=lrsf->basisFunctionsBegin();
	   bspl!=lrsf2->basisFunctionsEnd(); ++bspl, ++bspl2)
	{
	  int blevel = bspl->second->getNestLevel();
	  if (blevel != level)
	    continue;

	  Point coef =
	    LRProjection::Polynomial2Coef2(pol, degree1, degree2,
					  tot_degree,
					   bspl->second.get(), bsurftp);

	  
	  if (blevel > 0)
	    {
	      //std::cout << "Level: " << blevel << std::endl;
		
	      bool OK = bspl->second->adaptProjCoef(coef);
	      if (!OK)
		std::cout << "Method2. Adapt coefficient failed" << std::endl;
	    }
	  //std::cout << coef << " " << bspl2->second->Coef() << std::endl;
	  lrsf2->setCoef(coef, bspl->second.get());
	  num_update++;
	}
      if (num_update == 0)
	break;
    }

  lrsf2->writeStandardHeader(of2);
  lrsf2->write(of2);

  double maxdist2 = 0.0, avdist2 = 0.0;
   for (kj=0, v=-1; kj<nsample; ++kj, v+=vdel)
     {
       v = std::min(v, 1.0);
       for (ki=0, u=-1; ki<nsample; ++ki, u+=udel)
	 {
	   u = std::min(u, 1.0);

	   polynomialTerms(u, v, degree1, degree2, degree1+degree2, &tmp[0]);
	   double val = 0;
	   for (int kc=0; kc<num_terms; ++kc)
	     val += pol[kc]*tmp[kc];
	   Point pos1 = bsurftp->ParamSurface::point(u, v);
	   Point pos2 = lrsf2->ParamSurface::point(u, v);
	   double dd = pos1.dist(pos2);
	   maxdist2 = std::max(maxdist2, dd);
	   avdist2 += fac*dd;
	 }
     }

   std::cout << "Dual points. maxdist2: " << maxdist2 << ", avdist2: " << avdist2 << std::endl;

   vector<Point> par(num_terms);
   double udel2 = (range[1] - range[0])/(double)degree1;
   double vdel2 = (range[3] - range[2])/(double)degree2;
   int kr;
   for (kj=0, kr=0, v=range[2]; kj<=degree2; ++kj, v+=vdel2)
     {
       v = std::min(v, range[3]);
       for (ki=0, u=range[0]; ki<=degree1; ++ki, ++kr, u+=udel2)
	 {
	   u = std::min(u, range[1]);
	   par[kr] = Point(u, v);
	 }
     }
 
  vector<double> ls(num_terms*num_terms, 0.0), rs(num_terms, 0.0);
  vector<int> piv(num_terms);
  for (int ka=0; ka<num_terms; ++ka)
    piv[ka] = ka;

  for (int ka=0; ka<num_terms; ++ka)
    {
      polynomialTerms(par[ka][0], par[ka][1], degree1, degree2,
		      tot_degree, &ls[ka*num_terms]);
      double val = 0.0;
      for (int kb=0; kb<num_terms; ++kb)
	    val += pol[kb]*ls[ka*num_terms+kb];
      rs[ka] = val;
    }

    // Solve
  int kstat = 0;
  s6lufacp(&ls[0], &piv[0], num_terms, &kstat);
  if (kstat < 0)
    {
      std::cout << "Singular equation system" << std::endl;
      exit(-1);
    }

  s6lusolp(&ls[0], &rs[0], &piv[0], num_terms, &kstat);
  if (kstat < 0)
    {
      std::cout << "Failed solving equation system" << std::endl;
      exit(-1);
    }
   
  for (int level=0; level<max_level; ++level)
    {
      int num_update = 0;
      for (auto bspl=lrsf3->basisFunctionsBegin();
	   bspl!=lrsf3->basisFunctionsEnd(); ++bspl)
	{
	  int blevel = bspl->second->getNestLevel();
	  if (blevel != level)
	    continue;

	  Point coef =
	    LRProjection::Polynomial2Coef(rs, degree1, degree2,
					  tot_degree,
					  range[0], range[1],
					  range[2], range[3],
					  bspl->second.get(), bsurftp);
	  if (blevel > 0)
	    {
	      //std::cout << "Level: " << blevel << std::endl;
	      
	      bool OK = bspl->second->adaptProjCoef(coef);
	    }
	  lrsf3->setCoef(coef, bspl->second.get());
	  num_update++;
	}
      if (num_update == 0)
	break;
    }
  
  lrsf3->writeStandardHeader(of2);
  lrsf3->write(of2);

  double maxdist3 = 0.0, avdist3 = 0.0;
  for (kj=0, v=-1; kj<nsample; ++kj, v+=vdel)
    {
      v = std::min(v, 1.0);
      for (ki=0, u=-1; ki<nsample; ++ki, u+=udel)
	{
	  u = std::min(u, 1.0);

	  polynomialTerms(u, v, degree1, degree2, degree1+degree2, &tmp[0]);
	  double val = 0;
	  for (int kc=0; kc<num_terms; ++kc)
	    val += pol[kc]*tmp[kc];
	  Point pos1 = bsurftp->ParamSurface::point(u, v);
	  Point pos2 = lrsf3->ParamSurface::point(u, v);
	  double dd = pos1.dist(pos2);
	  maxdist3 = std::max(maxdist3, dd);
	  avdist3 += fac*dd;
	}
    }

  std::cout << "Interpolation. maxdist3: " << maxdist3 << ", avdist3: " << avdist3 << std::endl;

  double maxcfdist1=0.0, avcfdist1=0.0, maxcfdist2=0.0, avcfdist2=0.0, maxcfdist3=0.0, avcfdist3=0.0;
  for (auto bspl=lrsurftp2->basisFunctionsBegin(), bspl1=lrsf->basisFunctionsBegin(),
	 bspl2=lrsf2->basisFunctionsBegin(), bspl3=lrsf3->basisFunctionsBegin();
       bspl!=lrsurftp2->basisFunctionsEnd(); ++bspl, ++bspl1, ++bspl2, ++bspl3)
    {
      Point cf = bspl->second->coefTimesGamma();
      Point cf1 = bspl1->second->coefTimesGamma();
      Point cf2 = bspl2->second->coefTimesGamma();
      Point cf3 = bspl3->second->coefTimesGamma();
      double dd1 = cf.dist(cf1);
      double dd2 = cf.dist(cf2);
      double dd3 = cf.dist(cf3);
      maxcfdist1 = std::max(maxcfdist1, dd1);
      avcfdist1 += fac*dd1;
      maxcfdist2 = std::max(maxcfdist2, dd2);
      avcfdist2 += fac*dd2;
      maxcfdist3 = std::max(maxcfdist3, dd3);
      avcfdist3 += fac*dd3;
    }

  std::cout << "Compare coefficients." << std::endl;
  std::cout << "Knot insertion. maxcfdist1: " << maxcfdist1 << ", avcfdist1: " << avcfdist1 << std::endl;
  std::cout << "Dual points. maxcfdist2: " << maxcfdist2 << ", avcfdist2: " << avcfdist2 << std::endl;
  std::cout << "Interpolation. maxcfdist3: " << maxcfdist3 << ", avcfdist3: " << avcfdist3 << std::endl;
 }

