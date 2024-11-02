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

#include "GoTools/utils/GeneralFunctionMinimizer.h"
#include <vector>
using std::vector;
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/ClosestPoint.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/utils/Values.h"          // MAXDOUBLE

//***************************************************************************
//
// Implementation file of the free function ClosestPoint::closestPtCurves in namespace
// ClosestPoints, defined in ClosestPoint.h/
//
//***************************************************************************

using namespace Go;

// Anonymous namespace for definition of class CrvDistFun
namespace
{
// Distance function between two curves.  Used by the minimization algorithm
// initiated by ClosestPoint::closestPtCurves.
class CrvDistFun {
public:
    CrvDistFun(const ParamCurve* cv1, 
	       const ParamCurve* cv2,
	       const double* const minpar = 0,
	       const double* const maxpar = 0);
    
    inline double operator()(const double* arg) const;
    inline double grad(const double* arg, double* res) const;
    inline double minPar(int pardir) const;
    inline double maxPar(int pardir) const;

private:
    double minpar_[4];
    double maxpar_[4];
    const ParamCurve * const cv1_;
    const ParamCurve * const cv2_;
    mutable Point p1_, p2_, d_;
    mutable vector<Point> p1vec_, p2vec_;
};
} // End of anonymous namespace for declaration of class CrvDistFun


// Anonymous namespace for helper functions declaration.
namespace
{
void computeSeedCvCv(const SplineCurve* pc1, const SplineCurve* pc2,
		     double& seed1, double& seed2);

void computeSeedCvCv2(const ParamCurve* cv1, const ParamCurve* cv2,
		     double& seed1, double& seed2);

void insideParamDomain(double& delta, double acoef, double astart,
		       double aend);

void nextStep(double& cdist, double& cdiff1, double& cdiff2,
	      std::vector<Point>& eval1, std::vector<Point>& eval2);

} // End of anonymous namespace for helper functions declaration.


namespace Go {

//===========================================================================
void ClosestPoint::closestPtCurves(const ParamCurve* cv1, const ParamCurve* cv2,
				   double& par1, double& par2, double& dist,
				   Point& ptc1, Point& ptc2, int max_passes)
//===========================================================================

// Compute the closest point between two curves.
{

  DEBUG_ERROR_IF(cv1->dimension()!=cv2->dimension(), "Dimension mismatch.");

  double seed1, seed2;
  double tmin1,tmax1,tmin2,tmax2;

  // Use all of the parameter domain.
  tmin1 = cv1->startparam();
  tmax1 = cv1->endparam();
  tmin2 = cv2->startparam();
  tmax2 = cv2->endparam();
 
  if (cv1->instanceType() == Class_SplineCurve &&
      cv2->instanceType() == Class_SplineCurve) {
      const SplineCurve *pc1 = dynamic_cast<const SplineCurve*>(cv1);
      const SplineCurve *pc2 = dynamic_cast<const SplineCurve*>(cv2);

    // Compute seed values
    computeSeedCvCv(pc1, pc2, seed1, seed2);
  }
  else {
    computeSeedCvCv2(cv1, cv2, seed1, seed2);
    // seed1 = 0.5*(tmin1+tmax1);
    // seed2 = 0.5*(tmin2+tmax2);
  }

  // Iterate for closest point
  ClosestPoint::closestPtCurves(cv1,cv2,tmin1,tmax1,tmin2,tmax2,seed1,seed2,par1,par2,
				dist,ptc1,ptc2, max_passes);

}

// (s1770)
//===========================================================================
void  ClosestPoint::closestPtCurves(const ParamCurve* cv1, const ParamCurve* cv2, double tmin1,
 		     double tmax1, double tmin2, double tmax2,
 		     double seed1, double seed2, double& par1, double& par2,
				    double& dist, Point& ptc1, Point& ptc2,
				    int max_passes)
//===========================================================================
{

    DEBUG_ERROR_IF(cv1->dimension()!=cv2->dimension(), "Dimension mismatch.");

    const double REL_COMP_RES = 0.000000000000001;
    double anext1 = seed1; // Estimated start values
    double anext2 = seed2;

    double tdelta1 = std::max(1.0, cv1->endparam() - cv1->startparam());
    double tdelta2 = std::max(1.0, cv2->endparam() - cv2->startparam());

    double tprev = MAXDOUBLE;

    // Evaluate 0-1.st derivatives of both curves.
    int nder=1;    // Order of derivatives to be calulated
    std::vector<Point> sval1(nder+1), sval2(nder+1);
    cv1->point(sval1, anext1, nder);
    cv2->point(sval2, anext2, nder);
    tprev = sval1[0].dist(sval2[0]);

    // Compute the distance vector and value and the new step.
    double tdist, cdiff1, cdiff2;
    double td[2],t1[2],tdn[2];  // Distances between old and new parameter-
    // value in the two parameter directions.  
    nextStep(tdist, cdiff1, cdiff2, sval1, sval2);
    td[0] = cdiff1;
    td[1] = cdiff2; 

    // Adjust if we are not inside the parameter interval.
    t1[0] = td[0];
    t1[1] = td[1];
    insideParamDomain(t1[0], anext1, tmin1, tmax1);
    insideParamDomain(t1[1], anext2, tmin2, tmax2);


    // Iterate for closest point
    int kdir;                  // Changing direction. 
    //  int stat = 0;
  
    for ( int knbit = 0; knbit < max_passes; knbit++) {
    
	// Evaluate 0-1.st derivatives of both curves
 
	cv1->point(sval1, anext1+t1[0], nder);
	cv2->point(sval2, anext2+t1[1], nder);      
 
	// Compute the distance vector and value and the new step.
	nextStep(tdist, cdiff1, cdiff2, sval1, sval2);
	tdn[0] = cdiff1;
	tdn[1] = cdiff2; 
        
	// Check if the direction of the step have changed.
      
	kdir = (td[0]*tdn[0]+td[1]*tdn[1] >= 0.0);     // 0 if changed.
      
	// Ordinary converging.
      
	if (tdist < tprev*0.9 || (kdir && tdist < tprev)) 
	  {
	    anext1 += t1[0];
	    anext2 += t1[1];
      
	    td[0] = tdn[0];
	    td[1] = tdn[1];
      
	    // Correct if we are not inside the parameter intervall.
      
	    t1[0] = td[0];
	    t1[1] = td[1];
	    insideParamDomain(t1[0], anext1, tmin1, tmax1);
	    insideParamDomain(t1[1], anext2, tmin2, tmax2);
      
	    if ( (fabs(t1[0]/tdelta1) <= REL_COMP_RES) &&
		 (fabs(t1[1]/tdelta2) <= REL_COMP_RES) ) 
	      break;
      
	    tprev = tdist;
	}
    
	// Not converging, adjust and try again.
    
	else {
	    if ( (fabs(t1[0]/tdelta1) <= REL_COMP_RES) &&
		 (fabs(t1[1]/tdelta2) <= REL_COMP_RES) ) 
	      break;
      
	    if (tprev+tdist < REL_COMP_RES)
	      break;

	    t1[0] = tprev*t1[0]/(tprev+tdist);
	    t1[1] = tprev*t1[1]/(tprev+tdist);

	}
    }
  
    par1 = anext1;
    par2 = anext2;

    // Compute the points and the distance between them.
    cv1->point(ptc1, par1);
    cv2->point(ptc2, par2);
    dist = ptc1.dist(ptc2);

}
} // End of namespace Go


namespace {  // Anonymous namespace for class CrvDistFun.
//===========================================================================
CrvDistFun::CrvDistFun(const ParamCurve* cv1, 
		       const ParamCurve* cv2,
		       const double* const minpar,
		       const double* const maxpar)
//===========================================================================
    : cv1_(cv1), cv2_(cv2), p1vec_(2), p2vec_(2)
{
    if (!minpar) {
	minpar_[0] = cv1_->startparam();
	minpar_[1] = cv2_->startparam();
    } else {
	minpar_[0] = minpar[0];
	minpar_[1] = minpar[1];
    }
    if (!maxpar) {
	maxpar_[0] = cv1_->endparam();
	maxpar_[1] = cv2_->endparam();
    } else {
	maxpar_[0] = maxpar[0];
	maxpar_[1] = maxpar[1];
    }
}

//===========================================================================    
double CrvDistFun::operator()(const double* arg) const
//===========================================================================
{
    cv1_->point(p1_, arg[0]);
    cv2_->point(p2_, arg[1]);
    return p1_.dist2(p2_);
}

//===========================================================================
double CrvDistFun::grad(const double* arg, double* res) const
//===========================================================================
{
    cv1_->point(p1vec_, arg[0], 1);
    cv2_->point(p2vec_, arg[1], 1);
    d_ = p1vec_[0] - p2vec_[0];
    
    res[0] = 2 * d_ * p1vec_[1];
    res[1] = -2 * d_ * p2vec_[1];
    
    return d_.length2();
}

//===========================================================================
double CrvDistFun::minPar(int pardir) const
//===========================================================================
{
    ASSERT(pardir == 0 || pardir == 1);
    return minpar_[pardir];
}

//===========================================================================
double CrvDistFun:: maxPar(int pardir) const
//===========================================================================
{
    ASSERT(pardir == 0 || pardir == 1);
    return maxpar_[pardir];
}
}  // End of anonymous namespace for class CrvDistFun.


// Anonymous namespace for helper functions definition.
namespace {
/** Computes initial start points for iteration along the curves.
   * \param pc1 Curve number one.
   * \param pc2 Curve number two.
   * \retval seed1 Start point for iteration along curve number one.
   * \retval seed2 Start point for iteration along curve number two.
   */
//***************************************************************************
void computeSeedCvCv(const SplineCurve* cv1, const SplineCurve* cv2,
		     double& seed1, double& seed2)
//***************************************************************************
{

  // Make guess point to the iteration.
  // Find position of closest vertices
  std::vector<double>::const_iterator co1 = cv1->coefs_begin();
  std::vector<double>::const_iterator co2 = cv2->coefs_begin();
  std::vector<double>::const_iterator co3;
  std::vector<double>::const_iterator co12 = cv1->coefs_end();
  std::vector<double>::const_iterator co22 = cv2->coefs_end();

  const int dim = cv1->dimension();
  DEBUG_ERROR_IF(dim!=cv2->dimension(), "Dimension mismatch.");
  double td, tmin=1.0e8;
  int minidx1=0, minidx2=0;
  int ki, k1, k2;
  for (k1=0; co1<co12; co1+=dim, k1++) {
    for (k2=0, co3=co2; co3<co22; co3+=dim, k2++) {
      for (td=0.0, ki=0; ki<dim; ki++)
	td += (co1[ki]-co3[ki])*(co1[ki]-co3[ki]);
      if (td < tmin) {
	tmin = td;
	minidx1 = k1;
	minidx2 = k2;
      }
    }
  }

  // Estimate parameter value of vertices
  std::vector<double>::const_iterator st;
  int kk = cv1->order();
  for (k1=minidx1+1, st=cv1->basis().begin(), seed1=0.0;
       k1<minidx1+kk; seed1+=st[k1], k1++);
  seed1 /= (double)(kk-1);
  kk = cv2->order();
  for (k1=minidx2+1, st=cv2->basis().begin(), seed2=0.0;
       k1<minidx2+kk; seed2+=st[k1], k1++);
  seed2 /= (double)(kk-1);

}

//***************************************************************************
void computeSeedCvCv2(const ParamCurve* cv1, const ParamCurve* cv2,
		      double& seed1, double& seed2)
//***************************************************************************
{
  const int dim = cv1->dimension();
  DEBUG_ERROR_IF(dim!=cv2->dimension(), "Dimension mismatch.");

  // Make guess point to the iteration.
  // Find position of closest vertices
  std::vector<double>::const_iterator co1;
  std::vector<double>::const_iterator co2;
  std::vector<double>::const_iterator co3;
  std::vector<double>::const_iterator co12;
  std::vector<double>::const_iterator co22;

  const SplineCurve *pc1 = dynamic_cast<const SplineCurve*>(cv1);
  const CurveOnSurface *sfc1 = dynamic_cast<const CurveOnSurface*>(cv1);
  if (sfc1)
    pc1 = dynamic_cast<const SplineCurve*>(sfc1->spaceCurve().get());
  const SplineCurve *pc2 = dynamic_cast<const SplineCurve*>(cv2);
  const CurveOnSurface *sfc2 = dynamic_cast<const CurveOnSurface*>(cv2);
  if (sfc2)
    pc2 = dynamic_cast<const SplineCurve*>(sfc2->spaceCurve().get());

  int num_sample = 5;
  vector<double> pts1;
  vector<double> pts2;
  vector<double> par1(num_sample);
  vector<double> par2(num_sample);
  if (pc1 && pc1->numCoefs() > 3)
    {
      co1 = pc1->coefs_begin();
      co12 = pc1->coefs_end();
    }
  else
    {
      double t1 = cv1->startparam();
      double t2 = cv1->endparam();
      double tdel = (t2 - t1)/(double)(num_sample+1);
      double tpar;
      int ka;
      for (ka=0, tpar=t1+0.5*tdel; ka<num_sample; ++ka, tpar+=tdel)
	{
	  Point pt = cv1->point(tpar);
	  pts1.insert(pts1.end(), pt.begin(), pt.end());
	  par1[ka] = tpar;
	}
      co1 = pts1.begin();
      co12 = pts1.end();
     }
  
  if (pc2 && pc2->numCoefs() > 3)
    {
      co1 = pc2->coefs_begin();
      co22 = pc2->coefs_end();
    }
  else
    {
      double t1 = cv2->startparam();
      double t2 = cv2->endparam();
      double tdel = (t2 - t1)/(double)(num_sample+1);
      double tpar;
      int ka;
      for (ka=0, tpar=t1+0.5*tdel; ka<num_sample; ++ka, tpar+=tdel)
	{
	  Point pt = cv2->point(tpar);
	  pts2.insert(pts2.end(), pt.begin(), pt.end());
	  par2[ka] = tpar;
	}
      co2 = pts2.begin();
      co22 = pts2.end();
   }
  
  double td, tmin=1.0e8;
  int minidx1=0, minidx2=0;
  int ki, k1, k2;
  for (k1=0; co1<co12; co1+=dim, k1++) {
    for (k2=0, co3=co2; co3<co22; co3+=dim, k2++) {
      for (td=0.0, ki=0; ki<dim; ki++)
	td += (co1[ki]-co3[ki])*(co1[ki]-co3[ki]);
      if (td < tmin) {
	tmin = td;
	minidx1 = k1;
	minidx2 = k2;
      }
    }
  }

  // Estimate parameter value of vertices
  if (pc1 && pc1->numCoefs() > 3)
    {
      std::vector<double>::const_iterator st;
      int kk = pc1->order();
      for (k1=minidx1+1, st=pc1->basis().begin(), seed1=0.0;
	   k1<minidx1+kk; seed1+=st[k1], k1++);
      seed1 /= (double)(kk-1);
    }
  else
    seed1 = par1[minidx1];
  
  if (pc2 && pc2->numCoefs() > 3)
    {
      std::vector<double>::const_iterator st;
      int kk = pc2->order();
      for (k1=minidx2+1, st=pc2->basis().begin(), seed2=0.0;
	   k1<minidx2+kk; seed2+=st[k1], k1++);
      seed2 /= (double)(kk-1);
    }
  else
    seed2 = par2[minidx2];

}

/** Adjust delta to satisfy: \f[ astart \leq acoef+delta \leq aend \f]
 * Ported from the sisl function s1770_s9corr.
 */
//***************************************************************************
// (s1770_s9corr)
void insideParamDomain(double& delta, double acoef, double astart,
		       double aend)
//***************************************************************************
{
  // Make sure that the corrected parameters still lies in the domain.
  //  astart <= acoef+delta <= aend

  if (acoef + delta < astart)
    delta = astart - acoef;
  else if (acoef + delta > aend)
    delta = aend - acoef;
}

/** Computes the distance vector and value beetween
 * a point on the first curve and a point on the second
 * curve. And computes a next step on both curves.
 * This is equivalent to the nearest way to the
 * parameter plane in the tangent plane from a point in the
 * distance surface between two curves.
 * Ported from the sisl function s1770_s9dir.
 * METHOD : The method is to compute the parameter distance to the points
 * on both tangents which is closest to each other.
 * The difference vector beetween these points are orthogonal
 * to both tangents. If the distance vector beetween the two
 * points on the curve is "diff" and the two derivative vectors
 * are "der1" and "der2", and the two wanted parameter distances
 * are "dt1" and "dt2", then we get the following system of 
 * equations:
 * @verbatim
 <dt1*der1+dist-dt2*der2,der2> = 0
 <dt1*der1+dist-dt2*der2,der1> = 0
 This is further:
 
 | -<der1,der2>   <der2,der2> |  | dt1 |   | <diff,der2> |
 |                            |  |     | = |             |
 | -<der1,der1>   <der1,der2> |  | dt2 |   | <diff,der1> |
 @endverbatim
 * 
 * The solution of this matrix equation dt1,dt2 are returned in the
 * parameters cdiff1,cdiff2.
 *
 */
//***************************************************************************
// (s1770_s9dir)
void nextStep(double& cdist, double& cdiff1, double& cdiff2,
	      std::vector<Point>& eval1, std::vector<Point>& eval2)
//***************************************************************************
{
  const double TOL = 1.0e-12;

  Point& p1 = eval1[0];  // Value
  Point& d1 = eval1[1];  // 1. derivative

  Point& p2 = eval2[0];
  Point& d2 = eval2[1]; 

  Point gdiff = p1 - p2;  // Distance vector
  cdist = gdiff.length(); // Length of distance vector

  double t1,t2,t3,t4,t5;   // Variables in equation system
  // scalar products
  t1 = d1*d1;
  t2 = d1*d2;
  t3 = d2*d2;
  t4 = gdiff*d1;
  t5 = gdiff*d2;

  double tdet = t2*t2 - t1*t3;  // Determinant

  //  double delta_t1, delta_t2;
  if (fabs(tdet) < TOL) {
    cdiff1 = 0.0;
    cdiff2 = 0.0;
  }
  else {   // Using Cramer's rule to find the solution of the system
    cdiff1 =  (t4*t3 - t5*t2)/tdet;
    cdiff2 =  (t2*t4 - t1*t5)/tdet;
  }
}

} // End of anonymous namespace for helper functions definition.

