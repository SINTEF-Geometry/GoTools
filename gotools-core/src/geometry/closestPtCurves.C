#include "GoTools/utils/GeneralFunctionMinimizer.h"
#include <vector>
using std::vector;
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/closestPtCurves.h"
#include "GoTools/utils/Values.h"          // MAXDOUBLE

//***************************************************************************
//
// Implementation file of the free function closestPtCurves defined in
// closestPtCurves.h/
//
//***************************************************************************

using namespace Go;

namespace { // anonymous namespace 

// distance function between two curves.  Used by the minimization algorithm
// initiated by closestPtCurves.
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


} // end anonymous namespace

namespace Go {


//===========================================================================
void closestPtCurves(const ParamCurve* cv1, const ParamCurve* cv2,
		     double& par1, double& par2, double& dist,
		     Point& ptc1, Point& ptc2)
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
    seed1 = 0.5*(tmin1+tmax1);
    seed2 = 0.5*(tmin2+tmax2);
  }

  // Iterate for closest point
  closestPtCurves(cv1,cv2,tmin1,tmax1,tmin2,tmax2,seed1,seed2,par1,par2,
		  dist,ptc1,ptc2);

}

// //===========================================================================
// void closestPtCurves(const ParamCurve* cv1,  // curve number one
// 		     const ParamCurve* cv2,  // curve number two
// 		     double umin,           // min. parameter value for cv1
// 		     double umax,           // max. parameter value for cv1
// 		     double vmin,           // min. parameter value for cv2
// 		     double vmax,           // max. parameter value for cv2
// 		     double u_seed,         // start point for iter. along cv1
// 		     double v_seed,         // start point for iter. along cv2
// 		     double& par1,           // param. val. of found closest pt. in cv1
// 		     double& par2,           // param. val. of found closest pt. in cv2
// 		     double& dist,           // distance between found closest points
// 		     Point& ptc1,            // found closest pt. in cv1
// 		     Point& ptc2)            // found closest pt. in cv2
// //===========================================================================
// {
//     const double TOL = 1.0e-8;
//     double seed[2];
//     seed[0] = u_seed;
//     seed[1] = v_seed;
//     double minpar[2];
//     minpar[0] = umin;
//     minpar[1] = vmin;
//     double maxpar[2];
//     maxpar[0] = umax;
//     maxpar[1] = vmax;

//     // establing distance function to minimize
//     CrvDistFun distfun(cv1, cv2, minpar, maxpar);

//     // minimize the distance function
//     FunctionMinimizer<CrvDistFun> funmin(2, distfun, seed, TOL);
//     minimise_conjugated_gradient(funmin);//, 3); // number of iterations in each cycle

//     // calculate and copy results
//     par1 = funmin.getPar(0);
//     par2 = funmin.getPar(1);
//     dist = sqrt(funmin.fval());
//     ptc1 = cv1->point(par1);
//     ptc2 = cv2->point(par2);
// }

// //===========================================================================
// void closestPtCurves(const ParamCurve* cv1,  // curve number one
// 		     const ParamCurve* cv2,  // curve number two
// 		     double umin,           // min. parameter value for cv1
// 		     double umax,           // max. parameter value for cv1
// 		     double vmin,           // min. parameter value for cv2
// 		     double vmax,           // max. parameter value for cv2
// 		     double u_seed,         // start point for iter. along cv1
// 		     double v_seed,         // start point for iter. along cv2
// 		     double& par1,           // param. val. of found closest pt. in cv1
// 		     double& par2,           // param. val. of found closest pt. in cv2
// 		     double& dist,           // distance between found closest points
// 		     Point& ptc1,            // found closest pt. in cv1
// 		     Point& ptc2)            // found closest pt. in cv2
// //===========================================================================
// {
//     // defining distance function
//     const double TOL = 1.0e-8;
//     const double EPS = 1.0e-10;
//     DistanceFunctionMinimizer dfmin(cv1, cv2, umin, umax, vmin, vmax, u_seed, v_seed, TOL);

//     // minimizing distance function using conjugated gradients
//     double gradient[2], old_gradient[2];
//     dfmin.grad(old_gradient);
//     double dir[2];
//     dir[0] = -old_gradient[0];
//     dir[1] = -old_gradient[1]; // using negative gradient as first step direction
//     double dir_norm_2 = (dir[0] * dir[0]) + (dir[1] * dir[1]);
//     int i;

//     int num_minimizations = 0; // @@ debug purposes
//     double old_val = dfmin.fval();

//     while (true) { 
	
// 	// make sure direction is not uphill (is this already guaranteed??), 
// 	// and truncating if at border of domain
// 	if (dir[0] * old_gradient[0] + dir[1] * old_gradient[1] > 0) {
// 	    dir[0] *= -1;
// 	    dir[1] *= -1;
// 	}
// 	if ((dfmin.atUmin() && dir[0] < 0 ) || (dfmin.atUmax() && dir[0] > 0)) {
// 	    dir[0] = 0;
// 	}
// 	if ((dfmin.atVmin() && dir[1] < 0 ) || (dfmin.atVmax() && dir[1] > 0)) {
// 	    dir[1] = 0;
// 	}

// 	dir_norm_2 = dir[0] * dir[0] + dir[1] * dir[1];
	
// 	if (dir_norm_2 < EPS) {
// 	    // we believe we have reached a minimum
// 	    break;
// 	}
	
// 	// minimize along this direction
// 	num_minimizations++;
// 	bool hit_domain_edge = false; // 'dir' has not (yet) been truncated
// 	const bool allow_premature_end = true; // let the minimization algorithm end
// 	                                       // prematurely if it decides the direction 
// 	                                       // is not optimal enough.
// 	double new_val = dfmin.minimize(dir, hit_domain_edge, allow_premature_end); 

// 	if (2.0 * fabs(new_val - old_val) <= TOL * (fabs(new_val) + fabs(old_val) + EPS)) {
// 	    // we have reached a minimum
// 	    break;
// 	} else {
// 	    old_val = new_val;
// 	}

// 	// choose new direction using conjugated gradients (Polak-Ribiere variant)
// 	dfmin.grad(gradient);
// 	double factor = 0;
// 	double old_grad_norm_2 = 0;
// 	if (!hit_domain_edge) {
// 	    // we reached a non-border minimum on the last iteration, which makes it 
// 	    // worthwhile to seek a conjugate direction.  We must calculate a nonzero 
// 	    // factor
// 	    for (i = 0; i < 2; ++i) {
// 		factor += gradient[i] * (gradient[i] - old_gradient[i]);
// 		old_grad_norm_2 += old_gradient[i] * old_gradient[i];
// 	    }
// 	    factor /= old_grad_norm_2;
// 	}
// 	for (i = 0; i < 2; ++i) {
// 	    old_gradient[i] = gradient[i];
// 	    dir[i] = dir[i] * factor - gradient[i];
// 	}



//     }    
//     // copying result variables
//     par1 = dfmin.curU();
//     par2 = dfmin.curV();
//     dist = sqrt(dfmin.fval());
//     dfmin.points(ptc1, ptc2);

//     //std::cout << "Number of directions tried: " << num_minimizations << std::endl; // @@ debug purposes
// }



// (s1770)
//===========================================================================
void closestPtCurves(const ParamCurve* cv1, const ParamCurve* cv2, double tmin1,
 		     double tmax1, double tmin2, double tmax2,
 		     double seed1, double seed2, double& par1, double& par2,
 		     double& dist, Point& ptc1, Point& ptc2)
//===========================================================================
{

    DEBUG_ERROR_IF(cv1->dimension()!=cv2->dimension(), "Dimension mismatch.");

    const double REL_COMP_RES = 0.000000000000001;
    double anext1 = seed1; // Estimated start values
    double anext2 = seed2;

    double tdelta1 = cv1->endparam() - cv1->startparam();
    double tdelta2 = cv2->endparam() - cv2->startparam();

    double tprev = MAXDOUBLE;

    // Evaluate 0-1.st derivatives of both curves.
    int nder=1;    // Order of derivatives to be calulated
    std::vector<Point> sval1(nder+1), sval2(nder+1);
    cv1->point(sval1, anext1, nder);
    cv2->point(sval2, anext2, nder);

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
    const int max_passes = 30;
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
      
	if (tdist < tprev*0.9 || kdir) {
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
		 (fabs(t1[1]/tdelta2) <= REL_COMP_RES) ) break;
      
	    tprev = tdist;
	}
    
	// Not converging, adjust and try again.
    
	else {
	    //      t1[0] *= 0.5;
	    //      t1[1] *= 0.5; 
	    t1[0] = tprev*t1[0]/(tprev+tdist);
	    t1[1] = tprev*t1[1]/(tprev+tdist);
	    /* knbit--; */
	}
    }
  
    par1 = anext1;
    par2 = anext2;

    // Compute the points and the distance between them.
    cv1->point(ptc1, par1);
    cv2->point(ptc2, par2);
    dist = ptc1.dist(ptc2);

}

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

} // namespace Go  

namespace {

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

};
