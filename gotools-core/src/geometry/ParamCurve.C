//===========================================================================
//                                                                           
// File: ParamCurve.C                                                      
//                                                                           
// Created: Fri Oct 13 16:36:54 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: ParamCurve.C,v 1.40 2008-11-28 08:12:41 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/utils/GeneralFunctionMinimizer.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/utils/errormacros.h"

//using namespace std;
using std::vector;
using std::min;
using std::max;
using namespace Go;
using std::shared_ptr;

namespace {

// object used to compute the squared distance of a point on a curve and a fixed 
// point in space.
//===========================================================================
class SquaredDistFunctor 
//===========================================================================
{
public:
    SquaredDistFunctor(const Go::ParamCurve* cv, 
		       const Go::Point& pt, 
		       double minpar, 
		       double maxpar)
	: cv_(cv), pt_(pt), minpar_(minpar), maxpar_(maxpar), ptvec_(2) {}
    double operator()(const double* arg) const
    {
	cv_->point(ptvec_[0], *arg);
	return ptvec_[0].dist2(pt_);
    }

    double grad(const double* arg, double* grad) const
    {
	cv_->point(ptvec_, *arg, 1);
	*grad = 2 * (ptvec_[0] - pt_) * ptvec_[1];
	return ptvec_[0].dist2(pt_);
    }

    double minPar(int n) const {return minpar_;}
    double maxPar(int n) const {return maxpar_;}
private:
    const Go::ParamCurve* cv_;
    const Go::Point& pt_;
    double minpar_;
    double maxpar_;
    mutable std::vector<Go::Point> ptvec_;
};

}; // end anonymous namespace


namespace Go
{


//===========================================================================
ParamCurve::~ParamCurve()
//===========================================================================
{
}


//===========================================================================
CompositeBox ParamCurve::compositeBox() const
//===========================================================================
{
    BoundingBox temp = boundingBox();
    return CompositeBox(temp.low(), temp.high());
}


//===========================================================================
Point ParamCurve::point(double tpar) const
//===========================================================================
{
    Point pt(dimension());
    point(pt, tpar);
    return pt;
}


//===========================================================================
std::vector<Point> 
ParamCurve::point(double tpar,
		    int derivs, bool from_right) const
//===========================================================================
{
    std::vector<Point> pts(derivs+1, Point(dimension()));
    point(pts, tpar, derivs, from_right);
    return pts;
}


//===========================================================================
bool ParamCurve::isClosed()
//===========================================================================
{
    // Returns true if the start and end points coincinde within a
    // tolerance. This function may be redefined in subclasses.

    const double space_epsilon = 1.0e-8;

    double startpar = startparam();
    double endpar = endparam();
    Point startpt, endpt;
    point(startpt, startpar);
    point(endpt, endpar);

    return (startpt - endpt).length() < space_epsilon;
}


//===========================================================================
void ParamCurve::uniformEvaluator(int num, vector<Point>& points, vector<double>& param) const
//===========================================================================
{
  param.resize(num);
  points.resize(num);

  double current_par = startparam();
  double step;
  if (num == 1)   // Only to avoid zero division
    step = 0.0;
  else
    step = (endparam() - current_par) / (double)(num - 1);

  for (int i = 0; i < num; ++i, current_par += step)
    {
      param[i] = current_par;
      point(points[i], current_par);
    }

}



//===========================================================================
double ParamCurve::estimatedCurveLength(int numpts) const
//===========================================================================
{
  Point pprev = point(startparam());
  Point pnext;
  double length = 0;
  const double tmp = double(1) / double(numpts - 1);
  for (int i = 1; i < numpts; ++i) {
      double fac = double(i) * tmp;
      pnext = point((1.0-fac)*startparam() + fac*endparam());
      length += pnext.dist(pprev);
      pprev = pnext;
  }
  return length;
}

//===========================================================================
double ParamCurve::estimatedCurveLength(double tmin, double tmax, int numpts) const
//===========================================================================
{
    const double pareps = 1.0e-10;
    bool outside_domain = (tmin < startparam() - pareps ||
			   tmax > endparam() + pareps);
    ALWAYS_ERROR_IF(outside_domain,
		    "Input interval not in curve domain.");

    Point pprev = point(tmin);
    Point pnext;
    double length = 0;
    const double tmp = double(1) / double(numpts - 1);
    for (int i = 1; i < numpts; ++i) {
	double fac = double(i) * tmp;
	pnext = point((1.0-fac)*tmin + fac*tmax);
	length += pnext.dist(pprev);
	pprev = pnext;
    }
    return length;
}



//===========================================================================
void ParamCurve::closestPoint(const Point&   pt,
				double&   clo_t,
				Point&  clo_pt,
				double&   clo_dist) const
//===========================================================================
{
    // Calls the virtual closestPoint function
    closestPoint(pt,
		 startparam(),
		 endparam(),
		 clo_t,
		 clo_pt,
		 clo_dist);
}


//===========================================================================
double ParamCurve::nextSegmentVal(double par, bool forward, double tol) const
//===========================================================================
{
  if (forward)
    return endparam();
  else
    return startparam();
}


//===========================================================================
void ParamCurve::closestPointGeneric(const Point&   pt,
				       double    tmin,
				       double    tmax,
				       double guess_param,
				       double&   clo_t,
				       Point&  clo_pt,
				       double&   clo_dist) const
//===========================================================================
{
    const double TOL = 1.0e-8;
    static bool use_s1771 = true;

    // Make sure that the start parameter lies within the legal interval
    guess_param = std::max(guess_param, tmin);
    guess_param = std::min(guess_param, tmax);

    if (use_s1771)
    {
	int kstat = 0;
	s1771(pt, TOL, tmin, tmax, guess_param, clo_t, &kstat);
	clo_pt = this->point(clo_t);
	clo_dist = pt.dist(clo_pt);
    }
    else
    {
    SquaredDistFunctor distfun(this, pt, tmin, tmax);
    FunctionMinimizer<SquaredDistFunctor> funmin(1, distfun, &guess_param, TOL);

    Point dir(1);
    distfun.grad(&guess_param, &(dir[0])); // determining direction

    bool hit_domain_edge = false;
    funmin.minimize(dir, hit_domain_edge);

    clo_t = funmin.getPar(0);
    clo_pt = this->point(clo_t);
    clo_dist = sqrt(funmin.fval());
    }
}


//===========================================================================
double ParamCurve::length(double tol, double tstart, double tend)
//===========================================================================
{
  vector <vector<double> > len;
  int termination_depth = 15;
  vector<Point> pts(2);

  point(pts[0], tstart);
  point(pts[1], tend);

  len.resize(1);
  len[0].resize(1);
  len[0][0] = pts[0].dist(pts[1]);

  for (int i = 1; i < termination_depth; ++i)
    {
      // Expand the array of points
	int oldn = (int)pts.size();
      vector<Go::Point> newpts(2*oldn - 1);
      double piecelen = (tend - tstart)/double(2*oldn - 2);
      for(int j = 0; j < oldn-1; ++j) {
	// Copy the old point
	newpts[2*j] = pts[j];
	// Compute the new point
	double t = tstart + (2*j + 1)*piecelen;
	point(newpts[2*j+1],t);
      }
      // Copy the last point
      newpts[2*oldn - 2] = pts[oldn - 1];
      // Swap new point vector to output vector
      pts.swap(newpts);

      // Expand the array and compute the leftmost array of the new row.
      len.resize(i+1);
      len[i].resize(i+1);
      //      len[i][0] = compositeChordLength(pts);
      vector<double> cdist;
      cdist.resize(pts.size()-1);
      for (int j = 0; j < int(pts.size()) - 1; ++j)
	cdist[j] = pts[j].dist(pts[j+1]);
      for (int j = ((int)pts.size() - 1) >> 1; j > 0; j >>= 1)
	for (int k = 0; k < j; ++k)
	  cdist[k] = cdist[k<<1] + cdist[(k<<1) | 1];
      len[i][0] = cdist[0];

      // Compute rest of Romberg table
      double fac = 0.0;
      for (int j = 1; j < i+1; ++j) {
	fac *= 4.0;
	fac += 3.0;
	len[i][j] = len[i][j-1] + (len[i][j-1] - len[i-1][j-1])/fac;
      }

      double err = fabs(len[i][i] - len[i-1][i-1]);
      if (err < tol * len[i][i]) return len[i][i];

    }

  return len[termination_depth-1][termination_depth-1];

}

//===========================================================================
//
// Default implementation. Should be specialized in sub classes if it 
// is expected to be used
  vector<shared_ptr<ParamCurve> > 
ParamCurve::split(double param, double fuzzy) const
//===========================================================================
{
  vector<shared_ptr<ParamCurve> > cvs(2);
  cvs[0] = shared_ptr<ParamCurve>(subCurve(startparam(), param, fuzzy));
  cvs[1] = shared_ptr<ParamCurve>(subCurve(param, endparam(), fuzzy));
  return cvs;
}



} // namespace Go

// Below is commented-out code which belongs to an earlier implementation
// of the closest point calculation.  It is no longer used, but kept here
// for future reference.

// // an internally used functor.  Doxygen documentation unnecessary.
// class dist_deriv_functor
// {
// public:
//     dist_deriv_functor(const ParamCurve* cv, const Point& pt)
// 	: cv_(cv), pt_(pt), x_(2,Point(cv_->dimension())) {}
//   ~dist_deriv_functor()
//   {}
//     double operator() (double t)
// 	{
// 	    cv_->point(x_, t, 1);
// 	    x_[0] -= pt_;
// 	    return x_[0]*x_[1];
// 	}
// private:
//     const ParamCurve* cv_;
//     const Point& pt_;
//     std::vector<Point> x_;
// };

// // an internally used functor.  Doxygen documentation unnecessary.
// class dist_functor
// {
// public:
//     dist_functor(const ParamCurve* cv, const Point& pt)
// 	: cv_(cv), pt_(pt) {}
//     ~dist_functor()
//   {}
//     double operator() (double t) const
// 	{
// 	    return (cv_->point(t) - pt_).length2(); // No need to extract root.
// 	}
// private:
//     const ParamCurve* cv_;
//     const Point& pt_;
// };

/*
//-----------------------------------------------------------------------------
void ParamCurve::closestPointGeneric(const Point&   pt,
				       double    tmin,
 				       double    tmax,
 				       double guess_param,
 				       double&   clo_t,
 				       Point&  clo_pt,
 				       double&   clo_dist) const
//-----------------------------------------------------------------------------
{
    // Solving the equation min{f(t)}, with f being the dist_functor:
    // f(t) = (x(t)-pt)*(x(t)-p(t)) (inner product, hence dist^2).
    dist_functor f(this, pt);

    double ax = guess_param;
    double step = (tmax - tmin)*1e-05;
    double bx = (guess_param - step < tmin)
        ? guess_param + step : guess_param - step;
    double cx, fa, fb, fc;
    bracket_minimum(f, ax, bx, cx, fa, fb, fc); // We try to bracket solution.
    // Function may have iterated outside boundaries.
    ax = (ax < tmin ? tmin : ax);
    ax = (tmax < ax ? tmax : ax);
    bx = (bx < tmin ? tmin : bx);
    bx = (tmax < bx ? tmax : bx);
    cx = (cx < tmin ? tmin : cx);
    cx = (tmax < cx ? tmax : cx);
//      ax = max(tmin, min(tmax, ax));
//      bx = max(tmin, min(tmax, bx));
//      cx = max(tmin, min(tmax, cx));
    double tolerance = 1e-14;
    if ((fabs(cx-ax) < tolerance)) { //tmin-bx) < tolerance) || (fabs(bx-tmax) < tolerance)) {
// 	// Minimum is an end point, we return that value as the closest point.
	// We have already reached the solution within satisfactory tolerance.
	clo_t = 0.5*(ax + cx); //bx;
    } else { // Solution is bracketed, we locate the minimum.
	if (bx == ax || bx == cx) {
	    double dx = 0.5*(ax + cx);
	    double fd = f(dx);
	    fa = f(ax);
	    fc = f(cx);
	    if (fd < min(fa, fc)) {
		bx = dx;
	    }
	}
	dist_deriv_functor df(this, pt);
	try {
	    derivative_brent_minimize(&f, &df, ax, bx, cx, clo_t, tolerance);
	} catch (const IterationCountExceeded& ob) {
	    MESSAGE("Error in solving equation. Iterationcount exceeded. "
		       "Returning best guess.");
	    clo_t = ob.guess_par;
	} catch (...) {
	    GO_ERROR("Unexpected exception occured!", UnknownError());
	}
    }

    clo_pt = point(clo_t);
    clo_dist = (clo_pt - pt).length();
}
*/

void
ParamCurve::s1771(Point pt,double aepsge,
		  double astart,double aend,double anext,double &cpos,int *jstat) const
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Newton iteration on the distance function between
*              a curve and a point, to find a closest point or an
*              intersection point.
*              If a bad choice for the guess parameter is given in, the
*              iteration may end at a local, not global closest point.
*
*
* INPUT      : ppoint  - The point in the closest point problem.
*              pcurve  - The curve in the closest point problem.
*              aepsge  - Geometrical resolution.
*              astart  - Curve parameter giving the start of the search
*                        interval.
*              aend    - Curve parameter giving the end of the search
*                        interval.
*              anext   - Curve guess parameter for the closest point
*                        iteration.
*
*
*
* OUTPUT     : cpos    - Resulting curve parameter from the iteration.
*              jstat   - status messages
*                                = 2   : A minimum distanse found.
*                                = 1   : Intersection found.
*                                < 0   : error.
*
*
* METHOD     : Newton iteration.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, Feb 1989
* Revised by : Christophe Rene Birkeland, SINTEF Oslo, May 1993.
*              Error message corrected
*
*********************************************************************
*/
{
    int kstat = 0;
  int kdim;              /* Dimension of space the curves lie in             */
  double tdelta;         /* Parameter interval of the curve.                 */
  double tdist;          /* Distance between position and origo.             */
  double td;             /* Distances between old and new parameter value in
			    the two parameter directions.                    */
  double tprev;          /* Previous difference between the curves.          */
  vector<Point> val(3);     /* Value ,first and second derivatie on curve 1    */
  Point diff;         /* Difference between the curves                    */
  int quick = (*jstat);  /* Indicates if the exactness requirement is
                            relaxed.                                         */
  int max_it = 20;       /* Maximum number of iterations.                    */
  vector<Point> pnt(3);

  if (quick) max_it = 10;

  /* Test input.  */

  kdim = dimension();

  /* Fetch endpoints and the intervals of parameter interval of curves.  */

  tdelta = endparam() - startparam();

  /* Initiate variable.  */

  tprev = (double)1.0e10;

  /* Evaluate 0-2.st derivatives of the curve. */

  val = point(anext, 2); 

  diff = pt - val[0];

  tprev = tdist = pt.dist(val[0]);

  td = s1771_s9del(diff.begin(), val[1].begin(), val[2].begin(), kdim);

  /* Correct if we are not inside the parameter intervall. */

  if (anext + td < astart) td = astart - anext;
  else if (anext + td > aend) td = aend - anext;

  s1771_s9point(pt, val,diff,astart,aend,max_it,&anext,
	  &td,tdelta,&tdist,tprev, &kstat);

  /* Iteration stopped, test if point found is within resolution */

  if (tdist <= aepsge)
    *jstat = 1;
  else
    *jstat = 2;

  /* Iteration completed.  */
  cpos = anext;

}

void
ParamCurve::s1771_s9point(Point pt, vector<Point> val, Point diff,
			  double astart,double aend,int max_it,double *cnext,double *ad,
			  double adel,double *cdist,double aprev,int *jstat) const
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Newton iteration on the distance function between
*              a point and a curve to find a closest point or an
*              intersection point.
*
*
* INPUT      : pcurve  - Pointer to the curve.
*              eval1[] - Array containing the coefisients to the point.
*              eval2[] - Array containing the coefisients to a start
*                        iteration point on the curve.
*              ediff[] - The vector between the point and the start
*                        point on the curve.
*              astart  - Start parameter value on the curve.
*              aend    - End parameter value on the curve.
*              max_it  - Maximum number of iterations.
*              ad      - A first computed parametric step.
*              adel    - The differrence between start and end
*                        parameter value.
*              aprev    -A prevous distance between the curve and
*                        the point.
*              ileft   - An estimate of the knot number whice is
*                        closest to the point on the curve.
*
*
* OUTPUT     : cnext   - Parameter value of curve in the closest point.
*              cdist    -The current distance between the curve and
*                        the point.
*              jstat   - status messages
*                                >= 0   : OK.
*                                <  0   : error.
*
*
* METHOD     :
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, May 1989
* MODIFIED BY : Vibeke Skytt, SINTEF SI, 10.93. Stop earlier when divergence.
*
*********************************************************************
*/
{
  // int kstat = 0;            /* Local status variable.                      */
  int kdim = dimension();  /* Dimension of space the curves lie in        */
  int knbit;                /* Number of iterations                        */
  int kdiv = 0;             /* Counts number of diverging steps.           */
  int kdiv2 = 0;            /* Counts number of almost divergence.         */
  double REL_COMP_RES = 1.0e-15;

  /* Correct if we are not inside the parameter intervall. */

  if (*cnext + *ad < astart)    *ad = astart - *cnext;
  else if (*cnext + *ad > aend) *ad = aend - *cnext;

  for (knbit=0;knbit<max_it;knbit++)
    {
      /* Evaluate 0-2.st derivatives of the curve. */

	val = point(*cnext + *ad, 2); 

      diff = pt - val[0];

      *cdist = pt.dist(val[0]);

      if (*cdist -aprev  <= REL_COMP_RES)
	{
	   if (kdiv2 > 4) break;
	   if (*cdist -aprev >= 0.0) kdiv2++;

	   kdiv = 0;
	  aprev = *cdist;
	  *cnext += *ad;

	  *ad = s1771_s9del(diff.begin(), val[1].begin(), val[2].begin(), kdim);

	  /* Correct if we are not inside the parameter intervall. */

	  if (*cnext + *ad < astart)    *ad = astart - *cnext;
	  else if (*cnext + *ad > aend) *ad = aend - *cnext;
	}
      else
	{
	   kdiv++;
	   if (kdiv > 3) break;
	  (*ad) /= (double)2;
	  /** knbit--; */
	}
      if (fabs((*ad)/std::max(fabs(*cnext),adel)) <= REL_COMP_RES) break;
    }

}

double
ParamCurve::s1771_s9del(double *eco,double *eco1,double *eco2,int idim) const
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To compute the distance on the parameter line to a point
*            on the curve on which the tangent is orthogonal to the
*            difference vector from this point on the curve to the
*            point in the space.
*
*
* INPUT      : eco   - The differens vector beetween the point and the
*                      current posision on the curve.
*              eco1  - The first derevative vector in the  current posision
*                      on the curve.
*              eco2  - The second derevative vector in the  current posision
*                      on the curve.
*              idim  - Dimension of space the vectors lie in.
*
*
* OUTPUT     : s1771_s9del - The computed parameter distance.
*
*
* METHOD     : We have to find the parameter distance "dt" from
*              the equation:
*                <ecoef-dt*ecoef1,ecoef1+dt*ecoef2> = 0.
*              This may be written:
*                  <ecoef,ecoef1> + <ecoef,ecoef2>*dt
*                - <ecoef1,ecoef1>*dt + <ecoef1,ecoef2>*dt*dt = 0
*              The following function is the solution of this second
*              degree equation. We are always using the solution
*              with the smallest absolute value.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, Mar 1989
*
*********************************************************************
*/
{
  double t1,t2,t3,t4,t5,t6;   /* Constants in equation.                    */
  double tmax,tmax1;          /* Max values in equation.                   */
  double ttol=(double)1e-10;  /* Relative tolerance in equation.           */
  double REL_PAR_RES = 1.0e-12;

  t1 =  inner(eco, eco+idim, eco1);
  t3 =  inner(eco1, eco1+idim, eco1);
  t2 =  t3 - inner(eco, eco+idim, eco2);
  t4 =  -(double)2 * inner(eco1, eco1+idim, eco2);

  tmax  = max(fabs(t1),fabs(t2));
  tmax1 = max(fabs(t3),fabs(t4));
  tmax  = max(tmax1,tmax);

  if (fabs(tmax) < REL_PAR_RES)                    return 0.0;

  else if (fabs(t4)/tmax < ttol) /* The second degree part is degenerated. */
    {
      if (fabs(t2)/tmax < ttol)
	{
          if (fabs(t3)/tmax < ttol)        return 0.0;
          else                             return (t1/t3);
	}
      else                                  return (t1/t2);
    }
  else  /* An ordinary second degree equation.    */
    {
      t5 = t2*t2 - (double)2*t4*t1;
      if (t5 < 0.0)                       return (t1/t3);
      else
	{
          t6 = sqrt(t5);
          t5 = (t2 + t6)/t4;
          t6 = (t2 - t6)/t4;
	  t1 *= t3;


          /* We have two solutions and we want to use the one
	     with the same sign as we get while using an other
	     metode t1/t3. If both solutions have the same
	     sign we use the one with smallest value. */

          if (t1 < 0.0)
	    {
	      if (t5 <= 0.0 && t6 <= 0.0)
		{
		  if (t5 > t6)             return t5;
	          else                     return t6;
		}
	      else if (t5 <= 0.0)        return t5;
	      else if (t6 <= 0.0)        return t6;
              else                         return min(t5,t6);
	    }
	  else if (t1 > 0.0)
	    {
	      if (t5 >= 0.0 && t6 >= 0.0)
		{
		  if (t5 < t6)             return t5;
	          else                     return t6;
		}
	      else if (t5 >= 0.0)        return t5;
	      else if (t6 >= 0.0)        return t6;
              else                         return max(t5,t6);
	    }
	  else                             return min(fabs(t5),fabs(t6));
	}
    }
}
