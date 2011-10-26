//===========================================================================
//                                                                           
// File: new_closestPtSurfSurfPlane.C                                        
//                                                                           
// Created: Mon Apr 25 12:45:04 2005                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: closestPtSurfSurfPlane_functional.C,v 1.7 2005-08-15 07:52:46 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include <vector>
#include <limits>
#include "GoTools/geometry/closestPtSurfSurfPlane.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/utils/Point.h"
#include "GoTools/utils/GeneralFunctionMinimizer.h"

using namespace Go;
using namespace std;

namespace {
double compute_constraint_multiplier(const Point& plane_normal, // plane normal
				     const vector<Point>& point1, // 0-1 derivs in p1
				     const vector<Point>& point2, // 0-1 derivs in p2
				     const double aepsge); // tolerance

double compute_distance_to_plane(const Point& p1, 
				 const Point& p2, 
				 const Point& plane_point,
				 const Point& plane_normal);

// double compute_distance_to_plane(const double* p, 
// 				 const Point& plane_point, 
// 				 const Point& plane_normal);


// function object class for passing the function f(x) + c * P(x) to the 
// DistanceFunctionMinimizer class.  ( f(x) is the distance between the two points
// on the two surfaces, and P(x) is the distance from the points to the defined plane).
class SurfDistFun {
public:
    SurfDistFun(const ParamSurface* psurf1, 
		const ParamSurface* psurf2,
		double C,
		const Point& plane_point,
		const Point& plane_normal);
    
    inline double operator()(const double* arg) const;
    inline void grad(const double* arg, double* res) const;
    inline double minPar(int pardir) const;
    inline double maxPar(int pardir) const;

    static void compute_distance_grad(const vector<Point>& p1,
				      const vector<Point>& p2,
				      double* result);

private:
    const ParamSurface * const ps1_;
    const ParamSurface * const ps2_;
    const double C_;
    const Point& pl_normal_;
    const double plane_affinity_;
    double min_par_[4];
    double max_par_[4];

    mutable Point pt1_, pt2_; // scratch
    mutable vector<Point> ptvec1_, ptvec2_; // scratch

    inline double signedDistFromPlane(const Point& p) const;
};

}; // end anonymous namespace

namespace Go
{

//===========================================================================
void 
closestPtSurfSurfPlaneFunctional(const std::vector<Point>& epoint, //plane description
				 const std::vector<Point>& epnt1, // start pt. in surf. 1
				 const std::vector<Point>& epnt2, // start pt. in surf. 2
				 const Point& epar1, // parameter start pt. in surf. 1
				 const Point& epar2, // parameter start pt. in surf. 2
				 const ParamSurface* psurf1, // ptr. to surf. 1
				 const ParamSurface* psurf2, // ptr. to surf. 2
				 double aepsge, // absolute tolerance
				 std::vector<Point>& gpnt1, // result of iter. in surf. 1
				 std::vector<Point>& gpnt2, // result of iter. in surf. 2
				 Point& gpar1, Point& gpar2, int& jstat) // results of param.
//===========================================================================
{
    // we formulate the problem as minimizing f(x) + c_k P(x), where f(x) is the distance
    // of the points in surface 1 and surface 2, and P(x) is the distance to the plane.
    // The criteria that the points must be lying in the plane is approached as the 
    // multiplicative constant c_k increases.  By choosing a sufficiently high value
    // for c_k, we expect to converge within the specified tolerance. 

    gpar1 = gpar2 = Point(2);
    gpnt1.resize(3); gpnt2.resize(3);

    const Point plane_normal = epoint[1] / epoint[1].length();

    double C = 1; // = compute_constraint_multiplier(plane_normal, epnt1, epnt2, aepsge);

    double dist_to_plane, old_dist_to_plane;
    double seed[4];
    seed[0] = epar1[0]; seed[1] = epar1[1];
    seed[2] = epar2[0]; seed[3] = epar2[1];

    dist_to_plane = numeric_limits<double>::max();

    while(true) {
	old_dist_to_plane = dist_to_plane;
	
	// make function object here
	SurfDistFun distfun(psurf1, psurf2, C, epoint[0], plane_normal);
	
	// specify function minimizer
	FunctionMinimizer<SurfDistFun> dfmin(4, distfun, seed, aepsge);
	
	// use conjugated gradient algorithm to minimize this function
	minimise_conjugated_gradient(dfmin);//, 5); // number of iterations in each cycle
	
	dist_to_plane = 
	    compute_distance_to_plane(psurf1->point(dfmin.getPar(0), dfmin.getPar(1)),
				      psurf2->point(dfmin.getPar(2), dfmin.getPar(3)),
				      epoint[0],
				      plane_normal);
	
	if ( (dist_to_plane < aepsge) || (dist_to_plane >= old_dist_to_plane)) {
	    // we will quit the loop at this point
	    gpar1[0] = dfmin.getPar(0);
	    gpar1[1] = dfmin.getPar(1);
	    gpar2[0] = dfmin.getPar(2);
	    gpar2[1] = dfmin.getPar(3);
	    
	    gpnt1.resize(3);
	    gpnt2.resize(3);
	    psurf1->point(gpnt1, gpar1[0], gpar1[1], 1);
	    psurf2->point(gpnt2, gpar2[0], gpar2[1], 1);
	    
	    if (dist_to_plane < aepsge) {
		jstat = (gpnt1[0].dist(gpnt2[0]) < aepsge) ? 1 : 2;
	    } else {
		// we are not able to converge properly to plane
		jstat = 3;
	    }

	    // to fulfill the contract that second derivatives and normal should also
	    // be computed, we will do that here.  Note, however, that these values
	    // were used nowhere else in this algorithm.
	    gpnt1.resize(7);
	    gpnt2.resize(7);
	    psurf1->point(gpnt1, gpar1[0], gpar1[1], 2);
	    psurf2->point(gpnt2, gpar2[0], gpar2[1], 2);
	    gpnt1[6] = gpnt1[1].cross(gpnt1[2]); // nb: not normalized here!
	    gpnt2[6] = gpnt2[1].cross(gpnt2[2]); // nb: not normalized here!
	    return;
	}
	
	cout << "multiplying C!  New C is : " << C * 10 << endl;
	C *= 10; // if we redo the loop, we need to augment the value of C
	copy(dfmin.getPar(), dfmin.getPar() + 4, seed);
    } 
}

}; // end namespace Go

namespace {

//===========================================================================
double compute_constraint_multiplier(const Point& plane_normal, // plane normal, normalised
				     const vector<Point>& p1, // 0-1 derivs in p1
				     const vector<Point>& p2, // 0-1 derivs in p2
				     const double eps) // tolerance
//===========================================================================
{
    // we compute the constant C from the formula:
    // C = ||f(x)|| / (eps * min_i ||grad(g_i(x))||)
    // And then augment it somewhat to be on the "safe side".
    // f(x) is here the 'distance function' between the two points.
    const double AUGMENTATIVE_CONSTANT = 2; // should be >= 1
    
    const Point diff = p1[0] - p2[0]; // distance vector between points

    Point grad_f(4);

    SurfDistFun::compute_distance_grad(p1, p2, grad_f.begin());
    const double norm_f = grad_f.length();
    
    // calculating the minimum norm gradient function for distance to plane
    Point grad_g1(2), grad_g2(2);
    grad_g1[0] = plane_normal * p1[1];
    grad_g1[1] = plane_normal * p1[2];
    grad_g2[0] = plane_normal * p2[1];
    grad_g2[1] = plane_normal * p2[2];
    const double grad_g1_norm = grad_g1.length();
    const double grad_g2_norm = grad_g2.length();

    double min_norm = (grad_g1_norm < grad_g2_norm) ? grad_g1_norm : grad_g2_norm;
    min_norm  = (min_norm > eps) ? min_norm : eps;
    
    return AUGMENTATIVE_CONSTANT * norm_f / (eps * min_norm);
}


//===========================================================================
double compute_distance_to_plane(const Point& p1, 
				 const Point& p2, 
				 const Point& plane_point,
				 const Point& plane_normal)
//===========================================================================
{
    const double d = - plane_point * plane_normal;
    Point dist(2);
    dist[0] = plane_normal[0] * p1[0] + plane_normal[1] * p1[1] + plane_normal[2] * p1[2] + d;
    dist[1] = plane_normal[0] * p2[0] + plane_normal[1] * p2[1] + plane_normal[2] * p2[2] + d;
    return dist.length();
}


//===========================================================================
SurfDistFun::SurfDistFun(const ParamSurface* psurf1, 
			 const ParamSurface* psurf2,
			 double C,
			 const Point& plane_point,
 			 const Point& plane_normal)
//===========================================================================
    : ps1_(psurf1), ps2_(psurf2), C_(C), pl_normal_(plane_normal),
      plane_affinity_(-pl_normal_ * plane_point),
      ptvec1_(3), ptvec2_(3)
{
    RectDomain rdom = ps1_->containingDomain();
    
    min_par_[0] = rdom.umin(); max_par_[0] = rdom.umax();
    min_par_[1] = rdom.vmin(); max_par_[1] = rdom.vmax();

    rdom = ps2_->containingDomain();
    
    min_par_[2] = rdom.umin(); max_par_[2] = rdom.umax();
    min_par_[3] = rdom.vmin(); max_par_[3] = rdom.vmax();
}

//===========================================================================
inline double SurfDistFun::operator()(const double* arg) const
//===========================================================================
{
    ps1_->point(pt1_, arg[0], arg[1]);
    ps2_->point(pt2_, arg[2], arg[3]);

    const double g1 = signedDistFromPlane(pt1_);
    const double g2 = signedDistFromPlane(pt2_);

    return pt1_.dist2(pt2_) + C_ * (g1 * g1 + g2 * g2);
}

//===========================================================================
inline double SurfDistFun::signedDistFromPlane(const Point& p) const
//===========================================================================
{
    return p * pl_normal_ + plane_affinity_;
}

//===========================================================================
inline void SurfDistFun::grad(const double* arg, double* res) const
//===========================================================================
{
    ps1_->point(ptvec1_, arg[0], arg[1], 1);
    ps2_->point(ptvec2_, arg[2], arg[3], 1);

    compute_distance_grad(ptvec1_, ptvec2_, res);
    
    const double g1 = signedDistFromPlane(ptvec1_[0]);
    const double g2 = signedDistFromPlane(ptvec2_[0]);

    res[0] += C_ * 2 * g1 * pl_normal_ * ptvec1_[1];
    res[1] += C_ * 2 * g1 * pl_normal_ * ptvec1_[2];
    res[2] += C_ * 2 * g2 * pl_normal_ * ptvec2_[1];
    res[3] += C_ * 2 * g2 * pl_normal_ * ptvec2_[2];

    //return ptvec1_[0].dist2(ptvec2_[0]) + C_ * (g1 * g1 + g2 * g2);
}

//===========================================================================
inline double SurfDistFun::minPar(int pardir) const
//===========================================================================
{
    return min_par_[pardir];
}

//===========================================================================
inline double SurfDistFun::maxPar(int pardir) const
//===========================================================================
{
    return max_par_[pardir];
}

//===========================================================================
inline void SurfDistFun::compute_distance_grad(const vector<Point>& p1,
					       const vector<Point>& p2,
					       double* result)
//===========================================================================
{
    const Point diff = p1[0] - p2[0]; // distance vector between points

    result[0] = 2 * diff * p1[1]; // derivative with respect to first param.
    result[1] = 2 * diff * p1[2]; // derivative with respect to second param.
    result[2] = - 2 * diff * p2[1]; // derivative with respect to third param.
    result[3] = - 2 * diff * p2[2]; // derivative with respect to fourth param.
}

}; // end anonymous namespace

