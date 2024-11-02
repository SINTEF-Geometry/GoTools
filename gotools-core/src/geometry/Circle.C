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

#include "GoTools/geometry/Circle.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/GeometryTools.h"
#include "sisl.h"
#include <vector>
#include <fstream>

//#define DEBUG

using std::vector;
using std::istream;
using std::ostream;
using std::streamsize;
using std::cout;
using std::endl;
using std::swap;


namespace Go {


// Default constructor
//===========================================================================
  Circle::Circle()
    : ElementaryCurve()
//===========================================================================
  {
  }

// Constructor
//===========================================================================
Circle::Circle(double radius,
               Point centre, Point normal, Point x_axis,
               bool isReversed)
  : ElementaryCurve(), radius_(radius), centre_(centre),
      normal_(normal), vec1_(x_axis),
      parbound1_(0.0), parbound2_(2.0*M_PI),
      startparam_(0.0), endparam_(2.0 * M_PI)
//===========================================================================
{
    int dim = centre.dimension();
    if (dim != 2 && dim != 3)
        THROW("Dimension must be 2 or 3");

    if (dim == 3)
        normal_.normalize();
    setSpanningVectors();

    if (isReversed)
        reverseParameterDirection();
}

  // Copy constructor
//===========================================================================
Circle& Circle::operator=(const Circle& other)
//===========================================================================
{
  if (&other == this)
    return *this;
  else
    {
      radius_ = other.radius_;
      centre_ = other.centre_;
      normal_ = other.normal_;
      vec1_ = other.vec1_;
      vec2_ = other.vec2_;
      parbound1_ = other.parbound1_;
      parbound2_ = other.parbound2_;
      startparam_ = other.startparam_;
      endparam_ = other.endparam_;
      isReversed_ = other.isReversed_;
      return *this;
    }
}

// Destructor
//===========================================================================
Circle::~Circle()
//===========================================================================
{
}

//===========================================================================
void Circle::read(std::istream& is)
//===========================================================================
{
    bool is_good = is.good();
    if (!is_good) {
        THROW("Invalid geometry file!");
    }

    int dim;
    is >> dim;
    centre_.resize(dim);
    normal_.resize(dim);
    vec1_.resize(dim);
    is >> radius_
       >> centre_
       >> normal_
       >> vec1_;

    if(dim == 3)
        normal_.normalize();
    setSpanningVectors();

    is >> parbound1_ >> parbound2_;

#if 0
    // This hack is dangerous, will typically fail if this cirle is part of a CurveOnSurface.
    // Handle on the outside through setParamBounds() if necessary.
    //
    // Need to take care of rounding errors: If pars are "roughly"
    // (0, 2*M_PI) it is probably meant *exactly* (0, 2*M_PI).
    if (fabs(parbound1_) < ptol_) 
      parbound1_ = 0.0;
    if (fabs(parbound2_ - 2.0*M_PI) < ptol_)        
      parbound2_ = 2.0 * M_PI;
#endif
    startparam_ = parbound1_;
    endparam_ = parbound2_;

    // "Reset" reversion
    isReversed_ = false;

    // Swapped flag
    int isReversed; // 0 or 1
    is >> isReversed;
    bool has_param_int = (isReversed >= 10);
    if (has_param_int)
      {
	is >> startparam_ >> endparam_;
      }
    isReversed = isReversed % 10;
    if (isReversed == 0) {
        // Do nothing
    }
    else if (isReversed == 1) {
        reverseParameterDirection();
    }
    else {
        THROW("Swapped flag must be 0 or 1");
    }

}


//===========================================================================
void Circle::write(std::ostream& os) const
//===========================================================================
{
    streamsize prev = os.precision(15);
    int dim = dimension();
    os << dim << endl
       << radius_ << endl
       << centre_ << endl
       << normal_ << endl
       << vec1_ << endl
       << parbound1_ << " " << parbound2_ << endl;

    if (!isReversed()) {
        os << "10" << endl;
    }
    else {
        os << "11" << endl;
    }
    os << startparam_ << " " << endparam_ << endl;

    os.precision(prev);   // Reset precision to it's previous value    
}


//===========================================================================
BoundingBox Circle::boundingBox() const
//===========================================================================
{
    // A rather unefficient hack...
    Circle* circ = const_cast<Circle*>(this);
    SplineCurve* tmp = circ->geometryCurve();
    BoundingBox box = tmp->boundingBox();
    delete tmp;
    return box;
}


//===========================================================================
int Circle::dimension() const
//===========================================================================
{
    // Should be 2 or 3
    return centre_.dimension();
}


//===========================================================================
ClassType Circle::instanceType() const
//===========================================================================
{
    return classType();
}


//===========================================================================
ClassType Circle::classType()
//===========================================================================
{
    return Class_Circle;
}


//===========================================================================
Circle* Circle::clone() const
//===========================================================================
{
    Circle* circle = new Circle(radius_, centre_, normal_, vec1_,
        isReversed_);
    circle->setParamBounds(parbound1_, parbound2_);
    circle->setParameterInterval(startparam_, endparam_);
    return circle;
}


//===========================================================================
void Circle::point(Point& pt, double tpar) const
//===========================================================================
{
  getReversedParameter(tpar);
  tpar = parbound1_ + 
    (tpar-startparam_)*(parbound2_-parbound1_)/(endparam_-startparam_);
  pt = centre_ + radius_ * (cos(tpar) * vec1_ + sin(tpar) * vec2_);
}



//===========================================================================
void Circle::point(vector<Point>& pts,
                   double tpar,
                   int derivs,
                   bool from_right) const
//===========================================================================
{
    DEBUG_ERROR_IF(derivs < 0, 
                   "Negative number of derivatives makes no sense.");
    int totpts = (derivs + 1);
    int ptsz = (int)pts.size();
    DEBUG_ERROR_IF(ptsz < totpts, 
                   "The vector of points must have sufficient size.");

    int dim = dimension();
    for (int i = 0; i < totpts; ++i) {
        if (pts[i].dimension() != dim) {
            pts[i].resize(dim);
        }
        pts[i].setValue(0.0);
    }

    point(pts[0], tpar);
    if (derivs == 0)
        return;

    // // We use a trick that holds for a circle C(t) at the origin: The
    // // n'th derivative of C equals C(t + n*pi/2). This should work also
    // // for reversed parameters.
    getReversedParameter(tpar);
    double fac = (parbound2_-parbound1_)/(endparam_-startparam_);
    tpar = parbound1_ + fac*(tpar - startparam_);
		
    int sgn = (isReversed_) ? -1 : 1;
    for (int i = 1; i <= derivs; ++i) {
      double tpar2 = tpar + sgn*i*0.5*M_PI;
      pts[i] = radius_ * (cos(tpar2) * vec1_ + sin(tpar2) * vec2_);
      pts[i] *= fac;
      fac *= ((parbound2_-parbound1_)/(endparam_-startparam_));
    }

    return;

}


//===========================================================================
double Circle::startparam() const
//===========================================================================
{
    return startparam_;
}


//===========================================================================
double Circle::endparam() const
//===========================================================================
{
    return endparam_;
}


//===========================================================================
void Circle::swapParameters2D()
//===========================================================================
{
    if (dimension() == 2) {
        swap(centre_[0], centre_[1]);
        swap(vec1_[0], vec1_[1]);
        swap(vec2_[0], vec2_[1]);
    }
}


//===========================================================================
void Circle::setParameterInterval(double t1, double t2)
//===========================================================================
{
  startparam_ = t1;
  endparam_ = t2;
}


//===========================================================================
SplineCurve* Circle::geometryCurve()
//===========================================================================
{
    return createSplineCurve();
}


//===========================================================================
SplineCurve* Circle::createSplineCurve() const
//===========================================================================
{
    // Based on SISL function s1522.

    double tworoot = sqrt ((double) 2.0);
    double weight  = (double) 1.0 / tworoot;
    double factor = 2.0 * M_PI;

    // Knot vector
    double et[12];
    et[0] = 0.0;
    int i;
    for ( i=1;  i < 3;  i++ ) {
        et[i]     = 0.0;
        et[2 + i] = factor * 0.25;
        et[4 + i] = factor * 0.5;
        et[6 + i] = factor * 0.75;
        et[8 + i] = factor;
    }
    et[11] = factor;

    // Vertices
    double coef[36];
    int dim = dimension();
    Point axis1 = radius_ * vec1_;
    Point axis2 = radius_ * vec2_;
    if (dim == 2) {
        for ( i=0;  i < 2;  i++ ) {
            coef[     i] = centre_[i] + axis1[i];
            coef[3 +  i] = weight*(centre_[i] + axis1[i] + axis2[i]);
            coef[6 +  i] = centre_[i] + axis2[i];
            coef[9 + i] = weight*(centre_[i] - axis1[i] + axis2[i]);
            coef[12 + i] = centre_[i] - axis1[i];
            coef[15 + i] = weight*(centre_[i] - axis1[i] - axis2[i]);
            coef[18 + i] = centre_[i] - axis2[i];
            coef[21 + i] = weight*(centre_[i] + axis1[i] - axis2[i]);
            coef[24 + i] = centre_[i] + axis1[i];
        }
        // The rational weights.
        coef[2] = 1.0;
        coef[5] = weight;
        coef[8] = 1.0;
        coef[11] = weight;
        coef[14] = 1.0;
        coef[17] = weight;
        coef[20] = 1.0;
        coef[23] = weight;
        coef[26] = 1.0;
    }
    else {
        for ( i=0;  i < 3;  i++ ) {
            coef[     i] = centre_[i] + axis1[i];
            coef[4 +  i] = weight*(centre_[i] + axis1[i] + axis2[i]);
            coef[8 +  i] = centre_[i] + axis2[i];
            coef[12 + i] = weight*(centre_[i] - axis1[i] + axis2[i]);
            coef[16 + i] = centre_[i] - axis1[i];
            coef[20 + i] = weight*(centre_[i] - axis1[i] - axis2[i]);
            coef[24 + i] = centre_[i] - axis2[i];
            coef[28 + i] = weight*(centre_[i] + axis1[i] - axis2[i]);
            coef[32 + i] = centre_[i] + axis1[i];
        }
        // The rational weights.
        coef[3] = 1.0;
        coef[7] = weight;
        coef[11] = 1.0;
        coef[15] = weight;
        coef[19] = 1.0;
        coef[23] = weight;
        coef[27] = 1.0;
        coef[31] = weight;
        coef[35] = 1.0;
    }

    int ncoefs = 9;
    int order = 3;
    bool rational = true;
    SplineCurve curve(ncoefs, order, et, coef, dim, rational);

    // Extract segment. We need all this because 'curve' is not an
    // arc-length parametrized circle.
    double tmppar = parbound2_ - parbound1_;
    Point pt = centre_ + radius_*(cos(tmppar)*vec1_ + sin(tmppar)*vec2_);
    Point tmppt;
    double tmpt, tmpdist;
    double tmin = 0.0;
    double tmax = 2.0 * M_PI;
    double seed = parbound2_ - parbound1_;
    curve.closestPoint(pt, tmin, tmax, tmpt, tmppt, tmpdist, &seed);
    if (tmpt < ptol_ && parbound2_-parbound1_ == 2.0 * M_PI) {
        tmpt = 2.0 * M_PI;
    }
    SplineCurve* segment = curve.subCurve(0.0, tmpt);
    segment->basis().rescale(startparam_, endparam_);
    GeometryTools::translateSplineCurve(-centre_, *segment);
    GeometryTools::rotateSplineCurve(normal_, parbound1_, *segment);
    GeometryTools::translateSplineCurve(centre_, *segment);

    if (isReversed())
        segment->reverseParameterDirection();
#ifdef DEBUG
    std::ofstream of("circle_seg.g2");
    writeStandardHeader(of);
    write(of);
    segment->writeStandardHeader(of);
    segment->write(of);
#endif
    return segment;
}


//===========================================================================
bool Circle::isDegenerate(double degenerate_epsilon)
//===========================================================================
{
    // We consider a Circle as degenerate if the radius is smaller
    // than the epsilon.

    return radius_ < degenerate_epsilon;
}


//===========================================================================
Circle* Circle::subCurve(double from_par, double to_par,
                         double fuzzy) const
//===========================================================================
{
  Circle* circle = clone();
  if (isReversed())
    {
      double start = endparam_ - (to_par - startparam_);
      double end = startparam_ + (endparam_ - from_par);
      if (start > end)
	{
	  std::swap(start, end);
	}
      double bound1 = parbound1_ + 
	(start-startparam_)*(parbound2_-parbound1_)/(endparam_-startparam_);
      double bound2 = parbound1_ + 
	(end-startparam_)*(parbound2_-parbound1_)/(endparam_-startparam_);
      circle->setParamBounds(bound1, bound2);
      if (from_par > to_par)
	std::swap(from_par, to_par);
      circle->setParameterInterval(from_par, to_par);
    }
  else
    {
      double bound1 = parbound1_ + 
	(from_par-startparam_)*(parbound2_-parbound1_)/(endparam_-startparam_);
      double bound2 = parbound1_ + 
	(to_par-startparam_)*(parbound2_-parbound1_)/(endparam_-startparam_);
      circle->setParamBounds(bound1, bound2);
      circle->setParameterInterval(from_par, to_par);
    }
  return circle;
}

//===========================================================================
SplineCurve* Circle::createNonRationalSpline(double eps) const
//===========================================================================
{						
  // Create sisl curve
  int status;
  SISLCurve *sislcrv = NULL;
  Point startpt = centre_ + radius_ * (cos(parbound1_) * vec1_ + 
				       sin(parbound1_) * vec2_);
  double *pt1 = startpt.begin();
  double *mid = const_cast<double*>(centre_.begin());
  double *axis = const_cast<double*>(normal_.begin());
  s1303(pt1, eps, parbound2_-parbound1_, mid,
	axis, centre_.dimension(), &sislcrv, &status);
  if (status < 0 || sislcrv == NULL)
    return NULL;  // Approximation failed

  SplineCurve *crv = new SplineCurve(sislcrv->in, sislcrv->ik, sislcrv->et,
				     sislcrv->ecoef, centre_.dimension());
  crv->setParameterInterval(startparam_, endparam_);
  if (isReversed_)
    crv->reverseParameterDirection();

  freeCurve(sislcrv);
  return crv;
}

//===========================================================================
DirectionCone Circle::directionCone() const
//===========================================================================
{
    double tmin = startparam();
    double tmax = endparam();
    vector<Point> pts(2);
    point(pts, 0.5*(tmin+tmax), 1);
    return DirectionCone(pts[1], 0.5*(tmax-tmin));
}


//===========================================================================
void Circle::appendCurve(ParamCurve* cv, bool reparam)
//===========================================================================
{
  double dist;
  appendCurve(cv, 0, dist, reparam);
}


//===========================================================================
void Circle::appendCurve(ParamCurve* cv, int continuity, double& dist,
			 bool reparam, double tol)
//===========================================================================
{
  // Check input
  if (cv->instanceType() != Class_Circle)
    THROW("Inconsistency in curve types for appendCurve"); 
  
  double eps = 1.0e-5;
  double angtol = 0.01;
  Circle *other = dynamic_cast<Circle*>(cv);
  double cdist = centre_.dist(other->getCentre());
  double ang = normal_.angle(other->getNormal());
  double drad = fabs(radius_-other->getRadius());
  if (cdist > eps || drad > eps || ang > angtol ||
      (isReversed_ && !other->isReversed()) ||
      (!isReversed_ && other->isReversed()))
    THROW("AppendCurve: Circle descriptions not compatible");

  // Check also consistency of circle segment
  double bound1 = other->parbound1_;
  double bound2 = other->parbound2_;
  double axis_ang = vec1_.angle(other->getXAxis());
  if (axis_ang > angtol)
    {
      THROW("AppendCurve, Circle: Non-compatible axes");
      // Adjust parameter bounds of the other curve
      if (true)  // Too simple. The sign depends on the direction of
	// the rotation and whether or not the curves are reversed
	{
	  bound1 += axis_ang;
	  bound2 += axis_ang;
	}
    }
  double fac = (endparam_ - startparam_)/(parbound2_ - parbound1_);
  if (isReversed())
    {
      if (bound2 > parbound1_ + ptol_)
	{
	  bound1 -= 2*M_PI;
	  bound2 -= 2*M_PI;
	}
      if (fabs(bound2-parbound1_) > eps)
	THROW("AppendCurve: Circle segments not continuous");

      parbound1_ -= (bound2 - bound1);
      startparam_ -= fac*(bound2 - bound1);
      if (parbound2_ - parbound1_ > 2*M_PI)
	{
	  parbound2_ = parbound1_ + 2*M_PI;
	  endparam_ = startparam_ + 2*fac*M_PI;
	}
     }
  else
    {
      if (bound1 < parbound2_ - ptol_)
	{
	  bound1 += 2*M_PI;
	  bound2 += 2*M_PI;
	}
      else if (bound1 > parbound2_ + eps)
	{
	  bound1 -= 2*M_PI;
	  bound2 += 2*M_PI;
	}
      if (fabs(bound1-parbound2_) > eps)
	THROW("AppendCurve: Circle segments not continuous");

      parbound2_ += (bound2 - bound1);
      endparam_ += fac*(bound2 - bound1);
      if (parbound2_ - parbound1_ > 2*M_PI)
	{
	  parbound2_ = parbound1_ + 2*M_PI;
	  endparam_ = startparam_ + 2*fac*M_PI;
	}
    }

  if (parbound1_ < -2.0 * M_PI)
    {
      parbound1_ += 2.0*M_PI;
      parbound2_ += 2.0*M_PI;
    }
  if (parbound2_ > 2.0*M_PI)
    {
      parbound1_ -= 2.0*M_PI;
      parbound2_ -= 2.0*M_PI;
    }

  dist = cdist + drad;
}


//===========================================================================
void Circle::closestPoint(const Point& pt,
                        double tmin,
                        double tmax,
                        double& clo_t,
                        Point& clo_pt,
                        double& clo_dist,
                        double const *seed) const
//===========================================================================
{
    double clo_t2;

    // Check and fix the parameter bounds
    if (tmin < startparam_) {
        MESSAGE("tmin too small. Using startparam_.");
        tmin = startparam_;
    }
    if (tmax > endparam_) {
        MESSAGE("tmax too large. Using endparam_.");
        tmax = endparam_;
    }

    tmin = parbound1_ + (tmin-startparam_)*(parbound2_-parbound1_)/(endparam_-startparam_);
    tmax = parbound1_ + (tmax-startparam_)*(parbound2_-parbound1_)/(endparam_-startparam_);

    // If input is on the "centre line", we arbitrarily assign the
    // point with t = tmin.
    Point vec = pt - centre_;
    Point tmp = (dimension() == 3) ? vec.cross(normal_) : vec;
    if (tmp.length() == 0.0) {
        clo_t = (seed != NULL) ? *seed : tmin;
	clo_pt = centre_ + radius_*(cos(clo_t)*vec1_ + sin(clo_t)*vec2_);
        clo_dist = radius_;
        //MESSAGE("Input to Circle::closestPoint() is the centre.");
        return;
    }

    Point proj;
    if (dimension() == 2)
        proj = vec;
    else if (dimension() == 3)
        proj = vec - (vec * normal_) * normal_;
    else
        THROW("Dimension must be 2 or 3");
    double x = proj * vec1_;
    double y = proj * vec2_;
    if (x == 0.0) {
        if (y > 0.0) {
            clo_t = 0.5 * M_PI;
            clo_pt = centre_ + radius_ * vec2_;
        }
        else {
            clo_t = 1.5 * M_PI;
            clo_pt = centre_ - radius_ * vec2_;
        }
	clo_t2 = startparam_ + (clo_t-parbound1_)*(endparam_-startparam_)/(parbound2_-parbound1_);
        getReversedParameter(clo_t2);
    }
    else {
        clo_t = atan(y / x);
        // We need to correct the angle when we are in quadrants II,
        // III and IV
        if (x < 0.0)
            clo_t += M_PI; // II + III
        if (x > 0.0 && y < 0.0)
            clo_t += 2.0 * M_PI; // IV

	// If we are epsilon-close to the seam and were given a seed, we may want to move to the other side of the seam.
	if ((seed != NULL) && ((clo_t < ptol_) || (fabs(2.0*M_PI - clo_t) < ptol_)))
	{
	    clo_t = (*seed < M_PI) ? 0.0 : 2.0*M_PI;
	}

	clo_t2 = startparam_ + (clo_t-parbound1_)*(endparam_-startparam_)/(parbound2_-parbound1_);
        getReversedParameter(clo_t2);
        point(clo_pt, clo_t2);
    }
    clo_dist = (clo_pt - pt).length();
    // if (seed) // We do not want to move the point to the other side of the seam.
    // {
    //     return;
    // }
    
    // We must handle the case of a proper circle segment. Meaning a subcurve of the whole curve.
    double tlen = tmax - tmin;
    double tmp_t = clo_t - tmin;
    // We first make sure that we are inside the valid range.
    if (tmp_t > 2.0 * M_PI + ptol_)
        tmp_t -= 2.0 * M_PI;
    else if (tmp_t < - ptol_)
        tmp_t += 2.0 * M_PI;
    bool is_subcv = (tlen < 2.0*M_PI - ptol_);
    if (is_subcv)
    {
        if (tmp_t >= 0.5 * tlen + M_PI) {
            // Start of segment is closest
            clo_t = tmin;
	    clo_t2 = startparam_ + 
	      (clo_t-parbound1_)*(endparam_-startparam_)/(parbound2_-parbound1_);
	    getReversedParameter(clo_t2);                
            point(clo_pt, clo_t2);
            clo_dist = (clo_pt - pt).length();
	    clo_t = clo_t2;
            return;
        }
        if (tmp_t >= tlen) {
                // End of segment is closest
                clo_t = tmax;
		clo_t2 = startparam_ + 
		  (clo_t-parbound1_)*(endparam_-startparam_)/(parbound2_-parbound1_);
		getReversedParameter(clo_t2);                
		point(clo_pt, clo_t2);
                clo_dist = (clo_pt - pt).length();
		clo_t = clo_t2;
                return;
        }
    }
    // If we get here, point on segment is closest
    clo_t = tmp_t + tmin;
    clo_t2 = 
      startparam_ + (clo_t-parbound1_)*(endparam_-startparam_)/(parbound2_-parbound1_);
    getReversedParameter(clo_t2);
    point(clo_pt, clo_t2);
    clo_t = clo_t2;
    clo_dist = (clo_pt - pt).length();

}


//===========================================================================
double Circle::length(double tol)
//===========================================================================
{
    return (parbound2_ - parbound1_) * radius_;
}


//===========================================================================
void Circle::setParamBounds(double startpar, double endpar)
//===========================================================================
{
  if (fabs(startpar) < ptol_)
    startpar = 0.0;
  else if (fabs(2.0*M_PI-startpar) < ptol_)
    startpar = 2.0*M_PI;
  if (fabs(endpar) < ptol_)
    endpar = 0.0;
  else if (fabs(2.0*M_PI-endpar) < ptol_)
    endpar = 2.0*M_PI;
  
  if (startpar > -2.0 * M_PI - ptol_ && startpar < -2.0 * M_PI)
    startpar = -2.0 * M_PI;
  if (endpar < 2.0 * M_PI + ptol_ && endpar >2.0 * M_PI)
    endpar = 2.0 * M_PI;
  if (startpar >= endpar)
    THROW("First parameter must be strictly less than second.");
  if (startpar < -2.0 * M_PI || endpar > 2.0 * M_PI)
    THROW("Parameters must be in [-2pi, 2pi].");
  if (endpar - startpar > 2.0 * M_PI)
    THROW("(endpar - startpar) must not exceed 2pi.");

  double start =  
    parbound1_ + (startpar-startparam_)*(parbound2_-parbound1_)/(endparam_-startparam_);
  double end =  
    parbound1_ + (endpar-startparam_)*(parbound2_-parbound1_)/(endparam_-startparam_);
  parbound1_ = startpar;
  parbound2_ = endpar;
  startparam_ = start;
  endparam_ = end;
}

//===========================================================================
bool Circle::isClosed() const
//===========================================================================
{
  return (parbound2_ - parbound1_ >= 2.0*M_PI - ptol_ &&
	  parbound2_ - parbound1_ <= 2.0*M_PI + ptol_);
}


//===========================================================================
bool Circle::isAxisRotational(Point& centre, Point& axis, Point& vec,
			      double& angle)
//===========================================================================
{
  double radius;
  return isAxisRotational(centre, axis, vec, angle, radius);
}

//===========================================================================
bool Circle::isAxisRotational(Point& centre, Point& axis, Point& vec,
			      double& angle, double& radius)
//===========================================================================
{
  centre = centre_;
  axis = normal_;
  radius = radius_;
  double tol = 1.0e-10;
  if (isClosed() || radius_ < tol)
    {
      vec = vec1_;
      angle = 2.0*M_PI;
    }
  else
    {
      Point pt;
      if (isReversed())
	pt = ParamCurve::point(endparam_);
      else
	pt = ParamCurve::point(startparam_);
      vec = pt - centre_;
      vec.normalize();
      angle = endparam_ - startparam_;
    }
  return true;
}

//===========================================================================
void Circle::translateCurve(const Point& dir)
//===========================================================================
{
  centre_ += dir;
}

//===========================================================================
void Circle::setSpanningVectors()
//===========================================================================
{
    // In 3D, the spanning vectors vec1_, vec2_, and the vector
    // normal_ defines a right-handed coordinate system. Similar to an
    // axis2_placement_3d entity in STEP.

    int dim = centre_.dimension();
    if (dim == 2) {
        vec2_.resize(2);
        vec2_[0] = -vec1_[1];
        vec2_[1] = vec1_[0];
    }
    else if (dim ==3) {
        Point tmp = vec1_ - (vec1_ * normal_) * normal_;
        if (tmp.length() == 0.0) 
            THROW("X-axis parallel to normal.");
        vec1_ = tmp;
        vec2_ = normal_.cross(vec1_);
    }
    else {
        THROW("Dimension must be 2 or 3");
    }
    vec1_.normalize();
    vec2_.normalize();

}


//===========================================================================
bool Circle::isInPlane(const Point& loc, const Point& axis,
		       double eps, Point& normal) const
//===========================================================================
{
  normal = normal_;

  Point vec = normal.cross(axis);
  if (vec.length() < eps)
    return false;  // The circle lies in a plane orthogonal to the axis

  vec.normalize();
  double dist = (centre_ - loc)*normal;
  return (dist < eps);
}

//===========================================================================
bool Circle::isInPlane(const Point& norm,
		       double eps, Point& pos) const

//===========================================================================
{
  double ang = norm.angle(normal_);
  pos = centre_;

  return (ang <= eps || fabs(M_PI-ang) <= eps);
}

} // namespace Go
