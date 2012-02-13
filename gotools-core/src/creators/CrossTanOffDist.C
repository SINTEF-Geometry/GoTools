//===========================================================================
//                                                                           
// File: CrossTanOffDist.C                                                      
//                                                                           
// Created: 02.05.13
//                                                                           
// Author: Vibeke Skytt, SINTEF
//                                                                           
// Revision: 
//                                                                           
// Description: Blending between a set of curves seen as an evaluator
//              based curve.
//===========================================================================
//

#include "GoTools/creators/CrossTanOffDist.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/utils/errormacros.h"

using namespace Go;
using std::vector;
using std::max;
using std::min;

//===========================================================================

CrossTanOffDist::CrossTanOffDist(shared_ptr<SplineCurve>& poscurve,
				 shared_ptr<SplineCurve>& tangcv1,
				 shared_ptr<SplineCurve>& tangcv2,
				 shared_ptr<SplineCurve>& blend1,
				 shared_ptr<SplineCurve>& blend2,
				 shared_ptr<SplineCurve>& opposite1,
				 shared_ptr<SplineCurve>& opposite2,
				 double factor)
  : poscurve_(poscurve), epstol_(1.0e-15)
  //===========================================================================
{
  // Test input
  ALWAYS_ERROR_IF(poscurve.get() == 0 || tangcv1.get() == 0 || 
		  tangcv2.get() == 0 || blend1.get() == 0 || 
		  blend2.get() == 0 || opposite1.get() == 0 ||
		  opposite2.get() == 0, 
		  "Missing curve");
  ALWAYS_ERROR_IF(poscurve->dimension() != 3 ||
		  poscurve->dimension() != tangcv1->dimension() ||
		  poscurve->dimension() != tangcv2->dimension() ||
		  poscurve->dimension() != opposite1->dimension() ||
		  poscurve->dimension() != opposite2->dimension(),
		  "Dimension mismatch");
  ALWAYS_ERROR_IF(blend1->dimension() != 1 || blend2->dimension() != 1,
		  "Blending function of dimension different from 1");

  //poscurve_ = poscurve;
  tangcurves_.push_back(tangcv1);
  tangcurves_.push_back(tangcv2);
  blends_.push_back(blend1);
  blends_.push_back(blend2);

  double startpar = poscurve->startparam();
  double endpar = poscurve->endparam();
  double midpar = 0.5*(startpar+endpar);
  double dist;
  double minlen = 0.0001;

  // Create the opposite position curve.
  if (opposite1.get() == opposite2.get() && factor == 1.0)
    {
      // The same opposite curve for both halves of the cross tangent
      // curve
 #if ((_MSC_VER > 0) && (_MSC_VER < 1300))
     oppositepos_ = shared_ptr<SplineCurve>
       (dynamic_cast<SplineCurve*>(opposite1->clone()));
#else
     oppositepos_ = shared_ptr<SplineCurve>(opposite1->clone());
#endif
      oppositepos_->reverseParameterDirection();
    }
  else
    {
      // The opposite curve must be computed.
      // Fetch curve pieces
      shared_ptr<SplineCurve> sub1, sub2, sub3, sub4, sum2;
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
      sub1 = shared_ptr<SplineCurve>(dynamic_cast<SplineCurve*>
				      (opposite1->subCurve(midpar, endpar)));
      sub2 = shared_ptr<SplineCurve>(dynamic_cast<SplineCurve*>
				      (opposite2->subCurve(startpar, midpar)));
      sub3 = shared_ptr<SplineCurve>(dynamic_cast<SplineCurve*>
				      (poscurve->subCurve(startpar, midpar)));
      sub4 = shared_ptr<SplineCurve>(dynamic_cast<SplineCurve*>
				      (poscurve->subCurve(midpar, endpar)));
#else
      sub1 = shared_ptr<SplineCurve>(opposite1->subCurve(midpar, endpar));
      sub2 = shared_ptr<SplineCurve>(opposite2->subCurve(startpar, midpar));
      sub3 = shared_ptr<SplineCurve>(poscurve->subCurve(startpar, midpar));
      sub4 = shared_ptr<SplineCurve>(poscurve->subCurve(midpar, endpar));
#endif
      sub1->reverseParameterDirection();
      sub1->setParameterInterval(startpar, midpar);
      sub2->reverseParameterDirection();
      sub2->setParameterInterval(midpar, endpar);
      
      try {
	  oppositepos_ = shared_ptr<SplineCurve>(GeometryTools::curveSum(*sub1, 0.5,
							  *sub3, 0.5));
	  sum2 = shared_ptr<SplineCurve>(GeometryTools::curveSum(*sub2, 0.5,
						  *sub4, 0.5));
      } catch (...) {
	  THROW("Failed adding curves.");
      }
      oppositepos_->appendCurve(sum2.get(), 0, dist);
    }

  // Make length factor curve. First compute the fraction between the
  // tangent length in the start of the curve and the distance to the
  // opposite curve in the start.
  Point crtan = evalblend(startpar);
  Point diffvec = evaldiff(startpar);
  double len = std::max(minlen, diffvec.length());
  double fac1 = crtan.length()/len;

  // Compute fraction in the endpoint.
  crtan = evalblend(endpar);
  diffvec = evaldiff(endpar);
  len = std::max(minlen, diffvec.length());
  double fac2 = crtan.length()/len;
  
  // Make curve
  Point pt1(1), pt2(1);
  pt1.setValue(fac1);
  pt2.setValue(fac2);
  lengthfac_ = shared_ptr<SplineCurve>(new SplineCurve(pt1, startpar, pt2, endpar));

  // Make linear curve interpolating the blended cross tangent curve
  // in the endpoints.
  Point cross1, cross2;
  cross1 = evalblend(startpar);
  cross2 = evalblend(endpar);
  avcross_ = shared_ptr<SplineCurve>
    (new SplineCurve(cross1, startpar, cross2, endpar));
}

//===========================================================================

CrossTanOffDist::~CrossTanOffDist()
//===========================================================================
{
}

//===========================================================================
Point CrossTanOffDist::eval( double t) const
//===========================================================================
{
  int dim = poscurve_->dimension();
  Point pos(dim), cross(dim), len(1), proj(dim), diffvec(dim);

  poscurve_->point(pos, t);
  evalcrtan(t, cross, proj);
  diffvec = evaldiff(t);

  // Modify the cross tangent vector. 
  double startpar = poscurve_->startparam();
  double endpar = poscurve_->endparam();
  double del = 0.1*(endpar - startpar);
  double fac = (t < startpar+del || t > endpar-del) ?
    1.0 - std::min(t-startpar,endpar-t)/del : 0.0;
  cross = fac*cross + (1.0-fac)*proj;

//    cross = proj;

  lengthfac_->point(len,t);

  if (cross.length() < epstol_)
    return pos;
  else
    return pos + cross*len[0]*diffvec.length()/cross.length();
}

//===========================================================================
void CrossTanOffDist::eval(double t, int n, Point der[]) const
//===========================================================================
{
  MESSAGE_IF(n > 1, "Only one derivative will be computed");

  if (n == 0)
    der[0] = eval(t);
  else
    {
      Point blendder[2], projder[2], cross[2], hc[2], diffvec[2];
      vector<Point> pos(2), len(2);
      poscurve_->point(pos, t, 1);
      evalcrtan(t, 1, blendder, projder);
      lengthfac_->point(len, t, 1);
      evaldiff(t, 1, diffvec);

      double startpar = poscurve_->startparam();
      double endpar = poscurve_->endparam();
      double del = 0.1*(endpar - startpar);
      double fac = (t < startpar+del || t > endpar-del) ?
	1.0 - std::min(t-startpar,endpar-t)/del : 0.0;
      
      cross[0] = fac*blendder[0] + (1.0-fac)*projder[0];

      double dfac = (t < startpar+del || t > endpar-del) ?
	((t < startpar+del) ? -1.0/del : 1.0/del) : 0.0;

      cross[1] = dfac*blendder[0] + fac*blendder[1]
	- dfac*projder[0] + (1-fac)*projder[1];

      double clen = cross[0].length();
      double dlen = diffvec[0].length();
      if (clen < epstol_ || dlen < epstol_)
	{
	  der[0] = pos[0];
	  der[1] = pos[1];
	}
      else
	{
	  hc[0] = cross[0]/clen;
	  hc[1] = (1.0 - hc[0]*hc[0])*cross[1]/clen;

	  der[0] = pos[0] + hc[0]*len[0][0]*dlen;
	  der[1] = pos[1] + hc[1]*len[0][0]*dlen + hc[0]*len[1][0]*dlen
	    + hc[0]*len[0][0]*diffvec[0]*diffvec[1]/dlen;
	}
    }
}


//===========================================================================
double CrossTanOffDist::start() const
//===========================================================================
{
  return poscurve_->startparam();
}


//===========================================================================

double CrossTanOffDist::end() const
//===========================================================================
{
  return poscurve_->endparam();
}

//===========================================================================

int CrossTanOffDist::dim() const
//===========================================================================
{
  return poscurve_->dimension();
}

//===========================================================================

bool CrossTanOffDist::approximationOK(double par, Point approxpos,
				      double tol1, double tol2) const
//===========================================================================
{
  double tol3 = 0.000001*tol1;
  Point pos = eval(par);
  double dist = pos.dist(approxpos);  // Distance between original
                                      // curve and approximation
  if (dist > tol1)
    return false;   // Approximation not good enough

  if (dist < tol3)
    return true;

  Point pospt;
  poscurve_->point(pospt, par);
  Point diff = approxpos - pospt;
  Point pt1, pt2;
  tangcurves_[0]->point(pt1, par);
  tangcurves_[1]->point(pt2, par);
  Point normal = pt1 % pt2;
  double ang = normal.angle(diff);
  double pihalf = 3.141592653589793/2.0;
  if (fabs(pihalf-ang) > tol2)
    return false;   // Approximated cross tangent not in tangent plane

  return true;
}


//===========================================================================

void CrossTanOffDist::evalcrtan(double t, Point& blend, 
				Point& projdiff) const
//===========================================================================
{
  int dim = poscurve_->dimension();
  Point point1(dim), point2(dim);
  Point fac1, fac2;

  tangcurves_[0]->point(point1, t);
  tangcurves_[1]->point(point2, t);
  blends_[0]->point(fac1, t);
  blends_[1]->point(fac2, t);
  blend = fac1[0]*point1 + fac2[0]*point2;

  Point avvec(dim);
  avcross_->point(avvec, t);
  point1.normalize();
  point2.normalize();
  double ab = point1*point2;
  double tdiv = 1.0 - ab*ab;
  double t1 = (point1*avvec - ab*(point2*avvec))/tdiv;
  double t2 = (point2*avvec - ab*(point1*avvec))/tdiv;
  
  projdiff = t1*point1 + t2*point2;
}

//===========================================================================

Point CrossTanOffDist::evalblend( double t) const
//===========================================================================
{
  int dim = poscurve_->dimension();
  Point point1(dim), point2(dim);
  Point fac;
  point1.setValue(0.0);

  for (size_t ki=0; ki<tangcurves_.size(); ki++)
    {
      tangcurves_[ki]->point(point2, t);
      blends_[ki]->point(fac, t);
      point1 += fac[0]*point2;
    }

  return point1;
}

//===========================================================================

Point CrossTanOffDist::evaldiff( double t) const
//===========================================================================
{
  int dim = poscurve_->dimension();
  Point point1(dim), point2(dim);
  poscurve_->point(point1, t);
  oppositepos_->point(point2, t);

  return point2 - point1;
}

//===========================================================================

void CrossTanOffDist::evaldiff( double t, int n, Point der[]) const
//===========================================================================
{
  MESSAGE_IF(n > 1, "Only one derivative will be computed");

  if (n == 0)
    der[0] = evaldiff(t);
  else
    {
      vector<Point> point1(2), point2(2);
      poscurve_->point(point1, t, 1);
      oppositepos_->point(point2, t, 1);

      der[0] = point2[0] - point1[0];
      der[1] = point2[1] - point1[1];
    }
}

//===========================================================================

void CrossTanOffDist::evalcrtan(double t, int n, Point derblend[],
				  Point derproj[]) const
//===========================================================================
{
  MESSAGE_IF(n > 1, "Only one derivative will be computed");

  if (n == 0)
    evalcrtan(t, derblend[0], derproj[0]);
  else
    {
//       int dim = poscurve_->dimension();
      vector<Point> point1(2), point2(2), fac1(2), fac2(2);
      vector<Point> len(2), avvec(2); 
      
      tangcurves_[0]->point(point1, t, 1);
      tangcurves_[1]->point(point2, t, 1);
      blends_[0]->point(fac1, t, 1);
      blends_[1]->point(fac2, t, 1);
      derblend[0] = fac1[0][0]*point1[0] + fac2[0][0]*point2[0];
      derblend[1] = fac1[1][0]*point1[0] + fac1[0][0]*point1[1] +
	fac2[1][0]*point2[0] + fac2[0][0]*point2[1];

      avcross_->point(avvec, t, 1);
      double ab = point1[0]*point2[0];
      double aa = point1[0]*point1[0];
      double bb = point2[0]*point2[0];
      double ac = point1[0]*avvec[0];
      double bc = point2[0]*avvec[0];
      double daa = 2.0*point1[0]*point1[1];
      double dbb = 2.0*point2[0]*point2[1];
      double dab = point1[0]*point2[1] + point1[1]*point2[0];
      double dac = point1[1]*avvec[0]+point1[0]*avvec[1];
      double dbc = point2[1]*avvec[0]+point2[0]*avvec[1];
      double tdiv = aa*bb - ab*ab;
      double t1 = (bb*ac - ab*bc)/tdiv;
      double t2 = (aa*bc - ab*ac)/tdiv;

      derproj[0] = t1*point1[0] + t2*point2[0];

      double dtdiv = daa*bb + aa*dbb - 2.0*ab*dab;
      double du = dbb*ac + bb*dac - dab*bc - ab*dbc;
      double dv = daa*bc + aa*dbc - dab*ac - ab*dac;

      derproj[1] = point1[0]*(du-t1*dtdiv)/tdiv + t1*point1[1] +
	point2[0]*(dv-t2*dtdiv)/tdiv + t2*point2[1];
    }
}


