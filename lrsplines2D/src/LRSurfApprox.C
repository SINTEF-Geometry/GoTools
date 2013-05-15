// ----------------------------------------------------------------
//       Implementation file for class LRSurfApprox
// ----------------------------------------------------------------

#include "GoTools/lrsplines2D/LRSurfApprox.h"
#include "GoTools/lrsplines2D/LRSurfSmoothLS.h"
#include "GoTools/lrsplines2D/LRBSpline2DUtils.h"
#include "GoTools/creators/ApproxSurf.h"
#include "GoTools/geometry/PointCloud.h"
#include <iostream>
#include <fstream>

//#define DEBUG

using std::vector;

using namespace Go;

//==============================================================================
LRSurfApprox::LRSurfApprox(vector<double>& points, 
			   int dim, double epsge, bool closest_dist,
			   bool repar)
  : points_(points), maxdist_(-10000.0), avdist_(0.0), outsideeps_(0), aepsge_(epsge),
    smoothweight_(1.0e-3), repar_(repar), check_close_(closest_dist)
//==============================================================================
{
  edge_derivs_[0] = edge_derivs_[1] = edge_derivs_[2] = edge_derivs_[3] = 0;

  // Create an LR B-spline surface with the domain given by the 
  // parameter domain of the points. Only one element will be
  // created
  makeInitSurf(dim);
}

//==============================================================================
LRSurfApprox::LRSurfApprox(shared_ptr<SplineSurface>& srf,
			   vector<double>& points, 
			   double epsge, bool closest_dist,
			   bool repar)
  : points_(points), maxdist_(-10000.0), avdist_(0.0), outsideeps_(0), aepsge_(epsge),
    smoothweight_(1.0e-3), repar_(repar), check_close_(closest_dist)
//==============================================================================
{
  edge_derivs_[0] = edge_derivs_[1] = edge_derivs_[2] = edge_derivs_[3] = 0;

  // Create an LR B-spline surface based on the given spline surface
  makeInitSurf(srf);
}

//==============================================================================
LRSurfApprox::LRSurfApprox(shared_ptr<LRSplineSurface>& srf,
			   vector<double>& points, 
			   double epsge, bool closest_dist,
			   bool repar)
//==============================================================================
  : points_(points), maxdist_(-10000.0), avdist_(0.0), outsideeps_(0), aepsge_(epsge),
    smoothweight_(1.0e-3), repar_(repar), check_close_(closest_dist)
{
  edge_derivs_[0] = edge_derivs_[1] = edge_derivs_[2] = edge_derivs_[3] = 0;
  srf_ = srf;
}

//==============================================================================
LRSurfApprox::LRSurfApprox(int ncoef_u, int order_u, int ncoef_v, int order_v,
			   vector<double>& points, 
			   int dim, double epsge, bool closest_dist,
			   bool repar)
//==============================================================================
  : points_(points), maxdist_(-10000.0), avdist_(0.0), outsideeps_(0), aepsge_(epsge),
    smoothweight_(1.0e-3), repar_(repar), check_close_(closest_dist)
{
  edge_derivs_[0] = edge_derivs_[1] = edge_derivs_[2] = edge_derivs_[3] = 0;

  // Create an LR B-spline surface with unset coefficients and the domain
  // given by the parameter domain of the points. The size of the spline
  // space is given. The knots will be equally spaced
  makeInitSurf(dim, ncoef_u, order_u, ncoef_v, order_v);
}

//==============================================================================
LRSurfApprox::~LRSurfApprox()
//==============================================================================
{
}

//==============================================================================
 shared_ptr<LRSplineSurface> LRSurfApprox::getApproxSurf(double& maxdist, 
							 double& avdist,
							 int& nmb_out_eps, 
							 int max_iter)
//==============================================================================
{
  // Initiate approximation engine
  LRSurfSmoothLS LSapprox(srf_, coef_known_);

  // Initiate with data points
  LSapprox.addDataPoints(points_);

  // Initial smoothing of LR B-spline surface
  performSmooth(&LSapprox);

#ifdef DEBUG
  std::ofstream of1("init_sf.g2");
  shared_ptr<LRSplineSurface> tmp(srf_->clone());
  tmp->to3D();
  tmp->writeStandardHeader(of1);
  tmp->write(of1);
  of1 << std::endl;
#endif

  // Compute accuracy in data points
  computeAccuracy();

  for (int ki=0; ki<max_iter; ++ki)
    {
      // Check if the requested accuracy is reached
      if (maxdist_ <= aepsge_ || outsideeps_ == 0)
	break;

      // Refine surface
      refineSurf();
#ifdef DEBUG
  std::ofstream of2("refined_sf.g2");
  shared_ptr<LRSplineSurface> tmp2(srf_->clone());
  tmp2->to3D();
  tmp2->writeStandardHeader(of2);
  tmp2->write(of2);
  of2 << std::endl;

  std::ofstream of3("point_clouds.g2");
  for (LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
       it != srf_->elementsEnd(); ++it)
    {
      vector<double>& elem_data = it->second->getDataPoints();
      if (elem_data.size() > 0)
	{
	  PointCloud3D cloud(elem_data.begin(), elem_data.size()/3);
	  cloud.writeStandardHeader(of3);
	  cloud.write(of3);
	}
    }
#endif

 
      // Update surface
      LSapprox.updateLocals();
      performSmooth(&LSapprox);

#ifdef DEBUG
  std::ofstream of4("updated_sf.g2");
  shared_ptr<LRSplineSurface> tmp3(srf_->clone());
  tmp3->to3D();
  tmp3->writeStandardHeader(of4);
  tmp3->write(of4);
  of4 << std::endl;
#endif

       computeAccuracy();
    }

  // Set accuracy information
  maxdist = maxdist_;
  avdist = avdist_;
  nmb_out_eps = outsideeps_;
  
  return srf_;
}

//==============================================================================
void LRSurfApprox::performSmooth(LRSurfSmoothLS *LSapprox)
//==============================================================================
{
  double weight = 1.0;  // For the time being (until smoothing is implemented)
  LSapprox->setLeastSquares(weight);

  shared_ptr<LRSplineSurface> lrsf_out;
  int isOK = LSapprox->equationSolve(lrsf_out);
#ifdef DEBUG
  std::cout << "isOK: " << isOK << std::endl;
#endif

  srf_ = lrsf_out;
}

//==============================================================================
void LRSurfApprox::computeAccuracy()
//==============================================================================
{
  // Check the accuracy of all data points, element by element
  // Note that only points more distant from the surface than the tolerance
  // are considered in avdist_

  // Initiate accuracy information
  maxdist_ = 0.0;
  avdist_ = 0.0;
  outsideeps_ = 0;

  int dim = srf_->dimension();
  int del = 2 + dim;
  for (LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
       it != srf_->elementsEnd(); ++it)
    {
      if (!it->second->hasDataPoints())
	continue;   // No points in which to check accuracy
      
      vector<double> points = it->second->getDataPoints();
      int nmb = it->second->nmbDataPoints();

      // Local error information
      double max_err = 0.0;
      double av_err = 0.0;
      int outside = 0;

      int ki;
      double *curr;
      for (ki=0, curr=&points[0]; ki<nmb; ++ki, curr+=del)
	{
	  double dist;
	  if (check_close_ && dim == 3)
	    {
	      // Compute closest point
	      // VSK. 052013. Should be changed to make use of the fact that we know
	      // the element (at least initially)
	      double upar, vpar;
	      Point close_pt;
	      srf_->closestPoint(Point(curr+2, curr+del), upar, vpar, close_pt,
				 dist, aepsge_, NULL, curr);
	      if (repar_)
		{
		  curr[0] = upar;
		  curr[1] = vpar;
		}
	    }
	  else
	    {
	      // Evaluate
	      Point pos;
	      srf_->point(pos, curr[0], curr[1]);
	      dist = pos.dist(Point(curr+2, curr+del));
	    }

	  // Accumulate approximation error
	  maxdist_ = std::max(maxdist_, dist);
	  max_err = std::max(max_err, dist);
	  if (dist > aepsge_)
	    {
	      avdist_ += dist;
	      outsideeps_++;
	      av_err += dist;
	      outside++;
	    }
	}
      if (outside > 0)
	av_err /= (double)outside;

      // Store accuracy information in the element
      it->second->setAccuracyInfo(av_err, max_err, outside);
    }
  if (outsideeps_ > 0)
    avdist_ /= (double)outsideeps_;
}

//==============================================================================
void LRSurfApprox::refineSurf()
//==============================================================================
{
  // Prelimenary implementation. Write accuracy information to standard output
  // and get refinement informatiation
  // For each element, check if any data points are stored
  int idx=0;
  for (LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
       it != srf_->elementsEnd(); ++it)
    {
      if (it->second->hasAccuracyInfo())
	{
	  double av_err, max_err;
	  int nmb_out;
	  it->second->getAccuracyInfo(av_err, max_err, nmb_out);
	  std::cout << "Element " << idx << ", domain: [(" << it->second->umin() << ",";
	  std::cout << it->second->umax() << ")x(" << it->second->vmin() << ",";
	  std::cout << it->second->vmax() << ")]" << std::endl;
	  std::cout << "Average error: " << av_err << std::endl;
	  std::cout << "Maximume error: " << max_err << std::endl;
	  std::cout << "Number of points: " << it->second->nmbDataPoints();
	  std::cout << ", number of error points: " << nmb_out << std::endl;
	}
      ++idx;
    }

  // Specify refinement
  while (true)
    {
      int more;
      std::cout << "Refine more? " << std::endl;
      std::cin >> more;
      if (!more)
	break;
      double parval, start, end;
      int dir;
      int mult = 1;
      std::cout << "New knot, give pardir (0,1), fixed value, start and end: ";
      std::cout << std::endl;
      std::cin >> dir;
      std::cin >> parval;
      std::cin >> start;
      std::cin >> end;

      // The refinement must be performed inserting one knot segment at the
      // time in order to maintain information stored in the element
      srf_->refine((dir==0) ? YFIXED : XFIXED, parval, start, end, mult);

      // Update coef_known from information in LR B-splines
      updateCoefKnown();
    }
  
}

//==============================================================================
void LRSurfApprox::makeInitSurf(int dim)
//==============================================================================
{
  // Create a cubic Bezier surface represented as an LR B-spline surface
  // Compute domain
  double umin, umax, vmin, vmax;
  computeParDomain(dim, umin, umax, vmin, vmax);

  vector<double> knots_u(8);
  vector<double> knots_v(8);

  int kj;
  for (kj=0; kj<4; ++kj)
    {
      knots_u[kj] = umin;
      knots_v[kj] = vmin;
      knots_u[kj+4] = umax;
      knots_v[kj+4] = vmax;
    }

  makeInitSurf(dim, 4, 4, 4, 4, &knots_u[0], &knots_v[0]);
}

//==============================================================================
void LRSurfApprox::makeInitSurf(int dim, int ncoef_u, int order_u, int ncoef_v, 
				int order_v)
//==============================================================================
{
  // Compute domain
  double umin, umax, vmin, vmax;
  computeParDomain(dim, umin, umax, vmin, vmax);

  vector<double> knots_u(ncoef_u+order_u);
  vector<double> knots_v(ncoef_v+order_v);

  int kj;
  double del_u = (umax - umin)/(double)(ncoef_u-order_u+1);
  double del_v = (vmax - vmin)/(double)(ncoef_v-order_v+1);
  for (kj=0; kj<order_u; ++kj)
    knots_u[kj] = umin;
  for (kj=0; kj<order_v; ++kj)
    knots_v[kj] = vmin;
  for (; kj<ncoef_u; ++kj)
    knots_u[kj] = umin + (kj-order_u+1)*del_u;
  for (; kj<ncoef_v; ++kj)
    knots_v[kj] = vmin + (kj-order_v+1)*del_v;
  for (; kj<ncoef_u+order_u; ++kj)
    knots_u[kj] = umax;
  for (; kj<ncoef_v+order_v; ++kj)
    knots_v[kj] = vmax;

  makeInitSurf(dim, ncoef_u, order_u, ncoef_v, order_v, &knots_u[0], &knots_v[0]);
}

//==============================================================================
void LRSurfApprox::makeInitSurf(int dim, int ncoef_u, int order_u, int ncoef_v, 
				int order_v, double *knots_u, double *knots_v)
//==============================================================================
{
  vector<double> coefs(ncoef_u*ncoef_v*dim, 0.0);
  shared_ptr<SplineSurface> sf1(new SplineSurface(ncoef_u, ncoef_v, 
						  order_u, order_v,  
						  knots_u,     // ptr to knots in u
						  knots_v,     // ptr to knots in v
						  &coefs[0], // ptr to coefs
						  dim));

  // Reformatting data to the format that ApproxSurf wants (separating parameter values and 
  // data points in two distinct vectors).
  int del = 2 + dim;  // Number of entities for each point
  int nmb_pts = (int)points_.size()/del;
  vector<double> pts, param;
  pts.reserve(dim*nmb_pts);
  param.reserve(2 * nmb_pts);
  int ki;
  for (vector<double>::iterator it = points_.begin(); it != points_.end(); ) 
    {
      for (ki=0; ki<2; ++ki)
	param.push_back(*it++);
      for (ki=0; ki<dim; ++ki)
	pts.push_back(*it++);
    }
    
  // Approximate data points to create initial surface
  ApproxSurf asurf(sf1,
		   pts, //  scattered data points
		   param, // parameter values of scattered data points
		   dim, // dimension
		   aepsge_, // geometric tolerance
		   0, // 'constdir' - doesn't matter as we are not going to reparameterize anyway
		   false, // 'approx_orig' 
		   false, // 'close_belt' 
		   0,     // 'nmb_stabil' 
		   false); // 'repar' 

  asurf.setSmoothingWeight(smoothweight_);
  asurf.setDoRefine(false);
  asurf.setFixBoundary(false);
  int max_iter = 1;
  shared_ptr<SplineSurface> result_surf = asurf.getApproxSurf(maxdist_, avdist_, outsideeps_, 
							      max_iter);

  // Make LR spline surface
  double knot_tol = 1.0e-6;

  srf_ = shared_ptr<LRSplineSurface>(new LRSplineSurface(result_surf.get(), knot_tol));
  coef_known_.assign(srf_->numBasisFunctions(), 0.0);  // Initially nothing is fixed
}

//==============================================================================
void LRSurfApprox::makeInitSurf(shared_ptr<SplineSurface> surf)
//==============================================================================
{
  // Make LR spline surface
  double knot_tol = 1.0e-6;

  srf_ = shared_ptr<LRSplineSurface>(new LRSplineSurface(surf.get(), knot_tol));
  coef_known_.assign(srf_->numBasisFunctions(), 0.0);  // Initially nothing is fixed
}

//==============================================================================
void LRSurfApprox::computeParDomain(int dim, double& umin, double& umax, 
				    double& vmin, double& vmax)
//==============================================================================
{
  // Compute domain
  umin = umax = points_[0];
  vmin = vmax = points_[1];
  int del = 2 + dim;
  for (size_t ki=del; ki<points_.size(); ki+=del)
    {
      umin = std::min(umin, points_[ki]);
      umax = std::max(umax, points_[ki]);
      vmin = std::min(vmin, points_[ki+1]);
      vmax = std::max(vmax, points_[ki+1]);
    }
}

//==============================================================================
void LRSurfApprox::setCoefKnown()
//==============================================================================
{
  // Set coef_known from boundary information
  // Note that k-multiple knots at boundaries are expected, otherwise no
  // coefficients will be fixed
  
  for (size_t ki=0; ki<4; ++ki)
    {
      if (edge_derivs_[ki] == 0)
	continue;  // No coefficients to fix
      
      // Set boundary characteristica
      Direction2D d = (ki == 1 || ki == 3) ? XFIXED : YFIXED;  // Orthogonal to the curve
      Direction2D d2 = (ki == 1 || ki == 3) ? YFIXED : XFIXED; // Along the curve
      bool atstart = (ki == 0 || ki == 3);  // Whether the constant
      // parameter curve is in the start related to the opposite parameter direction
      // of the surface

      // Define coefficients as fixed in the associated b-splines
      int fixcoef = 1;
      setCoefKnown(d, d2, atstart, fixcoef);
    }

  // Transfer information to the global vector
  updateCoefKnown();
}

//==============================================================================
void LRSurfApprox::setCoefKnown(Direction2D d, Direction2D d2, 
				bool atstart, int fixcoef)
//==============================================================================
{
  // Check if the basis has k-tupple knots in the current direction
  const Mesh2D& mesh = srf_->mesh();
  int ix = (atstart) ? mesh.firstMeshVecIx(d) : 
    mesh.lastMeshVecIx(d);  
  int mult = mesh.minMultInLine(d, ix);
  int deg = srf_->degree(d);  // Degree orthogonal to the constant parameter curve
  if (mult < deg)
    return;  // No fixed coefficients are set

  // Fetch LRBSplines. Traverse the knot vector to make keys 
  // for the bspline map.
  int dir = d;
  int startmult[2], endmult[2];
  double startval[2], endval[2];
  int num = mesh.numDistinctKnots(d);  // Number of knots orthogonal to the
  // constant parameter curve

  // Set parameter value at the constant parameter curve
  if (atstart)
    {
      startval[dir] = mesh.kval(d, 0);
    }
  else
    {
      endval[dir] = mesh.kval(d, num-1);
    }
  int deg2 = srf_->degree(d2);  // Degree along the constant parameter curve

  vector<int> knot_idx =  LRBSpline2DUtils::derive_knots(mesh, d2, 
							 mesh.firstMeshVecIx(d2),
							 mesh.lastMeshVecIx(d2),
							 atstart ? ix : ix-1,
							 atstart ? ix+1 : ix);

  // Since we have multiple knots in the surface boundary and we are only interested
  // in the basis functions being non-zero along the boundary, we know that these
  // basis functions must have a multiplicity equal to the order at the boundary
  // and one in the other end.
  startmult[dir] = atstart ? deg+1 : 1;
  endmult[dir] = atstart ? 1 : deg+1;
  size_t k1, k2;
  vector<double> coefs;
  for (k1=0, k2=deg2+1; k2<knot_idx.size(); ++k1, ++k2)
    {
      // Fetch domain of first minimal LR B-spline along the constant
      // parameter curve
      // First orthogonal to the curve
      if (atstart)
	{
	  int k_idx = 
	    Mesh2DUtils::search_upwards_for_nonzero_multiplicity(mesh, d, 1, 
								 knot_idx[k1], knot_idx[k2]);
	  endval[dir] = mesh.kval(d, k_idx);
	}
      else
	{
	  int k_idx = 
	    Mesh2DUtils::search_downwards_for_nonzero_multiplicity(mesh, d, num-2, 
								   knot_idx[k1], knot_idx[k2]);
	  startval[dir] = mesh.kval(d, k_idx);
	}

      // Along the curve
      startval[1-dir] = mesh.kval(d2, knot_idx[k1]);
      endval[1-dir] = mesh.kval(d2, knot_idx[k2]);

      // Count multiplicitity along the curve
      int km = 1;
      for (; km<=deg2; ++km)
	if (knot_idx[k1+km] > knot_idx[k1])
	  break;
      startmult[1-dir] = km;

      km = 1;
      for (; km<=deg2; ++km)
	if (knot_idx[k2-km] < knot_idx[k2])
	  break;
      endmult[1-dir] = km;

      // Fetch the associated LR B-spline
      LRSplineSurface::BSplineMap::iterator bm = srf_->bsplineFromDomain(startval[0], startval[1], 
									 endval[0], endval[1], 
									 startmult[0], startmult[1], 
									 endmult[0], endmult[1]);

      // Set fixed coefficient information in the b-spline
      bm->second->setFixCoef(fixcoef);
    }
}

//==============================================================================
void LRSurfApprox::updateCoefKnown()
//==============================================================================
{
  coef_known_.resize(srf_->numBasisFunctions());

  LRSplineSurface::BSplineMap::const_iterator it_bs;			
  size_t ki;
  for (it_bs=srf_->basisFunctionsBegin(), ki=0; it_bs!=srf_->basisFunctionsEnd(); 
       ++it_bs, ++ki)
    coef_known_[ki] = it_bs->second->coefFixed();
}


