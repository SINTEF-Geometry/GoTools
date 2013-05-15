#include "GoTools/lrsplines2D/LRSurfSmoothLS.h"
#include "GoTools/lrsplines2D/Element2D.h"
#include "GoTools/creators/SolveCG.h"

using namespace Go;
using std::vector;

  //============================================================================
  // Constructs a map from a B-spline to a simple Index.
  // ============================================================================
  BsplineIndexMap construct_bsplineindex_map( const LRSplineSurface & lrs )
  {
    BsplineIndexMap result;
    size_t i = 0;
    for (LRSplineSurface::BSplineMap::const_iterator it_bs=lrs.basisFunctionsBegin(); 
        it_bs!=lrs.basisFunctionsEnd(); ++it_bs)
      {
	if (it_bs->second->coefFixed())
	  continue;
	
	result[it_bs->second.get()] = i++;
      }
    return result;
  }


//==============================================================================
LRSurfSmoothLS::LRSurfSmoothLS(shared_ptr<LRSplineSurface> surf, vector<int>& coef_known)
//==============================================================================
  : srf_(surf), coef_known_(coef_known)
{
  // Distribute information about fixed coefficients to the B-splines
  ncond_ = 0;
  LRSplineSurface::BSplineMap::const_iterator it_bs;
  size_t ki;
  for (it_bs=srf_->basisFunctionsBegin(), ki=0; it_bs!=srf_->basisFunctionsEnd(); 
       ++it_bs, ++ki)
    {
      it_bs->second->setFixCoef(coef_known[ki]);
      if (coef_known[ki] == 0)
	ncond_++;
    }
  
  // Construct index map
  BSmap_ = construct_bsplineindex_map(*srf_);

  // Allocate scratch for equation system
  gmat_.assign(ncond_*ncond_, 0.0);
  gright_.assign(srf_->dimension()*ncond_, 0.0);
  
}

//==============================================================================
LRSurfSmoothLS::~LRSurfSmoothLS()
//==============================================================================
{
}

//==============================================================================
void LRSurfSmoothLS::updateLocals()
//==============================================================================
{
  ncond_ = srf_->numBasisFunctions();
  BSmap_ = construct_bsplineindex_map(*srf_);

  gmat_.assign(ncond_*ncond_, 0.0);
  gright_.assign(srf_->dimension()*ncond_, 0.0);
}

//==============================================================================
bool LRSurfSmoothLS::hasDataPoints() const
//==============================================================================
{
  // For each element, check if any data points are stored
  for (LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
       it != srf_->elementsEnd(); ++it)
    {
      if (it->second->hasDataPoints())
	return true;   // Data points are found associated to the current element
    }
  return false;  // No data points are found
}

int compare_u_par(const void* el1, const void* el2)
{
  if (((double*)el1)[0] < ((double*)el2)[0])
    return -1;
  else if (((double*)el1)[0] > ((double*)el2)[0])
    return 1;
  else
    return 0;
}

int compare_v_par(const void* el1, const void* el2)
{
  if (((double*)el1)[1] < ((double*)el2)[1])
    return -1;
  else if (((double*)el1)[1] > ((double*)el2)[1])
    return 1;
  else
    return 0;
}

//==============================================================================
void LRSurfSmoothLS::addDataPoints(vector<double>& points) 
//==============================================================================
{
  int dim = srf_->dimension();
  int del = dim+2;                   // Number of entries for each point
  int nmb = (int)points.size()/del;  // Number of data points

  // Erase point information in the elements
  for (LRSplineSurface::ElementMap::const_iterator it = srf_->elementsBegin();
       it != srf_->elementsEnd(); ++it)
    it->second->eraseDataPoints();

  // Sort the points according to the u-parameter
  qsort(&points[0], nmb, del*sizeof(double), compare_u_par);

  // Get all knot values in the u-direction
  const double* const uknots_begin = srf_->mesh().knotsBegin(XFIXED);
  const double* const uknots_end = srf_->mesh().knotsEnd(XFIXED);
  const double* knotu;

  // Get all knot values in the v-direction
  const double* const vknots_begin = srf_->mesh().knotsBegin(YFIXED);
  const double* const vknots_end = srf_->mesh().knotsEnd(YFIXED);
  const double* knotv;

  // Traverse points and divide them according to their position in the
  // u direction
  int pp0, pp1;
  for (pp0=0, knotu=uknots_begin, ++knotu; knotu!= uknots_end; ++knotu)
    {
      
      for (pp1=pp0; pp1<(int)points.size() && points[pp1] < (*knotu); pp1+=del);
      if (knotu+1 == uknots_end)
	pp1 = (int)points.size();

      // Sort the current sub set of points according to the v-parameter
      qsort(&points[0]+pp0, (pp1-pp0)/del, del*sizeof(double), compare_v_par);

      // Traverse the relevant points and store them in the associated element
      int pp2, pp3;
      for (pp2=pp0, knotv=vknots_begin, ++knotv; knotv!=vknots_end; ++knotv)
	{
	  for (pp3=pp2; pp3<pp1 && points[pp3+1] < (*knotv); pp3 += del);
	  if (knotv+1 == vknots_end)
	    pp3 = pp1;
	  
	  // Fetch associated element
	   Element2D* elem = 
	     const_cast<Element2D*>(srf_->coveringElement(0.5*(knotu[-1]+knotu[0]), 
							  0.5*(knotv[-1]+knotv[0])));
	  elem->addDataPoints(points.begin()+pp2, points.begin()+pp3);

	  pp2 = pp3;
	}
      pp0 = pp1;
    }
}

//==============================================================================
void LRSurfSmoothLS::setOptimize(const double weight1, const double weight2,
				 const double weight3)
//==============================================================================
{
  // Not implemented yet
}

//==============================================================================
void LRSurfSmoothLS::setLeastSquares(const double weight)
//==============================================================================
{
  int dim = srf_->dimension();

  // For each element
  for (LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
       it != srf_->elementsEnd(); ++it)
    {
      // Check if the element contains an associated least squares matrix
      bool has_LS_mat = it->second->hasLSMatrix();

      // Check if the element is changed
      bool is_modified = it->second->isModified();

      // Fetch B-splines
      const vector<LRBSpline2D*>& bsplines = it->second->getSupport();
      size_t nmb = bsplines.size();

      if (!has_LS_mat || is_modified)
	{
	  // Either no pre-computed least squares matrix exists or 
	  // the element or an associated B-spline is changed.
	  // Compute the least squares matrix associated to the 
	  // element
	  // First fetch data points
	  vector<double>& elem_data = it->second->getDataPoints();

	  // Compute sub matrix
	  // First get access to storage in the element
	  double *subLSmat, *subLSright;
	  int kcond;
 	  it->second->setLSMatrix();
	  it->second->getLSMatrix(subLSmat, subLSright, kcond);
 
	  localLeastSquares(elem_data, bsplines, subLSmat, subLSright, kcond);
	  int stop_break;
	}

      // Assemble stiffness matrix and right hand side based on the local least 
      // squares matrix
      // The size of the stiffness matrix is the squared number of LR B-splines
      // with a free coefficient. The size of the right hand side is equal to
      // the number of free coefficients times the dimension of the data points
      double *subLSmat, *subLSright;
      int kcond;
      it->second->getLSMatrix(subLSmat, subLSright, kcond);

      vector<size_t> in_bs(kcond);
      size_t ki, kj, kr, kh;
      for (ki=0, kj=0; ki<nmb; ++ki)
	{
	  if (bsplines[ki]->coefFixed())
	    continue;

	  // Fetch index in the stiffness matrix
	  size_t inb = BSmap_.at(bsplines[kj]);
	  in_bs[kj++] = inb;
	}

      int kk;
      for (ki=0, kr=0; ki<nmb; ++ki)
	{
	  if (bsplines[ki]->coefFixed())
	    continue;
	  size_t inb1 = in_bs[ki];
	  for (kk=0; kk<dim; ++kk)
	    gright_[kk*ncond_+inb1] += weight*subLSright[ki];
	  for (kj=0, kh=0; kj<nmb; ++kj)
	    {
	      if (bsplines[kj]->coefFixed())
		continue;
	      gmat_[inb1*ncond_+in_bs[kh]] += weight*subLSmat[kr*nmb+kh];
	      kh++;
	    }
	  kr++;
	}
    }

}

//==============================================================================
void LRSurfSmoothLS::setLeastSquares(vector<double>& points, const double weight)
//==============================================================================
{
  addDataPoints(points);
  setLeastSquares(weight);
}

//==============================================================================
int
LRSurfSmoothLS::equationSolve(shared_ptr<LRSplineSurface>& surf)
//==============================================================================
{
  int kstat = 0;
  int dim = srf_->dimension();
  int ki, kk;

  // Solve the equation system by Conjugate Gradient Method.

  vector<double> eb(dim*ncond_, 0.0);
  for (ki=0; ki<dim*ncond_; ki++)
    eb[ki] = gright_[ki];

  // Copy coefficients to array of unknowns
  LRSplineSurface::BSplineMap::const_iterator it_bs;
  for (it_bs=srf_->basisFunctionsBegin(), ki=0; 
       it_bs!=srf_->basisFunctionsEnd(); ++it_bs)
    {
      if (it_bs->second->coefFixed())
	continue;
      
      Point cf = it_bs->second->Coef();
      for (kk=0; kk<dim; kk++)
	gright_[kk*ncond_+ki] = cf[kk]; 
      ki++;
    }
       
  // Set up CG-object

  SolveCG solveCg;

  // Create sparse matrix.

  ASSERT(gmat_.size() > 0);
  solveCg.attachMatrix(&gmat_[0], ncond_);

  // Attach parameters.

  solveCg.setTolerance(0.00000001);
  // Preconditioning
  // @@sbr When using the side-constraints our system is not longer guaranteed
  //       to be symmetric and positive definite!
  //       We should then use a different solver.
  int precond = 1;
  int nmb_iter = precond ? ncond_ : 2*ncond_;
  solveCg.setMaxIterations(std::min(nmb_iter, 1000));
  if (precond) {
    double omega = 0.1;
    // 	   printf("Omega = ");
    // 	   scanf("%lf",&omega);
    solveCg.precondRILU(omega);
  }

  // Solve equation systems.
       
  for (kk=0; kk<dim; kk++)
    {
      kstat = solveCg.solve(&gright_[kk*ncond_], &eb[kk*ncond_],
			    ncond_);
      //	       printf("solveCg.solve status %d \n", kstat);
      if (kstat < 0)
	return kstat;
      if (kstat == 1)
	THROW("Failed solving system (within tolerance)!");
    }

  // Update coefficients
  for (it_bs=srf_->basisFunctionsBegin(), ki=0; 
       it_bs!=srf_->basisFunctionsEnd(); ++it_bs)
    {
      if (it_bs->second->coefFixed())
	continue;
      
      Point cf(dim);
      for (kk=0; kk<dim; kk++)
	cf[kk] = gright_[kk*ncond_+ki];
      srf_->setCoef(cf, it_bs->second.get());
      ki++;
    }


  // Output surface
  surf = srf_;

  return 0;
}

//==============================================================================
void LRSurfSmoothLS::localLeastSquares(vector<double>& points,
				       const vector<LRBSpline2D*>& bsplines,
				       double* mat, double* right, int ncond)
//==============================================================================
{
  size_t nmbb = bsplines.size();
  int dim = srf_->dimension();
  int del = dim+2;
  int nmbp = (int)points.size()/del;

  size_t ki, kj, kp, kq, kr, kk;
  double *pp;
  for (kr=0, pp=&points[0]; kr<nmbp; ++kr, pp+=del)
    {
      vector<double> sb = getBasisValues(bsplines, pp);
      for (ki=0, kj=0; ki<nmbb; ++ki)
	{
	  if (bsplines[ki]->coefFixed())
	    continue;
	  for (kk=0; kk<dim; ++kk)
	    right[kk*ncond+kj] += pp[2+kk]*sb[kj];
	  for (kp=0, kq=0; kp<nmbb; kp++)
	    {
	      int fixed = bsplines[kp]->coefFixed();
	      if (fixed == 2)
		  continue;

	      double val = sb[kj]*sb[kq];
	      if (fixed == 1)
		{
		  // Move contribution to the right hand side
		  const Point coef = bsplines[kp]->Coef();
		  for (kk=0; kk<dim; ++kk)
		    right[kk*ncond+kj] -= coef[kk]*val;
		}
	      else
		{
		  mat[kp*ncond+kj] += val;
		}
	      kq++;
	    }
	  kj++;
	}
    }
}

//==============================================================================
vector<double> LRSurfSmoothLS::getBasisValues(const vector<LRBSpline2D*>& bsplines,
					      double *par)
//==============================================================================
{
  vector<double> bs(bsplines.size());
  for (size_t ki=0; ki<bsplines.size(); ++ki)
    {
      const bool u_on_end = (par[0] == bsplines[ki]->umax());
      const bool v_on_end = (par[1] == bsplines[ki]->vmax());
      bs[ki] = bsplines[ki]->evalBasisFunction(par[0], par[1], 0, 0, 
					       u_on_end, v_on_end);
    }
  return bs;
}
