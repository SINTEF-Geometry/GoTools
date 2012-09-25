//===========================================================================
//                                                                           
// File: VolBoundaryCondition.C                                              
//                                                                           
// Created: Tue Nov  1 10:17:18 2011                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/isogeometric_model/VolBoundaryCondition.h"
#include "GoTools/isogeometric_model/VolSolution.h"
#include "GoTools/isogeometric_model/EvalFunctorSurface.h"
#if 0
#include "GoTools/compositemodel/AdaptEvalSurface.h"
#endif
#include <assert.h>


using std::pair;

namespace Go
{

  //===========================================================================
  VolBoundaryCondition::VolBoundaryCondition(int face_nmb, BdConditionType type,
					     BdCondFunctor *fbd,
					     std::vector<std::pair<double, double> >& domain,
					     VolSolution *solution)
    : BlockBoundaryCondition(type),
      parent_(solution),
      facenmb_(face_nmb),
      fbd_(fbd),
      domain_(domain),
      approx_err_(-1.0)
  //===========================================================================
  {
  }


  //===========================================================================
  VolBoundaryCondition::VolBoundaryCondition(int face_nmb, BdConditionType type,
					     const Point& const_val,
					     std::vector<std::pair<double, double> >& domain,
					     VolSolution *solution)
    : BlockBoundaryCondition(type),
      parent_(solution),
      facenmb_(face_nmb),
      fbd_(0),
      const_val_(const_val),
      domain_(domain),
      approx_err_(-1.0)
  //===========================================================================
  {
  }

  //===========================================================================
  VolBoundaryCondition::~VolBoundaryCondition()
  //===========================================================================
  {
  }

  //===========================================================================
  void 
  VolBoundaryCondition::getCoefficientsEnumeration(std::vector<int>& local_enumeration)
  //===========================================================================
  {
    vector<int> full_enumeration;
    parent_->getBoundaryCoefficients(facenmb_, full_enumeration);
    int dir1 = (facenmb_ < 2) ? 1 : 0;
    int dir2 = (facenmb_ < 2) ? 2 : 1;
    BsplineBasis bas1 = parent_->basis(dir1);
    BsplineBasis bas2 = parent_->basis(dir2);
    // @@sbr Can we expect the surfaces to have parametrization (u,v) or (u,w) etc? Or (w,u)?
    // I.e. which basises are the values in domain_ referring to?
    double min1 = domain_[0].first;
    double max1 = domain_[0].first;
    double min2 = domain_[0].second;
    double max2 = domain_[0].second;
    for (size_t ki = 1; ki < domain_.size() - 1; ++ki)
    {
	if (domain_[ki].first < min1)
	    min1 = domain_[ki].first;
	if (max1 < domain_[ki].first)
	    max1 = domain_[ki].first;

	if (domain_[ki].second < min2)
	    min2 = domain_[ki].second;
	if (max2 < domain_[ki].second)
	    max2 = domain_[ki].second;
    }

    int pos1_1 = bas1.knotIntervalFuzzy(min1);
    pos1_1 -= bas1.order() - 1;
    int pos1_2 = bas1.knotIntervalFuzzy(max1);
    int num_coef_1 = pos1_2 - pos1_1 + 1;

    int pos2_1 = bas2.knotIntervalFuzzy(min2);
    pos2_1 -= bas2.order() - 1;
    int pos2_2 = bas2.knotIntervalFuzzy(max2);
    int num_coef_2 = pos2_2 - pos2_1 + 1;

    int in1 = bas1.numCoefs();
    local_enumeration.resize(num_coef_1*num_coef_2);
    int cntr = 0;
    for (int kj = 0; kj < num_coef_2; ++kj)
    {
	int c2 = pos2_1 + kj;
	for (int ki = 0; ki < num_coef_1; ++ki)
	{
	    int c1 = pos1_1 + ki;
	    int ci = c2*in1 + c1;
	    local_enumeration[cntr] = full_enumeration[ci];
	    ++cntr;
	}
    }
 }


  //===========================================================================
  void 
  VolBoundaryCondition::getCoefficientsEnumeration(std::vector<int>& local_enumeration_bd,
						   std::vector<int>& local_enumeration_bd2)
  //===========================================================================
  {
    vector<int> full_enumeration_bd, full_enumeration_bd2;
    parent_->getBoundaryCoefficients(facenmb_,
				     full_enumeration_bd, full_enumeration_bd2);
    int dir1 = (facenmb_ < 2) ? 1 : 0;
    int dir2 = (facenmb_ < 2) ? 2 : 1;
    BsplineBasis bas1 = parent_->basis(dir1);
    BsplineBasis bas2 = parent_->basis(dir2);
    // @@sbr Can we expect the surfaces to have parametrization (u,v) or (u,w) etc? Or (w,u)?
    // I.e. which basises are the values in domain_ referring to? Guessing the former.
    double min1 = domain_[0].first;
    double max1 = domain_[0].first;
    double min2 = domain_[0].second;
    double max2 = domain_[0].second;
    for (size_t ki = 1; ki < domain_.size() - 1; ++ki)
    {
	if (domain_[ki].first < min1)
	    min1 = domain_[ki].first;
	if (max1 < domain_[ki].first)
	    max1 = domain_[ki].first;

	if (domain_[ki].second < min2)
	    min2 = domain_[ki].second;
	if (max2 < domain_[ki].second)
	    max2 = domain_[ki].second;
    }

    int pos1_1 = bas1.knotIntervalFuzzy(min1);
    pos1_1 -= bas1.order() - 1;
    int pos1_2 = bas1.knotIntervalFuzzy(max1);
    int num_coef_1 = pos1_2 - pos1_1 + 1;

    int pos2_1 = bas2.knotIntervalFuzzy(min2);
    pos2_1 -= bas2.order() - 1;
    int pos2_2 = bas2.knotIntervalFuzzy(max2);
    int num_coef_2 = pos2_2 - pos2_1 + 1;

    int in1 = bas1.numCoefs();
    local_enumeration_bd.resize(num_coef_1*num_coef_2);
    local_enumeration_bd2.resize(num_coef_1*num_coef_2);
    int cntr = 0;
    for (int kj = 0; kj < num_coef_2; ++kj)
    {
	int c2 = pos2_1 + kj;
	for (int ki = 0; ki < num_coef_1; ++ki)
	{
	    int c1 = pos1_1 + ki;
	    int ci = c2*in1 + c1;
	    local_enumeration_bd[cntr] = full_enumeration_bd[ci];
	    local_enumeration_bd2[cntr] = full_enumeration_bd2[ci];
	    ++cntr;
	}
    }

  }

  //===========================================================================
  void
  VolBoundaryCondition::getBdCoefficients(vector<pair<int, Point> >& coefs)
  //===========================================================================
  {
    if (!isDirichlet())
      return;

    const shared_ptr<SplineVolume> vol = parent_->getSolutionVolume();
    vector<int> coefs_enum;
    getCoefficientsEnumeration(coefs_enum);
    int coefs_size = (int)coefs_enum.size();
    int dim = vol->dimension();
    bool rational = vol->rational();
    int kdim = dim + (rational ? 1 : 0);
    coefs.resize(coefs_size);

    for (int i = 0; i < coefs_size; ++i)
      {
	Point p(dim);
	vector<double>::const_iterator it = vol->ctrl_begin() + coefs_enum[i] * kdim;
	if (rational)
	  for (int j = 0; j < dim; ++j)
	    p[j] = it[j] / it[dim];
	else
	  for (int j = 0; j < dim; ++j)
	    p[j] = it[j];
	coefs[i] = pair<int, Point>(coefs_enum[i], p);
      }
  }


  //===========================================================================
  void
  VolBoundaryCondition::getBdCoefficients(std::vector<std::pair<int, Point> >& coefs_bd,
					  std::vector<std::pair<int, Point> >& coefs_bd2)
  //===========================================================================
  {
    if (!isDirichlet())
      return;

    const shared_ptr<SplineVolume> vol = parent_->getSolutionVolume();
    vector<int> coefs_enum_bd, coefs_enum_bd2;
    getCoefficientsEnumeration(coefs_enum_bd, coefs_enum_bd2);
    int coefs_size = (int)coefs_enum_bd.size();
    int dim = vol->dimension();
    bool rational = vol->rational();
    int kdim = dim + (rational ? 1 : 0);
    coefs_bd.resize(coefs_size);
    coefs_bd2.resize(coefs_size);

    for (int i = 0; i < coefs_size; ++i)
      {
	Point p_bd(dim), p_bd2(dim);
	vector<double>::const_iterator it_bd = vol->ctrl_begin() + coefs_enum_bd[i] * kdim;
	vector<double>::const_iterator it_bd2 = vol->ctrl_begin() + coefs_enum_bd2[i] * kdim;
	if (rational)
	  for (int j = 0; j < dim; ++j)
	    {
	      p_bd[j] = it_bd[j] / it_bd[dim];
	      p_bd2[j] = it_bd2[j] / it_bd2[dim];
	    }
	else
	  for (int j = 0; j < dim; ++j)
	    {
	      p_bd[j] = it_bd[j];
	      p_bd2[j] = it_bd2[j];
	    }
	coefs_bd[i] = pair<int, Point>(coefs_enum_bd[i], p_bd);
	coefs_bd2[i] = pair<int, Point>(coefs_enum_bd2[i], p_bd2);
      }
  }


  //===========================================================================
  void VolBoundaryCondition::update()
  //===========================================================================
  {
#ifndef NDEBUG
    MESSAGE("update() under construction");
#endif

    if (!isDirichlet())
      return;    // Only Dirichlet conditions might update the approximation surface.

    double tol = getTolerances().gap;
    if (approx_err_ >= 0.0 && approx_err_ < tol)
      return;    // No update is done if approximation within tolerance has allready occured

    // The domain is described by a polygon in the parameter domain.
#if 1

    assert(domain_.size() > 0);

    // For now we assume that the domain is rectangular and axis-aligned.
#if 1
    double umin = domain_[0].first;
    double umax = umin;
    double vmin = domain_[0].second;
    double vmax = vmin;
    for (size_t ki = 1; ki < domain_.size(); ++ki)
    {
	// First the u-par.
	if (domain_[ki].first < umin)
	    umin = domain_[ki].first;
	if (domain_[ki].first > umax)
	    umax = domain_[ki].first;
	// Then the v-par.
	if (domain_[ki].second < vmin)
	    vmin = domain_[ki].second;
	if (domain_[ki].second > vmax)
	    vmax = domain_[ki].second;

    }
    // @@sbr201209 We should test that domain is indeed rectangular.
#endif

    // @@sbr201209 It seems that the approximated surface is aligned
    // with the direction of the volume, hence no need to flip or
    // reverse.
    int orientation = 1;
    bool same_dir = true;

    const shared_ptr<SplineVolume> sol_vol = parent_->getSolutionVolume();

    shared_ptr<SplineSurface> face_srf = shared_ptr<SplineSurface>(sol_vol->getBoundarySurface(facenmb_));
    shared_ptr<SplineSurface> bas_srf = shared_ptr<SplineSurface>(face_srf->subSurface(umin, vmin, umax, vmax));
    int ncoefs = bas_srf->numCoefs_u()*bas_srf->numCoefs_v();
    int ncoefs_u = bas_srf->numCoefs_u();
    int ncoefs_v = bas_srf->numCoefs_v();
    int order_u = bas_srf->order_u();
    int order_v = bas_srf->order_v();
    vector<double> knots_u(ncoefs_u + order_u);
    vector<double> knots_v(ncoefs_v + order_v);
    copy(bas_srf->basis_u().begin(), bas_srf->basis_u().end(), knots_u.begin());
    copy(bas_srf->basis_v().begin(), bas_srf->basis_v().end(), knots_v.begin());

    if (getBdConditionType() == DIRICHLET)
      {
	shared_ptr<SplineSurface> geo_face_srf = shared_ptr<SplineSurface>
	  (parent_->getGeometryVolume()->getBoundarySurface(facenmb_));
	shared_ptr<SplineSurface> geo_bas_srf = shared_ptr<SplineSurface>(geo_face_srf->subSurface(umin, vmin, umax, vmax));
	shared_ptr<EvalFunctorSurface> efc(new EvalFunctorSurface(fbd_, geo_bas_srf, bas_srf->dimension()));

#if 1
//	MESSAGE("Missing call to AdaptEvalSurface!");
	// @@sbr201209 We create an approximating surface which is
	// 0.0, as that is the current case. Not sure if we will need
	// to approximate this surface with non-zero values.  Since
	// our surface is zero a linear space is sufficient.
	vector<double> knots_u_appr(4), knots_v_appr(4);
	knots_u_appr[0] = knots_u_appr[1] = umin;
	knots_u_appr[2] = knots_u_appr[3] = umax;
	knots_v_appr[0] = knots_v_appr[1] = vmin;
	knots_v_appr[2] = knots_v_appr[3] = vmax;
	int order_u_appr = 2;
	int order_v_appr = 2;
	int num_coefs_u_appr = 2;
	int num_coefs_v_appr = 2;
	int dim = bas_srf->dimension();
	vector<double> coefs_appr(num_coefs_u_appr*num_coefs_v_appr*dim, 0.0);
	bdsrf_cond_ = shared_ptr<SplineSurface>(new SplineSurface(num_coefs_u_appr, num_coefs_v_appr,
								  order_u_appr, order_v_appr,
								  knots_u_appr.begin(), knots_v_appr.begin(),
								  coefs_appr.begin(), dim));
	

#else
	shared_ptr<AdaptEvalSurface> adap_srf;
	adap_srf =
	    shared_ptr<AdaptEvalSurface>(new AdaptEvalSurface(efc.get(), tol,
							      ncoefs_u, ncoefs_v, order_u, order_v,
							      knots_u, knots_v));

	// adap_srf->approximate(1);
	const int max_iter = 1;
	adap_srf->doApprox(geo_bas_srf, max_iter, );
    // shared_ptr<SplineSurface>
    //   doApprox(shared_ptr<SplineSurface> init_surf, int max_iter,
    // 	       shared_ptr<ftPointSet> points, double tol,
    // 	       double& max_error, double& mean_error);

			   double maxdist, avdist;
	bdsrf_cond_ = adap_srf->getAdaptSurface(maxdist, avdist);
	approx_err_ = avdist;
#endif

      }
    else
      {
	int dim = bas_srf->dimension();
	int ncoefs = ncoefs_u*ncoefs_v;
	vector<double> coefs(dim*ncoefs);
	for (int i = 0, pos = 0; i < ncoefs; ++i)
	  for (int j = 0; j < dim; ++j, ++pos)
	    coefs[pos] = const_val_[j];
	bdsrf_cond_ = shared_ptr<SplineSurface>(new SplineSurface(ncoefs_u, ncoefs_v, order_u, order_v,
								  knots_u.begin(), knots_v.begin(),
								  coefs.begin(), dim));
	approx_err_ = 0.0;
      }

    // Update coefficients according to boundary surface
    double start_par_u = knots_u[order_u - 1];
    int edge_srf_start_u;
    for (edge_srf_start_u = 0; edge_srf_start_u < (int)knots_u.size() && knots_u[edge_srf_start_u] < start_par_u; ++edge_srf_start_u);
    if (edge_srf_start_u == (int)knots_u.size())
      return;
    int cond_srf_start_u = order_u - face_srf->basis_u().knotMultiplicity(knots_u[edge_srf_start_u]);

    double start_par_v = knots_v[order_v - 1];
    int edge_srf_start_v;
    for (edge_srf_start_v = 0; edge_srf_start_v < (int)knots_v.size() && knots_v[edge_srf_start_v] < start_par_v; ++edge_srf_start_v);
    if (edge_srf_start_v == (int)knots_v.size())
      return;
    int cond_srf_start_v = order_v - face_srf->basis_v().knotMultiplicity(knots_v[edge_srf_start_v]);

    vector<int> coefs_enum;
    getCoefficientsEnumeration(coefs_enum);
    int coefs_size = (int)coefs_enum.size();
    int dim = sol_vol->dimension();
    bool rational = sol_vol->rational();
    int kdim = dim + (rational ? 1 : 0);

    int cond_srf_start = cond_srf_start_v*ncoefs_u + cond_srf_start_u;

    bool srf_rational = bdsrf_cond_->rational();
    int srf_kdim = dim + (srf_rational ? 1 : 0);
    vector<double>::const_iterator srf_it = srf_rational ? bdsrf_cond_->rcoefs_begin() : bdsrf_cond_->coefs_begin();
    srf_it += cond_srf_start * srf_kdim;

    int srf_it_pos = same_dir ? 0 : (coefs_size - 1) * kdim;
    for (int i = 0; i < coefs_size; ++i)
      {
	vector<double>::iterator surf_it = sol_vol->ctrl_begin() + coefs_enum[i] * kdim;
	for (int j = 0; j < dim; ++j)
	  surf_it[j] = srf_it[srf_it_pos + j];
	if (same_dir)
	  srf_it_pos += srf_kdim;
	else
	  srf_it_pos -= srf_kdim;
      }
#endif
  }

  //===========================================================================
  void VolBoundaryCondition::getBasisFunctions(int index_of_Gauss_point1,
					       int index_of_Gauss_point2,
					       int& const_dir,
					       vector<double>& basisValues,
					       vector<double>& basisDerivs_u,
					       vector<double>& basisDerivs_v) const
  //===========================================================================
  {
    double par1, par2;
    parent_->getGaussParameter(index_of_Gauss_point1, index_of_Gauss_point2,
			       const_dir, par1, par2);

    double param[2];
    param[0] = par1;
    param[1] = par2;
    shared_ptr<SplineSurface> bd_surf =
	shared_ptr<SplineSurface>(parent_->getSolutionVolume()->getBoundarySurface(facenmb_));
    bd_surf->computeBasis(param, basisValues, basisDerivs_u, basisDerivs_v);

    // shared_ptr<SplineSurface> bd_srf = parent_->getSolutionVolume()->getBoundarySurface(facenmb_);
    // bd_srf->computeBasis(param, basisValues, basisDerivs);
  }

  //===========================================================================
  shared_ptr<SplineSurface> VolBoundaryCondition::getSplineApproximation() const
  //===========================================================================
  {
    return bdsrf_cond_;
  }

  //===========================================================================
  void VolBoundaryCondition::updateBoundaryValue(BdCondFunctor* fbd)
  //===========================================================================
  {
    MESSAGE("updateBoundaryValue() not implemented");
  }

  //===========================================================================
  int VolBoundaryCondition::faceNumber() const
  //===========================================================================
  {
    return facenmb_;
  }

  //===========================================================================
  tpTolerances VolBoundaryCondition::getTolerances() const
  //===========================================================================
  {
      return parent_->getTolerances();
  }

}

