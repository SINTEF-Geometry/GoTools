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
    MESSAGE("update() under construction");


    if (!isDirichlet())
      return;    // Only Dirichlet conditions might update the approximation surface.

    double tol = getTolerances().gap;
    if (approx_err_ >= 0.0 && approx_err_ < tol)
      return;    // No update is done if approximation within tolerance has allready occured

    // The domain is described by a polygon in the parameter domain.
#if 0
    bool same_dir = domain_.first <= domain_.second;
    double start_par = same_dir ? domain_.first : domain_.second;
    double end_par = same_dir ? domain_.second : domain_.first;
    const shared_ptr<SplineVolume> sol_vol = parent_->getSolutionVolume();

    shared_ptr<SplineSurface> face_srf = shared_ptr<SplineSurface>(sol_vol->getBoundarySurface(facenmb_));
    shared_ptr<SplineSurface> bas_srf = shared_ptr<SplineSurface>(face_srf->subCurve(start_par, end_par));
    int ncoefs = bas_curve->numCoefs();
    int order = bas_curve->order();
    vector<double> knots(ncoefs + order);
    copy(bas_curve->knotsBegin(), bas_curve->knotsEnd(), knots.begin());

    if (getBdConditionType() == DIRICHLET)
      {
	shared_ptr<SplineCurve> geo_edge_curve = shared_ptr<SplineCurve>(parent_->getGeometrySurface()->edgeCurve(ccw_edge_number));
	shared_ptr<SplineCurve> geo_bas_curve = shared_ptr<SplineCurve>(geo_edge_curve->subCurve(start_par, end_par));
	shared_ptr<EvalFunctorCurve> efc(new EvalFunctorCurve(fbd_, geo_bas_curve, bas_curve->dimension()));

	shared_ptr<AdaptCurve> adap_crv = shared_ptr<AdaptCurve>(new AdaptCurve(efc.get(), tol, ncoefs, order, knots));
	adap_crv->approximate(1);
	double maxdist, avdist;
	bdcrv_cond_ = adap_crv->getAdaptCurve(maxdist, avdist);
	approx_err_ = avdist;
      }
    else
      {
	int dim = bas_curve->dimension();
	vector<double> coefs(dim * ncoefs);
	for (int i = 0, pos = 0; i < ncoefs; ++i)
	  for (int j = 0; j < dim; ++j, ++pos)
	    coefs[pos] = const_val_[j];
	bdcrv_cond_ = shared_ptr<SplineCurve>(new SplineCurve(ncoefs, order, knots.begin(), coefs.begin(), dim));
	approx_err_ = 0.0;
      }

    // Update coefficients according to boundary curve
    int edge_crv_start;
    for (edge_crv_start = 0; edge_crv_start < (int)knots.size() && knots[edge_crv_start] < start_par; ++edge_crv_start);
    if (edge_crv_start == (int)knots.size())
      return;
    int cond_crv_start = order - edge_curve->basis().knotMultiplicity(knots[edge_crv_start]);
    vector<int> local_enumeration;

    vector<int> coefs_enum;
    getCoefficientsEnumeration(coefs_enum);
    int coefs_size = (int)coefs_enum.size();
    int dim = sol_surf->dimension();
    bool rational = sol_surf->rational();
    int kdim = dim + (rational ? 1 : 0);

    bool crv_rational = bdcrv_cond_->rational();
    int crv_kdim = dim + (crv_rational ? 1 : 0);
    vector<double>::const_iterator crv_it = crv_rational ? bdcrv_cond_->rcoefs_begin() : bdcrv_cond_->coefs_begin();
    crv_it += cond_crv_start * crv_kdim;
    int crv_it_pos = same_dir ? 0 : (coefs_size - 1) * kdim;
    for (int i = 0; i < coefs_size; ++i)
      {
	vector<double>::iterator surf_it = sol_surf->ctrl_begin() + coefs_enum[i] * kdim;
	for (int j = 0; j < dim; ++j)
	  surf_it[j] = crv_it[crv_it_pos + j];
	if (same_dir)
	  crv_it_pos += crv_kdim;
	else
	  crv_it_pos -= crv_kdim;
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

