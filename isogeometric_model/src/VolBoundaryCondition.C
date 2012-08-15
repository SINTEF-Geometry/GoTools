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
    int dir1 = (facenmb_ << 2) ? 1 : 0;
    int dir2 = (facenmb_ == 2 || facenmb_ == 3) ? 2 : 1;
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
    int dir1 = (facenmb_ << 2) ? 1 : 0;
    int dir2 = (facenmb_ == 2 || facenmb_ == 3) ? 2 : 1;
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
    MESSAGE("getBdCoefficients() under construction");

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
    MESSAGE("update() not implemented");
  }

  //===========================================================================
  void VolBoundaryCondition::getBasisFunctions(int index_of_Gauss_point1,
					       int index_of_Gauss_point2,
					       shared_ptr<BasisDerivs> result,
					       int solutionspace_idx) const
  //===========================================================================
  {
    MESSAGE("getBasisFunctions() not implemented");

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

