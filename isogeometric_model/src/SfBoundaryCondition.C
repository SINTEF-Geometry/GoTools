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

#include "GoTools/geometry/PointCloud.h"
#include "GoTools/isogeometric_model/SfBoundaryCondition.h"
#include "GoTools/creators/AdaptCurve.h"
#include <algorithm>
#include "GoTools/isogeometric_model/SfSolution.h"
#include "GoTools/isogeometric_model/EvalFunctorCurve.h"
#include <assert.h>

using std::vector;
using std::pair;

namespace Go
{

  //===========================================================================
  SfBoundaryCondition::SfBoundaryCondition(int edge_nmb, BdConditionType type, BdCondFunctor *fbd,
					   pair<double, double> end_par, SfSolution *solution):
    BlockBoundaryCondition(type),
    parent_(solution),
    edgenmb_(edge_nmb),
    fbd_(fbd),
    domain_(end_par),
    approx_err_(-1.0)
  //===========================================================================
  {
  }


  //===========================================================================
  SfBoundaryCondition::SfBoundaryCondition(int edge_nmb, BdConditionType type, const Point& const_val,
					   pair<double, double> end_par, SfSolution *solution):
    BlockBoundaryCondition(type),
    parent_(solution),
    edgenmb_(edge_nmb),
    fbd_(0),
    const_val_(const_val),
    domain_(end_par),
    approx_err_(-1.0)
  //===========================================================================
  {
  }


  //===========================================================================
  SfBoundaryCondition::~SfBoundaryCondition()
  //===========================================================================
  {
  }


  //===========================================================================
  void SfBoundaryCondition::getCoefficientsEnumeration(std::vector<int>& local_enumeration)
  //===========================================================================
  {
    vector<int> full_enumeration;
    parent_->getBoundaryCoefficients(edgenmb_, full_enumeration);
    BsplineBasis bas = parent_->basis(1 - (edgenmb_ >> 1));    // Basis for v-dir if we are at umin, umax, or u-dir for vmin, vmax
    if (domain_.first <= domain_.second)
      {
	// Increasing order
	int pos = bas.knotIntervalFuzzy(domain_.first);
	pos -= bas.order() - 1;
	int end_pos = bas.knotIntervalFuzzy(domain_.second);
	local_enumeration.resize(end_pos + 1 - pos);
	for (int i = 0; pos <= end_pos; ++i, ++pos)
	  local_enumeration[i] = full_enumeration[pos];
      }
    else
      {
	// Decreasing order
	int end_pos = bas.knotIntervalFuzzy(domain_.second);
	end_pos -= bas.order() - 1;
	int pos = bas.knotIntervalFuzzy(domain_.first);
	local_enumeration.resize(pos + 1 - end_pos);
	for (int i = 0; pos >= end_pos; ++i, --pos)
	  local_enumeration[i] = full_enumeration[pos];
      }
  }


  //===========================================================================
  void SfBoundaryCondition::getCoefficientsEnumeration(std::vector<int>& local_enumeration_bd,
						       std::vector<int>& local_enumeration_bd2)
  //===========================================================================
  {
    vector<int> full_enumeration_bd, full_enumeration_bd2;
    parent_->getBoundaryCoefficients(edgenmb_,
				     full_enumeration_bd, full_enumeration_bd2);
    // Basis for v-dir if we are at umin, umax, or u-dir for vmin, vmax
    BsplineBasis bas = parent_->basis(1 - (edgenmb_ >> 1));
    if (domain_.first <= domain_.second)
      {
	// Increasing order
	int pos = bas.knotIntervalFuzzy(domain_.first);
	pos -= bas.order() - 1;
	int end_pos = bas.knotIntervalFuzzy(domain_.second);
	local_enumeration_bd.resize(end_pos + 1 - pos);
	local_enumeration_bd2.resize(end_pos + 1 - pos);
	for (int i = 0; pos <= end_pos; ++i, ++pos)
	  {
	    local_enumeration_bd[i] = full_enumeration_bd[pos];
	    local_enumeration_bd2[i] = full_enumeration_bd2[pos];
	  }
      }
    else
      {
	// Decreasing order
	int end_pos = bas.knotIntervalFuzzy(domain_.second);
	end_pos -= bas.order() - 1;
	int pos = bas.knotIntervalFuzzy(domain_.first);
	local_enumeration_bd.resize(pos + 1 - end_pos);
	local_enumeration_bd2.resize(pos + 1 - end_pos);
	for (int i = 0; pos >= end_pos; ++i, --pos)
	  {
	    local_enumeration_bd[i] = full_enumeration_bd[pos];
	    local_enumeration_bd2[i] = full_enumeration_bd2[pos];
	  }
      }
  }


  //===========================================================================
  void SfBoundaryCondition::getBdCoefficients(std::vector<std::pair<int, Point> >& coefs)
  //===========================================================================
  {
    if (!isDirichlet())
      return;

    const shared_ptr<SplineSurface> surf = parent_->getSolutionSurface();
    vector<int> coefs_enum;
    getCoefficientsEnumeration(coefs_enum);
    int coefs_size = (int)coefs_enum.size();
    int dim = surf->dimension();
    bool rational = surf->rational();
    int kdim = dim + (rational ? 1 : 0);
    coefs.resize(coefs_size);

    for (int i = 0; i < coefs_size; ++i)
      {
	Point p(dim);
	vector<double>::const_iterator it = surf->ctrl_begin() + coefs_enum[i] * kdim;
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
  SfBoundaryCondition::getBdCoefficients(vector<pair<int, Point> >& coefs_bd,
					 vector<pair<int, Point> >& coefs_bd2)
  //===========================================================================
  {
    if (!isDirichlet())
      return;

    const shared_ptr<SplineSurface> surf = parent_->getSolutionSurface();
    vector<int> coefs_enum_bd, coefs_enum_bd2;
    getCoefficientsEnumeration(coefs_enum_bd, coefs_enum_bd2);
    int coefs_size = (int)coefs_enum_bd.size();
    int dim = surf->dimension();
    bool rational = surf->rational();
    int kdim = dim + (rational ? 1 : 0);
    coefs_bd.resize(coefs_size);
    coefs_bd2.resize(coefs_size);

    for (int i = 0; i < coefs_size; ++i)
      {
	Point p_bd(dim), p_bd2(dim);
	vector<double>::const_iterator it_bd = surf->ctrl_begin() + coefs_enum_bd[i] * kdim;
	vector<double>::const_iterator it_bd2 = surf->ctrl_begin() + coefs_enum_bd2[i] * kdim;
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

#ifndef NDEBUG
#if 0
    // We write to file the sampled bd pts.
    vector<double> pts_bd(coefs_bd.size()*dim), pts_bd2(coefs_bd.size()*dim);
    assert(coefs_bd.size() == coefs_bd2.size());
    for (size_t ki = 0; ki < coefs_bd.size(); ++ki)
      {
	copy(coefs_bd[ki].second.begin(), coefs_bd[ki].second.end(), pts_bd.begin() + ki*dim);
	copy(coefs_bd2[ki].second.begin(), coefs_bd2[ki].second.end(), pts_bd2.begin() + ki*dim);
      }
    assert(dim == 3 || dim == 1);
    std::ofstream fileout("tmp/bd_pts.g2");
    if (dim == 1)
      {
	Go::PointCloud<1> bd_cloud(pts_bd.begin(), pts_bd.size()/dim);
	Go::PointCloud<1> bd_cloud2(pts_bd.begin(), pts_bd.size()/dim);
	bd_cloud.writeStandardHeader(fileout);
	bd_cloud.write(fileout);
	bd_cloud2.writeStandardHeader(fileout);
	bd_cloud2.write(fileout);
      }
    else if (dim == 3)
      {
	Go::PointCloud<3> bd_cloud(pts_bd.begin(), pts_bd.size()/dim);
	Go::PointCloud<3> bd_cloud2(pts_bd.begin(), pts_bd.size()/dim);
	bd_cloud.writeStandardHeader(fileout);
	bd_cloud.write(fileout);
	bd_cloud2.writeStandardHeader(fileout);
	bd_cloud2.write(fileout);
      }
    else
      {
      MESSAGE("Dim not supported!");
      }
#endif
#endif

  }


  //===========================================================================
  void SfBoundaryCondition::update()
  //===========================================================================
  {
    if (!isDirichlet())
      return;    // Only Dirichlet conditions might update the approximation curve

    double tol = getTolerances().gap;
    if (approx_err_ >= 0.0 && approx_err_ < tol)
      return;    // No update is done if approximation within tolerance has allready occured

    int ccw_edge_number;  // The input to SplineSurface::edgeCurve. Remeber that this value
                          // order is vmin-umax-vmax-umin while edgenmb_ uses umin-umax-vmin-vmax
    switch(edgenmb_)
      {
      case 0:   // umin
	ccw_edge_number = 3;
	break;
      case 1:   // umax
	ccw_edge_number = 1;
	break;
      case 2:   // vmin
	ccw_edge_number = 0;
	break;
      case 3:   // vmax
	ccw_edge_number = 2;
	break;
      }

    bool same_dir = domain_.first <= domain_.second;
    double start_par = same_dir ? domain_.first : domain_.second;
    double end_par = same_dir ? domain_.second : domain_.first;
    const shared_ptr<SplineSurface> sol_surf = parent_->getSolutionSurface();

    shared_ptr<SplineCurve> edge_curve = shared_ptr<SplineCurve>(sol_surf->edgeCurve(ccw_edge_number));
    shared_ptr<SplineCurve> bas_curve = shared_ptr<SplineCurve>(edge_curve->subCurve(start_par, end_par));
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
//    vector<int> local_enumeration;

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
  }


  //===========================================================================
  void SfBoundaryCondition::getBasisFunctions(int index_of_Gauss_point,
					      bool& u_dir,
					      vector<double>& basisValues,
					      vector<double>& basisDerivs) const
  //===========================================================================
  {
    u_dir = edgenmb_ == 2 || edgenmb_ == 3;
    double param = parent_->getGaussParameter(index_of_Gauss_point, u_dir ? 0 : 1);

    int ccw_edge_number;  // The input to SplineSurface::edgeCurve. Remeber that this value
                          // order is vmin-umax-vmax-umin while edgenmb_ uses umin-umax-vmin-vmax
    switch(edgenmb_)
      {
      case 0:   // umin
	ccw_edge_number = 3;
	break;
      case 1:   // umax
	ccw_edge_number = 1;
	break;
      case 2:   // vmin
	ccw_edge_number = 0;
	break;
      case 3:   // vmax
	ccw_edge_number = 2;
	break;
      }

    shared_ptr<SplineCurve> edge_curve = shared_ptr<SplineCurve>(parent_->getSolutionSurface()->edgeCurve(ccw_edge_number));
    edge_curve->computeBasis(param, basisValues, basisDerivs);
  }


  //===========================================================================
  shared_ptr<SplineCurve> SfBoundaryCondition::getSplineApproximation() const
  //===========================================================================
  {
    return bdcrv_cond_;
  }


  //===========================================================================
  void SfBoundaryCondition::updateBoundaryValue(BdCondFunctor* fbd)
  //===========================================================================
  {
    MESSAGE("updateBoundaryValue() not implemented");
  }


  //===========================================================================
  int SfBoundaryCondition::edgeNumber() const
  //===========================================================================
  {
    return edgenmb_;
  }


  //===========================================================================
  tpTolerances SfBoundaryCondition::getTolerances() const
  //===========================================================================
  {
    return parent_->getTolerances();
  }


}   // namespace Go
