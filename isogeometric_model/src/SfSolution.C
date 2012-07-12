//===========================================================================
//
// File : SfSolution.C
//
// Created: Wed Mar  3 09:52:55 2010
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================




#include "GoTools/isogeometric_model/SfSolution.h"
#include "GoTools/isogeometric_model/IsogeometricSfBlock.h"
#include "GoTools/geometry/GapRemoval.h"
#include "GoTools/geometry/SurfaceTools.h"
#include <algorithm>


using std::vector;
using std::pair;
using std::max;

namespace Go
{

  //===========================================================================
  SfSolution::SfSolution(IsogeometricSfBlock* parent, shared_ptr<SplineSurface> sol_surf):
    solution_(sol_surf),
    parent_(parent)
  //===========================================================================
  {
  }

  //===========================================================================
  SfSolution::~SfSolution()
  //===========================================================================
  {
  }

  //===========================================================================
  SfSolution* SfSolution::asSfSolution()
  //===========================================================================
  {
    return this;
  }


  //===========================================================================
  void SfSolution::addBoundaryCondition(int edge_nmb, BdConditionType type, BdCondFunctor *fbd,
					pair<double, double> end_par)
  //===========================================================================
  {
    boundary_conditions_.push_back(shared_ptr<SfBoundaryCondition>(new SfBoundaryCondition(edge_nmb, type, fbd, end_par, this)));
  }


  //===========================================================================
  void SfSolution::addBoundaryCondition(int edge_nmb, BdConditionType type, const Point& const_val,
					pair<double, double> end_par)
  //===========================================================================
  {
    boundary_conditions_.push_back(shared_ptr<SfBoundaryCondition>(new SfBoundaryCondition(edge_nmb, type, const_val, end_par, this)));
  }


  //===========================================================================
  void SfSolution::addDirichletPointBdCond(double param[],
					   Point& condition_value)
  //===========================================================================
  {
    MESSAGE("addDirichletPointBdCond() not implemented");
  }


  //===========================================================================
  int SfSolution::getNmbOfBoundaryConditions() const
  //===========================================================================
  {
      return (int)boundary_conditions_.size();
  }


  //===========================================================================
  shared_ptr<SfBoundaryCondition> SfSolution::getBoundaryCondition(int index) const
  //===========================================================================
  {
      if (index < 0 || index >= (int)boundary_conditions_.size())
      {
	shared_ptr<SfBoundaryCondition> bd_cond;
	return bd_cond;
      }

    return boundary_conditions_[index];
  }


  //===========================================================================
  void SfSolution::getEdgeBoundaryConditions(int edge_number, 
					     vector<shared_ptr<SfBoundaryCondition> >& bd_cond) const
  //===========================================================================
  {
    for (int i = 0; boundary_conditions_.size(); ++i)
      if (boundary_conditions_[i]->edgeNumber() == edge_number)
	bd_cond.push_back(boundary_conditions_[i]);
  }


  //===========================================================================
  void SfSolution::getEdgeBoundaryConditions(vector<shared_ptr<SfBoundaryCondition> >& bd_cond) const
  //===========================================================================
  {
    for (int i = 0; boundary_conditions_.size(); ++i)
      bd_cond.push_back(boundary_conditions_[i]);
  }


  //===========================================================================
  int SfSolution::getNmbOfPointBdConditions() const
  //===========================================================================
  {
      return (int)point_bd_cond_.size();
  }


  //===========================================================================
  shared_ptr<SfPointBdCond> SfSolution::getPointBdCondition(int index) const
  //===========================================================================
  {
      if (index < 0 || index >= (int)point_bd_cond_.size())
      {
	shared_ptr<SfPointBdCond> bd_cond;
	return bd_cond;
      }

    return point_bd_cond_[index];
  }


  //===========================================================================
  void SfSolution::getEdgePointBdConditions(int edge_number, 
					    vector<shared_ptr<SfPointBdCond> >& bd_cond) const
  //===========================================================================
  {
    for (int i = 0; point_bd_cond_.size(); ++i)
      if (point_bd_cond_[i]->edgeNumber() == edge_number)
	bd_cond.push_back(point_bd_cond_[i]);
  }


  //===========================================================================
  void SfSolution::getPointBdCond(vector<shared_ptr<SfPointBdCond> >& bd_cond) const
  //===========================================================================
  {
    for (int i = 0; point_bd_cond_.size(); ++i)
      bd_cond.push_back(point_bd_cond_[i]);
  }


  //===========================================================================
  bool SfSolution::matchingSplineSpace(BlockSolution* other) const
  //===========================================================================
  {
    vector<int> edges, edges_other;
    vector<bool> equal_orient;
    vector<bool> space_matches;
    neighbourInfo(other, edges, edges_other, equal_orient, space_matches);

    for (int i = 0; i < (int)edges.size(); ++i)
      if (!space_matches[i])
	return false;

    return true;
  }

  //===========================================================================
  void SfSolution::getMatchingCoefficients(BlockSolution* other, vector<pair<int,int> >& enumeration, int match_pos) const
  //===========================================================================
  {
    vector<int> edges, edges_other;
    vector<bool> equal_orient;
    vector<bool> space_matches;
    neighbourInfo(other, edges, edges_other, equal_orient, space_matches);

    if ((int)edges.size() <= match_pos || match_pos < 0)
      return;
    if (!space_matches[match_pos])
      return;

    SfSolution* sf_other = other->asSfSolution();
    vector<int> coefs_this, coefs_other;

    SurfaceTools::getCoefEnumeration(solution_, edges[match_pos], coefs_this);
    SurfaceTools::getCoefEnumeration(sf_other->solution_, edges_other[match_pos], coefs_other);

    int coefs_size = (int)coefs_this.size();
    if (equal_orient[match_pos])
      for (int i = 0; i < coefs_size; ++i)
	enumeration.push_back(pair<int, int>(coefs_this[i], coefs_other[i]));
    else
      for (int i = 0; i < coefs_size; ++i)
	enumeration.push_back(pair<int, int>(coefs_this[i], coefs_other[coefs_size-1-i]));
  }

  //===========================================================================
  void SfSolution::getBoundaryCoefficients(int boundary, vector<int>& enumeration) const
  //===========================================================================
  {
    SurfaceTools::getCoefEnumeration(solution_, boundary, enumeration);
  }


  //===========================================================================
  void
  SfSolution::getBoundaryCoefficients(int boundary,
				      std::vector<int>& enumeration_bd,
				      std::vector<int>& enumeration_bd2) const
  //===========================================================================
  {
    SurfaceTools::getCoefEnumeration(solution_, boundary,
				     enumeration_bd, enumeration_bd2);
  }


  //===========================================================================
  void SfSolution::makeMatchingSplineSpace(BlockSolution* other)
  //===========================================================================
  {
    double tol = getTolerances().gap;
    bool changed = false;

    vector<int> edges, edges_other;
    vector<bool> equal_orient;
    vector<bool> space_matches;
    neighbourInfo(other, edges, edges_other, equal_orient, space_matches);

    shared_ptr<SplineSurface> solution_other = other->asSfSolution()->solution_;

    for (int i = 0; i < (int)edges.size(); ++i)
      if (!space_matches[i])
	{
	  // Spline spaces are not equal

	  // Get boundary curve on first solution surface
	  shared_ptr<SplineCurve> bd_crv_this;
	  switch (edges[i])
	    {
	    case 0:
	      bd_crv_this = shared_ptr<SplineCurve>(solution_->constParamCurve(solution_->startparam_v(), false));
	    break;
	    case 1:
	      bd_crv_this = shared_ptr<SplineCurve>(solution_->constParamCurve(solution_->endparam_v(), false));
	      break;
	    case 2:
	      bd_crv_this = shared_ptr<SplineCurve>(solution_->constParamCurve(solution_->startparam_u(), true));
	      break;
	    case 3:
	      bd_crv_this = shared_ptr<SplineCurve>(solution_->constParamCurve(solution_->endparam_u(), true));
	      break;
	    }
	  shared_ptr<CurveOnSurface> curve_this(new CurveOnSurface(solution_, bd_crv_this, false));

	  // Get boundary curve on second surface
	  shared_ptr<SplineCurve> bd_crv_neighbour;
	  switch (edges_other[i])
	    {
	    case 0:
	      bd_crv_neighbour = shared_ptr<SplineCurve>(solution_other->constParamCurve(solution_other->startparam_v(), false));
	      break;
	    case 1:
	      bd_crv_neighbour = shared_ptr<SplineCurve>(solution_other->constParamCurve(solution_other->endparam_v(), false));
	      break;
	    case 2:
	      bd_crv_neighbour = shared_ptr<SplineCurve>(solution_other->constParamCurve(solution_other->startparam_u(), true));
	      break;
	    case 3:
	      bd_crv_neighbour = shared_ptr<SplineCurve>(solution_other->constParamCurve(solution_other->endparam_u(), true));
	      break;
	    }
	  shared_ptr<CurveOnSurface> curve_neighbour(new CurveOnSurface(solution_other, bd_crv_neighbour, false));

	  // Get limiting parameters
	  double start1 = (i < 2) ? solution_->startparam_v() : solution_->startparam_u();
	  double end1 = (i < 2) ? solution_->endparam_v() : solution_->endparam_u();
	  double start2 = (edges_other[i] < 2) ? solution_other->startparam_v() : solution_other->startparam_u();
	  double end2 = (edges_other[i] < 2) ? solution_other->endparam_v() : solution_other->endparam_u();

	  // Get corner points
	  Point p_start, p_end;
	  bd_crv_this->point(p_start, start1);
	  bd_crv_this->point(p_end, end1);

	  // Make uniform
	  bool eq_or = equal_orient[i];
	  GapRemoval::removeGapSpline(solution_, curve_this, start1, end1,
				      solution_other, curve_neighbour, start2, end2,
				      p_start, p_end, tol, &eq_or);
	  changed = true;
	}

    if (changed)
      updateConditions();
  }

  //===========================================================================
  void SfSolution::increaseDegree(int new_degree, int pardir)
  //===========================================================================
  {
    bool changed = false;

    int new_order = new_degree + 1;
    if (pardir == 0 && new_order > solution_->order_u())
      {
	solution_->raiseOrder(new_order - solution_->order_u(), 0);
	changed = true;
      }
    else if (pardir == 1 && new_order > solution_->order_v())
      {
	solution_->raiseOrder(0, new_order - solution_->order_v());
	changed = true;
      }

    if (changed)
      updateConditions();
  }

  //===========================================================================
  void SfSolution::insertKnots(const vector<int>& knot_intervals, int pardir)
  //===========================================================================
  {
    vector<double> old_knots;
    solution_->basis(pardir).knotsSimple(old_knots);
    vector<double> new_knots;

    for (vector<int>::const_iterator it = knot_intervals.begin(); it != knot_intervals.end(); ++it)
	if ((*it) < (int)old_knots.size() - 1)
	new_knots.push_back(0.5 * (old_knots[*it] + old_knots[(*it) + 1]));

    if (new_knots.size() > 0)
      {
	insertKnots(new_knots, pardir);
	updateConditions();
      }
  }

  //===========================================================================
  void SfSolution::insertKnots(const vector<double>& knots, int pardir)
  //===========================================================================
  {
    if (pardir == 0)
      solution_->insertKnot_u(knots);
    else
      solution_->insertKnot_v(knots);
  }


  //===========================================================================
  void SfSolution::erasePreEvaluatedBasisFunctions()
  //===========================================================================
  {
    shared_ptr<preEvaluationSf> empty;
    evaluated_grid_ = empty;
  }

  //===========================================================================
  void SfSolution::performPreEvaluation(vector<vector<double> >& Gauss_par)
  //===========================================================================
  {
    ASSERT (Gauss_par.size() == 2);

    erasePreEvaluatedBasisFunctions();
    evaluated_grid_ = shared_ptr<preEvaluationSf>(new preEvaluationSf);

    vector<double> par_u = Gauss_par[0];
    vector<double> par_v = Gauss_par[1];
    int nmb_par_u = (int)par_u.size();
    int nmb_par_v = (int)par_v.size();

    // int num_u = solution_->numCoefs_u();
    // int num_v = solution_->numCoefs_v();
    int ord_u = solution_->order_u();
    int ord_v = solution_->order_v();

    evaluated_grid_->gauss_par1_.resize(nmb_par_u);
    copy(par_u.begin(), par_u.end(), evaluated_grid_->gauss_par1_.begin());
    evaluated_grid_->gauss_par2_.resize(nmb_par_v);
    copy(par_v.begin(), par_v.end(), evaluated_grid_->gauss_par2_.begin());

    evaluated_grid_->basisvals_u_.resize(nmb_par_u * ord_u * 2);
    evaluated_grid_->basisvals_v_.resize(nmb_par_v * ord_v * 2);
    evaluated_grid_->left_u_.resize(nmb_par_u);
    evaluated_grid_->left_v_.resize(nmb_par_v);

    solution_->basis_u().computeBasisValues(&par_u[0], &par_u[0]+nmb_par_u,
					    &(evaluated_grid_->basisvals_u_[0]),
					    &(evaluated_grid_->left_u_[0]), 1);
    solution_->basis_v().computeBasisValues(&par_v[0], &par_v[0]+nmb_par_v,
					    &(evaluated_grid_->basisvals_v_[0]),
					    &(evaluated_grid_->left_v_[0]), 1);
    getGeometrySurface()->gridEvaluator(par_u, par_v,
					evaluated_grid_->points_,
					evaluated_grid_->deriv_u_,
					evaluated_grid_->deriv_v_);
  }


  //===========================================================================
  void SfSolution::getBasisFunctions(int index_of_Gauss_point1,
				     int index_of_Gauss_point2,
				     vector<double>& basisValues,
				     vector<double>& basisDerivs_u,
				     vector<double>& basisDerivs_v) const
  //===========================================================================
  {
    if (evaluated_grid_.get() == NULL)
      return;
    if (index_of_Gauss_point1 < 0 ||
	index_of_Gauss_point1 >= (int)evaluated_grid_->gauss_par1_.size() ||
	index_of_Gauss_point2 < 0 ||
	index_of_Gauss_point2 >= (int)evaluated_grid_->gauss_par2_.size())
      return;

    solution_->computeBasis(evaluated_grid_->basisvals_u_.begin() 
			    + 2 * index_of_Gauss_point1 * solution_->order_u(),
			    evaluated_grid_->basisvals_v_.begin() 
			    + 2 * index_of_Gauss_point2 * solution_->order_v(),
			    evaluated_grid_->left_u_[index_of_Gauss_point1],
			    evaluated_grid_->left_v_[index_of_Gauss_point2],
			    basisValues,
			    basisDerivs_u,
			    basisDerivs_v);
  }


  //===========================================================================
  void SfSolution::getBasisFunctions(double param1,
				     double param2,
				     vector<double>& basisValues,
				     vector<double>& basisDerivs_u,
				     vector<double>& basisDerivs_v) const
  //===========================================================================
  {
    double param[2];
    param[0] = param1;
    param[1] = param2;
    solution_->computeBasis(param, basisValues, basisDerivs_u, basisDerivs_v);
  }


  //===========================================================================
  void SfSolution::getBasisFunctionValues(int basis_func_id_u, int basis_func_id_v,
					  vector<int>& index_of_Gauss_points1,
					  vector<int>& index_of_Gauss_points2,
					  vector<double>& basisValues,
					  vector<double>& basisDerivs_u,
					  vector<double>& basisDerivs_v) const
  //===========================================================================
  {
    const int order_u = solution_->order_u();
    const int order_v = solution_->order_v();
    const int deg_u = order_u - 1;
    const int deg_v = order_v - 1;

    const int dim = solution_->dimension();

    // We run through the evaluated_grid_ and compute basis values for
    // the Gauss points in the support of our basis function.
    for (size_t kj = 0; kj < evaluated_grid_->left_v_.size(); ++kj)
	if (evaluated_grid_->left_v_[kj] - deg_v <= basis_func_id_v &&
	    basis_func_id_v < evaluated_grid_->left_v_[kj] + 1)
	{
	    int local_ind_v = basis_func_id_v + deg_v - evaluated_grid_->left_v_[kj];
	    for (size_t ki = 0; ki < evaluated_grid_->left_u_.size(); ++ki)
		if (evaluated_grid_->left_u_[ki] - deg_u <= basis_func_id_u &&
		    basis_func_id_u < evaluated_grid_->left_u_[ki] + 1)
		{
		    // We have found a Gauss point in the support of the function.
		    int local_ind_u = basis_func_id_u + deg_u - evaluated_grid_->left_u_[ki];

		    // We add the contribution from the sf coef (and
		    // weight for rational case).
		    vector<double> local_basisValues;
		    vector<double> local_basisDerivs_u;
		    vector<double> local_basisDerivs_v;
		    solution_->computeBasis(evaluated_grid_->basisvals_u_.begin()
					    + 2 * ki * order_u,
					    evaluated_grid_->basisvals_v_.begin() 
					    + 2 * kj * order_v,
					    evaluated_grid_->left_u_[ki],
					    evaluated_grid_->left_v_[kj],
					    local_basisValues,
					    local_basisDerivs_u,
					    local_basisDerivs_v);
		    basisValues.insert(basisValues.end(),
				       local_basisValues.begin() + (local_ind_v*order_u + local_ind_u)*dim,
				       local_basisValues.begin() + (local_ind_v*order_u + local_ind_u + 1)*dim);
		    basisDerivs_u.insert(basisDerivs_u.end(),
					 local_basisDerivs_u.begin() + (local_ind_v*order_u + local_ind_u)*dim,
					 local_basisDerivs_u.begin() + (local_ind_v*order_u + local_ind_u + 1)*dim);
		    basisDerivs_v.insert(basisDerivs_v.end(),
					 local_basisDerivs_v.begin() + (local_ind_v*order_u + local_ind_u)*dim,
					 local_basisDerivs_v.begin() + (local_ind_v*order_u + local_ind_u + 1)*dim);

		    // Storing the index of the Gauss points.
		    index_of_Gauss_points1.push_back((int)ki);
		    index_of_Gauss_points2.push_back((int)kj);

		}
	}
  }


  //===========================================================================
  double SfSolution::getJacobian(vector<int>& index_of_Gauss_point) const
  //===========================================================================
  {
    ASSERT (index_of_Gauss_point.size() == 2);
    if (evaluated_grid_.get() == NULL)
      return 0.0;

    int dim = getGeometrySurface()->dimension();
    ASSERT (dim == 2);
    int pos = dim * (index_of_Gauss_point[1] * ((int)evaluated_grid_->gauss_par1_.size())
		     + index_of_Gauss_point[0]);
    return (evaluated_grid_->deriv_u_[pos] * evaluated_grid_->deriv_v_[pos+1] -
	    evaluated_grid_->deriv_u_[pos+1] * evaluated_grid_->deriv_v_[pos]);
  }

  //===========================================================================
  void SfSolution::valuesInGaussPoint(const vector<int>& index_of_Gauss_point, vector<Point>& derivs) const
  //===========================================================================
  {
    if (evaluated_grid_.get() == NULL)
      return;

    int dim = getGeometrySurface()->dimension();
    int pos = dim * (index_of_Gauss_point[0] + index_of_Gauss_point[1] * (int)evaluated_grid_->gauss_par1_.size());
    derivs.resize(3);
    derivs[0] = Point(evaluated_grid_->points_.begin() + pos,
		      evaluated_grid_->points_.begin() + pos + dim);
    derivs[1] = Point(evaluated_grid_->deriv_u_.begin() + pos,
		      evaluated_grid_->deriv_u_.begin() + pos + dim);
    derivs[2] = Point(evaluated_grid_->deriv_v_.begin() + pos,
		      evaluated_grid_->deriv_v_.begin() + pos + dim);
  }

  //===========================================================================
  void SfSolution::setSolutionCoefficients(const vector<double>& coefs)
  //===========================================================================
  {
    int ncoefs = solution_->numCoefs_u() * solution_->numCoefs_v();
    int dim = solution_->dimension();
    copy(coefs.begin(), coefs.begin() + ncoefs * dim, solution_->coefs_begin());
    if (solution_->rational())
      {
	vector<double>::const_iterator c_it = solution_->coefs_begin();
	vector<double>::iterator r_it = solution_->rcoefs_begin();
	for (int i = 0; i < ncoefs; ++i, ++r_it)
	  {
	    double w = r_it[dim];
	    for (int j = 0; j < dim; ++j, ++r_it, ++c_it)
	      (*r_it) = w * (*c_it);
	  }
      }
  }


  //===========================================================================
  shared_ptr<SplineSurface> SfSolution::getSolutionSurface() const
  //===========================================================================
  {
    return solution_;
  }


  //===========================================================================
  shared_ptr<SplineSurface> SfSolution::getGeometrySurface() const
  //===========================================================================
  {
    return parent_->surface();
  }


  //===========================================================================
  void SfSolution::setMinimumDegree(int degree)
  //===========================================================================
  {
    int order = degree + 1;

    int raise_u = max(parent_->surface()->order_u(), order) - solution_->order_u();
    int raise_v = max(parent_->surface()->order_v(), order) - solution_->order_v();

    raise_u = max(raise_u, 0);
    raise_v = max(raise_v, 0);

    if (raise_u > 0 || raise_v > 0)
      {
	solution_->raiseOrder(raise_u, raise_v);
	updateConditions();
      }
  }


  //===========================================================================
  void SfSolution::refineToGeometry(int pardir)
  //===========================================================================
  {
    bool changed = false;

    BsplineBasis base_solution = solution_->basis(pardir);
    BsplineBasis base_geometry = parent_->surface()->basis(pardir);

    int order_solution = base_solution.order();
    int order_geometry = base_geometry.order();
    if (order_solution < order_geometry)
      {
	if (pardir == 0)
	  solution_->raiseOrder(order_geometry - order_solution, 0);
	else
	  solution_->raiseOrder(0, order_geometry - order_solution);
	order_solution = order_geometry;
	changed = true;
      }

    vector<double> geo_knots;
    base_geometry.knotsSimple(geo_knots);
    int order_diff = order_solution - order_geometry;
    vector<double> new_knots;
    for (int i = 0; i < (int)geo_knots.size(); ++i)
      {
	double knot_val = geo_knots[i];
	int knots_needed = order_diff + base_geometry.knotMultiplicity(knot_val) - base_solution.knotMultiplicity(knot_val);
	if (knots_needed > 0)
	  for (int j = 0; j < knots_needed; ++j)
	    new_knots.push_back(knot_val);
      }

    if (new_knots.size() > 0)
      {
	if (pardir == 0)
	  solution_->insertKnot_u(new_knots);
	else
	  solution_->insertKnot_v(new_knots);
	changed = true;
      }

    if(changed)
      updateConditions();
  }


  //===========================================================================
  tpTolerances SfSolution::getTolerances() const
  //===========================================================================
  {
    return parent_->getTolerances();
  }


  //===========================================================================
  int SfSolution::nmbCoefs() const
  //===========================================================================
  {
    return solution_->numCoefs_u() * solution_->numCoefs_v();
  }


  //===========================================================================
  int SfSolution::nmbCoefs(int pardir) const
  //===========================================================================
  {
    if (pardir == 0)
      return solution_->numCoefs_u();
    else
      return solution_->numCoefs_v();
  }


  //===========================================================================
  int SfSolution::degree(int pardir) const
  //===========================================================================
  {
    if (pardir == 0)
      return solution_->order_u() - 1;
    else
      return solution_->order_v() - 1;
  }


  //===========================================================================
  vector<double> SfSolution::knots(int pardir) const
  //===========================================================================
  {
    const BsplineBasis bas = (pardir == 0) ? solution_->basis_u() : solution_->basis_v();
    vector<double> result(bas.order() + bas.numCoefs());
    vector<double>::const_iterator bas_begin = bas.begin();
    vector<double>::const_iterator bas_end = bas.end();
    copy (bas_begin, bas_end, result.begin());
    return result;
  }


  //===========================================================================
  vector<double> SfSolution::distinctKnots(int pardir) const
  //===========================================================================
  {
    const BsplineBasis bas = (pardir == 0) ? solution_->basis_u() : solution_->basis_v();
    vector<double> result;
    bas.knotsSimple(result);
    return result;
  }



  //===========================================================================
  BsplineBasis SfSolution::basis(int pardir) const
  //===========================================================================
  {
    return solution_->basis(pardir);
  }


  //===========================================================================
  int SfSolution::dimension() const
  //===========================================================================
  {
    return solution_->dimension();
  }


  //===========================================================================
  void SfSolution::updateConditions()
  //===========================================================================
  {
      for (int i = 0; i < (int)boundary_conditions_.size(); ++i)
      boundary_conditions_[i]->update();
  }


  //===========================================================================
  double SfSolution::getGaussParameter(int index_of_Gauss_point, int pardir) const
  //===========================================================================
  {
    if (evaluated_grid_.get() == NULL)
      return 0.0;

    if (pardir == 0)
      {
	  if (index_of_Gauss_point < 0 || index_of_Gauss_point >= (int)evaluated_grid_->gauss_par1_.size())
	  return 0.0;
	else
	  return evaluated_grid_->gauss_par1_[index_of_Gauss_point];
      }
    else
      {
	  if (index_of_Gauss_point < 0 || index_of_Gauss_point >= (int)evaluated_grid_->gauss_par2_.size())
	  return 0.0;
	else
	  return evaluated_grid_->gauss_par2_[index_of_Gauss_point];
      }
  }


  //===========================================================================
  void SfSolution::neighbourInfo(BlockSolution* other, vector<int>& edges, vector<int>& edges_other,
				 vector<bool>& equal_oriented, vector<bool>& space_matches) const
  //===========================================================================
  {
    double tol = getTolerances().gap;

    SfSolution* sf_other = other->asSfSolution();
    parent_->getNeighbourInfo(sf_other->parent_, edges, edges_other, equal_oriented);

    space_matches.resize(edges.size());
    for (int i = 0; i < (int)edges.size(); ++i)
      {
	BsplineBasis basis1 = (edges[i] < 2) ? (solution_->basis_v()) : (solution_->basis_u());
	BsplineBasis basis2_original = (edges_other[i] < 2) ? (sf_other->solution_->basis_v()) : (sf_other->solution_->basis_u());
	BsplineBasis basis2(basis2_original);   // Make copy so we do not manipulate original basis
	if (!equal_oriented[i])
	  basis2.reverseParameterDirection();
	basis2.rescale(basis1.startparam(), basis1.endparam());

	space_matches[i] = basis1.sameSplineSpace(basis2, tol);
      }
  }

}   // namespace Go
