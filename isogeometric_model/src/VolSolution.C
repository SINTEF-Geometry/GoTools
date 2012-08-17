//===========================================================================
//
// File : VolSolution.C
//
// Created: Mon Mar  8 14:10:32 2010
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================



#include "GoTools/isogeometric_model/VolSolution.h"
#include "GoTools/isogeometric_model/IsogeometricVolBlock.h"
#include "GoTools/trivariate/VolumeTools.h"


using std::max;
using std::pair;


namespace Go
{

  //===========================================================================
  VolSolution::VolSolution(IsogeometricVolBlock* parent, shared_ptr<SplineVolume> sol_vol):
    parent_(parent),
    solution_(sol_vol)
  //===========================================================================
  {
  }

  //===========================================================================
  VolSolution::~VolSolution()
  //===========================================================================
  {
  }

  //===========================================================================
  VolSolution* VolSolution::asVolSolution()
  //===========================================================================
  {
    return this;
  }

  //===========================================================================
  void VolSolution::addBoundaryCondition(int face_nmb, BdConditionType type, BdCondFunctor *fbd,
					 vector<pair<double, double> >& domain)
					 // vector<pair<Point, Point> >& polygon)
  //===========================================================================
  {
    boundary_conditions_.push_back(shared_ptr<VolBoundaryCondition>
				   (new VolBoundaryCondition(face_nmb, type, fbd, domain, this)));
  }

  //===========================================================================
  void VolSolution::addDirichletPointBdCond(double param[],
			       Point& condition_value)
  //===========================================================================
  {
    MESSAGE("addDirichletPointBdCond() not implemented");
  }

  //===========================================================================
  int VolSolution::getNmbOfBoundaryConditions() const
  //===========================================================================
  {
    return (int)boundary_conditions_.size();
  }

  //===========================================================================
  shared_ptr<VolBoundaryCondition> VolSolution::getBoundaryCondition(int index) const
  //===========================================================================
  {
      if (index < 0 || index >= (int)boundary_conditions_.size())
      {
	shared_ptr<VolBoundaryCondition> bd_cond;
	return bd_cond;
      }

    return boundary_conditions_[index];
  }

  //===========================================================================
  void VolSolution::getFaceBoundaryConditions(int face_number, 
				 vector<shared_ptr<VolBoundaryCondition> >& bd_cond) const
  //===========================================================================
  {
    for (int i = 0; boundary_conditions_.size(); ++i)
      if (boundary_conditions_[i]->faceNumber() == face_number)
	bd_cond.push_back(boundary_conditions_[i]);
  }

  //===========================================================================
  void VolSolution::getFaceBoundaryConditions(vector<shared_ptr<VolBoundaryCondition> >& bd_cond) const
  //===========================================================================
  {
    for (int i = 0; boundary_conditions_.size(); ++i)
      bd_cond.push_back(boundary_conditions_[i]);
  }

  //===========================================================================
  int VolSolution::getNmbOfPointBdConditions() const
  //===========================================================================
  {
    return (int)point_bd_cond_.size();
  }

  //===========================================================================
  shared_ptr<VolPointBdCond> VolSolution::getPointBdCondition(int index) const
  //===========================================================================
  {
    if (index < 0 || index >= (int)point_bd_cond_.size())
      {
	shared_ptr<VolPointBdCond> bd_cond;
	return bd_cond;
      }

    return point_bd_cond_[index];
  }

  //===========================================================================
  void VolSolution::getFacePointBdConditions(int face_number, 
				vector<shared_ptr<VolPointBdCond> >& bd_cond) const
  //===========================================================================
  {
    for (int i = 0; point_bd_cond_.size(); ++i)
      if (point_bd_cond_[i]->faceNumber() == face_number)
	bd_cond.push_back(point_bd_cond_[i]);
  }

  //===========================================================================
  void VolSolution::getPointBdCond(vector<shared_ptr<VolPointBdCond> >& bd_cond) const
  //===========================================================================
  {
    for (int i = 0; point_bd_cond_.size(); ++i)
      bd_cond.push_back(point_bd_cond_[i]);
  }

  //===========================================================================
  bool VolSolution::matchingSplineSpace(BlockSolution* other) const
  //===========================================================================
  {
    vector<int> faces, faces_other;
    vector<int> orientation;
    vector<bool> same_dir_order;
    vector<bool> space_matches;
    neighbourInfo(other, faces, faces_other, orientation, same_dir_order, space_matches);

    for (int i = 0; i < (int)faces.size(); ++i)
      if (!space_matches[i])
	return false;

    return true;
  }

  //===========================================================================
  void VolSolution::getMatchingCoefficients(BlockSolution* other,
					    vector<pair<int,int> >& enumeration,
					    int match_pos) const
  //===========================================================================
  {
    MESSAGE("getMatchingCoefficients() under construction");

    vector<int> faces, faces_other;
    vector<int> orientation;
    vector<bool> same_dir_order;
    vector<bool> space_matches;
    neighbourInfo(other, faces, faces_other, orientation, same_dir_order, space_matches);

    if ((int)faces.size() <= match_pos || match_pos < 0)
      return;
    if (!space_matches[match_pos])
      return;

    VolSolution* vol_other = other->asVolSolution();
    vector<int> coefs_this, coefs_other;

    // @@sbr But can we really use the same match_pos for both vol blocks?
    VolumeTools::getVolCoefEnumeration(solution_, faces[match_pos], coefs_this);
    VolumeTools::getVolCoefEnumeration(vol_other->solution_, faces_other[match_pos], coefs_other);
    // Expecting that the coef ordering is u-dir first, then v-dir, finally w-dir.

    // We check if the sfs (in the volume intersection) have corr u- and v-dir.
    bool dir_order_ok = parent_->sameDirOrder(match_pos);
    int neighb_const_dir = faces[match_pos]/2;
    bool u_rev = false;//orientation[match_pos];
    bool v_rev = false;

//    int nb_face = faces_other[match_pos];

#if 0
    // We make sure the coefs for faces_other match those of faces.
    int coefs_size = (int)coefs_this.size();
    if (equal_orient[match_pos])
      for (int i = 0; i < coefs_size; ++i)
	enumeration.push_back(pair<int, int>(coefs_this[i], coefs_other[i]));
    else
      for (int i = 0; i < coefs_size; ++i)
	enumeration.push_back(pair<int, int>(coefs_this[i], coefs_other[coefs_size-1-i]));
#endif
  }

  //===========================================================================
  void VolSolution::getBoundaryCoefficients(int boundary,
					    vector<int>& enumeration) const
  //===========================================================================
  {
    MESSAGE("getBoundaryCoefficients() not implemented");
  }


  //===========================================================================
  void
  VolSolution::getBoundaryCoefficients(int boundary,
				       std::vector<int>& enumeration_bd,
				       std::vector<int>& enumeration_bd2) const
  //===========================================================================
  {
    MESSAGE("getBoundaryCoefficients() not implemented");
  }


  //===========================================================================
  void VolSolution::makeMatchingSplineSpace(BlockSolution* other)
  //===========================================================================
  {
    MESSAGE("makeMatchingSplineSpace() not implemented");
  }

  //===========================================================================
  void VolSolution::increaseDegree(int new_degree, int pardir)
  //===========================================================================
  {
    bool changed = false;

    int curr_order = solution_->order(pardir);
    int new_order = new_degree + 1;
    int raise_order = new_order - curr_order;
    if (pardir == 0 && raise_order > 0)
      {
	solution_->raiseOrder(raise_order, 0, 0);
	changed = true;
      }
    else if (pardir == 1 && raise_order > 0)
      {
	solution_->raiseOrder(0, raise_order, 0);
	changed = true;
      }
    else if (pardir == 2 && raise_order > 0)
      {
	solution_->raiseOrder(0, 0, raise_order);
	changed = true;
      }

    if (changed)
      updateConditions();
  }

  //===========================================================================
  void VolSolution::insertKnots(const vector<int>& knot_intervals, int pardir)
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
  void VolSolution::insertKnots(const vector<double>& knots, int pardir)
  //===========================================================================
  {
    solution_->insertKnot(pardir, knots);
  }

  //===========================================================================
  void VolSolution::erasePreEvaluatedBasisFunctions()
  //===========================================================================
  {
    shared_ptr<preEvaluationVol> empty;
    evaluated_grid_ = empty;
  }

  //===========================================================================
  void VolSolution::performPreEvaluation(vector<vector<double> >& Gauss_par)
  //===========================================================================
  {
    MESSAGE("performPreEvaluation() not implemented");
  }

  //===========================================================================
  void VolSolution::getBasisFunctions(int index_of_Gauss_point1,
				      int index_of_Gauss_point2,
				      int index_of_Gauss_point3,
				      shared_ptr<BasisDerivs> result) const
  //===========================================================================
  {
    MESSAGE("getBasisFunctions() under construction");

    if (evaluated_grid_.get() == NULL)
      return;
    if (index_of_Gauss_point1 < 0 ||
	index_of_Gauss_point1 >= (int)evaluated_grid_->gauss_par1_.size() ||
	index_of_Gauss_point2 < 0 ||
	index_of_Gauss_point2 >= (int)evaluated_grid_->gauss_par2_.size() ||
	index_of_Gauss_point3 < 0 ||
	index_of_Gauss_point3 >= (int)evaluated_grid_->gauss_par3_.size())
      return;

    solution_->computeBasis(evaluated_grid_->basisvals_u_.begin() 
			    [2 * index_of_Gauss_point1 * solution_->order(0)],
			    evaluated_grid_->basisvals_v_.begin()
			    [2 * index_of_Gauss_point2 * solution_->order(1)],
			    evaluated_grid_->basisvals_w_.begin() 
			    [2 * index_of_Gauss_point2 * solution_->order(2)],
			    *result);
			    // evaluated_grid_->left_u_[index_of_Gauss_point1],
			    // evaluated_grid_->left_v_[index_of_Gauss_point2],
			    // evaluated_grid_->left_w_[index_of_Gauss_point3],
			    // basisValues,
			    // basisDerivs_u,
			    // basisDerivs_v);

  }

  //===========================================================================
  void VolSolution::getBasisFunctions(double param1,
			 double param2,
			 double param3,
			 shared_ptr<BasisDerivs> result) const
  //===========================================================================
  {
    MESSAGE("getBasisFunctions() not implemented");
  }

  //===========================================================================
  void VolSolution::getBasisFunctionValues(int basis_func_id_u,
					   int basis_func_id_v,
					   int basis_func_id_w,
					   std::vector<int>& index_of_Gauss_points1,
					   std::vector<int>& index_of_Gauss_points2,
					   std::vector<int>& index_of_Gauss_points3,
					   shared_ptr<BasisDerivs> result) const
  //===========================================================================
  {
    MESSAGE("getBasisFunctions() not implemented");
  }

  //===========================================================================
  double VolSolution::getJacobian(vector<int>& index_of_Gauss_point) const
  //===========================================================================
  {
    MESSAGE("getBasisFunctions() not implemented");
    return 0.0;
#if 0
    ASSERT (index_of_Gauss_point.size() == 2);
    if (evaluated_grid_.get() == NULL)
      return 0.0;

    int dim = getGeometryVolume()->dimension();
    ASSERT (dim == 2);
    int pos = dim * (index_of_Gauss_point[1] * ((int)evaluated_grid_->gauss_par1_.size())
		     + index_of_Gauss_point[0]);

    // We compute the determinant of the Jacobian matrix.
    double det = 
    return (evaluated_grid_->deriv_u_[pos] * evaluated_grid_->deriv_v_[pos+1] -
	    evaluated_grid_->deriv_u_[pos+1] * evaluated_grid_->deriv_v_[pos]);
#endif
  }

  //===========================================================================
  void VolSolution::valuesInGaussPoint(const vector<int>& index_of_Gauss_point,
			  vector<Point>& derivs) const
  //===========================================================================
  {
    if (evaluated_grid_.get() == NULL)
      return;

    int dim = getGeometryVolume()->dimension();
    int pos = dim * (index_of_Gauss_point[0] +
		     index_of_Gauss_point[1] * (int)evaluated_grid_->gauss_par1_.size() +
		     index_of_Gauss_point[2] * (int)evaluated_grid_->gauss_par1_.size() + (int)evaluated_grid_->gauss_par1_.size() );
    derivs.resize(3);
    derivs[0] = Point(evaluated_grid_->points_.begin() + pos,
		      evaluated_grid_->points_.begin() + pos + dim);
    derivs[1] = Point(evaluated_grid_->deriv_u_.begin() + pos,
		      evaluated_grid_->deriv_u_.begin() + pos + dim);
    derivs[2] = Point(evaluated_grid_->deriv_v_.begin() + pos,
		      evaluated_grid_->deriv_v_.begin() + pos + dim);
  }

  //===========================================================================
  void VolSolution::setSolutionCoefficients(const vector<double>& coefs)
  //===========================================================================
  {
    int ncoefs = solution_->numCoefs(0) * solution_->numCoefs(1) * solution_->numCoefs(2);
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
  shared_ptr<SplineVolume> VolSolution::getSolutionVolume() const
  //===========================================================================
  {
    return solution_;
  }

  //===========================================================================
  shared_ptr<SplineVolume> VolSolution::getGeometryVolume() const
  //===========================================================================
  {
    return parent_->volume();
  }


  //===========================================================================
  void VolSolution::setMinimumDegree(int degree)
  //===========================================================================
  {
    int order = degree + 1;

    int raise_u = max(parent_->volume()->order(0), order) - solution_->order(0);
    int raise_v = max(parent_->volume()->order(1), order) - solution_->order(1);
    int raise_w = max(parent_->volume()->order(2), order) - solution_->order(2);

    raise_u = max(raise_u, 0);
    raise_v = max(raise_v, 0);
    raise_w = max(raise_w, 0);

    if (raise_u > 0 || raise_v > 0 || raise_w > 0)
      {
	solution_->raiseOrder(raise_u, raise_v, raise_w);
	updateConditions();
      }
  }

  //===========================================================================
  void refineToGeometry(int pardir)
  //===========================================================================
  {
    MESSAGE("refineToGeometry() not implemented");
  }

  //===========================================================================
  int VolSolution::nmbCoefs() const
  //===========================================================================
  {
    return solution_->numCoefs(0) * solution_->numCoefs(1) *
      solution_->numCoefs(2);
  }

  //===========================================================================
  int VolSolution::nmbCoefs(int pardir) const
  //===========================================================================
  {
    return solution_->numCoefs(pardir);
  }

  //===========================================================================
  int VolSolution::degree(int pardir) const
  //===========================================================================
  {
    int order = solution_->order(pardir);
    return order - 1;
  }

  //===========================================================================
  vector<double> VolSolution::knots(int pardir) const
  //===========================================================================
  {
    const BsplineBasis bas = solution_->basis(pardir);
    vector<double> result(bas.order() + bas.numCoefs());
    vector<double>::const_iterator bas_begin = bas.begin();
    vector<double>::const_iterator bas_end = bas.end();
    copy (bas_begin, bas_end, result.begin());
    return result;
  }

  //===========================================================================
  vector<double> VolSolution::distinctKnots(int pardir) const
  //===========================================================================
  {
    const BsplineBasis bas = solution_->basis(pardir);
    vector<double> result;
    bas.knotsSimple(result);
    return result;
  }

  //===========================================================================
  BsplineBasis VolSolution::basis(int pardir) const
  //===========================================================================
  {
    return solution_->basis(pardir);
  }

  //===========================================================================
  int VolSolution::dimension() const
  //===========================================================================
  {
    return solution_->dimension();
  }

  //===========================================================================
  void VolSolution::updateConditions()
  //===========================================================================
  {
      for (int i = 0; i < (int)boundary_conditions_.size(); ++i)
	boundary_conditions_[i]->update();
  }

  //===========================================================================
  double VolSolution::getGaussParameter(int index_of_Gauss_point, int pardir) const
  //===========================================================================
  {
    MESSAGE("getGaussParameter() not implemented");
  }

  //===========================================================================
  tpTolerances VolSolution::getTolerances() const
  //===========================================================================
  {
    return parent_->getTolerances();
  }

  //===========================================================================
  void VolSolution::neighbourInfo(BlockSolution* other, vector<int>& faces, vector<int>& faces_other,
				  vector<int>& orientation, vector<bool>& same_dir_order,
				  vector<bool>& space_matches) const
  //===========================================================================
  {
    MESSAGE("neighbourInfo() under construction");

    double tol = getTolerances().gap;

    VolSolution* vol_other = other->asVolSolution();
    parent_->getNeighbourInfo(vol_other->parent_, faces, faces_other, orientation, same_dir_order);

    space_matches.resize(faces.size());
    for (int i = 0; i < (int)faces.size(); ++i)
      {
	int bas_u = (faces[i] < 2) ? 1 : 0;
	int bas_v = (faces[i] < 2) ? 2 : 1;
	BsplineBasis basis_1 = solution_->basis(bas_u);
	BsplineBasis basis_2 = solution_->basis(bas_v);

	int bas_o_u = (faces_other[i] < 2) ? 1 : 0;
	int bas_o_v = (faces_other[i] < 2) ? 2 : 1;
	BsplineBasis basis_o_1 = vol_other->basis(bas_o_u);
	BsplineBasis basis_o_2 = vol_other->basis(bas_o_v);

	// If the basis are flipped we swap.
	if (!same_dir_order[i])
	  std::swap(basis_o_1, basis_o_2);

	int orient = orientation[i];
	// We then fix the direction of the basises.
	if ((bas_u == 0 && (orient == 1) || (orient == 4) || (orient == 5) || (orient == 7)) ||
	    (bas_u == 1 && (orient == 2) || (orient == 4) || (orient == 6) || (orient == 7)))
	  {
	    basis_o_1.reverseParameterDirection();
	  }
	if ((bas_v == 1 && (orient == 2) || (orient == 4) || (orient == 6) || (orient == 7)) ||
	    (bas_v == 2 && (orient == 3) || (orient == 5) || (orient == 6) || (orient == 7)))
	  {
	    basis_o_2.reverseParameterDirection();
	  }

	// And the domain.
	basis_o_1.rescale(basis_1.startparam(), basis_1.endparam());
	basis_o_2.rescale(basis_2.startparam(), basis_2.endparam());

	bool match1 = basis_1.sameSplineSpace(basis_o_1, tol);
	bool match2 = basis_2.sameSplineSpace(basis_o_2, tol);

	space_matches[i] = (match1 && match2);
      }
  }

}   // namespace Go
