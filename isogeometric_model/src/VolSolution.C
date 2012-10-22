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
#include "GoTools/trivariate/SurfaceOnVolume.h"
#include "GoTools/trivariate/GapRemovalVolume.h"


using std::max;
using std::pair;

class InsideInterval
{
public:
  InsideInterval(int deg, int basis_func_id)
    : deg_(deg), basis_func_id_(basis_func_id)
    {}

  bool operator()(int knot_ind) const
    {
      return (knot_ind - deg_ <= basis_func_id_ && basis_func_id_ < knot_ind + 1);
    }

private:
  int deg_;
  int basis_func_id_;
};

static bool
inside_interval(int deg, int basis_func_id, int knot_ind)
{
  return (knot_ind - deg <= basis_func_id && basis_func_id < knot_ind + 1);
}


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
  void VolSolution::addBoundaryCondition(int face_nmb, BdConditionType type,Point &const_val,
					 vector<pair<double, double> >& domain)
  //===========================================================================
  {
    boundary_conditions_.push_back(shared_ptr<VolBoundaryCondition>
				   (new VolBoundaryCondition(face_nmb, type, const_val, domain, this)));
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

    VolumeTools::getVolCoefEnumeration(solution_, faces[match_pos], coefs_this);
    VolumeTools::getVolCoefEnumeration(vol_other->solution_, faces_other[match_pos], coefs_other);
    // Expecting that the coef ordering is u-dir first, then v-dir, finally w-dir.

    int this_in1 = (faces[match_pos] < 2) ? solution_->numCoefs(1) : solution_->numCoefs(0);
    int this_in2 = (faces[match_pos] < 4) ? solution_->numCoefs(2) : solution_->numCoefs(1);
    // We check if the sfs (in the volume intersection) have corr u- and v-dir.
    bool dir_order_ok = parent_->sameDirOrder(faces[match_pos]);
    if (!dir_order_ok)
      {
	// We know that the two faces have matching coefs.
	// We transpose the matrix for the other volume.
	vector<int> coefs_other_tr(coefs_other.size());
	for (int ki = 0; ki < this_in2; ++ki)
	  for (int kj = 0; kj < this_in1; ++kj)
	    coefs_other_tr[ki*this_in1+kj] = coefs_other[kj*this_in2+ki];
	coefs_other = coefs_other_tr;
      }
    int const_dir = faces[match_pos]/2;
    int const_dir_other = faces_other[match_pos]/2;
    int orient = orientation[match_pos];
    bool u_rev = ((const_dir == 0 && (orient == 2 || orient == 4 || orient == 6 || orient == 7)) ||
		  (const_dir == 1 && (orient == 1 || orient == 4 || orient == 5 || orient == 7)) ||
		  (const_dir == 2 && (orient == 1 || orient == 4 || orient == 5 || orient == 7)));
    bool v_rev = ((const_dir == 1 && (orient == 3 || orient == 5 || orient == 6 || orient == 7)) ||
		  (const_dir == 0 && (orient == 3 || orient == 5 || orient == 6 || orient == 7)) ||
		  (const_dir == 2 && (orient == 2 || orient == 4 || orient == 6 || orient == 7)));

    // We make sure the coefs for faces_other match those of faces.
    int in2_start = (v_rev) ? this_in2 - 1 : 0;
    int in2_step = (v_rev) ? -1 : 1;
    int in1_start = (u_rev) ? this_in1 - 1 : 0;
    int in1_step = (u_rev) ? -1 : 1;
    for (int ki = 0; ki < this_in2; ++ki)
      {
	for (int kj = 0; kj < this_in1; ++kj)
	  {
	    int in = ki*this_in1 + kj;
	    int in_other = (in2_start + ki*in2_step)*this_in1 + in1_start + kj*in1_step;
	    enumeration.push_back(pair<int, int>(coefs_this[in], coefs_other[in_other]));
	  }
      }
  }

  //===========================================================================
  void VolSolution::getBoundaryCoefficients(int boundary,
					    vector<int>& enumeration) const
  //===========================================================================
  {
    VolumeTools::getVolCoefEnumeration(solution_, boundary, enumeration);
  }


  //===========================================================================
  void
  VolSolution::getBoundaryCoefficients(int boundary,
				       std::vector<int>& enumeration_bd,
				       std::vector<int>& enumeration_bd2) const
  //===========================================================================
  {
    VolumeTools::getVolCoefEnumeration(solution_, boundary,
				       enumeration_bd, enumeration_bd2);
  }


  //===========================================================================
  void VolSolution::makeMatchingSplineSpace(BlockSolution* other)
  //===========================================================================
  {
    double tol = getTolerances().gap;
    bool changed = false;

    vector<int> faces, faces_other;
    vector<int> orientation;
    vector<bool> same_dir_order;
    vector<bool> space_matches;
    neighbourInfo(other, faces, faces_other, orientation, same_dir_order, space_matches);

    shared_ptr<SplineVolume> solution_other = other->asVolSolution()->solution_;

    for (int i = 0; i < (int)faces.size(); ++i)
      if (!space_matches[i])
	{
	  // Spline spaces are not equal
	  // Get boundary curve on first solution surface
	  int const_dir = (faces[i]/2);
	  double const_par = (faces[i]%2 == 0) ?
	    solution_->startparam(const_dir) : solution_->endparam(const_dir);
	  shared_ptr<SplineSurface> bd_srf_this =
	    shared_ptr<SplineSurface>(solution_->constParamSurface(const_par, const_dir));
	  shared_ptr<SurfaceOnVolume> surface_this(new SurfaceOnVolume(solution_, bd_srf_this,
								       const_dir, const_par, faces[i], false));

	  int const_dir_other = (faces_other[i]/2);
	  double const_par_other = (faces_other[i]%2 == 0) ?
	    solution_other->startparam(const_dir_other) : solution_other->endparam(const_dir_other);
	  shared_ptr<SplineSurface> bd_srf_other =
	    shared_ptr<SplineSurface>(solution_other->constParamSurface(const_par_other, const_dir_other));
	  shared_ptr<SurfaceOnVolume> surface_other(new SurfaceOnVolume(solution_, bd_srf_this,
									const_dir_other, const_par_other, faces_other[i], false));

	  RectDomain sf_rec_domain = surface_this->containingDomain();
	  double sf1_start1 = sf_rec_domain.umin(); // const_u_this ? surface_this->startparam_v() :surface_this->startparam_u();
	  double sf1_end1 = sf_rec_domain.umax();// const_u_this ? surface_this->endparam_v() :
//	    surface_this->endparam_u();
	  double sf1_start2 =sf_rec_domain.vmin();// const_u_this ? surface_this->startparam_v() :
//	    surface_this->startparam_v();
	  double sf1_end2 =sf_rec_domain.vmax();// const_u_this ? surface_this->endparam_v() :
	  //surface_this->endparam_v();

	  RectDomain sf2_rec_domain = surface_other->containingDomain();
	  double sf2_start1 = sf2_rec_domain.umin();//const_u_neighbour ? surface_neighbour->startparam_v() :
	  //surface_neighbour->startparam_u();
	  double sf2_end1 = sf2_rec_domain.umax();//const_u_neighbour ? surface_neighbour->endparam_v() :
	  //surface_neighbour->endparam_u();
	  double sf2_start2 = sf2_rec_domain.vmin();//const_u_neighbour ? surface_neighbour->startparam_v() :
	  //  surface_neighbour->startparam_v();
	  double sf2_end2 = sf2_rec_domain.vmax();//const_u_neighbour ? surface_neighbour->endparam_v() :
	  //surface_neighbour->endparam_v();


	  // Get corner points
	  Point vertex_ll = surface_this->ParamSurface::point(sf1_start1, sf1_start2);
	  Point vertex_ur = surface_this->ParamSurface::point(sf1_end1, sf1_end2);
	  // Make uniform
	  GapRemoval::removeGapSpline(solution_, surface_this,
				      sf1_start1, sf1_end1, sf1_start2, sf1_end2,
				      solution_other, surface_other,
				      sf2_start1, sf2_end1, sf2_start2, sf2_end2,
				      vertex_ll, vertex_ur, tol, orientation[i]);

	  changed = true;
	}

    if (changed)
      updateConditions();
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
    ASSERT (Gauss_par.size() == 3);

    erasePreEvaluatedBasisFunctions();
    evaluated_grid_ = shared_ptr<preEvaluationVol>(new preEvaluationVol);

    vector<double> par_u = Gauss_par[0];
    vector<double> par_v = Gauss_par[1];
    vector<double> par_w = Gauss_par[2];
    int nmb_par_u = (int)par_u.size();
    int nmb_par_v = (int)par_v.size();
    int nmb_par_w = (int)par_w.size();

    // int num_u = solution_->numCoefs_u();
    // int num_v = solution_->numCoefs_v();
    int ord_u = solution_->order(0);
    int ord_v = solution_->order(1);
    int ord_w = solution_->order(2);

    evaluated_grid_->gauss_par1_.resize(nmb_par_u);
    copy(par_u.begin(), par_u.end(), evaluated_grid_->gauss_par1_.begin());
    evaluated_grid_->gauss_par2_.resize(nmb_par_v);
    copy(par_v.begin(), par_v.end(), evaluated_grid_->gauss_par2_.begin());
    evaluated_grid_->gauss_par3_.resize(nmb_par_w);
    copy(par_w.begin(), par_w.end(), evaluated_grid_->gauss_par3_.begin());

    evaluated_grid_->basisvals_u_.resize(nmb_par_u * ord_u * 2);
    evaluated_grid_->basisvals_v_.resize(nmb_par_v * ord_v * 2);
    evaluated_grid_->basisvals_w_.resize(nmb_par_w * ord_w * 2);
    evaluated_grid_->left_u_.resize(nmb_par_u);
    evaluated_grid_->left_v_.resize(nmb_par_v);
    evaluated_grid_->left_w_.resize(nmb_par_w);

    solution_->basis(0).computeBasisValues(&par_u[0], &par_u[0]+nmb_par_u,
					    &(evaluated_grid_->basisvals_u_[0]),
					    &(evaluated_grid_->left_u_[0]), 1);
    solution_->basis(1).computeBasisValues(&par_v[0], &par_v[0]+nmb_par_v,
					    &(evaluated_grid_->basisvals_v_[0]),
					    &(evaluated_grid_->left_v_[0]), 1);
    solution_->basis(2).computeBasisValues(&par_w[0], &par_w[0]+nmb_par_w,
					    &(evaluated_grid_->basisvals_w_[0]),
					    &(evaluated_grid_->left_w_[0]), 1);
    getGeometryVolume()->gridEvaluator(par_u, par_v, par_w,
				       evaluated_grid_->points_,
				       evaluated_grid_->deriv_u_,
				       evaluated_grid_->deriv_v_,
				       evaluated_grid_->deriv_w_);
  }

  //===========================================================================
  void VolSolution::getBasisFunctions(int index_of_Gauss_point1,
				      int index_of_Gauss_point2,
				      int index_of_Gauss_point3,
				      vector<double>& basisValues,
				      vector<double>& basisDerivs_u,
				      vector<double>& basisDerivs_v,
				      vector<double>& basisDerivs_w) const
				      // shared_ptr<BasisDerivs> result) const
  //===========================================================================
  {
    if (evaluated_grid_.get() == NULL)
      return;
    if (index_of_Gauss_point1 < 0 ||
	index_of_Gauss_point1 >= (int)evaluated_grid_->gauss_par1_.size() ||
	index_of_Gauss_point2 < 0 ||
	index_of_Gauss_point2 >= (int)evaluated_grid_->gauss_par2_.size() ||
	index_of_Gauss_point3 < 0 ||
	index_of_Gauss_point3 >= (int)evaluated_grid_->gauss_par3_.size())
      return;

    // The basis values are already computed, we skip ahead to accumulation.
    solution_->computeBasis(evaluated_grid_->basisvals_u_.begin() 
			    + 2 * index_of_Gauss_point1 * solution_->order(0),
			    evaluated_grid_->basisvals_v_.begin() 
			    + 2 * index_of_Gauss_point2 * solution_->order(1),
			    evaluated_grid_->basisvals_w_.begin() 
			    + 2 * index_of_Gauss_point3 * solution_->order(2),
			    evaluated_grid_->left_u_[index_of_Gauss_point1],
			    evaluated_grid_->left_v_[index_of_Gauss_point2],
			    evaluated_grid_->left_w_[index_of_Gauss_point3],
			    basisValues,
			    basisDerivs_u,
			    basisDerivs_v,
			    basisDerivs_w);
  }

  //===========================================================================
  void VolSolution::getBasisFunctions(double param1,
				      double param2,
				      double param3,
//				      shared_ptr<BasisDerivs> result,
				      vector<double>& basisValues,				      
				      vector<double>& basisDerivs_u,
				      vector<double>& basisDerivs_v,
				      vector<double>& basisDerivs_w) const
  //===========================================================================
  {
    double param[3];
    param[0] = param1;
    param[1] = param2;
    param[3] = param3;
    // All the basis values must be computed, before they are accumulated.
    solution_->computeBasis(param, basisValues,
			    basisDerivs_u, basisDerivs_v, basisDerivs_w);
  }

  //===========================================================================
  void VolSolution::getBasisFunctionValues(int basis_func_id_u,
					   int basis_func_id_v,
					   int basis_func_id_w,
					   std::vector<int>& index_of_Gauss_points1,
					   std::vector<int>& index_of_Gauss_points2,
					   std::vector<int>& index_of_Gauss_points3,
					   vector<double>& basisValues,
					   vector<double>& basisDerivs_u,
					   vector<double>& basisDerivs_v,
					   vector<double>& basisDerivs_w) const
  //					   shared_ptr<BasisDerivs> result) const
  //===========================================================================
  {
    const int order_u = solution_->order(0);
    const int order_v = solution_->order(1);
    const int order_w = solution_->order(2);
    const int deg_u = order_u - 1;
    const int deg_v = order_v - 1;
    const int deg_w = order_w - 1;

    const int dim = solution_->dimension();

    // To speed things up we locate the first and last occurence of
    // the index points (using the fact that the elements are sorted).

    // Since the c++11 standard is not yet fully supported we create a
    // class (for the predicate) instead of using a lambda function.
    vector<int>::const_iterator first_u =
      std::find_if(evaluated_grid_->left_u_.begin(), evaluated_grid_->left_u_.end(),
    		   InsideInterval(deg_u, basis_func_id_u));
    // vector<int>::const_iterator first_u =
    //   std::find_if(evaluated_grid_->left_u_.begin(), evaluated_grid_->left_u_.end(),
    // 		   [deg_u, basis_func_id_u] (int knot_ind_u)
    // 		   { return (knot_ind_u - deg_u <= basis_func_id_u && basis_func_id_u < knot_ind_u + 1); }
    // 	);
    vector<int>::const_iterator last_u = first_u;
    while ((*last_u - deg_u <= basis_func_id_u) && (last_u < evaluated_grid_->left_u_.end()))
      ++last_u;
    int first_u_ind = first_u - evaluated_grid_->left_u_.begin();
    int last_u_ind = last_u - evaluated_grid_->left_u_.begin(); // I.e. one passed the last index.
    // We split based on ind values.
    vector<vector<int> > u_ind;



    // vector<int>::const_iterator first_v =
    //   std::find_if(evaluated_grid_->left_v_.begin(), evaluated_grid_->left_v_.end(),
    // 		   [deg_v, basis_func_id_v] (int knot_ind_v)
    // 		   { return (knot_ind_v - deg_v <= basis_func_id_v && basis_func_id_v < knot_ind_v + 1); }
    // 	);
    vector<int>::const_iterator first_v =
      std::find_if(evaluated_grid_->left_v_.begin(), evaluated_grid_->left_v_.end(),
    		   InsideInterval(deg_v, basis_func_id_v));
    vector<int>::const_iterator last_v = first_v;
    while ((*last_v - deg_v <= basis_func_id_v) && (last_v < evaluated_grid_->left_v_.end()))
      ++last_v;
    int first_v_ind = first_v - evaluated_grid_->left_v_.begin();
    int last_v_ind = last_v - evaluated_grid_->left_v_.begin();

    // vector<int>::const_iterator first_w =
    //   std::find_if(evaluated_grid_->left_w_.begin(), evaluated_grid_->left_w_.end(),
    // 		   [deg_w, basis_func_id_w] (int knot_ind_w)
    // 		   { return (knot_ind_w - deg_w <= basis_func_id_w && basis_func_id_w < knot_ind_w + 1); }
    // 	);
    vector<int>::const_iterator first_w =
      std::find_if(evaluated_grid_->left_w_.begin(), evaluated_grid_->left_w_.end(),
    		   InsideInterval(deg_w, basis_func_id_w));
    vector<int>::const_iterator last_w = first_w;
    while ((*last_w - deg_w <= basis_func_id_w) && (last_w < evaluated_grid_->left_w_.end()))
      ++last_w;
    int first_w_ind = first_w - evaluated_grid_->left_w_.begin();
    int last_w_ind = last_w - evaluated_grid_->left_w_.begin();

#if 0
    std::cout << "first_u_ind: " << first_u_ind << ", first_v_ind: " << first_v_ind << ", first_w_ind: " << first_w_ind << std::endl;
    std::cout << "last_u_ind: " << last_u_ind << ", last_v_ind: " << last_v_ind << ", last_w_ind: " << last_w_ind << std::endl;
#endif

    // We run through the evaluated_grid_ and compute basis values for
    // the Gauss points in the support of our basis function.
    for (size_t kk = first_w_ind; kk < last_w_ind; ++kk)
      {
	int local_ind_w = basis_func_id_w + deg_w - evaluated_grid_->left_w_[kk];
	for (size_t kj = first_v_ind; kj < last_v_ind; ++kj)
	  {
	    int local_ind_v = basis_func_id_v + deg_v - evaluated_grid_->left_v_[kj];
	    for (size_t ki = first_u_ind; ki < last_u_ind; ++ki)
	      {
		int local_ind_u = basis_func_id_u + deg_u - evaluated_grid_->left_u_[ki];

		// We add the contribution from the sf coef (and
		// weight for rational case).
		vector<double> local_basisValues;
		vector<double> local_basisDerivs_u;
		vector<double> local_basisDerivs_v;
		vector<double> local_basisDerivs_w;
		// Size of returned vectors local_...: kk1*kk2*kk3, where kk1 is order_u etc.
		// The basis values are already computed, function skips directly to accumulation.
		solution_->computeBasis(evaluated_grid_->basisvals_u_.begin()
					+ 2 * ki * order_u,
					evaluated_grid_->basisvals_v_.begin() 
					+ 2 * kj * order_v,
					evaluated_grid_->basisvals_w_.begin() 
					+ 2 * kk * order_w,
					evaluated_grid_->left_u_[ki],
					evaluated_grid_->left_v_[kj],
					evaluated_grid_->left_w_[kk],
					local_basisValues,
					local_basisDerivs_u,
					local_basisDerivs_v,
					local_basisDerivs_w);
		basisValues.insert(basisValues.end(),
				   local_basisValues.begin() +
				   (local_ind_w*order_v*order_u + local_ind_v*order_u + local_ind_u)*dim,
				   local_basisValues.begin() +
				   (local_ind_w*order_v*order_u + local_ind_v*order_u + local_ind_u + 1)*dim);
		basisDerivs_u.insert(basisDerivs_u.end(),
				     local_basisDerivs_u.begin() +
				     (local_ind_w*order_v*order_u + local_ind_v*order_u + local_ind_u)*dim,
				     local_basisDerivs_u.begin() +
				     (local_ind_w*order_v*order_u + local_ind_v*order_u + local_ind_u + 1)*dim);
		basisDerivs_v.insert(basisDerivs_v.end(),
				     local_basisDerivs_v.begin() +
				     (local_ind_w*order_v*order_u + local_ind_v*order_u + local_ind_u)*dim,
				     local_basisDerivs_v.begin() +
				     (local_ind_w*order_v*order_u + local_ind_v*order_u + local_ind_u + 1)*dim);
		basisDerivs_w.insert(basisDerivs_w.end(),
				     local_basisDerivs_w.begin() +
				     (local_ind_w*order_v*order_u + local_ind_w*order_u + local_ind_u)*dim,
				     local_basisDerivs_w.begin() +
				     (local_ind_w*order_v*order_u + local_ind_w*order_u + local_ind_u + 1)*dim);

		// Storing the index of the Gauss points.
		index_of_Gauss_points1.push_back((int)ki);
		index_of_Gauss_points2.push_back((int)kj);
		index_of_Gauss_points3.push_back((int)kk);

	      }
	  }
      }
  }


  //===========================================================================
  void VolSolution::getBasisFunctionValues(int basis_func_id_u,
					   int basis_func_id_v,
					   int basis_func_id_w,
					   int elem_ind_u,
					   int elem_ind_v,
					   int elem_ind_w,
					   std::vector<int>& index_of_Gauss_points1,
					   std::vector<int>& index_of_Gauss_points2,
					   std::vector<int>& index_of_Gauss_points3,
					   std::vector<double>& basisValues,
					   std::vector<double>& basisDerivs_u,
					   std::vector<double>& basisDerivs_v,
					   std::vector<double>& basisDerivs_w) const
  //===========================================================================
  {
    MESSAGE("getBasisFunctionValues(): under construction.");

    return;

    const int order_u = solution_->order(0);
    const int order_v = solution_->order(1);
    const int order_w = solution_->order(2);
    const int deg_u = order_u - 1;
    const int deg_v = order_v - 1;
    const int deg_w = order_w - 1;

    const int dim = solution_->dimension();

    // To speed things up we locate the first and last occurence of
    // the index points (using the fact that the elements are sorted).

    // Since the c++11 standard is not yet fully supported we create a
    // class (for the predicate) instead of using a lambda function.
    vector<int>::const_iterator first_u =
      std::find_if(evaluated_grid_->left_u_.begin(), evaluated_grid_->left_u_.end(),
    		   InsideInterval(deg_u, basis_func_id_u));
    vector<int>::const_iterator last_u = first_u;
    while ((*last_u - deg_u <= basis_func_id_u) && (last_u < evaluated_grid_->left_u_.end()))
      ++last_u;
    int first_u_ind = first_u - evaluated_grid_->left_u_.begin();
    int last_u_ind = last_u - evaluated_grid_->left_u_.begin(); // I.e. one passed the last index.

    vector<int>::const_iterator first_v =
      std::find_if(evaluated_grid_->left_v_.begin(), evaluated_grid_->left_v_.end(),
    		   InsideInterval(deg_v, basis_func_id_v));
    vector<int>::const_iterator last_v = first_v;
    while ((*last_v - deg_v <= basis_func_id_v) && (last_v < evaluated_grid_->left_v_.end()))
      ++last_v;
    int first_v_ind = first_v - evaluated_grid_->left_v_.begin();
    int last_v_ind = last_v - evaluated_grid_->left_v_.begin();

    vector<int>::const_iterator first_w =
      std::find_if(evaluated_grid_->left_w_.begin(), evaluated_grid_->left_w_.end(),
    		   InsideInterval(deg_w, basis_func_id_w));
    vector<int>::const_iterator last_w = first_w;
    while ((*last_w - deg_w <= basis_func_id_w) && (last_w < evaluated_grid_->left_w_.end()))
      ++last_w;
    int first_w_ind = first_w - evaluated_grid_->left_w_.begin();
    int last_w_ind = last_w - evaluated_grid_->left_w_.begin();

#if 0
    std::cout << "first_u_ind: " << first_u_ind << ", first_v_ind: " << first_v_ind << ", first_w_ind: " << first_w_ind << std::endl;
    std::cout << "last_u_ind: " << last_u_ind << ", last_v_ind: " << last_v_ind << ", last_w_ind: " << last_w_ind << std::endl;
#endif

    // We run through the evaluated_grid_ and compute basis values for
    // the Gauss points in the support of our basis function.
    for (size_t kk = first_w_ind; kk < last_w_ind; ++kk)
      {
	int local_ind_w = basis_func_id_w + deg_w - evaluated_grid_->left_w_[kk];
	for (size_t kj = first_v_ind; kj < last_v_ind; ++kj)
	  {
	    int local_ind_v = basis_func_id_v + deg_v - evaluated_grid_->left_v_[kj];
	    for (size_t ki = first_u_ind; ki < last_u_ind; ++ki)
	      {
		int local_ind_u = basis_func_id_u + deg_u - evaluated_grid_->left_u_[ki];

		// We add the contribution from the sf coef (and
		// weight for rational case).
		vector<double> local_basisValues;
		vector<double> local_basisDerivs_u;
		vector<double> local_basisDerivs_v;
		vector<double> local_basisDerivs_w;
		// Size of returned vectors local_...: kk1*kk2*kk3, where kk1 is order_u etc.
		// The basis values are already computed, function skips directly to accumulation.
		solution_->computeBasis(evaluated_grid_->basisvals_u_.begin()
					+ 2 * ki * order_u,
					evaluated_grid_->basisvals_v_.begin() 
					+ 2 * kj * order_v,
					evaluated_grid_->basisvals_w_.begin() 
					+ 2 * kk * order_w,
					evaluated_grid_->left_u_[ki],
					evaluated_grid_->left_v_[kj],
					evaluated_grid_->left_w_[kj],
					local_basisValues,
					local_basisDerivs_u,
					local_basisDerivs_v,
					local_basisDerivs_w);
		basisValues.insert(basisValues.end(),
				   local_basisValues.begin() +
				   (local_ind_w*order_v*order_u + local_ind_v*order_u + local_ind_u)*dim,
				   local_basisValues.begin() +
				   (local_ind_w*order_v*order_u + local_ind_v*order_u + local_ind_u + 1)*dim);
		basisDerivs_u.insert(basisDerivs_u.end(),
				     local_basisDerivs_u.begin() +
				     (local_ind_w*order_v*order_u + local_ind_v*order_u + local_ind_u)*dim,
				     local_basisDerivs_u.begin() +
				     (local_ind_w*order_v*order_u + local_ind_v*order_u + local_ind_u + 1)*dim);
		basisDerivs_v.insert(basisDerivs_v.end(),
				     local_basisDerivs_v.begin() +
				     (local_ind_w*order_v*order_u + local_ind_v*order_u + local_ind_u)*dim,
				     local_basisDerivs_v.begin() +
				     (local_ind_w*order_v*order_u + local_ind_v*order_u + local_ind_u + 1)*dim);
		basisDerivs_w.insert(basisDerivs_w.end(),
				     local_basisDerivs_w.begin() +
				     (local_ind_w*order_v*order_u + local_ind_w*order_u + local_ind_u)*dim,
				     local_basisDerivs_w.begin() +
				     (local_ind_w*order_v*order_u + local_ind_w*order_u + local_ind_u + 1)*dim);

		// Storing the index of the Gauss points.
		index_of_Gauss_points1.push_back((int)ki);
		index_of_Gauss_points2.push_back((int)kj);
		index_of_Gauss_points3.push_back((int)kk);

	      }
	  }
      }
  }


  //===========================================================================
  double VolSolution::getJacobian(vector<int>& index_of_Gauss_point) const
  //===========================================================================
  {
    ASSERT (index_of_Gauss_point.size() == 3);
    if (evaluated_grid_.get() == NULL)
      return 0.0;

    int dim = getGeometryVolume()->dimension();
    ASSERT (dim == 3);

    int pos = dim * (index_of_Gauss_point[2] * ((int)evaluated_grid_->gauss_par1_.size()*
						(int)evaluated_grid_->gauss_par2_.size()) +
		     index_of_Gauss_point[1] * ((int)evaluated_grid_->gauss_par1_.size()) +
		     index_of_Gauss_point[0]);

    // We first create the Jacobian matrix.
    double jac_mat[3][3];
    for (int ki = 0; ki < 3; ++ki)
    {
	jac_mat[0][ki] = evaluated_grid_->deriv_u_[pos+ki];
	jac_mat[1][ki] = evaluated_grid_->deriv_v_[pos+ki];
	jac_mat[2][ki] = evaluated_grid_->deriv_w_[pos+ki];
    }

    // We then compute the determinant.
    double det = jac_mat[0][0]*(jac_mat[1][1]*jac_mat[2][2] - jac_mat[1][2]*jac_mat[2][1]) -
	jac_mat[1][0]*(jac_mat[0][1]*jac_mat[2][2] - jac_mat[0][2]*jac_mat[2][1]) +
	jac_mat[2][0]*(jac_mat[0][1]*jac_mat[1][2] - jac_mat[0][2]*jac_mat[1][1]);

    return det;
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
		     index_of_Gauss_point[2] * (int)evaluated_grid_->gauss_par1_.size() * (int)evaluated_grid_->gauss_par2_.size() );
    derivs.resize(4);
    derivs[0] = Point(evaluated_grid_->points_.begin() + pos,
		      evaluated_grid_->points_.begin() + pos + dim);
    derivs[1] = Point(evaluated_grid_->deriv_u_.begin() + pos,
		      evaluated_grid_->deriv_u_.begin() + pos + dim);
    derivs[2] = Point(evaluated_grid_->deriv_v_.begin() + pos,
		      evaluated_grid_->deriv_v_.begin() + pos + dim);
    derivs[3] = Point(evaluated_grid_->deriv_w_.begin() + pos,
		      evaluated_grid_->deriv_w_.begin() + pos + dim);
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
  void VolSolution::refineToGeometry(int pardir)
  //===========================================================================
  {
    bool changed = false;

    BsplineBasis base_solution = solution_->basis(pardir);
    BsplineBasis base_geometry = parent_->volume()->basis(pardir);

    int order_solution = base_solution.order();
    int order_geometry = base_geometry.order();
    if (order_solution < order_geometry)
      {
	if (pardir == 0)
	  solution_->raiseOrder(order_geometry - order_solution, 0, 0);
	else if (pardir == 1)
	  solution_->raiseOrder(0, order_geometry - order_solution, 0);
	else
	  solution_->raiseOrder(0, 0, order_geometry - order_solution);
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
	solution_->insertKnot(pardir, new_knots);
	changed = true;
      }

    if(changed)
      updateConditions();
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
  void VolSolution::getGaussParameter(int index_of_Gauss_point1,
				      int index_of_Gauss_point2,
				      int const_dir,
				      double& par1, double& par2) const
  //===========================================================================
  {
    if (index_of_Gauss_point1 < 0 || index_of_Gauss_point2 < 0)
	return;
    if (const_dir < 0 || const_dir > 2)
	return;

    if (evaluated_grid_.get() == NULL)
      return;

    if (const_dir == 0)
    {
	if ((index_of_Gauss_point1 >= (int)evaluated_grid_->gauss_par2_.size()) ||
	    (index_of_Gauss_point2 >= (int)evaluated_grid_->gauss_par3_.size()))
	    return;
	else
	{
	    par1 = evaluated_grid_->gauss_par2_[index_of_Gauss_point1];
	    par2 = evaluated_grid_->gauss_par3_[index_of_Gauss_point2];
	    return;
	}
      }
    else if (const_dir == 1)
    {
	if ((index_of_Gauss_point1 >= (int)evaluated_grid_->gauss_par1_.size()) ||
	    (index_of_Gauss_point2 >= (int)evaluated_grid_->gauss_par3_.size()))
	    return;
	else
	{
	    par1 = evaluated_grid_->gauss_par1_[index_of_Gauss_point1];
	    par2 = evaluated_grid_->gauss_par3_[index_of_Gauss_point2];
	    return;
	}
      }
    else
    {
	if ((index_of_Gauss_point1 >= (int)evaluated_grid_->gauss_par1_.size()) ||
	    (index_of_Gauss_point2 >= (int)evaluated_grid_->gauss_par2_.size()))
	    return;
	else
	{
	    par1 = evaluated_grid_->gauss_par1_[index_of_Gauss_point1];
	    par2 = evaluated_grid_->gauss_par2_[index_of_Gauss_point2];
	    return;
	}
      }
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
    double tol = getTolerances().gap;

    VolSolution* vol_other = other->asVolSolution();
    parent_->getNeighbourInfo(vol_other->parent_, faces, faces_other, orientation, same_dir_order);

    space_matches.resize(faces.size());
    for (int i = 0; i < (int)faces.size(); ++i)
      {
	int bas_u = (faces[i] < 2) ? 1 : 0;
	int bas_v = (faces[i] < 4) ? 2 : 1;
	BsplineBasis basis_1 = solution_->basis(bas_u);
	BsplineBasis basis_2 = solution_->basis(bas_v);

	int bas_o_u = (faces_other[i] < 2) ? 1 : 0;
	int bas_o_v = (faces_other[i] < 4) ? 2 : 1;
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
