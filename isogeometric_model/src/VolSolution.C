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
#include "GoTools/isogeometric_model/VolBoundaryCondition.h"


using std::max;
using std::shared_ptr;


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
			    vector<pair<Point, Point> >& polygon)
  //===========================================================================
  {
    MESSAGE("addBoundaryCondition() not implemented");
    // The polygon defines the domain.
#if 0
    boundary_conditions_.push_back(shared_ptr<VolBoundaryCondition>
				   (new VolBoundaryCondition(face_nmb, type, fbd, polygon, this)));
#endif
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
    MESSAGE("getNmbOfBoundaryConditions() not implemented");
    return 0;
  }

  //===========================================================================
  shared_ptr<VolBoundaryCondition> VolSolution::getBoundaryCondition(int index) const
  //===========================================================================
  {
    MESSAGE("getBoundaryCondition() not implemented");
    shared_ptr<VolBoundaryCondition> vbc;
    return vbc;
  }

  //===========================================================================
  void VolSolution::getFaceBoundaryConditions(int face_number, 
				 vector<shared_ptr<VolBoundaryCondition> >& bd_cond) const
  //===========================================================================
  {
    MESSAGE("getFaceBoundaryConditions() not implemented");
    // MESSAGE("getEdgeBoundaryConditions() not implemented");
  }

  //===========================================================================
  void VolSolution::getFaceBoundaryConditions(vector<shared_ptr<VolBoundaryCondition> >& bd_cond) const
  //===========================================================================
  {
    MESSAGE("getFaceBoundaryConditions() not implemented");
    // MESSAGE("getEdgeBoundaryConditions() not implemented");
  }

  //===========================================================================
  int VolSolution::getNmbOfPointBdConditions() const
  //===========================================================================
  {
    MESSAGE("getNmbOfPointBdConditions() not implemented");
    return 0;
  }

  //===========================================================================
  shared_ptr<VolPointBdCond> VolSolution::getPointBdCondition(int index) const
  //===========================================================================
  {
    MESSAGE("getPointBdCondition() not implemented");
    shared_ptr<VolPointBdCond> vpbc;
    return vpbc;
  }

  //===========================================================================
  void VolSolution::getFacePointBdConditions(int face_number, 
				vector<shared_ptr<VolPointBdCond> >& bd_cond) const
  //===========================================================================
  {
    MESSAGE("getFacePointBdConditions() not implemented");
  }

  //===========================================================================
  void VolSolution::getPointBdCond(vector<shared_ptr<VolPointBdCond> >& bd_cond) const
  //===========================================================================
  {
    MESSAGE("getPointBdCond() not implemented");
  }

  //===========================================================================
  bool VolSolution::matchingSplineSpace(BlockSolution* other) const
  //===========================================================================
  {
    MESSAGE("matchingSplineSpace() not implemented");
    return false;
  }

  //===========================================================================
  void VolSolution::getMatchingCoefficients(BlockSolution* other,
			       vector<pair<int,int> >& enumeration,
			       int match_pos) const
  //===========================================================================
  {
    MESSAGE("getMatchingCoefficients() not implemented");
  }

  //===========================================================================
  void VolSolution::getBoundaryCoefficients(int boundary,
					    vector<int>& enumeration) const
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
    MESSAGE("erasePreEvaluatedBasisFunctions() not implemented");
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
    MESSAGE("getBasisFunctions() not implemented");
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
  double VolSolution::getJacobian(vector<int>& index_of_Gauss_point) const
  //===========================================================================
  {
    MESSAGE("getJacobian() not implemented");
    return 0.0;
  }

  //===========================================================================
  void VolSolution::valuesInGaussPoint(const vector<int>& index_of_Gauss_point,
			  vector<Point>& derivs) const
  //===========================================================================
  {
    MESSAGE("valuesInGaussPoint() not implemented");
  }

  //===========================================================================
  void VolSolution::setSolutionCoefficients(const vector<double>& coefs)
  //===========================================================================
  {
    MESSAGE("setSolutionCoefficients() not implemented");
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
  shared_ptr<SplineVolume> VolSolution::getSolutionVolume() const
  //===========================================================================
  {
    MESSAGE("getSolutionVolume() not implemented");
    shared_ptr<SplineVolume> sv;
    return sv;
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
  tpTolerances VolSolution::getTolerances() const
  //===========================================================================
  {
    MESSAGE("getTolerances() not implemented");
    tpTolerances tt(0.0, 0.0, 0.0, 0.0);
    return tt;
  }

}   // namespace Go
