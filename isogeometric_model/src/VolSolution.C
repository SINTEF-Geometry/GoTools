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


namespace Go
{

  //===========================================================================
  VolSolution::VolSolution(IsogeometricVolBlock* parent, std::shared_ptr<SplineVolume> sol_vol):
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
    MESSAGE("increaseDegree() not implemented");
  }

  //===========================================================================
  void VolSolution::insertKnots(const vector<int>& knot_intervals, int pardir)
  //===========================================================================
  {
    MESSAGE("insertKnots() not implemented");
  }

  //===========================================================================
  void VolSolution::insertKnots(const vector<double>& knots, int pardir)
  //===========================================================================
  {
    MESSAGE("insertKnots() not implemented");
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
    MESSAGE("nmbCoefs() not implemented");
    return 0;
  }

  //===========================================================================
  int VolSolution::nmbCoefs(int pardir) const
  //===========================================================================
  {
    MESSAGE("nmbCoefs() not implemented");
    return 0;
  }

  //===========================================================================
  int VolSolution::degree(int pardir) const
  //===========================================================================
  {
    MESSAGE("degree() not implemented");
    return 0;
  }

  //===========================================================================
  vector<double> VolSolution::knots(int pardir) const
  //===========================================================================
  {
    MESSAGE("knots() not implemented");
    vector<double> k;
    return k;
  }

  //===========================================================================
  vector<double> VolSolution::distinctKnots(int pardir) const
  //===========================================================================
  {
    MESSAGE("distinctKnots() not implemented");
    vector<double> dk;
    return dk;
  }

  //===========================================================================
  BsplineBasis VolSolution::basis(int pardir) const
  //===========================================================================
  {
    MESSAGE("basis() not implemented");
    BsplineBasis b;
    return b;
  }

  //===========================================================================
  int VolSolution::dimension() const
  //===========================================================================
  {
    MESSAGE("dimension() not implemented");
    return 0;
  }

  //===========================================================================
  void VolSolution::updateConditions()
  //===========================================================================
  {
    MESSAGE("updateConditions() not implemented");
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
  tpTolerances VolSolution::getTolerances() const
  //===========================================================================
  {
    MESSAGE("getTolerances() not implemented");
    tpTolerances tt(0.0, 0.0, 0.0, 0.0);
    return tt;
  }

}   // namespace Go
