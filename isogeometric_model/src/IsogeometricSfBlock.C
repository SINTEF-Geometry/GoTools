//===========================================================================
//
// File : IsogeometricSfBlock.C
//
// Created: Tue Mar  2 14:56:02 2010
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================




#include "GoTools/isogeometric_model/IsogeometricSfBlock.h"
#include "GoTools/geometry/GapRemoval.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/SurfaceTools.h"



using std::vector;
using std::pair;

namespace Go
{

  //===========================================================================
  IsogeometricSfBlock::IsogeometricSfBlock(IsogeometricModel* model,
					   shared_ptr<SplineSurface> geom_sf,
					   vector<int> solution_space_dimension,
					   int index):
    IsogeometricBlock(model),
    surface_(geom_sf),
    index_(index)
   //===========================================================================
  {
    int ncoefs = geom_sf->numCoefs_u() * geom_sf->numCoefs_v();
    vector<double> weights;
    bool rational = geom_sf->rational();
    if (rational)
      {
	int dim = geom_sf->dimension();
	weights.resize(ncoefs);
	vector<double>::const_iterator it = geom_sf->rcoefs_begin();
	for (int i = 0, pos = dim; i < ncoefs; ++i, pos += dim + 1)
	  weights[i] = it[pos];
      }

    vector<double> coefs;
    int nmb_sol = (int)solution_space_dimension.size();
    solution_.resize(nmb_sol);
    BsplineBasis bas_u = geom_sf->basis_u();
    BsplineBasis bas_v = geom_sf->basis_v();

    for (int i = 0; i < nmb_sol; ++i)
      {
	int dim = solution_space_dimension[i];
	shared_ptr<SplineSurface> sol_surface;
	if (rational)
	  {
	    coefs.resize((dim+1) * ncoefs, 0.0);
	    for (int i = 0, pos = dim; i < ncoefs; ++i, pos += dim + 1)
	      coefs[pos] = weights[i];
	  }
	else
	  coefs.resize(dim * ncoefs, 0.0);
	sol_surface =
	  shared_ptr<SplineSurface>(new SplineSurface(bas_u, bas_v, coefs.begin(), dim, rational));
	shared_ptr<SfSolution> sol(new SfSolution(this, sol_surface));
	solution_[i] = sol;
      }
  }

  //===========================================================================
  IsogeometricSfBlock::~IsogeometricSfBlock()
  //===========================================================================
  {
  }


  //===========================================================================
  IsogeometricSfBlock* IsogeometricSfBlock::asIsogeometricSfBlock()
  //===========================================================================
  {
    return this;
  }


  //===========================================================================
  void IsogeometricSfBlock::addNeighbour(shared_ptr<IsogeometricSfBlock> neighbour,
					 int edge_nmb_this, int edge_nmb_other, bool equ_orient)
  //===========================================================================
  {
    neighbours_[edge_nmb_this] = neighbour;
    neighb_edge_[edge_nmb_this] = edge_nmb_other;
    equal_orientation_[edge_nmb_this] = equ_orient;
  }


  //===========================================================================
  int IsogeometricSfBlock::nmbOfNeighbours() const
  //===========================================================================
  {
    int ngb_count = 0;
    for (int i = 0; i < 4; ++i)
      if (neighbours_[i].get() != 0)
	++ngb_count;
    return ngb_count;
  }


  //===========================================================================
  IsogeometricSfBlock* IsogeometricSfBlock::getNeighbour(int edge_nmb) const
  //===========================================================================
  {
    return neighbours_[edge_nmb].get();
  }


  //===========================================================================
  bool IsogeometricSfBlock::isNeighbour(IsogeometricBlock* other) const
  //===========================================================================
  {
    IsogeometricSfBlock* other_sf = other->asIsogeometricSfBlock();
    if (other_sf == 0)
      return false;

    for (int i = 0; i < 4; ++i)
      if (neighbours_[i].get() == other_sf)
	return true;

    return false;
  }


  //===========================================================================
  int IsogeometricSfBlock::nmbCoefs() const
  //===========================================================================
  {
    return surface_->numCoefs_u() * surface_->numCoefs_v();
  }


  //===========================================================================
  BsplineBasis IsogeometricSfBlock::basis(int pardir) const
  //===========================================================================
  {
    return surface_->basis(pardir);
  }


  //===========================================================================
  shared_ptr<SplineCurve> IsogeometricSfBlock::getGeomBoundaryCurve(int edge_number) const
  //===========================================================================
  {
    switch (edge_number)
      {
      case 0:
	return shared_ptr<SplineCurve>(surface_->constParamCurve(surface_->startparam_u(), false));
      case 1:
	return shared_ptr<SplineCurve>(surface_->constParamCurve(surface_->endparam_u(), false));
      case 2:
	return shared_ptr<SplineCurve>(surface_->constParamCurve(surface_->startparam_v(), true));
      case 3:
	return shared_ptr<SplineCurve>(surface_->constParamCurve(surface_->endparam_v(), true));
      }

    shared_ptr<SplineCurve> empty_curve;
    return empty_curve;
  }


  //===========================================================================
  double IsogeometricSfBlock::getParamOnBdCurve(int edge_number, const Point& position) const
  //===========================================================================
  {
    shared_ptr<SplineCurve> crv = getGeomBoundaryCurve(edge_number);
    double param, dist;
    Point pt;
    crv->closestPoint(position, crv->startparam(), crv->endparam(), param, pt, dist);
    if (dist >= getTolerances().gap)
      THROW("Point is not on boundary curve");
    return param;
  }


  //===========================================================================
  bool IsogeometricSfBlock::geomIsDegenerate(vector<int>& degen_bd, double epsge)
  //===========================================================================
  {
    bool found = false;

    for (int i = 0; i < 4; ++i)
      {
	shared_ptr<SplineCurve> crv = getGeomBoundaryCurve(i);
	if (crv->isDegenerate(epsge))
	  {
	    degen_bd.push_back(i);
	    found = true;
	  }
      }

    return found;
  }


  //===========================================================================
  void IsogeometricSfBlock::getDegenEnumeration(vector<int>& degen_bd, 
						vector<vector<int> >& enumeration,
						double epsge)
  //===========================================================================
  {
    degen_bd.resize(0);
    bool is_degen = geomIsDegenerate(degen_bd, epsge);
    enumeration.resize(degen_bd.size());
    if (!is_degen)
      return;

    for (int i = 0; i < (int)degen_bd.size(); ++i)
      SurfaceTools::getCoefEnumeration(surface_, degen_bd[i], enumeration[i]);
  }


  //===========================================================================
  bool IsogeometricSfBlock::geomIsPeriodic(int per[], double epsge)
  //===========================================================================
  {
    bool is_periodic = false;
    SplineSurface *surf = surface_.get();
    for (int i = 0; i < 2; ++i)
      {
	per[i] = GeometryTools::analyzePeriodicity(*surf, i, epsge);
	if (per[i] >= 0)
	  is_periodic = true;
      }

    return is_periodic;
  }


  //===========================================================================
  bool IsogeometricSfBlock::getPeriodicEnumeration(int pardir, vector<pair<int, int> >& enumeration)
  //===========================================================================
  {
    SplineSurface *surf = surface_.get();
    if (GeometryTools::analyzePeriodicity(*surf, pardir, getTolerances().gap) == -1)
      return false;

    vector<int> coefs_min, coefs_max;
    if (pardir == 0)
      {
	SurfaceTools::getCoefEnumeration(surface_, 2, coefs_min);
	SurfaceTools::getCoefEnumeration(surface_, 3, coefs_max);
      }
    else
      {
	SurfaceTools::getCoefEnumeration(surface_, 0, coefs_min);
	SurfaceTools::getCoefEnumeration(surface_, 1, coefs_max);
      }

    enumeration.resize(coefs_min.size());
    for (int i = 0; i < (int)coefs_min.size(); ++i)
      enumeration[i] = pair<int, int>(coefs_min[i], coefs_max[i]);

    return true;
  }


  //===========================================================================
  void IsogeometricSfBlock::refineGeometry(vector<double> newknots, int pardir)
  //===========================================================================
  {
    if (pardir == 0)
      surface_->insertKnot_u(newknots);
    else if (pardir == 1)
      surface_->insertKnot_v(newknots);
    else
      return;  // Bad parameter direction
    for (int i = 0; i < (int)solution_.size(); ++i)
      solution_[i]->refineToGeometry(pardir);
  }


  //===========================================================================
  void IsogeometricSfBlock::refineGeometry(const BsplineBasis& other_basis, int pardir)
  //===========================================================================
  {
    if (pardir < 0 || pardir > 1)
      return;  // Bad parameter direction

    bool order_changed = false;
    BsplineBasis geo_basis = surface_->basis(pardir);

    int order_geo = geo_basis.order();
    int order_other = other_basis.order();
    if (order_geo < order_other)
      {
	if (pardir == 0)
	  surface_->raiseOrder(order_other - order_geo, 0);
	else
	  surface_->raiseOrder(0, order_other - order_geo);
	order_changed = true;
	order_geo = order_other;
      }

    vector<double> knots_other;
    other_basis.knotsSimple(knots_other);
    int order_diff = order_geo - order_other;
    vector<double> new_knots;
    for (int i = 0; i < (int)knots_other.size(); ++i)
      {
	double knot_val = knots_other[i];
	int knots_needed = order_diff + other_basis.knotMultiplicity(knot_val) -
	    geo_basis.knotMultiplicity(knot_val);
	if (knots_needed > 0)
	  for (int j = 0; j < knots_needed; ++j)
	    new_knots.push_back(knot_val);
      }

    if (new_knots.size() > 0)
      refineGeometry(new_knots, pardir);
    else if (order_changed)
	for (int i = 0; i < (int)solution_.size(); ++i)
	solution_[i]->refineToGeometry(pardir);
  }


  //===========================================================================
  void IsogeometricSfBlock::increaseGeometryDegree(int new_degree, int pardir)
  //===========================================================================
  {
    int new_order = new_degree + 1;
    if (pardir == 0)
      {
	if (surface_->order_u() >= new_order)
	  return;
	surface_->raiseOrder(new_order - surface_->order_u(), 0);
      }
    else if (pardir == 1)
      {
	if (surface_->order_v() >= new_order)
	  return;
	surface_->raiseOrder(0, new_order - surface_->order_v());
      }
    else
      return;  // Bad parameter direction

    for (int i = 0; i < (int)solution_.size(); ++i)
      solution_[i]->refineToGeometry(pardir);
  }


  //===========================================================================
  void IsogeometricSfBlock::updateGeometry(shared_ptr<SplineCurve> new_boundary, int edge_number)
  //===========================================================================
  {
    // Find out what to do. Remember to update solution spaces to inlude the geometry space
    MESSAGE("updateGeometry() not implemented");
  }


  //===========================================================================
  void IsogeometricSfBlock::erasePreEvaluatedBasisFunctions()
  //===========================================================================
  {
      for (int i = 0; i < (int)solution_.size(); ++i)
      solution_[i]->erasePreEvaluatedBasisFunctions();
  }


  //===========================================================================
  int IsogeometricSfBlock::getNmbOfBoundaryConditions() const
  //===========================================================================
  {
    if (solution_.size() == 0)
      return 0;
    else
      return solution_[0]->getNmbOfBoundaryConditions();
  }


  //===========================================================================
  void IsogeometricSfBlock::getEdgeBoundaryConditions(int edge_number, 
						      vector<shared_ptr<SfBoundaryCondition> >& bd_cond) const
  //===========================================================================
  {
      for (int i = 0; i < (int)solution_.size(); ++i)
      solution_[i]->getEdgeBoundaryConditions(edge_number, bd_cond);
  }


  //===========================================================================
  int IsogeometricSfBlock::getNmbOfPointBdConditions() const
  //===========================================================================
  {
    if (solution_.size() == 0)
      return 0;
    else
      return solution_[0]->getNmbOfPointBdConditions();
  }


  //===========================================================================
  void IsogeometricSfBlock::getEdgePointBdConditions(int edge_number, 
						     vector<shared_ptr<SfPointBdCond> >& bd_cond) const
  //===========================================================================
  {
      for (int i = 0; i < (int)solution_.size(); ++i)
      solution_[i]->getEdgePointBdConditions(edge_number, bd_cond);
  }


  //===========================================================================
  shared_ptr<SfSolution> IsogeometricSfBlock::getSolutionSpace(int solution_index)
  //===========================================================================
  {
    return solution_[solution_index];
  }


  //===========================================================================
  shared_ptr<SplineSurface> IsogeometricSfBlock::surface() const
  //===========================================================================
  {
    return surface_;
  }


  //===========================================================================
  void IsogeometricSfBlock::setMinimumDegree(int degree, int solutionspace_idx)
  //===========================================================================
  {
      if (solutionspace_idx >= 0 && solutionspace_idx < (int)solution_.size())
      solution_[solutionspace_idx]->setMinimumDegree(degree);
  }


  //===========================================================================
  bool IsogeometricSfBlock::updateSolutionSplineSpace(int solutionspace_idx)
  //===========================================================================
  {
    double tol = getTolerances().gap;

    shared_ptr<SplineSurface> surface_this = getSolutionSpace(solutionspace_idx)->getSolutionSurface();

    // Test at each edge
    for (int i = 0; i < 4; ++i) // umin, umax, vmin, vmax
      {
	shared_ptr<IsogeometricSfBlock> neighbour = neighbours_[i];
	if (neighbour.get() == NULL)
	  continue;  // No neighbour at this edge

	bool const_u_this = i < 2;
	BsplineBasis basis_this_edge; // Basis along common edge.
	BsplineBasis basis_this_const; // Basis in const par dir.
	if (const_u_this)
	  {
	    basis_this_edge = surface_this->basis_v();
	    basis_this_const = surface_this->basis_u();
	  }
	else
	  {
	    basis_this_edge = surface_this->basis_u();
	    basis_this_const = surface_this->basis_v();
	  }

	double const_par_this;
	if (i == 0 || i == 2)
	  const_par_this = basis_this_const.startparam();
	else
	  const_par_this = basis_this_const.endparam();

	shared_ptr<SplineSurface> surface_neighbour = neighbour->getSolutionSpace(solutionspace_idx)->getSolutionSurface();

	bool const_u_neighbour = neighb_edge_[i] < 2;
	BsplineBasis basis_neighbour_edge_pre;
	BsplineBasis basis_neighbour_const;
	if (const_u_neighbour)
	  {
	    basis_neighbour_edge_pre = surface_neighbour->basis_v();
	    basis_neighbour_const = surface_neighbour->basis_u();
	  }
	else
	  {
	    basis_neighbour_edge_pre = surface_neighbour->basis_u();
	    basis_neighbour_const = surface_neighbour->basis_v();
	  }

	double const_par_neighbour;
	if (neighb_edge_[i] == 0 || neighb_edge_[i] == 2)
	  const_par_neighbour = basis_neighbour_const.startparam();
	else
	  const_par_neighbour = basis_neighbour_const.endparam();

	// Make copy, to avoid manipulation of original basis
	BsplineBasis basis_neighbour_edge = basis_neighbour_edge_pre;
	if (!equal_orientation_[i])
	  basis_neighbour_edge.reverseParameterDirection();
	basis_neighbour_edge.rescale(basis_this_edge.startparam(),
				     basis_this_edge.endparam());

	// Test if spline spaces are equal
	if (basis_this_edge.sameSplineSpace(basis_neighbour_edge, tol))
	  continue;

	// Spline spaces are not equal. Make them equal, and return

	// Get boundary curve on first surface
	shared_ptr<SplineCurve> bd_crv_this = shared_ptr<SplineCurve>(surface_this->constParamCurve(const_par_this, !const_u_this));
	shared_ptr<CurveOnSurface> curve_this(new CurveOnSurface(surface_this, bd_crv_this, const_u_this ? 1 : 2, const_par_this, i));

	// Get boundary curve on second surface
	shared_ptr<SplineCurve> bd_crv_neighbour =
	  shared_ptr<SplineCurve>(surface_neighbour->constParamCurve(const_par_neighbour, !const_u_neighbour));
	shared_ptr<CurveOnSurface> curve_neighbour(new CurveOnSurface(surface_neighbour, bd_crv_neighbour,
								      const_u_neighbour ? 1 : 2, const_par_neighbour, neighb_edge_[i]));

	// Get limiting parameters
	double start1 = const_u_this ? surface_this->startparam_v() : surface_this->startparam_u();
	double end1 = const_u_this ? surface_this->endparam_v() : surface_this->endparam_u();
	double start2 = const_u_neighbour ? surface_neighbour->startparam_v() : surface_neighbour->startparam_u();
	double end2 = const_u_neighbour ? surface_neighbour->endparam_v() : surface_neighbour->endparam_u();

	// Get corner points
	Point p_start, p_end;
	bd_crv_this->point(p_start, start1);
	bd_crv_this->point(p_end, end1);

	// Make uniform
	GapRemoval::removeGapSpline(surface_this, curve_this, start1, end1,
				    surface_neighbour, curve_neighbour, start2, end2,
				    p_start, p_end, tol, &(equal_orientation_[i]));
	return true;
      }

    return false;
  }


  //===========================================================================
  int IsogeometricSfBlock::nmbSolutionSpaces() const
  //===========================================================================
  {
      return (int)solution_.size();
  }


  //===========================================================================
  int IsogeometricSfBlock::getEdgeOrientation(shared_ptr<ParamCurve> crv, double tol)
  //===========================================================================
  {
    vector<Point> sf_corners(4);
    vector<Point> curve_corners(2);
    crv->point(curve_corners[0], crv->startparam());
    crv->point(curve_corners[1], crv->endparam());

    double u_min = surface_->startparam_u();
    double u_max = surface_->endparam_u();
    double v_min = surface_->startparam_v();
    double v_max = surface_->endparam_v();
    surface_->point(sf_corners[0], u_min, v_min);
    surface_->point(sf_corners[1], u_max, v_min);
    surface_->point(sf_corners[2], u_min, v_max);
    surface_->point(sf_corners[3], u_max, v_max);

    double tol2 = tol * tol;
    bool close[4][2];
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 2; ++j)
	close[i][j] = sf_corners[i].dist2(curve_corners[j]) < tol2;

    if (close[0][0] && close[2][1])
      return 0;
    if (close[2][0] && close[0][1])
      return 1;
    if (close[1][0] && close[3][1])
      return 2;
    if (close[3][0] && close[1][1])
      return 3;
    if (close[0][0] && close[1][1])
      return 4;
    if (close[1][0] && close[0][1])
      return 5;
    if (close[2][0] && close[3][1])
      return 6;
    if (close[3][0] && close[2][1])
      return 7;

    return -1;
  }

  //===========================================================================
  void IsogeometricSfBlock::getNeighbourInfo(IsogeometricSfBlock* other,
					     std::vector<int>& edges,
					     std::vector<int>& edges_other,
					     std::vector<bool>& equal_oriented)
  //===========================================================================
  {
    edges.resize(0);
    edges_other.resize(0);
    equal_oriented.resize(0);
    for (int i = 0; i < 4; ++i)
      if (neighbours_[i].get() == other)
	{
	  edges.push_back(i);
	  edges_other.push_back(neighb_edge_[i]);
	  equal_oriented.push_back(equal_orientation_[i]);
	}
  }


} // end namespace Go
