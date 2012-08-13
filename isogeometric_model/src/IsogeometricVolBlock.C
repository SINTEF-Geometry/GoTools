//===========================================================================
//
// File : IsogeometricVolBlock.C
//
// Created: Thu Mar  4 11:03:39 2010
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================



#include "GoTools/isogeometric_model/IsogeometricVolBlock.h"
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/trivariate/VolumeTools.h"
#include "GoTools/trivariate/SurfaceOnVolume.h"
#include "GoTools/trivariate/GapRemovalVolume.h"



using std::vector;
using std::pair;

namespace Go
{

  //===========================================================================
  IsogeometricVolBlock::IsogeometricVolBlock(IsogeometricModel* model,
					     shared_ptr<SplineVolume> geom_vol,
					     vector<int> solution_space_dimension,
					     int index)
    : IsogeometricBlock(model),
      volume_(geom_vol),
      index_(index)
  //===========================================================================
  {
    int ncoefs = geom_vol->numCoefs(0) * geom_vol->numCoefs(1) * geom_vol->numCoefs(2);
    vector<double> weights;
    bool rational = geom_vol->rational();
    if (rational)
      {
	int dim = geom_vol->dimension();
	weights.resize(ncoefs);
	vector<double>::const_iterator it = geom_vol->rcoefs_begin();
	for (int i = 0, pos = dim; i < ncoefs; ++i, pos += dim + 1)
	  weights[i] = it[pos];
      }

    vector<double> coefs;
    int nmb_sol = (int)solution_space_dimension.size();
    solution_.resize(nmb_sol);
    BsplineBasis bas_u = geom_vol->basis(0);
    BsplineBasis bas_v = geom_vol->basis(1);
    BsplineBasis bas_w = geom_vol->basis(2);

    for (int i = 0; i < nmb_sol; ++i)
      {
	int dim = solution_space_dimension[i];
	shared_ptr<SplineVolume> sol_volume;
	if (rational)
	  {
	    coefs.resize((dim+1) * ncoefs, 0.0);
	    for (int i = 0, pos = dim; i < ncoefs; ++i, pos += dim + 1)
	      coefs[pos] = weights[i];
	  }
	else
	  coefs.resize(dim * ncoefs, 0.0);
	sol_volume =
	  shared_ptr<SplineVolume>(new SplineVolume(bas_u, bas_v, bas_w,
						    coefs.begin(), dim, rational));
	shared_ptr<VolSolution> sol(new VolSolution(this, sol_volume));
	solution_[i] = sol;
      }
  }

  //===========================================================================
  IsogeometricVolBlock::~IsogeometricVolBlock()
  //===========================================================================
  {
  }


  //===========================================================================
  IsogeometricVolBlock* IsogeometricVolBlock::asIsogeometricVolBlock()
  //===========================================================================
  {
    return this;
  }

  //===========================================================================
  void IsogeometricVolBlock::addNeighbour(shared_ptr<IsogeometricVolBlock> neighbour,
					  int face_nmb_this,
					  int face_nmb_other,
					  int orientation,
					  bool same_dir_order)
  //===========================================================================
  {
    neighbours_[face_nmb_this] = neighbour;
    neighb_face_[face_nmb_this] = face_nmb_other;
    orientation_[face_nmb_this] = orientation;
    same_dir_order_[face_nmb_this] = same_dir_order;
  }

  //===========================================================================
  int IsogeometricVolBlock::nmbOfNeighbours() const
  //===========================================================================
  {
    int ngb_count = 0;
    for (int i = 0; i < 6; ++i)
      if (neighbours_[i].get() != 0)
	++ngb_count;
    return ngb_count;
  }

  //===========================================================================
  IsogeometricVolBlock* IsogeometricVolBlock::getNeighbour(int bd_nmb) const
  //===========================================================================
  {
    return neighbours_[bd_nmb].get();
  }

  //===========================================================================
  bool IsogeometricVolBlock::isNeighbour(IsogeometricBlock* other) const
  //===========================================================================
  {

    IsogeometricVolBlock* other_vol = other->asIsogeometricVolBlock();
    if (other_vol == 0)
      return false;

    for (int i = 0; i < 6; ++i) {
      if (neighbours_[i].get() == other_vol)
	return true;
    }
    return false;
  }

  //===========================================================================
  int IsogeometricVolBlock::nmbCoefs() const
  //===========================================================================
  {
    return volume_->numCoefs(0) * volume_->numCoefs(1) * volume_->numCoefs(2);
  }

  //===========================================================================
  BsplineBasis IsogeometricVolBlock::basis(int pardir) const
  //===========================================================================
  {
    return volume_->basis(pardir);
  }

  //===========================================================================
  shared_ptr<SplineSurface> IsogeometricVolBlock::getGeomBoundarySurface(int face_number) const
  //===========================================================================
  {
    /// Fetch one boundary surface
    shared_ptr<SplineSurface> ss = volume_->getBoundarySurface(face_number);

    return ss;
  }

  //===========================================================================
  vector<double> IsogeometricVolBlock::getParamOnBdSurf(int face_number, const Point& position) const
  //===========================================================================
  {
    shared_ptr<SplineSurface> srf = getGeomBoundarySurface(face_number);
    double clo_u, clo_v, clo_dist;
    Point clo_pt;
    double eps_geo = getTolerances().gap;

    srf->closestPoint(position, clo_u, clo_v, clo_pt, clo_dist, eps_geo);
    if (clo_dist >= getTolerances().gap)
      THROW("Point is not on boundary surface");

    vector<double> clo_par(2);
    clo_par[0] = clo_u;
    clo_par[1] = clo_v;

    return clo_par;
  }

  //===========================================================================
  bool IsogeometricVolBlock::geomIsDegenerate(vector<pair<int,int> >& degen_bd, 
					      double epsge)
  //===========================================================================
  {
    bool found = false;
    for (int i = 0; i < 6; ++i)
      {
	int type = 0;
	bool b, r, t, l;
	volume_->isDegenerate(i, type, b, r, t, l, epsge);
	if (type > 0)
	  {
	    degen_bd.push_back(std::make_pair(i, type));
	    found = true;
	  }
      }

    return found;
  }

  //===========================================================================
  void IsogeometricVolBlock::getDegenEnumeration(vector<pair<int,int> >& degen_bd, 
						 vector<vector<int> >& enumeration,
						 double epsge)
  //===========================================================================
  {
    degen_bd.resize(0);
    bool is_degen = geomIsDegenerate(degen_bd, epsge);
    enumeration.resize(degen_bd.size());
    if (!is_degen)
      return;

    for (int ki = 0; ki < (int)degen_bd.size(); ++ki)
    {
      shared_ptr<SplineSurface> ssrf =
	volume_->getBoundarySurface(degen_bd[ki].first);
      SurfaceTools::getCoefEnumeration(ssrf, degen_bd[ki].second, enumeration[ki]);
    }
  }

  //===========================================================================
  bool IsogeometricVolBlock::geomIsPeriodic(int per[], double epsge)
  //===========================================================================
  {
    bool is_periodic = false;
    SplineVolume *vol = volume_.get();
    for (int i = 0; i < 2; ++i)
      {
	per[i] = VolumeTools::analyzePeriodicity(*vol, i, epsge);

	if (per[i] >= 0)
	  is_periodic = true;
      }

    return is_periodic;
  }

  //===========================================================================
  bool IsogeometricVolBlock::getPeriodicEnumeration(int pardir,
						    vector<pair<int, int> >& enumeration)
  //===========================================================================
  {
    if (pardir < 0 || pardir > 2)
      THROW("Bad parameter direction."); //return false;  // Bad parameter direction

    SplineVolume *vol = volume_.get();
    if (VolumeTools::analyzePeriodicity(*vol, pardir, getTolerances().gap) == -1)
      return false;

    vector<int> coefs_min, coefs_max;
    int bd_min, bd_max;
    if (pardir == 0)
      {
	bd_min = 0;
	bd_max = 1;
      }
    if (pardir == 1)
      {
	bd_min = 2;
	bd_max = 3;
      }
    else
      {
	bd_min = 4;
	bd_max = 5;
      }

    VolumeTools::getVolCoefEnumeration(volume_, bd_min, coefs_min);
    VolumeTools::getVolCoefEnumeration(volume_, bd_max, coefs_max);

    enumeration.resize(coefs_min.size());
    for (int i = 0; i < (int)coefs_min.size(); ++i)
      enumeration[i] = pair<int, int>(coefs_min[i], coefs_max[i]);

    return true;
  }

  //===========================================================================
  void IsogeometricVolBlock::refineGeometry(vector<double> newknots, int pardir)
  //===========================================================================
  {
    if (pardir >= 0 && pardir <= 2)
      volume_->insertKnot(pardir, newknots);
    else
      THROW("Bad parameter direction."); //return;  // Bad parameter direction

    for (int i = 0; i < (int)solution_.size(); ++i)
      solution_[i]->insertKnots(newknots, pardir);
  }

  //===========================================================================
  void IsogeometricVolBlock::refineGeometry(const BsplineBasis& other_basis, int pardir)
  //===========================================================================
  {
    if (pardir < 0 || pardir > 2)
      THROW("Bad parameter direction.");//return;  // Bad parameter direction

    bool order_changed = false;
    BsplineBasis geo_basis = volume_->basis(pardir);

    int order_geo = geo_basis.order();
    int order_other = other_basis.order();
    if (order_geo < order_other)
      {
	int raise_order = order_other - order_geo;
	int raise_u = (pardir == 0) ? raise_order : 0;
	int raise_v = (pardir == 1) ? raise_order : 0;
	int raise_w = (pardir == 2) ? raise_order : 0;
	volume_->raiseOrder(raise_u, raise_v, raise_w);
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
	solution_[i]->increaseDegree(order_other - 1, pardir);
  }

  //===========================================================================
  void IsogeometricVolBlock::increaseGeometryDegree(int new_degree, int pardir)
  //===========================================================================
  {
    int new_order = new_degree + 1;
    if (pardir >= 0 && pardir <= 2)
      {
	if (volume_->order(pardir) >= new_order)
	  return;
	int raise_order = new_order - volume_->order(pardir);
	int raise_u = (pardir == 0) ? raise_order : 0;
	int raise_v = (pardir == 1) ? raise_order : 0;
	int raise_w = (pardir == 2) ? raise_order : 0;
	volume_->raiseOrder(raise_u, raise_v, raise_w);
      }
    else
      THROW("Bad parameter direction.");//return;  // Bad parameter direction

    for (int i = 0; i < (int)solution_.size(); ++i)
      solution_[i]->increaseDegree(new_order, pardir);
  }

  //===========================================================================
  void IsogeometricVolBlock::updateGeometry(shared_ptr<SplineSurface> new_boundary,
					    int face_number)
  //===========================================================================
  {
    // Solution space must be updated to include the geometry space.
    MESSAGE("updateGeometry() not implemented");
  }

  //===========================================================================
  void IsogeometricVolBlock::erasePreEvaluatedBasisFunctions()
  //===========================================================================
  {
      for (int i = 0; i < (int)solution_.size(); ++i)
      solution_[i]->erasePreEvaluatedBasisFunctions();
  }

  //===========================================================================
  int IsogeometricVolBlock::getNmbOfBoundaryConditions() const
  //===========================================================================
  {
    if (solution_.size() == 0)
      return 0;
    else
      return solution_[0]->getNmbOfBoundaryConditions();
  }

  //===========================================================================
  shared_ptr<VolBoundaryCondition> IsogeometricVolBlock::getBoundaryCondition(int index) const
  //===========================================================================
  {
    if (solution_.size() == 0)
    {
      shared_ptr<VolBoundaryCondition> vbc;
      return vbc;
    }
    else
      return solution_[0]->getBoundaryCondition(index);
  }

  //===========================================================================
  void IsogeometricVolBlock::getFaceBoundaryConditions(int face_number, 
						       vector<shared_ptr<VolBoundaryCondition> >& bd_cond) const
  //===========================================================================
  {
      for (int i = 0; i < (int)solution_.size(); ++i)
      solution_[i]->getFaceBoundaryConditions(face_number, bd_cond);
  }

  //===========================================================================
  int IsogeometricVolBlock::getNmbOfPointBdConditions() const
  //===========================================================================
  {
    if (solution_.size() == 0)
      return 0;
    else
      return solution_[0]->getNmbOfPointBdConditions();
  }
  
  //===========================================================================
  shared_ptr<VolPointBdCond> IsogeometricVolBlock::getPointBdCondition(int index) const
  //===========================================================================
  {
    if (solution_.size() == 0)
    {
      shared_ptr<VolPointBdCond> vpbc;
      return vpbc;
    }
    else
      return solution_[0]->getPointBdCondition(index);
  }

  //===========================================================================
  void
  IsogeometricVolBlock::getFacePointBdConditions(int face_number, 
						 vector<shared_ptr<VolPointBdCond> >& bd_cond) const
  //===========================================================================
  {
      for (int i = 0; i < (int)solution_.size(); ++i)
      solution_[i]->getFacePointBdConditions(face_number, bd_cond);
  }


  //===========================================================================
  shared_ptr<VolSolution> IsogeometricVolBlock::getSolutionSpace(int solution_index)
  //===========================================================================
  {
      if (solution_index < (int)solution_.size())
      return solution_[solution_index];
    else
      {
	MESSAGE("No solution exists!");
	shared_ptr<VolSolution> vs_dummy;
	return vs_dummy;
      }
  }


  //===========================================================================
  shared_ptr<SplineVolume> IsogeometricVolBlock::volume() const
  //===========================================================================
  {
    return volume_;
  }


  //===========================================================================
  void IsogeometricVolBlock::setMinimumDegree(int degree, int solutionspace_idx)
  //===========================================================================
  {
    if (solutionspace_idx >= 0 && solutionspace_idx < (int)solution_.size())
      solution_[solutionspace_idx]->setMinimumDegree(degree);
  }


  //===========================================================================
  bool IsogeometricVolBlock::updateSolutionSplineSpace(int solutionspace_idx)
  //===========================================================================
  {
    MESSAGE("updateSolutionSplineSpace() not yet tested!");
    return false;

    double tol = getTolerances().gap;

    shared_ptr<SplineVolume> volume_this =
      getSolutionSpace(solutionspace_idx)->getSolutionVolume();

    // Test at each face
    for (int i = 0; i < 6; ++i) // umin, umax, vmin, vmax, wmin, wmax
      {
	shared_ptr<IsogeometricVolBlock> neighbour = neighbours_[i];
	if (neighbour.get() == NULL)
	  continue;  // No neighbour at this face

	// We locate the two pairs of matching basises (u- and v- for both sfs).

	int const_dir = i/2;
	BsplineBasis basis_this_1 = (const_dir == 0) ?
	  volume_this->basis(1) : volume_this->basis(0);
	BsplineBasis basis_this_2 = (const_dir == 2) ?
	  volume_this->basis(1) : volume_this->basis(2);
	BsplineBasis basis_this_const = volume_this->basis(const_dir);

	double const_par_this = (i%2 == 0) ?
	  basis_this_const.startparam() : basis_this_const.endparam();

	shared_ptr<SplineVolume> volume_neighbour =
	  neighbour->getSolutionSpace(solutionspace_idx)->getSolutionVolume();
	int const_dir_neighbour = neighb_face_[i]/2;
	// We must also map the varying basises between the volumes.
	int neighb_dir_1 = (const_dir_neighbour == 0) ? 1 : 0;
	BsplineBasis basis_neighbour_1_pre = volume_neighbour->basis(neighb_dir_1);
	// BsplineBasis basis_neighbour_2_pre = (const_dir_neighbour == 2) ?
	//   volume_neighbour->basis(1) : volume_neighbour->basis(2);
	int neighb_dir_2 = (const_dir_neighbour == 2) ? 1 : 2;
	BsplineBasis basis_neighbour_2_pre = volume_neighbour->basis(neighb_dir_2);
	BsplineBasis basis_neighbour_const = volume_this->basis(const_dir);

	double const_par_neighbour = (neighb_face_[i]%2 == 0) ?
	  basis_neighbour_const.startparam() : basis_neighbour_const.endparam();

	// If u-dir in a sf corr to v-dir in the other sf, we swap.
	// Note that u- and v- are referring to the parametrization of the boundary
	// surface, not the volume.
	if (!same_dir_order_[i])
	{
	    std::swap(basis_neighbour_1_pre, basis_neighbour_2_pre);
	    std::swap(neighb_dir_1, neighb_dir_2);
	}

	// Value in orientation.
	// Make copies, to avoid manipulation of original basis
	BsplineBasis basis_neighbour_1 = basis_neighbour_1_pre;
	bool opp_orientation_1 = (((neighb_dir_1 == 0) &&
				   ((orientation_[i] == 1) || (orientation_[i] == 4) ||
				    (orientation_[i] == 5) || (orientation_[i] == 7))) ||
				  ((neighb_dir_1 == 1) &&
				   ((orientation_[i] == 2) || (orientation_[i] == 4) ||
				    (orientation_[i] == 6) || (orientation_[i] == 7))) ||
				  ((neighb_dir_1 == 2) &&
				   ((orientation_[i] == 3) || (orientation_[i] == 5) ||
				    (orientation_[i] == 6) || (orientation_[i] == 7))));
	if (opp_orientation_1)
	    basis_neighbour_1.reverseParameterDirection();
	basis_neighbour_1.rescale(basis_this_1.startparam(), basis_this_1.endparam());
	BsplineBasis basis_neighbour_2 = basis_neighbour_2_pre;
	bool opp_orientation_2 = (((neighb_dir_2 == 0) &&
				   ((orientation_[i] == 1) || (orientation_[i] == 4) ||
				    (orientation_[i] == 5) || (orientation_[i] == 7))) ||
				  ((neighb_dir_2 == 1) &&
				   ((orientation_[i] == 2) || (orientation_[i] == 4) ||
				    (orientation_[i] == 6) || (orientation_[i] == 7))) ||
				  ((neighb_dir_2 == 2) &&
				   ((orientation_[i] == 3) || (orientation_[i] == 5) ||
				    (orientation_[i] == 6) || (orientation_[i] == 7))));
	if (opp_orientation_2)
	  basis_neighbour_2.reverseParameterDirection();
	basis_neighbour_2.rescale(basis_this_2.startparam(), basis_this_2.endparam());

	// Test if spline spaces are equal
	bool common_basis_1 = basis_this_1.sameSplineSpace(basis_neighbour_1, tol);
	bool common_basis_2 = basis_this_2.sameSplineSpace(basis_neighbour_2, tol);
	if (common_basis_1 && common_basis_2)
	  continue;

	// Spline spaces are not equal. Make them equal, and return

	// Get boundary surface on first volume
	shared_ptr<SplineSurface> bd_srf_this =
	  shared_ptr<SplineSurface>(volume_this->constParamSurface(const_par_this,
								   const_dir));
	shared_ptr<SurfaceOnVolume> surface_this
	  (new SurfaceOnVolume(volume_this, bd_srf_this,
			       const_dir, const_par_this, i, false));

	// Get boundary surface on second volume
	shared_ptr<SplineSurface> bd_srf_neighbour =
	  shared_ptr<SplineSurface>(volume_neighbour->constParamSurface
				  (const_par_neighbour, const_dir_neighbour));
	shared_ptr<SurfaceOnVolume> surface_neighbour
	  (new SurfaceOnVolume(volume_neighbour, bd_srf_neighbour,
			       const_dir_neighbour, const_par_neighbour,
			       neighb_face_[i], (opp_orientation_1 || opp_orientation_2)));

	// Get limiting parameters for the matching face parts.
	// @@sbr Here we are expecting rectangular match between the
	// faces. May be generalized in the future.
	RectDomain sf_rec_domain = surface_this->containingDomain();
	double sf1_start1 = sf_rec_domain.umin(); // const_u_this ? surface_this->startparam_v() :surface_this->startparam_u();
	double sf1_end1 = sf_rec_domain.umax();// const_u_this ? surface_this->endparam_v() :
//	    surface_this->endparam_u();
	double sf1_start2 =sf_rec_domain.vmin();// const_u_this ? surface_this->startparam_v() :
//	    surface_this->startparam_v();
	double sf1_end2 =sf_rec_domain.vmax();// const_u_this ? surface_this->endparam_v() :
	//surface_this->endparam_v();

	RectDomain sf2_rec_domain = surface_neighbour->containingDomain();
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
	GapRemoval::removeGapSpline(volume_this, surface_this,
				    sf1_start1, sf1_end1, sf1_start2, sf1_end2,
				    volume_neighbour, surface_neighbour,
				    sf2_start1, sf2_end1, sf2_start2, sf2_end2,
				    vertex_ll, vertex_ur, tol, orientation_[i]);
	return true;
      }

    return false;
  }


  //===========================================================================
  int IsogeometricVolBlock::nmbSolutionSpaces() const
  //===========================================================================
  {
      return (int)solution_.size();
  }


  //===========================================================================
  int IsogeometricVolBlock::getFaceOrientation(shared_ptr<ParamSurface> srf,
					       double tol)
  //===========================================================================
  {

#ifndef NDEBUG
      // We write to file the vol shell and input srf.
      std::ofstream debug("tmp/debug.g2");
      // Since we are lacking support for SplineVolume in viewer we
      // write the shell to file.
      // Sequence: u_min, u_max, v_min, v_max, w_min, w_max
      vector<shared_ptr<ParamSurface> > bd_faces =
	  volume_->getAllBoundarySurfaces();
      for (size_t ki = 0; ki < bd_faces.size(); ++ki)
      {
	  bd_faces[ki]->writeStandardHeader(debug);
	  bd_faces[ki]->write(debug);
      }
      srf->writeStandardHeader(debug);
      srf->write(debug);
#endif

    vector<pair<Point, Point> > srf_corners;
    srf->getCornerPoints(srf_corners); // (umin, vmin) and then ccw.

    vector<Point> vol_corners(8);
    double u_min = volume_->startparam(0);
    double u_max = volume_->endparam(0);
    double v_min = volume_->startparam(1);
    double v_max = volume_->endparam(1);
    double w_min = volume_->startparam(2);
    double w_max = volume_->endparam(2);
    volume_->point(vol_corners[0], u_min, v_min, w_min);
    volume_->point(vol_corners[1], u_max, v_min, w_min);
    volume_->point(vol_corners[2], u_min, v_max, w_min);
    volume_->point(vol_corners[3], u_max, v_max, w_min);
    volume_->point(vol_corners[4], u_min, v_min, w_max);
    volume_->point(vol_corners[5], u_max, v_min, w_max);
    volume_->point(vol_corners[6], u_min, v_max, w_max);
    volume_->point(vol_corners[7], u_max, v_max, w_max);

    double tol2 = tol * tol;
    bool close[8][4]; // In total we should get only 4 hits, for each srf corner.
    for (int i = 0; i < 8; ++i)
      for (int j = 0; j < 4; ++j)
	close[i][j] = vol_corners[i].dist2(srf_corners[j].first) < tol2;

    bool umin = (close[0][0] || close[0][1] || close[0][2] || close[0][3]) &&
      (close[2][0] || close[2][1] || close[2][2] || close[2][3]) &&
      (close[4][0] || close[4][1] || close[4][2] || close[4][3]) &&
      (close[6][0] || close[6][1] || close[6][2] || close[6][3]);
    bool umax = (close[1][0] || close[1][1] || close[4][2] || close[4][3]) &&
      (close[3][0] || close[3][1] || close[3][2] || close[3][3]) &&
      (close[5][0] || close[5][1] || close[5][2] || close[5][3]) &&
      (close[7][0] || close[7][1] || close[7][2] || close[7][3]);
    bool vmin = (close[0][0] || close[0][1] || close[0][2] || close[0][3]) &&
      (close[1][0] || close[1][1] || close[1][2] || close[1][3]) &&
      (close[4][0] || close[4][1] || close[4][2] || close[4][3]) &&
      (close[5][0] || close[5][1] || close[5][2] || close[5][3]);
    bool vmax = (close[2][0] || close[2][1] || close[2][2] || close[2][3]) &&
      (close[3][0] || close[3][1] || close[3][2] || close[3][3]) &&
      (close[6][0] || close[6][1] || close[6][2] || close[6][3]) &&
      (close[7][0] || close[7][1] || close[7][2] || close[7][3]);
    bool wmin = (close[0][0] || close[0][1] || close[0][2] || close[0][3]) &&
      (close[1][0] || close[1][1] || close[1][2] || close[1][3]) &&
      (close[2][0] || close[2][1] || close[2][2] || close[2][3]) &&
      (close[3][0] || close[3][1] || close[3][2] || close[3][3]);
    bool wmax = (close[4][0] || close[4][1] || close[4][2] || close[4][3]) &&
      (close[5][0] || close[5][1] || close[5][2] || close[5][3]) &&
      (close[6][0] || close[6][1] || close[6][2] || close[6][3]) &&
      (close[7][0] || close[7][1] || close[7][2] || close[7][3]);

    // We first set value denoting constant parameter.
    int return_val = -1;
    if (umin || umax)
      return_val = 0;
    else if (vmin || vmax)
      return_val = 4;
    else if (wmin || wmax)
      return_val = 8;

    // If we failed matching against volume edges something is wrong.
    if (return_val == -1)
    {
	MESSAGE("No face<->block match."); // @@sbr Message useful for debugging.
	return return_val; // Wrong input (or bug).
    }

    // We then check if we are at a start or end parameter.
    bool at_max = false;
    if (umax || vmax || wmax)
    {
	at_max = true;
	return_val += 2;
    }

    // Finally we check if the srf orientation matches that of the volume.
    // The orientation of the surface is given by the cross product of partial
    // derivatives, whilst the VolBlock boundary normals points into the block
    // for min iso values, outwards for max iso values.
    // u x v, v x w, w x u defines block normals (with u etc partial derivs),
    // i.e. the system (SplineVolume) is assumed to be right-handed.
    bool is_left_handed = volume_->isLeftHanded();
    if (is_left_handed) // @@sbr201110 Not expecting left-handed system to be ok.
	MESSAGE("System left handed, make sure it is supported!");
    vector<Point> pts(4);
    // We test in the (umin, umax) pt of srf.
    double upar = (close[0][0] || close[2][0] || close[4][0] || close[6][0]) ?
	volume_->startparam(0) : volume_->endparam(0);
    double vpar = (close[0][0] || close[1][0] || close[4][0] || close[5][0]) ?
	volume_->startparam(0) : volume_->endparam(0);
    double wpar = (close[0][0] || close[1][0] || close[2][0] || close[3][0]) ?
	volume_->startparam(0) : volume_->endparam(0);
    volume_->point(pts, upar, vpar, wpar, 1);
    Point vol_normal = (umin || umax) ? pts[1] : ((vmin || vmax) ? pts[2] : pts[3]);
    // We then compute the srf normal.
    Point sf_normal;
    srf->normal(sf_normal, srf_corners[0].second[0], srf_corners[0].second[1]);

    if (sf_normal*vol_normal < 0.0)
	return_val += 1; // Opposite normal.

    return return_val;
  }


  //===========================================================================
  void IsogeometricVolBlock::getNeighbourInfo(IsogeometricVolBlock* other,
					      std::vector<int>& faces,
					      std::vector<int>& faces_other,
					      std::vector<int>& orientation)
  //===========================================================================
  {
      MESSAGE("getNeighbourInfo() not implemented");
  }

} // end namespace Go
