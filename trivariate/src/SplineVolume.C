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

#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineUtils.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/trivariate/VolumeTools.h"
#include "GoTools/geometry/Utils.h"

#include <iomanip>
#include <fstream>


using std::vector;


namespace Go
{


//===========================================================================
SplineVolume::~SplineVolume()

//===========================================================================
{
}

//===========================================================================
void SplineVolume::read (std::istream& is)
//===========================================================================
{
    // We verify that the object is valid.
    bool is_good = is.good();
    if (!is_good) {
	THROW("Invalid geometry file!");
    }
    // Canonical data
    is >> dim_;
    is >> rational_;
    is >> basis_u_;
    is >> basis_v_;
    is >> basis_w_;
    int nc = basis_u_.numCoefs()*basis_v_.numCoefs()*basis_w_.numCoefs();
    is_good = is.good();
    if (!is_good) {
	THROW("Invalid geometry file!");
    }
    if (rational_) {
	int n = nc * (dim_ + 1);
	rcoefs_.resize(n);
	for (int i = 0; i < n; ++i)
	    is >> rcoefs_[i];
	coefs_.resize(nc*dim_);
	updateCoefsFromRcoefs();
    } else {
	int n = nc*dim_;
	coefs_.resize(n);
	for (int i = 0; i < n; ++i)
	    is >> coefs_[i];
    }

    is_good = is.good();
    if (!is_good) {
	THROW("Invalid geometry file!");
    }
}


//===========================================================================
void SplineVolume::write (std::ostream& os) const
//===========================================================================
{
    os << std::setprecision(15);

    os << dim_ << ' ' << rational_ << '\n';
    os << basis_u_;
    os << basis_v_;
    os << basis_w_;
    int n = basis_u_.numCoefs()*basis_v_.numCoefs()*basis_w_.numCoefs();
    int kdim = dim_ + (rational_ ? 1 : 0);
    const vector<double>& co = rational_ ? rcoefs_ : coefs_;
    for (int i = 0; i < n; ++i) {
	os << co[i*kdim];
	for (int d = 1; d < kdim; ++d) {
	    os << ' ' << co[i*kdim + d];
	}
	os << '\n';
    }
    os << std::endl;
}


//===========================================================================
BoundingBox SplineVolume::boundingBox() const
//===========================================================================
{
    BoundingBox box;
    box.setFromArray(&coefs_[0], &coefs_[0] + coefs_.size(), dim_);
    return box;
}

//===========================================================================
int SplineVolume::dimension() const
//===========================================================================
{
    return dim_;
}


//===========================================================================
ClassType SplineVolume::instanceType() const
//===========================================================================
{
    return classType();
}


//===========================================================================
double SplineVolume::startparam(int i) const
//===========================================================================
{
  BsplineBasis b = basis(i);
  std::vector<double>::const_iterator knot = b.begin();
  return knot[b.order() - 1];
}


//===========================================================================
double SplineVolume::endparam(int i) const
//===========================================================================
{
  BsplineBasis b = basis(i);
  std::vector<double>::const_iterator knot = b.begin();
  return knot[b.numCoefs()];
}


//===========================================================================
const Array<double,6> SplineVolume::parameterSpan() const
//===========================================================================
{
  Array<double,6> pSpan;
  for (int i = 0; i < 3; ++i)
    {
      BsplineBasis b = basis(i);
      pSpan[2*i] = b.startparam();
      pSpan[2*i+1] = b.endparam();
    }
  return pSpan;
}


//===========================================================================
double SplineVolume::nextSegmentVal(int dir, double par, bool forward, double tol) const
//===========================================================================
{
  double startpar = startparam(dir);
  double endpar = endparam(dir);
  const BsplineBasis& bas = basis(dir);

  if (!forward && par <= startpar)
    return startpar;
  
  if (forward && par >= endpar)
    return endpar;

  std::vector<double>::const_iterator knot;
  if (forward) {
    par += fabs(tol);
    knot = std::upper_bound(bas.begin(),bas.end(),par);
    if (knot == bas.end())
      return endpar;
    else
      return *knot;
  }
  else {
    par -= fabs(tol);
    for (knot=bas.end()-1; knot>bas.begin(); --knot) {
      if (*knot < par)
	return *knot;
    }
    return *(bas.begin());
  }
}


//===========================================================================
void SplineVolume::replaceCoefficient(int ix, Point coef)
//===========================================================================
{
  ASSERT(dim_ == coef.dimension());
  vector<double>::iterator c1 = coefs_begin() + ix*dim_;
  for (int ki=0; ki<dim_; ++ki)
    c1[ki] = coef[ki];

  if (rational_)
    {
      vector<double>::iterator c2 = rcoefs_begin() + ix*dim_;
      for (int ki=0; ki<dim_; ++ki)
	c2[ki] = c1[ki]*c2[dim_];
    }
}

////===========================================================================
void SplineVolume::getWeights(std::vector<double>& weights) const
//===========================================================================
{
    int ncoefs = basis_u_.numCoefs()*basis_v_.numCoefs()*basis_w_.numCoefs();
    if ((int)weights.size() != ncoefs)
	weights.resize(ncoefs);
    if (rational_)
    {
	int ki, ki2;
	for (ki=0, ki2=0; ki<ncoefs; ++ki, ki2+=(dim_+1))
	    weights[ki] = rcoefs_[ki2+dim_];
    }
    else
	std::fill(weights.begin(), weights.end(), 1.0);
    
}


//===========================================================================
void SplineVolume::setParameterDomain(double u1, double u2, 
				      double v1, double v2,
				      double w1, double w2)
//===========================================================================
{
  basis_u_.rescale(u1, u2);
  basis_v_.rescale(v1, v2);
  basis_w_.rescale(w1, w2);
} 


//===========================================================================
void SplineVolume::removeKnot(int pardir, double tpar)
//===========================================================================
{
  int kdim = rational_ ? dim_+1 : dim_;

  if (pardir == 0)
    swapParameterDirection(0,2);
  else if (pardir == 1)
    swapParameterDirection(1,2);

  // Make a hypercurve from this volume
  SplineCurve cv(numCoefs(2), order(2), basis_w_.begin(),
		 activeCoefs().begin(), kdim*numCoefs(1)*numCoefs(0), false);

  // Remove the knot from the curve
  cv.removeKnot(tpar);
  // Remove the knot from the basis
  basis_w_.removeKnot(tpar);

  // Copy back the data
  if (!rational_)
    {
      coefs_.resize(numCoefs(0)*numCoefs(1)*numCoefs(2)*dim_);
      std::copy(cv.coefs_begin(), cv.coefs_end(), activeCoefs().begin());
    }
  else
    {
      rcoefs_.resize(numCoefs(0)*numCoefs(1)*numCoefs(2)*(dim_+1));
      std::copy(cv.coefs_begin(), cv.coefs_end(), activeCoefs().begin());
      updateCoefsFromRcoefs();
    }

  if (pardir == 0)
    swapParameterDirection(0,2);
  else if (pardir == 1)
    swapParameterDirection(1,2);

}


//===========================================================================
DirectionCone SplineVolume::tangentCone(int pardir) const
//===========================================================================
{
  ALWAYS_ERROR_IF(dim_ != 3, "Normal only defined in 3D");

  int nu = numCoefs(0);
  int nv = numCoefs(1);
  int nw = numCoefs(2);

  const double* start = &coefs_[0];
  Point tmp(dim_);

  int end_u = nu;
  int end_v = nv;
  int end_w = nw;
  int step;

  if (pardir == 0)
    {
      step = dim_;
      --end_u;
    }
  else if (pardir == 1)
    {
      step = nu*dim_;
      --end_v;
    }
  else
    {
      step = nu*nv*dim_;
      --end_w;
    }

  vector<double> coefs;
  for (int i = 0; i < end_u; ++i)
    for (int j = 0; j < end_v; ++j)
      for (int k = 0; k < end_w; ++k)
	{
	  int from = ((k*nv + j)*nu + i) * dim_;
	  for (int dd = 0; dd < dim_; ++dd, ++from) {
	    tmp[dd] = *(start + from) - *(start + from + step);
	  }
	  if (tmp.length2() > 0.0) // Ignore vectors of zero length
	    for (int dd = 0; dd < dim_; ++dd)
	      coefs.push_back(tmp[dd]);
	}

  DirectionCone cone;
  cone.setFromArray(&coefs[0], &coefs[0] + coefs.size(), dim_);
  return cone;
}



//===========================================================================
void SplineVolume::appendVolume(SplineVolume* vol, int join_dir,
				int cont, bool repar)
//===========================================================================
{
  ALWAYS_ERROR_IF(rational_,
		  "Trying to append to a rational volume. Not implemented yet");
  ALWAYS_ERROR_IF(vol->rational(),
		  "Trying to append a rational volume. Not implemented yet");

  int dir1 = (join_dir == 0) ? 1 : 0;
  int dir2 = (join_dir == 2) ? 1 : 2;

  // Represent as surfaces, to unify the other two basis
  vector<shared_ptr<SplineSurface> > volsAsSurface;
  volsAsSurface.push_back(VolumeTools::representVolumeAsSurface(*this,dir1,dir2));
  volsAsSurface.push_back(VolumeTools::representVolumeAsSurface(*vol,dir1,dir2));
  GeometryTools::unifySurfaceSplineSpace(volsAsSurface, DEFAULT_PARAMETER_EPSILON);

  // Represent as curves, to append
  shared_ptr<SplineVolume> new_volume
    = VolumeTools::representSurfaceAsVolume(*(volsAsSurface[0]), dir1, dir2, basis(join_dir), 0); 
  shared_ptr<SplineCurve> curve1 = VolumeTools::representVolumeAsCurve(*new_volume, join_dir);
  new_volume = VolumeTools::representSurfaceAsVolume(*(volsAsSurface[1]), dir1, dir2, vol->basis(join_dir), 0); 
  shared_ptr<SplineCurve> curve2 = VolumeTools::representVolumeAsCurve(*new_volume, join_dir);
  double dist;
  curve1->appendCurve(curve2.get(), cont, dist, repar);

  new_volume = VolumeTools::representCurveAsVolume(*curve1, join_dir,
						   volsAsSurface[0]->basis_u(),
						   volsAsSurface[0]->basis_v(),
						   0);
  basis_u_ = new_volume->basis(0);
  basis_v_ = new_volume->basis(1);
  basis_w_ = new_volume->basis(2);
  int ncoefs = (new_volume->numCoefs(0)) * (new_volume->numCoefs(1)) * (new_volume->numCoefs(2)) * dim_;
  coefs_.resize(ncoefs);
  copy(new_volume->coefs_begin(), new_volume->coefs_begin()+ncoefs, coefs_.begin());
}



//===========================================================================
void SplineVolume::reverseParameterDirection(int pardir)
//===========================================================================
{
  if (pardir == 0)
    swapParameterDirection(0,2);
  else if (pardir == 1)
    swapParameterDirection(1,2);

  // we will treat the volume coefficients as curve coefficients 
  // on a spline curve in high dimensionality
  int kdim = dim_ + (rational_ ? 1 : 0) ;
  int n = numCoefs(2);
  int high_dimension = numCoefs(0) * numCoefs(1) * kdim;
  std::vector<double>& co = rational_ ? rcoefs_ : coefs_;
  for (int i = 0; i < n/2; ++i) {
    for (int dd = 0; dd < high_dimension; ++dd) {
      std::swap(co[i * high_dimension + dd], co[(n - 1 - i) * high_dimension + dd]);
    }
  }
  basis_w_.reverseParameterDirection();

  if (pardir == 0)
    swapParameterDirection(0,2);
  else if (pardir == 1)
    swapParameterDirection(1,2);

  if (rational_) {
    updateCoefsFromRcoefs();
  }

  getBoundarySurfaces(true);

}


//===========================================================================
void SplineVolume::swapParameterDirection(int pardir1, int pardir2)
//===========================================================================
{
  if (pardir1 < 0 || pardir2 < 0 || pardir1 > 2 || pardir2 > 2 || pardir1 == pardir2)
    return;

  int kdim = dim_ + (rational_ ? 1 : 0);
  int old_i_coefs = numCoefs(0);
  int old_j_coefs = numCoefs(1);
  int old_k_coefs = numCoefs(2);
  double* actCoefs = &(activeCoefs()[0]);
  vector<double> oldCoefs(actCoefs, actCoefs + kdim * old_i_coefs * old_j_coefs * old_k_coefs);
  int old_pos = 0;

  if ((pardir1 == 0 && pardir2 == 1) || (pardir1 == 1 && pardir2 == 0))
    {    // Swap u and v
      for (int k = 0; k < old_k_coefs; ++k)
	for (int j = 0; j < old_j_coefs; ++j)
	  for (int i = 0; i < old_i_coefs; ++i)
	    for (int d = 0; d < kdim; d++)
	      actCoefs[((k*old_i_coefs + i)*old_j_coefs + j)*kdim + d] = oldCoefs[old_pos++];
      basis_u_.swap(basis_v_);
      std::swap(bd_sfs_[0], bd_sfs_[2]);
      std::swap(bd_sfs_[1], bd_sfs_[3]);
      std::swap(periodicity_info_[0], periodicity_info_[1]);
      std::swap(degen_[0], degen_[2]);
      std::swap(degen_[1], degen_[3]);
    }

  else if ((pardir1 == 0 && pardir2 == 2) || (pardir1 == 2 && pardir2 == 0))
    {    // Swap u and w
      for (int k = 0; k < old_k_coefs; ++k)
	for (int j = 0; j < old_j_coefs; ++j)
	  for (int i = 0; i < old_i_coefs; ++i)
	    for (int d = 0; d < kdim; d++)
	      actCoefs[((i*old_j_coefs + j)*old_k_coefs + k)*kdim + d] = oldCoefs[old_pos++];
      basis_u_.swap(basis_w_);
      std::swap(bd_sfs_[0], bd_sfs_[4]);
      std::swap(bd_sfs_[1], bd_sfs_[5]);
      std::swap(periodicity_info_[0], periodicity_info_[2]);
      std::swap(degen_[0], degen_[4]);
      std::swap(degen_[1], degen_[5]);
    }

  else if ((pardir1 == 1 && pardir2 == 2) || (pardir1 == 2 && pardir2 == 1))
    {    // Swap v and w
      for (int k = 0; k < old_k_coefs; ++k)
	for (int j = 0; j < old_j_coefs; ++j)
	  for (int i = 0; i < old_i_coefs; ++i)
	    for (int d = 0; d < kdim; d++)
	      actCoefs[((j*old_k_coefs + k)*old_i_coefs + i)*kdim + d] = oldCoefs[old_pos++];
      basis_v_.swap(basis_w_);
      std::swap(bd_sfs_[2], bd_sfs_[4]);
      std::swap(bd_sfs_[3], bd_sfs_[5]);
      std::swap(periodicity_info_[1], periodicity_info_[2]);
      std::swap(degen_[2], degen_[4]);
      std::swap(degen_[3], degen_[5]);
    }

  if (rational_)
    updateCoefsFromRcoefs();
}


//===========================================================================
double SplineVolume::knotSpan(int pardir, int iknot) const
{
  vector<double>::const_iterator knots = basis(pardir).begin();
  if (iknot < 0 || iknot >= numCoefs(pardir) + order(pardir) - 1)
    return 0.0;
  return knots[iknot+1] - knots[iknot];
}
//===========================================================================


//===========================================================================
SplineSurface* SplineVolume::constParamSurface(double parameter,
					       int pardir) const
//===========================================================================
{
  int kdim = dim_ + (rational_ ? 1 : 0);
  const std::vector<double>& co = rational_ ? rcoefs_ : coefs_;
  std::vector<double> huge_curve_coefs;
  int bas0_idx; // The volume parameter direction of the first parameter direction on plane
  int bas1_idx; // The volume parameter direction of the sexond parameter direction on plane
  int const_bas_idx; // The volume parameter direction where the surface is constant

  // Set up coefficients and basis in right order depending on constant parameter direction
  if (pardir == 0)
    {
      // Coefficients come in parameter directions 1 (innermost), 2, 0
      int n0 = numCoefs(0);
      int n1 = numCoefs(1);
      int n2 = numCoefs(2);
      huge_curve_coefs.reserve(n0*n1*n2*kdim);
      for (int i = 0; i < n0; ++i)
	for (int k = 0; k < n2; ++k)
	  for (int j = 0; j < n1; ++j)
	    for (int l = 0; l < kdim; ++l)
	      huge_curve_coefs.push_back(co[((k*n1+j)*n0+i)*kdim + l]);
      bas0_idx = 1;
      bas1_idx = 2;
      const_bas_idx = 0;
    }
  else if (pardir == 1)
    {
      // Coefficients come in parameter directions 0, 2, 1
      int n0 = numCoefs(0);
      int n1 = numCoefs(1);
      int n2 = numCoefs(2);
      huge_curve_coefs.reserve(n0*n1*n2*kdim);
      for (int j = 0; j < n1; ++j)
	for (int k = 0; k < n2; ++k)
	  for (int i = 0; i < n0; ++i)
	    for (int l = 0; l < kdim; ++l)
	      huge_curve_coefs.push_back(co[((k*n1+j)*n0+i)*kdim + l]);
      bas0_idx = 0;
      bas1_idx = 2;
      const_bas_idx = 1;
    }
  else
    {
      // pardir == 2
      // Coefficients come in parameter directions 0, 1, 2
      huge_curve_coefs = co;
      bas0_idx = 0;
      bas1_idx = 1;
      const_bas_idx = 2;
    }

  const BsplineBasis& bas0 = basis(bas0_idx); // Bspline basis along first parameter direction on plane
  const BsplineBasis& bas1 = basis(bas1_idx); // Bspline basis along second parameter direction on plane
  const BsplineBasis& const_bas = basis(const_bas_idx); // Bspline basis for the parameter where the surface is constant

  // Place coefficients in a curve in a high dimensional space
  int num = const_bas.numCoefs();
  int order = const_bas.order();
  std::vector<double>::const_iterator knotstart = const_bas.begin();
  std::vector<double>::const_iterator coefstart = huge_curve_coefs.begin();
  int dim = bas0.numCoefs() * bas1.numCoefs() * kdim;
  SplineCurve huge_curve(num, order, knotstart, coefstart, dim, false);

  // Evaluate to get new coefficients in the final surface
  std::vector<double> coefs_wanted(dim);
  Point p(&coefs_wanted[0], &coefs_wanted[0] + dim, false);
  huge_curve.point(p, parameter);

  // Create final surface
  int num0 = bas0.numCoefs();
  int order0 = bas0.order();
  std::vector<double>::const_iterator knotstart0 = bas0.begin();
  int num1 = bas1.numCoefs();
  int order1 = bas1.order();
  std::vector<double>::const_iterator knotstart1 = bas1.begin();

  SplineSurface* ss = new SplineSurface (num0, num1, order0, order1,
					 knotstart0, knotstart1, coefs_wanted.begin(),
					 dim_, rational_);
  return ss;
}


//===========================================================================
void SplineVolume::swap(SplineVolume& other)
//===========================================================================
{
    std::swap(dim_, other.dim_);
    std::swap(rational_, other.rational_);
    basis_u_.swap(other.basis_u_);
    basis_v_.swap(other.basis_v_);
    basis_w_.swap(other.basis_w_);
    coefs_.swap(other.coefs_);
    rcoefs_.swap(other.rcoefs_);
    for (int i = 0; i < SPLINE_VOLUME_BD_SFS_SIZE; ++i)
      std::swap(bd_sfs_[i], other.bd_sfs_[i]);
    for (int i = 0; i < SPLINE_VOLUME_PERIOD_INFO_SIZE; ++i)
      std::swap(periodicity_info_[i], other.periodicity_info_[i]);
    for (int i = 0; i < SPLINE_VOLUME_DEGEN_SIZE; ++i)
      std::swap(degen_[i], other.degen_[i]);
}



//===========================================================================
int SplineVolume::volumePeriodicity(int pardir, double epsilon) const
//===========================================================================
{
    int cont = GeometryTools::analyzePeriodicity(basis(pardir));
    int nmb = (cont < 0) ? 0 : cont;  // Number of rows of coefficients to check for equality
    int count;
    int numc[2];

    // Set indexes related to parameter direction
    int kj, kr;
    int start_idx, end_idx, del1, del2, del3;
    if (pardir == 0)
    {
	del1 = numCoefs(0);
	del2 = numCoefs(0)*numCoefs(1);
	del3 = 1;
	start_idx = 0;
	end_idx = numCoefs(0) - nmb - 1;
	numc[0] = numCoefs(1);
	numc[1] = numCoefs(2);
    }
    else if (pardir == 1)
    {
	del1 = 1;
	del2 = numCoefs(0)*numCoefs(1);
	del3 = numCoefs(0);
	start_idx = 0;
	end_idx = numCoefs(0)*(numCoefs(1)-nmb-1);
	numc[0] = numCoefs(0);
	numc[1] = numCoefs(2);
    }
    else if (pardir == 2)
    {
	del1 = 1;
	del2 = numCoefs(0);
	del3 = numCoefs(0)*numCoefs(1);;
	start_idx = 0;
	end_idx = numCoefs(0)*numCoefs(1)*(numCoefs(2)-nmb-1);
	numc[0] = numCoefs(0);
	numc[1] = numCoefs(1);
    }
    else 
	return -1; 

    vector<double>::const_iterator vol_coefs1 = coefs_begin();  // Coefficients of volume
    vector<double>::const_iterator vol_coefs2 = vol_coefs1 + end_idx*dim_;
    vector<double>::const_iterator vc1;  // Running coefficient pointer
    vector<double>::const_iterator vc2;  // Running coefficient pointer
    double dist2;
    double eps2 = epsilon*epsilon;
    for (count=0; count<=nmb; ++count, vol_coefs1+=del3*dim_, vol_coefs2+=del3*dim_)
    {
	vc1 = vol_coefs1;
	vc2 = vol_coefs2;
	for (kr=0; kr<numc[1]; ++kr, vc1=vol_coefs1+kr*del2*dim_, vc2=vol_coefs2+kr*del2*dim_)
	{
	    for (kj=0; kj<numc[0]; ++kj, vc1+=del1*dim_, vc2+=del1*dim_)
	    {
		dist2 = Utils::distance_squared(&vc1[0], &vc1[dim_], &vc2[0]);
		if (dist2 > eps2)
		    break;
	    }
	    if (kj < numc[0])
		break;
	}
	if (kr < numc[1])
	    break;
	
    }

    // The volume is open in this parameter direction
    int result;
    if (count == 0)
	result = -1;
    else if (cont < 0 && count == 1)
	result = 0;
    else
	result = count;

    periodicity_info_[pardir] = std::make_pair(result, epsilon);

    return periodicity_info_[pardir].first;
}

//===========================================================================
vector<shared_ptr<ParamSurface> > SplineVolume::getAllBoundarySurfaces() const
//===========================================================================
{
  vector<shared_ptr<SplineSurface> > sfs = getBoundarySurfaces(false);
  vector<shared_ptr<ParamSurface> > sfs2(sfs.begin(), sfs.end());
  return sfs2;
}

 //===========================================================================
vector<shared_ptr<SplineSurface> > SplineVolume::getBoundarySurfaces(bool do_clear) const
//===========================================================================
{
    
    vector<shared_ptr<SplineSurface> > boundary_sfs(6);
    if (bd_sfs_[0].get() != NULL && !do_clear)
    {
	// Surfaces already computed
	std::copy(bd_sfs_, bd_sfs_+6, boundary_sfs.begin());
	return boundary_sfs;
    }

    // For each parameter direction, fetch the boundary surfaces
    int dim = dim_;  // Array dimension of coefficient
    if (rational_)
	dim++;   

    int delta1[3], delta2[3], idx_start[3], idx_end[3];
    delta1[0] = numCoefs(0);
    delta2[0] = numCoefs(0)*numCoefs(1);
    idx_start[0] = 0;
    idx_end[0] = numCoefs(0) - 1;
    delta1[1] = 1;
    delta2[1] = numCoefs(0)*numCoefs(1);
    idx_start[1] = 0;
    idx_end[1] = numCoefs(0)*(numCoefs(1)-1);
    delta1[2] = 1;
    delta2[2] = numCoefs(0);
    idx_start[2] = 0;
    idx_end[2] = numCoefs(0)*numCoefs(1)*(numCoefs(2)-1);
    for (int ki=0; ki<3; ++ki)
    {
	int cont = GeometryTools::analyzePeriodicity(basis(ki));
	if (cont >= 0)
	{
	    // Periodic boundary conditions, use constParamSurface
	    bd_sfs_[2*ki] = shared_ptr<SplineSurface>(constParamSurface(startparam(ki), ki));
	    bd_sfs_[2*ki+1] = shared_ptr<SplineSurface>(constParamSurface(endparam(ki), ki));
	}
	else
	{
	    // Pick coefficients of boundary surfaces
	    int kj, kr, kh, del1, del2, nmbc=1;
	    int nmb[2], bdir[2];

	    // Compute number of coefficients
	    for (kj=0, kh=0; kj<3; ++kj)
		if (kj != ki)
		{
		    bdir[kh] = kj;
		    nmb[kh++] = numCoefs(kj);
		    nmbc *= numCoefs(kj);
		}

	    vector<double> coefs(nmbc*dim);  // Scratch for coefficient array
	    vector<double>::const_iterator vol_coefs = ctrl_begin();  // Coefficients of volume
	    vector<double>::const_iterator vc = vol_coefs;  // Running coefficient pointer

	    // Front surface
	    del1 = delta1[ki];
	    del2 = delta2[ki];
	    for (kr=0, kh=0; kr<nmb[1]; ++kr)
	      {
		vc = vol_coefs+kr*del2*dim;
		for (kj=0; kj<nmb[0]; ++kj, kh+=dim)
		{
		  if (kj > 0)
		    vc+=del1*dim;
		  std::copy(vc, vc+dim, coefs.begin()+kh);
		}
	      }

	    bd_sfs_[2*ki] = shared_ptr<SplineSurface>(new SplineSurface(basis(bdir[0]),
									 basis(bdir[1]),
									 &coefs[0],
									 dim_,
									 rational_));
	    // Back surface
	    for (kr=0, kh=0, vol_coefs+=idx_end[ki]*dim, vc=vol_coefs; kr<nmb[1]; 
		 ++kr)
	      {
		vc = vol_coefs+kr*del2*dim;
		for (kj=0; kj<nmb[0]; ++kj, kh+=dim)
		{
		  if (kj > 0)
		    vc+=del1*dim;
		  std::copy(vc, vc+dim, coefs.begin()+kh);
		}
	      }

	    bd_sfs_[2*ki+1] = shared_ptr<SplineSurface>(new SplineSurface(basis(bdir[0]),
									   basis(bdir[1]),
									   &coefs[0],
									   dim_,
									   rational_));
	}
    }
    std::copy(bd_sfs_, bd_sfs_+6, boundary_sfs.begin());
    return boundary_sfs;
}

//===========================================================================
void SplineVolume::checkDegeneracy(double tol, int is_degenerate[]) const
//===========================================================================
{
    bool b, r, t, l;
    int type;
    std::ofstream of("out_deg.g2");
    for (int ki=0; ki<6; ki++)
    {
      // bool is_degen = isDegenerate(ki, type, b, r, t, l, tol);
        isDegenerate(ki, type, b, r, t, l, tol);
	is_degenerate[ki] = type;
	bd_sfs_[ki]->writeStandardHeader(of);
	bd_sfs_[ki]->write(of);
    }
}

//===========================================================================
  void SplineVolume::getElementBdPar(int elem_ix, double elem_par[]) const
//===========================================================================
{
  // Fetch number of patches in all parameter directions
  int nu = numberOfPatches(0);
  int nv = numberOfPatches(1);
  int nw = numberOfPatches(2);

  if (elem_ix < 0 || elem_ix >= nu*nv*nw)
    {
      elem_par[0] = elem_par[1] = elem_par[2] = elem_par[3] = elem_par[4] = elem_par[5] = 0.0;
      return;
    }

  // 3-variate index
  int iw = elem_ix/(nu*nv);
  int iv = (elem_ix - iw*nu*nv)/nu;
  int iu = elem_ix - iw*nu*nv - iv*nu;

  // Parameter values
  vector<double> knots_u;
  vector<double> knots_v;
  vector<double> knots_w;
  const BsplineBasis basis_u = basis(0);
  const BsplineBasis basis_v = basis(1);
  const BsplineBasis basis_w = basis(2);
  basis_u.knotsSimple(knots_u);
  basis_v.knotsSimple(knots_v);
  basis_w.knotsSimple(knots_w);

  double u1 = knots_u[iu];
  double u2 = knots_u[iu+1];
  double v1 = knots_v[iv];
  double v2 = knots_v[iv+1];
  double w1 = knots_w[iw];
  double w2 = knots_w[iw+1];

  elem_par[0] = u1;
  elem_par[1] = u2;
  elem_par[2] = v1;
  elem_par[3] = v2;
  elem_par[4] = w1;
  elem_par[5] = w2;
}

//===========================================================================
  vector<shared_ptr<SplineSurface> > SplineVolume::getElementBdSfs(int elem_ix,
								   double elem_par[]) const
//===========================================================================
{
  vector<shared_ptr<SplineSurface> > result;

  // Fetch number of patches in all parameter directions
  int nu = numberOfPatches(0);
  int nv = numberOfPatches(1);
  int nw = numberOfPatches(2);

  if (elem_ix < 0 || elem_ix >= nu*nv*nw)
    return result;

  // 3-variate index
  int iw = elem_ix/(nu*nv);
  int iv = (elem_ix - iw*nu*nv)/nu;
  int iu = elem_ix - iw*nu*nv - iv*nu;

  // Parameter values
  vector<double> knots_u;
  vector<double> knots_v;
  vector<double> knots_w;
  const BsplineBasis basis_u = basis(0);
  const BsplineBasis basis_v = basis(1);
  const BsplineBasis basis_w = basis(2);
  basis_u.knotsSimple(knots_u);
  basis_v.knotsSimple(knots_v);
  basis_w.knotsSimple(knots_w);

  double u1 = knots_u[iu];
  double u2 = knots_u[iu+1];
  double v1 = knots_v[iv];
  double v2 = knots_v[iv+1];
  double w1 = knots_w[iw];
  double w2 = knots_w[iw+1];

  elem_par[0] = u1;
  elem_par[1] = u2;
  elem_par[2] = v1;
  elem_par[3] = v2;
  elem_par[4] = w1;
  elem_par[5] = w2;

  // Fetch full side surfaces
  shared_ptr<SplineSurface> tmp_u1(constParamSurface(u1, 0));
  shared_ptr<SplineSurface> tmp_u2(constParamSurface(u2, 0));
  shared_ptr<SplineSurface> tmp_v1(constParamSurface(v1, 1));
  shared_ptr<SplineSurface> tmp_v2(constParamSurface(v2, 1));
  shared_ptr<SplineSurface> tmp_w1(constParamSurface(w1, 2));
  shared_ptr<SplineSurface> tmp_w2(constParamSurface(w2, 2));

  // Restrict surfaces
  result.resize(6);
  result[0] = shared_ptr<SplineSurface>(tmp_u1->subSurface(v1, w1, v2, w2));
  result[1] = shared_ptr<SplineSurface>(tmp_u2->subSurface(v1, w1, v2, w2));
  result[2] = shared_ptr<SplineSurface>(tmp_v1->subSurface(u1, w1, u2, w2));
  result[3] = shared_ptr<SplineSurface>(tmp_v2->subSurface(u1, w1, u2, w2));
  result[4] = shared_ptr<SplineSurface>(tmp_w1->subSurface(u1, v1, u2, v2));
  result[5] = shared_ptr<SplineSurface>(tmp_w2->subSurface(u1, v1, u2, v2));

  return result;
}

//===========================================================================
bool SplineVolume::isDegenerate(int which_sf, int& type, bool& b, bool& r,
				bool& t, bool& l, double tol) const
//===========================================================================
{
    if (degen_[which_sf].is_set_ &&
	fabs(degen_[which_sf].tol_ - tol) < 1.0e-15)
    {
	type = degen_[which_sf].type_;
	b = degen_[which_sf].b_;
	r = degen_[which_sf].r_;
	t = degen_[which_sf].t_;
	l = degen_[which_sf].l_;
	return (type != 0);
    }

    // Make sure that all boundary surfaces are computed
    vector<shared_ptr<SplineSurface> > bd_sf = getBoundarySurfaces();

    bool deg = bd_sf[which_sf]->isDegenerate(b, r, t, l, tol);
    double dist2;
    double tol2 = tol*tol;
    if (deg)
    {
	type = 0;
	if (b && r && t && l)
	{
	    // Check if the surface degenerates to a point
	    vector<double>::const_iterator coef = bd_sf[which_sf]->coefs_begin();
	    int nn = bd_sf[which_sf]->numCoefs_u()*bd_sf[which_sf]->numCoefs_v();
	    int ki;
	    for (ki=1; ki<nn; ++ki)
	    {
		dist2 = Utils::distance_squared(&coef[0], &coef[dim_], &coef[ki*dim_]);
		if (dist2 > tol2)
		    break;
	    }
	    if (ki == nn)
		type = 3;
		
	}
	if (type==0 && b && t)
	{
	    // Check if the surface degenerates to a line
	    int nn1 = bd_sf[which_sf]->numCoefs_u();
	    int nn2 = bd_sf[which_sf]->numCoefs_v();
	    vector<double>::const_iterator coef = bd_sf[which_sf]->coefs_begin();
	    int ki;
	    for (ki=0; ki<nn2; ki++, coef+=nn1*dim_)
	    {
		dist2 = Utils::distance_squared(&coef[0], &coef[dim_], &coef[(nn1-1)*dim_]);
		if (dist2 > tol2)
		    break;
	    }
	    if (ki == nn2)
		type = 2;
	}
	if (type==0 && r && l)
	{
	    // Check if the surface degenerates to a line
	    int nn1 = bd_sf[which_sf]->numCoefs_u();
	    int nn2 = bd_sf[which_sf]->numCoefs_v();
	    vector<double>::const_iterator coef = bd_sf[which_sf]->coefs_begin();
	    vector<double>::const_iterator coef2 = coef + (nn2-1)*nn1*dim_;
	    int ki;
	    for (ki=0; ki<nn1; ki++, coef+=dim_, coef2+=dim_)
	    {
		dist2 = Utils::distance_squared(&coef[0], &coef[dim_], &coef2[0]);
		if (dist2 > tol2)
		    break;
	    }
	    if (ki == nn1)
		type = 2;
	}
	if (type == 0)
	    type = 1;
    }
    // VSK, May 2010. Maybe a good idea, but doesn't work as planned
//     else if (dim_ == 3)
//     {
// 	// Check if the surface degenerates to a line anyway
// 	vector<double>::iterator coef = bd_sf[which_sf]->coefs_begin();
// 	Point pnt(&coef[0], &coef[dim_], true);
// 	Point pnt2;
// 	Point diff; 
// 	int nn = bd_sf[which_sf]->numCoefs_u()*bd_sf[which_sf]->numCoefs_v();
// 	int ki;
// 	for (ki=1; ki<nn; ++ki)
// 	{
// 	    pnt2 = Point(&coef[ki*dim_], &coef[(ki+1)*dim_], true);
// 	    diff = pnt2 - pnt;
// 	    double len = diff.normalize_checked();
// 	    if (len != 0.0)
// 		break;
// 	}
	
// 	// A line is found, check if all remaining coefficients lie on this line
// 	Point norm1, norm2;
// 	getPlaneNormals(pnt, diff, norm1, norm2);
// 	for (; ki<nn; ki++)
// 	{
// 	    double d1 = (pnt - Point(&coef[ki*dim_], &coef[(ki+1)*dim_], true))*norm1;
// 	    double d2 = (pnt - Point(&coef[ki*dim_], &coef[(ki+1)*dim_], true))*norm2;
// 	    if (d1 > tol || d2 > tol)
// 		break;
// 	}
// 	if (ki == nn)
// 	    type = 2;
// 	else
// 	    type = 0;
//     }
    else
	type = 0;

    return deg;
}


//===========================================================================
void SplineVolume::translate(const Point& vec)
//===========================================================================
{
  ALWAYS_ERROR_IF(dim_ != vec.dimension(), "Volume and translation vector of different dimension");

  // Change rcoefs_ if rational
  if (rational_)
    for (vector<double>::iterator it = rcoefs_.begin(); it != rcoefs_.end(); )
      {
	double w = it[dim_];
	for (int i = 0; i < dim_; ++i, ++it)
	  (*it) += vec[i] * w;
	++it;
      }

  // Change coefs, both when rational and not rational
  for (vector<double>::iterator it = coefs_.begin(); it != coefs_.end(); )
    {
      for (int i = 0; i < dim_; ++i, ++it)
	(*it) += vec[i];
    }
}


//===========================================================================
void SplineVolume::scale(double fac)
//===========================================================================
{
  // Scale coefficients
  for (size_t ki=0; ki<coefs_.size(); ++ki)
    coefs_[ki] *= fac;

  if (rational_)
    {
      for (size_t ki=0; ki<rcoefs_.size(); ki+=(dim_+1))
	{
	  for (int kj=0; kj<dim_; ++kj)
	    rcoefs_[ki+kj] *= fac;
	}
    }
}

// Added by KMO for ICADA usage.
//===========================================================================
void SplineVolume::deform(const std::vector<double>& vec, int vdim)
//===========================================================================
{
  int i, j;
  vector<double>::iterator it;
  if (vdim == 0) vdim = dim_;

  // Change rcoefs_ if rational
  if (rational_)
    for (it = rcoefs_.begin(), j = 0; it != rcoefs_.end(); ++it)
      {
	double w = it[dim_];
	for (i = 0; i < dim_ && i < vdim; ++i)
	  it[i] += vec[j+i] * w;
	it += dim_;
	j += vdim;
      }

  // Change coefs, both when rational and not rational
  for (it = coefs_.begin(), j = 0; it != coefs_.end(); )
    {
      for (i = 0; i < dim_; ++i)
	it[i] += vec[j+i];
      it += dim_;
      j += vdim;
    }
}

//===========================================================================
void SplineVolume::add(const SplineVolume* other, double tol)
//===========================================================================
{
  int ord_u = basis_u_.order();
  int ord_v = basis_v_.order();
  int ord_w = basis_w_.order();
  int ncoefs_u = basis_u_.numCoefs();
  int ncoefs_v = basis_v_.numCoefs();
  int ncoefs_w = basis_w_.numCoefs();
  int ncoefs = ncoefs_u * ncoefs_v * ncoefs_w;
  int kdim = rational_ ? dim_+1 : dim_;
  ALWAYS_ERROR_IF(ord_u != other->basis_u_.order(), "Volumes have different order in first parameter direction");
  ALWAYS_ERROR_IF(ord_v != other->basis_v_.order(), "Volumes have different order in second parameter direction");
  ALWAYS_ERROR_IF(ord_w != other->basis_w_.order(), "Volumes have different order in third parameter direction");
  ALWAYS_ERROR_IF(ncoefs_u != other->basis_u_.numCoefs(), "Volumes have different number of coefficients in first parameter direction");
  ALWAYS_ERROR_IF(ncoefs_v != other->basis_v_.numCoefs(), "Volumes have different number of coefficients in second parameter direction");
  ALWAYS_ERROR_IF(ncoefs_w != other->basis_w_.numCoefs(), "Volumes have different number of coefficients in third parameter direction");
  ALWAYS_ERROR_IF(dim_ != other->dim_, "Volumes have different geometry space dimension");
  ALWAYS_ERROR_IF(rational_ != other->rational_, "Can not add rational to non-rational volume");

  vector<double>::const_iterator knots, knots_other;
  knots = basis_u_.begin();
  knots_other = other->basis_u_.begin();
  for (int i = 0; i < ncoefs_u + ord_u; ++i)
    ALWAYS_ERROR_IF(knots[i] != knots_other[i], "Volumes have different knot vector in first parameter direction");
  knots = basis_v_.begin();
  knots_other = other->basis_v_.begin();
  for (int i = 0; i < ncoefs_v + ord_v; ++i, ++knots, ++knots_other)
    ALWAYS_ERROR_IF((*knots) != (*knots_other), "Volumes have different knot vector in second parameter direction");
  knots = basis_w_.begin();
  knots_other = other->basis_w_.begin();
  for (int i = 0; i < ncoefs_w + ord_w; ++i, ++knots, ++knots_other)
    ALWAYS_ERROR_IF((*knots) != (*knots_other), "Volumes have different knot vector in third parameter direction");

  if (rational_)
    {
      vector<double>::const_iterator weights = rcoefs_.begin() + dim_;
      vector<double>::const_iterator weights_other = other->rcoefs_.begin() + dim_;
      for (int i = 0; i < ncoefs; ++i, weights += kdim, weights_other += kdim)
	ALWAYS_ERROR_IF(fabs ((*weights) - (*weights_other)) >= tol, "Volumes have different weights");
    }

  // All tests passed - we can now add coefficients

  vector<double>::iterator coefs_it = coefs_.begin();
  vector<double>::const_iterator other_coefs_it = other->coefs_.begin();
  for (int i = 0; i < ncoefs * dim_; ++i, ++coefs_it, ++other_coefs_it)
    (*coefs_it) += (*other_coefs_it);

  // Update rcoefs_ in rational case

  if (rational_)
    {
      coefs_it = coefs_.begin();
      vector<double>::iterator rcoefs_it = rcoefs_.begin();
      for (int i = 0; i < ncoefs; ++i, coefs_it += dim_, rcoefs_it += kdim)
	for (int k = 0; k < dim_; ++k)
	  rcoefs_it[k] = rcoefs_it[dim_] * coefs_it[k];
    }

}


//===========================================================================
bool SplineVolume::isLeftHanded()
//===========================================================================
{
  if (dim_ != 3)
    return false;

  vector<Point> pts(4);
  point(pts,
	0.5 * (startparam(0) + endparam(0)),
	0.5 * (startparam(1) + endparam(1)),
	0.5 * (startparam(2) + endparam(2)),
	1);

  return ((pts[1] % pts[2]) * pts[3]) < 0.0;
}


//===========================================================================
//
//                 Private helper functions
//
//===========================================================================



//===========================================================================
void SplineVolume::updateCoefsFromRcoefs()
//===========================================================================
{
    int nmb = numCoefs(0)*numCoefs(1)*numCoefs(2);
    coefs_.resize(nmb*dim_);
    SplineUtils::make_coef_array_from_rational_coefs(&rcoefs_[0],
					&coefs_[0],
					nmb,
					dim_);
}

//===========================================================================
void SplineVolume::getPlaneNormals(Point pnt, Point vec, Point& norm1, Point& norm2) const
//===========================================================================
{
    Point dum(0.0, 0.0, 0.0);
    if (fabs(vec[0]) < fabs(vec[1]) && fabs(vec[0]) < fabs(vec[2]))
	dum[0] = 1.0;
    else if (fabs(vec[1]) < fabs(vec[2]))
	dum[1] = 1.0;
    else 
	dum[2] = 1.0;

    norm1 = vec.cross(dum);
    norm1.normalize();

    norm2 = vec.cross(norm1);
    norm2.normalize();
}


} // namespace Go
