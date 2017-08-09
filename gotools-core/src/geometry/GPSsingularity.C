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

#include "GoTools/utils/GeneralFunctionMinimizer.h"
#include "GoTools/geometry/ParamSurface.h"
#include <fstream>

using namespace Go;
using namespace std;

namespace {
//===========================================================================
// squared distance function between a point and a surface. Used by the
// minimization algorithm initiated by SplineSurface::closestPoint
class SingDist2 {
//===========================================================================
public:
    SingDist2(const ParamSurface *sf,
	      const RectDomain* rd) 
	: sf_(sf), tmp_ptvec_(6) {
	if (rd) {
	    ll_[0] = rd->umin(); 
	    ll_[1] = rd->vmin();	    
	    ur_[0] = rd->umax(); 
	    ur_[1] = rd->vmax();
	} else {
	    RectDomain domain = sf->containingDomain();
	    ll_[0] = domain.umin();
	    ll_[1] = domain.vmin();
	    ur_[0] = domain.umax();
	    ur_[1] = domain.vmax();
	}
    }
    inline double operator()(const double* arg) const;
    inline void grad(const double* arg, double* res) const;
    inline double minPar(int pardir) const { 
	return ll_[pardir];
	//return (pardir == 0) ? sf_.startparam_u() : sf_.startparam_v();
    }
    inline double maxPar(int pardir) const {
	return ur_[pardir];
	//return (pardir == 0) ? sf_.endparam_u() : sf_.endparam_v();
    }
private:
    const ParamSurface *sf_;
    double ll_[2]; // lower left corner of domain
    double ur_[2]; // upper right corner of domain
    mutable Point tmp_pt_, tmp_pt_du_, tmp_pt_dv_;
    mutable vector<Point> tmp_ptvec_;
};

//===========================================================================
double SingDist2::operator()(const double* arg) const
//===========================================================================
{
    sf_->point(tmp_ptvec_, arg[0], arg[1], 1);
    tmp_pt_ = tmp_ptvec_[1]%tmp_ptvec_[2];
    return tmp_pt_.length2();
}

//===========================================================================
void SingDist2::grad(const double* arg, double* res) const
//===========================================================================
{
    sf_->point(tmp_ptvec_, arg[0], arg[1], 2);
    tmp_pt_ = tmp_ptvec_[1]%tmp_ptvec_[2];
    tmp_pt_du_ = tmp_ptvec_[3]%tmp_ptvec_[2] + tmp_ptvec_[1]%tmp_ptvec_[4];
    tmp_pt_dv_ = tmp_ptvec_[4]%tmp_ptvec_[2] + tmp_ptvec_[1]%tmp_ptvec_[5];
    res[0] = 2 * tmp_pt_ * tmp_pt_du_;
    res[1] = 2 * tmp_pt_ * tmp_pt_dv_;
}

//===========================================================================
// find a good seed for closest point computation
void seedfind(const ParamSurface& sf, 
	      const RectDomain* rd,
	      double& u,
	      double& v)
//===========================================================================
{
    // Use midpoint
    if (rd)
    {
	u = 0.5*(rd->umin() + rd->umax());
	v = 0.5*(rd->vmin() + rd->vmax());
    }
    else
    {
	RectDomain domain = sf.containingDomain();
	u = 0.5*(domain.umin() + domain.umax());
	v = 0.5*(domain.vmin() + domain.vmax());
    }
}

}; // end anonymous namespace 


namespace Go {

//===========================================================================
void ParamSurface::singularity(double& sing_u,
			       double& sing_v, 
			       Point& sing_pt,
			       double& sing_dist,
			       double epsilon,
			       const RectDomain* rd,
			       double *seed) const
//===========================================================================
{
    SingDist2 dist_fun(this, rd);
    
    double seed_buf[2];
    if (!seed) {
	// no seed given, we must compute one
	seed = seed_buf;
	seedfind(*this, rd, seed[0], seed[1]);
    }
    
    // define distance function
    FunctionMinimizer<SingDist2> funmin(2, dist_fun, seed, epsilon);
    
    // minimise distance function
    minimise_conjugated_gradient(funmin); //, 3) ; // number of iterations in each cycle

    // return results
    sing_u = funmin.getPar(0);
    sing_v = funmin.getPar(1);
    sing_dist = sqrt(funmin.fval());
    point(sing_pt, sing_u, sing_v);
}

} // namespace Go

