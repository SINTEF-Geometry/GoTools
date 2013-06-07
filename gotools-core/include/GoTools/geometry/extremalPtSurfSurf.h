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

#ifndef _EXTREMALPTSURFSURF_H
#define _EXTREMALPTSURFSURF_H

/** @file extremalPtSurfSurf.h
 * Declaration file for an intersection function operating on 
 * objects of class ParamSurface belonging to geometry.
 * Ported from the sisl function shsing_ext.
 */


#include "GoTools/geometry/ParamSurface.h"

namespace Go {
  /** Finds one point in each surface where the two normals
   * are parallel to each other and to the difference vector
   * between the two points.
   * The edge info makes it possible to constrain the search
   * to an ISO-curve in one of or both surfaces.
   * Ported from the sisl-function \c shsing_ext.
   @verbatim
METHOD :  - Start with a guess value (u,v) in domain of surface 1 (S(u,v))
      (a) - Find domain value (r,t) of closest point (to S(u,v) in surface 2 (Q(r,t))
          - If vf1(u,v) = <Su,Normal(Q> and vf2(u,v)= <Sv,Normal(Q> is small enough  stop
                 (<,> means scalar prod.)
          - Find du and dv by taylorizing vf1 and vf2.
            This include finding the derivatives of the closest point function (r(u,v),t(u,v))
            with respect to u and v. (called h(u,v) in article, see comments in shsing_s9dir)
          - u:= u+du v:= v+dv, goto (a)

   
REFERENCES : Solutions of tangential surface and curve intersections.
             R P Markot and R L Magedson
             Computer-Aided Design; vol. 21, no 7 sept. 1989, page 421-429
   @endverbatim
   * \param psurf1 Pointer to the first surface.
   * \param psurf2 Pointer to the second surface
   * \param constraints Constraints flag. 
   * - constraints
   *     -# constraints[0] Constraint 1. surface
   *     -# constraints[1] Constraint 2. surface
   * - constraints values
   *     -# -1 : No constraints.
   *     -#  0 : Constant in 1. par direction
   *     -#  1 : Constant in 2. par direction
   *
   * ------------------ or ----------------
   *  \li \c constraints[0] Constraint 1. surface
   *  \li \c constraints[1] Constraint 2. surface
   *  \li \c -1 No constraints.
   *  \li \c  0 Constant in 1. par direction
   *  \li \c  1 Constant in 2. par direction
   *
   * \param constraints_par Constraints parameters.
   *   -# constraints_par[0] Constraints parameter 1. surface.
   *   -# constraints_par[1] Constraints parameter 2. surface.
   * \param limit  - Parameter borders of both surfaces.
   *  -# limit[0] - limit[1] Parameter interval. 1. surface  1. direction.
   *  -# limit[2] - limit[3] Parameter interval. 1. surface  2. direction.
   *  -# limit[4] - limit[5] Parameter interval. 2. surface  1. direction.
   *  -# limit[6] - limit[7] Parameter interval. 2. surface  2. direction.
   *
   * \param enext Parameter start value for iteration(4 values).
   * \param angle_tol The angular tolerance for success.
   * \param gpos Parameter values of the found singularity(4 values).
   * \return Status message. =1 Extremum found.  = 0 Extremum NOT found.
   */
int extremalPtSurfSurf(ParamSurface* psurf1, ParamSurface* psurf2,
		      int    constraints[2], double constraints_par[2],
		      double limit[], double enext[],
		      double gpos[], double angle_tol);


// Anonymous namespace
namespace {


  /** To be sure that we are inside the border of the parameter plane.
   * If we are outside clipping is used to adjust the step value.
   * Ported from the sisl function \c shsing_s9corr.
   *
   * \param gd Proposed delta values.
   * \param coef Current position.
   * \param limit Parameter borders of both surfaces.
   */
void insideParamDomain(double gd[], double coef[],double limit[]);


  /** Calculates the increments in the first surface domain.
   * Ported from the sisl function \c shsing_s9dir.
   * \param evals Value and derivatives on first surface.
   * \param evalq Value and derivatives on second surface.
   * \param cdiff1 Parameter increments in two directions.
   @verbatim
   METHOD  : See comments in function shsing.
             The only thing missing in the article is the derivation of the
             function h(u,v). Calling this function (r(u,v), t(u,v))
             we know that
             For all (u,v) (<,> meaning scalar product)
             <S(u,v)-Q(r(u,v),t(u,v)),Qt(r(u,v),t(u,v))> = 0
             <S(u,v)-Q(r(u,v),t(u,v)),Qr(r(u,v),t(u,v))> = 0
             This means that (derivation by u and v)
             <Su-[Qt*tu+Qr*ru],Qt> + <S-Q,Qtt*tu + Qtr*ru> = 0
             <Su-[Qt*tu+Qr*ru],Qr> + <S-Q,Qrt*tu + Qrr*ru> = 0
             <Sv-[Qt*tv+Qr*rv],Qt> + <S-Q,Qtt*tv + Qtr*rv> = 0
             <Sv-[Qt*tv+Qr*rv],Qr> + <S-Q,Qrt*tv + Qrr*rv> = 0
             Solving these four equations gives us ru,rv,tu,tv.
   @endverbatim
   */
void nextStep(double cdiff[],std::vector<Point> evals,
	      std::vector<Point> evalq);


} // Anonymous namespace

} // namespace Go  

#endif // _EXTREMALPTSURFSURF_H
