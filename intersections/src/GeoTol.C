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

#include "GoTools/intersections/GeoTol.h"
#include <iostream>

using namespace std;
namespace Go
{
//===========================================================================
GeoTol::GeoTol(double epsge, double rel_par_res, double numerical_tol)
    : epsge_(epsge), 
      eps_bracket_(0.1), 
      ang_tol_(0.01), 
      rel_par_res_(rel_par_res),
      numerical_tol_(numerical_tol)
//===========================================================================
  {
      if (epsge <= 10e-5)
	  ref_ang_ = 0.01;
      else if (epsge <= 10e-3)
	  ref_ang_ = 0.1;
      else
	  ref_ang_ = 0.5;
  }

//===========================================================================
GeoTol::GeoTol(GeoTol *epsge)
{
    epsge_ = epsge->epsge_;
    eps_bracket_ = epsge->eps_bracket_;
    ref_ang_ = epsge->ref_ang_;
    ang_tol_ = epsge->ang_tol_;
    rel_par_res_ = epsge->rel_par_res_;
    numerical_tol_ = epsge->numerical_tol_;
}

//===========================================================================
  GeoTol::~GeoTol()
//===========================================================================
{
}

//===========================================================================
void GeoTol::write(ostream& os) const
//===========================================================================
{
    os << epsge_ << ' ';
    os << eps_bracket_ << ' ';
    os << ref_ang_ << ' ';
    os << ang_tol_ << ' ';

    os << rel_par_res_ << ' ';
    os << numerical_tol_ << ' ';
}

//===========================================================================
void GeoTol::read(istream& is)
//===========================================================================
{
    is >> epsge_;
    is >> eps_bracket_;
    is >> ref_ang_;
    is >> ang_tol_;
    is >> rel_par_res_;
    is >> numerical_tol_;
}

};
