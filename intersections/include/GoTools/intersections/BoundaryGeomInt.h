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

#ifndef _BOUNDARYGEOMINT_H
#define _BOUNDARYGEOMINT_H


#include <memory>


namespace Go {


class ParamGeomInt;


/// This struct is a helper struct that bundles boundary information
/// for an object of type ParamGeomInt.

struct BoundaryGeomInt {

    shared_ptr<ParamGeomInt> bd_obj_;
    int pardir_; // 2dim case: 0 || 1 (dir of par_, i.e. opp that of
		 // bd_obj_).
    double par_;

    /// Constructor
    /// \param bd_obj shared pointer to the object of interest
    /// \param dir parameter direction that specifies the boundary
    /// \param par the value of the relevant parameter at the boundary
    BoundaryGeomInt(shared_ptr<ParamGeomInt> bd_obj,
		    int dir, double par)
    {
	bd_obj_ =  bd_obj;
	pardir_ = dir;
	par_ = par;
    }

    /// Destructor
    ~BoundaryGeomInt() { }

    /// Get the parameter direction that specifies the boundary
    /// \return the parameter direction
    int getDir()
    { return pardir_; }
     
    /// Get the value of the relevant parameter at the boundary
    /// \return the parameter value
    double getPar()
    { return par_; }

    /// Get the object which the boundary belongs to
    /// \return shared pointer to the object
    shared_ptr<ParamGeomInt> getObject()
    { return bd_obj_; }
};


} // namespace Go


#endif // _BOUNDARYGEOMINT_H

