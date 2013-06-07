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

#ifndef __VOLUMEMODELCREATOR_H
#define __VOLUMEMODELCREATOR_H

#include "GoTools/utils/config.h"

namespace Go
{
  class SurfaceModel;
  class VolumeModel;

  namespace VolumeModelCreator
  {
    // Given a boundary represented solid, check if it is possible to create
    // a block structured volume model using rotation. In that case, create 
    // the model
    bool createRotationalModel(shared_ptr<SurfaceModel>& sfmodel,
			       shared_ptr<VolumeModel>& volmodel);

    // This function recognizes a planar surface set swept along a line
    // in a boundary represented model and creates the associated
    // volume model. 
    // Non-planar surface sets in a similar configuration is not recognized and
    // neither is linear sweep where the curve along which to sweep is not linear
    bool linearSweptModel(shared_ptr<SurfaceModel>& sfmodel,
			  shared_ptr<VolumeModel>& volmodel);

  };  // VolumeModelCreator
}  // Go
#endif // __VOLUMEMODELCREATOR_H
