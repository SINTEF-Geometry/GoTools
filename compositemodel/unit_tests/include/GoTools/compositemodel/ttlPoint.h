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

#ifndef _TTLPOINT_H_
#define _TTLPOINT_H_

#include "GoTools/compositemodel/ftPointSet.h"

namespace Go
{

  /** Utility class for representing a ftSamplePoint with additional parameter value.
   * Useful when dealing with sets of surfaces treated as one.
   */
  class ttlPoint
  {

  public:
    /// Constructor
  ttlPoint(PointIter pnt_iter, double x, double y, double z = 0.0)
    : par_val_(Vector3D(x, y, z)), pnt_iter_(pnt_iter) {}
    /// Destructor
    ~ttlPoint() {}
  
    inline const PointIter& pnt_iter() const { return pnt_iter_; };
    inline double x() const { return par_val_[0]; } ;
    inline double y() const { return par_val_[1]; } ;
    inline double z() const { return par_val_[2]; } ;


  private:

    // As point is to be included in a larger set, we must set new parameter value.
    Vector3D par_val_; // z coordinate is by default set to 0.
    PointIter pnt_iter_; // We need access to the actual object.

  };


} // namespace Go

#endif
