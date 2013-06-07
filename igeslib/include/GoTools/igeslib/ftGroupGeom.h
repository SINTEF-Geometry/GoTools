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

#ifndef _FTGROUPGEOM_H
#define _FTGROUPGEOM_H


#include <vector>
#include "GoTools/igeslib/ftTangPriority.h"
#include "GoTools/utils/config.h"


namespace Go {


  /// A group of geometrical objects
class ftGroupGeom
{
public:
  // Constructor
       ftGroupGeom()
       : type_(ftNoType)
    {}

  // Destructor
       ~ftGroupGeom()
    {}

       /// Fetch object number idx
  shared_ptr<GeomObject> operator[] (int idx) const
    { return geomobj_[idx]; }

  /// The number of objects in the group
  int size() const
    { return (int)geomobj_.size(); }

  /// Whether the group of objects (surfaces) are a master, a
  /// slave or not specified. Related to tangent plane continuity
  /// between groups of surfaces
  ftTangPriority getType() const
    { return type_; }

  /// Set whether the group of objects (surfaces) are a master, a
  /// slave or not specified. Related to tangent plane continuity
  /// between groups of surfaces
  void setType(ftTangPriority type)
    { type_ = type; }

  /// Add a new geometry entitiy to the group
  void addGeomObj(shared_ptr<GeomObject> obj)
    { geomobj_.push_back(obj); }
      

protected:
  std::vector<shared_ptr<GeomObject> > geomobj_;
  ftTangPriority type_;
  
};

} // namespace Go

#endif // _FTGROUPGEOM_H
