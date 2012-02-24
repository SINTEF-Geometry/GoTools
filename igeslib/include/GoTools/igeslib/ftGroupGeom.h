//===========================================================================
//                                                                           
// File: ftGroupGeom.h                                                   
//                                                                           
// Created: Wed Mar 21 2001                                         
//                                                                           
// Author: Vibeke Skytt 
//                                                                           
// Revision: 
//                                                                           
// Description: 
//
// Implementation files: ftGroupGeom.C
//                                                                           
//===========================================================================

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
