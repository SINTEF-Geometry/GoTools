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


namespace Go {
    class GeomObject;
}
using Go::GeomObject;



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

  std::shared_ptr<GeomObject> operator[] (int idx) const
    { return geomobj_[idx]; }

  int size() const
    { return (int)geomobj_.size(); }

  ftTangPriority getType() const
    { return type_; }

  void setType(ftTangPriority type)
    { type_ = type; }

  void addGeomObj(std::shared_ptr<GeomObject> obj)
    { geomobj_.push_back(obj); }
      

protected:
  std::vector<std::shared_ptr<GeomObject> > geomobj_;
  ftTangPriority type_;
  
};

#endif // _FTGROUPGEOM_H
