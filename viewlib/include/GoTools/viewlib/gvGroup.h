//===========================================================================
//                                                                           
// File: gvGroup.h                                                          
//                                                                           
// Created: Wed Oct  9 09:28:14 2002                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@math.sintef.no>
//                                                                           
// Revision: $Id: gvGroup.h,v 1.2 2007-05-02 14:39:24 sbr Exp $
//                                                                           
// Description: Structure used to store name and index of members of a group
//              of objects.                                                                           
//
//===========================================================================

#ifndef _GVGROUP_H
#define _GVGROUP_H

#include <vector>
#include <QLineEdit>

/// Structure used to store name and index of members of a group of objects.   
class gvGroup
{
 public: 

  gvGroup(std::vector<int>& members, QString name);

  ~gvGroup();

  QString name() const;

  int size() const;

  int operator[] (int index) const;


 private:

  std::vector<int> members_;
  QString name_;

};


#endif // _GVGROUP_H
