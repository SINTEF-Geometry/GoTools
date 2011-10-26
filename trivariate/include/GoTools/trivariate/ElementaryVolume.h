//===========================================================================
//
// File : ElementaryVolume.h
//
// Created: Fri Nov  6 13:38:38 2009
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================



#ifndef __ELEMENTARYVOLUME_H
#define __ELEMENTARYVOLUME_H


#include "GoTools/trivariate/ParamVolume.h"


namespace Go
{

  class SplineVolume;


  /// \brief ElementaryVolume is a base class for elementary volumes
  /// like boxes and solid cylinders. Such volumes have natural
  /// parametrizations and ElementartVolume is therefore a subclass of
  /// ParamVolume. These volumes are non-self-intersecting.

  class ElementaryVolume : public ParamVolume
  {

  public:

    /// Virtual destructor, enables safe inheritance.
    virtual ~ElementaryVolume() { }

    // --- Functions inherited from GeomObject ---
    virtual ElementaryVolume* clone() const = 0;


    virtual SplineVolume* geometryVolume() const = 0;

  };    // Class ElementaryVolume


} // namespace Go



#endif    // #ifndef __ELEMENTARYVOLUME_H
