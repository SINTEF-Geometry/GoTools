//===========================================================================
//
// File : VolumeAdjacency
//
// Created: January 5, 2010
//
// Author: Vibeke Skytt, SINTEF
//
// Revision: 
//
// Description:
//
//===========================================================================

#ifndef _VOLUMEADJACENCY_H
#define _VOLUMEADJACENCY_H

#include "GoTools/compositemodel/Body.h"

namespace Go
{
  class ftSurface;
  class ParamSurface;

  /// \brief Adjacency analysis of volume models

  class VolumeAdjacency
  {
  public:
    /// Constructor giving tolerances for adjacency analysis
    VolumeAdjacency(double gap, double neighbour);

    /// Destructor
    ~VolumeAdjacency();

    /// Add the solids of the volume model
    void setAdjacency(std::vector<shared_ptr<Body> >& solids);

    /// Add the solids of the volume model, indicates the solids for which 
    /// topology analysis is not performed already
    void setAdjacency(std::vector<shared_ptr<Body> >& solids, 
		      int new_solid_pos);

    private:
    /// Gap between volumes
    double gap_;

    /// Tolerance to check if two volumes can be adjacent at all, 
    /// neighbour_ > gap_
    double neighbour_;

    /// Check for adjacency between two solids
    void setAdjacency(shared_ptr<Body> solid1, shared_ptr<Body> solid2);

    /// Check for adjacency between two boundary faces, split faces in case
    /// of partial adjacency (one face embedded in the other, general partial
    /// coincidence is not handled)
    int faceAdjacency(shared_ptr<ftSurface> face1, 
		      shared_ptr<ftSurface> face2,
		      std::vector<shared_ptr<ftSurface> >& new_faces1,
		      std::vector<shared_ptr<ftSurface> >& new_faces2);

    /// Handle embedded boundary surfaces
    void splitSurface(shared_ptr<ParamSurface> srf1,
		      shared_ptr<ParamSurface> srf2,
		      std::vector<shared_ptr<ParamSurface> >& new_sfs);
  };

} // namespace Go


#endif // _VOLUMEADJACENCY_H
