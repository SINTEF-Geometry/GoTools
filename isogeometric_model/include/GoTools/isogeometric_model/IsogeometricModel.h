//===========================================================================
//
// File : IsogeometricModel.h
//
// Created: 30th of October, 2009
//
// Author: Vibeke Skytt, SINTEF
//
// Revision: 
//
// Description:
//
//===========================================================================


#ifndef __ISOGEOMETRICMODEL_H
#define __ISOGEOMETRICMODEL_H


#include "GoTools/topology/tpTopologyTable.h"



namespace Go
{

  class IsogeometricModel
  {
  public:
    // Constructor
    IsogeometricModel(tpTolerances toptol)
      : toptol_(toptol)
    {
    }


    // Destructor
    virtual ~IsogeometricModel();

    virtual int getNmbOfBoundaries() const = 0;

    // Ensure minimum degree of solution space
    virtual void setMinimumDegree(int degree, int solutionspace_idx) = 0;

    // Update spline spaces of the solution to ensure consistence
    virtual void updateSolutionSplineSpace() = 0;

    // Get tolerances
    virtual tpTolerances getTolerances() const
    { return toptol_; }

    // Set tolerances
    virtual void setTolerances(tpTolerances t)
    { toptol_ = t; }

  private:

    tpTolerances toptol_;

  };   // end class IsogeometricModel

} // end namespace Go


#endif    // #ifndef __ISOGEOMETRICMODEL_H
