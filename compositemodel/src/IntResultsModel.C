//===========================================================================
//
// File : IntResultsModel.C
//
// Created: July 9, 2009
//
// Author: Vibeke Skytt, SINTEF
//
// Revision: 
//
// Description:
//
//===========================================================================


#include "GoTools/compositemodel/IntResultsModel.h"

using namespace std;

namespace Go
{
  //===========================================================================
  // Constructor
  IntResultsModel::IntResultsModel(IntersectionType type)
  //===========================================================================
    : type_(type)
  {
    numpar_ = 0;
    if (type_ == SurfaceModel_SurfaceModel)
      numpar_ = 4;
    else if (type_ == SurfaceModel_Plane || type_ == SurfaceModel_Line ||
	     type_ == CompositeCurve_CompositeCurve)
      numpar_ = 2;
    else
      numpar_ = 1;
  }

} // namespace Go
