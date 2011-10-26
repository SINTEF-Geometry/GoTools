//===========================================================================
//
// File : CompositeModel.C
//
// Created: Thu Feb 21 09:40:05 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision:
//
// Description:
//
//===========================================================================


#include "GoTools/compositemodel/CompositeModel.h"

namespace Go
{

  //===========================================================================
  CompositeModel::CompositeModel(double gap, double neighbour, double kink, double bend)
    : toptol_(tpTolerances(gap, neighbour, kink, bend))
      //===========================================================================
  {
      closest_idx_ = -1;
  }



  //===========================================================================
  CompositeModel::~CompositeModel()
  //===========================================================================
  {
  }


  //===========================================================================
  void CompositeModel::setTolerances(double gap, double neighbour,
				     double kink, double bend)
  //===========================================================================
  {
    toptol_ = tpTolerances(gap, neighbour, kink, bend);
  }


  //===========================================================================
  double CompositeModel::boxVecDist(const BoundingBox& box, 
				  const Point& vec) const
  //===========================================================================
  {
    Point b = box.low() - vec;
    Point t = vec - box.high();
    double dist2 = 0;
    double d = std::max(0.0, std::max(b[0], t[0]));
    dist2 += d*d;
    d = std::max(0.0, std::max(b[1], t[1]));
    dist2 += d*d;
    d = std::max(0.0, std::max(b[2], t[2]));
    dist2 += d*d;
    return sqrt(dist2);
  }

//===========================================================================
bool
CompositeModel::boxExtreme(const BoundingBox& box, const Point& dir, 
			 const Point& curr_pnt) const
//===========================================================================
{
  // Fetch box corners
  Point corners[8];
  corners[0] = corners[1] = corners[2] = corners[3] = box.low();
  corners[4] = corners[5] = corners[6] = corners[7] = box.high();
  corners[1][0] = corners[7][0];
  corners[2][1] = corners[7][1];
  corners[3][0] = corners[7][0];
  corners[3][1] = corners[7][1];
  corners[4][0] = corners[0][0];
  corners[4][1] = corners[0][1];
  corners[5][1] = corners[0][1];
  corners[6][0] = corners[0][0];

  // For all corners, check if it lies further in the given
  // direction than the input point
  for (int ki=0; ki<8; ki++)
    {
      if (corners[ki]*dir > curr_pnt*dir)
	return true;
    }

  return false;
}

} // namespace Go
