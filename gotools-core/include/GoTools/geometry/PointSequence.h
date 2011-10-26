//===========================================================================
//
// File : PointSequence.h
//
// Created: Thu Oct 29 14:23:42 2009
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================



#ifndef __POINTSEQUENCE_H
#define __POINTSEQUENCE_H


#include <vector>
#include <algorithm>


namespace Go
{


  enum PointSequenceType
  {
    PSTPoint,         // Each grid position holds a point
    PSTPointTangent,  // Each grid position holds a point, then a tangent
    PSTScattered      // Regard this an unordered point cloud
  };


  class PointSequence
  {

  public:

    // Constructors
    PointSequence()
      : dim_(-1), type_(PSTPoint)
    {
      grid_length_.resize(0);
    }


    template <typename RandomIterator> PointSequence(int dim,
						     int nmb_pts,
						     RandomIterator coefs_start,
						     PointSequenceType pst = PSTPoint)
      : dim_(dim), type_(pst)
    {
      grid_length_.resize(1);
      grid_length_[0] = nmb_pts;
      int coefs_size;
      switch (pst)
	{
	case PSTPoint:
	  coefs_size = nmb_pts * dim;
	  break;
	case PSTPointTangent:
	  coefs_size = 2 * nmb_pts * dim;
	  break;
	}
      coefs_.resize(coefs_size);
      copy(coefs_start, coefs_start + coefs_size, coefs_.begin());
    }



    template <typename RandomIterator> PointSequence(int dim,
						     int nmb_pts_1,
						     int nmb_pts_2,
						     RandomIterator coefs_start,
						     PointSequenceType pst = PSTPoint)
      : dim_(dim), type_(pst)
    {
      grid_length_.resize(2);
      grid_length_[0] = nmb_pts_1;
      grid_length_[1] = nmb_pts_2;
      int coefs_size;
      switch (pst)
	{
	case PSTPoint:
	  coefs_size = nmb_pts_1 * nmb_pts_2 * dim;
	  break;
	case PSTPointTangent:
	  coefs_size = 2 * nmb_pts_1 * nmb_pts_2 * dim;
	  break;
	}
      coefs_.resize(coefs_size);
      copy(coefs_start, coefs_start + coefs_size, coefs_.begin());
    }




    template <typename RandomIterator> PointSequence(int dim,
						     int nmb_pts_1,
						     int nmb_pts_2,
						     int nmb_pts_3,
						     RandomIterator coefs_start,
						     PointSequenceType pst = PSTPoint)
      : dim_(dim), type_(pst)
    {
      grid_length_.resize(3);
      grid_length_[0] = nmb_pts_1;
      grid_length_[1] = nmb_pts_2;
      grid_length_[2] = nmb_pts_3;
      int coefs_size;
      switch (pst)
	{
	case PSTPoint:
	  coefs_size = nmb_pts_1 * nmb_pts_2 * nmb_pts_3 * dim;
	  break;
	case PSTPointTangent:
	  coefs_size = 2 * nmb_pts_1 * nmb_pts_2 * nmb_pts_3 * dim;
	  break;
	}
      coefs_.resize(coefs_size);
      copy(coefs_start, coefs_start + coefs_size, coefs_.begin());
    }


    // Destructor
    virtual ~PointSequence() { }


    int dimension() const;

    int grid_dimension() const;

    int grid_length(int i) const;

    PointSequenceType type() const;

    std::vector<double>::iterator coefs_begin()
    { return coefs_.begin(); }

    std::vector<double>::iterator coefs_end()
    { return coefs_.end(); }

    std::vector<double>::const_iterator coefs_begin() const
    { return coefs_.begin(); }

    std::vector<double>::const_iterator coefs_end() const
    { return coefs_.end(); }


  private:

    std::vector<double> coefs_;      // The coordinates of the points
    int dim_;                        // The dimension of the space of points
    std::vector<int> grid_length_;   // Number of points in each direction in grid. grid_length_.size()
                                     // gives the grid dimension (flat, 2D, 3D, etc)
    PointSequenceType type_;

  };    // Class PointSequence


} // namespace Go


#endif    // #ifndef __POINTSEQUENCE_H
