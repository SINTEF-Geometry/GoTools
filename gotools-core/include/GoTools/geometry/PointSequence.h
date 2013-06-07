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

#ifndef __POINTSEQUENCE_H
#define __POINTSEQUENCE_H


#include <vector>
#include <algorithm>


namespace Go
{


  /// Type of point
  enum PointSequenceType
  {
    /// Each grid position holds a point
    PSTPoint,         
    /// Each grid position holds a point, then a tangent
    PSTPointTangent,  
    /// Regard this an unordered point cloud
    PSTScattered      
  };


  class PointSequence
  {

  public:

    /// Default constructor
    PointSequence()
      : dim_(-1), type_(PSTPoint)
    {
      grid_length_.resize(0);
    }


    /// Constructor given a sequence of points of a specified type
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



    /// Constructor given a number of sequences of points of a specified type,
    /// Each point sequence must contain the same number of points, i.e.
    /// two dimensional grid
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



    /// Constructor given a 3 dimensional grid of points of a specified type
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


    /// Point dimension
    int dimension() const;

    /// Dimension of grid
    int grid_dimension() const;

    /// Length of grid in the given direction
    int grid_length(int i) const;

    /// Point type
    PointSequenceType type() const;

    /// Iterator to grid start
    std::vector<double>::iterator coefs_begin()
    { return coefs_.begin(); }

     /// Iterator to grid end
   std::vector<double>::iterator coefs_end()
    { return coefs_.end(); }

    /// Constant iterator to grid start
    std::vector<double>::const_iterator coefs_begin() const
    { return coefs_.begin(); }

    /// Constant iterator to grid end
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
