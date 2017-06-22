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

#ifndef _FACECONNECTIVITY_H_
#define _FACECONNECTIVITY_H_

#include <vector>

using std::vector;
using std::pair;

namespace Go
{

/// A structure storing the connectivity information between
/// two adjacent faces
template <class edgeType>
struct FaceConnectivity
{
public:
  /// Edge of first face along the common boundary
  edgeType *e1_;
  /// Edge of second face along the common boundary
  edgeType *e2_;
  /// The status is:
  /// 0 : edges join smoothly. \f$G^1\f$.
  /// 1 : edges join, but normals are slightly discontinous. A kink.
  /// 2 : edges join, but the normals are discontinous. \f$G^0\f$.
  /// 3 : edges almost join. A gap.
  /// 4 : edges are totally discontinous.
  /// The minimal tpTopologicalInfo has a one-element status vector
  /// and a two-element parameters vector
  vector<int> status_;
  /// Parameter intervals limiting the areas of the found state of continuity 
  vector< pair<double, double> > parameters_; // The pair refers to parameter of e1_ & e2_.

  /// Constructor
  FaceConnectivity(edgeType* e1, edgeType *e2)
  {
    e1_ = e1;
    e2_ = e2;
  }

  /// Reset edge info
  void setEdges(edgeType* e1, edgeType *e2)
  {
    e1_ = e1;
    e2_ = e2;
  }

  /// The highest continuity between the two faces
  int BestStatus() const
  {
    int s = 4;
    for (size_t i = 0; i < status_.size(); ++i)
      if (status_[i] < s)
	s = status_[i];
    return s;
  }
  /// The lowest continuity between the two faces
  int WorstStatus() const
  {
    int s = 0;
    for (size_t i = 0; i < status_.size(); ++i)
      if (status_[i] > s) s = status_[i];
    return s;
  }
};

} // namespace Go
#endif  //_FACECONNECTIVITY_H_
