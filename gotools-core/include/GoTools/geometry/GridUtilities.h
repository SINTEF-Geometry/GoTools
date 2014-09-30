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

#ifndef _GRIDUTILITIES_H
#define _GRIDUTILITIES_H

#include "GoTools/geometry/PointCloud.h"
#include <vector>
#include <inttypes.h>

using namespace Go;

namespace GridUtilities
{
  void grid2PointCloud(double* corners, double nodata_val,
		       uint16_t **valZ, uint cols, uint rows,
		       PointCloud3D& points);

  void grid2PointCloud(double* corners, double nodata_val,
		       short **valZ, uint cols, uint rows,
		       PointCloud3D& points);

  void grid2PointCloud(double* corners, double nodata_val,
  		       float **valZ, uint cols, uint rows,
  		       PointCloud3D& points);

  void fillGridVals(float nodata_val, float **valZ, 
		    uint cols, uint rows, int nmb);

 void extractLimitIxs(double nodata_val,
		       uint16_t **valZ, uint cols, uint rows,
		       std::vector<std::vector<uint> >& limits);

 void extractLimitIxs(double nodata_val,
		       short **valZ, uint cols, uint rows,
		       std::vector<std::vector<uint> >& limits);

  void extractLimitIxs(double nodata_val,
		       float **valZ, uint cols, uint rows,
		       std::vector<std::vector<uint> >& limits);

  void getLimitSeqs(std::vector<std::vector<uint> >& limits,
		    std::vector<std::vector<uint> >& seqs);

  void seq2Points(double* corners, 
		  uint16_t **valZ, uint cols, uint rows,
		  std::vector<std::vector<uint> >& seqs,
		  std::vector<std::vector<double> >& point_seqs);

  void seq2Points(double* corners, 
		  short **valZ, uint cols, uint rows,
		  std::vector<std::vector<uint> >& seqs,
		  std::vector<std::vector<double> >& point_seqs);

  void seq2Points(double* corners, 
		  float **valZ, uint cols, uint rows,
		  std::vector<std::vector<uint> >& seqs,
		  std::vector<std::vector<double> >& point_seqs);
}


#endif
