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

#include "GoTools/geometry/GridUtilities.h"
#include "GoTools/geometry/PointCloud.h"
#include <iostream>
#include <string.h>
#include <fstream>

using namespace Go;
using std::vector;
using std::pair;

void GridUtilities::grid2PointCloud(double* corners, double nodata_val,
				   uint16_t **valZ, uint cols, uint rows,
				   PointCloud3D& points)
{
  // corners are given as lower left, upper right
  double del1 = (corners[2]-corners[0])/cols; //deltaX
  double del2 = (corners[3]-corners[1])/rows; //deltaY

  uint ki, kj;
  double x, y, z;
  vector<double> data;
  for (ki=0, y=corners[1]; ki<rows; ++ki, y+=del2)
    for (kj=0, x=corners[0]; kj<cols; ++kj, x+=del1)
      {
	z = (double)valZ[ki][kj];
	if (z == nodata_val)
	  continue;
	data.push_back(x);
	data.push_back(y);
	data.push_back(z);
      }

  points = PointCloud3D(data.begin(), data.size()/3);
}

void GridUtilities::grid2PointCloud(double* corners, double nodata_val,
				   short **valZ, uint cols, uint rows,
				   PointCloud3D& points)
{
  // corners are given as lower left, upper right
  double del1 = (corners[2]-corners[0])/cols; //deltaX
  double del2 = (corners[3]-corners[1])/rows; //deltaY

  uint ki, kj;
  double x, y, z;
  vector<double> data;
  for (ki=0, y=corners[1]; ki<rows; ++ki, y+=del2)
    for (kj=0, x=corners[0]; kj<cols; ++kj, x+=del1)
      {
	z = (double)valZ[ki][kj];
	if (z == nodata_val)
	  continue;
	data.push_back(x);
	data.push_back(y);
	data.push_back(z);
      }

  points = PointCloud3D(data.begin(), data.size()/3);
}

void GridUtilities::fillGridVals(float nodata_val, float **valZ, 
				uint cols, uint rows, int nmb)
{
  uint ki, kj, kr, kh;
  float z0, z1, z2;

  // Row wise
  for (ki=0; ki<rows; ++ki)
    {
      z0 = nodata_val;
      for (kj=0; kj<cols; kj=kr)
	{
	  kr = kj + 1;
	  z1 = valZ[ki][kj];
	  if (z0 != nodata_val && z1 == nodata_val)
	    {
	      // Count number of nodata values
	      z2 = nodata_val;
	      for (kr=kj+1; kr<cols; ++kr)
		{
		  if (kr-kj > nmb)
		    break;
		  z2 = (double)valZ[ki][kr];
		  if (z2 != nodata_val)
		    break;
		}
	      if (z2 != nodata_val && kr-kj <= nmb)
		{
		  // Construct z-values
		  float del = (float)(kr-kj+1);
		  for (kh=kj; kh<kr; ++kh)
		    valZ[ki][kh] = ((kr-kh)*z0 + (kh-kj+1)*z2)/del;
		}
	      z1 = z2;
	    }
	  z0 = z1;
	}
    }

  // Column wise
  for (ki=0; ki<cols; ++ki)
    {
      z0 = nodata_val;
      for (kj=0; kj<rows; kj=kr)
	{
	  kr = kj + 1;
	  z1 = valZ[kj][ki];
	  if (z0 != nodata_val && z1 == nodata_val)
	    {
	      // Count number of nodata values
	      z2 = nodata_val;
	      for (kr=kj+1; kr<rows; ++kr)
		{
		  if (kr-kj > nmb)
		    break;
		  z2 = (double)valZ[kr][ki];
		  if (z2 != nodata_val)
		    break;
		}
	      if (z2 != nodata_val && kr-kj <= nmb)
		{
		  // Construct z-values
		  float del = (float)(kr-kj+1);
		  for (kh=kj; kh<kr; ++kh)
		    valZ[kh][ki] = ((kr-kh)*z0 + (kh-kj+1)*z2)/del;
		}
	      z1 = z2;
	    }
	  z0 = z1;
	}
    }

}

void GridUtilities::grid2PointCloud(double* corners, double nodata_val,
				   float **valZ, uint cols, uint rows,
				   PointCloud3D& points)
{
  double del1 = (corners[2]-corners[0])/cols; //deltaX
  double del2 = (corners[3]-corners[1])/rows; //deltaY

  uint ki, kj;
  double x, y, z;
  vector<double> data;
  for (ki=0, y=corners[3]-del2; ki<rows; ++ki, y-=del2)
    for (kj=0, x=corners[0]; kj<cols; ++kj, x+=del1)
      {
	z = (double)valZ[ki][kj];
	if (z == nodata_val)
	  continue;
	data.push_back(x);
	data.push_back(y);
	data.push_back(z);
      }

  points = PointCloud3D(data.begin(), data.size()/3);
}

void GridUtilities::extractLimitIxs(double nodata_val,
				   uint16_t **valZ, uint cols, uint rows,
				   vector<vector<uint> >& limits)
{
  uint ki, kj;
  uint16_t z;
  for (ki=0; ki<rows; ++ki)
    {
      bool prev_ix = false;
      vector<uint> lim_ix;
      for (kj=0; kj<cols; ++kj)
	{
	  z = valZ[ki][kj];
	  if ((double)z != nodata_val)
	    {
	      // Check neighbouring values
	      if (ki == 0 || ki == rows-1)
		lim_ix.push_back(kj);
	      else if (prev_ix == false)
		{
		  lim_ix.push_back(kj);
		  prev_ix = true;
		}
	      else if (kj == cols-1)
		lim_ix.push_back(kj);
	      else if ((double)valZ[ki][kj+1] == nodata_val)
		lim_ix.push_back(kj);
	      else if ((double)valZ[ki-1][kj] == nodata_val)
		lim_ix.push_back(kj);
	      else if ((double)valZ[ki+1][kj] == nodata_val)
		lim_ix.push_back(kj);
	    }
	  else
	    prev_ix = false;
	}
	  limits.push_back(lim_ix);
    }
}

void GridUtilities::extractLimitIxs(double nodata_val,
				   short **valZ, uint cols, uint rows,
				   vector<vector<uint> >& limits)
{
  uint ki, kj;
  short z;
  for (ki=0; ki<rows; ++ki)
    {
      bool prev_ix = false;
      vector<uint> lim_ix;
      for (kj=0; kj<cols; ++kj)
	{
	  z = valZ[ki][kj];
	  if ((double)z != nodata_val)
	    {
	      // Check neighbouring values
	      if (ki == 0 || ki == rows-1)
		lim_ix.push_back(kj);
	      else if (prev_ix == false)
		{
		  lim_ix.push_back(kj);
		  prev_ix = true;
		}
	      else if (kj == cols-1)
		lim_ix.push_back(kj);
	      else if ((double)valZ[ki][kj+1] == nodata_val)
		lim_ix.push_back(kj);
	      else if ((double)valZ[ki-1][kj] == nodata_val)
		lim_ix.push_back(kj);
	      else if ((double)valZ[ki+1][kj] == nodata_val)
		lim_ix.push_back(kj);
	    }
	  else
	    prev_ix = false;
	}
	  limits.push_back(lim_ix);
    }
}

void GridUtilities::extractLimitIxs(double nodata_val,
				   float **valZ, uint cols, uint rows,
				   vector<vector<uint> >& limits)
{
  uint ki, kj;
  float z;
  for (ki=0; ki<rows; ++ki)
    {
      bool prev_ix = false;
      vector<uint> lim_ix;
      for (kj=0; kj<cols; ++kj)
	{
	  z = valZ[ki][kj];
	  if ((double)z != nodata_val)
	    {
	      // Check neighbouring values
	      if (ki == 0 || ki == rows-1)
		lim_ix.push_back(kj);
	      else if (prev_ix == false)
		{
		  lim_ix.push_back(kj);
		  prev_ix = true;
		}
	      else if (kj == cols-1)
		lim_ix.push_back(kj);
	      else if ((double)valZ[ki][kj+1] == nodata_val)
		lim_ix.push_back(kj);
	      else if ((double)valZ[ki-1][kj] == nodata_val)
		lim_ix.push_back(kj);
	      else if ((double)valZ[ki+1][kj] == nodata_val)
		lim_ix.push_back(kj);
	    }
	  else
	    prev_ix = false;
	}
	  limits.push_back(lim_ix);
    }
}

void GridUtilities::getLimitSeqs(vector<vector<uint> >& limits,
				vector<vector<uint> >& seqs)
{
  uint ki, kj, kr;
  size_t kh, kv;

  // Make intial sequences
  // First row
  for (kj=0; kj<limits[0].size(); kj=kr)
    {
      vector<uint> curr;
      for (kr=kj; kr<limits[0].size(); ++kr)
	{
	  if (kr > kj && limits[0][kr]-limits[0][kr-1] > 1)
	    break;
	  curr.push_back(0);
	  curr.push_back(limits[0][kr]);
	}
      seqs.push_back(curr);
    }

  // Remaining rows
  for (ki=1; ki<limits.size(); ++ki)
    {
      for (kj=0; kj<limits[ki].size(); kj=kr)
	{
	  vector<uint> curr;
	  for (kr=kj; kr<limits[ki].size(); ++kr)
	    {
	      if (kr > kj && limits[ki][kr]-limits[ki][kr-1] > 1)
		break;
	      curr.push_back(ki);
	      curr.push_back(limits[ki][kr]);
	    }

	  // Check if the new sequence can be attached one or two
	  // existing ones
	  for (int dir=0; dir<2; ++dir)
	    {
	      // Try both orientations of the current vector
	      if (curr.size() > 0)
		{
		  for (kh=0; kh<seqs.size(); ++kh)
		    if (ki-seqs[kh][seqs[kh].size()-2] == 1 &&
			abs(curr[1]-seqs[kh][seqs[kh].size()-1]) <= 1)
		      {
			if (curr.size() == 4 && 
			    abs(curr[1]-seqs[kh][seqs[kh].size()-1]) == 1 &&
			    abs(curr[3]-seqs[kh][seqs[kh].size()-1]) == 0)
			  std::swap(curr[1], curr[3]);
			seqs[kh].insert(seqs[kh].end(), curr.begin(), curr.end());
			curr.clear();
			break;
		      }
		}
	      if (curr.size() > 0)
		{
		  for (kh=0; kh<seqs.size(); ++kh)
		    if (ki-seqs[kh][0] == 1 &&
			abs(curr[curr.size()-1]-seqs[kh][1]) <= 1)
		      {
			if (curr.size() == 4 && 
			    abs(curr[3]-seqs[kh][1]) == 1 &&
			    abs(curr[1]-seqs[kh][1]) == 0)
			  std::swap(curr[1], curr[3]);
			curr.insert(curr.end(), seqs[kh].begin(), seqs[kh].end());
			seqs[kh].swap(curr);
			curr.clear();
			break;
		      }
		}
	      if (curr.size() <= 2)
		break;

	      // Turn the sequence of the points in curr
	      int size = (int)curr.size();
	      int nmb = size/4;
	      for (int ix=0; ix<nmb; ++ix)
		std::swap(curr[2*ix+1], curr[size-2*(ix+1)+1]);
	    }

	  if (curr.size() > 0)
	    seqs.push_back(curr);
	}

      // Check if any sequences can be joined
      for (kh=0; kh<seqs.size();)
      	{
      	  for (kv=kh+1; kv<seqs.size();)
      	    {
      	      if (abs(seqs[kh][seqs[kh].size()-2]-seqs[kv][0]) == 1 &&
      		       abs(seqs[kh][seqs[kh].size()-1]-seqs[kv][1]) <= 1)
      		{
      		  seqs[kh].insert(seqs[kh].end(), seqs[kv].begin(),
      				  seqs[kv].end());
      		  seqs.erase(seqs.begin()+kv);
      		  break;
      		}
      	      else if (abs(seqs[kv][seqs[kv].size()-2]-seqs[kh][0]) == 1 &&
      		       abs(seqs[kv][seqs[kv].size()-1]-seqs[kh][1]) <= 1)
      		{
      		  seqs[kv].insert(seqs[kv].end(), seqs[kh].begin(),
      				  seqs[kh].end());
      		  std::swap(seqs[kh], seqs[kv]);
      		  seqs.erase(seqs.begin()+kv);
      		  break;
      		}
      	      else if (abs(seqs[kh][0]-seqs[kv][0]) == 1 &&
      		  abs(seqs[kh][1]-seqs[kv][1]) <= 1)
      		{
      		  int size = (int)seqs[kh].size();
      		  int nmb = size/4;
      		  for (int ix=0; ix<nmb; ++ix)
      		    {
      		      std::swap(seqs[kh][2*ix], seqs[kh][size-2*(ix+1)]);
      		      std::swap(seqs[kh][2*ix+1], seqs[kh][size-2*(ix+1)+1]);
      		    }
      		  seqs[kh].insert(seqs[kh].end(), seqs[kv].begin(),
      				  seqs[kv].end());
      		  seqs.erase(seqs.begin()+kv);
      		  break;
      		}
      	      else if (abs(seqs[kh][seqs[kh].size()-2]-
      			   seqs[kv][seqs[kv].size()-2]) == 1 &&
      		       abs(seqs[kh][seqs[kh].size()-1]-
      			   seqs[kv][seqs[kv].size()-1]) <= 1)
      		{
      		  int size = (int)seqs[kv].size();
      		  int nmb = size/4;
      		  for (int ix=0; ix<nmb; ++ix)
      		    {
      		      std::swap(seqs[kv][2*ix], seqs[kv][size-2*(ix+1)]);
      		      std::swap(seqs[kv][2*ix+1], seqs[kv][size-2*(ix+1)+1]);
      		    }
      		  seqs[kh].insert(seqs[kh].end(), seqs[kv].begin(),
      				  seqs[kv].end());
      		  seqs.erase(seqs.begin()+kv);
      		  break;
      		}
      	      else 
      		++kv;
      	    }
      	  if (kv == seqs.size())
      	    kh++;
      	}
    }

  // One or more closed loops are expected. Make a final attempt to join
  // sequences. 
  // Identify closest endpoints
  vector<pair<size_t, size_t> > end_idx;
  for (kh=0; kh<seqs.size(); ++kh)
    {
      for (ki=0, kj=0; ki<2; ++ki, kj=seqs[kh].size()-2)
	{
	  for (kr=0; kr<end_idx.size(); ++kr)
	    if (end_idx[kr].first == 2*kh+ki || end_idx[kr].second == 2*kh+ki)
	      break;
	  if (kr < end_idx.size())
	    continue;

	  size_t ix = 2*kh+1;
	  int mind = abs(seqs[kh][0] - seqs[kh][seqs[kh].size()-2]) +
	    abs(seqs[kh][1] - seqs[kh][seqs[kh].size()-1]);
	  for (kv=kh+1; kv<seqs.size(); ++kv)
	    {
	      int d1 = abs(seqs[kh][kj] - seqs[kv][0]) +
		abs(seqs[kh][kj+1] - seqs[kv][1]);
	      int d2 = abs(seqs[kh][kj] - seqs[kv][seqs[kv].size()-2]) +
		abs(seqs[kh][kj+1] - seqs[kv][seqs[kv].size()-1]);

	      for (kr=0; kr<end_idx.size(); ++kr)
		{
		  if (end_idx[kr].first == 2*kv || end_idx[kr].second == 2*kv)
		    d1 = 2*mind;
		  if (end_idx[kr].first == 2*kv+1 || end_idx[kr].second == 2*kv+1)
		    d2 = 2*mind;
		}

	      if (d1 < mind)
		{
		  mind = d1;
		  ix = 2*kv;
		}
	      if (d2 < mind)
		{
		  mind = d2;
		  ix = 2*kv + 1;
		}
	    }
	  end_idx.push_back(std::make_pair(2*kh+ki, ix));
	  if (ix - 2*kh == 1)
	    break;
	}
    }

  // Reorganize points towards the seqence ends if necessary
  for (ki=0; ki<end_idx.size(); ++ki)
    {
      for (int iend=0; iend<2; ++iend)
	{
	  // Check if the seqence of the last points must be swapped
	  size_t ix1 = (iend == 0) ? end_idx[ki].first : end_idx[ki].second;
	  size_t ix2 = (iend == 0) ? end_idx[ki].second : end_idx[ki].first;
	  size_t start = (ix2%2 == 0) ? 0 : seqs[ix2/2].size()-2;
	  size_t stop =  (ix2%2 == 0) ? seqs[ix2/2].size()-4 : 2;
	  int sgn = (ix2%2 == 0) ? 1 : -1;
	  kr = (ix1%2 == 0) ? 0 : seqs[ix1/2].size()-2;
	  ix1 = ix1/2;
	  ix2 = ix2/2;
	  if (seqs[ix2].size() < 6)
	    continue;    // No room for reorganization
	  int d2 = abs(seqs[ix1][kr] - seqs[ix2][start]) +
	    abs((seqs[ix1][kr+1] - seqs[ix2][start+1]));
	  for (kj=start+2*sgn; sgn*(int)kj<sgn*(int)stop; kj+=2*sgn)
	    {
	      int d1 = abs(seqs[ix1][kr] - seqs[ix2][kj]) +
		abs((seqs[ix1][kr+1] - seqs[ix2][kj+1]));
	      if (d1 < d2)
		{
		  std::swap(seqs[ix2][kj], seqs[ix2][start]);
		  std::swap(seqs[ix2][kj+1], seqs[ix2][start+1]);
		  d2 = d1;
		}
	      else
		{
		  int size = abs((int)kj - (int)start - 2*sgn);
		  int nmb = size/4;
		  int ka;
		  for (kh=start+2*sgn, ka=0; ka<nmb; ++ka, kh+=2*sgn)
		    {
		      std::swap(seqs[ix2][kh], seqs[ix2][kj-sgn*2*(ka+1)]);
		      std::swap(seqs[ix2][kh+1], seqs[ix2][kj-sgn*2*(ka+1)+1]);
		    } 
		  break;
		}
	    }
	}
    }

  // Join sequences
  for (ki=0; ki<end_idx.size(); ++ki)
    {
      if (abs(end_idx[ki].second - end_idx[ki].first) > 1)
	{
	  size_t ix1 = end_idx[ki].first/2;
	  size_t ix2 = end_idx[ki].second/2;
	  size_t ix3, endix;
	  if (end_idx[ki].first%2 == 1 && end_idx[ki].second%2 == 0)
	    {
	      seqs[ix1].insert(seqs[ix1].end(), seqs[ix2].begin(), seqs[ix2].end());
	      seqs.erase(seqs.begin()+ix2);
	      ix3 = ix2;
	      endix = end_idx[ki].first;
	    }
	  else if (end_idx[ki].first%2 == 0 && end_idx[ki].second%2 == 1)
	    {
	      seqs[ix2].insert(seqs[ix2].end(), seqs[ix1].begin(), seqs[ix1].end());
	      seqs.erase(seqs.begin()+ix1);
	      ix3 = ix1;
	      endix = end_idx[ki].second - 2;
	    }
	  else if (end_idx[ki].first%2 == 0 && end_idx[ki].second%2 == 0)
	    {
	      int size = (int)seqs[ix1].size();
	      int nmb = size/4;
	      for (int ix=0; ix<nmb; ++ix)
		{
		  std::swap(seqs[ix1][2*ix], seqs[ix1][size-2*(ix+1)]);
		  std::swap(seqs[ix1][2*ix+1], seqs[ix1][size-2*(ix+1)+1]);
		}
  	      seqs[ix1].insert(seqs[ix1].end(), seqs[ix2].begin(), seqs[ix2].end());
	      seqs.erase(seqs.begin()+ix2);
	      end_idx[ki].first = 1;
	      ix3 = ix2;
	      endix = end_idx[ki].first + 1;
	    }
	  else if (end_idx[ki].first%2 == 1 && end_idx[ki].second%2 == 1)
	    {
	      int size = (int)seqs[ix2].size();
	      int nmb = size/4;
	      for (int ix=0; ix<nmb; ++ix)
		{
		  std::swap(seqs[ix2][2*ix], seqs[ix2][size-2*(ix+1)]);
		  std::swap(seqs[ix2][2*ix+1], seqs[ix2][size-2*(ix+1)+1]);
		}
  	      seqs[ix2].insert(seqs[ix2].end(), seqs[ix1].begin(), seqs[ix1].end());
	      seqs.erase(seqs.begin()+ix1);
	      end_idx[ki].second = 0;
	      ix3 = ix1;
	      endix = end_idx[ki].first;
	    }

	  for (kj=ki+1; kj<end_idx.size(); ++kj)
	    {
	      if (end_idx[kj].first > 2*ix3+1)
		end_idx[kj].first -= 2;
	      if (end_idx[kj].second > 2*ix3+1)
		end_idx[kj].second -= 2;
	      if (end_idx[kj].first == 2*ix3 || end_idx[kj].first == 2*ix3+1)
		{
		  end_idx[kj].first = endix;
		}
	      if (end_idx[kj].second == 2*ix3 || end_idx[kj].second == 2*ix3+1)
		{
		  end_idx[kj].second = endix;
		}
	    }
	}
    }
  
#ifdef DEBUG
  std::ofstream of("trim_seq.g2");
  for (ki=0; ki<seqs.size(); ++ki)
    {
      of << "400 1 0 4 255 0 0 255" << std::endl;
      of << 1 << std::endl;
      of << Point(seqs[ki][0], seqs[ki][1], 0) << std::endl << std::endl;

      of << "410 1 0 4 0 255 0 255" << std::endl;
      of << seqs[ki].size()/2-1 << std::endl;
      for (kj=2; kj<seqs[ki].size(); kj+=2)
	{
	  of << Point(seqs[ki][kj-2], seqs[ki][kj-1], 0) << " ";
	  of << Point(seqs[ki][kj], seqs[ki][kj+1], 0) << std::endl;
	}
      of << std::endl;
    }
#endif
}

void GridUtilities::seq2Points(double* corners, 
			      uint16_t **valZ, uint cols, uint rows,
			      vector<vector<uint> >& seqs,
			      vector<vector<double> >& point_seqs)
{
  // corners are given as lower left, upper right
  double del1 = (corners[2]-corners[0])/cols; //deltaX
  double del2 = (corners[3]-corners[1])/rows; //deltaY

  for (size_t kr=0; kr<seqs.size(); ++kr)
    {
      uint ki, kj;
      double x, y, z;
      vector<double> pts;
      for (size_t kh=0; kh<seqs[kr].size(); kh+=2)
	{
	  ki = seqs[kr][kh];
	  kj = seqs[kr][kh+1];
	  x = corners[0] + kj*del1;
	  y = corners[1] + ki*del2;
	  z = (double)valZ[ki][kj];

	  pts.push_back(x);
	  pts.push_back(y);
	  pts.push_back(z);
	}

      point_seqs.push_back(pts);
    }
}

void GridUtilities::seq2Points(double* corners, 
			      short **valZ, uint cols, uint rows,
			      vector<vector<uint> >& seqs,
			      vector<vector<double> >& point_seqs)
{
  // corners are given as lower left, upper right
  double del1 = (corners[2]-corners[0])/cols; //deltaX
  double del2 = (corners[3]-corners[1])/rows; //deltaY

  for (size_t kr=0; kr<seqs.size(); ++kr)
    {
      uint ki, kj;
      double x, y, z;
      vector<double> pts;
      for (size_t kh=0; kh<seqs[kr].size(); kh+=2)
	{
	  ki = seqs[kr][kh];
	  kj = seqs[kr][kh+1];
	  x = corners[0] + kj*del1;
	  y = corners[1] + ki*del2;
	  z = (double)valZ[ki][kj];

	  pts.push_back(x);
	  pts.push_back(y);
	  pts.push_back(z);
	}

      point_seqs.push_back(pts);
    }
}

void GridUtilities::seq2Points(double* corners, 
			      float **valZ, uint cols, uint rows,
			      vector<vector<uint> >& seqs,
			      vector<vector<double> >& point_seqs)
{
  // corners are given as lower left, upper right
  double del1 = (corners[2]-corners[0])/cols; //deltaX
  double del2 = (corners[3]-corners[1])/rows; //deltaY

  for (size_t kr=0; kr<seqs.size(); ++kr)
    {
      uint ki, kj;
      double x, y, z;
      vector<double> pts;
      for (size_t kh=0; kh<seqs[kr].size(); kh+=2)
	{
	  ki = seqs[kr][kh];
	  kj = seqs[kr][kh+1];
	  x = corners[0] + kj*del1;
	  y = corners[1] + ki*del2;
	  z = (double)valZ[ki][kj];

	  pts.push_back(x);
	  pts.push_back(y);
	  pts.push_back(z);
	}

      point_seqs.push_back(pts);
    }
}
