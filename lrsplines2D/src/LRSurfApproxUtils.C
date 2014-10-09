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


#include "GoTools/lrsplines2D/LRSurfApproxUtils.h"
#include "GoTools/geometry/Utils.h"
#include <iostream>
#include <fstream>

using namespace Go;
using std::vector;

#define DEBUG

int compare_u_par(const void* el1, const void* el2)
{
  if (((double*)el1)[0] < ((double*)el2)[0])
    return -1;
  else if (((double*)el1)[0] > ((double*)el2)[0])
    return 1;
  else
    return 0;
}

int compare_v_par(const void* el1, const void* el2)
{
  if (((double*)el1)[1] < ((double*)el2)[1])
    return -1;
  else if (((double*)el1)[1] > ((double*)el2)[1])
    return 1;
  else
    return 0;
}

//==============================================================================
void LRSurfApproxUtils::computeTrimInfo(vector<double>& points, int dim,
					int max_level, int nmb_u, int nmb_v,
					vector<double>& seq)
//==============================================================================
{
  int del = dim+2;
  int nmb = (int)points.size()/del;  // Number of data points

  // Compute the domain corresponding to the point set
  int ki, kj, kh;
  double minmax[4];
  for (ki=0; ki<2; ++ki)
    minmax[2*ki] = minmax[2*ki+1] = points[ki];
  for (kj=1; kj<nmb; ++kj)
    for (ki=0; ki<2; ++ki)
      {
	minmax[2*ki] = std::min(minmax[2*ki], points[kj*del+ki]);
	minmax[2*ki+1] = std::max(minmax[2*ki+1], points[kj*del+ki]);
      }
  
  vector<vector<double> > seqs2;
  int nmb_sub = 0;
  computeTrimInfo2(&points[0], dim, nmb, minmax, max_level, nmb_u, nmb_v, 
		   33, seqs2, nmb_sub);

#ifdef DEBUG
  std::ofstream of0("trim_seqs0.g2");
  (void)of0.precision(15);
  for (ki=0; ki<(int)seqs2.size(); ++ki)
    {
      of0 << "410 1 0 0" << std::endl;
      of0 << seqs2[ki].size()/2-1 << std::endl;
      for (kj=0; kj<=(int)seqs2[ki].size()-4; kj+=2)
	{
	  of0 << seqs2[ki][kj] << " " << seqs2[ki][kj+1] << " " << 0 << " ";
	  of0 << seqs2[ki][kj+2] << " " << seqs2[ki][kj+3] << " " << 0 << std::endl;
	}
    }
#endif

  // Merge sequences with "exact" match
  double eps2 = 1.0e-12; // A small number
  double limit = std::min((minmax[1] - minmax[0])/(double)(nmb_u),
			  (minmax[3] - minmax[2])/(double)(nmb_v));
  limit = eps2;
  for (ki=0; ki<(int)seqs2.size(); )
    {
      for (kj=ki+1; kj<(int)seqs2.size(); ++kj)
	{
	  double dist2 = Utils::distance_squared(seqs2[kj].begin(),
						 seqs2[kj].begin()+2,
						 seqs2[ki].begin()+seqs2[ki].size()-2);
	  if (dist2 < limit)
	    {
	      seqs2[ki].insert(seqs2[ki].end(), seqs2[kj].begin(), seqs2[kj].end());
	      seqs2.erase(seqs2.begin()+kj);
	      break;
	    }

	  dist2 = Utils::distance_squared(seqs2[ki].begin(),
					  seqs2[ki].begin()+2,
					  seqs2[kj].begin()+seqs2[kj].size()-2);
	  if (dist2 < limit)
	    {
	      seqs2[kj].insert(seqs2[kj].end(), seqs2[ki].begin(), seqs2[ki].end());
	      std::swap(seqs2[ki], seqs2[kj]);
	      seqs2.erase(seqs2.begin()+kj);
	      break;
	    }
	}
      if (kj == (int)seqs2.size())
	++ki;
    }

  // Remove duplicate parameter points
  for (kh=0; kh<(int)seqs2.size(); ++kh)
    {
      for (kj=0, ki=2; ki<(int)seqs2[kh].size(); )
	{
	  double dist2 = Utils::distance_squared(seqs2[kh].begin()+kj,
						 seqs2[kh].begin()+ki,
						 seqs2[kh].begin()+ki);
	  if (dist2 < eps2)
	    seqs2[kh].erase(seqs2[kh].begin()+ki, seqs2[kh].begin()+ki+2);
	  else
	    {
	      ki += 2;
	      kj += 2;
	    }
	}
    }

  // Remove fragments
  for (ki=0; ki<(int)seqs2.size(); )
    {
      if (seqs2[ki].size() <= 6)
	seqs2.erase(seqs2.begin()+ki);
      else
	++ki;
    }

#ifdef DEBUG
  std::ofstream of("trim_seqs.g2");
  (void)of.precision(15);
  for (ki=0; ki<(int)seqs2.size(); ++ki)
    {
      of << "410 1 0 0" << std::endl;
      of << seqs2[ki].size()/2-1 << std::endl;
      for (kj=0; kj<=(int)seqs2[ki].size()-4; kj+=2)
	{
	  of << seqs2[ki][kj] << " " << seqs2[ki][kj+1] << " " << 0 << " ";
	  of << seqs2[ki][kj+2] << " " << seqs2[ki][kj+3] << " " << 0 << std::endl;
	}
    }
#endif

  // Merge sequences with closest match
  limit = 2.0*std::min((minmax[1] - minmax[0])/(double)(nmb_u),
		       (minmax[3] - minmax[2])/(double)(nmb_v));
  limit *= limit;
  for (ki=0; ki<(int)seqs2.size()-1; )
    {
      double min1 = limit, min2 = limit;
      int ix1 = -1, ix2 = -1;
      for (kj=ki+1; kj<(int)seqs2.size(); ++kj)
	{
	  double dist0 = (seqs2[kj].size() < 8) ? limit :
	    Utils::distance_squared(seqs2[kj].begin(), seqs2[kj].begin()+2,
				    seqs2[kj].begin()+seqs2[kj].size()-2);

	  double dist2 = Utils::distance_squared(seqs2[kj].begin(),
						 seqs2[kj].begin()+2,
						 seqs2[ki].begin()+seqs2[ki].size()-2);
	  double dist3 = Utils::distance_squared(seqs2[ki].begin(),
					  seqs2[ki].begin()+2,
					  seqs2[kj].begin()+seqs2[kj].size()-2);

	  double dista =  Utils::distance_squared(seqs2[kj].begin(),
						 seqs2[kj].begin()+2,
						 seqs2[ki].begin());
	  double distb = Utils::distance_squared(seqs2[kj].begin()+seqs2[kj].size()-2,
						 seqs2[kj].begin()+seqs2[kj].size(),
						 seqs2[ki].begin()+seqs2[ki].size()-2);

	  if (dist2 < min1 && dist2 < dist0)
	    {
	      min1 = dist2;
	      ix1 = kj;
	    }

	  if (dist3 < min2 && dist3 < dist0)
	    {
	      min2 = dist3;
	      ix2 = kj;
	    }
	}
      if (ix1 >= 0 && min1 < min2)
	{
	  seqs2[ki].insert(seqs2[ki].end(), seqs2[ix1].begin(), seqs2[ix1].end());
	  seqs2.erase(seqs2.begin()+ix1);
	}
      else if (ix2 >= 0)
	{
	  seqs2[ix2].insert(seqs2[ix2].end(), seqs2[ki].begin(), seqs2[ki].end());
	  std::swap(seqs2[ki], seqs2[ix2]);
	  seqs2.erase(seqs2.begin()+ix2);
	}
      else
	++ki;
    }
  seq = seqs2[0];
  
  for (ki=0; ki<(int)seqs2.size()-1; )
    {
      double min1 = HUGE, min2 = HUGE;
      int ix1 = -1, ix2 = -1;
      for (kj=ki+1; kj<(int)seqs2.size(); ++kj)
	{
	  double dist0 = (seqs2[kj].size() < 8) ? limit :
	    Utils::distance_squared(seqs2[kj].begin(), seqs2[kj].begin()+2,
				    seqs2[kj].begin()+seqs2[kj].size()-2);

	  double dist2 = Utils::distance_squared(seqs2[kj].begin(),
						 seqs2[kj].begin()+2,
						 seqs2[ki].begin()+seqs2[ki].size()-2);
	  double dist3 = Utils::distance_squared(seqs2[ki].begin(),
					  seqs2[ki].begin()+2,
					  seqs2[kj].begin()+seqs2[kj].size()-2);

	  double dista =  Utils::distance_squared(seqs2[kj].begin(),
						 seqs2[kj].begin()+2,
						 seqs2[ki].begin());
	  double distb = Utils::distance_squared(seqs2[kj].begin()+seqs2[kj].size()-2,
						 seqs2[kj].begin()+seqs2[kj].size(),
						 seqs2[ki].begin()+seqs2[ki].size()-2);

	  if (dist2 < min1 && dist2 < dist0)
	    {
	      min1 = dist2;
	      ix1 = kj;
	    }

	  if (dist3 < min2 && dist3 < dist0)
	    {
	      min2 = dist3;
	      ix2 = kj;
	    }
	}
      if (ix1 >= 0 && min1 < min2)
	{
	  seqs2[ki].insert(seqs2[ki].end(), seqs2[ix1].begin(), seqs2[ix1].end());
	  seqs2.erase(seqs2.begin()+ix1);
	}
      else if (ix2 >= 0)
	{
	  seqs2[ix2].insert(seqs2[ix2].end(), seqs2[ki].begin(), seqs2[ki].end());
	  std::swap(seqs2[ki], seqs2[ix2]);
	  seqs2.erase(seqs2.begin()+ix2);
	}
      else
	++ki;
    }
  seq = seqs2[0];
  
#ifdef DEBUG
  std::ofstream of2("trim_seqs2.g2");
  (void)of2.precision(15);
  for (ki=0; ki<(int)seqs2.size(); ++ki)
    {
      of2 << "410 1 0 0" << std::endl;
      of2 << seqs2[ki].size()/2-1 << std::endl;
      for (kj=0; kj<=(int)seqs2[ki].size()-4; kj+=2)
	{
	  of2 << seqs2[ki][kj] << " " << seqs2[ki][kj+1] << " " << 0 << " ";
	  of2 << seqs2[ki][kj+2] << " " << seqs2[ki][kj+3] << " " << 0 << std::endl;
	}
    }
#endif

  
}


//==============================================================================
void mergeTrimSeqs(vector<vector<double> >& seqs,
		   vector<vector<double> >& seqs2,
		   int side2)
//==============================================================================
{
  double eps2 = 1.0e-12; // A small number
  for (size_t ka=0; ka<seqs2.size(); ++ka)
    {
      if (seqs2[ka].size() == 0)
	continue;

      int s3 = side2 % 10;
      int s4 = side2 - s3;
      if (seqs.size() == 0)
	seqs.push_back(seqs2[ka]);
      else if (s3 == 1 || (s4 == 20 && s3 != 2))
	{

	  size_t kb;
	  for (kb=0; kb<seqs.size(); ++kb)
	    {
	      double dist2 = 
		Utils::distance_squared(seqs[kb].begin(),
					seqs[kb].begin()+2,
					seqs2[ka].begin()+seqs2[ka].size()-2);
	      if (dist2 < eps2)
		{
		  seqs2[ka].insert(seqs2[ka].end(), 
				   seqs[kb].begin(), 
				   seqs[kb].end());
		  std::swap(seqs[kb], seqs2[ka]);
		  break;
		}
	    }
	  if (kb == seqs.size())
	    seqs.push_back(seqs2[ka]);
	}
      else
	{
	  size_t kb;
	  for (kb=0; kb<seqs.size(); ++kb)
	    {
	      double dist2 = 
		Utils::distance_squared(seqs2[ka].begin(),
					seqs2[ka].begin()+2,
					seqs[kb].begin()+seqs[kb].size()-2);
	      if (dist2 < eps2)
		{
		  seqs[kb].insert(seqs[kb].end(), 
				  seqs2[ka].begin(), 
				  seqs2[ka].end());
		  break;
		}
	    }
	  if (kb == seqs.size())
	    seqs.push_back(seqs2[ka]);
	}
    }
}

//==============================================================================
void LRSurfApproxUtils::computeTrimInfo2(double *points, int dim,
					 int nmb_pts, double minmax[], 
					 int max_level, int nmb_u, 
					 int nmb_v, int side,
					 vector<vector<double> >& seqs,
					 int& nmb_sub)
//==============================================================================
{
  int del = dim+2;

  int s1 = side % 10;
  int s2 = side - s1;
  if (max_level == 0)
    {
#ifdef DEBUG
      std::ofstream of01("sub_cloud1.g2");
      (void)of01.precision(15);
      of01 << "400 1 0 0" << std::endl;
      of01 << nmb_pts << std::endl;
      for (int kr=0; kr<nmb_pts; ++kr)
	{
	  for (int kh=0; kh<del-1; ++kh)
	    of01 << points[kr*del+kh] << " ";
	  of01 << 0 << std::endl;
	}
#endif

      double corner_vec[8];
      corner_vec[0] = minmax[0];
      corner_vec[1] = minmax[3];
      corner_vec[2] = minmax[0];
      corner_vec[3] = minmax[2];
      corner_vec[4] = minmax[1];
      corner_vec[5] = minmax[2];
      corner_vec[6] = minmax[1];
      corner_vec[7] = minmax[3];

      // Fetch corners of domain
      vector<double> curr_seq;
      if (s2 == 0)
	{
	  if (s1 == 1)  // Left edge, top to bottom
	    curr_seq.insert(curr_seq.end(), corner_vec, corner_vec+4);
	  else if (s1 == 2)  // Right edge, bottom to top
	    curr_seq.insert(curr_seq.end(), corner_vec+4, corner_vec+8);
	}
      else if (s1 == 0)
	{
	  if (s2 == 10)  // Bottom edge, left to right
	    curr_seq.insert(curr_seq.end(), corner_vec+2, corner_vec+6);
	  else if (s2 == 20)  // Top edge, right to left
	    {
	      curr_seq.insert(curr_seq.end(), corner_vec+6, corner_vec+8);
	      curr_seq.insert(curr_seq.end(), corner_vec, corner_vec+2);
	    }
	}
      else if ((s2 == 20 || s2 == 30) && s1 == 1)
	{
	  // Top edge, left edge and possibly bottom edge
	  curr_seq.insert(curr_seq.end(), corner_vec+6, corner_vec+8);
	  curr_seq.insert(curr_seq.end(), corner_vec, 
		     corner_vec+((s2==20)?4:6));
	}
      else if (s2 == 10 && (s1 == 1 || s1 == 3))  
	// Left edge, bottom edge and possibly right edge
	curr_seq.insert(curr_seq.end(), corner_vec, corner_vec+((s1==1)?6:8));
      else if (s1 == 2 && s2 == 10) 
	// Bottom edge and right edge
	curr_seq.insert(curr_seq.end(), corner_vec+2, corner_vec+8);
      else if (s1 == 2 && s2 == 30)
	{
	  // Bottom edge, right edge and top edge
	  curr_seq.insert(curr_seq.end(), corner_vec+2, corner_vec+8);
	  curr_seq.insert(curr_seq.end(), corner_vec, corner_vec+2);
	}
      else if (s2 == 20 && (s1 == 2 || s1 == 3))
	{
	  // Right edge, top edge and possibly left edge
	  curr_seq.insert(curr_seq.end(), corner_vec+4, corner_vec+8);
	  curr_seq.insert(curr_seq.end(), corner_vec, corner_vec+((s1==2)?2:4));
	}
      seqs.push_back(curr_seq);
      nmb_sub = 1;
    }
  else
    {
#ifdef DEBUG
      std::ofstream of("sub_cloud2.g2");
      (void)of.precision(15);
      of << "400 1 0 0" << std::endl;
      of << nmb_pts << std::endl;
      for (int kr=0; kr<nmb_pts; ++kr)
	{
	  for (int kh=0; kh<del-1; ++kh)
	    of << points[kr*del+kh] << " ";
	  of << 0 << std::endl;
	}
#endif

      // Reduce the number of divisions if the number of points is small
      int nmb2 = (int)sqrt((double)nmb_pts);
      nmb_u = std::min(nmb_u, std::max(3, nmb2));
      nmb_v = std::min(nmb_v, std::max(3, nmb2));

      nmb_sub = 0;
      double u_del = (minmax[1] - minmax[0])/(double)(nmb_u);
      double v_del = (minmax[3] - minmax[2])/(double)(nmb_v);
      vector<vector<int> > end_ix;

      // Sort points in v-direction
      qsort(points, nmb_pts, del*sizeof(double), compare_v_par);
      
      // Distribute points into strips
      int ki, kj;
      int pp0, pp1, pp2, pp3;
      int ix0, ix1, ix2, ix3;
      double upar, vpar;
      for (ki=0, pp0=0, ix0=0, vpar=minmax[2]+v_del; ki<nmb_v; 
	   ++ki, vpar+=v_del, ix0=ix1, pp0=pp1)
	{
	  vector<int> ixs;
	  if (ki == nmb_v-1)
	    {
	      ix1 = nmb_pts;
	      pp1 = ix1*del;
	    }
	  else
	    {
	      for (pp1=pp0, ix1=ix0; ix1<nmb_pts && points[pp1+1]<vpar; 
		   pp1+=del, ++ix1);
	    }

	  // Sort according to the u-parameter
	  qsort(points+pp0, (pp1-pp0)/del, del*sizeof(double), compare_u_par);
	   
	  ixs.push_back(ix0);
	  for (kj=0, pp2=pp0, upar=minmax[0]+u_del, ix2=ix0; kj<nmb_u;
	       ++kj, upar+=u_del, ix2=ix3, pp2=pp3)
	    {
	      if (kj == nmb_v-1)
		{
		  pp3 = pp1;
		  ix3 = pp3/del;
		}
	      else
		{
		  for (pp3=pp2, ix3=ix2; pp3<pp1 && points[pp3]<upar; 
		       pp3+=del, ++ix3);
		}
	      ixs.push_back(ix3);
	    }
	  end_ix.push_back(ixs);
	}
    
	
      // Recursively compute trim parameters
      double domain[4];
      domain[2] = minmax[2];
      for (ki=0; ki<nmb_v; ++ki, domain[2]=domain[3])
	{
	  domain[3] = std::min(domain[2]+v_del, minmax[3]);

	  domain[0] = minmax[0];
	  int prev, curr, next;
	  int nmb_sub2 = 0;
	  int side2_prev = 0;
	  for (kj=0, prev=0, curr=end_ix[ki][1]-end_ix[ki][0];
	       kj<nmb_u; 
	       ++kj, domain[0]=domain[1], prev=curr, curr=next)
	    {
	      domain[1] = std::min(domain[0]+u_del, minmax[1]);

	      // Number of points in sub cloud
	      next = (kj == nmb_u-1) ? 0 : 
		end_ix[ki][kj+2] - end_ix[ki][kj+1];

	      if (curr == 0)
		continue;   // Not inside the trimming loop
	      nmb_sub++;

	      // Check if the sub cloud is to be divided further
	      int v1=0, v2=0;
	      for (int kr=0; kr<ki; ++kr)
		v1 += (end_ix[kr][kj+1]-end_ix[kr][kj]);
	      for (int kr=ki+1; kr<nmb_v; ++kr)
		v2 += (end_ix[kr][kj+1]-end_ix[kr][kj]);

	      int side2 = 0;
	      if ((s1 == 1 || s1 == 3 || kj > 0) && 
		  end_ix[ki][kj]-end_ix[ki][0]==0) //prev == 0)
		side2 += 1;
	      if ((s1 == 2 || s1 == 3 || kj < nmb_u-1) && 
		  end_ix[ki][nmb_u]-end_ix[ki][kj+1] == 0) //next == 0)
		side2 += 2;
	      if ((s2 == 10 || s2 == 30 || ki > 0) && v1 == 0)
		//(ki == 0 || end_ix[ki-1][kj+1]-end_ix[ki-1][kj] == 0))
		side2 += 10;
	      if ((s2 == 20 || s2 == 30 || ki < nmb_v-1) && v2 == 0)
		//(ki == nmb_v-1 || end_ix[ki+1][kj+1]-end_ix[ki+1][kj] == 0))
		side2 += 20;

	      if (side2 > 0 || 
		  (max_level > 1 && side2_prev > 0 && nmb_sub2 < (nmb_u-1)*(nmb_v-1)))
		{
		  int first_ix = end_ix[ki][kj]*del;
#ifdef DEBUG
		  std::ofstream ofa("sub_cloud2_curr.g2");
		  (void)ofa.precision(15);
		  ofa << "400 1 0 0" << std::endl;
		  ofa << curr << std::endl;
		  for (int kr=0; kr<curr; ++kr)
		    {
		      for (int kh=0; kh<del-1; ++kh)
			ofa << points[first_ix+kr*del+kh] << " ";
		      ofa << 0 << std::endl;
		    }
#endif
		  // Compute trim parameters related to sub cloud
		  vector<vector<double> > seqs2;
		  nmb_sub2 = 0;
		  // if (side2 == 0)
		  //   side2 = side2_prev;
		  computeTrimInfo2(points+first_ix, dim,
				   curr, domain, 
				   max_level-1, nmb_u, nmb_v,
				   side2, seqs2, nmb_sub2);
		  
#ifdef DEBUG
		  std::ofstream ofc("trim_seqs_curr.g2");
		  (void)ofc.precision(15);
		  for (int kr=0; kr<(int)seqs2.size(); ++kr)
		    {
		      ofc << "410 1 0 0" << std::endl;
		      ofc << seqs2[kr].size()/2-1 << std::endl;
		      for (int kh=0; kh<=(int)seqs2[kr].size()-4; kh+=2)
			{
			  ofc << seqs2[kr][kh] << " " << seqs2[kr][kh+1] << " " << 0 << " ";
			  ofc << seqs2[kr][kh+2] << " " << seqs2[kr][kh+3] << " " << 0 << std::endl;
			}
		    }
#endif		 

		  // Merge trim parameter info
		  mergeTrimSeqs(seqs, seqs2, side2);

		  // Check if the previous block should be included in
		  // the next recursion level
		  if (max_level > 1 && side2_prev == 0 && prev > 0 && 
		      nmb_sub2 < (nmb_u-1)*(nmb_v-1))
		    {
		      first_ix = end_ix[ki][kj-1]*del;
#ifdef DEBUG
		      std::ofstream ofb("sub_cloud2_curr2.g2");
		      (void)ofb.precision(15);
		      ofb << "400 1 0 0" << std::endl;
		      ofb << prev << std::endl;
		      for (int kr=0; kr<prev; ++kr)
			{
			  for (int kh=0; kh<del-1; ++kh)
			    ofb << points[first_ix+kr*del+kh] << " ";
			  ofb << 0 << std::endl;
			}
#endif
		  
		      seqs2.clear();
		      int nmb_sub3 = 0;
		      double domain3[4];
		      domain3[0] = domain[0] - u_del;
		      domain3[1] = domain[0];
		      domain3[2] = domain[2];
		      domain3[3] = domain[3];
		      computeTrimInfo2(points+first_ix, dim,
				       prev, domain3, 
				       max_level-1, nmb_u, nmb_v,
				       0/*side2*/, seqs2, nmb_sub3);
		      
#ifdef DEBUG
		  std::ofstream ofd("trim_seqs_curr.g2");
		  (void)ofd.precision(15);
		  for (int kr=0; kr<(int)seqs2.size(); ++kr)
		    {
		      ofd << "410 1 0 0" << std::endl;
		      ofd << seqs2[kr].size()/2-1 << std::endl;
		      for (int kh=0; kh<=(int)seqs2[kr].size()-4; kh+=2)
			{
			  ofd << seqs2[kr][kh] << " " << seqs2[kr][kh+1] << " " << 0 << " ";
			  ofd << seqs2[kr][kh+2] << " " << seqs2[kr][kh+3] << " " << 0 << std::endl;
			}
		    }
#endif		 
		      mergeTrimSeqs(seqs, seqs2, 0);
		    }
		}
	      side2_prev = side2;
	    }
	}
#ifdef DEBUG
  std::ofstream of2("trim_seqs_loc.g2");
  (void)of2.precision(15);
  for (ki=0; ki<(int)seqs.size(); ++ki)
    {
      of2 << "410 1 0 0" << std::endl;
      of2 << seqs[ki].size()/2-1 << std::endl;
      for (kj=0; kj<=(int)seqs[ki].size()-4; kj+=2)
	{
	  of2 << seqs[ki][kj] << " " << seqs[ki][kj+1] << " " << 0 << " ";
	  of2 << seqs[ki][kj+2] << " " << seqs[ki][kj+3] << " " << 0 << std::endl;
	}
    }
  int stop_break = 1;
#endif
    }
}
