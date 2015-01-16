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
#include <string.h>

using namespace Go;
using std::vector;
using std::string;

#define DEBUG


//==============================================================================
void LRSurfApproxUtils::readBlockMeta(std::ifstream& is,
				      vector<string>& inblock,
				      vector<int>& nmb_points,
				      vector<double>& domain)
//==============================================================================
{
<<<<<<< HEAD
  int del = dim+2;
  int nmb = (int)points.size()/del;  // Number of data points

  // Compute the domain corresponding to the point set
  int ki, kj;
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

  // Merge sequences with "exact" match
  double eps2 = 1.0e-12; // A small number
  double limit = std::min((minmax[1] - minmax[0])/(double)(nmb_u),
			  (minmax[3] - minmax[2])/(double)(nmb_v));
  limit = eps2;
  for (ki=0; ki<(int)seqs2.size(); )
=======
  char cc;
  is >> cc;
  int nmb;
  double xx;
  while (cc != ']')
>>>>>>> Tiling of data sets
    {
      is >> cc;
      while (cc == '{')
	{
    	  is >> cc;  // Expects "
    	  is >> cc;  // Expects f
    	  while (cc != '"')
    	    is >> cc;  // Expects ile"
    	  is >> cc;    // Expects :
    	  is >> cc;    // Expects "
    	  is >> cc;    // Expects first character in filename
	  int ki=0;
	  char filename[80];
    	  while (cc != '"')
	    {
	      filename[ki++] = cc;
	      is >> cc;   // Expects remaining file name + "
	    }
	  inblock.push_back(std::string(filename, filename+ki));
    	  while (cc != ':')
    	    is >> cc;                  // Expects "number of points"
	  is >> nmb;
	  nmb_points.push_back(nmb);   // Number of points in file
    	  is >> cc;
	  while (cc != '[')
	    is >> cc;            // Expects "domain" : [
	  for (int ki=0; ki<4; ++ki)
	    {
	      is >> xx;              // Parameter domain: xmin, xmax, ymin, ymax
	      domain.push_back(xx);
	      is >> cc;
	    }
	  is >> cc;  // Expects ]
	  is >> cc;  // Expects }
	  is >> cc;  // Expects , or ]
	}
    }
}

<<<<<<< HEAD
#ifdef DEBUG
  for (ki=0; ki<(int)seqs2.size(); ++ki)
    {
      of << "410 1 0 4 255 0 0 255" << std::endl;
      of << seqs2[ki].size()/2-1 << std::endl;
      for (kj=0; kj<=(int)seqs2[ki].size()-4; kj+=2)
	{
	  of << seqs2[ki][kj] << " " << seqs2[ki][kj+1] << " " << 0 << " ";
	  of << seqs2[ki][kj+2] << " " << seqs2[ki][kj+3] << " " << 0 << std::endl;
	}
=======
//==============================================================================
void LRSurfApproxUtils::writeBlockMeta(std::ofstream& os,
				       vector<int>& nmb_points,
				       vector<double>& domain,
				       vector<string>& file_name)
//==============================================================================
{
  // For each set of output entities, create filename based on the given
  // root and write filename and related information to the given output
  // stream. Return constructed filenames.

  // Check input
  if (nmb_points.size() != file_name.size() ||
      4*nmb_points.size() != domain.size())
    return;

  (void)os.precision(15);

  int ki, kj;
  int nmb_blocks = (int)nmb_points.size();
  os << "[" << std::endl;
  for (ki=0; ki<nmb_blocks; ++ki)
    {
      os << "{" << std::endl;
      os << "\"file\": ";
      os << "\"" << file_name[ki] << "\"," << std::endl;
      os << "\"number of points\": " << nmb_points[ki] <<"," << std::endl;
      os << "\"domain\": [" << domain[4*ki];
      for (kj=1; kj<4; ++kj)
	os << "," << domain[ki*4+kj];
      os << "]" << std::endl << "}";
      if (ki < nmb_blocks-1)
	os << ",";
      os << std::endl;
>>>>>>> Tiling of data sets
    }
  os << "]" << std::endl;
}

//==============================================================================
void LRSurfApproxUtils::readTileMeta(std::ifstream& is,
				     double total_domain[],
				     int& nmb_u, int& nmb_v,
				     double& u_overlap, double& v_overlap,
				     vector<string>& intile,
				     vector<int>& nmb_points,
				     vector<double>& domain)
//==============================================================================
{
  char cc;
  double xx;
  int nmb;
  is >> cc;  // Expects {
  is >> cc;  // Expects "
  while (cc != '{')
    is >> cc;  // Expects Meta":
  is >> cc;
  while (cc != '[')
    is >> cc;            // Expects "Total domain" : [
  for (int ki=0; ki<4; ++ki)
    {
      is >> xx;              // Parameter domain: xmin, xmax, ymin, ymax
      total_domain[ki] = xx;
      is >> cc;
    }
  while (cc != ':')
    is >> cc;                  // Expects "Nmb u"
  is >> nmb_u;
  is >> cc;
  while (cc != ':')
    is >> cc;                  // Expects "Nmb v"
  is >> nmb_v;
  is >> cc;
  while (cc != ':')
    is >> cc;                  // Expects "Overlap u"
  is >> u_overlap;
  is >> cc;
  while (cc != ':')
    is >> cc;                  // Expects "Overlap v"
  is >> v_overlap;
  is >> cc;
  while (cc != '[')
    is >> cc;                  // Expects "}, "Detail": "

  while (cc != ']')
    {
      is >> cc;
      while (cc == '{')
	{
    	  is >> cc;  // Expects "
    	  is >> cc;  // Expects f
    	  while (cc != '"')
    	    is >> cc;  // Expects ile"
    	  is >> cc;    // Expects :
    	  is >> cc;    // Expects "
    	  is >> cc;    // Expects first character in filename
	  int ki=0;
	  char filename[80];
    	  while (cc != '"')
	    {
	      filename[ki++] = cc;
	      is >> cc;   // Expects remaining file name + "
	    }
	  intile.push_back(std::string(filename, filename+ki));
    	  while (cc != ':')
    	    is >> cc;                  // Expects "number of points"
	  is >> nmb;
	  nmb_points.push_back(nmb);   // Number of points in file
    	  is >> cc;
	  while (cc != '[')
	    is >> cc;            // Expects "domain" : [
	  for (int ki=0; ki<4; ++ki)
	    {
	      is >> xx;              // Parameter domain: xmin, xmax, ymin, ymax
	      domain.push_back(xx);
	      is >> cc;
	    }
<<<<<<< HEAD
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
  
  // Remove duplicate parameter points
  for (kj=0, ki=2; ki<(int)seq.size(); )
    {
      double dist2 = Utils::distance_squared(seq.begin()+kj,
					     seq.begin()+ki,
					     seq.begin()+ki);
      if (dist2 < eps2)
	seq.erase(seq.begin()+ki, seq.begin()+ki+2);
      else
	{
	  ki += 2;
	  kj += 2;
	}
    }
  
=======
	  is >> cc;  // Expects ]
	  is >> cc;  // Expects }
	  is >> cc;  // Expects , or ]
	}
    }
>>>>>>> Tiling of data sets
}

//==============================================================================
void LRSurfApproxUtils::readSurfMeta(std::ifstream& is,
				     double total_domain[],
				     int& nmb_u, int& nmb_v,
				     double& eps, int& max_iter,
				     vector<int>& nmb_points,
				     vector<double>& max_dists,
				     vector<double>& av_dists,
				     vector<int>& nmb_outside,
				     vector<string>& file_name)
//==============================================================================
{
  char cc;
  double xx;
  int nmb;
  double dist;
  is >> cc;  // Expects {
  is >> cc;  // Expects "
  while (cc != '{')
    is >> cc;  // Expects Meta":
  is >> cc;
  while (cc != '[')
    is >> cc;            // Expects "Total domain" : [
  for (int ki=0; ki<4; ++ki)
    {
      is >> xx;              // Parameter domain: xmin, xmax, ymin, ymax
      total_domain[ki] = xx;
      is >> cc;
    }
  while (cc != ':')
    is >> cc;                  // Expects "Nmb u"
  is >> nmb_u;
  is >> cc;
  while (cc != ':')
    is >> cc;                  // Expects "Nmb v"
  is >> nmb_v;
  is >> cc;
  while (cc != ':')
    is >> cc;                  // Expects "Tolerance"
  is >> eps;
  is >> cc;
  while (cc != ':')
    is >> cc;                  // Expects "Maximum number of iterations"
  is >> max_iter;
  is >> cc;
  while (cc != '[')
    is >> cc;                  // Expects "}, "Detail": "

  while (cc != ']')
    {
      is >> cc;
      while (cc == '{')
	{
    	  is >> cc;  // Expects "
    	  is >> cc;  // Expects F
    	  while (cc != '"')
    	    is >> cc;  // Expects ile"
    	  is >> cc;    // Expects :
    	  is >> cc;    // Expects "
    	  is >> cc;    // Expects first character in filename
	  int ki=0;
	  char filename[80];
    	  while (cc != '"')
	    {
	      filename[ki++] = cc;
	      is >> cc;   // Expects remaining file name + "
	    }
	  file_name.push_back(std::string(filename, filename+ki));
    	  while (cc != ':')
    	    is >> cc;                  // Expects "Number of points"
	  is >> nmb;
	  nmb_points.push_back(nmb);   // Number of points in file
    	  is >> cc;
    	  while (cc != ':')
    	    is >> cc;                  // Expects "Maximum distance
	  is >> dist;
	  max_dists.push_back(dist);
    	  is >> cc;
    	  while (cc != ':')
    	    is >> cc;                  // Expects "Average distance
	  is >> dist;
	  av_dists.push_back(dist);
    	  is >> cc;
    	  while (cc != ':')
    	    is >> cc;                  // Expects "Number of points outside tolerance"
	  is >> nmb;
	  nmb_outside.push_back(nmb);
	  is >> cc;  // Expects }
	  is >> cc;  // Expects , or ]
	}
    }
}

//==============================================================================
void LRSurfApproxUtils::writeTileMeta(std::ofstream& os,
				      double total_domain[],
				      int nmb_u, int nmb_v,
				      double u_overlap, double v_overlap,
				      vector<int>& nmb_points,
				      vector<double>& domain,
				      vector<string>& file_name)
//==============================================================================
{
  // For each set of output entities, write filename and related 
  // information to the given output stream.

  // Check input
  int nmb_tiles = (int)nmb_points.size();
  if (nmb_tiles != nmb_u*nmb_v)
    return;
  if (nmb_tiles != (int)file_name.size() ||
      4*nmb_tiles != (int)domain.size())
    return;

  (void)os.precision(15);

  int ki, kj;
  os << "{" << std::endl;
  os << "\"Meta\": {" << std::endl;
  os << "\"Total domain\": [" << total_domain[0];
  for (kj=1; kj<4; ++kj)
    os << "," << total_domain[kj];
  os << "]," << std::endl;
  os <<"\"Nmb u\": " << nmb_u << "," << std::endl;
  os <<"\"Nmb v\": " << nmb_v << "," << std::endl;
  os <<"\"Overlap u\": " << u_overlap << "," << std::endl;
  os <<"\"Overlap v\": " << v_overlap << std::endl;
  os <<"}," << std::endl;
  os << "\"Detail\": [" << std::endl;
  for (ki=0; ki<nmb_tiles; ++ki)
    {
      os << "{" << std::endl;
      os << "\"file\": ";
      os << "\"" << file_name[ki] << "\"," << std::endl;
      os << "\"number of points\": " << nmb_points[ki] <<"," << std::endl;
      os << "\"domain\": [" << domain[4*ki];
      for (kj=1; kj<4; ++kj)
	os << "," << domain[ki*4+kj];
      os << "]" << std::endl << "}";
      if (ki < nmb_tiles-1)
	os << ",";
      os << std::endl;
    }
<<<<<<< HEAD
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
=======
  os << "]" << std::endl;
}
>>>>>>> Tiling of data sets

//==============================================================================
void LRSurfApproxUtils::writeSurfMeta(std::ofstream& os,
				      double total_domain[],
				      int nmb_u, int nmb_v,
				      double eps, int max_iter,
				      vector<int>& nmb_points,
				      vector<double>& max_dists,
				      vector<double>& av_dists,
				      vector<int>& nmb_outside,
				      vector<string>& file_name)
//==============================================================================
{
  // Check input
  int nmb_tiles = (int)nmb_points.size();
  if (nmb_tiles != (int)file_name.size() ||
      nmb_tiles != (int)max_dists.size() ||
      nmb_tiles != (int)av_dists.size() ||
      nmb_tiles != (int)nmb_outside.size())
    return;

  (void)os.precision(15);

  int ki, kj;
  os << "{" << std::endl;
  os << "\"Meta\": {" << std::endl;
  os << "\"Total domain\": [" << total_domain[0];
  for (kj=1; kj<4; ++kj)
    os << "," << total_domain[kj];
  os << "]," << std::endl;
  os <<"\"Nmb u\": " << nmb_u << "," << std::endl;
  os <<"\"Nmb v\": " << nmb_v << "," << std::endl;
  os <<"\"Tolerance\": " << eps << "," << std::endl;
  os <<"\"Maximum number of iterations\": " << max_iter << "," << std::endl;
  os <<"}," << std::endl;
  os << "\"Detail\": [" << std::endl;
  for (ki=0; ki<nmb_tiles; ++ki)
    {
      os << "{" << std::endl;
      os << "\"File\": ";
      os << "\"" << file_name[ki] << "\"," << std::endl;
      os << "\"Number of points\": " << nmb_points[ki] <<"," << std::endl;
      os << "\"Maximum distance\": " << max_dists[ki] <<"," << std::endl;
      os << "\"Average distance\": " << av_dists[ki] <<"," << std::endl;
      os << "\"Number of points outside tolerance\": " << nmb_outside[ki] << std::endl;
      os << std::endl << "}";
      if (ki < nmb_tiles-1)
	os << ",";
      os << std::endl;
    }
  os << "]" << std::endl;
}

 //==============================================================================
void LRSurfApproxUtils::fetchFileNames(const char* file_root,
				       int extension_type,
				       int nmb_files,
				       vector<string>& file_name)
//==============================================================================
{
  // For each set of output entities, create filename based on the given
  // root 
  for (int ki=0; ki<nmb_files; ++ki)
    {
      char outfile[90];
      strcpy(outfile, file_root);
      char tmp[5];
      sprintf(tmp, "_%d", ki+1);
      strncat(outfile, tmp, 4);
      if (extension_type == 1)
	strncat(outfile, ".g2", 3);
      else
	strncat(outfile, ".txt", 4);  // For the time being, plans also las
 	
      file_name.push_back(std::string(outfile));
    }
}
