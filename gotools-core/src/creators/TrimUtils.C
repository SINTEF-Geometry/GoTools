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


#include "GoTools/creators/TrimUtils.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/GoIntersections.h"
#include "GoTools/geometry/SplineCurve.h"
#include <iostream>
#include <fstream>

using namespace Go;
using std::vector;
using std::pair;

#define DEBUG

static int compare_u_par(const void* el1, const void* el2)
{
  if (((double*)el1)[0] < ((double*)el2)[0])
    return -1;
  else if (((double*)el1)[0] > ((double*)el2)[0])
    return 1;
  else
    return 0;
}

static int compare_v_par(const void* el1, const void* el2)
{
  if (((double*)el1)[1] < ((double*)el2)[1])
    return -1;
  else if (((double*)el1)[1] > ((double*)el2)[1])
    return 1;
  else
    return 0;
}

//==============================================================================
TrimUtils::TrimUtils(double *points, int nmb_pts, int dim)
  : eps2_(1.0e-12), points_(points), nmb_pts_(nmb_pts), dim_(dim), 
    del_u_(0.0), del_v_(0.0)
//==============================================================================
{
  // Compute the domain corresponding to the point set
  int ki, kj;
  int del = dim_ + 2;
  for (ki=0; ki<2; ++ki)
    domain_[2*ki] = domain_[2*ki+1] = points_[ki];
  for (kj=1; kj<nmb_pts_; ++kj)
    for (ki=0; ki<2; ++ki)
      {
	domain_[2*ki] = std::min(domain_[2*ki], points_[kj*del+ki]);
	domain_[2*ki+1] = std::max(domain_[2*ki+1], points_[kj*del+ki]);
      }
}

//==============================================================================
TrimUtils::TrimUtils(double *points, int nmb_pts, int dim, double domain[])
  : eps2_(1.0e-12), points_(points), nmb_pts_(nmb_pts), dim_(dim), 
    del_u_(0.0), del_v_(0.0)
//==============================================================================
{
  // Set domain. The domain is assumed to be equal to or larger than
  // the domain corresponding to the point set
  for (int ki=0; ki<4; ++ki)
    domain_[ki] = domain[ki];
}

//==============================================================================
TrimUtils::~TrimUtils()
//==============================================================================
{
}

//==============================================================================
void TrimUtils::computeTrimSeqs(int max_level, int nmb_div,
				  vector<vector<double> >& seqs,
				  bool outer_only)
//==============================================================================
{ 
  int del = dim_ + 2;

  SubCloud cloud;
  cloud.setInfo(nmb_pts_, 0, nmb_pts_*del, domain_, domain_);
  computeTrimInfo(cloud, max_level, nmb_div, seqs);

  // Remove fragments
  int limitsize = 4;
  for (int ki=1; ki<max_level; ++ki)
    limitsize += 4;
  limitsize *= 2;
  cleanTrimResults(limitsize, seqs);

  // Organize sequences and possibly remove inner trimming sequences
  if (seqs.size() > 0)
    reOrganizeSeqs(seqs, outer_only);

}

//==============================================================================
void TrimUtils::computeAllTrimSeqs(int max_level, int nmb_div,
				   vector<vector<vector<double> > >& allseqs)
//==============================================================================
{ 
  int del = dim_ + 2;

  SubCloud cloud;
  cloud.setInfo(nmb_pts_, 0, nmb_pts_*del, domain_, domain_);
  vector<vector<double> > seqs;
  computeTrimInfo(cloud, max_level, nmb_div, seqs);

  // Remove fragments
  int limitsize = 4;
  for (int ki=1; ki<max_level; ++ki)
    limitsize += 4;
  limitsize *= 2;
  cleanTrimResults(limitsize, seqs);

  // Organize sequences and possibly remove inner trimming sequences
  if (seqs.size() > 0)
    {
      vector<int> outer(seqs.size(), 0);
      reOrganizeSeqs(seqs, outer);
      size_t ki, kj=0;
      for (ki=0; ki<seqs.size(); ki=kj)
	{
	  for (kj=ki+1; kj<seqs.size() && outer[kj]==0; ++kj);
	  vector<vector<double> > loop(seqs.begin()+ki, seqs.begin()+kj);
	  allseqs.push_back(loop);
	}
    }
  else
    allseqs.push_back(seqs);

}

//==============================================================================
void TrimUtils::cleanTrimResults(int limitsize,
				   vector<vector<double> >& seqs)
//==============================================================================
{
  // Merge sequences with exact match
  int ki, kj, kh;
  for (ki=0; ki<(int)seqs.size(); )
    {
      for (kj=ki+1; kj<(int)seqs.size(); ++kj)
	{
	  double dist2 = Utils::distance_squared(seqs[kj].begin(),
						 seqs[kj].begin()+2,
						 seqs[ki].begin()+seqs[ki].size()-2);
	  if (dist2 < eps2_)
	    {
	      seqs[ki].insert(seqs[ki].end(), seqs[kj].begin(), seqs[kj].end());
	      seqs.erase(seqs.begin()+kj);
	      break;
	    }

	  dist2 = Utils::distance_squared(seqs[ki].begin(),
					  seqs[ki].begin()+2,
					  seqs[kj].begin()+seqs[kj].size()-2);
	  if (dist2 < eps2_)
	    {
	      seqs[kj].insert(seqs[kj].end(), seqs[ki].begin(), seqs[ki].end());
	      std::swap(seqs[ki], seqs[kj]);
	      seqs.erase(seqs.begin()+kj);
	      break;
	    }
	}
      if (kj == (int)seqs.size())
	++ki;
    }

#ifdef DEBUG
  std::ofstream of1("trim_seq1.g2");
  for (ki=0; ki<seqs.size(); ++ki)
    {
      of1 << "410 1 0 0" << std::endl;
      of1 << seqs[ki].size()/2 - 1 << std::endl;
      for (kj=2; kj<seqs[ki].size(); kj+=2)
	{
	  of1 << seqs[ki][kj-2] << " " << seqs[ki][kj-1] << " " << 0 << " ";
	  of1 << seqs[ki][kj] << " " << seqs[ki][kj+1] << " " << 0 << std::endl;
	}
    }
#endif

  // Remove fragments
  for (ki=0; ki<(int)seqs.size(); )
    {
      // double dist2 = 
      // 	Utils::distance_squared(seqs[ki].begin(),
      // 				seqs[ki].begin()+2,
      // 				seqs[ki].begin()+seqs[ki].size()-2);
      if (/*dist2 <eps2_ && */(int)seqs[ki].size() <= limitsize)
	seqs.erase(seqs.begin()+ki);
      else
	++ki;
    }
  
#ifdef DEBUG
  std::ofstream of2("trim_seq2.g2");
  for (ki=0; ki<seqs.size(); ++ki)
    {
      of2 << "410 1 0 0" << std::endl;
      of2 << seqs[ki].size()/2 - 1 << std::endl;
      for (kj=2; kj<seqs[ki].size(); kj+=2)
	{
	  of2 << seqs[ki][kj-2] << " " << seqs[ki][kj-1] << " " << 0 << " ";
	  of2 << seqs[ki][kj] << " " << seqs[ki][kj+1] << " " << 0 << std::endl;
	}
    }
#endif

  // Merge segments with closest match
  for (ki=0; ki<(int)seqs.size()-1; )
    {
      double min1 = std::numeric_limits<double>::max();
      double min2 = std::numeric_limits<double>::max();
      int ix1 = -1, ix2 = -1;
      for (kj=ki+1; kj<(int)seqs.size(); ++kj)
	{
	  double dist0 = /*(seqs[kj].size() < limitsize) ? eps2_ :*/
	    Utils::distance_squared(seqs[kj].begin(), seqs[kj].begin()+2,
				    seqs[kj].begin()+seqs[kj].size()-2);

	  double dist2 = Utils::distance_squared(seqs[kj].begin(),
						 seqs[kj].begin()+2,
						 seqs[ki].begin()+seqs[ki].size()-2);
	  double dist3 = Utils::distance_squared(seqs[ki].begin(),
					  seqs[ki].begin()+2,
					  seqs[kj].begin()+seqs[kj].size()-2);

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
	  seqs[ki].insert(seqs[ki].end(), seqs[ix1].begin(), seqs[ix1].end());
	  seqs.erase(seqs.begin()+ix1);
	}
      else if (ix2 >= 0)
	{
	  seqs[ix2].insert(seqs[ix2].end(), seqs[ki].begin(), seqs[ki].end());
	  std::swap(seqs[ki], seqs[ix2]);
	  seqs.erase(seqs.begin()+ix2);
	}
      else
	++ki;
    }

  // Remove duplicate parameter points
  for (kh=0; kh<(int)seqs.size(); ++kh)
    {
      for (kj=0, ki=2; ki<(int)seqs[kh].size(); )
	{
	  double dist2 = Utils::distance_squared(seqs[kh].begin()+kj,
						 seqs[kh].begin()+ki,
						 seqs[kh].begin()+ki);
	  if (dist2 < eps2_)
	    seqs[kh].erase(seqs[kh].begin()+ki, seqs[kh].begin()+ki+2);
	  else
	    {
	      ki += 2;
	      kj += 2;
	    }
	}
    }

#ifdef DEBUG
  std::ofstream of3("trim_seq3.g2");
  for (ki=0; ki<seqs.size(); ++ki)
    {
      of3 << "410 1 0 0" << std::endl;
      of3 << seqs[ki].size()/2 - 1 << std::endl;
      for (kj=2; kj<seqs[ki].size(); kj+=2)
	{
	  of3 << seqs[ki][kj-2] << " " << seqs[ki][kj-1] << " " << 0 << " ";
	  of3 << seqs[ki][kj] << " " << seqs[ki][kj+1] << " " << 0 << std::endl;
	}
    }
#endif

  // Remove loops in the sequences
  for (kh=0; kh<(int)seqs.size(); ++kh)
    {
      for (ki=2; ki<(int)seqs[kh].size()-2; ki+=2)
	{
	  for (kj=ki+2; kj<(int)seqs[kh].size()-2; kj+=2)
	    {
	      double dist2 = Utils::distance_squared(seqs[kh].begin()+ki,
						     seqs[kh].begin()+ki+2,
						     seqs[kh].begin()+kj);
	      if (dist2 < eps2_)
		{
		  // A loop is found. 
		  // Distinguish between small and large loops
		  if (kj - ki <= limitsize)
		    {
		      // Remove loop
		      seqs[kh].erase(seqs[kh].begin()+ki+2,
				     seqs[kh].begin()+kj);
		    }
		  else
		    {
		      // Split loop
		      vector<double> curr_seq(seqs[kh].begin()+ki,
					      seqs[kh].begin()+kj+2);
		      seqs[kh].erase(seqs[kh].begin()+ki+2,
				     seqs[kh].begin()+kj);
		      seqs.push_back(curr_seq);
		    }
		  break;
		}
	    }
	  // if (kj >= (int)seqs[kh].size())
	  //   ki += 2;
	}
    }

#ifdef DEBUG
  std::ofstream of4("trim_seq4.g2");
  for (ki=0; ki<seqs.size(); ++ki)
    {
      of4 << "410 1 0 0" << std::endl;
      of4 << seqs[ki].size()/2 - 1 << std::endl;
      for (kj=2; kj<seqs[ki].size(); kj+=2)
	{
	  of4 << seqs[ki][kj-2] << " " << seqs[ki][kj-1] << " " << 0 << " ";
	  of4 << seqs[ki][kj] << " " << seqs[ki][kj+1] << " " << 0 << std::endl;
	}
    }
#endif

  // Remove fragments
  for (ki=0; ki<(int)seqs.size(); )
    {
      // double dist2 = 
      // 	Utils::distance_squared(seqs[ki].begin(),
      // 				seqs[ki].begin()+2,
      // 				seqs[ki].begin()+seqs[ki].size()-2);
      if (/*dist2 <eps2_ && */(int)seqs[ki].size() <= limitsize)
	seqs.erase(seqs.begin()+ki);
      else
	++ki;
    }
  
  // Close loops
  for (ki=0; ki<seqs.size(); ++ki)
    {
      double dist2 = 
      	Utils::distance_squared(seqs[ki].begin(),
      				seqs[ki].begin()+2,
      				seqs[ki].begin()+seqs[ki].size()-2);
      if (dist2 > eps2_)
	seqs[ki].insert(seqs[ki].end(), seqs[ki].begin(), seqs[ki].begin()+2);
    }

 #ifdef DEBUG
  std::ofstream of5("trim_seq5.g2");
  for (ki=0; ki<seqs.size(); ++ki)
    {
      of5 << "410 1 0 0" << std::endl;
      of5 << seqs[ki].size()/2 - 1 << std::endl;
      for (kj=2; kj<seqs[ki].size(); kj+=2)
	{
	  of5 << seqs[ki][kj-2] << " " << seqs[ki][kj-1] << " " << 0 << " ";
	  of5 << seqs[ki][kj] << " " << seqs[ki][kj+1] << " " << 0 << std::endl;
	}
    }
#endif
}

//==============================================================================
void TrimUtils::distributePointCloud(int ix1, int ix2,
				     double domain[4],
				     int nmb_u, int nmb_v,
				     vector<SubCloud>& sub_clouds)
//==============================================================================
{
  sub_clouds.resize(nmb_u*nmb_v);
  int del = dim_ + 2;

  double u_del = (domain[1] - domain[0])/(double)(nmb_u);
  double v_del = (domain[3] - domain[2])/(double)(nmb_v);

  // Sort points in v-direction
  int nmb_pts = (ix2 - ix1)/del;
  double *points = points_+ix1;
  qsort(points, nmb_pts, del*sizeof(double), compare_v_par);

#ifdef DEBUG
  std::ofstream of("division_lines.g2");
  of << "410 1 0 4 255 0 0 255" << std::endl;
  of << nmb_u+nmb_v+2 << std::endl;
  for (int ki=0; ki<=nmb_u; ++ki)
    {
      Point p1(domain[0]+ki*u_del, domain[2], 0.0);
      Point p2(domain[0]+ki*u_del, domain[3], 0.0);
      of << p1 << " " << p2 << std::endl;
    }
  for (int ki=0; ki<=nmb_v; ++ki)
    {
      Point p1(domain[0], domain[2]+ki*v_del, 0.0);
      Point p2(domain[1], domain[2]+ki*v_del, 0.0);
      of << p1 << " " << p2 << std::endl;
    }

#endif
  
  // Distribute points into strips
  int ki, kj, kr;
  int pp0, pp1, pp2, pp3;
  int ppmax = ix2 - ix1;
  double upar, vpar;
  double subdomain[4];
  for (kr=0, ki=0, pp0=0, vpar=domain[2]+v_del; ki<nmb_v; 
       ++ki, vpar+=v_del, pp0=pp1)
    {
      subdomain[2] = vpar - v_del;
      subdomain[3] = vpar;

      if (ki == nmb_v-1)
	{
	  pp1 = ppmax;
	}
      else
	{
	  for (pp1=pp0; pp1<ppmax && points[pp1+1]<vpar; pp1+=del);
	}

      // Sort according to the u-parameter
      qsort(points+pp0, (pp1-pp0)/del, del*sizeof(double), compare_u_par);
	   
      for (kj=0, pp2=pp0, upar=domain[0]+u_del; kj<nmb_u;
	   ++kj, upar+=u_del, pp2=pp3, ++kr)
	{
	  subdomain[0] = upar - u_del;
	  subdomain[1] = upar;

	  double bb[4];
	  bb[0] = subdomain[1];
	  bb[1] = subdomain[0];
	  bb[2] = subdomain[3];
	  bb[3] = subdomain[2];
	  
	  for (pp3=pp2; pp3<pp1 && points[pp3]<upar; pp3+=del)
	    {
	      bb[0] = std::min(bb[0], points[pp3]);
	      bb[1] = std::max(bb[1], points[pp3]);
	      bb[2] = std::min(bb[2], points[pp3+1]);
	      bb[3] = std::max(bb[3], points[pp3+1]);
	    }
	  if (kj == nmb_u-1)
	    {
	      for (; pp3<pp1; pp3+=del)
		{
		  bb[0] = std::min(bb[0], points[pp3]);
		  bb[1] = std::max(bb[1], points[pp3]);
		  bb[2] = std::min(bb[2], points[pp3+1]);
		  bb[3] = std::max(bb[3], points[pp3+1]);
		}
	    }
	  sub_clouds[kr].setInfo((pp3-pp2)/del, pp2+ix1, pp3+ix1, 
				 subdomain, bb);
	}
    }
}

//==============================================================================
void TrimUtils::setSubSeq(SubCloud& cloud,
			    vector<double>& seq)
//==============================================================================
{ 
  seq.resize(10);
  seq[0] = cloud.dom_[0];
  seq[1] = cloud.dom_[3];
  seq[2] = cloud.dom_[0];
  seq[3] = cloud.dom_[2];
  seq[4] = cloud.dom_[1];
  seq[5] = cloud.dom_[2];
  seq[6] = cloud.dom_[1];
  seq[7] = cloud.dom_[3];
  seq[8] = cloud.dom_[0];
  seq[9] = cloud.dom_[3];

  del_u_ = std::max(del_u_, cloud.dom_[1]-cloud.dom_[0]);
  del_v_ = std::max(del_v_, cloud.dom_[3]-cloud.dom_[2]);
}

//==============================================================================
void TrimUtils::computeTrimInfo(SubCloud& cloud,
				int max_level, int nmb_div,
				vector<vector<double> >& seqs)
//==============================================================================
{ 
  double fac = 0.66;   // Reduction factor in the number of sub_clouds
  if (max_level == 0)
    {
      vector<double> curr_seq;
      setSubSeq(cloud, curr_seq);
      seqs.push_back(curr_seq);
      return;
    }

  // Compute number of sub clouds in each direction
  int nmb_u, nmb_v;
  if (cloud.dom_[3]-cloud.dom_[2] > cloud.dom_[1]-cloud.dom_[0])
    {
      nmb_u = nmb_div;
      nmb_v = 
	(int)(nmb_div*(cloud.dom_[3]-cloud.dom_[2])/(cloud.dom_[1]-cloud.dom_[0]));
    }
  else
    {
      nmb_v = nmb_div;
      nmb_u = 
	(int)(nmb_div*(cloud.dom_[1]-cloud.dom_[0])/(cloud.dom_[3]-cloud.dom_[2]));
    }
  nmb_div = (int)(fac*nmb_div);
  
  // Distribute current point cloud into sub clouds
  vector<SubCloud> sub_clouds;
  distributePointCloud(cloud.ix1_, cloud.ix2_, cloud.dom_,
		       nmb_u, nmb_v, sub_clouds);

  // Traverse the sub clouds and check whether they need to be processed
  // further
  int ki, kj, kr;
  double frac = 0.9;
  for (kr=0, kj=0; kj<nmb_v; ++kj)
    {
      for (ki=0; ki<nmb_u; ++ki, ++kr)
	{
	  if (sub_clouds[kr].nmb_pts_ == 0)
	    continue;  // No points in sub cloud. Not inside trimming loop
	  
#ifdef DEBUG
	  int del = dim_ + 2;
	  std::ofstream ofa("sub_cloud_curr.g2");
	  (void)ofa.precision(15);
	  ofa << "400 1 0 0" << std::endl;
	  ofa << sub_clouds[kr].nmb_pts_ << std::endl;
	  for (int kr1=0; kr1<sub_clouds[kr].nmb_pts_; ++kr1)
	    {
	      for (int kh1=0; kh1<del-1; ++kh1)
		ofa << points_[sub_clouds[kr].ix1_+kr1*del+kh1] << " ";
	      ofa << 0 << std::endl;
	    }
#endif

	  vector<double> limitseq;  // Point sequences towards non-processed
	                            // sub clouds

	  // For each neighbour, check if it exists or if it contains no 
	  // points or contains points only in a part of the domain
	  int nmb_nolimits = 0;
	  if (true /*sub_clouds[kr].limitedSupport(frac)*/)
	    nmb_nolimits = 4;
	  else
	    {
	      if (ki==0 || sub_clouds[kr-1].nmb_pts_ == 0 ||
		  /*sub_clouds[kr-1].processed_ == true*/ sub_clouds[kr-1].limitedSupport(frac))
		{
		  nmb_nolimits++;
		  if (ki > 0 && sub_clouds[kr-1].processed_ && sub_clouds[kr].limit_[0])
		    {
		      vector<double> bdseq;
		      sub_clouds[kr].leftBd(bdseq);
		      limitseq.insert(limitseq.end(), bdseq.begin(), bdseq.end());
		      sub_clouds[kr-1].limit_[1] = 1;
		    }
		}
	      else
		{
		  vector<double> bdseq;
		  sub_clouds[kr].leftBd(bdseq);
		  limitseq.insert(limitseq.end(), bdseq.begin(), bdseq.end());
		  sub_clouds[kr-1].limit_[1] = 1;
		}
	      if (ki==nmb_u-1 || sub_clouds[kr+1].nmb_pts_ == 0 ||
		  sub_clouds[kr+1].limitedSupport(frac))
		{
		  nmb_nolimits++;
		}
	      else
		{
		  vector<double> bdseq;
		  sub_clouds[kr].rightBd(bdseq);
		  limitseq.insert(limitseq.end(), bdseq.begin(), bdseq.end());
		  sub_clouds[kr+1].limit_[0] = 1;
		}
	      if (kj==0 || sub_clouds[kr-nmb_u].nmb_pts_ == 0 ||
		  /*sub_clouds[kr-nmb_u].processed_ == true*/ sub_clouds[kr-nmb_u].limitedSupport(frac))
		{
		  nmb_nolimits++;
		  if (kj > 0 && sub_clouds[kr-nmb_u].processed_ && 
		      sub_clouds[kr].limit_[2] == 1)
		    {
		      vector<double> bdseq;
		      sub_clouds[kr].lowerBd(bdseq);
		      limitseq.insert(limitseq.end(), bdseq.begin(), bdseq.end());
		      sub_clouds[kr-nmb_u].limit_[3] = 1;
		    }
		}
	      else
		{
		  vector<double> bdseq;
		  sub_clouds[kr].lowerBd(bdseq);
		  limitseq.insert(limitseq.end(), bdseq.begin(), bdseq.end());
		  sub_clouds[kr-nmb_u].limit_[3] = 1;
		}
	      if (kj==nmb_v-1 || sub_clouds[kr+nmb_u].nmb_pts_ == 0 ||
		  sub_clouds[kr+nmb_u].limitedSupport(frac))
		{
		  nmb_nolimits++;
		}
	      else
		{
		  vector<double> bdseq;
		  sub_clouds[kr].upperBd(bdseq);
		  limitseq.insert(limitseq.end(), bdseq.begin(), bdseq.end());
		  sub_clouds[kr+nmb_u].limit_[2] = 1;
		}
	    }
	      
	  if (nmb_nolimits == 0)
	      continue;  // Not a boundary cloud. Do not process further

	  // Recursively compute seqences of the trimming loops
	  vector<vector<double> > seqs2;
	  computeTrimInfo(sub_clouds[kr], max_level-1, nmb_div, seqs2);
	  sub_clouds[kr].processed_ = true;

#ifdef DEBUG
	  std::ofstream ofd0("trim_seqs_curr0.g2");
	  (void)ofd0.precision(15);
	  for (int kr1=0; kr1<(int)seqs2.size(); ++kr1)
	    {
	      ofd0 << "410 1 0 0" << std::endl;
	      ofd0 << seqs2[kr1].size()/2-1 << std::endl;
	      for (int kh1=0; kh1<=(int)seqs2[kr1].size()-4; kh1+=2)
		{
		  ofd0 << seqs2[kr1][kh1] << " " << seqs2[kr1][kh1+1] << " " << 0 << " ";
		  ofd0 << seqs2[kr1][kh1+2] << " " << seqs2[kr1][kh1+3] << " " << 0 << std::endl;
		}
	    }
	  for (int kr1=0; kr1<(int)limitseq.size(); kr1+=4)
	    {
	      ofd0 << "410 1 0 4 0 255 0 255 " << std::endl;
	      ofd0 << 1 << std::endl;
	       ofd0 << limitseq[kr1] << " " << limitseq[kr1+1] << " " << 0 << " ";
		  ofd0 << limitseq[kr1+2] << " " << limitseq[kr1+3] << " " << 0 << std::endl;
	    }
#endif		 

	  // Remove trim seqences adjacent to non processed sub clouds
	  if (limitseq.size() > 0)
	    removeFalseTrimSeqs(limitseq, seqs2);

#ifdef DEBUG
	  std::ofstream ofd("trim_seqs_curr.g2");
	  (void)ofd.precision(15);
	  for (int kr1=0; kr1<(int)seqs2.size(); ++kr1)
	    {
	      ofd << "410 1 0 0" << std::endl;
	      ofd << seqs2[kr1].size()/2-1 << std::endl;
	      for (int kh1=0; kh1<=(int)seqs2[kr1].size()-4; kh1+=2)
		{
		  ofd << seqs2[kr1][kh1] << " " << seqs2[kr1][kh1+1] << " " << 0 << " ";
		  ofd << seqs2[kr1][kh1+2] << " " << seqs2[kr1][kh1+3] << " " << 0 << std::endl;
		}
	    }
#endif		 
	  // Merge trimming sequences
	  mergeTrimSeqs(seqs, seqs2);
#ifdef DEBUG
	  std::ofstream off("trim_seqs.g2");
	  (void)off.precision(15);
	  for (int kr1=0; kr1<(int)seqs.size(); ++kr1)
	    {
	      off << "410 1 0 0" << std::endl;
	      off << seqs[kr1].size()/2-1 << std::endl;
	      for (int kh1=0; kh1<=(int)seqs[kr1].size()-4; kh1+=2)
		{
		  off << seqs[kr1][kh1] << " " << seqs[kr1][kh1+1] << " " << 0 << " ";
		  off << seqs[kr1][kh1+2] << " " << seqs[kr1][kh1+3] << " " << 0 << std::endl;
		}
	    }
	  int stop_break = 1;
#endif		 
	}
    }

  // Simplify output
  // Merge sequences with exact match
  for (ki=0; ki<(int)seqs.size(); )
    {
      for (kj=ki+1; kj<(int)seqs.size(); ++kj)
	{
	  double dist2 = Utils::distance_squared(seqs[kj].begin(),
						 seqs[kj].begin()+2,
						 seqs[ki].begin()+seqs[ki].size()-2);
	  if (dist2 < eps2_)
	    {
	      seqs[ki].insert(seqs[ki].end(), seqs[kj].begin()+2, 
			      seqs[kj].end());
	      seqs.erase(seqs.begin()+kj);
	      break;
	    }

	  dist2 = Utils::distance_squared(seqs[ki].begin(),
					  seqs[ki].begin()+2,
					  seqs[kj].begin()+seqs[kj].size()-2);
	  if (dist2 < eps2_)
	    {
	      seqs[kj].insert(seqs[kj].end(), seqs[ki].begin()+2, 
			      seqs[ki].end());
	      std::swap(seqs[ki], seqs[kj]);
	      seqs.erase(seqs.begin()+kj);
	      break;
	    }
	}
      if (kj == (int)seqs.size())
	++ki;
    }

}

//==============================================================================
void TrimUtils::removeFalseTrimSeqs(vector<double>& limitseqs,
				      vector<vector<double> >& seqs)
//==============================================================================
{ 
  if (limitseqs.size() == 0)
    return;  // No need to remove segments due to limit sequences towards
  // non-processed point cloudes

  double eps = 1.0e-12;  // A small number used as numeric tolerance
  size_t ki, kj, kr, kh;
  for (ki=0; ki<seqs.size(); ++ki)
    {
      for (kj=0; kj<seqs[ki].size(); )
	{
	  bool found = false;
	  bool merge = false;

	  // Compare parameter value with end point of limit sequence
	  for (kr=0; kr<limitseqs.size(); kr+=4)
	    {
	      double d1 = Utils::distance_squared(limitseqs.begin()+kr,
						  limitseqs.begin()+kr+2,
						  seqs[ki].begin()+kj);
	      double d2 = Utils::distance_squared(limitseqs.begin()+kr+2,
						  limitseqs.begin()+kr+4,
						  seqs[ki].begin()+kj);
	      Point v1(limitseqs[kr+2]-limitseqs[kr],
		       limitseqs[kr+3]-limitseqs[kr+1]);
	      Point v2(seqs[ki][kj]-limitseqs[kr],
		       seqs[ki][kj+1]-limitseqs[kr+1]);
	      v1.normalize();
	      double dd = (v2 - (v1*v2)*v1).length2();
	      if (d1 <= eps || d2 <= eps || dd <= eps)
		{
		  // A match is found. Compute the extent
		  double d3, dd2;
		  size_t ix = (d2 < eps) ? kr : kr+2;

		  for (kh=kj+2; kh<seqs[ki].size(); kh+=2)
		    {
		      d3 = Utils::distance_squared(limitseqs.begin()+ix,
						   limitseqs.begin()+ix+2,
						   seqs[ki].begin()+kh);
		      Point vec1(limitseqs[kr+2]-limitseqs[kr],
				 limitseqs[kr+3]-limitseqs[kr+1]);
		      Point vec2(seqs[ki][kh]-limitseqs[kr],
				 seqs[ki][kh+1]-limitseqs[kr+1]);
		      vec1.normalize();
		      dd2 = (vec2 - (vec1*vec2)*vec1).length2();
		      if (d3 <= eps)
			break;
		      if (dd2 > eps)
			break;
		    }

		  if (d3 > eps || kh==seqs[ki].size())
		    {
		      // Not a complete match. Modify extent
		      kh -= 2;
		    }

		  if (kh > kj)
		    {
		      // A match is found. Remove piece and split
		      // the trim seqence if necessary
		      if (d3 <= eps && d1 <= eps)
			limitseqs.erase(limitseqs.begin()+kr,
					limitseqs.begin()+kr+4);
		      else
			{
			  // Partial match. Modify extent of limit seqence
			  if (d1 <= eps)
			    {
			      limitseqs[kr] = seqs[ki][kh];
			      limitseqs[kr+1] = seqs[ki][kh+1];
			    }
			  else if (d3 <= eps)
			    {
			      limitseqs[kr+2] = seqs[ki][kj];
			      limitseqs[kr+3] = seqs[ki][kj+1];
			    }
			  else
			    {
			      limitseqs.push_back(seqs[ki][kh]);
			      limitseqs.push_back(seqs[ki][kh+1]);
			      limitseqs.push_back(limitseqs[kr+2]);
			      limitseqs.push_back(limitseqs[kr+3]);
			      limitseqs[kr+2] = seqs[ki][kj];
			      limitseqs[kr+3] = seqs[ki][kj+1];
			    }
			}

		      if (kj == 0)
			seqs[ki].erase(seqs[ki].begin(), 
				       seqs[ki].begin()+kh);
		      else if (kh == seqs[ki].size() - 2)
			seqs[ki].erase(seqs[ki].begin()+kj+2, seqs[ki].end());
		      else
			{
			  vector<double> nseq(seqs[ki].begin()+kh,
					      seqs[ki].end());
			  seqs[ki].erase(seqs[ki].begin()+kj+2, 
					 seqs[ki].end());
			  dd = Utils::distance_squared(nseq.end()-2,
						       nseq.end(),
						       seqs[ki].begin());
			  if (dd < eps)
			    {
			      nseq.insert(nseq.end(), seqs[ki].begin()+2, 
					  seqs[ki].end());
			      std::swap(nseq, seqs[ki]);
			      merge = true;
			    }
			  else
			    seqs.push_back(nseq);
			}

		      found = true;
		      break;  // loop: for (kr=0; kr<limitseqs.size(); ...
		    }
		}
	    }
	  if (!found)
	    kj += 2; // No match found. Check next parameter pair
	  else if (merge)
	    kj = 0;  // Reiterate sequence
	}
    }
  for (ki=0; ki<seqs.size(); )
    {
      if (seqs[ki].size() == 0)
	seqs.erase(seqs.begin()+ki);
      else
	++ki;
    }
}


//==============================================================================
void TrimUtils::mergeTrimSeqs(vector<vector<double> >& seqs,
				vector<vector<double> >& seqs2)
//==============================================================================
{ 
  double eps = 1.0e-12; // A small number used as numeric tolerance

  int ki, kj, kr, kh;
  int kv, kw, ks, kt;

  // Remove double segments
  for (ki=0; ki<(int)seqs2.size(); ++ki)
    {
      for (kj=0; kj<(int)seqs2[ki].size(); kj+=2)
	{
	  for (kr=0; kr<(int)seqs.size(); ++kr)
	    {
	      for (kh=0; kh<(int)seqs[kr].size(); kh+=2)
		{
		  double dd = 
		    Utils::distance_squared(seqs[kr].begin()+kh,
					    seqs[kr].begin()+kh+2,
					    seqs2[ki].begin()+kj);
		  if (dd < eps)
		    {
		      // Equality in one point. Search for an interval
		      // equality. Adjacent sub group boundaries will
		      // have opposite direction
		      // Forward for seqs2
		      for (kv=kj+2, kw=kh-2; 
			   kv<(int)seqs2[ki].size() && kw>=0;
			   kv+=2, kw-=2)
			{
			  double dd2 =
			    Utils::distance_squared(seqs[kr].begin()+kw,
						    seqs[kr].begin()+kw+2,
						    seqs2[ki].begin()+kv);
			  if (dd2 >= eps)
			    break;
			}

		      // Backwards for seqs2
		      for (ks=kj-2, kt=kh+2;
			   ks>=0 && kt<(int)seqs[kr].size();
			   ks-=2, kt+=2)
			{
			  double dd2 =
			    Utils::distance_squared(seqs[kr].begin()+kt,
						    seqs[kr].begin()+kt+2,
						    seqs2[ki].begin()+ks);
			  if (dd2 >= eps)
			    break;
			}
		      
		      if (kv - ks > 4)
			{
			  // Remove double segment
			  if (ks < 0)
			    seqs2[ki].erase(seqs2[ki].begin(), 
					    seqs2[ki].begin()+kv-2);
			  else if (kv > (int)seqs2[ki].size())
			    seqs2[ki].erase(seqs2[ki].begin()+ks+4,
					    seqs2[ki].end());
			  else
			    {
			      vector<double> nseq(seqs2[ki].begin()+kv-2,
						  seqs2[ki].end());
			      seqs2.push_back(nseq);
			      seqs2[ki].erase(seqs2[ki].begin()+ks+4, 
					      seqs2[ki].end());
			    }
			}
		      if (kt - kw > 4)
			{
			  // Remove double segment
			  if (kw < 0)
			    seqs[kr].erase(seqs[kr].begin(), 
					   seqs[kr].begin()+kt-2);
			  else if (kt > (int)seqs[kr].size())
			    seqs[kr].erase(seqs[kr].begin()+kw+4,
					    seqs[kr].end());
			  else
			    {
			      vector<double> nseq(seqs[kr].begin()+kt-2,
						  seqs[kr].end());
			      seqs.push_back(nseq);
			      seqs[kr].erase(seqs[kr].begin()+kw+4, 
					     seqs[kr].end());
			    }
			}
		    }
		}
	    }
	}
    }

  // Clean up
  for (ki=0; ki<(int)seqs.size(); )
    {
      if (seqs[ki].size() <= 2) //== 0)
	seqs.erase(seqs.begin()+ki);
      else
	++ki;
    }

  for (ki=0; ki<(int)seqs2.size(); )
    {
      if (seqs2[ki].size() <= 2) //== 0)
	seqs2.erase(seqs2.begin()+ki);
      else
	++ki;
    }

#ifdef DEBUG
	  std::ofstream of1("trim_seqs_merge.g2");
	  (void)of1.precision(15);
	  for (int kr1=0; kr1<(int)seqs.size(); ++kr1)
	    {
	      of1 << "410 1 0 4 255 0 0 255" << std::endl;
	      of1 << seqs[kr1].size()/2-1 << std::endl;
	      for (int kh1=0; kh1<=(int)seqs[kr1].size()-4; kh1+=2)
		{
		  of1 << seqs[kr1][kh1] << " " << seqs[kr1][kh1+1] << " " << 0 << " ";
		  of1 << seqs[kr1][kh1+2] << " " << seqs[kr1][kh1+3] << " " << 0 << std::endl;
		}
	    }

	  for (int kr1=0; kr1<(int)seqs2.size(); ++kr1)
	    {
	      of1 << "410 1 0 4 0 255 0 255" << std::endl;
	      of1 << seqs2[kr1].size()/2-1 << std::endl;
	      for (int kh1=0; kh1<=(int)seqs2[kr1].size()-4; kh1+=2)
		{
		  of1 << seqs2[kr1][kh1] << " " << seqs2[kr1][kh1+1] << " " << 0 << " ";
		  of1 << seqs2[kr1][kh1+2] << " " << seqs2[kr1][kh1+3] << " " << 0 << std::endl;
		}
	    }
#endif

  // Join sequences
  if (seqs.size() == 0)
    seqs = seqs2;
  else
    {
      for (ki=0; ki<(int)seqs2.size(); ++ki)
	{

	  for (kr=0; kr<(int)seqs.size(); ++kr)
	    {
	      double dd1 =
		Utils::distance_squared(seqs[kr].begin(),
					seqs[kr].begin()+2,
					seqs2[ki].begin()+seqs2[ki].size()-2);
	      double dd2 =
		Utils::distance_squared(seqs2[ki].begin(),
					seqs2[ki].begin()+2,
					seqs[kr].begin()+seqs[kr].size()-2);
	      if (dd1 < eps)
		{
		  seqs2[ki].insert(seqs2[ki].end(), seqs[kr].begin()+2,
				   seqs[kr].end());
		  std::swap(seqs[kr], seqs2[ki]);
		  break;
		}
	      else if (dd2 < eps)
		{
		  seqs[kr].insert(seqs[kr].end(), seqs2[ki].begin()+2,
				  seqs2[ki].end());
		  break;
		}
	    }
	  if (kr == (int)seqs.size())
	    seqs.push_back(seqs2[ki]);
	}
    } 

}
#if 0
//==============================================================================
void TrimUtils::reOrganizeSeqs(vector<vector<double> >& seqs, bool outer_only)
//==============================================================================
{ 
  double tol = 1.0e-6;  // Intersection tolerance

  // Assume that the seqence with most elements is the outer one
  size_t ki, kj, kh, kr;
  for (ki=0; ki<seqs.size(); ++ki)
    for (kj=ki+1; kj<seqs.size(); ++kj)
      if (seqs[ki].size() < seqs[kj].size())
	std::swap(seqs[ki], seqs[kj]);

  if (outer_only)
    {
      // Close largest loop
      seqs[0].insert(seqs[0].end(), seqs[0].begin(), seqs[0].begin()+2);

      // Check if the remaining sequences lie internally to the first one
      for (ki=1; ki<seqs.size(); )
	{
	  // Select an arbitrary point and make a horizontal line through
	  // this point
	  Point p1(seqs[ki][0], seqs[ki][1]);
	  double xlen = domain_[1] - domain_[0];
	  shared_ptr<SplineCurve> 
	    line1(new SplineCurve(Point(p1[0]-xlen, p1[1]), -1,
				  Point(p1[0]+xlen, p1[1]), 1));

	  // For each segment in the first loop, check if it intersects the
	  // line and count the number of intersections to the left and to
	  // the right
	  vector<vector<double> > val_side(2);

	  for (kj=2; kj<seqs[0].size(); kj+=2)
	    {
	      if (std::min(seqs[0][kj-1], seqs[0][kj+1]) > p1[1] ||
		  std::max(seqs[0][kj-1], seqs[0][kj+1]) < p1[1])
		continue; // Not an intersection

	      double u1 = std::min(seqs[0][kj-2], seqs[0][kj]);
	      double u2 = std::max(seqs[0][kj-2], seqs[0][kj]);
	      int ix = -1;
	      if (u1 > p1[0])
		ix = 1;
	      else if (u2 < p1[0])
		ix = 0;
	      else
		{
		  // Intersect
		  vector<pair<double,double> > int_pts;
		  shared_ptr<SplineCurve> 
		    line2(new SplineCurve(Point(seqs[0][kj-2], seqs[0][kj-1]),
					  Point(seqs[0][kj], seqs[0][kj+1])));
		  intersectcurves(line1.get(), line2.get(), tol, int_pts);
		  for (kh=0; kh<int_pts.size(); ++kh)
		    {
		      if (int_pts[kh].first < -tol)
			ix = 0;
		      else if (int_pts[kh].first > tol)
			ix = 1;
		    }
		}
	      for (kr=0; kr<val_side[ix].size(); kr+=2)
		{
		  if (u1 > val_side[ix][kr]-eps2_ && u1 < val_side[ix][kr+1]+eps2_)
		    {
		      if (u2 > val_side[ix][kr+1]+eps2_)
			val_side[ix][kr+1] = u2;
		      break;
		    }
		  if (u2 > val_side[ix][kr]-eps2_ && u2 < val_side[ix][kr+1]+eps2_)
		    {
		      if (u1 < val_side[ix][kr]-eps2_)
			val_side[ix][kr] = u1;
		      break;
		    }
		}
	      if (kr == val_side[ix].size())
		{
		  val_side[ix].push_back(u1);
		  val_side[ix].push_back(u2);
		}
	    }

	  for (int ka=0; ka<2; ++ka)
	    {
	      if (val_side[ka].size() > 0)
		{
		  std::sort(val_side[ka].begin(), val_side[ka].end());
		  for (kr=2; kr<val_side[ka].size()-1; )
		    {
		      if (val_side[ka][kr]-val_side[ka][kr-1] < eps2_)
			val_side[ka].erase(val_side[ka].begin()+kr-1, val_side[ka].begin()+kr+1);
		    }
		}
	    }
	  int nmb_left = val_side[0].size()/2;
	  int nmb_right = val_side[1].size()/2;
	  if (nmb_left % 2 == 1 || nmb_right % 2 == 1)
	    {
	      // If one of the remainders are different from 1, there is an
	      // ambiguity, probably caused by a corner hit. Corner hits on 
	      // both sides may also cause false results. Worth a closer
	      // investigation
	      // Remove
	      seqs.erase(seqs.begin()+ki);
	    }
	  else
	    ++ki;
	}

      // Change largest sequence back to the original
      seqs[0].erase(seqs[0].end()-2, seqs[0].end());

      // If there is more than one loop left, the same exercise should be done 
      // for other combination of loops
    }
}
#endif
//==============================================================================
void TrimUtils::reOrganizeSeqs(vector<vector<double> >& seqs, bool outer_only)
//==============================================================================
{ 
  double tol = 1.0e-6;  // Intersection tolerance

  // Assume that the seqence with most elements is the outer one
  size_t ki, kj, kh, kr;
  for (ki=0; ki<seqs.size(); ++ki)
    for (kj=ki+1; kj<seqs.size(); ++kj)
      if (seqs[ki].size() < seqs[kj].size())
	std::swap(seqs[ki], seqs[kj]);

  if (outer_only)
    {
      // Close largest loop
      seqs[0].insert(seqs[0].end(), seqs[0].begin(), seqs[0].begin()+2);

      // Check if the remaining sequences lie internally to the first one
      for (ki=1; ki<seqs.size(); )
	{
	  // Select an arbitrary point and make a horizontal line through
	  // this point
	  Point p1(seqs[ki][0], seqs[ki][1]);
	  double xlen = domain_[1] - domain_[0];
	  shared_ptr<SplineCurve> 
	    line1(new SplineCurve(Point(p1[0]-xlen, p1[1]), -1,
				  Point(p1[0]+xlen, p1[1]), 1));

	  // For each segment in the first loop, check if it intersects the
	  // line and count the number of intersections to the left and to
	  // the right
	  int nmb_left = 0, nmb_right = 0;
	  for (kj=2; kj<seqs[0].size(); kj+=2)
	    {
	      if (std::min(seqs[0][kj-1], seqs[0][kj+1]) > p1[1] ||
		  std::max(seqs[0][kj-1], seqs[0][kj+1]) < p1[1])
		continue; // Not an intersection

	      double u1 = std::min(seqs[0][kj-2], seqs[0][kj]);
	      double u2 = std::max(seqs[0][kj-2], seqs[0][kj]);
	      if (u1 > p1[0])
		nmb_left += 2;
	      else if (u2 < p1[0])
		nmb_right += 2;
	      else if (fabs(p1[1]-seqs[0][kj-1]) < tol &&
		       fabs(p1[1]-seqs[0][kj+1]) < tol)
		{
		  nmb_left++;
		  nmb_right++;
		}
	      else
		{
		  // Intersect
		  vector<pair<double,double> > int_pts;
		  shared_ptr<SplineCurve> 
		    line2(new SplineCurve(Point(seqs[0][kj-2], seqs[0][kj-1]),
					  Point(seqs[0][kj], seqs[0][kj+1])));
		  intersectcurves(line1.get(), line2.get(), tol, int_pts);
		  for (kh=0; kh<int_pts.size(); ++kh)
		    {
		      if (int_pts[kh].first < -tol)
			nmb_left += 2;
		      else if (int_pts[kh].first > tol)
			nmb_right += 2;
		    }
		}
	    }


	  if (nmb_left % 4 != 0 || nmb_right % 4 != 0)
	    {
	      // If one of the remainders are different from 1, there is an
	      // ambiguity, probably caused by a corner hit. Corner hits on 
	      // both sides may also cause false results. Worth a closer
	      // investigation
	      // Remove
	      seqs.erase(seqs.begin()+ki);
	    }
	  else
	    ++ki;
	}

      // Change largest sequence back to the original
      seqs[0].erase(seqs[0].end()-2, seqs[0].end());

      // If there is more than one loop left, the same exercise should be done 
      // for other combination of loops
    }
}

#if 0
//==============================================================================
void TrimUtils::reOrganizeSeqs(vector<vector<double> >& seqs, vector<int>& outer)
//==============================================================================
{ 
  double tol = 1.0e-6;  // Intersection tolerance

  outer.resize(seqs.size());
  std::fill(outer.begin(), outer.end(), 0);

  // Assume that the seqence with most elements is the outer one
  size_t last = seqs.size();
  for (size_t kr=0; kr<seqs.size(); kr=last)
    {
      last = seqs.size();
      size_t ki, kj, kh;
      for (ki=kr; ki<seqs.size(); ++ki)
	for (kj=ki+1; kj<seqs.size(); ++kj)
	  if (seqs[ki].size() < seqs[kj].size())
	    std::swap(seqs[ki], seqs[kj]);


 
      // Close largest loop
      seqs[kr].insert(seqs[kr].end(), seqs[kr].begin(), seqs[kr].begin()+2);
      outer[kr] = 1;

      // Check if the remaining sequences lie internally to the first one
      for (ki=kr+1; ki<last; )
	{
	  // Select an arbitrary point and make a horizontal line through
	  // this point
	  Point p1(0.5*(seqs[ki][0]+seqs[ki][2]), 0.5*(seqs[ki][1]+seqs[ki][3]));
	  double xlen = domain_[1] - domain_[0];
	  shared_ptr<SplineCurve> 
	    line1(new SplineCurve(Point(p1[0]-xlen, p1[1]), -1,
				  Point(p1[0]+xlen, p1[1]), 1));

	  // For each segment in the first loop, check if it intersects the
	  // line and count the number of intersections to the left and to
	  // the right
	  vector<vector<double> > val_side(2);

	  for (kj=2; kj<seqs[kr].size(); kj+=2)
	    {
	      if (std::min(seqs[kr][kj-1], seqs[kr][kj+1]) > p1[1] ||
		  std::max(seqs[kr][kj-1], seqs[kr][kj+1]) < p1[1])
		continue; // Not an intersection

	      double u1 = std::min(seqs[kr][kj-2], seqs[kr][kj]);
	      double u2 = std::max(seqs[kr][kj-2], seqs[kr][kj]);
	      int ix = -1;
	      if (std::min(seqs[kr][kj-2], seqs[kr][kj]) > p1[0])
		ix = 1;
	      else if (std::max(seqs[kr][kj-2], seqs[kr][kj]) < p1[0])
		ix = 0;
	      else
		{
		  // Intersect
		  vector<pair<double,double> > int_pts;
		  shared_ptr<SplineCurve> 
		    line2(new SplineCurve(Point(seqs[kr][kj-2], seqs[kr][kj-1]),
					  Point(seqs[kr][kj], seqs[kr][kj+1])));
		  intersectcurves(line1.get(), line2.get(), tol, int_pts);
		  for (kh=0; kh<int_pts.size(); ++kh)
		    {
		      if (int_pts[kh].first < -tol)
			ix = 0;
		      else if (int_pts[kh].first > tol)
			ix = 1;
		    }
		}
	      for (kh=0; kh<val_side[ix].size(); kh+=2)
		{
		  if (u1 > val_side[ix][kh]-eps2_ && u1 < val_side[ix][kh+1]+eps2_)
		    {
		      if (u2 > val_side[ix][kh+1]+eps2_)
			val_side[ix][kh+1] = u2;
		      break;
		    }
		  if (u2 > val_side[ix][kh]-eps2_ && u2 < val_side[ix][kh+1]+eps2_)
		    {
		      if (u1 < val_side[ix][kh]-eps2_)
			val_side[ix][kh] = u1;
		      break;
		    }
		}
	      if (kh == val_side[ix].size())
		{
		  val_side[ix].push_back(u1);
		  val_side[ix].push_back(u2);
		}
	    }

	  for (int ka=0; ka<2; ++ka)
	    {
	      if (val_side[ka].size() > 0)
		{
		  std::sort(val_side[ka].begin(), val_side[ka].end());
		  for (kh=2; kh<val_side[ka].size()-1; )
		    {
		      if (val_side[ka][kh]-val_side[ka][kh-1] < eps2_)
			val_side[ka].erase(val_side[ka].begin()+kh-1, val_side[ka].begin()+kh+1);
		      else
			kh += 2;
		    }
		}
	    }
	  int nmb_left = val_side[0].size()/2;
	  int nmb_right = val_side[1].size()/2;
	  if (nmb_left % 2 == 0 || nmb_right % 2 == 0)
	    {
	      // Loop outside the current outer loop. Move to end
	      std::swap(seqs[ki], seqs[last-1]);
	      --last;
	    }
	  else
	    ++ki;
	}

    }
}
#endif
//==============================================================================
void TrimUtils::reOrganizeSeqs(vector<vector<double> >& seqs, vector<int>& outer)
//==============================================================================
{ 
  double tol = 1.0e-6;  // Intersection tolerance

  outer.resize(seqs.size());
  std::fill(outer.begin(), outer.end(), 0);

  // Assume that the seqence with most elements is the outer one
  size_t last = seqs.size();
  for (size_t kr=0; kr<seqs.size(); kr=last)
    {
      last = seqs.size();
      size_t ki, kj, kh;
      for (ki=kr; ki<seqs.size(); ++ki)
	for (kj=ki+1; kj<seqs.size(); ++kj)
	  if (seqs[ki].size() < seqs[kj].size())
	    std::swap(seqs[ki], seqs[kj]);

#ifdef DEBUG
      std::ofstream ofout("outer_seq.g2");
      (void)ofout.precision(15);
      ofout << "410 1 0 4 255 0 0 255" << std::endl;
      ofout << seqs[kr].size()/2-1 << std::endl;
      for (int kh1=0; kh1<=(int)seqs[kr].size()-4; kh1+=2)
	{
	  ofout << seqs[kr][kh1] << " " << seqs[kr][kh1+1] << " " << 0 << " ";
	  ofout << seqs[kr][kh1+2] << " " << seqs[kr][kh1+3] << " " << 0 << std::endl;
	}
#endif
 
      // Close largest loop
      seqs[kr].insert(seqs[kr].end(), seqs[kr].begin(), seqs[kr].begin()+2);
      outer[kr] = 1;

      // Check if the remaining sequences lie internally to the first one
      for (ki=kr+1; ki<last; )
	{
#ifdef DEBUG
	  std::ofstream ofc("curr_seq.g2");
	  (void)ofc.precision(15);
	  ofc << "410 1 0 4 0 255 0 255" << std::endl;
	  ofc << seqs[ki].size()/2-1 << std::endl;
	  for (int kh1=0; kh1<=(int)seqs[ki].size()-4; kh1+=2)
	    {
	      ofc << seqs[ki][kh1] << " " << seqs[ki][kh1+1] << " " << 0 << " ";
	      ofc << seqs[ki][kh1+2] << " " << seqs[ki][kh1+3] << " " << 0 << std::endl;
	    }
#endif
 
	  // Select a point and make a horizontal line through
	  // this point
	  for (kj=0; kj<seqs[ki].size(); kj+=4)
	    if (fabs(seqs[ki][kj+1]-seqs[ki][kj+3]) > tol)
	      break;
	  if (kj >= seqs[ki].size())
	    kj = 0;
	  Point p1(0.5*(seqs[ki][kj]+seqs[ki][kj+2]), 0.5*(seqs[ki][kj+1]+seqs[ki][kj+3]));
	  double xlen = domain_[1] - domain_[0];
	  shared_ptr<SplineCurve> 
	    line1(new SplineCurve(Point(p1[0]-xlen, p1[1]), -1,
				  Point(p1[0]+xlen, p1[1]), 1));

	  // For each segment in the first loop, check if it intersects the
	  // line and count the number of intersections to the left and to
	  // the right

	  int nmb_left = 0, nmb_right = 0;
	  for (kj=2; kj<seqs[kr].size(); kj+=2)
	    {
	      if (std::min(seqs[kr][kj-1], seqs[kr][kj+1]) > p1[1] ||
		  std::max(seqs[kr][kj-1], seqs[kr][kj+1]) < p1[1])
		continue; // Not an intersection

	      double u1 = std::min(seqs[0][kj-2], seqs[0][kj]);
	      double u2 = std::max(seqs[0][kj-2], seqs[0][kj]);
	      if (u1 > p1[0])
		nmb_left += 2;
	      else if (u2 < p1[0])
		nmb_right += 2;
	      else if (fabs(p1[1]-seqs[0][kj-1]) < tol &&
		       fabs(p1[1]-seqs[0][kj+1]) < tol)
		{
		  nmb_left++;
		  nmb_right++;
		}
	      else
		{
		  // Intersect
		  vector<pair<double,double> > int_pts;
		  shared_ptr<SplineCurve> 
		    line2(new SplineCurve(Point(seqs[0][kj-2], seqs[0][kj-1]),
					  Point(seqs[0][kj], seqs[0][kj+1])));
		  intersectcurves(line1.get(), line2.get(), tol, int_pts);
		  for (kh=0; kh<int_pts.size(); ++kh)
		    {
		      if (int_pts[kh].first < -tol)
			nmb_left += 2;
		      else if (int_pts[kh].first > tol)
			nmb_right += 2;
		    }
		}
	    }

	  if (nmb_left % 4 == 0 || nmb_right % 4 == 0)
	    {
	      // Loop outside the current outer loop. Move to end
	      std::swap(seqs[ki], seqs[last-1]);
	      --last;
	    }
	  else
	    ++ki;
	}

    }
}


//==============================================================================
void TrimUtils::extractOutsidePoint(shared_ptr<BoundedSurface>& surf,
				    vector<double>& outpoints,
				    int max_nmb)
//==============================================================================
{
  outpoints.clear();
  int del = surf->dimension() + 2;

  int ki, kj;
  double dom[4];
  for (ki=0; ki<2; ++ki)
    dom[2*ki] = dom[2*ki+1] = points_[ki];
  for (kj=1; kj<nmb_pts_; ++kj)
    for (ki=0; ki<2; ++ki)
      {
	dom[2*ki] = std::min(dom[2*ki], points_[kj*del+ki]);
	dom[2*ki+1] = std::max(dom[2*ki+1], points_[kj*del+ki]);
      }

  SubCloud cloud;
  cloud.setInfo(nmb_pts_, 0, nmb_pts_*del, dom, dom);
  extractOutsidePoint(surf, cloud, outpoints, max_nmb);
}

//==============================================================================
void TrimUtils::extractOutsidePoint(shared_ptr<BoundedSurface>& surf,
				    SubCloud& cloud,
				    vector<double>& outpoints,
				    int max_nmb)
//==============================================================================
{
  int del = surf->dimension() + 2;
  if (cloud.nmb_pts_ <= max_nmb)
    {
      //std::cout << "Bottom \n";
      // Test point for inclusion
      for (int ki=0; ki<cloud.nmb_pts_; ++ki)
	{
	  int ix = cloud.ix1_;
	  if (!surf->inDomain(points_[ix+ki*del], points_[ix+ki*del+1]))
	    outpoints.insert(outpoints.end(), points_+ix+ki*del, 
			     points_+ix+ki*del+del);
	}
    }
  else
    {
      // Distribute current point cloud into sub clouds
      int nmb_u = std::max(3, std::min(cloud.nmb_pts_/500, 10));
      int nmb_v = std::max(3, std::min(cloud.nmb_pts_/500, 10));
      vector<SubCloud> sub_clouds;
      distributePointCloud(cloud.ix1_, cloud.ix2_, cloud.dom_,
			   nmb_u, nmb_v, sub_clouds);

      // Traverse the sub clouds and check whether they need to be processed
      // further
      //std::cout << "Nmb pts: " << cloud.nmb_pts_ << ", nmb sub clouds: " << sub_clouds.size() << std::endl;
      int nmb_processed = 0;
      for (size_t kj=0; kj<sub_clouds.size(); ++kj)
	{
	  int ki;
	  double par[10];
	  par[0] = par[4] = sub_clouds[kj].bb_[0];
	  par[2] = par[6] = sub_clouds[kj].bb_[1];
 	  par[1] = par[3] = sub_clouds[kj].bb_[2];
	  par[5] = par[7] = sub_clouds[kj].bb_[3];
	  // par[8] = 0.5*(sub_clouds[kj].bb_[0]+sub_clouds[kj].bb_[1]);
	  // par[9] = 0.5*(sub_clouds[kj].bb_[2]+sub_clouds[kj].bb_[3]);
	  int nmb_in=0, nmb_out=0;
	  for (ki=0; ki<8; ki+=2)
	    {
	      if (surf->inDomain(par[ki], par[ki+1]))
		++nmb_in;
	      else
		++nmb_out;
	    }

	  /*if (nmb_out == 5)
	    outpoints.insert(outpoints.end(), points_.begin()+sub_clouds[kj].ix1_,
			  points_.begin()+sub_clouds[kj].ix2_);
			  else*/ if (nmb_out > 0)
	    {
	      //std::cout << "Sub cloud " << kj << ", nmb points " << sub_clouds[kj].nmb_pts_;
	      //std::cout << " Nmb outside pts: " << outpoints.size()/3 << std::endl;
	      extractOutsidePoint(surf, sub_clouds[kj], outpoints, max_nmb);
	      nmb_processed++;
	    }
	  // else
	  //   std::cout << "sub cloud not processed \n";
	}
#ifdef DEBUG
      std::cout << "Nmb pts: " << cloud.nmb_pts_ << ", nmb sub clouds: " << sub_clouds.size();
      std::cout << ", nmb processed. " << nmb_processed << ", nmb out: " <<  outpoints.size()/3 << std::endl;
#endif
    }
}

 
