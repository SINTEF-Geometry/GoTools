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

#ifndef TRIMUTILS_H
#define TRIMUTILS_H

#include "GoTools/geometry/BoundedSurface.h"
#include <vector>
#include <math.h>

namespace Go
{
  class TrimUtils
  {
  public:
    // Note that the sequence of the points will be changed
    TrimUtils(double* points, int nmb_pts, int dim);

    TrimUtils(double* points, int nmb_pts, int dim, double domain[]);

    ~TrimUtils();

    void computeTrimSeqs(int max_level, int nmb_div,
			 std::vector<std::vector<double> >& seqs,
			 bool outer_only = true);

    void computeAllTrimSeqs(int max_level, int nmb_div,
			    std::vector<std::vector<std::vector<double> > >& seqs);

    double getDomainDiag()
    {
      return sqrt(del_u_*del_u_ + del_v_*del_v_);
    }

    void getDomainLengths(double& del_u, double& del_v)
    {
      del_u = del_u_;
      del_v = del_v_;
    }

    void extractOutsidePoint(shared_ptr<BoundedSurface>& surf,
			     std::vector<double>& outpoints,
			     int max_nmb);

  private:
    struct SubCloud
    {
      SubCloud()
      {
	nmb_pts_ = 0;
	ix1_ = 0; 
	ix2_ = 0;
	processed_ = false;
	limit_[0] = limit_[1] = limit_[2] = limit_[3] = 0;
      }

      void setInfo(int nmb_pts, int ix1, int ix2, double dom[], double bb[])
      {
	nmb_pts_ = nmb_pts;
	ix1_ = ix1;
	ix2_ = ix2;
	for (int ki=0; ki<4; ++ki)
	  {
	    dom_[ki] = dom[ki];
	    bb_[ki] = bb[ki];
	  }
      }

      void leftBd(std::vector<double>& bd)
      {
	bd.push_back(dom_[0]);
	bd.push_back(dom_[3]);
	bd.push_back(dom_[0]);
	bd.push_back(dom_[2]);
      }
    
      void rightBd(std::vector<double>& bd)
      {
	bd.push_back(dom_[1]);
	bd.push_back(dom_[2]);
	bd.push_back(dom_[1]);
	bd.push_back(dom_[3]);
      }

      void lowerBd(std::vector<double>& bd)
      {
	bd.push_back(dom_[0]);
	bd.push_back(dom_[2]);
	bd.push_back(dom_[1]);
	bd.push_back(dom_[2]);
      }

      void upperBd(std::vector<double>& bd)
      {
	bd.push_back(dom_[1]);
	bd.push_back(dom_[3]);
	bd.push_back(dom_[0]);
	bd.push_back(dom_[3]);
      }

      bool limitedSupport(double frac)
      {
	if ((bb_[1]-bb_[0])/(dom_[1]-dom_[0]) < frac)
	  return true;
	if ((bb_[3]-bb_[2])/(dom_[3]-dom_[2]) < frac)
	  return true;
	return false;
      }

      double dom_[4];
      double bb_[4];
      int nmb_pts_;
      int ix1_;
      int ix2_;
      bool processed_;
      short limit_[4];   // Sequence: left, right, lower, upper
    };

    double eps2_;   // Equality tolerance
    double *points_;
    int nmb_pts_;
    int dim_;
    double domain_[4];
    double del_u_;
    double del_v_;

    void computeTrimInfo(SubCloud& cloud,
			 int max_level, int nmb_div,
			 std::vector<std::vector<double> >& seqs);

    void distributePointCloud(int ix1, int ix2,
			      double domain[4],
			      int nmb_u, int nmb_v,
			      std::vector<SubCloud>& sub_clouds);

    void setSubSeq(SubCloud& cloud,
		   std::vector<double>& seq);
    
    void removeFalseTrimSeqs(std::vector<double>& limitseqs,
			     std::vector<std::vector<double> >& seqs);

    void mergeTrimSeqs(std::vector<std::vector<double> >& seqs,
		       std::vector<std::vector<double> >& seqs2);

    void cleanTrimResults(int limitsize,
			  std::vector<std::vector<double> >& seqs);

    void reOrganizeSeqs(std::vector<std::vector<double> >& seqs,
			bool outer_only);

    void reOrganizeSeqs(std::vector<std::vector<double> >& seqs,
			std::vector<int>& outer);

   void extractOutsidePoint(shared_ptr<BoundedSurface>& surf,
			    SubCloud& cloud,
			    std::vector<double>& outpoints,
			    int max_nmb);

   };
};

#endif
