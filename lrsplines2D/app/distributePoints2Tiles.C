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
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/utils/BoundingBox.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace Go;
using std::vector;
using std::string;

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

int main(int argc, char *argv[])
{
  if (argc != 6) {
    std::cout << "Usage: input meta file, output meta file, file root, nmb, overlap" << std::endl;
    return -1;
  }

  std::ifstream filein(argv[1]);
  std::ofstream fileout(argv[2]);
  const char* file_root(argv[3]);
  int nmb = atoi(argv[4]);
  double overlap = atof(argv[5]);
  overlap = std::min(overlap, 0.5);
  overlap = std::max(overlap, 0.0);

  // Read input meta data
  vector<string> inblock;
  vector<int> nmb_points;
  vector<double> domain;
  LRSurfApproxUtils::readBlockMeta(filein, inblock, nmb_points, domain);
  
  // Total domain of surface set
  double tot_domain[4];
  tot_domain[0] = tot_domain[2] = HUGE;
  tot_domain[1] = tot_domain[3] = -HUGE;
  for (size_t ki=0; ki<domain.size(); ki+=4)
    {
      tot_domain[0] = std::min(tot_domain[0], domain[ki]);
      tot_domain[1] = std::max(tot_domain[1], domain[ki+1]);
      tot_domain[2] = std::min(tot_domain[2], domain[ki+2]);
      tot_domain[3] = std::max(tot_domain[3], domain[ki+3]);
    }

  // Estimate number of surfaces assuming that the density of the 
  // compound point set is fairly uniform
  int nmb_blocks = nmb_points.size();
  int all_points = 0;  // Total number of points
  int ki, kj, kr;
  for (ki=0; ki<nmb_blocks; ++ki)
    all_points += nmb_points[ki];

  double fac = (tot_domain[1]-tot_domain[0])/(tot_domain[3]-tot_domain[2]);
  double nmb2 = (double)all_points/(double)nmb;
  double len = sqrt((double)nmb2/fac);
  int u_nmb = std::max(1, (int)(fac*len));
  int v_nmb = std::max(1, (int)len);

  // Fetch names of files containing the tiles. Make text file due to 
  // lack of knowledge of the number of points. This information has
  // to be added later
  vector<string> tiles;
  LRSurfApproxUtils::fetchFileNames(file_root, 2, u_nmb*v_nmb, tiles);
  for (ki=0; ki<(int)tiles.size(); ++ki)
    {
      FILE *fp = fopen(tiles[ki].c_str(), "w");
      fclose(fp);
    }
  vector<int> nmb_points_tiles(tiles.size(), 0);
  vector<double> domain_tiles(4*tiles.size(), 0.0);

  // Create output files

  // Make tile raster
  double udel = (tot_domain[1] - tot_domain[0])/(double)u_nmb;
  double vdel = (tot_domain[3] - tot_domain[2])/(double)v_nmb;
  double overlap_u = overlap*udel;
  double overlap_v = overlap*vdel;
  vector<double> u_val(3*u_nmb+1); 
  vector<double> v_val(3*v_nmb+1); 
  double upar, vpar;
  for (ki=0, kj=0, upar=tot_domain[0]; ki<u_nmb; ++ki, upar+=udel)
    {
      u_val[kj++] = upar;
      u_val[kj++] = upar + overlap_u;
      u_val[kj++] = upar + udel - overlap_u;
    }
  u_val[kj++] = tot_domain[1];

  for (ki=0, kj=0, vpar=tot_domain[2]; ki<v_nmb; ++ki, vpar+=vdel)
    {
      v_val[kj++] = vpar;
      v_val[kj++] = vpar + overlap_v;
      v_val[kj++] = vpar + vdel - overlap_v;
    }
  v_val[kj++] = tot_domain[3];

  int uel = (int)u_val.size()-1;
  int vel = (int)v_val.size()-1;
  vector<int> raster(uel*vel);

  // Distribute points
  // Perform for each input block
  int dim = 3;   // This might need to be reconsidered if working
  // with enhanced point clouds
  int pp0, pp1, pp2, pp3;
  for (ki=0; ki<nmb_blocks; ++ki)
    {
      std::ofstream of("pnt_grp.g2");
      (void)of.precision(15);

      // Read block of points. Currently the points are expected to be
      // stored in g2-format
      std::ifstream infile(inblock[ki].c_str());
      ObjectHeader header;
      header.read(infile);
      PointCloud3D points;
      points.read(infile);
      vector<double> data(points.rawData(), points.rawData()+3*nmb_points[ki]);

      // Identify relevant point groups
      size_t ux1, ux2, vx1, vx2;
      for (ux1=0; ux1<u_val.size()-1; ++ux1)
	if (domain[4*ki] < u_val[ux1+1])
	  break;
      for (ux2=ux1+1; ux2<u_val.size(); ++ux2)
	if (domain[4*ki+1] < u_val[ux2])
	  break;
      ux2 = std::min(ux2, u_val.size()-1);

      for (vx1=0; vx1<v_val.size()-1; ++vx1)
	if (domain[4*ki+2] < v_val[vx1+1])
	  break;
      for (vx2=vx1+1; vx2<v_val.size(); ++vx2)
	if (domain[4*ki+3] < v_val[vx2])
	  break;
      vx2 = std::min(vx2, v_val.size()-1);

      printf("Start sort block nr %d \n", ki);

     // Sort points in the y-direction
      qsort(&data[0], nmb_points[ki], 3*sizeof(double), compare_v_par);

      // For each slot in the y-direction, identify the points that belongs
      // to a slot in the x-direction.
      std::fill(raster.begin(), raster.end(), -1);
      for (pp0=0, kj=vx1+1, vpar=v_val[kj]; kj<=vx2; ++kj)
	{
	  vpar=v_val[kj];

	  // Identify points belonging to slot in y-direction
	  for (pp1=pp0; pp1<(int)data.size() && data[pp1+1]<vpar; pp1+=3);

	  // Sort points in x-direction
	  qsort(&data[0]+pp0, (pp1-pp0)/3, 3*sizeof(double), compare_u_par);
	  for (pp2=pp0, kr=ux1+1, upar=u_val[kr]; kr<=ux2; ++kr)
	    {
	      upar=u_val[kr];

	      // Identify points belonging to slot in y-direction
	      for (pp3=pp2; pp3<pp1 && data[pp3]<upar; pp3+=3);
	      if (kr == ux2-1)
	      	for (; pp3<pp1 && data[pp3] <= upar; pp3+=3);

	      // Set raster pointer
	      raster[(kj-1)*uel+kr-1] = pp2;

	      // if (pp3 > pp2)
	      // 	{
	      // 	  of << "400 1 0 0" << std::endl;
	      // 	  of << (pp3-pp2)/3 << std::endl;
	      // 	  for (int pp=pp2; pp<pp3; pp+=3)
	      // 	    of << data[pp] << " " << data[pp+1] << " " << data[pp+2] << std::endl;
	      // 	}
	      pp2 = pp3;
	    }
	  pp0 = pp1;
	} 
      printf("Read and sorted block nr %d \n", ki);

      // Write tiles to files
      int uu = std::max((int)ux1-1,0);
      int ux = 3*(uu/3)-1;
      int umax = 3*(ux2/3)+1;
      umax = std::min(umax, uel-1);
      int vv = std::max((int)vx1-1,0);
      int vx = 3*(vv/3)-1;
      int vmax = 3*(vx2/3)+1;
      vmax = std::min(vmax, vel-1);
      int file_ix = 0;
      int ix;
      int prev;
      for(; vx<vmax; vx+=3)
      	for (ux=3*(ux1/3)-1; ux<umax; ux+=3)
      	  {
      	    // Open output file
      	    file_ix = ((vx+1)/3)*u_nmb + ((ux+1)/3);
      	    FILE *fp = fopen(tiles[file_ix].c_str(), "a");

      	    prev = 0;
      	    for (kj=std::max(0, vx); kj<std::min(vel, vx+5); ++kj)
      	      for (kr=std::max(0, ux); kr<std::min(uel, ux+5); ++kr)
      		{
      		  // Write sequence of points to file
      		  ix = kj*uel+kr;
      		  int pmin = raster[ix];
      		  int pmax;
      		  if (pmin < 0 || pmin < prev)
      		    continue;
      		  if (ix+1 < (int)raster.size() && raster[ix+1] >= pmin)
      		    pmax = raster[ix+1];
      		  else 
      		    {
      		      int ix1;
      		      for (ix1=ix+1; ix1<(int)raster.size(); ++ix1)
      			if (raster[ix1] >= pmin)
      			  {
      			    pmax = raster[ix1];
      			    break;
      			  }
      		      if (ix1 == (int)raster.size())
      			pmax = data.size();
      		    }

      		  for (int pp=pmin; pp<pmax; pp+=3)
      		    fprintf(fp,"%7.12f %7.12f %7.12f \n", 
      			    data[pp], data[pp+1], data[pp+2]);
		  
      		  prev = pmin;
      		  nmb_points_tiles[file_ix] += (pmax-pmin)/3;
      		  int stop1 = 1;

      		}

      	    fclose(fp);
      	  }
      printf("Tiles from block nr %d to file\n", ki);
    }

  vector<string> tiles2;
  LRSurfApproxUtils::fetchFileNames(file_root, 1, u_nmb*v_nmb, tiles2);
  
  // Create g2 files
  for (ki=0; ki<(int)tiles.size(); ++ki)
    {
      std::ofstream tmp("tmp.txt");
      tmp << "400 1 0 0" << std::endl;
      tmp << nmb_points_tiles[ki] << std::endl;

      std::ifstream if_a("tmp.txt");
      std::ifstream if_b(tiles[ki].c_str());
      std::ofstream of_c(tiles2[ki].c_str());
      of_c << if_a.rdbuf() << if_b.rdbuf();
    }
  printf("g2-files created \n");

  // Local domains
  for (ki=0; ki<u_nmb; ++ki)
    domain_tiles[4*ki+2] = v_val[0];
  for (kj=1, kr=3; kj<v_nmb; ++kj, kr+=3)
    {
      for (ki=0; ki<u_nmb; ++ki)
  	{
  	  domain_tiles[4*((kj-1)*u_nmb+ki)+3] = v_val[kr+1];
  	  domain_tiles[4*(kj*u_nmb+ki)+2] = v_val[kr-1];
  	}
    }
  for (ki=0; ki<u_nmb; ++ki)
    domain_tiles[4*(u_nmb*(v_nmb-1)+ki)+3] = v_val[vel];

  for (kj=0; kj<v_nmb; ++kj)
    domain_tiles[4*kj*u_nmb] = u_val[0];
  for (ki=1, kr=3; ki<u_nmb; ++ki, kr+=3)
    {
      for (kj=0; kj<v_nmb; ++kj)
  	{
  	  domain_tiles[4*(kj*u_nmb+ki-1)+1] = u_val[kr+1];
  	  domain_tiles[4*(kj*u_nmb+ki)] = u_val[kr-1];
  	}
    }
  for (kj=0; kj<v_nmb; ++kj)
    domain_tiles[4*(kj*u_nmb+u_nmb-1)+1] = u_val[uel];

  printf("Domain tiles set \n");

  // Write metadata
  LRSurfApproxUtils::writeTileMeta(fileout, tot_domain, u_nmb, v_nmb,
  				   overlap_u, overlap_v, nmb_points_tiles,
  				   domain_tiles, tiles2);
  printf("Finished \n");
  
}

