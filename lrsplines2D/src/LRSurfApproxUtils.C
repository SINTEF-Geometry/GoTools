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
  char cc;
  is >> cc;
  int nmb;
  double xx;
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
	  is >> cc;  // Expects ]
	  is >> cc;  // Expects }
	  is >> cc;  // Expects , or ]
	}
    }
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
  os << "]" << std::endl;
}

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
