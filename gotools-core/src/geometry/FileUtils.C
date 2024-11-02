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


#include "GoTools/geometry/FileUtils.h"
#include "GoTools/geometry/Utils.h"
#include <fstream>
#include <string.h>

using namespace Go;

int compare(const char *str1, char str2[][8], int nmb)
{
  for (int ki=0; ki<nmb; ++ki)
    if (strcmp(str1, str2[ki]) == 0)
      return ki;
  return -1;
}

//==============================================================================
int FileUtils::fileType(char *file, char keys[][8], int nmb_keys)
//==============================================================================
{
  //  Find file extension
  char *loc;
  char *last = 0;
  loc = strchr(file, '.');
  while (loc != NULL)
    {
      last = loc;
      loc = strchr(loc+1, '.');
    }
  if (last == NULL)
    {
      return -1;
    }
  char *input_type = last+1;

  // Check type
  int type;
  try {
    type = compare(input_type, keys, nmb_keys);
  }
  catch (...)
    {
      return -1;
    }

  return type;
}

//==============================================================================
void FileUtils::readTxtPointFile(std::ifstream& is, int del,
				 std::vector<double>& data, int& nmb_pts,
				 std::vector<double>& extent)
//==============================================================================
{
   if (!is.good())
    THROW("Invalid file!");

  // Read points
  nmb_pts = 0;
  char xx;
  extent.resize(2*del);
  for (int ki=0; ki<del; ++ki)
    {
      extent[2*ki] = std::numeric_limits<double>::max();
      extent[2*ki+1] = std::numeric_limits<double>::lowest();
    }

  // Check if the file has got a header
  char firstline[80];
  is >> xx;
  if (isdigit(xx) || xx == '-')
    {
      is.putback(xx);
    }
  else
    is >> firstline;

  while (!is.eof())
    {
      double tmp;
      is >> xx;
      if (!(isdigit(xx) || xx == '-'))
	{
	  is >> firstline;
	  Utils::eatwhite(is);
	}
      else
	{
	  is.putback(xx);
	  is >> tmp;
	  data.push_back(tmp);
	  extent[0] = std::min(tmp, extent[0]);
	  extent[1] = std::max(tmp, extent[1]);
	  for (int ki=1; ki<del; ++ki)
	    {
	      is >> xx;
	      if (xx != ',' && xx != ';')
		is.putback(xx);
	      is >> tmp;
	      data.push_back(tmp);
	      extent[2*ki] = std::min(tmp, extent[2*ki]);
	      extent[2*ki+1] = std::max(tmp, extent[2*ki+1]);
	    }
	  nmb_pts++;
	  Utils::eatwhite(is);
	}
    }
}
 
