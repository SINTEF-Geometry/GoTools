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

#include <vector>
#include <fstream>
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/errormacros.h"
#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/geometry/Utils.h"

using namespace Go;
using namespace std;

#define MAX_COLORS 12
int colors[MAX_COLORS][3] = {
  {255, 0, 0},
  {0, 255, 0},
  {0, 0, 255},
  {255, 255, 0},
  {255, 0, 255},
  {0, 255, 255},
  {128, 255, 0},
  {255, 128, 0},
  {128, 0, 255},
  {255, 0, 128},
  {0, 128, 255},
  {0, 255, 128},
};

int main(int argc, char* argv[] )
{

  if (argc < 3)
    {
      cout << "Usage: " << argv[0] << " infile outfile [color_info]" << endl;
      exit(-1);
    }
  // Color info, lists are cyclic:
  //   (empty) - no colors
  //   0 - Use all colors
  //   1 n - Use n first colors
  //   2 c1 c2 c3 ... cn - Use given list of colors
  //   3 n r1 g1 b1 r2 g2 b2 ... rn bn gn c1 c2 ... cm - Use list c1...cn of colors, where color_i has rgb-values ri gi bi (range 0...255)

  bool use_colors;
  int cycle_size = 0;
  int col_table_size;
  vector<int> col_idx_table;
  vector<int> col_r_table, col_g_table, col_b_table;

  use_colors = (argc > 3);
  if (use_colors)
    {
      col_r_table.resize(MAX_COLORS);
      col_g_table.resize(MAX_COLORS);
      col_b_table.resize(MAX_COLORS);
      col_idx_table.resize(MAX_COLORS);
      for (int i = 0; i < MAX_COLORS; ++i)
	{
	  col_r_table[i] = colors[i][0];
	  col_g_table[i] = colors[i][1];
	  col_b_table[i] = colors[i][2];
	  col_idx_table[i] = i;
	}

      int cycle_info_pos;
      switch(atoi(argv[3]))
	{
	case 0:
	  cycle_size = MAX_COLORS;
	  break;

	case 1:
	  if (argc < 5)
	    {
	      use_colors = false;
	      break;
	    }
	  cycle_size = min(MAX_COLORS, max(0, atoi(argv[4])));
	  break;

	case 2:
	  if (argc < 5)
	    {
	      use_colors = false;
	      break;
	    }
	  cycle_size = argc - 4;
	  col_idx_table.resize(cycle_size);
	  for (int i = 0; i < cycle_size; ++i)
	    col_idx_table[i] = min(MAX_COLORS - 1, max(0, atoi(argv[i+4])));
	  break;

	case 3:
	  if (argc < 5)
	    {
	      use_colors = false;
	      break;
	    }
	  col_table_size = atoi(argv[4]);
	  if (argc - 5 <= 3 * col_table_size)
	    {
	      use_colors = false;
	      break;
	    }
	  col_r_table.resize(col_table_size);
	  col_g_table.resize(col_table_size);
	  col_b_table.resize(col_table_size);
	  for (int i = 0; i < col_table_size; ++i)
	    {
	      col_r_table[i] = atoi(argv[5+3*i]);
	      col_g_table[i] = atoi(argv[6+3*i]);
	      col_b_table[i] = atoi(argv[7+3*i]);
	    }
	  cycle_info_pos = 5 + 3 * col_table_size;
	  cycle_size = argc - cycle_info_pos;
	  col_idx_table.resize(cycle_size);
	  for (int i = 0; i < cycle_size; ++i)
	    col_idx_table[i] = min(MAX_COLORS, max(0, atoi(argv[i+cycle_info_pos])));
	  break;

	default:
	  use_colors = false;
	}
    }

  // Open input volume file
  ifstream is(argv[1]);
  ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

  // Open outfile
  ofstream os(argv[2]);
  ALWAYS_ERROR_IF(os.bad(), "Bad or no output filename");

  int col_pos = 0;

  while (!is.eof())
    {

      // Read volume from file
      ObjectHeader head;
      SplineVolume vol;
      is >> head >> vol;

      SplineSurface* ss;

      int col = 0;
      if (use_colors)
	{
	  col = col_idx_table[col_pos];
	  ++col_pos;
	  if (col_pos == cycle_size)
	    col_pos = 0;
	}

      for (int i = 0; i < 3; ++i)
	for (int j = 0; j < 2; ++j)
	  {
	    if (j == 0)
	      ss = vol.constParamSurface(vol.startparam(i), i);
	    else
	      ss = vol.constParamSurface(vol.endparam(i), i);

	    if (use_colors == 0)
	      ss->writeStandardHeader(os);
	    else
	      os << "200 1 0 4 " << col_r_table[col] << " " << col_g_table[col] << " " << col_b_table[col] << " 255" << endl;
	    ss->write(os);

	    delete ss;
	  }

      Utils::eatwhite(is);
    }
}
