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

#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/utils/Array.h"
#include "GoTools/geometry/sisl_file_io.h"
#include "sislP.h"
#include <fstream>
#include <stdlib.h>
#include <stdio.h>


#ifndef FALSE
#define FALSE 0
#endif




int main(int argc, char** argv)
{
  if (argc!=3)
    {
      std::cerr << "usage: <infile name> <outfile name>\n";
      return -1;
    }

  FILE* from_file = fopen(argv[1], "r" );
  if (!from_file)
    {
      std::cout << "File does not exist!" << std::endl;
      return -1;
    }
  std::ofstream outfile(argv[2]);

  int nmb_objs = 1;
  // char string[80];       /* character to read $ which identifies a comment */
  //read_non_comment(from_file,string);
  //sscanf(string,"%d",&nmb_objs);
  int kstat = 0;
  int ki;
  for (ki = 0; ki < nmb_objs; ++ki)
    {
      SISLObject* ob=NULL;
      file_to_obj(from_file, &ob, &kstat);
      if (kstat<0)
	{
	  std::cerr << "file_to_obj returned :" << kstat << std::endl;
	  return -1;
	}
      if (ob && ob->iobj==2)
	{
	  Go::SplineSurface* gs=Go::SISLSurf2Go(ob->s1);
	  gs->writeStandardHeader(outfile);
	  gs->write (outfile);
	}
      else if (ob && ob->iobj==1)
	{
	  Go::SplineCurve* gc=Go::SISLCurve2Go(ob->c1);
	  gc->writeStandardHeader(outfile);
	  gc->write (outfile);
	}
      else
	{
	  std::cerr << "Nothing to write" << std::endl;	  
	}
      if (ob) delete ob;
    }

  if (from_file) fclose(from_file);
}
