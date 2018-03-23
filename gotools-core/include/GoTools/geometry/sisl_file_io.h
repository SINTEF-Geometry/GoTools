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

#ifndef _SISL_FILE_IO_H
#define _SISL_FILE_IO_H


#include <memory>
#include <vector>

#include <fstream>

#include "sislP.h"
#include "GoTools/utils/config.h"

#ifdef __BORLANDC__
#include <cstdio>
#define STD_FILE std::FILE
#else
#define STD_FILE FILE
#endif

//struct SISLObject;

void read_non_comment(STD_FILE* fp, char* string);

void file_to_obj(STD_FILE* fp, SISLObject** wo, int* jstat);

void curve_to_file(STD_FILE *f, struct SISLCurve *c1);

void surface_to_file(STD_FILE *f, struct SISLSurf *surf);

int get_next_surface(STD_FILE *fp, SISLSurf **qc);

int get_sisl_surfaces(STD_FILE* fp,
		      std::vector<shared_ptr<SISLSurf> >& sisl_sfs);


#endif // _SISL_FILE_IO_H

