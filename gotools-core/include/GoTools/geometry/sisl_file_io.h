//===========================================================================
//                                                                           
// File: sisl_file_io.h                                                      
//                                                                           
// Created: Thu Mar 25 15:26:54 2004                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: sisl_file_io.h,v 1.4 2009-05-13 07:32:30 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _SISL_FILE_IO_H
#define _SISL_FILE_IO_H


#include <memory>
#include <vector>

#include <fstream>

#include "GoTools/utils/config.h"

#ifdef __BORLANDC__
#include <cstdio>
#define STD_FILE std::FILE
#else
#define STD_FILE FILE
#endif

struct SISLObject;

void read_non_comment(STD_FILE* fp, char* string);

void file_to_obj(STD_FILE* fp, SISLObject** wo, int* jstat);

void curve_to_file(STD_FILE *f, struct SISLCurve *c1);

void surface_to_file(STD_FILE *f, struct SISLSurf *surf);

int get_next_surface(STD_FILE *fp, SISLSurf **qc);

int get_sisl_surfaces(STD_FILE* fp,
		      std::vector<shared_ptr<SISLSurf> >& sisl_sfs);


#endif // _SISL_FILE_IO_H

