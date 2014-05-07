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

#include "GoTools/compositemodel/PointSetApp.h"
#include "rply.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h> // For atof()

#define DEBUG

using namespace std;
using namespace Go;

/* ----------------------------------------------------------------------
 * Make sure we get our integer types right
 * ---------------------------------------------------------------------- */
#if defined(_MSC_VER) && (_MSC_VER < 1600)
/* C99 stdint.h only supported in MSVC++ 10.0 and up */
typedef __int8 t_ply_int8;
typedef __int16 t_ply_int16;
typedef __int32 t_ply_int32;
typedef unsigned __int8 t_ply_uint8;
typedef unsigned __int16 t_ply_uint16;
typedef unsigned __int32 t_ply_uint32;
#define PLY_INT8_MAX (127)
#define PLY_INT8_MIN (-PLY_INT8_MAX-1)
#define PLY_INT16_MAX (32767)
#define PLY_INT16_MIN (-PLY_INT16_MAX-1)
#define PLY_INT32_MAX (2147483647)
#define PLY_INT32_MIN (-PLY_INT32_MAX-1)
#define PLY_UINT8_MAX (255)
#define PLY_UINT16_MAX (65535)
#define PLY_UINT32_MAX  (4294967295)
#else
#include <stdint.h>
typedef int8_t t_ply_int8;
typedef int16_t t_ply_int16;
typedef int32_t t_ply_int32;
typedef uint8_t t_ply_uint8;
typedef uint16_t t_ply_uint16;
typedef uint32_t t_ply_uint32;
#define PLY_INT8_MIN INT8_MIN
#define PLY_INT8_MAX INT8_MAX
#define PLY_INT16_MIN INT16_MIN
#define PLY_INT16_MAX INT16_MAX
#define PLY_INT32_MIN INT32_MIN
#define PLY_INT32_MAX INT32_MAX
#define PLY_UINT8_MAX UINT8_MAX
#define PLY_UINT16_MAX UINT16_MAX
#define PLY_UINT32_MAX UINT32_MAX
#endif

/* ----------------------------------------------------------------------
 * Constants 
 * ---------------------------------------------------------------------- */
#define WORDSIZE 256
#define LINESIZE 1024
#define BUFFERSIZE (8*1024)

typedef enum e_ply_io_mode_ {
    PLY_READ, 
    PLY_WRITE
} e_ply_io_mode;

static const char *const ply_storage_mode_list[] = {
    "binary_big_endian", "binary_little_endian", "ascii", NULL
};     /* order matches e_ply_storage_mode enum */

static const char *const ply_type_list[] = {
    "int8", "uint8", "int16", "uint16", 
    "int32", "uint32", "float32", "float64",
    "char", "uchar", "short", "ushort", 
    "int", "uint", "float", "double",
    "list", NULL
};     /* order matches e_ply_type enum */

/* ----------------------------------------------------------------------
 * Property reading callback argument
 *
 * element: name of element being processed
 * property: name of property being processed
 * nelements: number of elements of this kind in file
 * instance_index: index current element of this kind being processed
 * length: number of values in current list (or 1 for scalars)
 * value_index: index of current value int this list (or 0 for scalars)
 * value: value of property
 * pdata/idata: user data defined with ply_set_cb
 *
 * Returns handle to PLY file if succesful, NULL otherwise.
 * ---------------------------------------------------------------------- */
typedef struct t_ply_argument_ {
    p_ply_element element;
    long instance_index;
    p_ply_property property;
    long length, value_index;
    double value;
    void *pdata;
    long idata;
} t_ply_argument;

/* ----------------------------------------------------------------------
 * Property information
 *
 * name: name of this property
 * type: type of this property (list or type of scalar value)
 * length_type, value_type: type of list property count and values
 * read_cb: function to be called when this property is called
 *
 * Returns 1 if should continue processing file, 0 if should abort.
 * ---------------------------------------------------------------------- */
typedef struct t_ply_property_ {
    char name[WORDSIZE];
    e_ply_type type, value_type, length_type;
    p_ply_read_cb read_cb;
    void *pdata;
    long idata;
} t_ply_property; 

/* ----------------------------------------------------------------------
 * Element information
 *
 * name: name of this property
 * ninstances: number of elements of this type in file
 * property: property descriptions for this element
 * nproperty: number of properties in this element
 *
 * Returns 1 if should continue processing file, 0 if should abort.
 * ---------------------------------------------------------------------- */
typedef struct t_ply_element_ {
    char name[WORDSIZE];
    long ninstances;
    p_ply_property property;
    long nproperties;
} t_ply_element;

/* ----------------------------------------------------------------------
 * Input/output driver
 *
 * Depending on file mode, different functions are used to read/write 
 * property fields. The drivers make it transparent to read/write in ascii, 
 * big endian or little endian cases.
 * ---------------------------------------------------------------------- */
typedef int (*p_ply_ihandler)(p_ply ply, double *value);
typedef int (*p_ply_ichunk)(p_ply ply, void *anydata, size_t size);
typedef struct t_ply_idriver_ {
    p_ply_ihandler ihandler[16];
    p_ply_ichunk ichunk;
    const char *name;
} t_ply_idriver;
typedef t_ply_idriver *p_ply_idriver;

typedef int (*p_ply_ohandler)(p_ply ply, double value);
typedef int (*p_ply_ochunk)(p_ply ply, void *anydata, size_t size);
typedef struct t_ply_odriver_ {
    p_ply_ohandler ohandler[16];
    p_ply_ochunk ochunk;
    const char *name;
} t_ply_odriver;
typedef t_ply_odriver *p_ply_odriver;

/* ----------------------------------------------------------------------
 * Ply file handle. 
 *
 * io_mode: read or write (from e_ply_io_mode)
 * storage_mode: mode of file associated with handle (from e_ply_storage_mode)
 * element: elements description for this file
 * nelement: number of different elements in file
 * comment: comments for this file
 * ncomments: number of comments in file
 * obj_info: obj_info items for this file
 * nobj_infos: number of obj_info items in file
 * fp: file pointer associated with ply file
 * rn: skip extra char after end_header? 
 * buffer: last word/chunck of data read from ply file
 * buffer_first, buffer_last: interval of untouched good data in buffer
 * buffer_token: start of parsed token (line or word) in buffer
 * idriver, odriver: input driver used to get property fields from file 
 * argument: storage space for callback arguments
 * welement, wproperty: element/property type being written
 * winstance_index: index of instance of current element being written
 * wvalue_index: index of list property value being written 
 * wlength: number of values in list property being written
 * error_cb: error callback
 * pdata/idata: user data defined with ply_open/ply_create
 * ---------------------------------------------------------------------- */
typedef struct t_ply_ {
    e_ply_io_mode io_mode;
    e_ply_storage_mode storage_mode;
    p_ply_element element;
    long nelements;
    char *comment;
    long ncomments;
    char *obj_info;
    long nobj_infos;
    FILE *fp;
    int rn;
    char buffer[BUFFERSIZE];
    size_t buffer_first, buffer_token, buffer_last;
    p_ply_idriver idriver;
    p_ply_odriver odriver;
    t_ply_argument argument;
    long welement, wproperty;
    long winstance_index, wvalue_index, wlength;
    p_ply_error_cb error_cb;
    void *pdata;
    long idata;
} t_ply;



vector<double> vertices;
vector<int> triangles;

  static int vertex_cb(p_ply_argument argument) {
    long eol;
    ply_get_argument_user_data(argument, NULL, &eol);
    //sprintf(s1, "%g", ply_get_argument_value(argument));
    //vertices[v1++] =  ply_get_argument_value(argument);
    vertices.push_back(ply_get_argument_value(argument));
    //if (eol) sprintf(s1,"\n");
    //else sprintf(s1," ");
    return 1;
  }

  static int face_cb(p_ply_argument argument) {
    long length, value_index;
    ply_get_argument_property(argument, NULL, &length, &value_index);
    switch (value_index) {
    case 0:
    case 1: 
      triangles.push_back(ply_get_argument_value(argument));
      break;
    case 2:
      triangles.push_back(ply_get_argument_value(argument));
      break;
    default: 
      break;
    }
    return 1;
  }

int main( int argc, char* argv[] )
{
  if (argc != 3) {
    std::cout << "Input parameters : Input triangulation file(.ply), output parameterized points (.txt)"  << std::endl;
    exit(-1);
  }

  // Read input arguments
  char *input = argv[1];
  std::ofstream os(argv[2]);

  p_ply ply = ply_open(input, NULL, 0, NULL);

  if (!ply)
    return -1;
  if (!ply_read_header(ply))
    return -1;
  
  // Read input file
  long nvertices = ply_set_read_cb(ply, "vertex", "x", vertex_cb, NULL, 0);
  ply_set_read_cb(ply, "vertex", "y", vertex_cb, NULL, 0);
  ply_set_read_cb(ply, "vertex", "z", vertex_cb, NULL, 1);
  long ntriangles = ply_set_read_cb(ply, "face", "vertex_indices", face_cb, 
				    NULL, 0);
  vertices.reserve(nvertices*3);
  triangles.reserve(ntriangles*3);
  if (!ply_read(ply)) 
    return -1;

  int ki;

  streamsize prev = os.precision(15);
#ifdef DEBUG
  std::ofstream of("triangulation.g2");
  of << "400 1 0 0" << std::endl;
  of << nvertices << std::endl;
  for (ki=0; ki<nvertices; ++ki)
    
  for (ki=0; ki<nvertices; ++ki)
      of << vertices[3*ki] << " " << vertices[3*ki+1] << " " << vertices[3*ki+2] << std::endl;
  of << "410 1 0 0" << std::endl;
  of << 3*ntriangles << std::endl;
  for (ki=0; ki<ntriangles; ++ki)
    {
      int kj = triangles[3*ki];
      int kr = triangles[3*ki+1];
      int kh = triangles[3*ki+2];
      of << vertices[3*kj] << " " << vertices[3*kj+1] << " " << vertices[3*kj+2] << " ";
      of << vertices[3*kr] << " " << vertices[3*kr+1] << " " << vertices[3*kr+2] << std::endl;
      of << vertices[3*kj] << " " << vertices[3*kj+1] << " " << vertices[3*kj+2] << " ";
      of << vertices[3*kh] << " " << vertices[3*kh+1] << " " << vertices[3*kh+2] << std::endl;
      of << vertices[3*kh] << " " << vertices[3*kh+1] << " " << vertices[3*kh+2] << " ";
      of << vertices[3*kr] << " " << vertices[3*kr+1] << " " << vertices[3*kr+2] << std::endl;
    }
#endif
  
  // Parameterize
  vector<double> param;
  PointSetApp::parameterizeTriang(&vertices[0], nvertices, 
				  &triangles[0], ntriangles,
  				  param);

  for (ki=0; ki<nvertices; ++ki)
    {
      os << param[2*ki] << ", " << param[2*ki+1] << ", ";
      os << vertices[3*ki] << ", " << vertices[3*ki+1] << ", " << vertices[3*ki+2];
      os << std::endl;
    }
}

