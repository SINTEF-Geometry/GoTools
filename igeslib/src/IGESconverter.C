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

#define IGESLIB_DEBUG

#include "GoTools/igeslib/IGESconverter.h"

#include <stdlib.h>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <sstream>
#include <vector>
#include <memory>
// #include "errno.h"

//#ifdef __BORLANDC__
#include <iterator>
//#endif

#include "sislP.h"
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/orientCurves.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/utils/Values.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/SplineDebugUtils.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/Ellipse.h"
#include "GoTools/geometry/Line.h"
#include "GoTools/utils/RotatedBox.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/creators/CurveCreators.h"


using namespace Go;
using std::vector;
using std::string;
using std::ostringstream;
using std::back_inserter;
using std::ostream;
using std::istream;
using std::map;


//*****************************************************************************
//
// NOTE: In IGES docs etc. columns are referred to from 1 to 80.
//       Here, in comments we refer to columns 0..79 instead.
//
//*****************************************************************************


// Global utility fctn
inline void pad(string& s, int line_length = 64, char filler = ' ')
{
    int sl = (int)s.length();
    int l = sl/line_length;
    if (l*line_length != sl)
	s += string((l+1)*line_length - sl, filler);
}


//-----------------------------------------------------------------------------
IGESheader::IGESheader()
    : pardel(','),
      recdel(';'),
      sending_system_prodid("Unknown"),
      original_filename("Unknown"),
      creator_system_prodid("IGES converter 1.1"),
      preprocessor_version_number("1.1"),
      integer_bits(32),
      single_prec_magn(38),
      single_prec_sign(6),
      double_prec_magn(308),
      double_prec_sign(15),
      receiving_system_prodid("Unknown"),
      model_space_scale(1.0),
      // Default unit is mm
      unit_flag(2),
      unit_description("MM"),
      max_weight_grads(1),
      max_line_width(0.01),
      timestamp_file_changed("800101.120000"),
      min_resolution(0.01),
      max_coordinate(1.0), // Should be changed when actually used!
      author("Unknown author"),
      organisation("SINTEF Applied Mathematics"),
      IGES_version(8),
      drafting_standard(0),
      timestamp_model_changed("800101.120000"),
      num_parameters(25)
//-----------------------------------------------------------------------------
{
}


//-----------------------------------------------------------------------------
bool EntityList::validEntity(int ent)
//-----------------------------------------------------------------------------
{
    for (size_t ki = 0; ki < entities_.size(); ++ki)
	if (entities_[ki] == ent)
	    return true;

    return false;
}


//-----------------------------------------------------------------------------
IGESconverter::IGESconverter()
    : filled_with_data_(false), geom_(), group_()
//-----------------------------------------------------------------------------
{
    GoTools::init();
}


//-----------------------------------------------------------------------------
IGESconverter::~IGESconverter()
//-----------------------------------------------------------------------------
{
}


//-----------------------------------------------------------------------------
void IGESconverter::readIGES(istream& is)
//-----------------------------------------------------------------------------
{

    // An IGES file consists of five sections. We read the content of each
    // section into a string, while checking that the line numbers are correct.

    char line_buffer[300]; // Only 81 chars will be filled, but we're safing...
    int line_number = 0;
    IGESSection sect;
    string sbufs[5];
    sbufs[0]=sbufs[1]=sbufs[2]=sbufs[3]=sbufs[4]="";
    num_lines_[0]=num_lines_[1]=num_lines_[2]=num_lines_[3]=num_lines_[4]=0;
    int Pcurr;
    // The P section is usually quite large, so we reserve a megabyte
    // of memory for it here.
    sbufs[P].reserve(1000000);
    while (readSingleIGESLine(is, line_buffer, line_number, sect) &&
	   sect < E) 
      {
	// Special treatment of P section throws away object indexing
	// (odd numbers in columns 64..71). First remember the number.
      if (sect == P)
	{
	  Pcurr = atoi(line_buffer+64);
	  if (Pnumber_.size() == 0 || Pnumber_[Pnumber_.size()-1] < Pcurr)
	    Pnumber_.push_back(Pcurr);
	  line_buffer[64] = 0;
	}
      sbufs[sect] += line_buffer;
      ++num_lines_[sect];
      DEBUG_ERROR_IF(num_lines_[sect] != line_number,
	       "Error in line numbers detected in IGES file (count vs. read line number): "
	       << num_lines_[sect] << " != " << line_number);
    }

    // Now we verify that the terminating section claims the same number of
    // lines that we counted for every section:

    sect = T;
    const char* cs = sbufs[sect].c_str();
    char numbuf[8];
    for (int i=0; i<4; i++) {
	strncpy(numbuf, cs + 1 + i*8, 8);
	numbuf[7] = 0;
	DEBUG_ERROR_IF(atoi(numbuf) != num_lines_[i],
		       "Terminate section line counts disagree with my counts.\n"
		       "Section count " << i << " is " << atoi(numbuf) <<
		       ", I counted " << num_lines_[i]);
    }

    // The start section data are not treated currently. All data there
    // are replaced by data from this converter (but it's just a comment
    // section, so it is not significant).

    // The global section is read into the header_ variable:
    readIGESheader(sbufs[G]);

    // Now, we need to scan the directory section
    // They are (as opposed to the global section) NOT stored in
    // a class variable for later use. afr: Now they are, too.
    // As we scan the D section we discover entities of type 128 and
    // type 126. They are immediately read into surf_ and crv_.
        // The P-numbers are also remembered.
    int sz = (int)sbufs[D].length();
    //cout << sz << endl << sbufs[D] << endl;
    int num_entries = sz/144;
    DEBUG_ERROR_IF (num_entries*144 != sz,
		"Directory section string has length != N*144,"
		"that is: a whole number times two lines");
    direntries_.resize(num_entries);
    const char* posD = sbufs[D].c_str();
    const char* posP0 = sbufs[P].c_str();
    start_of_P_section_ = posP0;
    const char* posP = posP0;
    //char pd = ',';
//     char rd = ';';
    // First we read all entities that may be included as part of
    // other entities (such as curve segments and surfaces, used for
    // composite curves and trimmed surfaces).
    // @@sbr We really should read all parts that are not created
    // using other entities.
    for (int i=0; i<num_entries; ++i) {
	direntries_[i] = readIGESdirentry(posD + i*144);
// 	std::cout << i << ' ' << direntries_[i].entity_type_number << ' '
// 	     << direntries_[i].param_data_start << ' '
// 	     << direntries_[i].line_count << std::endl;
// 	std::cout << (posP0 + 72*(direntries_[i].param_data_start-1)) << std::en//dl;
        posP = posP0 + 64*(direntries_[i].param_data_start-1);

	int entity_number = direntries_[i].entity_type_number;
	if (!supp_ent_.validEntity(entity_number))
	{
	    MESSAGE("Unknown entity-type (" << entity_number <<
		    ") in file! Object neglected.");
	}
	else if (entity_number == 128)
        {
          local_geom_.push_back(readIGESsurface(posP,
                                                direntries_[i].line_count));
	  local_colour_.push_back(direntries_[i].color);
          geom_id_.push_back(Pnumber_[i]);
          geom_used_.push_back(0);
        }
	else if (entity_number == 126)
        {
          shared_ptr<SplineCurve> crv =
	    readIGEScurve(posP, direntries_[i].line_count, i);
          local_geom_.push_back(crv);
	  local_colour_.push_back(direntries_[i].color);
          geom_id_.push_back(Pnumber_[i]);
          geom_used_.push_back(0);
	  pnumber_to_plane_normal_index_[Pnumber_[i]] = (int)plane_normal_.size()-1;
        }
	else if (entity_number == 110)
        {
          shared_ptr<BoundedCurve> crv =
	      readIGESline(posP, direntries_[i].line_count, i);
          local_geom_.push_back(crv);
	  local_colour_.push_back(direntries_[i].color);
          geom_id_.push_back(Pnumber_[i]);
          geom_used_.push_back(0);
 	  pnumber_to_plane_normal_index_[Pnumber_[i]] = (int)plane_normal_.size()-1;
	}
	else if (entity_number == 116)
        {
          shared_ptr<PointCloud3D> pt_cl =
	      readIGESpointCloud(posP, direntries_[i].line_count);
          local_geom_.push_back(pt_cl);
	  local_colour_.push_back(direntries_[i].color);
          geom_id_.push_back(Pnumber_[i]);
          geom_used_.push_back(0);
        }
	else if (entity_number == 123)
        {
          shared_ptr<PointCloud3D> pt_cl =
	      readIGESdirection(posP, direntries_[i].line_count);
          local_geom_.push_back(pt_cl);
	  local_colour_.push_back(direntries_[i].color);
          geom_id_.push_back(Pnumber_[i]);
          geom_used_.push_back(0);
        }
	else if (entity_number == 124)
        {
	    shared_ptr< CoordinateSystem<3> > cs
		= readIGEStransformation(posP, direntries_[i].line_count);
	    coordsystems_[Pnumber_[i]] = *cs;
	}
	else if  (entity_number == 108)
	{
	    // Planar surface
          shared_ptr<Plane> plane =
	      readIGESplane(posP, direntries_[i].line_count,
			    direntries_[i].form);
          local_geom_.push_back(plane);
	  local_colour_.push_back(direntries_[i].color);
          geom_id_.push_back(Pnumber_[i]);
          geom_used_.push_back(0);
	}
	else if (entity_number == 104)
	{
	  // Conic arc (parabola, ellipse, hyperbola)
	  shared_ptr<ElementaryCurve> conic_arc =
	      readIGESconicArc(posP, direntries_[i].line_count,
			       direntries_[i].form);
	  local_geom_.push_back(conic_arc);
 	  local_colour_.push_back(direntries_[i].color);
	  geom_id_.push_back(Pnumber_[i]);
	  geom_used_.push_back(0);
	}
    }
    
    // We scan the directory, looking for entities of type 100,
    // circular segment
    posP = posP0;
    for (int i=0; i<num_entries; ++i) {
	if (direntries_[i].entity_type_number == 100) {
	    posP = posP0 + 64*(direntries_[i].param_data_start-1);
	    shared_ptr<SplineCurve> crv =
		readIGEScircularsegment(posP, direntries_[i].line_count, i);
	    local_geom_.push_back(crv);
	    local_colour_.push_back(direntries_[i].color);
	    geom_id_.push_back(Pnumber_[i]);
	    geom_used_.push_back(0);
	}
    }

    // @@@ For the moment this does not handle general cases.
    // (Entities of type 'point' and 'connected point' not supported.)
    // We assume input entities have been represented as spline curves.
    // We scan the directory, looking for entities of type 102,
    // composite curve.
    posP = posP0;
    for (int i=0; i<num_entries; ++i) {
	if (direntries_[i].entity_type_number == 102) {
	    posP = posP0 + 64*(direntries_[i].param_data_start-1);
            vector<int> crv_vec;
	    readIGEScompositeCurve(posP, direntries_[i].line_count, i, crv_vec);
	    // Both 141 & 143 need modifications to handle composite curves.
            local_comp_curve_.push_back(crv_vec);
            comp_curve_id_.push_back(Pnumber_[i]);
	}
    }

    // We scan the directory, looking for entities of type 106 (form 12),
    // linear path entity.
    posP = posP0;
    for (int i=0; i<num_entries; ++i) {
	if (direntries_[i].entity_type_number == 106 &&
	    direntries_[i].form == 12){
	    posP = posP0 + 64*(direntries_[i].param_data_start-1);
	    shared_ptr<SplineCurve> crv =
		readIGESlinearPath(posP, direntries_[i].line_count,
				   direntries_[i].form);
	    local_geom_.push_back(crv);
	    local_colour_.push_back(direntries_[i].color);
	    geom_id_.push_back(Pnumber_[i]);
	    geom_used_.push_back(0);
	}
    }

    // We scan the directory, looking for entities of type 118,
    // ruled surface
    posP = posP0;
    for (int i=0; i<num_entries; ++i) {
	if (direntries_[i].entity_type_number == 118) {
	    posP = posP0 + 64*(direntries_[i].param_data_start-1);
	    local_geom_.push_back
                (readIGESruledSurface(posP, direntries_[i].line_count,
                                      direntries_[i].form));
	    geom_id_.push_back(Pnumber_[i]);
	    geom_used_.push_back(0);
	    local_colour_.push_back(direntries_[i].color);
	}
    }

    // We scan the directory, looking for entities of type 120,
    // surface of revolution
    posP = posP0;
    for (int i=0; i<num_entries; ++i) {
	if (direntries_[i].entity_type_number == 120) {
	    posP = posP0 + 64*(direntries_[i].param_data_start-1);
	    local_geom_.push_back
	      (readIGESsurfOfRevolution(posP, direntries_[i].line_count));
	    geom_id_.push_back(Pnumber_[i]);
	    geom_used_.push_back(0);
	    local_colour_.push_back(direntries_[i].color);
	}
    }


    // We scan the directory, looking for entities of type 122,
    // tabulated cylinder.
    posP = posP0;
    for (int i=0; i<num_entries; ++i) {
	if (direntries_[i].entity_type_number == 122) {
	    posP = posP0 + 64*(direntries_[i].param_data_start-1);
	    local_geom_.push_back
	      (readIGEStabulatedCylinder(posP, direntries_[i].line_count));
	    geom_id_.push_back(Pnumber_[i]);
	    geom_used_.push_back(0);
	    local_colour_.push_back(direntries_[i].color);
	}
    }

    // Now, we need to scan the directory section again.
    // This time we look for entities of type 141, boundary entity.
    posP = posP0;
    for (int i=0; i<num_entries; ++i) {
	if (direntries_[i].entity_type_number == 141) {
	    posP = posP0 + 64*(direntries_[i].param_data_start-1);
            vector<shared_ptr<CurveOnSurface> > crv_vec;
	    try {
	      readIGESboundary(posP,direntries_[i].line_count,crv_vec);
	    } catch (...) {
		MESSAGE("Failed reading entity number 141! Trying to continue anyway.");
	    }
            local_loop_.push_back(crv_vec);
            loop_id_.push_back(Pnumber_[i]);
	}
    }
    // We scan the directory section again.
    // We look for entities of type 143, boundary surface.
    posP = posP0;
    for (int i=0; i<num_entries; ++i) {
	if (direntries_[i].entity_type_number == 143) {
	    posP = posP0 + 64*(direntries_[i].param_data_start-1);
	    shared_ptr<BoundedSurface> bd_sf;
	    try {
	      bd_sf = readIGESboundedSurf(posP,
					  direntries_[i].line_count);
	    } catch (...) {
		MESSAGE("Failed reading entity number 143! Trying to continue anyway.");
	    }
	    local_geom_.push_back(bd_sf);
	    geom_id_.push_back(Pnumber_[i]);
	    geom_used_.push_back(0);
	    local_colour_.push_back(direntries_[i].color);
	}
    }
    
    // Now, we need to scan the directory section again.
    // This time we look for entities of type 142, curve on parametric surface.
    // @@ Warning! This entity is currently saved in local_loop_, although it
    // may not be a loop at all! Change variable name?
    posP = posP0;
    for (int i=0; i<num_entries; ++i) {
	if (direntries_[i].entity_type_number == 142) {
	    posP = posP0 + 64*(direntries_[i].param_data_start-1);
	    vector<shared_ptr<CurveOnSurface> > crv_vec;
	    try {
		readIGEScurveOnSurf(posP, direntries_[i].line_count, crv_vec);
	    } catch (...) {
		// @@sbr Possibly throw here.
		MESSAGE("Failed reading entity number 142! "
			"Trying to continue anyway.");
	    }
	    local_loop_.push_back(crv_vec);
	    loop_id_.push_back(Pnumber_[i]);
	}
    }

    // We scan the directory section again.
    // We look for entities of type 144, trimmed (parametric) surface.
    posP = posP0;
    for (int i=0; i<num_entries; ++i) {
	if (direntries_[i].entity_type_number == 144) {
	    posP = posP0 + 64*(direntries_[i].param_data_start-1);
 	    shared_ptr<BoundedSurface> bd_sf;
	    try {
 		bd_sf = readIGEStrimmedSurf(posP, direntries_[i].line_count);
 	    } catch (...) {
	      MESSAGE("Failed reading trimmed sf (144)! Continuing anyway.");
	    }
 	    local_geom_.push_back(bd_sf);
            geom_id_.push_back(Pnumber_[i]);
            geom_used_.push_back(0);
	    local_colour_.push_back(direntries_[i].color);
	}
    }

    // We scan the directory section again.
    // We look for entities of type 190, plane surface.
    posP = posP0;
    for (int i=0; i<num_entries; ++i) {
      if (direntries_[i].entity_type_number == 190) {
	posP = posP0 + 64*(direntries_[i].param_data_start-1);
	vector<double> colour;
	string cname;
	shared_ptr<Plane> plane_sf =
	  readIGESplaneSurface(posP, direntries_[i].line_count,
			       direntries_[i].form);
	local_geom_.push_back(plane_sf);
	local_colour_.push_back(direntries_[i].color);
	geom_id_.push_back(Pnumber_[i]);
	geom_used_.push_back(0);
      }
    }

    // We scan the directory section again.
    // We look for entities of type 192, right circular cylindrical
    // surfaces.
    posP = posP0;
    for (int i=0; i<num_entries; ++i) {
	if (direntries_[i].entity_type_number == 192) {
	    posP = posP0 + 64*(direntries_[i].param_data_start-1);
	    vector<double> colour;
	    string cname;
	  // Cylinder
	    shared_ptr<Cylinder> cylinder =
		readIGESrightCircularCylindricalSurface(posP,
							direntries_[i].line_count,
							direntries_[i].form);
	    local_geom_.push_back(cylinder);
	    local_colour_.push_back(direntries_[i].color);
	    geom_id_.push_back(Pnumber_[i]);
	    geom_used_.push_back(0);
	}
    }

    // We scan the directory section again.
    // We look for entities of type 314, colour description.
    posP = posP0;
    for (int i=0; i<num_entries; ++i) {
	if (direntries_[i].entity_type_number == 314) {
	    posP = posP0 + 64*(direntries_[i].param_data_start-1);
	    vector<double> colour;
	    string cname;
	    readIGEScolour(posP, direntries_[i].line_count, colour, cname);
	    colour_objects_.push_back(colour);
	    colour_name_.push_back(cname);
	    colour_id_.push_back(Pnumber_[i]);
	}
    }

    // We scan the directory section again.
    // This time we look for entities of type 402, group associativity.
    posP = posP0;
    for (int i=0; i<num_entries; ++i) {
        if (direntries_[i].entity_type_number == 402 &&
	    direntries_[i].form == 7){
	    posP = posP0 + 64*(direntries_[i].param_data_start-1);
	    group_.push_back(readIGESgroupAssembly(posP,
						   direntries_[i].line_count));
	    ftTangPriority prio_type = ftNoType;
	    if (direntries_[i].entity_number == 1)
	      prio_type = ftMaster;
	    else if (direntries_[i].entity_number == 2)
	      prio_type = ftSlave;
	    group_[group_.size()-1].setType(prio_type);
	}
    }

    // We scan the directory section again.
    // This time we look for entities of type 502, vertex list.
    posP = posP0;
    for (int i=0; i<num_entries; ++i) {
      if (direntries_[i].entity_type_number == 502) {
	posP = posP0 + 64*(direntries_[i].param_data_start-1);
	MESSAGE("Entity number 502 (vertex list) soon to be supported!"
		"Object currently neglected.");
      }
    }

    // We scan the directory section again.
    // This time we look for entities of type 504, edge.
    posP = posP0;
    for (int i=0; i<num_entries; ++i) {
      if (direntries_[i].entity_type_number == 504) {
	posP = posP0 + 64*(direntries_[i].param_data_start-1);
	MESSAGE("Entity number 504 (edge list) soon to be supported! "
		"Object currently neglected.");
      }
    }

    // We scan the directory section again.
    // This time we look for entities of type 508, loop.
    posP = posP0;
    for (int i=0; i<num_entries; ++i) {
      if (direntries_[i].entity_type_number == 508) {
	posP = posP0 + 64*(direntries_[i].param_data_start-1);
	MESSAGE("Entity number 508 (loop) soon to be supported! "
		"Object currently neglected.");
      }
    }

    // We scan the directory section again.
    // This time we look for entities of type 510, vertex list.
    posP = posP0;
    for (int i=0; i<num_entries; ++i) {
      if (direntries_[i].entity_type_number == 510) {
	posP = posP0 + 64*(direntries_[i].param_data_start-1);
	MESSAGE("Entity number 510 (face) soon to be supported! "
		"Object currently neglected.");
      }
    }

    // Collect all individual entities
    for (size_t i = 0; i < local_geom_.size(); ++i) {
        if (!geom_used_[i]) {
            geom_.push_back(local_geom_[i]);
            // We extract colour information.
            if (local_colour_[i] >= 0)
                colour_.push_back(transformFromIGEScolour(local_colour_[i]));
            else {
                size_t j;
                for (j = 0; j < colour_id_.size(); ++j)
                    if (colour_id_[j] == -local_colour_[i])
                        break;
                if (j == colour_id_.size()) {
                    THROW("Colour number " << j << " not found!");
                }
                else {
                    colour_.push_back(colour_objects_[j]);
                }
                // std::cout << local_geom_[i]->instanceType() << std::endl;
            }
        }
    }
    filled_with_data_ = true;
}





//-----------------------------------------------------------------------------
void IGESconverter::writeIGES(ostream& os)
//-----------------------------------------------------------------------------
{
    DEBUG_ERROR_IF(!filled_with_data_,
	       "Cannot write. I contain no data.\n" 
	       "Please call one of the read functions first");
    if (!filled_with_data_) return;

    // A line's worth of spaces. Used for padding strings.

    // An IGES file consists of five sections.
    // First is the start section, which consists of a comment:
    string comment
	= " Sintef Applied Mathematics IGES converter 1.1 IGES version 5.3";
    pad(comment,72);
    // The comment is assumed to be one line
    int num_lines[5];
    num_lines[S] = 1;
    writeSingleIGESLine(os, comment.c_str(), 1, S);

    // Then the global section follows
    string sec;
    // Write header section into sec
    writeIGESheader(sec);
    pad(sec,72);
    num_lines[G] = (int)sec.length()/72;
    for (int i=0; i<num_lines[G]; ++i)
	writeSingleIGESLine(os, sec.c_str() + 72*i, i+1, G);

    // Next is the directory section. BUT we need the line numbers from
    // the parameter section for each entity, so we must create that one
    // first.
    string parsect = "";
    vector<IGESdirentry> ent;
    /* ent.reserve(geom_.size());  // Can be larger due to bounded surfaces. */
    writeIGESparsect(parsect, ent);
    // Write dir section into sec
    writeIGESdirectory(sec, ent);
    // Write out D section
    for (size_t i=0; i<2*ent.size(); ++i)
	writeSingleIGESLine(os, sec.c_str() + 72*(int)i, (int)i+1, D);
    num_lines[D] = 2*(int)ent.size();
    // Write out P section
    num_lines[P] = (int)parsect.length()/64;
    char line72[73];
    int geom_num = 0;
    for (int i=0; i<num_lines[P]; ++i) {
	strncpy(line72, parsect.c_str() + 64*i, 64);
	if (geom_num < (int)ent.size()-1 &&
	    (i+1 >= ent[geom_num+1].param_data_start))
	    ++geom_num;
	sprintf(line72+64, "%8i", geom_num*2 + 1);
	writeSingleIGESLine(os, line72, i+1, P);
    }
    sprintf(line72, "S%7iG%7iD%7iP%7i", num_lines[S], num_lines[G],
	    num_lines[D], num_lines[P]);
    for (int i=32; i<72; ++i)
	line72[i] = ' ';
    writeSingleIGESLine(os, line72, 1, T);
}


//-----------------------------------------------------------------------------
void IGESconverter::readsislsrfs(istream& is)
//-----------------------------------------------------------------------------
{
    int nmb_sfs = 0;
    int value;             /* type of B-Spline and numbers */
    int kcuopen1;          /* Open/closed flag for surf  */
    int kcuopen2;          /* Open/closed flag for surf  */
    int kdim;              /* dimension of target space */
    int kdim1;             /* dimension + 1             */
    int kk1;               /* Order of curve or in 1st parameter direction of surface */
    int kk2;               /* Order in 2nd parameter direction of surface */
    int kn1;               /* # of vertices of curve or in 1st parameter direction of surf */
    int kn2;               /* # of vertices in 2nd parameter direction of surface */
    vector<double> coefs; //double *scoef=NULL;    /* Vertices of curve or surface */
    double dummy;          /* trashbin for fourth value in vertices */
    double matrix[16];     /* instance matrix */
    int ki, kh;            /* Counters */
    vector<double> knots_u, knots_v;
    double knot;
    shared_ptr<SplineSurface> sf;
    int rat = 0;

    SkipSislComments(is);
    is >> nmb_sfs;

    for (kh = 0; kh < nmb_sfs; ++kh)
      {
	knots_u.clear();
	knots_v.clear();

	SkipSislComments(is);
	is >> value; // Type of object, not used. Could of course verify that it matches that of Spline Sf.

	// orders of parameter directions
	SkipSislComments(is);
	is >> kk1;
	SkipSislComments(is);
	is >> kk2;

	// number of vertices
	SkipSislComments(is);
	is >> kn1;
	SkipSislComments(is);
	is >> kn2;

	// dimension
	SkipSislComments(is);
	is >> kdim;
	kdim1 = kdim;

	// surface open/closed
	SkipSislComments(is);
	is >> kcuopen1 >> kcuopen2;

	// rational B-Spline or not
	SkipSislComments(is);
	is >> rat;
	if (rat) kdim1++;

	// Allocate space for first knot vector
	knots_u.reserve(kk1+kn1);

	// Read in first knot vector.
	SkipSislComments(is);
	for (ki = 0; ki < kn1+kk1; ++ki)
	  {
	    is >> knot;
	    knots_u.push_back(knot);
	  }
	// Allocate space for second knot vector
	knots_v.reserve(kk2+kn2);

	// Read in second knot vector.
	SkipSislComments(is);
	for (ki = 0; ki < kn2+kk2; ++ki)
	  {
	    is >> knot;
	    knots_v.push_back(knot);
	  }

	// Allocate space for control vertices
	coefs.reserve(kn1*kn2*kdim1);
	for (ki = 0; ki < kn1*kn2*kdim1; ++ki)
	  coefs.push_back(0.0);

	// Read in control vertices.
	SkipSislComments(is);
	is >> coefs[0];
	if (kdim > 1) is >> coefs[1]; //fscanf(fp,"%lf",&coefs[1]);
	if (kdim > 2) is >> coefs[2]; //fscanf(fp,"%lf",&coefs[2]);
	is >> dummy;
	if (rat) coefs[kdim] = dummy;
	for (ki = kdim1; ki < kn1*kn2*kdim1; ki += kdim1)
	  {
	    is >> coefs[ki];
	    if (kdim > 1) is >> coefs[ki+1]; //fscanf(fp,"%lf",&coefs[ki+1]);
	    if (kdim > 2) is >> coefs[ki+2]; //fscanf(fp,"%lf",&coefs[ki+2]);
	    is >> dummy;
	    if (rat) coefs[ki+kdim] = dummy;
	  }

	// transformation data
	SkipSislComments(is);
	for (ki = 0; ki < (kdim+1)*(kdim+1); ki ++)
	  is >> matrix[ki];
	// assuming matrix is unit matrix, thus not being used.

	// Create surface.
	// @@sbr Not handling closed or periodic sfs.
	sf = shared_ptr<SplineSurface>(new SplineSurface(kn1, kn2, kk1, kk2, knots_u.begin(),
							 knots_v.begin(),
							 coefs.begin(), 3, rat!=0));

	if (kcuopen1 == -1 || kcuopen2 == -1)
	  {
	    // We want spline space to be open.
	    shared_ptr<SplineSurface> sub_sf(sf->subSurface(sf->startparam_u(), sf->startparam_v(),
							    sf->endparam_u(), sf->endparam_v()));
	    sf = sub_sf;
// 	    MESSAGE("Surface is periodic, currently not supported!");
	  }
// 	if (kcuopen1 == -1 || kcuopen2 == -1) {
// 	    MESSAGE("Surface is periodic, currently not supported!");
// 	}

	shared_ptr<GeomObject> obj(Factory::createObject(Class_SplineSurface));
	obj = sf;
	geom_.push_back(obj);
	vector<double> col;
	colour_.push_back(col);
	SkipSislComments(is);
      }

    filled_with_data_ = true;
}


//-----------------------------------------------------------------------------
void IGESconverter::readsislcrvs(istream& is)
//-----------------------------------------------------------------------------
{
    int nmb_cvs = 0;
    int value;             /* type of B-Spline and numbers */
    int kcuopen;          /* Open/closed flag for surf  */
    int kdim;              /* dimension of target space */
    int kdim1;             /* dimension + 1             */
    int kk;               /* Order of curve */
    int kn;               /* # of vertices of curve */
    vector<double> coefs; //double *scoef=NULL;    /* Vertices of curve */
    double dummy;          /* trashbin for fourth value in vertices */
    double matrix[16];     /* instance matrix */
    int ki, kh;            /* Counters */
    vector<double> knots;
    double knot;
    shared_ptr<SplineCurve> cv;
    int rat = 0;

    SkipSislComments(is);
    is >> nmb_cvs;

    for (kh = 0; kh < nmb_cvs; ++kh)
      {
	knots.clear();

	SkipSislComments(is);
	is >> value; // Type of object, not used. Could of course verify that it matches that of Spline Cv.

	// order of cv.
	SkipSislComments(is);
	is >> kk;

	// number of vertices
	SkipSislComments(is);
	is >> kn;

	// dimension
	SkipSislComments(is);
	is >> kdim;
	kdim1 = kdim;

	// curve open/closed
	SkipSislComments(is);
	is >> kcuopen;

	// rational B-Spline or not
	SkipSislComments(is);
	is >> rat;
	if (rat) kdim1++;

	// Allocate space for knot vector
	knots.reserve(kk+kn);

	// Read in knot vector.
	SkipSislComments(is);
	for (ki = 0; ki < kn+kk; ++ki)
	  {
	    is >> knot;
	    knots.push_back(knot);
	  }

	// Allocate space for control vertices
	coefs.reserve(kn*kdim1);
	for (ki = 0; ki < kn*kdim1; ++ki)
	  coefs.push_back(0.0);

	// Read in control vertices.
	SkipSislComments(is);
	is >> coefs[0];
	if (kdim > 1) is >> coefs[1]; //fscanf(fp,"%lf",&coefs[1]);
	if (kdim > 2) is >> coefs[2]; //fscanf(fp,"%lf",&coefs[2]);
	is >> dummy;
	if (rat) coefs[kdim] = dummy;
	for (ki = kdim1; ki < kn*kdim1; ki += kdim1)
	  {
	    is >> coefs[ki];
	    if (kdim > 1) is >> coefs[ki+1]; //fscanf(fp,"%lf",&coefs[ki+1]);
	    if (kdim > 2) is >> coefs[ki+2]; //fscanf(fp,"%lf",&coefs[ki+2]);
	    is >> dummy;
	    if (rat) coefs[ki+kdim] = dummy;
	  }

	// transformation data
	SkipSislComments(is);
	for (ki = 0; ki < (kdim+1)*(kdim+1); ki ++)
	  is >> matrix[ki];
	// assuming matrix is unit matrix, thus not being used.

	// Create surface.
	// @@sbr Not handling closed or periodic cvs.
	cv = shared_ptr<SplineCurve>(new SplineCurve(kn, kk, knots.begin(),
						     coefs.begin(), 3, rat!=0));

	if (kcuopen == -1)
	  {
	    // We want spline space to be open.
	    shared_ptr<SplineCurve> sub_cv(cv->subCurve(cv->startparam(), cv->endparam()));
	    cv = sub_cv;
// 	    MESSAGE("Curve is periodic, currently not supported!");
	  }
// 	if (kcuopen == -1) {
// 	    MESSAGE("Curve is periodic, currently not supported!");
// 	}

	shared_ptr<GeomObject> obj(Factory::createObject(Class_SplineCurve));
	obj = cv;
	geom_.push_back(obj);
	vector<double> col;
	colour_.push_back(col);
	SkipSislComments(is);
      }

    filled_with_data_ = true;
}


//-----------------------------------------------------------------------------
void IGESconverter::readgo(istream& is)
//-----------------------------------------------------------------------------
{
    ObjectHeader header;

    while (!is.eof()) {
	header.read(is);
	//Read(is, header);
	shared_ptr<GeomObject> obj(Factory::createObject(header.classType()));
	//Read(is, *obj);
	obj->read(is);
	geom_.push_back(obj);
	vector<double> col;
	if (header.auxdataSize() == 4) {
	    for (int ki = 0; ki < 4; ++ki)
		col.push_back((double)header.auxdata(ki)*100/255);
	}
	colour_.push_back(col);
	Utils::eatwhite(is);
	//SkipComments(is);  
    }

    filled_with_data_ = true;
}


//-----------------------------------------------------------------------------
void IGESconverter::writego(ostream& os)
//-----------------------------------------------------------------------------
{
     DEBUG_ERROR_IF(!filled_with_data_,
	       "Cannot write. I contain no data.\n" 
	       "Please call one of the read functions first");
    if (!filled_with_data_) return;
    for (size_t i = 0; i < geom_.size(); ++i) {
	// We introduce local pointer as we do not handle trimmed surfaces.
	shared_ptr<GeomObject> object = geom_[i];
	if (object.get() == 0) {
	    MESSAGE("Object missing! Moving on to next!");
	    continue;
	}
	int major = 1;
	int minor = 0;
	ClassType class_type = geom_[i]->instanceType();
	// @@@ Local hack as writing of trimmed surface is yet to be implemented.
	// afr: Removed it, as trimmed surfaces do read and write now, as long as
	// the underlying geometries are splines.
//  	if (class_type == 210) {
//  	    MESSAGE("Bounded surface is written as an untrimmed surface!");
//  	    class_type = Class_SplineSurface;
//  	    object = dynamic_pointer_cast<BoundedSurface, GeomObject>(object)->
//  		underlyingSurface();
//  	}
	vector<int> colour;
	if ((i < colour_.size()) && (int(colour_[i].size()) == 3)) {
	    for (int j = 0; j < 3; ++j)
		colour.push_back((int)(255.0*colour_[i][j]/100.0));
	    colour.push_back(255);
	} else {
	    // Default colour is blue
	    colour.push_back(0);
	    colour.push_back(0);
	    colour.push_back(255);
	    colour.push_back(255);
	}
	ObjectHeader local_header(class_type, major, minor, colour);
	local_header.write(os);
	object->write(os);
    }
}


//-----------------------------------------------------------------------------
void IGESconverter::addGeom(shared_ptr<GeomObject> sp)
//-----------------------------------------------------------------------------
{
    geom_.push_back(sp);
    vector<double> dummy_vec;
    colour_.push_back(dummy_vec);
    filled_with_data_ = true; // Needed when doing push_back onto a new object
}




//---------------------- Private members --------------------------




//-----------------------------------------------------------------------------
void IGESconverter::readIGESheader(string g)
//-----------------------------------------------------------------------------
{
    // Find parameter delimiter
    const char* whereami = g.c_str();
    string it;
    it = readIGESstring(whereami, ',');
    char pd = ',';
    if (it.size() != 0) {
	pd = it[0];
    }
    header_.pardel = pd;
    skipDelimiter(whereami, pd);

    // Find record delimiter
    it = readIGESstring(whereami, pd);
    char rd = ';';
    if (it.size() != 0) {
	rd = it[0];
    }
    header_.recdel = rd;
    skipDelimiter(whereami, pd);

    //    cout << pd << "  " << rd << endl;

    // Now whereami points to the third parameter
    header_.sending_system_prodid = readIGESstring(whereami, pd);
    skipDelimiter(whereami, pd);

    header_.original_filename = readIGESstring(whereami, pd);
    skipDelimiter(whereami, pd);

    header_.creator_system_prodid = readIGESstring(whereami, pd);
    skipDelimiter(whereami, pd);

    header_.preprocessor_version_number = readIGESstring(whereami, pd);
    skipDelimiter(whereami, pd);

    header_.integer_bits = readIGESint(whereami, pd, rd);
    skipDelimiter(whereami, pd);

    header_.single_prec_magn = readIGESint(whereami, pd, rd);
    skipDelimiter(whereami, pd);

    header_.single_prec_sign = readIGESint(whereami, pd, rd);
    skipDelimiter(whereami, pd);

    header_.double_prec_magn = readIGESint(whereami, pd, rd);
    skipDelimiter(whereami, pd);

    header_.double_prec_sign = readIGESint(whereami, pd, rd);
    skipDelimiter(whereami, pd);

    header_.receiving_system_prodid = readIGESstring(whereami, pd);
    skipDelimiter(whereami, pd);

    header_.model_space_scale = readIGESdouble(whereami, pd, rd);
    skipDelimiter(whereami, pd);

    header_.unit_flag = readIGESint(whereami, pd, rd);
    skipDelimiter(whereami, pd);

    header_.unit_description = readIGESstring(whereami, pd);
    skipDelimiter(whereami, pd);

    header_.max_weight_grads = readIGESint(whereami, pd, rd);
    skipDelimiter(whereami, pd);

    header_.max_line_width = readIGESdouble(whereami, pd, rd);
    skipDelimiter(whereami, pd);

    header_.timestamp_file_changed = readIGESstring(whereami, pd);
    skipDelimiter(whereami, pd);

    header_.min_resolution = readIGESdouble(whereami, pd, rd);
    skipDelimiter(whereami, pd);

    header_.max_coordinate = readIGESdouble(whereami, pd, rd);
    skipDelimiter(whereami, pd);

    header_.author = readIGESstring(whereami, pd);
    skipDelimiter(whereami, pd);

    header_.organisation = readIGESstring(whereami, pd);
    skipDelimiter(whereami, pd);

    header_.IGES_version = readIGESint(whereami, pd, rd);
    skipDelimiter(whereami, pd);

    header_.drafting_standard = readIGESint(whereami, pd, rd);


    if (*whereami != rd) {
	skipDelimiter(whereami, pd);
	header_.timestamp_model_changed = readIGESstring(whereami, pd);
	if (*whereami != rd) {
	  skipDelimiter(whereami, pd);
	  header_.application_protocol_description = readIGESstring(whereami, pd);
	  header_.num_parameters = 26;
	}
	else
	  header_.num_parameters = 25;
    } else {
	header_.timestamp_model_changed = header_.timestamp_file_changed;
	header_.num_parameters = 24;
    }
    skipDelimiter(whereami, rd);
}


//-----------------------------------------------------------------------------
void IGESconverter::writeIGESheader(string& g)
//-----------------------------------------------------------------------------
{
    char pd = header_.pardel;
    char rd = header_.recdel;

    g = "";
    g += "1H";
    g += pd;
    g += pd;
    g += "1H";
    g += rd;
    g += pd;

    g += writeIGESstring(header_.sending_system_prodid);
    g += pd;

    g += writeIGESstring(header_.original_filename);
    g += pd;

    g += writeIGESstring(header_.creator_system_prodid);
    g += pd;

    g += writeIGESstring(header_.preprocessor_version_number);
    g += pd;

    g += writeIGESint(header_.integer_bits);
    g += pd;

    g += writeIGESint(header_.single_prec_magn);
    g += pd;

    g += writeIGESint(header_.single_prec_sign);
    g += pd;

    g += writeIGESint(header_.double_prec_magn);
    g += pd;

    g += writeIGESint(header_.double_prec_sign);
    g += pd;

    g += writeIGESstring(header_.receiving_system_prodid);
    g += pd;

    g += writeIGESdouble(header_.model_space_scale);
    g += pd;

    g += writeIGESint(header_.unit_flag);
    g += pd;

    g += writeIGESstring(header_.unit_description);
    g += pd;

    g += writeIGESint(header_.max_weight_grads);
    g += pd;

    g += writeIGESdouble(header_.max_line_width);
    g += pd;

    g += writeIGESstring(header_.timestamp_file_changed);
    g += pd;

    g += writeIGESdouble(header_.min_resolution);
    g += pd;

    g += writeIGESdouble(header_.max_coordinate);
    g += pd;

    g += writeIGESstring(header_.author);
    g += pd;

    g += writeIGESstring(header_.organisation);
    g += pd;

    g += writeIGESint(header_.IGES_version);
    g += pd;

    g += writeIGESint(header_.drafting_standard);
    g += pd;

    g += writeIGESstring(header_.timestamp_model_changed);
    g += pd;

    g += writeIGESstring(header_.application_protocol_description);
    g += rd;
}

//-----------------------------------------------------------------------------
void IGESconverter::writeIGESdirectory(string& d, 
				       const vector<IGESdirentry>& dirent)
//-----------------------------------------------------------------------------
{
    d = "";
    char buffer[9];
    for (size_t i=0; i<dirent.size(); ++i) {
	const IGESdirentry& de = dirent[i];
	sprintf(buffer, "%8i", de.entity_type_number);
	d += buffer;
	sprintf(buffer, "%8i", de.param_data_start);
	d += buffer;
	sprintf(buffer, "%8i", de.structure);
	d += buffer;
	sprintf(buffer, "%8i", de.line_font_pattern);
	d += buffer;
	sprintf(buffer, "%8i", de.level);
	d += buffer;
	sprintf(buffer, "%8i", de.view);
	d += buffer;
	sprintf(buffer, "%8i", de.trans_matrix);
	d += buffer;	
	sprintf(buffer, "%8i", de.label_display);
	d += buffer;
	d += de.status;
	sprintf(buffer, "%8i", de.entity_type_number);
	d += buffer;	
	sprintf(buffer, "%8f", de.line_weight);
	d += buffer;
	sprintf(buffer, "%8i", de.color);
	d += buffer;
	sprintf(buffer, "%8i", de.line_count);
	d += buffer;
	sprintf(buffer, "%8i", de.form);
	d += buffer;
	d += "                "; // Empty reserved fields
	d += de.entity_label;
	sprintf(buffer, "%8i", de.entity_number);
	d += buffer;
    }




    
//     int sz = sbufs[D].length();
//     //cout << sz << endl << sbufs[D] << endl;
//     int num_entries = sz/144;
//     DEBUG_ERROR_IF (num_entries*144 != sz,
//                 "Directory section string has length != N*144");
//     vector<IGESdirentry> direntries_(num_entries);
//     vector<int> surf_in_dir;
//     char* posD = sbufs[D].c_str();
//     char* posP0 = sbufs[P].c_str();
//     const char* posP = posP0;
//     num_surfs_ = 0;
//     for (int i=0; i<num_entries; ++i) {
// 	direntries_[i] = readIGESdirentry(posD + i*144);
// 	//cout << i << ' ' << direntries_[i].entity_type_number << ' '
// 	//     << direntries_[i].param_data_start << ' '
// 	//     << direntries_[i].line_count << endl;
// 	//cout << (posP0 + 72*(direntries_[i].param_data_start-1)) << endl;
// 	if (direntries_[i].entity_type_number == 128) {
// 	    surf_in_dir.push_back(i);
// 	    posP = posP0 + 64*(direntries_[i].param_data_start-1);
// 	    surfs_.push_back(readIGESsurface(posP, direntries_[i].line_count));
// 	    ++num_surfs_;
// 	}
//     }
    

}


//-----------------------------------------------------------------------------
void IGESconverter::writeIGESparsect(string& g, vector<IGESdirentry>& dirent)
//-----------------------------------------------------------------------------
{
//     char pd = header_.pardel;
//     char rd = header_.recdel;
    g = "";

        // First write all geometries (curves and surfaces and
        // bounded surfaces).
    int Pcurr = -1;
    // But before that we need to know about colours.
    vector<vector<double> > unique_colours = uniqueColours(colour_);
    for (size_t i = 0; i < unique_colours.size(); ++i) {
	writeIGEScolour(unique_colours[i], g, dirent, Pcurr);
    }
    for (size_t i=0; i<geom_.size(); ++i) {
	// We must find index of colour_[i].
	int col = 0;
	if (colour_[i].size() != 0) { // We must transform colour information.
	    size_t j;
	    for (j = 0; j < unique_colours.size(); ++j) {
		vector<double> inters;
		set_intersection(colour_[i].begin(), colour_[i].begin() + 3,
				 unique_colours[j].begin(), unique_colours[j].begin() + 3,
				 back_inserter(inters));
		if (inters.size() == 3)
		    break;
	    }
	    ASSERT(j < unique_colours.size());
	    col = -(2*int(j) + 1); // @@ Using the fact that colour info is placed first in iges-file...
	}
	if (geom_[i]->instanceType() == Class_SplineSurface)
	    writeIGESsurface(dynamic_cast<SplineSurface*>(geom_[i].get()),
			     col, g, dirent, Pcurr);
	else if (geom_[i]->instanceType() == Class_SplineCurve)
	    writeIGEScurve(dynamic_cast<SplineCurve*>(geom_[i].get()),
			   col, g, dirent, Pcurr);
	else if (geom_[i]->instanceType() == Class_BoundedSurface)
	    writeIGESboundedSurf(dynamic_cast<BoundedSurface*>(geom_[i].get()),
				 col,g, dirent, Pcurr);
    }
}

//-----------------------------------------------------------------------------
void IGESconverter::writeIGEScolour(const vector<double>& colour, string& g,
				    vector<IGESdirentry>& dirent, int& Pcurr)
//-----------------------------------------------------------------------------
{
    ASSERT(colour.size() > 2); // Even if alpha information is present, we only use the first 3 values.

    char pd = header_.pardel;
    char rd = header_.recdel;

    IGESdirentry curr_dirent;
    curr_dirent.entity_type_number = 314;
    curr_dirent.structure = 0;
    curr_dirent.line_font_pattern = 0;
    curr_dirent.level = 0;
    curr_dirent.view = 0;
    curr_dirent.trans_matrix = 0;
    curr_dirent.label_display = 0;
    ostringstream os;
    os << "00000000";
    curr_dirent.status = os.str();
    curr_dirent.line_weight = 0;
    curr_dirent.color = 0;
    curr_dirent.form = 0;
    curr_dirent.entity_label = "  colour";
    curr_dirent.entity_number = 0;

    int sz = (int)g.length();
    int ln = sz/64;
    DEBUG_ERROR_IF (ln*64 != sz, "String length not a multiple of 64.");

    curr_dirent.param_data_start = ln+1;

    char buffer[100]; // Enough digits for most numbers, I hope.

    g += "314";
    g += pd;
    sprintf(buffer, "%.5G", colour[0]);
    g += buffer;
    g += pd;
    sprintf(buffer, "%.5G", colour[1]);
    g += buffer;
    g += pd;
    sprintf(buffer, "%.5G", colour[2]);
    g += buffer;
    g += rd;
    // We don't bother to include any name. Takne from where?
    pad(g);

    curr_dirent.line_count = (int)g.size()/64 - curr_dirent.param_data_start + 1;
    dirent.push_back(curr_dirent);
    Pcurr+=2;
}



//-----------------------------------------------------------------------------
void IGESconverter::writeIGESsurface(SplineSurface* surf, int colour, string& g,
                                     vector<IGESdirentry>& dirent, int& Pcurr,
				     int dependency)
//-----------------------------------------------------------------------------
{
    char pd = header_.pardel;
    char rd = header_.recdel;

    IGESdirentry curr_dirent;
    curr_dirent.entity_type_number = 128;
    curr_dirent.structure = 0;
    curr_dirent.line_font_pattern = 1;
    curr_dirent.level = 0;
    curr_dirent.view = 0;
    curr_dirent.trans_matrix = 0;
    curr_dirent.label_display = 0;
    ostringstream os;
    os << "000" << dependency << "0001";
    curr_dirent.status = os.str();
    curr_dirent.line_weight = 0;
    curr_dirent.color = colour;
    curr_dirent.form = 0;
    curr_dirent.entity_label = "   nurbs";
    curr_dirent.entity_number = 0;

    int sz = (int)g.length();
    int ln = sz/64;
    DEBUG_ERROR_IF (ln*64 != sz, "String length not a multiple of 64.");

    curr_dirent.param_data_start = ln+1;


    g += "128";
    g += pd;

    char buffer[100]; // Enough digits for most numbers, I hope.

    sprintf(buffer, "%i", surf->numCoefs_u() - 1);
    g += buffer;
    g += pd;

    sprintf(buffer, "%i", surf->numCoefs_v() - 1);
    g += buffer;
    g += pd;

    sprintf(buffer, "%i", surf->order_u() - 1);
    g += buffer;
    g += pd;

    sprintf(buffer, "%i", surf->order_v() - 1);
    g += buffer;
    g += pd;

//  	g += (surf->closedness_u() == GoOPEN) ? "0" : "1";
    g += "0";
    g += pd;

//  	g += (surf->closedness_v() == GoOPEN) ? "0" : "1";
    g += "0";
    g += pd;

//  	if (surf->splinetype() == GoRATIONAL_BSPLINE 
//  	    || surf->splinetype() == GoRATIONAL_BEZIER)
//  	    g += "0";
//  	else
//  	    g += "1";
    int polynomial = (surf->rational()) ? 0 : 1;
    sprintf(buffer, "%i", polynomial);
    g += buffer;
    g += pd;
	
//  	g += (surf->closedness_u() <= GoCLOSED) ? "0" : "1";
    g += "0";
    g += pd;

//  	g += (surf->closedness_v() <= GoCLOSED) ? "0" : "1";
    g += "0";
    g += pd;

	// Pad with spaces until a 64-char line is filled
    pad(g);

	// Knots
    std::vector<double>::const_iterator it;
    for (it=surf->basis_u().begin();
         it!=surf->basis_u().end(); ++it) {
      sprintf(buffer, "%.15G", *it);
      if ((strlen(buffer)+g.length()+1)/64 > g.length()/64) {
		// Pad it till EOL
        pad(g);
      }
      g += buffer;
      g += pd;
    }       
    for (it=surf->basis_v().begin();
         it!=surf->basis_v().end(); ++it) {
      sprintf(buffer, "%.15G", *it);
      if ((strlen(buffer)+g.length()+1)/64 > g.length()/64) {
		// Pad it till EOL
		pad(g);
      }   
      g += buffer;
      g += pd;
    }

	// Pad with spaces until a 64-char line is filled
    pad(g);

	// @afr: Rationals are not handled properly
	// For now, we write 1.0 for all weights.
    int n = surf->numCoefs_v()*surf->numCoefs_u();
    int n1 = surf->numCoefs_u();
    int n2 = surf->numCoefs_v();
    int j, k;
    if (surf->rational())
    {
      std::vector<double>::const_iterator co = surf->rcoefs_begin();
      for (j=0; j<n2; ++j) {
        for (int k=0; k<n1; ++k) {
          sprintf(buffer, "%.13G",co[(k+j*n1)*4+3]);
          if ((strlen(buffer)+g.length()+1)/64 > g.length()/64) {
		// Pad it till EOL
            pad(g);
          }
          g += buffer;
          g += pd;
        }
      }
      pad(g);

      for (j=0; j<n2; ++j) {
        for (k=0; k<n1; ++k) {
          sprintf(buffer, "%.13G",co[(k+j*n1)*4]/co[(k+j*n1)*4+3]);
          g += buffer;
          g += pd;
          sprintf(buffer, "%.13G",co[(k+j*n1)*4 + 1]/co[(k+j*n1)*4+3]);
          g += buffer;
          g += pd;
          sprintf(buffer, "%.13G",co[(k+j*n1)*4 + 2]/co[(k+j*n1)*4+3]);
          g += buffer;
          g += pd;
          pad(g);
        }
      }
    }
    else
    {
      for (j=0; j<n; ++j) {
        sprintf(buffer, "1.0");
        if ((strlen(buffer)+g.length()+1)/64 > g.length()/64) {
		// Pad it till EOL
          pad(g);
        }
        g += buffer;
        g += pd;
      }
      pad(g);

      std::vector<double>::const_iterator co = surf->coefs_begin();
      for (j=0; j<n2; ++j) {
        for (k=0; k<n1; ++k) {
          sprintf(buffer, "%.13G",co[(k+j*n1)*3]);
          g += buffer;
          g += pd;
          sprintf(buffer, "%.13G",co[(k+j*n1)*3 + 1]);
          g += buffer;
          g += pd;
          sprintf(buffer, "%.13G",co[(k+j*n1)*3 + 2]);
          g += buffer;
          g += pd;
          pad(g);
        }
      }
    }

// #ifdef _MSC_VER
//     const RectDomain& dom = static_cast<const RectDomain&>(surf->parameterDomain());
// #else
    const RectDomain& dom = surf->parameterDomain();
// #endif

	// Write u and v parameter ranges
    sprintf(buffer, "%.13G", dom.umin());
    g += buffer;
    g += pd;
    sprintf(buffer, "%.13G", dom.umax());
    g += buffer;
    g += pd;
    pad(g);
    sprintf(buffer, "%.13G", dom.vmin());
    g += buffer;
    g += pd;
    sprintf(buffer, "%.13G", dom.vmax());
    g += buffer;
    g += rd;
    pad(g);


    curr_dirent.line_count = (int)g.size()/64 - curr_dirent.param_data_start + 1;
    dirent.push_back(curr_dirent);
    Pcurr+=2;
}


//-----------------------------------------------------------------------------
void IGESconverter::writeIGEScurve(SplineCurve* curve, int colour, string& g,
                                   vector<IGESdirentry>& dirent, int& Pcurr,
				   int dependency)
//-----------------------------------------------------------------------------
{
    char pd = header_.pardel;
    char rd = header_.recdel;

    int dim = curve->dimension();
    int planar = (dim == 2) ? 1 : 0;

    IGESdirentry curr_dirent;
    curr_dirent.entity_type_number = 126;
    curr_dirent.structure = 0;
    curr_dirent.line_font_pattern = 1;
    curr_dirent.level = 0;
    curr_dirent.view = 0;
    curr_dirent.trans_matrix = 0;
    curr_dirent.label_display = 0;
    ostringstream os;
    os << "000" << dependency << (planar ? "0501" : "0001");
    curr_dirent.status = os.str();
    curr_dirent.line_weight = 0;
    curr_dirent.color = colour;
    curr_dirent.form = 0;
    curr_dirent.entity_label = "   nurbs";
    curr_dirent.entity_number = 0;

    int sz = (int)g.length();
    int ln = sz/64;
    DEBUG_ERROR_IF (ln*64 != sz, "String length not a multiple of 64.");

    curr_dirent.param_data_start = ln+1;


    g += "126";
    g += pd;

    char buffer[100]; // Enough digits for most numbers, I hope.

    sprintf(buffer, "%i", curve->numCoefs() - 1);
    g += buffer;
    g += pd;

    sprintf(buffer, "%i", curve->order() - 1);
    g += buffer;
    g += pd;

    sprintf(buffer, "%i", planar);
    g += buffer;
    g += pd;


//  	g += (curve->closedness() == GoOPEN) ? "0" : "1";
    g += "0";
    g += pd;

//  	if (curve->splinetype() == GoRATIONAL_BSPLINE 
//  	    || curve->splinetype() == GoRATIONAL_BEZIER)
//  	    g += "0";
//  	else
//  	    g += "1";
    int polynomial = (curve->rational()) ? 0 : 1;
    sprintf(buffer, "%i", polynomial);
    g += buffer;
    g += pd;
	
    // Periodicity
    g += "0";
    g += pd;

	// Pad with spaces until a 64-char line is filled
    pad(g);

	// Knots
    std::vector<double>::const_iterator it;
    for (it=curve->basis().begin();
         it!=curve->basis().end(); ++it) {
      sprintf(buffer, "%.15G", *it);
      if ((strlen(buffer)+g.length()+1)/64 > g.length()/64) {
		// Pad it till EOL
        pad(g);
      }
      g += buffer;
      g += pd;
    }       

	// Pad with spaces until a 64-char line is filled
    pad(g);

    int n = curve->numCoefs();
    int j;
    if (curve->rational())
    {
      std::vector<double>::const_iterator co = curve->rcoefs_begin();
      for (j=0; j<n; ++j) {
        sprintf(buffer, "%.13G",co[j*(dim+1)+dim]);
        if ((strlen(buffer)+g.length()+1)/64 > g.length()/64) {
              // Pad it till EOL
          pad(g);
        }
        g += buffer;
        g += pd;
      }
      pad(g);

      for (j=0; j<n; ++j) {
        sprintf(buffer, "%.13G",co[j*(dim+1)]/co[j*(dim+1)+dim]);
        g += buffer;
        g += pd;
        sprintf(buffer, "%.13G",co[j*(dim+1) + 1]/co[j*(dim+1)+dim]);
        g += buffer;
        g += pd;
        if (planar)
          sprintf(buffer, "0.0");
        else
          sprintf(buffer, "%.13G",co[j*(dim+1) + 2]/co[j*(dim+1)+dim]);
        g += buffer;
        g += pd;
        pad(g);
      }
    }
    else
    {
      for (j=0; j<n; ++j) {
        sprintf(buffer, "1.0");
        if ((strlen(buffer)+g.length()+1)/64 > g.length()/64) {
		// Pad it till EOL
          pad(g);
        }
        g += buffer;
        g += pd;
      }
      pad(g);

      std::vector<double>::const_iterator co = curve->coefs_begin();
      for (j=0; j<n; ++j) {
        sprintf(buffer, "%.13G",co[j*dim]);
        g += buffer;
        g += pd;
        sprintf(buffer, "%.13G",co[j*dim + 1]);
        g += buffer;
        g += pd;
        if (planar)
          sprintf(buffer, "0.0");
        else
          sprintf(buffer, "%.13G",co[j*dim + 2]);
        g += buffer;
        g += pd;
        pad(g);
      }
    }

    double tmin = curve->startparam();
    double tmax = curve->endparam();
	// Write parameter range
    sprintf(buffer, "%.13G", tmin);
    g += buffer;
    g += pd;
    sprintf(buffer, "%.13G", tmax);
    g += buffer;

    if (planar)
    {
      g += pd;
      g += "0.0";
      g += pd;
      g += "0.0";
      g += pd;
      g += "1.0";
    }
    g += rd;
    pad(g);


    curr_dirent.line_count = (int)g.size()/64 - curr_dirent.param_data_start + 1;
    dirent.push_back(curr_dirent);
    Pcurr+=2;
}


//-----------------------------------------------------------------------------
void IGESconverter::
writeIGESboundary(shared_ptr<CurveLoop>  loop,
                  string& g, vector<IGESdirentry>& dirent,
                  int Psurf, int& Pcurr,
		  int dependency)
//-----------------------------------------------------------------------------
{
    char pd = header_.pardel;
    char rd = header_.recdel;

    IGESdirentry curr_dirent;
    curr_dirent.entity_type_number = 141;
    curr_dirent.structure = 0;
    curr_dirent.line_font_pattern = 1;
    curr_dirent.level = 0;
    curr_dirent.view = 0;
    curr_dirent.trans_matrix = 0;
    curr_dirent.label_display = 0;
    ostringstream os;
    os << "000" << dependency << "0001";
    curr_dirent.status = os.str();
    curr_dirent.line_weight = 0;
    curr_dirent.color = 0;
    curr_dirent.form = 0;
    curr_dirent.entity_label = "boundary";
    curr_dirent.entity_number = 0;

        // Write curves
    int nmb_space = loop->size();
    vector<std::pair<int,int> > Ppointer;
    int p1, p2;
    SplineCurve *space, *param;
    shared_ptr<SplineCurve> tmp_space;
    CurveOnSurface *curr_crv;
    int col = 0;
    for (int i=0; i<nmb_space; i++)
    {
      curr_crv = dynamic_cast<CurveOnSurface*>((*loop)[i].get());
      space = dynamic_cast<SplineCurve*>(curr_crv->spaceCurve().get());
      param = dynamic_cast<SplineCurve*>(curr_crv->parameterCurve().get());
      //      ASSERT((curr_crv != NULL) && (space != NULL) && (param != NULL));
      ASSERT((curr_crv != NULL) && ((space != NULL) || (param != NULL)));
      if (curr_crv->spaceCurve().get() && !space)
	{
	  tmp_space = 
	    shared_ptr<SplineCurve>(curr_crv->spaceCurve()->geometryCurve());
	  space = tmp_space.get();
	}
//        if (!space)
// 	 {
// 	   space = dynamic_cast<SplineCurve*>(curr_crv->createSpaceCurve());
// 	   free_space = true;
// 	 }
       if (space) {
	  writeIGEScurve(space, col, g, dirent, Pcurr, 1); // Physically dependent
	  p1 = Pcurr;
      } else {
	  p1 = 0;
      }

      if (param)
      {
        writeIGEScurve(param, col, g, dirent, Pcurr, 1);
        p2 = Pcurr;
      }
      else
        p2 = 0;

      Ppointer.push_back(std::make_pair(p1,p2));
       }
    
    int sz = (int)g.length();
    int ln = sz/64;
    DEBUG_ERROR_IF (ln*64 != sz, "String length not a multiple of 64.");

    curr_dirent.param_data_start = ln+1;


    g += "141";
    g += pd;

    char buffer[100]; // Enough digits for most numbers, I hope.

    curr_crv = dynamic_cast<CurveOnSurface*>((*loop)[0].get());
    int parexist = (curr_crv->parameterCurve().get() != 0);
    sprintf(buffer, "%i", parexist);
    g += buffer;
    g += pd;


    int pref = (curr_crv->parPref()) ? 2 : 1;
    sprintf(buffer, "%i", pref);
    g += buffer;
    g += pd;

        // Write surface pointer
    sprintf(buffer, "%i", Psurf);
    g += buffer;
    g += pd;

    sprintf(buffer, "%i", nmb_space);
    g += buffer;

    for (int i=0; i<nmb_space; i++)
    {
      g += pd;
      
          // Fetch the P-number of the space curve
      sprintf(buffer, "%i", Ppointer[i].first);
      if ((strlen(buffer)+g.length()+2)/64 > g.length()/64) 
		// Pad it till EOL
	pad(g);
      g += buffer;
      g += pd;

          // Not reversed
      g += "1";
      g += pd;

      if (Ppointer[i].second > 0)
      {
          // One parameter curve
        g += "1";
        g += pd;
        sprintf(buffer, "%i", Ppointer[i].second);
	if ((strlen(buffer)+g.length()+2)/64 > g.length()/64) 
		// Pad it till EOL
	  pad(g);
        g += buffer;
      }
      else
      {
            // No parameter curve
        g += "0";
      }
    }
    
    g += rd;
    pad(g);
    
    curr_dirent.line_count = (int)g.size()/64 - curr_dirent.param_data_start + 1;
    dirent.push_back(curr_dirent);
    Pcurr+=2;
}

//-----------------------------------------------------------------------------
void IGESconverter::
writeIGESboundedSurf(BoundedSurface *surf, int colour, string& g, 
                     vector<IGESdirentry>& dirent, int& Pcurr)
//-----------------------------------------------------------------------------
{
    char pd = header_.pardel;
    char rd = header_.recdel;
    
    //    std::cout << surf << " underlying " << surf->underlyingSurface().get() << std::endl;

    IGESdirentry curr_dirent;
    curr_dirent.entity_type_number = 143;
    curr_dirent.structure = 0;
    curr_dirent.line_font_pattern = 1;
    curr_dirent.level = 0;
    curr_dirent.view = 0;
    curr_dirent.trans_matrix = 0;
    curr_dirent.label_display = 0;
    curr_dirent.status = "00000001";
    curr_dirent.line_weight = 0;
    curr_dirent.color = colour;
    curr_dirent.form = 0;
    curr_dirent.entity_label = " bounded";
    curr_dirent.entity_number = 0;

        // First write the surface
    int Psurf;
    shared_ptr<ParamSurface> srf = surf->underlyingSurface();
    int col = 0;
    //ASSERT(srf->instanceType() == Class_SplineSurface);
    shared_ptr<SplineSurface> spline_sf = 
      dynamic_pointer_cast<SplineSurface, ParamSurface>(srf);
    if (!spline_sf.get())
      {
	ElementarySurface *elem_sf = dynamic_cast<ElementarySurface*>(srf.get());
	if (elem_sf)
	  spline_sf = shared_ptr<SplineSurface>(elem_sf->geometrySurface());
      }
    ASSERT(spline_sf.get() != 0);

    writeIGESsurface(spline_sf.get(), col,
                     g, dirent, Pcurr, 1); // Physically dependent
    Psurf = Pcurr;
    
        // Write all loops
    int i;
    int nmbloop = surf->numberOfLoops();
    vector<int> Ppointer;
    for (i=0; i<nmbloop; i++)
    {
      shared_ptr<CurveLoop>  curr_loop = surf->loop(i);
      writeIGESboundary(curr_loop, g, dirent, Psurf, Pcurr, 1); // Physically dependent
      Ppointer.push_back(Pcurr);
    }
    
    int sz = (int)g.length();
    int ln = sz/64;
    DEBUG_ERROR_IF (ln*64 != sz, "String length not a multiple of 64.");

    curr_dirent.param_data_start = ln+1;


    g += "143";
    g += pd;

    char buffer[100]; // Enough digits for most numbers, I hope.

    shared_ptr<ParamCurve> cv = (*(surf->loop(0).get()))[0];
    int parexist =
	(dynamic_pointer_cast<CurveOnSurface, ParamCurve>(cv).get() != 0);

    sprintf(buffer, "%i", parexist);
    g += buffer;
    g += pd;

        // Pointer to surface
    sprintf(buffer, "%i", Psurf);
    g += buffer;
    g += pd;

        // Number of loops
    sprintf(buffer, "%i", nmbloop);
    g += buffer;

        // Pointer to all loops
    for (i=0; i<nmbloop; i++)
    {
      g += pd;
      
      sprintf(buffer, "%i", Ppointer[i]);
      if ((strlen(buffer)+g.length()+2)/64 > g.length()/64) 
		// Pad it till EOL
	pad(g);
      g += buffer;
    }
    
    g += rd;
    pad(g);
    
    curr_dirent.line_count = (int)g.size()/64 - curr_dirent.param_data_start + 1;
    dirent.push_back(curr_dirent);
    Pcurr+=2;
}


//-----------------------------------------------------------------------------
IGESdirentry IGESconverter::readIGESdirentry(const char* start)
//-----------------------------------------------------------------------------
{
    char buffer[9];
    IGESdirentry ent;

    strncpy(buffer, start + 0, 9);
    ent.entity_type_number = atoi(buffer);

    strncpy(buffer, start + 8, 9);
    ent.param_data_start = atoi(buffer);

    strncpy(buffer, start + 16, 9);
    ent.structure = atoi(buffer);

    strncpy(buffer, start + 24, 9);
    ent.line_font_pattern = atoi(buffer);

    strncpy(buffer, start + 32, 9);
    ent.level = atoi(buffer);

    strncpy(buffer, start + 40, 9);
    ent.view = atoi(buffer);

    strncpy(buffer, start + 48, 9);
    ent.trans_matrix = atoi(buffer);

    strncpy(buffer, start + 56, 9);
    buffer[8] = '\0'; // As the next element is not a space (which
		      // atoi expects).
    ent.label_display = atoi(buffer);

    strncpy(buffer, start + 64, 9);
    buffer[8] = '\0';
    ent.status = buffer;

    strncpy(buffer, start + 80, 9);
    ent.line_weight = atof(buffer);  // @ DANGER, possible D-notation trouble

    strncpy(buffer, start + 88, 9);
    ent.color = atoi(buffer);

    strncpy(buffer, start + 96, 9);
    ent.line_count = atoi(buffer);

    strncpy(buffer, start + 104, 9);
    ent.form = atoi(buffer);

    strncpy(buffer, start + 128, 9);
    buffer[8] = '\0';
    ent.entity_label = buffer;

    strncpy(buffer, start + 136, 9);
    ent.entity_number = atoi(buffer);

    return ent;
}


//-----------------------------------------------------------------------------
int IGESconverter::whichPLine(ccp whereami)
//-----------------------------------------------------------------------------
{
    int offset = (int)(whereami - start_of_P_section_)/64;
    return offset + 1;
}


//-----------------------------------------------------------------------------
void IGESconverter::skipDelimiter(ccp& whereami, char pd)
//-----------------------------------------------------------------------------
{
    //  cout << "*whereami = " << *whereami << endl
    // << "pd        = " << pd << endl;

    // Skip whitespace
    while (isspace(*whereami)) ++whereami;
    if ((*whereami) == pd)
	++whereami;
    else {
	//char buf[101];
	//cout << endl << pd << endl;
	//cout << "-------------------------------------" << endl;
	//strncpy(buf,whereami,100);
	//cout << buf << endl;
	int line = whichPLine(whereami);
	MESSAGE ("Parsing anomaly, could not locate separator (" << pd << ") on p-line " << line);
    }
}


//-----------------------------------------------------------------------------
bool IGESconverter::checkDelimiter(ccp& whereami, char wanted_delimiter,
				   char alternative_delimiter)
//-----------------------------------------------------------------------------
{
    //  cout << "*whereami = " << *whereami << endl
    // << "pd        = " << pd << endl;

    // Skip whitespace
    while (isspace(*whereami)) ++whereami;
    if ((*whereami) == wanted_delimiter)
      {
	++whereami;
	return true;
      }
    else if ((*whereami) == alternative_delimiter)
      {
	++whereami;
	return false;
      }
    else 
      {
	//char buf[101];
	//cout << endl << pd << endl;
	//cout << "-------------------------------------" << endl;
	//strncpy(buf,whereami,100);
	//cout << buf << endl;
	int line = whichPLine(whereami);
	MESSAGE ("Parsing anomaly, could not locate separator (" << wanted_delimiter 
		 << " or " << alternative_delimiter<< ") on p-line " << line);
	return false;
      }
}


//-----------------------------------------------------------------------------
void IGESconverter::skipOptionalTrailingArguments(ccp& whereami,
						  char pd, char rd,
						  int max_to_skip)
//-----------------------------------------------------------------------------
{
    for (int j = 0; j < max_to_skip; ++j) {
	bool ismore = checkDelimiter(whereami, pd, rd);
	if (ismore) {
	    // Read past the extra pointers in the file
	    int npointer = readIGESint(whereami, pd, rd);
	    int pdummy;
	    for (int i = 0; i < npointer; i++) {
		skipDelimiter(whereami, pd);
		pdummy = readIGESint(whereami, pd, rd);
	    }
	} else {
	    break;
	}
    }
}


//-----------------------------------------------------------------------------
string IGESconverter::readIGESstring(ccp& start, char delim, char delim2)
//-----------------------------------------------------------------------------
{
    while (isspace(*start)) ++start;
    char numbuf[8];
    int numdig = 0;
    if (start[0] == delim || start[0] == delim2) {
	return string();
    }
    for(int i=0; i<7; i++) {
	if (start[i] == 'H') {
	    numdig = i;
	    break;
	} else {
	    numbuf[i] = start[i];
	}
    }
    numbuf[numdig] = 0; // Terminate numbuf
    //        cout << numdig << " " << numbuf << endl;
    int numchars = atoi(numbuf);
    if (numchars < 1) {
	MESSAGE("Less than one character in string: " << numchars);
	return string();
    }
    start += numdig + 1 + numchars;
    //    cout << "[[[    " << numchars << "    ]]]" << endl;
    return string(start-numchars, numchars);
}

//-----------------------------------------------------------------------------
string IGESconverter::writeIGESstring(const string& instring)
//-----------------------------------------------------------------------------
{
    char digits[10];
    sprintf(digits, "%i", (int)instring.length());
    string out = "";
    if (atoi(digits) == 0)
      return out;
    out += digits;
    out += 'H';
    out += instring;
    return out;
}


//-----------------------------------------------------------------------------
double IGESconverter::readIGESdouble(ccp& start, char pd, char rd)
//-----------------------------------------------------------------------------
{
    while (isspace(*start)) ++start;
    char numbuf[32];
    int numdig = 0;
    for(int i=0; i<31; i++) {
	if (start[i] == pd || start[i] == rd) {
	    numdig = i;
	    break;
	} else {
	    numbuf[i] = start[i];
	    if (numbuf[i] == 'D') numbuf[i] = 'E'; // FP notation...
	}
    }
    start += numdig; // Next value.
    // We must also disregard any trailing white spaces in numbuf.
    int nmb_trailing_spaces = 0;
    while (isspace(start[-1-nmb_trailing_spaces]))
	++nmb_trailing_spaces;
    numbuf[numdig-nmb_trailing_spaces] = 0; // Terminate numbuf

    std::stringstream ss ( numbuf );
    double val = -1.0;
    ss>> val;

    return val;
}





//-----------------------------------------------------------------------------
string IGESconverter::writeIGESdouble(double d)
//-----------------------------------------------------------------------------
{
    char number[30];
    sprintf(number, "%.15G", d);
    return string(number);
}





//-----------------------------------------------------------------------------
int IGESconverter::readIGESint(ccp& start, char pd, char rd)
//-----------------------------------------------------------------------------
{
    while (isspace(*start)) ++start;
    char numbuf[32];
    int numdig = 0;
    for(int i=0; i<31; i++) {
	if (start[i] == pd || start[i] == rd) {
	    numdig = i;
	    break;
	} else {
	    numbuf[i] = start[i];
	}
    }
    numbuf[numdig] = 0; // Terminate numbuf
    start += numdig;
    return atoi(numbuf);
}

//-----------------------------------------------------------------------------
string IGESconverter::writeIGESint(int i)
//-----------------------------------------------------------------------------
{
    char number[30];
    sprintf(number, "%i", i);
    return string(number);
}

//-----------------------------------------------------------------------------
vector<double> IGESconverter::transformFromIGEScolour(int i)
//-----------------------------------------------------------------------------
{
    vector<double> colour(3, 0.0);
    if (i == 0) // No colour
	colour.clear();
    else if (i == 1) // Black
	; // Black, we do nothing.
    else if (i == 2) // Red
	colour[0] = 100.0;
    else if (i == 3) // Green
	colour[1] = 100.0;
    else if (i == 4) // Blue
	colour[2] = 100.0;
    else if (i == 5) { // Yellow
	colour[0] = 100.0;
	colour[1] = 100.0;
    } else if (i == 6) { // Magenta
	colour[0] = 100.0;
	colour[2] = 100.0;
    } else if (i == 7) { // Cyan
	colour[1] = 100.0;
	colour[2] = 100.0;
    } else if (i == 8) { // White
	colour[0] = 100.0;
	colour[1] = 100.0;
	colour[2] = 100.0;
    }
    return colour;
}

//-----------------------------------------------------------------------------
vector<vector<double> > IGESconverter::uniqueColours(vector<vector<double> >& colours)
//-----------------------------------------------------------------------------
{
    vector<vector<double> > unique_colours; // We're only interested in the first 3 vals (skipping alpha).
    for (size_t ki = 0; ki < colours.size(); ++ki) {
	if (colours[ki].size() == 0)
	    continue;
	size_t kj;
	for (kj = 0; kj < unique_colours.size(); ++kj) {
	    vector<double> set_int;
	    set_intersection(colours[ki].begin(), colours[ki].begin() + 3,
			     unique_colours[kj].begin(), unique_colours[kj].begin() + 3,
			     std::back_inserter(set_int));
	    if (set_int.size() == 3)
		break;
	}
	if (kj == unique_colours.size()) { // We couldn't find colour triplet.
	    vector<double> new_colour(colours[ki].begin(), colours[ki].begin() + 3);
	    unique_colours.push_back(new_colour);
	}
    }

    return unique_colours;
}

//-----------------------------------------------------------------------------
shared_ptr<SplineSurface> IGESconverter::readIGESsurface(const char* start,
							   int num_lines)
//-----------------------------------------------------------------------------
{
    char pd = header_.pardel;
    char rd = header_.recdel;

    // Verify that start points to '128,.....' warn if not
    int type = readIGESint(start, pd, rd);
    DEBUG_ERROR_IF(type!=128, "Entity is not of type 128.");
    skipDelimiter(start, pd);

    // Read parameters
    int n1 = readIGESint(start, pd, rd) + 1;
    skipDelimiter(start, pd);
    int n2 = readIGESint(start, pd, rd) + 1;
    skipDelimiter(start, pd);
    int k1 = readIGESint(start, pd, rd) + 1;
    skipDelimiter(start, pd);
    int k2 = readIGESint(start, pd, rd) + 1;
    skipDelimiter(start, pd);
    int clo1, clo2;
    clo1 = readIGESint(start, pd, rd);
    skipDelimiter(start, pd);
    clo2 = readIGESint(start, pd, rd);
    skipDelimiter(start, pd);
    int polynomial = readIGESint(start, pd, rd);
    int rational = (polynomial == 0);
    skipDelimiter(start, pd);
    int per1, per2;
    per1 = readIGESint(start, pd, rd);
    skipDelimiter(start, pd);
    per2 = readIGESint(start, pd, rd);
    skipDelimiter(start, pd);

    int N = n1+k1;
    int i=0, j=0;
    std::vector<double> knot1(N);
    //cout << "\nknot1: \n";
    for (i=0; i<N; ++i) {
	knot1[i] = readIGESdouble(start, pd, rd);
	skipDelimiter(start, pd);
	//cout << knot1[i] << endl;
    }

    N = n2+k2;
    //cout << "\nknot2: \n";
    std::vector<double> knot2(N);
    for (i=0; i<N; ++i) {
	knot2[i] = readIGESdouble(start, pd, rd);
	skipDelimiter(start, pd);
	//cout << knot2[i] << endl;
    }

    bool all_weights_are_one = true;
    N = n1*n2;
    std::vector<double> weights(N);
    for (i=0; i<N; ++i) {
	weights[i] = readIGESdouble(start, pd, rd);
	if (weights[i] != 1.0) all_weights_are_one = false;
	skipDelimiter(start, pd);
    }

    N = n1*n2*3;
    std::vector<double> coefs(N);
    for (i=0; i<N; ++i) {
	coefs[i] = readIGESdouble(start, pd, rd);
	skipDelimiter(start, pd);
    }

    double u1, u2, v1, v2;
    u1 = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    u2 = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    v1 = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    v2 = readIGESdouble(start, pd, rd);

    skipOptionalTrailingArguments(start, pd, rd);
  

//      GoClosedness c1 = GoOPEN;
//      if (clo1==1) c1 = GoCLOSED;
//      if (per1==1 && k1<=5) c1 = static_cast<GoClosedness>(k1-2);

//      GoClosedness c2 = GoOPEN;
//      if (clo2==1) c2 = GoCLOSED;
//      if (per2==1 && k2<=5) c2 = static_cast<GoClosedness>(k2-2);

//     shared_ptr<SplineSurface> srf;
    SplineSurface *sfptr;
    if (rational == 1 && !all_weights_are_one)
    {
      N = n1*n2;
      vector<double> coefs2(4*N);
      for (i=0; i<N; ++i) {
        for (j=0; j<3; j++)
          coefs2[4*i+j] = coefs[3*i+j]*weights[i];
        coefs2[4*i+3] = weights[i];
      }
      sfptr = new SplineSurface( n1, n2, k1, k2, knot1.begin(),
                                   knot2.begin(),
                                   coefs2.begin(), 3, true);
    }
    else
      sfptr = new SplineSurface( n1, n2, k1, k2, knot1.begin(),
                                   knot2.begin(),
                                   coefs.begin(), 3);
    shared_ptr<SplineSurface> srf(sfptr);

    return srf;
}


//-----------------------------------------------------------------------------
shared_ptr<BoundedCurve> IGESconverter::readIGESline(const char* start,
						   int num_lines,
						   int direntry_index)
//-----------------------------------------------------------------------------
{
    char pd = header_.pardel;
    char rd = header_.recdel;

    // Verify that start points to '110,.....' warn if not
    int type = readIGESint(start, pd, rd);
    DEBUG_ERROR_IF(type!=110, "Entity is not of type 110.");
    skipDelimiter(start, pd);

    Array<double, 3> pnt1;
    Array<double, 3> pnt2;
    // Read parameters
    pnt1[0] = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    pnt1[1] = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    pnt1[2] = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    pnt2[0] = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    pnt2[1] = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    pnt2[2] = readIGESdouble(start, pd, rd);

    skipOptionalTrailingArguments(start, pd, rd);

    // Extract the coordinate system
    int csentry = direntries_[direntry_index].trans_matrix;
    CoordinateSystem<3> cs; 
    if (csentry != 0) { // If value of directory entry is 0, we should use identity.
	map< int, CoordinateSystem<3> >::iterator it
	    = coordsystems_.find(csentry);
	if (it == coordsystems_.end()) {
	    MESSAGE("Could not find the referred coordinate system ("
		    << csentry << ") in the file. Using identity.");
	} else {
	    cs = it->second;
	}
    }

    pnt1 = cs*pnt1;
    pnt2 = cs*pnt2;
    Point p1(pnt1.begin(), pnt1.end());
    Point p2(pnt2.begin(), pnt2.end());
    Point dir = p2 - p1;

//     shared_ptr<SplineCurve> crv(new SplineCurve(p1, 0.0, p2, 1.0));
    shared_ptr<Line> crv(new Line(p1, dir));
    crv->setParameterInterval(0.0, 1.0);
    plane_normal_.push_back(Point());

    shared_ptr<BoundedCurve> bd_cv(new BoundedCurve(crv, p1, p2));

    return bd_cv;
}


//-----------------------------------------------------------------------------
shared_ptr<PointCloud3D> IGESconverter::readIGESpointCloud(const char* start,
							   int num_lines)
//-----------------------------------------------------------------------------
{
    char pd = header_.pardel;
    char rd = header_.recdel;

    // Verify that start points to '116,.....' warn if not
    int type = readIGESint(start, pd, rd);
    DEBUG_ERROR_IF(type!=116, "Entity is not of type 116.");
    skipDelimiter(start, pd);

    // Read parameters
    Point pt(3);

    pt[0] = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    pt[1] = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    pt[2] = readIGESdouble(start, pd, rd);

    skipOptionalTrailingArguments(start, pd, rd);

    shared_ptr<PointCloud3D> pt_cl(new PointCloud3D(pt.begin(), 1));
    return pt_cl;
}


//-----------------------------------------------------------------------------
shared_ptr<PointCloud3D> IGESconverter::readIGESdirection(const char* start,
							  int num_lines)
//-----------------------------------------------------------------------------
{
    char pd = header_.pardel;
    char rd = header_.recdel;

    // Verify that start points to '123,.....' warn if not
    int type = readIGESint(start, pd, rd);
    DEBUG_ERROR_IF(type!=123, "Entity is not of type 123.");
    skipDelimiter(start, pd);

    // Read parameters
    Point pt(3);

    pt[0] = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    pt[1] = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    pt[2] = readIGESdouble(start, pd, rd);

    skipOptionalTrailingArguments(start, pd, rd);

    shared_ptr<PointCloud3D> pt_cl(new PointCloud3D(pt.begin(), 1));
    return pt_cl;
}


//-----------------------------------------------------------------------------
shared_ptr<SplineCurve> IGESconverter::readIGEScurve(const char* start,
						     int num_lines,
						     int direntry_index)
//-----------------------------------------------------------------------------
{
    char pd = header_.pardel;
    char rd = header_.recdel;

    // Verify that start points to '126,.....' warn if not
    int type = readIGESint(start, pd, rd);
    DEBUG_ERROR_IF(type!=126, "Entity is not of type 126.");
    skipDelimiter(start, pd);

    // Read parameters
    int n1 = readIGESint(start, pd, rd) + 1;
    skipDelimiter(start, pd);
    int k1 = readIGESint(start, pd, rd) + 1;
    skipDelimiter(start, pd);
    int planar = readIGESint(start, pd, rd);
    skipDelimiter(start, pd);
    int clo1;
    clo1 = readIGESint(start, pd, rd);
    skipDelimiter(start, pd);
    int polynomial = readIGESint(start, pd, rd);
    int rational = (polynomial == 0);
    skipDelimiter(start, pd);
    int per1;
    per1 = readIGESint(start, pd, rd);
    skipDelimiter(start, pd);

    int N = n1+k1;
    int i=0, j=0;
    std::vector<double> knot1(N);
    //cout << "\nknot1: \n";
    for (i=0; i<N; ++i) {
	knot1[i] = readIGESdouble(start, pd, rd);
	skipDelimiter(start, pd);
	//cout << knot1[i] << endl;
    }


    bool all_weights_are_one = true;
    N = n1;
    std::vector<double> weights(N);
    for (i=0; i<N; ++i) {
	weights[i] = readIGESdouble(start, pd, rd);
	if (weights[i] != 1.0) all_weights_are_one = false;
	skipDelimiter(start, pd);
    }

    N = n1*3;
    std::vector<double> coefs(N);
    for (i=0; i<N; ++i) {
	coefs[i] = readIGESdouble(start, pd, rd);
	skipDelimiter(start, pd);
    }

    double u1 = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    double u2 = readIGESdouble(start, pd, rd);

    double norm[3];
    int dim = 3;
    double *coefptr = &coefs[0];
    if (planar)
    {
      for (i=0; i<3; ++i) {
	skipDelimiter(start, pd);
        norm[i] = readIGESdouble(start, pd, rd);
      }
      //      skipDelimiter(start, rd);
	
      plane_normal_.push_back(Point(norm[0],norm[1],norm[2]));
    }
    else
      plane_normal_.push_back(Point());

    skipOptionalTrailingArguments(start, pd, rd);

    // Extract the coordinate system
    int csentry = direntries_[direntry_index].trans_matrix;
    CoordinateSystem<3> cs; 
    if (csentry != 0) { // If value of directory entry is 0, we should
			// use identity.
	map< int, CoordinateSystem<3> >::iterator it
	    = coordsystems_.find(csentry);
	if (it == coordsystems_.end()) {
	    MESSAGE("Could not find the referred coordinate system ("
		    << csentry << ") in the file. Using identity.");
	} else {
	    cs = it->second;
	    MESSAGE("Transformation matrix for spline curve object "
		    "is missing!");
	}
    }

//     // Need to skip ut to 5 arguments, as we can have
//     // both a normal (an error, but common) and extra
//     // property pointers.
//     //    skipOptionalTrailingArguments(start, pd, rd, 5);
//     // @afr: Instead, we just ignore everything until rd (usually semicolon).
//     while ((*start) != rd) {
// 	++start;
//     }
//     ++start;

    SplineCurve *cvptr;
    if (rational == 1 && !all_weights_are_one)
    {
      int dim1 = dim+1;
      vector<double> coefs2(dim1*n1);
      for (i=0; i<n1; ++i) {
        for (j=0; j<dim; j++)
          coefs2[dim1*i+j] = coefptr[dim*i+j]*weights[i];
        coefs2[dim1*i+3] = weights[i];
      }
      cvptr = new SplineCurve(n1, k1, knot1.begin(), coefs2.begin(),
                              dim, true);
    }
    else
      cvptr = new SplineCurve(n1, k1, knot1.begin(), coefptr, dim);
    shared_ptr<SplineCurve> crv = (shared_ptr<SplineCurve>)(cvptr);

    double fuzzypar = 1e-10;
    if ((fabs(u1-crv->startparam()) > fuzzypar) ||
	(fabs(crv->endparam() - u2) > fuzzypar))
    {
#ifndef NDEBUG
        MESSAGE("Extracting subcurve for entity 126!");
	if (u1 < crv->startparam() || u2 > crv->endparam())
	  MESSAGE("Parameter value(s) outside domain, moving them inside!");
#endif
// #if 0
//         MESSAGE("Should extract subcurve! Skipping for now ...");
// #endif
	if (u1 < crv->startparam())
	  u1 = crv->startparam();
	if (u2 > crv->endparam())
	  u2 = crv->endparam();
	shared_ptr<SplineCurve> sub_crv(crv->subCurve(u1, u2));
	crv = sub_crv;
    }

    return crv;
}


//-----------------------------------------------------------------------------
shared_ptr<Plane> IGESconverter::readIGESplane(const char* start,
					       int num_lines, int form)
//-----------------------------------------------------------------------------
{
    char pd = header_.pardel;
    char rd = header_.recdel;

    // Verify that start points to '108,.....' warn if not
    int type = readIGESint(start, pd, rd);
    DEBUG_ERROR_IF(type!=108, "Entity is not of type 108.");
    skipDelimiter(start, pd);

    // Read the data
    // Description of plane: ta*x + tb*y + tz+c = d
    double ta = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    double tb = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    double tc = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    double td = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);

    // Pointer to trimming loop. Not treated at this stage
    ASSERT(form == 0); // Unbounded plane.
    int crv_id;
    crv_id = readIGESint(start, pd, rd); 
    skipDelimiter(start, pd);
   
    // Location point and size for visualization
    double p1 = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    double p2 = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    double p3 = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    double tsize;
    tsize  = readIGESdouble(start, pd, rd);

    skipOptionalTrailingArguments(start, pd, rd);

    shared_ptr<Plane> planarsf = shared_ptr<Plane>(new Plane(ta, tb, tc, td));

    double dist;
    dist = planarsf->distance(Point(p1,p2,p3));
    return planarsf;
}


//-----------------------------------------------------------------------------
shared_ptr<ElementaryCurve>
IGESconverter::readIGESconicArc(const char* start,
				int num_lines, int form)
//-----------------------------------------------------------------------------
{
    shared_ptr<ElementaryCurve> conic_arc;

    DEBUG_ERROR_IF((form != 1) && (form != 2) && (form != 3),
		   "Form should be 1, 2 or 3, not " << form << "!");

    char pd = header_.pardel;
    char rd = header_.recdel;

    // Verify that start points to '104,.....' warn if not
    int type = readIGESint(start, pd, rd);
    DEBUG_ERROR_IF(type!=104, "Entity is not of type 104.");
    skipDelimiter(start, pd);

    // Read the data
    double ta = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    double tb = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    double tc = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    double td = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    double te = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    double tf = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    double zt = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    double x1 = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    double y1 = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    double x2 = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    double y2 = readIGESdouble(start, pd, rd);
    skipOptionalTrailingArguments(start, pd, rd);

    // We then decide what type of conic arc (ellipse, hyperbola or
    // parabola) we are dealing with, by computing some determinants
    // and such.
    double q1 = ta*(tc*tf - 0.25*te*te) -
	0.5*tb*(0.5*tb*tf - 0.25*td*te) + 0.5*td*(0.25*tb*te - 0.5*tc*td);
    double q2 = ta*tc - 0.25*tb*tb;
    double q3 = ta + tc;

    // The expressions in the IGES manual does not seem to be correct,
    // using mathworld.wolfram.com.
    if ((q1 != DZERO) && (q2 > 0.0) && (q1/q3 < 0.0)) {
	;//MESSAGE("Just found an ellipse!");
    } else if ((q2 < 0.0) && (q2 != DZERO)) {
	;//MESSAGE("Just found a hyperbola!");
    } else if ((q2 == DZERO) && (q1 != DZERO)) {
	;//MESSAGE("Just found a parabola!");
    } else {
	MESSAGE("Input values for conic arc do not match an "
		"ellipse, hyperbola or parabola! q1 = " <<
		q1 << ", q2 = " << q2 << ", q3 = " << q3);
    }

    if (form == 1) {
#ifdef SBR_DBG
	std::cout << "ta=" << ta << ", tb=" << tb << ", tc=" << tc <<
	    ", td=" << td << ", te=" << te << ", tf=" << tf << ", zt=" <<
	    zt << ", x1=" << x1 << ", y1=" << y1 << ", x2=" << x2 <<
	    ", y2=" << y2 << std::endl;
#endif

	double a = sqrt(-(tf/ta));
	double b = sqrt(-(tf/tc));
	double r1 = a/x2;
	double r2 = b/y2;
	Point pos(x1, y1, zt);
	Point ref_dir(x2, y2, 0.0);
	Point normal(0.0, 0.0, 1.0);
	conic_arc =
	    shared_ptr<Ellipse>(new Ellipse(pos, ref_dir, normal, r1, r2));
    } else if (form == 2) {
	MESSAGE("Form number for a hyperbola, not yet supported.");
    } else if (form == 3) {
	MESSAGE("Form number for a parabola, not yet supported.");
    }

    return conic_arc;
}


//-----------------------------------------------------------------------------
shared_ptr<Plane>
IGESconverter::readIGESplaneSurface(const char* start,
				    int num_lines, int form)
//-----------------------------------------------------------------------------
{
    shared_ptr<Plane> plane;

    char pd = header_.pardel;
    char rd = header_.recdel;

    DEBUG_ERROR_IF((form != 0) && (form != 1),
		   "Form should be 0 or 1, not " << form << "!");

    // Verify that start points to '190,.....' warn if not
    int type = readIGESint(start, pd, rd);
    DEBUG_ERROR_IF(type!=190, "Entity is not of type 190.");
    skipDelimiter(start, pd);
    int pos_id = readIGESint(start, pd, rd); // Pointer to location entity.
    skipDelimiter(start, pd);
    int normal_id = readIGESint(start, pd, rd); // Pointer to normal entity.
    // We do not have to worry about 
    int ref_dir_id = -1; // Pointer to reference dir entity.
    if (form == 1) {
	skipDelimiter(start, pd);
	ref_dir_id = readIGESint(start, pd, rd); 
    }
    // Not reading optional pointers (text etc).
    skipOptionalTrailingArguments(start, pd, rd);

    // We then fetch the location, normal and possibly reference dir.
    // Parametrization of sf:
    // x = <
    // s(u,v) = loc + u

    size_t ki;
    for (ki = 0; ki < local_geom_.size(); ++ki)
	if (geom_id_[ki] == pos_id)
	    break;
    DEBUG_ERROR_IF(ki==local_geom_.size(),
		   "Missing object. Could not find pos #" << pos_id);
    shared_ptr<PointCloud3D> pos_cl =
	dynamic_pointer_cast<PointCloud3D, GeomObject>(local_geom_[ki]);
    Point pos(pos_cl->rawData()[0], pos_cl->rawData()[1], pos_cl->rawData()[2]);
    geom_used_[ki] = 1;
    for (ki = 0; ki < local_geom_.size(); ++ki)
	if (geom_id_[ki] == normal_id)
	    break;
    DEBUG_ERROR_IF(ki==local_geom_.size(),
		   "Missing object. Could not find dir #" << normal_id);
    shared_ptr<PointCloud3D> normal_cl =
	dynamic_pointer_cast<PointCloud3D, GeomObject>(local_geom_[ki]);
    Point normal(normal_cl->rawData()[0],
		 normal_cl->rawData()[1], normal_cl->rawData()[2]);
    geom_used_[ki] = 1;

    if (form == 0) {
	// We do not have to worry about parametrizing the surface.	
	plane = shared_ptr<Plane>(new Plane(pos, normal));
    } else {
	for (ki = 0; ki < local_geom_.size(); ++ki)
	    if (geom_id_[ki] == ref_dir_id)
		break;
	DEBUG_ERROR_IF(ki==local_geom_.size(),
		       "Missing object. Could not find dir #" << ref_dir_id);
	shared_ptr<PointCloud3D> dir_cl =
	    dynamic_pointer_cast<PointCloud3D, GeomObject>(local_geom_[ki]);
	Point dir(dir_cl->rawData()[0],
		  dir_cl->rawData()[1], dir_cl->rawData()[2]);
	geom_used_[ki] = 1;
	plane = shared_ptr<Plane>(new Plane(pos, normal, dir));
    }

    return plane;
}


//-----------------------------------------------------------------------------
shared_ptr<Go::Cylinder>
IGESconverter::readIGESrightCircularCylindricalSurface(const char* start,
						       int num_lines, int form)
//-----------------------------------------------------------------------------
{
    shared_ptr<Cylinder> cylinder;

    DEBUG_ERROR_IF((form != 0) && (form != 1),
		   "Form should be 0 or 1, not " << form << ".");

    char pd = header_.pardel;
    char rd = header_.recdel;

    // Verify that start points to '192,.....' warn if not
    int type = readIGESint(start, pd, rd);
    DEBUG_ERROR_IF(type!=192, "Entity is not of type 192.");
    skipDelimiter(start, pd);

    int pos_id = readIGESint(start, pd, rd);
    skipDelimiter(start, pd);
    int axis_id = readIGESint(start, pd, rd);
    skipDelimiter(start, pd);
    double radius = readIGESdouble(start, pd, rd);
    int ref_dir_id = -1;
    if (form == 1) {
	skipDelimiter(start, pd);
	ref_dir_id = readIGESint(start, pd, rd);
    }
    // Not reading optional pointers (text etc).
    skipOptionalTrailingArguments(start, pd, rd);

    size_t ki;
    for (ki = 0; ki < local_geom_.size(); ++ki)
	if (geom_id_[ki] == pos_id)
	    break;
    DEBUG_ERROR_IF(ki==local_geom_.size(),
		   "Missing object. Could not find pos #" << pos_id);
    shared_ptr<PointCloud3D> pos_cl =
	dynamic_pointer_cast<PointCloud3D, GeomObject>(local_geom_[ki]);
    Point pos(pos_cl->rawData()[0], pos_cl->rawData()[1], pos_cl->rawData()[2]);
    geom_used_[ki] = 1;

    for (ki = 0; ki < local_geom_.size(); ++ki)
	if (geom_id_[ki] == axis_id)
	    break;
    DEBUG_ERROR_IF(ki==local_geom_.size(),
		   "Missing object. Could not find pos #" << axis_id);
    shared_ptr<PointCloud3D> axis_cl =
	dynamic_pointer_cast<PointCloud3D, GeomObject>(local_geom_[ki]);
    Point axis(axis_cl->rawData()[0],
	       axis_cl->rawData()[1], axis_cl->rawData()[2]);
    geom_used_[ki] = 1;

    if (form == 0) {
	// We need to set a dir vector, orthogonal to the axis of the cylinder.
	Point dir = ((axis[0] != 0.0) && (axis[1] != 0.0)) ?
	    Point(axis[1], -axis[0], 0.0) : Point(1.0, 0.0, 0.0);
	cylinder = shared_ptr<Cylinder>(new Cylinder(radius, pos, axis, dir));
    } else {
	for (ki = 0; ki < local_geom_.size(); ++ki)
	    if (geom_id_[ki] == ref_dir_id)
		break;
	DEBUG_ERROR_IF(ki==local_geom_.size(),
		       "Missing object. Could not find dir #" << ref_dir_id);
	shared_ptr<PointCloud3D> dir_cl =
	    dynamic_pointer_cast<PointCloud3D, GeomObject>(local_geom_[ki]);
	Point dir(dir_cl->rawData()[0],
		  dir_cl->rawData()[1], dir_cl->rawData()[2]);
	geom_used_[ki] = 1;
	cylinder = shared_ptr<Cylinder>(new Cylinder(radius, pos, axis, dir));
    }

    return cylinder;
}


//-----------------------------------------------------------------------------
shared_ptr<SplineCurve>
IGESconverter::readIGEScircularsegment(const char* start,
				       int num_lines,
				       int direntry_index)
//-----------------------------------------------------------------------------
{
    char pd = header_.pardel;
    char rd = header_.recdel;

    // Verify that start points to '100,.....' warn if not
    int type = readIGESint(start, pd, rd);
    DEBUG_ERROR_IF(type!=100, "Entity is not of type 100.");
    skipDelimiter(start, pd);

    // Read the data
    double zt;
    Array<double, 3> center;
    Array<double, 3> seg1;
    Array<double, 3> seg2;
    zt = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    center[0] = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    center[1] = readIGESdouble(start, pd, rd);
    center[2] = zt;
    skipDelimiter(start, pd);
    seg1[0] = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    seg1[1] = readIGESdouble(start, pd, rd);
    seg1[2] = zt;
    skipDelimiter(start, pd);
    seg2[0] = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    seg2[1] = readIGESdouble(start, pd, rd);
    //skipDelimiter(start, rd);
    seg2[2] = zt;

    skipOptionalTrailingArguments(start, pd, rd);

    // Extract the coordinate system
    int csentry = direntries_[direntry_index].trans_matrix;
    CoordinateSystem<3> cs; 
    if (csentry != 0) { // If value of directory entry is 0, we should
			// use identity.
	map< int, CoordinateSystem<3> >::iterator it
	    = coordsystems_.find(csentry);
	if (it == coordsystems_.end()) {
	    MESSAGE("Could not find the referred coordinate system ("
		    << csentry << ") in the file. Using identity.");
	} else {
	    cs = it->second;
	}
    }
    double angle = (seg1-center).angle(seg2-center);
    if (angle > -1e-6 && angle < 1e-6) angle = 2*M_PI;
    if (angle < 0) angle += 2*M_PI;
    Array<double, 3> c_tr = cs*center;
    Array<double, 3> s_tr = cs*seg1;
//     Array<double, 3> axis_tr = cs*Array<double, 3>(0.0, 0.0, 1.0);
    Array<double, 3> axis_tr = cs.rot()*Array<double, 3>(0.0, 0.0, 1.0);
    SISLCurve* sisl_cv;
    int stat;
    s1303(s_tr.begin(), 1e-6, angle, c_tr.begin(),
	  axis_tr.begin(), 3, &sisl_cv, &stat);
    if (stat < 0) {
	THROW("Could not make circular arc curve. s1303 failed with error " << stat << ".");
    }
    shared_ptr<SplineCurve> crv(SISLCurve2Go(sisl_cv));
    freeCurve(sisl_cv);

    return crv;
}


//-----------------------------------------------------------------------------
shared_ptr< CoordinateSystem<3> >
IGESconverter::readIGEStransformation(const char* start,
				      int num_lines)
//-----------------------------------------------------------------------------
{
    char pd = header_.pardel;
    char rd = header_.recdel;

    // Verify that start points to '124,.....' warn if not
    int type = readIGESint(start, pd, rd);
    DEBUG_ERROR_IF(type!=124, "Entity is not of type 124.");
    skipDelimiter(start, pd);

    shared_ptr< CoordinateSystem<3> > cs(new CoordinateSystem<3>);
    // Read the data
    cs->rot()(0, 0) = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    cs->rot()(0, 1) = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    cs->rot()(0, 2) = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    cs->tr()[0] = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    cs->rot()(1, 0) = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    cs->rot()(1, 1) = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    cs->rot()(1, 2) = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    cs->tr()[1] = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    cs->rot()(2, 0) = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    cs->rot()(2, 1) = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    cs->rot()(2, 2) = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    cs->tr()[2] = readIGESdouble(start, pd, rd);
    //skipDelimiter(start, rd);

    skipOptionalTrailingArguments(start, pd, rd);

    return cs;
}


//-----------------------------------------------------------------------------
void
IGESconverter::readIGEScompositeCurve(const char* start, int num_lines,
				      int direntry_index, vector<int>& crv_vec)
//-----------------------------------------------------------------------------
{
    char pd = header_.pardel;
    char rd = header_.recdel;

    // Verify that start points to '102,.....' warn if not
    int type = readIGESint(start, pd, rd);
    // This does not handle very general cases
    //    MESSAGE("Structure is yet to be determined.");
    DEBUG_ERROR_IF(type!=102, "Entity is not of type 102.");
    skipDelimiter(start, pd);

    // Read parameters
    int nmb_ent = readIGESint(start, pd, rd); // Number of entities

    vector<int> ent_id(nmb_ent);
    crv_vec.resize(nmb_ent);

    int i;
    for (i = 0; i < nmb_ent; i++) {
	skipDelimiter(start, pd);
	ent_id[i] = readIGESint(start, pd, rd);

	// Find curve
	size_t h;
	for (h=0; h<local_geom_.size(); ++h)
	    if (geom_id_[h] == ent_id[i])
		break;
	DEBUG_ERROR_IF(h==local_geom_.size(),
		       "Missing object. Could not find curve #" << ent_id[i]);
	geom_used_[h] = 1;
	crv_vec[i] = (int)h;
    }
    //skipDelimiter(start, rd);

    skipOptionalTrailingArguments(start, pd, rd);

    // Extract the coordinate system
    int csentry = direntries_[direntry_index].trans_matrix;
    CoordinateSystem<3> cs; 
    if (csentry != 0) { // If value of directory entry is 0, we should
			// use identity.
	map< int, CoordinateSystem<3> >::iterator it
	    = coordsystems_.find(csentry);
	if (it == coordsystems_.end()) {
	    MESSAGE("Could not find the referred coordinate system ("
		    << csentry << ") in the file. Using identity.");
	} else {
	  MESSAGE("Transformation matrix for composite curve "
		  "object is missing!");
	    cs = it->second;
	}
    }
}


//-----------------------------------------------------------------------------
shared_ptr<SplineCurve>
IGESconverter::readIGESlinearPath(const char* start, int num_lines, int form)
//-----------------------------------------------------------------------------
{
    char pd = header_.pardel;
    char rd = header_.recdel;

    ASSERT(form == 12);

    // Verify that start points to '106,.....' warn if not
    int type = readIGESint(start, pd, rd);
    DEBUG_ERROR_IF(type!=106, "Entity is not of type 106.");
    skipDelimiter(start, pd);

    int ip = readIGESint(start, pd, rd);
    DEBUG_ERROR_IF(ip!=2, "Expecting 3d-points as form is 12. Got 2d-points!");
    skipDelimiter(start, pd);

    int nmb_pts = readIGESint(start, pd, rd);

    int dim = 3;
    vector<double> knots(nmb_pts + 2);
    knots[0] = 0.0;
    vector<double> coefs(dim*nmb_pts);
    for (int ki = 0; ki < nmb_pts; ++ki) {
	skipDelimiter(start, pd);
	coefs[dim*ki] = readIGESdouble(start, pd, rd);
	skipDelimiter(start, pd);
	coefs[dim*ki+1] = readIGESdouble(start, pd, rd);
	skipDelimiter(start, pd);
	coefs[dim*ki+2] = readIGESdouble(start, pd, rd);
	knots[ki+1] = (double)ki;
    }
    knots[nmb_pts+1] = (double)nmb_pts - 1;
    //skipDelimiter(start, rd);

    skipOptionalTrailingArguments(start, pd, rd);

    return shared_ptr<SplineCurve>(new SplineCurve(nmb_pts, 2, knots.begin(),
						       coefs.begin(), dim));
}


//-----------------------------------------------------------------------------

void
IGESconverter::readIGESboundary(const char* start, int num_lines,
                                vector<shared_ptr<CurveOnSurface> >& crv_vec)
//-----------------------------------------------------------------------------
{
//   SplineCurve *tmpc1;
//   SplineCurve *tmpc2;
    char pd = header_.pardel;
    char rd = header_.recdel;

    // Verify that start points to '141,.....' warn if not
    int type = readIGESint(start, pd, rd);
    DEBUG_ERROR_IF(type!=141, "Entity is not of type 141.");
    skipDelimiter(start, pd);

    // Read parameters
    int parcrv_exist;
    parcrv_exist = readIGESint(start, pd, rd);
    skipDelimiter(start, pd);

    int pref = readIGESint(start, pd, rd);
    skipDelimiter(start, pd);

    int surf_id = readIGESint(start, pd, rd);  // Untrimmed surface
    skipDelimiter(start, pd);
        // Look for the surface
    shared_ptr<ParamSurface> surf;
    int i=0, j=0;
    size_t h=0;
    for (i=0; i<int(local_geom_.size()); i++)
      if (geom_id_[i] == surf_id)
        break;

    DEBUG_ERROR_IF(i==int(local_geom_.size()), "Missing object. Could not find surface #" << surf_id);
    shared_ptr<GeomObject> lg = local_geom_[i];
    surf = dynamic_pointer_cast<ParamSurface, GeomObject>(lg);
    //std::cout << surf.get() << std::endl;
    geom_used_[i] = 1;

    int nmb_crv = readIGESint(start, pd, rd);  // Number of boundary entities

//     shared_ptr<SplineCurve> curr_spacecrv;
    shared_ptr<ParamCurve> curr_spacecrv;
    for (i=0; i<nmb_crv; i++)
    {
      skipDelimiter(start, pd);
      int space_id = readIGESint(start, pd, rd);  // Space curve identifier
      skipDelimiter(start, pd);

          // Find space curve
      for (h=0; h<local_geom_.size(); h++)
        if (geom_id_[h] == space_id)
          break;
      DEBUG_ERROR_IF(h==local_geom_.size(), "Missing object. Could not find curve #" << space_id);
      lg = local_geom_[h];
      curr_spacecrv = dynamic_pointer_cast<ParamCurve, GeomObject>(lg);
      geom_used_[h] = 1;

          // Check if the space curve must be turned
      int sense =  readIGESint(start, pd, rd);  // Space curve identifier
      skipDelimiter(start, pd);
      if (sense == 2)
      {
	  curr_spacecrv = shared_ptr<ParamCurve>
	      (dynamic_cast<ParamCurve*>(local_geom_[h]->clone()));
        curr_spacecrv->reverseParameterDirection();
      }
      
      int nmb_par = readIGESint(start, pd, rd); // Number of parameter curves
      if (nmb_par == 0)
      {
	  CurveOnSurface *srf_crv = new CurveOnSurface(surf,
						       curr_spacecrv, 0);
        crv_vec.push_back(shared_ptr<CurveOnSurface>(srf_crv));
      }
      else
      {
	skipDelimiter(start, pd);
        int par_id = readIGESint(start, pd, rd);  // Parameter curve id.;
        double par1 = curr_spacecrv->startparam();
        double par3 = curr_spacecrv->endparam();
        double par2;
        shared_ptr<ParamCurve> spacecrv = curr_spacecrv;
// 	SplineCurve *pc;
        for (j=1; j<nmb_par; j++, par1=par2)
        {
	  skipDelimiter(start, pd);

	  // Find parameter curve
	  for (h=0; h<local_geom_.size(); h++)
	    if (geom_id_[h] == par_id)
	      break;
          
	  DEBUG_ERROR_IF(h==local_geom_.size(), "Missing object. Could not find curve #" << par_id);
	  geom_used_[h] = 1;

              // Fetch the appropriate piece of the corresponding space
              // curve. First find the endparameter of the segment.
              // Evaluate the surface in the endparameter of the
              // current parameter curve
          shared_ptr<GeomObject> lg = local_geom_[h];
          shared_ptr<ParamCurve> pcurve =
	      dynamic_pointer_cast<ParamCurve, GeomObject>(lg);

	  // Check if the curve is planar
// 	  int pind = pnumber_to_plane_normal_index_[Pnumber_[h]];
	  int pind = pnumber_to_plane_normal_index_[geom_id_[h]];
	  if (plane_normal_[pind].dimension() == 0)
	    // Assume that the curve lies in the xy-plane
	    plane_normal_[pind].setValue(0.0, 0.0, 1.0);

	  // Fetch the planar curve.
// 	  pc = dynamic_pointer_cast<SplineCurve, ParamCurve>(pcurve).get();
// 	  pcurve = GeometryTools::projectCurve(*pc, plane_normal_[pind], true);
	  pcurve = GeometryTools::projectCurve(pcurve, plane_normal_[pind], true);

          double ppar = pcurve->endparam();
          Point parpnt = pcurve->point(ppar);
          Point pntonsrf = surf->ParamSurface::point(parpnt[0], parpnt[1]);

              // Find closest point on the space curve.
          Point clopnt;
          double clodist;
          curr_spacecrv->closestPoint(pntonsrf, par1, par3, par2, clopnt,
                                      clodist);

	  // Fetch piece of space curve
  	  spacecrv = shared_ptr<ParamCurve>( 
  		      (curr_spacecrv->subCurve(par1, par2))); 
        
	  // Make curve-on-surface curve
          CurveOnSurface *srf_crv = new CurveOnSurface(surf,
						       pcurve,
						       spacecrv,
						       (pref != 1));
          crv_vec.push_back(shared_ptr<CurveOnSurface>(srf_crv));
          
          par_id = readIGESint(start, pd, rd);  // Next parameter curve id.
        }

	// Find parameter curve
	for (h=0; h<local_geom_.size(); h++)
	  if (geom_id_[h] == par_id)
	    break;
          
	DEBUG_ERROR_IF(h==local_geom_.size(), "Missing object. Could not find curve #" << par_id);
	geom_used_[h] = 1;
            // Fetch piece of space curve
// 	spacecrv = shared_ptr<SplineCurve>
// 	    (new SplineCurve(*(curr_spacecrv->subCurve(par1, par3))));
        spacecrv = curr_spacecrv;

            // Make curve-on-surface curve
        shared_ptr<GeomObject> lg = local_geom_[h];
        shared_ptr<ParamCurve> pcurve =
	    dynamic_pointer_cast<ParamCurve, GeomObject>(lg);

	// Check if the curve is planar
// 	int pind = pnumber_to_plane_normal_index_[Pnumber_[h]];
	int pind = pnumber_to_plane_normal_index_[geom_id_[h]];
	if (plane_normal_[pind].dimension() == 0)
	  // Assume that the curve lies in the xy-plane
	  plane_normal_[pind].setValue(0.0, 0.0, 1.0);

	// Fetch the planar curve.
// 	pc = dynamic_pointer_cast<SplineCurve, ParamCurve>(pcurve).get();
// 	pcurve = projectCurve(*pc, plane_normal_[pind], true);
	pcurve = GeometryTools::projectCurve(pcurve, plane_normal_[pind], true);

// 	tmpc1 = dynamic_cast<SplineCurve*>(pcurve.get());
// 	tmpc2 = dynamic_cast<SplineCurve*>(spacecrv.get());
        CurveOnSurface *srf_crv = new CurveOnSurface(surf,
                                                         pcurve,
                                                         spacecrv,
                                                         (pref != 1));
        crv_vec.push_back(shared_ptr<CurveOnSurface>(srf_crv));
          
      }
    }
    skipOptionalTrailingArguments(start, pd, rd);
}



// A helper function for readIGEScurveOnSurf().
//-----------------------------------------------------------------------------
shared_ptr<CurveOnSurface>
IGESconverter::makeCurveOnSurface(int pc_ind,
				  int sc_ind,
				  shared_ptr<ParamSurface> surf,
				  bool pref_param, int ccm)
//-----------------------------------------------------------------------------
{
    shared_ptr<ParamCurve> pc, sc;
    if (pc_ind != -1) {
	pc = dynamic_pointer_cast<ParamCurve, GeomObject>(local_geom_[pc_ind]);
	// Project the parametric curve (usually, the z-coord is zero, but we'd better
	// make it safe).
// 	int pind = pnumber_to_plane_normal_index_[Pnumber_[pc_ind]];
 	int pind = pnumber_to_plane_normal_index_[geom_id_[pc_ind]];
	if (plane_normal_[pind].dimension() == 0) {
	    // Assume that the curve lies in the xy-plane
	    plane_normal_[pind].setValue(0.0, 0.0, 1.0);
	}
// 	SplineCurve* pcspline
// 	    = dynamic_cast<SplineCurve*>(pc.get());
// 	pc = projectCurve(*pcspline, plane_normal_[pind], true);
	pc = GeometryTools::projectCurve(pc, plane_normal_[pind], true);
    }
    if (sc_ind != -1) {
	sc = dynamic_pointer_cast<ParamCurve, GeomObject>(local_geom_[sc_ind]);
    }

    shared_ptr<CurveOnSurface> cv;
    if (pc_ind != -1 && sc_ind != -1) {
	cv.reset(new CurveOnSurface(surf, pc, sc, pref_param, ccm));
    } else if (pc_ind != -1) {
	cv.reset(new CurveOnSurface(surf, pc, true));
    } else if (sc_ind != -1) {
	cv.reset(new CurveOnSurface(surf, sc, false));
    } else {
	THROW ("No valid indices for either param curve or space curve: "
	       << pc_ind << ", " << sc_ind);
    }
    return cv;
}

// Entity defining the curves may be given as composite curves.
// Our result is hence stored in a local vector.
//-----------------------------------------------------------------------------
void
IGESconverter::readIGEScurveOnSurf(const char* start, int num_lines,
				   vector<shared_ptr<CurveOnSurface> >& crv_vec)
//-----------------------------------------------------------------------------
{
    crv_vec.clear();
    bool pc_comp_curve = false; // By default we assume curve is a NURB.
    bool sc_comp_curve = false; // By default we assume curve is a NURB.

    char pd = header_.pardel;
    char rd = header_.recdel;

    // Verify that start points to '142,.....' warn if not
    int type = readIGESint(start, pd, rd);
    DEBUG_ERROR_IF(type!=142,
	       "Entity is not of type 142.");
    skipDelimiter(start, pd);

    // Read parameters.
    int ccm = readIGESint(start, pd, rd); // Curve creation method:
					  // 0=undefined, 1=proj,
					  // 2=sf-sf-int, 3=isopar-cv
    skipDelimiter(start, pd);

    // Look for the surface.
    int surf_id = readIGESint(start, pd, rd);  // Untrimmed surface
    skipDelimiter(start, pd);
    int i=0;
    size_t h=0, g=0;
    for (i=0; i<int(local_geom_.size()); i++)
	if (geom_id_[i] == surf_id)
	    break;
    if (i==int(local_geom_.size()))
	THROW("Missing object. Could not find surface #" << surf_id);
    shared_ptr<GeomObject> lg = local_geom_[i];
    geom_used_[i] = 1;
    shared_ptr<ParamSurface> surf =
	dynamic_pointer_cast<ParamSurface, GeomObject>(lg);

    // The surface might be a plane. In that case it is necessary to
    // make a planar representation
//     shared_ptr<Plane> plane;
//     // Check also for plane
//     plane = dynamic_pointer_cast<Plane, GeomObject>(lg);

    // Find parameter curve.
    int par_id = readIGESint(start, pd, rd);  // Parameter curve id.;
					      // may be 0 (cv not
					      // given)
    skipDelimiter(start, pd);
    int pc_ind_spline = -1;
    int pc_ind_comp = -1;
    if (par_id != 0) {
	for (h=0; h<local_geom_.size(); h++) {
	    if (geom_id_[h] == par_id) {
		pc_ind_spline = (int)h;
		break;
	    }
	}
	if (h==local_geom_.size()) {
	    for (g=0; g<local_comp_curve_.size(); g++) {
		if (comp_curve_id_[g] == par_id) {
		    pc_ind_comp = (int)g;
		    break;
		}
	    }
	    DEBUG_ERROR_IF(g==local_geom_.size(),
			   "Missing object. Could not find curve #" << par_id);
	    pc_comp_curve = true;
	}
    }
    // Now, either local_geom_[pc_ind_spline] or
    // local_comp_curve_[pc_ind_comp] is the parameter curve.
    // The bool comp_curve indicates which.

    // Now, we do the same for the space curve.
    int space_id = readIGESint(start, pd, rd);  // Space curve
						// identifier; may be
						// 0 (cv not given).
    skipDelimiter(start, pd);
    DEBUG_ERROR_IF((par_id == 0) && (space_id == 0),
	     "Missing both parameter curve and space curve!");
    int sc_ind_spline = -1;
    int sc_ind_comp = -1;
    if (space_id != 0) {
	for (h=0; h<local_geom_.size(); h++) {
	    if (geom_id_[h] == space_id) {
		sc_ind_spline = (int)h;
		break;
	    }
	}
	if (h==local_geom_.size()) {
	    for (g=0; g<local_comp_curve_.size(); g++) {
		if (comp_curve_id_[g] == space_id) {
		    sc_ind_comp = (int)g;
		    break;
		}
	    }
	    DEBUG_ERROR_IF(g==local_geom_.size(), "Missing object. Could not find curve #" << space_id);
	    sc_comp_curve = true;
	}
    }

//     if (plane.get() != 0)
//     {
// 	// Create trimmed plane. Use the space curves to limit the planar surface
// 	vector<shared_ptr<ParamCurve> > space_crvs;
// 	if (!sc_comp_curve)
// 	    space_crvs.push_back(dynamic_pointer_cast<ParamCurve, GeomObject>
// 				 (local_geom_[sc_ind_spline]));
// 	else
// 	{
// 	   int num_crv = local_comp_curve_[sc_ind_comp].size(); 
// 	   for (int i=0; i<num_crv; ++i)
// 	     space_crvs.push_back(dynamic_pointer_cast<ParamCurve, GeomObject>
// 				  (local_geom_[local_comp_curve_[sc_ind_comp][i]]));  
// 	}
// 	surf = makeTrimmedPlane(plane, space_crvs);
//     }

#ifdef SBR_DBG
    std::ofstream outfile("tmp/curr_bd_sf.g2");
    surf->writeStandardHeader(outfile);
    surf->write(outfile);

    vector<shared_ptr<ParamCurve> > space_crvs;
    if (space_id != 0)
	if (!sc_comp_curve)
	    space_crvs.push_back(dynamic_pointer_cast<ParamCurve, GeomObject>
				 (local_geom_[sc_ind_spline]));
	else
	{
	    int num_crv = local_comp_curve_[sc_ind_comp].size(); 
	    for (int i=0; i<num_crv; ++i)
		space_crvs.push_back(dynamic_pointer_cast<ParamCurve, GeomObject>
				     (local_geom_[local_comp_curve_[sc_ind_comp][i]]));  
	}

    // @@@ VSK 0209. This is a hack. It seems that either some iges
    // files have a misplaced plane or the plane definition is
    // buggy. Move the plane position to fit with the trimming curve
    if (surf->instanceType() == Class_Plane)
    {
	shared_ptr<Plane> plane = dynamic_pointer_cast<Plane>(surf);
	BoundedUtils::translatePlaneToCurves(plane, space_crvs);
    }

    for (size_t kj = 0; kj < space_crvs.size(); ++kj) {
	space_crvs[kj]->writeStandardHeader(outfile);
	space_crvs[kj]->write(outfile);
    }
    vector<shared_ptr<ParamCurve> > par_crvs;
    if (par_id != 0)
	if (!sc_comp_curve)
	    par_crvs.push_back(dynamic_pointer_cast<ParamCurve, GeomObject>
			       (local_geom_[pc_ind_spline]));
	else
	{
	    int num_crv = local_comp_curve_[pc_ind_comp].size();
	    for (int i=0; i<num_crv; ++i)
		par_crvs.push_back
		    (dynamic_pointer_cast<ParamCurve, GeomObject>
		     (local_geom_[local_comp_curve_[pc_ind_comp][i]]));
	}
    for (size_t kj = 0; kj < par_crvs.size(); ++kj)
	if (par_crvs[kj] != NULL)
	    if (par_crvs[kj]->dimension() == 3) {
		par_crvs[kj]->writeStandardHeader(outfile);
		par_crvs[kj]->write(outfile);
	    } else {
		shared_ptr<SplineCurve> spline_par_cv;
		if (par_crvs[kj]->instanceType() == Class_SplineCurve)
		    spline_par_cv = dynamic_pointer_cast<SplineCurve, ParamCurve>
			(par_crvs[kj]);
		else
		    spline_par_cv = shared_ptr<SplineCurve>
			(par_crvs[kj]->geometryCurve());
// 		writeSpaceParamCurve(*par_crvs[kj], outfile, 0.0);
		writeSpaceParamCurve(*spline_par_cv, outfile, 0.0);
	    }
	else
	    MESSAGE("Unsupported curve type!");
    double debug_val = 0.0;
#endif

    // Curve preference.
    int pref = readIGESint(start, pd, rd);  // Indicates preferance in
					    // send syst: 0=not_set,
					    // 1=par, 2=space, 3=equal.

    skipOptionalTrailingArguments(start, pd, rd);

    // Param and space curves must either both be composite, or both splines.
    if (par_id != 0 && space_id != 0 && pc_comp_curve != sc_comp_curve) {
	THROW("Only one of the (parametric and spatial) curves are composite.");
    }

    // If non-composite, deal with it and exit.
    bool prefer_par = (pref != 2);
    if (par_id == 0)
	prefer_par = false;
    else if (space_id == 0)
	prefer_par = true;
    if (!pc_comp_curve && !sc_comp_curve) {
	// We're dealing with non-composite curves.
	crv_vec.push_back(makeCurveOnSurface(pc_ind_spline,
					     sc_ind_spline,
					     surf,
					     prefer_par, ccm));
	return;
    }

    // If composite, we check if they have the same number of subcurves.
    int num_crv_pc = (pc_ind_comp < 0) ? 0 : (int)local_comp_curve_[pc_ind_comp].size();
    int num_crv_sc = (sc_ind_comp < 0) ? 0 : (int)local_comp_curve_[sc_ind_comp].size();
    if (num_crv_pc == num_crv_sc || num_crv_pc == 0 || num_crv_sc == 0) {
	int num_crv = (prefer_par) ? num_crv_pc : num_crv_sc;
	// We assume that space and param curves correspond.
	for (int i = 0; i < num_crv; ++i) {
	    crv_vec.push_back(makeCurveOnSurface((num_crv_pc == 0) ? -1 : 
						 local_comp_curve_[pc_ind_comp][i],
						 (num_crv_sc == 0) ? -1 : 
						 local_comp_curve_[sc_ind_comp][i],
						 surf,
						 prefer_par, ccm));
	}
	return;
    }
	

    // Since we got here, the parametric and spatial
    // curve chains do not correspond directly (different size).
    if (num_crv_pc < num_crv_sc) {
	THROW("More spatial than parametric curves: " << num_crv_sc 
	      << " spatial vs. " << num_crv_pc << " parametric.");
    }
    // Probably, the extra parametric curves have a degenerate
    // image on the surface.
    // Find which parametric curves are degenerated.
    int numdegen = 0;
    // pstart and pend are 3D because the curves are, too.
    Point pstart(3);
    Point pend(3);
    Point sstart(3);
    Point send(3);
    vector<bool> deg(num_crv_pc, false);
    for (i = 0; i < num_crv_pc; ++i) {
	shared_ptr<GeomObject> obj
	    = local_geom_[local_comp_curve_[pc_ind_comp][i]];
	const ParamCurve& cv
	    = dynamic_cast<const ParamCurve&>(*obj);
	cv.point(pstart, cv.startparam());
	cv.point(pend, cv.endparam());
	// @ Assuming that curves are given in the XY-plane.
	surf->point(sstart, pstart[0], pstart[1]);
	surf->point(send, pend[0], pend[1]);
	double d = sstart.dist(send);
	if (d < header_.min_resolution) {
	    // Degenerate
	    deg[i] = true;
	    ++numdegen;
	}
    }
    if (numdegen != num_crv_pc - num_crv_sc) {
	THROW("Found " << numdegen << " degenerate curves, should be "
	      << num_crv_pc - num_crv_sc);
    }
    int cur_p = 0;
    int cur_sc = 0;
    // @@sbr I whould prefer to add also the degenerate cvs (as they may define a loop
    // in the parameter domain, and hence removing them may yield a non-closed loop).
    int total_nmb_cvs = num_crv_sc + numdegen;
    for (i = 0; i < total_nmb_cvs; ++i, ++cur_p) {
	// Find next non-degenerate parametric curve
// 	while (deg[cur_p]) ++cur_p;
	// Make the curve-on-surface
	if (deg[cur_p]) {
	    crv_vec.push_back(makeCurveOnSurface(local_comp_curve_[pc_ind_comp][cur_p],
						 -1,
						 surf,
						 pref != 2, ccm));
	} else {
	    crv_vec.push_back(makeCurveOnSurface(local_comp_curve_[pc_ind_comp][cur_p],
						 local_comp_curve_[sc_ind_comp][cur_sc],
						 surf,
						 pref != 2, ccm));
	    ++cur_sc;
	}
    }
}


// @afr: Commented out because the code was replaced.
    /*
    for (i=0; i<nmb_curves; i++) {
	// We fetch the parameter curve, if it exists.
	shared_ptr<ParamCurve> pcurve;
	if (par_id != 0) {
	    if (comp_curve) {
		lg = local_comp_curve_[g][i];
		for (h=0; h<local_geom_.size(); h++)
		    {
			shared_ptr<GeomObject> lg = local_geom_[h];
			if (local_comp_curve_[g][i] ==
			    dynamic_pointer_cast<ParamCurve, GeomObject>(lg))
			    break;
		    }
		ALWAYS_ERROR_IF(h==local_geom_.size(),
			    "This should never happen. Seems like a bug.",
			    InputError());
	    } else {
		lg = local_geom_[h];
		geom_used_[h] = 1;
	    }
	    // Check if the curve is planar
	    int pind = pnumber_to_plane_normal_index_[Pnumber_[h]];
	    if (plane_normal_[pind].dimension() == 0)
		// Assume that the curve lies in the xy-plane
		plane_normal_[pind].setValue(0.0, 0.0, 1.0);

	    pcurve = dynamic_pointer_cast<ParamCurve, GeomObject>(lg);

	    // Fetch the planar curve.
	    SplineCurve *pc;
	    pc = dynamic_pointer_cast<SplineCurve, ParamCurve>(pcurve).get();
	    pcurve = projectCurve(*pc, plane_normal_[pind], true);
	}

	shared_ptr<SplineCurve> spacecrv;
	if (space_id != 0) {
	    // Make space curve
	    if (comp_curve)
		lg = local_comp_curve_[j][i];
	    if (!comp_curve) {
		lg = local_geom_[j];
		geom_used_[j] = 1;
	    }
	    spacecrv = dynamic_pointer_cast<SplineCurve, GeomObject>(lg);
	}

	CurveOnSurface *crv_on_surf;
	if ((par_id != 0) && (space_id != 0)) {
	    // Make curve-on-surface curve
	    crv_on_surf = new CurveOnSurface(surf, pcurve, spacecrv, 1);
	} else {
	    if (par_id != 0)
		crv_on_surf = new CurveOnSurface(surf, pcurve, true);
	    else crv_on_surf = new CurveOnSurface(surf, spacecrv, false);
	}

	shared_ptr<CurveOnSurface> crv;
	crv = (shared_ptr<CurveOnSurface>)(crv_on_surf);
	crv_vec.push_back(crv);
    }
    int pref = readIGESint(start, pd, rd);  // Indicates preferance in send syst.
    skipDelimiter(start, rd);
}
*/


//-----------------------------------------------------------------------------
shared_ptr<BoundedSurface>
IGESconverter::readIGESboundedSurf(const char* start, int num_lines)
//-----------------------------------------------------------------------------
{
    char pd = header_.pardel;
    char rd = header_.recdel;

    // Verify that start points to '143,.....' warn if not
    int type = readIGESint(start, pd, rd);
    DEBUG_ERROR_IF(type!=143, "Entity is not of type 143.");
    skipDelimiter(start, pd);

    // Read parameters
    int parcrv_exist;
    parcrv_exist = readIGESint(start, pd, rd);
    skipDelimiter(start, pd);

    int surf_id = readIGESint(start, pd, rd);  // Untrimmed surface
    skipDelimiter(start, pd);

    int nmb_bd = readIGESint(start, pd, rd);  // Number of boundary entities

    vector<int> bd_id(nmb_bd);
    for (int i=0; i<nmb_bd; i++) {
      skipDelimiter(start, pd);
      bd_id[i] = readIGESint(start, pd, rd);  // Current boundary entities
    }

    skipOptionalTrailingArguments(start, pd, rd);


        // Set pointers to the underlying geometry entities
    shared_ptr<ParamSurface> surf;
    vector<vector<shared_ptr<CurveOnSurface> > > boundaries(nmb_bd);

        // First look for the surface
    size_t isz;
    for (isz=0; isz<local_geom_.size(); isz++)
      if (geom_id_[isz] == surf_id)
        break;

    DEBUG_ERROR_IF(isz==local_geom_.size(),
		   "Missing object. Could not find surface #" << surf_id);
    shared_ptr<GeomObject> lg = local_geom_[isz];
    surf = dynamic_pointer_cast<ParamSurface, GeomObject>(lg);
    geom_used_[isz] = 1;

        // Fetch the boundary loops
    for (int j=0; j<nmb_bd; j++)
    {
      size_t i;
      for (i=0; i<local_loop_.size(); i++)
      {
        if (loop_id_[i] == bd_id[j])
          break;
      }
      DEBUG_ERROR_IF(i==local_loop_.size(),
		     "Missing object. Could not find boundary #" << bd_id[j]);
      for (size_t h=0; h<local_loop_[i].size(); h++)
	boundaries[j].push_back(local_loop_[i][h]);
//        boundaries[j] = local_loop_[i];
    }

// #ifdef TRY_FIX_INPUT
//     // Check that the loops don't have gaps in them
//     // (that are too large).
//     double tolerance = header_.min_resolution;
//     for (int j=0; j<nmb_bd; j++) {
// 	double maxgap = Go::computeLoopGap(boundaries[j]);
// 	if (maxgap > tolerance) {
// 	    // Try to fix by rearranging the segments.
// 	    vector<int> perm;
// 	    vector<bool> flip;
// 	    orientCurves::orientCurves(boundaries[j], perm, flip,
// 			 tolerance, false);
// 	    // Making the new boundary vector
// 	    vector< shared_ptr<CurveOnSurface> > new_boundary;
// 	    new_boundary.reserve(boundaries[j].size());
// 	    for (size_t bi = 0; bi < boundaries[j].size(); ++bi) {
// 		new_boundary.push_back(boundaries[j][perm[bi]]);
// 		if (flip[bi]) {
// 		    new_boundary[bi]->reverseParameterDirection();
// 		}
// 	    }
// 	    boundaries[j].swap(new_boundary);
// 	    // We check if that helped.
// 	    maxgap = Go::computeLoopGap(boundaries[j]);
// 	    if (maxgap > tolerance) {
// 		THROW("Cannot fix boundary that does not form a loop.\n"
// 		      "IGES file probably buggy.");
// 	    }
// 	}
//     }
// #endif // TRY_FIX_INPUT

    // Make bounded surface
    shared_ptr<BoundedSurface> bdsurf =
	shared_ptr<BoundedSurface>(new BoundedSurface(surf, boundaries,
						      header_.min_resolution));
    return bdsurf;
}

//-----------------------------------------------------------------------------
shared_ptr<BoundedSurface>
// shared_ptr<SplineSurface>
IGESconverter::readIGEStrimmedSurf(const char* start, int num_lines)
//-----------------------------------------------------------------------------
{
    char pd = header_.pardel;
    char rd = header_.recdel;

    // Verify that start points to '144,.....' warn if not
    int type = readIGESint(start, pd, rd);
    DEBUG_ERROR_IF(type!=144, "Entity is not of type 144.");
    skipDelimiter(start, pd);

    // Read parameters
    int surf_id = readIGESint(start, pd, rd);  // Untrimmed surface
    skipDelimiter(start, pd);
    // First look for the surface
    size_t i;
    for (i=0; i<local_geom_.size(); i++)
	if (geom_id_[i] == surf_id)
	    break;
    DEBUG_ERROR_IF(i==local_geom_.size(),
		   "Missing object. Could not find surface #" << surf_id);
    shared_ptr<GeomObject> lg = local_geom_[i];
    geom_used_[i] = 1;
    shared_ptr<ParamSurface> surf =
	dynamic_pointer_cast<ParamSurface, GeomObject>(lg);
//     shared_ptr<SplineSurface> surf =
// 	dynamic_pointer_cast<SplineSurface, GeomObject>(lg);

    // The surface might be a plane. In that case it is necessary to
    // make a planar representation
//     shared_ptr<Plane> plane;

//     // Check also for plane
//     plane = dynamic_pointer_cast<Plane, GeomObject>(lg);

    int obt = readIGESint(start, pd, rd); // Outer boundary type. 0 if
                                          // outer bnd is parameter
                                          // domain, otherwise 1.
    if (obt == 0) {
	THROW("Trimmed surface has no outer boundary defined,"
	      " we cannot handle this.");
    }

    skipDelimiter(start, pd);

    int nmb_i_bd = readIGESint(start, pd, rd);  // Number of inner
						// boundary entities

    vector<vector<shared_ptr<CurveOnSurface> > > boundaries(nmb_i_bd + obt);

    double loop_tol = header_.min_resolution;
    if (loop_tol == 0) {
	MESSAGE("Tolerance given was zero. Using 0.05 instead.");
	loop_tol = 0.05;
    }

    int o_bd_id = -1; // Outer boundary; yet to be initialized
    // Get o_bd_id, given that outer boundary is different from the
    // parameter domain
    if (obt != 0) {
	skipDelimiter(start, pd);
	o_bd_id = readIGESint(start, pd, rd);

	// Find outer boundary curve(s)
	for (i=0; i<local_loop_.size(); ++i)
	    if (loop_id_[i] == o_bd_id)
		break;
	DEBUG_ERROR_IF(i==local_loop_.size(),
		       "Missing object. Could not find boundary #" << o_bd_id);

	vector<shared_ptr<CurveOnSurface> > o_bd_crvs;
	for (size_t j=0; j<local_loop_[i].size(); ++j) {
	    lg = local_loop_[i][j];
	    shared_ptr<CurveOnSurface> o_bd_crv =
		dynamic_pointer_cast<CurveOnSurface, GeomObject>(lg);
	    o_bd_crvs.push_back(o_bd_crv);
	}

	boundaries[0] = o_bd_crvs;
// 	if (plane.get() != 0 && o_bd_crvs.size() > 0)
// 	{
// 	    // The underlying surface is a plane. Fetch surface from
// 	    // trimming curve
// 	    surf = dynamic_pointer_cast<SplineSurface, ParamSurface>
// 	      (o_bd_crvs[0]->underlyingSurface());
// 	}
	
    }

    vector<int> i_bd_id(nmb_i_bd);
    for (int j=0; j<nmb_i_bd; ++j) {
	skipDelimiter(start, pd);
      i_bd_id[j] = readIGESint(start, pd, rd);  // Current boundary entities
    }

    skipOptionalTrailingArguments(start, pd, rd);

    for (int j=0; j<nmb_i_bd; ++j)
    {
	// Find inner boundary curve(s)
	for (i = 0; i < local_loop_.size(); ++i)
	    if (loop_id_[i] == i_bd_id[j])
		break;
	DEBUG_ERROR_IF(i==local_loop_.size(),
		       "Missing object. Could not find boundary #" <<
		       i_bd_id[j]);
	vector<shared_ptr<CurveOnSurface> > i_bd_crvs; // Plural!
	for (size_t h=0; h<local_loop_[i].size(); ++h) {
	    lg = local_loop_[i][h];
	    shared_ptr<CurveOnSurface> i_bd_crv =
		dynamic_pointer_cast<CurveOnSurface, GeomObject>(lg);
// 	    if (plane.get() != 0)
// 	    {
// 		// The surface belonging to the curve-on-surface
// 		// curves of the inner loops is not the same as the
// 		// surface corresponding to the outer loop. Set
// 		// surface of inner loops
// 		i_bd_crv->setUnderlyingSurface(surf);
// 	    }
	    i_bd_crvs.push_back(i_bd_crv);
	}

	// i_bd_crvs should constitute a loop. Otherwise the BoundedSurface
	// constructor will fail.
	boundaries[obt+j] = i_bd_crvs;
    }

    // Make bounded surface
//     shared_ptr<BoundedSurface> trimsurf;
// //     shared_ptr<SplineSurface> under_surf;
// #ifdef TRY_FIX_INPUT
//     try {
// 	// Check that the loops don't have gaps in them
// 	// (that are too large).
// 	int nmb_bd = boundaries.size();
// 	double tolerance = loop_tol;
// 	for (int j=0; j<nmb_bd; j++) {
// 	    if (boundaries[j].size() == 0)
// 		continue;   // No curve in boundary loop

// 	    double maxgap = Go::computeLoopGap(boundaries[j]);
// 	    if (maxgap > tolerance) {
// 		// We first check if there exists a space curve with
// 		// legal definition (or parametric of we prefer space).
// 		// boundaries is a vector of CurveOnSurface.
// 		double maxgap_par = computeLoopGap(boundaries[j], true);
// 		double maxgap_space = computeLoopGap(boundaries[j], false);
// 		std::cout << "Loop number " << j << std::endl;
// 		std::cout << "Max gap par: " << maxgap_par << std::endl;
// 		std::cout << "Max gap space: " << maxgap_space << std::endl;
// 		if ((maxgap_par >= 0.0) &&
// 		    (maxgap_par < tolerance)) {
// 		    for (int k = 0; k < boundaries[j].size(); ++k)
// 		      {
// 			if (!boundaries[j][k]->makeCurvesConsistent(tolerance,
// 								    true))
// 			  boundaries[j][k] =
// 			    shared_ptr<CurveOnSurface>
// 			    (new CurveOnSurface
// 			     (boundaries[j][k]->underlyingSurface(),
// 			      boundaries[j][k]->parameterCurve(),
// 			      boundaries[j][k]->spaceCurve(),
// 			      true));
// 		      }
// 		    std::cout << "Loop rescued!" << std::endl;
// 		} else if ((maxgap_space >= 0.0) &&
// 			   (maxgap_space < tolerance)) {
// 		    for (int k = 0; k < boundaries[j].size(); ++k)
// 		      {
// 			if (!boundaries[j][k]->makeCurvesConsistent(tolerance,
// 								    false))
// 			  boundaries[j][k] =
// 			    shared_ptr<CurveOnSurface>
// 			    (new CurveOnSurface
// 			     (boundaries[j][k]->underlyingSurface(),
// 			      boundaries[j][k]->parameterCurve(),
// 			      boundaries[j][k]->spaceCurve(),
// 			      false));
// 		      }
// 		    std::cout << "Loop rescued!" << std::endl;
// 		} else {
// 		    // Try to fix by rearranging the segments.
// 		    vector<int> perm;
// 		    vector<bool> flip;
// 		    orientCurves::orientCurves(boundaries[j], perm, flip,
// 				 tolerance, false);
// 		    // Making the new boundary vector
// 		    vector< shared_ptr<CurveOnSurface> > new_boundary;
// 		    new_boundary.reserve(boundaries[j].size());
// 		    for (size_t bi = 0; bi < boundaries[j].size(); ++bi) {
// 			new_boundary.push_back(boundaries[j][perm[bi]]);
// 			if (flip[bi]) {
// 			    new_boundary[bi]->reverseParameterDirection();
// 			}
// 		    }
// 		    boundaries[j].swap(new_boundary);
// 		    // We check if that helped.
// 		    maxgap = Go::computeLoopGap(boundaries[j]);
// 		    if (maxgap > tolerance) {
// 			cerr << "Gap > tolerance: " << maxgap << " > "
// 			     << tolerance << endl;
// 			THROW("Cannot fix boundary that does not form a loop.\n"
// 			      "IGES file probably buggy.");
// 		    }
// 		}
// 	    }
// 	}
//  	trimsurf.reset(new BoundedSurface(surf, boundaries, loop_tol));
// // 	under_surf = surf;
//     }
//     catch (...) {
// #ifdef IGESLIB_DEBUG
// 	ofstream os("debug_trimsurf.g2");
// 	for (i = 0; i < boundaries[0].size(); ++i) {
// 	    if (boundaries[0][i]->spaceCurve().get() != 0) {
// 		boundaries[0][i]->spaceCurve()->writeStandardHeader(os);
// 		boundaries[0][i]->spaceCurve()->write(os);
// 	    }
// // 	    boundaries[0][i]->parameterCurve()->writeStandardHeader(os);
// // 	    boundaries[0][i]->parameterCurve()->write(os);
// 	}
// 	boundaries[0][0]->underlyingSurface()->writeStandardHeader(os);
// 	boundaries[0][0]->underlyingSurface()->write(os);
// #endif // IGESLIB_DEBUG
// 	std::cout << "Something went wrong with making the trimmed surface,"
// 	    " trying alternative approach." << std::endl;
// 	bool sf_rescued = true;
// 	int nmb_bd = boundaries.size();
// 	std::cout << "Number of loops: " << nmb_bd << std::endl;
// 	double tolerance = loop_tol;
// 	for (int j=0; j<nmb_bd; j++) {
// 	    if (boundaries[j].size() == 0)
// 		continue;   // No curve in boundary loop
// 	    double maxgap_par = computeLoopGap(boundaries[j], true);
// 	    double maxgap_space = computeLoopGap(boundaries[j], false);
// 	    std::cout << "Max gap par: " << maxgap_par << std::endl;
// 	    std::cout << "Max gap space: " << maxgap_space << std::endl;
// 	    if ((maxgap_space >= 0.0) &&
// 		(maxgap_space < tolerance)) {
// 		for (int k = 0; k < boundaries[j].size(); ++k)
// 		    boundaries[j][k] =
// 			shared_ptr<CurveOnSurface>
// 			(new CurveOnSurface
// 			 (boundaries[j][k]->underlyingSurface(),
// 			  boundaries[j][k]->parameterCurve(),
// 			  boundaries[j][k]->spaceCurve(),
// 			  false));
// 		std::cout << "Loop rescued!" << std::endl;
// 	    } else if ((maxgap_par >= 0.0) &&
// 		       (maxgap_par < tolerance)) {
// 		for (int k = 0; k < boundaries[j].size(); ++k)
// 		    boundaries[j][k] =
// 			shared_ptr<CurveOnSurface>
// 			(new CurveOnSurface
// 			 (boundaries[j][k]->underlyingSurface(),
// 			  boundaries[j][k]->parameterCurve(),
// 			  boundaries[j][k]->spaceCurve(),
// 			  true));
// 		std::cout << "Loop rescued!" << std::endl;
// 	    }  else {
// 		sf_rescued = false;
// 		break;
// 	    }
// 	}
// 	if (sf_rescued) {
// 	    std::cout << "Surface rescued!" << std::endl;
// 	    try {
// 		trimsurf.reset(new BoundedSurface(surf, boundaries, loop_tol));
// 	    } catch (...) {
// 		THROW("Something went wrong with making the trimmed surface.");
// 	    }
// 	} else {
// 	    THROW("Something went wrong with making the trimmed surface.");
// 	}
//     }
// #else
//     trimsurf.reset(new BoundedSurface(surf, boundaries, loop_tol));
// #endif // TRY_FIX_INPUT

//#ifndef NDEBUG
    // We write to file the surf and boundaries.
    std::ofstream debug_file("tmp/debug.g2");
    surf->writeStandardHeader(debug_file);
    surf->write(debug_file);
    for (size_t kj = 0; kj < boundaries.size(); ++kj)
	for (size_t ki = 0; ki < boundaries[kj].size(); ++ki)
	{
	    shared_ptr<ParamCurve> par_cv = boundaries[kj][ki]->parameterCurve();
	    shared_ptr<SplineCurve> par_scv = dynamic_pointer_cast<SplineCurve>(par_cv);
	    if (par_scv.get() != NULL)
		SplineDebugUtils::writeSpaceParamCurve(*par_scv, debug_file, 0.0);

	    shared_ptr<ParamCurve> space_cv = boundaries[kj][ki]->spaceCurve();
	    if (space_cv.get() != 0)
	    {
		space_cv->writeStandardHeader(debug_file);
		space_cv->write(debug_file);
	    }
	}
    //#endif

    shared_ptr<BoundedSurface> trimsurf
	(new BoundedSurface(surf, boundaries, loop_tol));

    std::ofstream deb2("tmp/bd_surf.g2");
    trimsurf->writeStandardHeader(deb2);
    trimsurf->write(deb2);

    return trimsurf;
//     return under_surf;
}

//-----------------------------------------------------------------------------
shared_ptr<SplineSurface>
IGESconverter::readIGESruledSurface(const char* start, int num_lines, int form)
//-----------------------------------------------------------------------------
{
    // MESSAGE("@jbt: Ruled surface - entity 118...");

    shared_ptr<SplineSurface> srf;

    if (form == 0) {
        MESSAGE("\nRuled Surface: Form is 0 (arc length parametrization)\n"
                << "Not implemented - using Form 1 (given parametrizations)"
                << "instead...");
    }
    else if (form != 1) {
        MESSAGE("What? Form not 0 or 1? Undefined!");
        return srf;
    }

    // The surface is defined as:
    // sf(u,v) = (1-v)*c1(t) + v*c2(s),
    // for 0 <= u,v <= 1.
    // We will not bother to define s and t here (see the documentation),
    // but t goes "along" the u direction and s goes either "along" or
    // "the opposite way", depending on the value of the flag DIRFLG (0 or
    // 1, respectively).
    // Another flag, DEVFLG, is 1 if sf is known to be developable, 0 if it
    // is not known.

    char pd = header_.pardel;
    char rd = header_.recdel;

    // Verify that start points to '118,.....' warn if not
    int type = readIGESint(start, pd, rd);
    DEBUG_ERROR_IF(type!=118, "Entity is not of type 118.");
    skipDelimiter(start, pd);

    // Read parameters
    int c1_id = readIGESint(start, pd, rd);  // First curve entity.
    skipDelimiter(start, pd);
    int c2_id = readIGESint(start, pd, rd);  // Second curve entity.
    skipDelimiter(start, pd);
    int dirflg = readIGESint(start, pd, rd);  // Direction flag
    skipDelimiter(start, pd);
    int devflg;
    devflg = readIGESint(start, pd, rd);  // Developable flag

    skipOptionalTrailingArguments(start, pd, rd);

    // Get the two curves
    shared_ptr<SplineCurve> sc1, sc2;
    bool found_c1 = false;
    bool found_c2 = false;
    for (size_t i = 0; i < local_geom_.size(); ++i) {
        if (geom_id_[i] == c1_id) {
            shared_ptr<ParamCurve> pc1 
                = dynamic_pointer_cast<ParamCurve>(local_geom_[i]);
            sc1 = shared_ptr<SplineCurve>(pc1->geometryCurve());
            if (sc1.get() == NULL) {
                MESSAGE("First curve does not have a spline representation!");
                return srf;
            }
            found_c1 = true;
        }
        if (geom_id_[i] == c2_id) {
            shared_ptr<ParamCurve> pc2 
                = dynamic_pointer_cast<ParamCurve>(local_geom_[i]);
            sc2 = shared_ptr<SplineCurve>(pc2->geometryCurve());
            if (sc2.get() == NULL) {
                MESSAGE("Second curve does not have a spline representation!");
                return srf;
            }
            found_c2 = true;
        }
        if (found_c1 && found_c2)
            break;
    }
    if (!found_c1) {
        MESSAGE("Could not find first curve");
        return srf;
    }
    if (!found_c2) {
        MESSAGE("Could not find second curve");
        return srf;
    }

    // If the Form number is 0, we should reparametrize at this point...
    if (form == 0) {
        // Do nothing...
    }

    // If direction flag is 1, we turn the direction of the second curve
    if (dirflg == 1) {
        shared_ptr<SplineCurve> tmp(sc2->clone());
        tmp->reverseParameterDirection();
        sc2 = tmp;
    }

    // We then use lofting to connect the curves.
    vector<shared_ptr<SplineCurve> > cvs(2);
    cvs[0] = sc1;
    cvs[1] = sc2;
    srf = shared_ptr<SplineSurface>(CoonsPatchGen::loftSurface(cvs.begin(), 2));

    // We set the parameter domain to the expected one according to
    // the IGES specification
    srf->setParameterDomain(0.0, 1.0, 0.0, 1.0);

// #ifndef NDEBUG
//     std::ofstream file_out("debug_ruled_surface.g2");
//     srf->writeStandardHeader(file_out);
//     srf->write(file_out);
//     sc1->writeStandardHeader(file_out);
//     sc1->write(file_out);
//     sc2->writeStandardHeader(file_out);
//     sc2->write(file_out);
// #endif // NDEBUG

    return srf;
}

//-----------------------------------------------------------------------------
shared_ptr<SplineSurface>
IGESconverter::readIGESsurfOfRevolution(const char* start, int num_lines)
//-----------------------------------------------------------------------------
{
    shared_ptr<SplineSurface> srf;

    char pd = header_.pardel;
    char rd = header_.recdel;

    // Verify that start points to '120,.....' warn if not
    int type = readIGESint(start, pd, rd);
    DEBUG_ERROR_IF(type!=120, "Entity is not of type 120.");
    skipDelimiter(start, pd);

    // Read parameters
    int axis_id = readIGESint(start, pd, rd);  // Untrimmed surface
    skipDelimiter(start, pd);
    int generatrix_id = readIGESint(start, pd, rd);
    skipDelimiter(start, pd);
    double start_angle = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    double end_angle = readIGESdouble(start, pd, rd);

    skipOptionalTrailingArguments(start, pd, rd);

    // Locate the axis
    size_t i;
    for (i=0; i<local_geom_.size(); i++)
	if (geom_id_[i] == axis_id)
	    break;
//     const SplineCurve& axis
// 	= dynamic_cast<const SplineCurve&>(*local_geom_[i]);
    shared_ptr<ParamCurve> par_cv =
	dynamic_pointer_cast<ParamCurve>(local_geom_[i]);
    shared_ptr<SplineCurve> axis(par_cv->geometryCurve());
//     if (local_geom_[i]->instanceType() == Class_SplineCurve)
//       axis = dynamic_pointer_cast<SplineCurve, GeomObject>(local_geom_[i]);
//     else if (local_geom_[i]->instanceType() == Class_Line)
//       {
// 	shared_ptr<Line> axis_line =
// 	  dynamic_pointer_cast<Line, GeomObject>(local_geom_[i]);
// 	axis = shared_ptr<SplineCurve>(axis_line->geometryCurve());
//       }
//     else
//       MESSAGE("Unknown curve type!");
    if (axis.get() == NULL) {
	MESSAGE("Curve does not have a spline representation!");
	return srf;
    }
    // Locate the generatrix
    for (i=0; i<local_geom_.size(); i++)
	if (geom_id_[i] == generatrix_id)
	    break;
//     const SplineCurve& generatrix
// 	= dynamic_cast<const SplineCurve&>(*local_geom_[i]);
    par_cv = dynamic_pointer_cast<ParamCurve>(local_geom_[i]);
    shared_ptr<SplineCurve> generatrix(par_cv->geometryCurve());
//     if (local_geom_[i]->instanceType() == Class_SplineCurve)
//       generatrix = dynamic_pointer_cast<SplineCurve, GeomObject>(local_geom_[i]);
//     else if (local_geom_[i]->instanceType() == Class_Line)
//       {
// 	shared_ptr<Line> gen_line =
// 	  dynamic_pointer_cast<Line, GeomObject>(local_geom_[i]);
// 	generatrix = shared_ptr<SplineCurve>(gen_line->geometryCurve());
//       }
//     else
//       MESSAGE("Unknown curve type!");
    if (generatrix.get() == NULL) {
	MESSAGE("Curve does not have a spline representation!");
	return srf;
    }

    // Rotate the generatrix about the axis to the start angle
    Point ab(3);
    axis->point(ab, axis->startparam());    
    Point ae(3);
    axis->point(ae, axis->endparam());
    Point adir = ae-ab;
    MatrixXD<double, 3> rot;
    rot.setToRotation(start_angle, adir[0], adir[1], adir[2]);
    // Go through the control points of the generatrix, transform.
    // Doing it on a copy of the generatrix.
    SplineCurve gen(*generatrix);
    int numctrl = gen.numCoefs();
    Array<double, 3> cp;
    Array<double, 3> axis_begin(ab.begin());
    for (int j = 0; j < numctrl; ++j) {
	cp[0] = gen.coefs_begin()[3*j];
	cp[1] = gen.coefs_begin()[3*j+1];
	cp[2] = gen.coefs_begin()[3*j+2];
	// Translate...
	cp -= axis_begin;
	// ... rotate...
	cp = rot*cp;
	// ...translate back.
	cp += axis_begin;
	gen.coefs_begin()[3*j] = cp[0];
	gen.coefs_begin()[3*j+1] = cp[1];
	gen.coefs_begin()[3*j+2] = cp[2];
    }

    // Call sisl to get the surface of revolution
    SISLCurve* sisl_cv = Curve2SISL(gen, false);
    SISLSurf* sisl_sf;
    int stat;
    s1302(sisl_cv,
	  1e-6,
	  end_angle-start_angle,
	  ab.begin(),
	  adir.begin(),
	  &sisl_sf,
	  &stat);
    freeCurve(sisl_cv);
    if (stat < 0) {
	THROW("Error in creating surface of revolution.");
    }

    srf = shared_ptr<SplineSurface>(SISLSurf2Go(sisl_sf));
    // SISL creates a surface where the generatrix becomes the
    // first parameter direction. IGES expects the opposite.
    srf->swapParameterDirection();

    // We also make sure that parameter domain in the rotational
    // direction is given by the angles (as for instance trim curves
    // will expect this).
    srf->setParameterDomain(srf->startparam_u(), srf->endparam_u(),
			    start_angle, end_angle);
#ifdef IGESLIB_DEBUG
    std::ofstream out("tmp/ruled_sf.g2");
    out << "400 1 0 4 255 255 0 255\n2\n"
	<< ab << "  " << ae << std::endl;
    srf->writeStandardHeader(out);
    srf->write(out);
    gen.writeStandardHeader(out);
    gen.write(out);
#endif // IGESLIB_DEBUG

    return srf;
}

//-----------------------------------------------------------------------------
shared_ptr<SplineSurface>
IGESconverter::readIGEStabulatedCylinder(const char* start, int num_lines)
//-----------------------------------------------------------------------------
{
    MESSAGE("@@sbr Under construction!");

    // The surface is defined as:
    // sf(u,v) = cv(u) + v*(pt-cv(0)),
    // for 1 <= u <= 1 & 1 <= v <= 1,
    // assuming cv is parametrized to the unit interval.
    // The 'pt' is the end point of the generatrix.

    shared_ptr<SplineSurface> srf;

    char pd = header_.pardel;
    char rd = header_.recdel;

    // Verify that start points to '122,.....' warn if not
    int type = readIGESint(start, pd, rd);
    DEBUG_ERROR_IF(type!=122, "Entity is not of type 122.");
    skipDelimiter(start, pd);

    // Read parameters
    int directrix_id = readIGESint(start, pd, rd);  // Directrix curve entity.
    skipDelimiter(start, pd);
    double gen_x = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    double gen_y = readIGESdouble(start, pd, rd);
    skipDelimiter(start, pd);
    double gen_z = readIGESdouble(start, pd, rd);

    skipOptionalTrailingArguments(start, pd, rd);

    // Locate the axis
    size_t i;
    for (i=0; i<local_geom_.size(); i++)
	if (geom_id_[i] == directrix_id)
	    break;
//     const SplineCurve& axis
// 	= dynamic_cast<const SplineCurve&>(*local_geom_[i]);
    shared_ptr<ParamCurve> dir_cv =
	dynamic_pointer_cast<ParamCurve>(local_geom_[i]);
    shared_ptr<SplineCurve> directrix_cv(dir_cv->geometryCurve());
//     if (local_geom_[i]->instanceType() == Class_SplineCurve)
//       axis = dynamic_pointer_cast<SplineCurve, GeomObject>(local_geom_[i]);
//     else if (local_geom_[i]->instanceType() == Class_Line)
//       {
// 	shared_ptr<Line> axis_line =
// 	  dynamic_pointer_cast<Line, GeomObject>(local_geom_[i]);
// 	axis = shared_ptr<SplineCurve>(axis_line->geometryCurve());
//       }
//     else
//       MESSAGE("Unknown curve type!");
    if (directrix_cv.get() == NULL) {
	MESSAGE("Curve does not have a spline representation!");
	return srf;
    }
    else
    {
	MESSAGE("Curve does indeed have a spline representation!");
    }

    // We compute the generatrix line segment (the offset value i.e.).
    Point gen_start_pt = directrix_cv->ParamCurve::point(directrix_cv->startparam());
    Point gen_end_pt(gen_x, gen_y, gen_z);
    Point gen_diff = gen_end_pt - gen_start_pt;
    // We then create the spline curve using this offset point.
    shared_ptr<SplineCurve> offset_cv(CurveCreators::offsetCurve(*directrix_cv, gen_diff));

    // We then use lofting to connect the curves.
    vector<shared_ptr<SplineCurve> > cvs(2);
    cvs[0] = directrix_cv;
    cvs[1] = offset_cv;
    srf = shared_ptr<SplineSurface>(CoonsPatchGen::loftSurface(cvs.begin(), 2));

    // // SISL creates a surface where the generatrix becomes the
    // // first parameter direction. IGES expects the opposite.
    // srf->swapParameterDirection();

    // We also make sure that parameter domain in the rotational
    // direction is given by the angles (as for instance trim curves
    // will expect this).
    srf->setParameterDomain(0.0, 1.0, 0.0, 1.0);

#ifndef NDEBUG
    std::ofstream file_out("tmp/debug.g2");
    srf->writeStandardHeader(file_out);
    srf->write(file_out);
    directrix_cv->writeStandardHeader(file_out);
    directrix_cv->write(file_out);
    offset_cv->writeStandardHeader(file_out);
    offset_cv->write(file_out);
#endif // NDEBUG

    return srf;
}

//-----------------------------------------------------------------------------
void
IGESconverter::readIGEScolour(const char* start, int num_lines,
			      vector<double>& colour, string& cname)
//-----------------------------------------------------------------------------
{
    colour.clear(); // Just to be on the safe side...
    char pd = header_.pardel;
    char rd = header_.recdel;

    // Verify that start points to '314,.....' warn if not
    int type = readIGESint(start, pd, rd);
    DEBUG_ERROR_IF(type!=314, "Entity is not of type 314.");

    // The colour object consists of three doubles (RGB, range 0.0 to 100.0).
    for (int i = 0; i < 3; ++i) {
	skipDelimiter(start, pd);
	double val = readIGESdouble(start, pd, rd);
	colour.push_back(val);
    }
    if (start[0] == pd) { // If name string is present it is read, otherwise cname is not set.
	skipDelimiter(start, pd);
	cname = readIGESstring(start, pd, rd);
    }

    skipOptionalTrailingArguments(start, pd, rd);
}


//-----------------------------------------------------------------------------
ftGroupGeom
IGESconverter::readIGESgroupAssembly(const char* start, int num_lines)
//-----------------------------------------------------------------------------
{
    char pd = header_.pardel;
    char rd = header_.recdel;

    // Verify that start points to '402,.....' warn if not
    int type = readIGESint(start, pd, rd);
    DEBUG_ERROR_IF(type!=402, "Entity is not of type 402.");
    skipDelimiter(start, pd);
    // The form number is 7.

    // Read parameters
    int nmb_entries = readIGESint(start, pd, rd);

    ftGroupGeom group;
    for (int i=0; i<nmb_entries; i++)
    {
      skipDelimiter(start, pd);
      int entry_id = readIGESint(start, pd, rd);  // Entity

          // First look for the entity
      size_t h;
      for (h=0; h<local_geom_.size(); h++)
        if (geom_id_[h] == entry_id)
          break;

      DEBUG_ERROR_IF(h==local_geom_.size(),
	       "Missing object for group assembly: #" << entry_id);
      group.addGeomObj(local_geom_[h]);
      geom_used_[h] = 1;
    }
    skipOptionalTrailingArguments(start, pd, rd);

    return group;
}


//-----------------------------------------------------------------------------
bool IGESconverter::readSingleIGESLine(istream& is, char line_terminated[81],
				       int& line_number, IGESSection& sect)
//-----------------------------------------------------------------------------
{
    // Read any lonely endlines
    char c;
    while(is.get(c)) {
	if (c != '\n') {
	    is.putback(c);
	    break;
	}
    }

    // If we have reached end of file, return false
    if (is.eof()) return false;

    // We set the section indicator character to '\000' so
    // our switch further down is guaranteed to work.
    // Usually, this is unneeded, because the same is
    // done after the switch, but doing it here, too, 
    // means that we no longer depend on having a proper
    // line at the top of the file, and we are no longer
    // required to repeatedly call this function with the
    // same buffer as argument.
    line_terminated[72] = 0;

    // Read a line
    is.get(line_terminated, 81);

    // Get rid of the trailing newline, which would otherwise
    // terminate next get operation
    is.get(c);

    switch (line_terminated[72])
	{
	case 'S':
	    {
		sect = S;
		break;
	    }
	case 'G':
	    {
		sect = G;
		break;
	    }
	case 'D':
	    {
		sect = D;
		break;
	    }
	case 'P':
	    {
		sect = P;
		break;
	    }
	case 'T':
	    {
		sect = T;
		break;
	    }
	case '\000':
	  {
	    sect = E;
	    break;
	  }
	default:
	    THROW("No valid section code for line.");
	}

    // Having read the section, we put a terminator there, so that
    // the line returned will be a string terminating just after the
    // content (including columns 0..71):
    line_terminated[72] = 0;

    // Read line numbers
    line_number = atoi(line_terminated + 73);

    return true;
}

//-----------------------------------------------------------------------------
void IGESconverter::writeSingleIGESLine(ostream& os,
					const char line_terminated[73],
					int line_number, IGESSection sect)
//-----------------------------------------------------------------------------
{
    // Write 72 chars from line_terminated
    for (int i=0; i<72; ++i)
	os << line_terminated[i];

    // Write section code
    switch (sect)
	{
	case S:
	    {
		os << 'S';
		break;
	    }
	case G:
	    {
		os << 'G';
		break;
	    }
	case D:
	    {
		os << 'D';
		break;
	    }
	case P:
	    {
		os << 'P';
		break;
	    }
	case T:
	    {
		os << 'T';
		break;
	    }
	default:
	    THROW("No recognized section code: " << sect);
	}

    // Write line number and endline
    os.width(7);
    os << line_number << std::endl;
}

 //-----------------------------------------------------------------------------
 void IGESconverter::SkipSislComments(std::istream& is)
 //-----------------------------------------------------------------------------
 {
    Utils::eatwhite(is);
 
    char c;
    bool eof_reached = false;
 
    while ((eof_reached = ((is.get(c)).eof())) && c == '$')//GoCOMMENT_START)
    {
      while ((eof_reached = ((is.get(c)).eof())) && c != '\n'); //GoCOMMENT_END );
 
      Utils::eatwhite(is);
    }
 
    if (!eof_reached)
       is.putback(c);
 }


// //===========================================================================
// /// Computes the largest gap in the loop specified by the vector of curves
// double
// IGESconverter::computeLoopGap(const std::vector<shared_ptr<CurveOnSurface> >& curves,
// 			      bool pref_par)
// //===========================================================================
// {

//     // Here, we should check that the given curves indeed are forming a
//     // loop, so every endpoint is within space_epsilon of the start of the
//     // next curve.
//     // Also, we make sure that all curves have the same dimension
//     // and that the curves are of the same type

//     int dim = curves[0]->dimension();
//     int n = curves.size();
//     int i;
//     for (i = 1; i < n; ++i) {
// 	if (curves[i]->dimension() != dim) {
// 	    THROW("Curves do not have the same dimension.");
// 	}
//     }
//     Point startp(dim);
//     Point endp(dim);
//     double maxdist = 0.0;
//     double dist;
//     for (i = 1; i < n; ++i) {
// 	shared_ptr<ParamCurve> prev_cv = (pref_par) ?
// 	    curves[i-1]->parameterCurve() : curves[i-1]->spaceCurve();
// 	shared_ptr<ParamCurve> curr_cv = (pref_par) ?
// 	    curves[i]->parameterCurve() : curves[i]->spaceCurve();
// 	if (!prev_cv || !curr_cv)
// 	    return -1.0;
// 	prev_cv->point(endp, prev_cv->endparam());
// 	curr_cv->point(startp, curr_cv->startparam());
// 	dist = endp.dist(startp);
// 	if (dist > maxdist) maxdist = dist;
//     }
//     shared_ptr<ParamCurve> last_cv = (pref_par) ?
// 	curves[n-1]->parameterCurve() : curves[n-1]->spaceCurve();
//     shared_ptr<ParamCurve> first_cv = (pref_par) ?
// 	curves[0]->parameterCurve() : curves[0]->spaceCurve();
//     last_cv->point(endp, last_cv->endparam());
//     first_cv->point(startp, first_cv->startparam());
//     dist = endp.dist(startp);
//     if (dist > maxdist) maxdist = dist;
//     return maxdist;
// }
