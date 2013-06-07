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

#ifndef IGESLIB_H
#define IGESLIB_H

#include "GoTools/utils/CoordinateSystem.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/LineCloud.h"
#include "GoTools/geometry/Line.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/BoundedCurve.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/utils/Point.h"
#include "GoTools/igeslib/ftGroupGeom.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/Plane.h"
#include "GoTools/utils/config.h"
#include "sislP.h"
#include <memory>
#include <map>
#include <vector>
#include <string>
#include <iostream>

#if 0
#include "Vertex.h"
#include "ftEdgeBase.h"
#include "Loop.h"
#include "ftFaceBase.h"
#endif

namespace Go
{

enum FileFormat { go, disp, IGES };
enum IGESSection { S, G, D, P, T, E };



/// Storage of all data contained in the IGES header 
struct IGESheader
{
public:
  /// Constructor. Needed to handle write to file.
    IGESheader(); 

    char pardel;
    char recdel;
    std::string sending_system_prodid;
    std::string original_filename;
    std::string creator_system_prodid;
    std::string preprocessor_version_number;
    int integer_bits;                       // Typically 32
    int single_prec_magn;                   // Typically 38
    int single_prec_sign;                   // Typically 6
    int double_prec_magn;                   // Typically 308
    int double_prec_sign;                   // Typically 15
    std::string receiving_system_prodid;         // Typically Unknown
    double model_space_scale;
    int unit_flag;
    std::string unit_description;
    int max_weight_grads;
    double max_line_width;
    std::string timestamp_file_changed;
    double min_resolution;
    double max_coordinate;
    std::string author;
    std::string organisation;
    int IGES_version;
    int drafting_standard;
    std::string timestamp_model_changed;
    std::string application_protocol_description;

    int num_parameters;
};


/// The entity number of all supported IGES entites
struct EntityList
{
public:

  EntityList()
  {
      // 20110815: Added entity 118. More than 24 entities supported...
    /// 24 entities currently supported (072009). Note that for
    /// entities with multiple forms some forms may be missing.
    /// Circular arc.
    entities_.push_back(100); 
    /// Composite curve.
    entities_.push_back(102); 
    /// Conic arc.
    entities_.push_back(104); 
    /// Linear path (form 12).
    entities_.push_back(106); 
    /// Plane.
    entities_.push_back(108); 
    /// Line.
    entities_.push_back(110); 
    /// Point.
    entities_.push_back(116); 
    /// Ruled surface.
    entities_.push_back(118); 
    /// Surface of revolution.
    entities_.push_back(120); 
    /// Tabulated cylinder.
    entities_.push_back(122); 
    /// Direction.
    entities_.push_back(123); 
    /// Transformation matrix.
    entities_.push_back(124); 
    /// Rational B-spline curve.
    entities_.push_back(126); 
    /// Rational B-spline surface.
    entities_.push_back(128); 
    /// Boundary.
    entities_.push_back(141); 
    /// Curve on a parametric surface.
    entities_.push_back(142); 
    /// Bounded surface.
    entities_.push_back(143); 
    /// Trimmed (parametric) surface.
    entities_.push_back(144); 
    //     entities_.push_back(186); // Manifold solid b-rep object.
    /// Plane surface.
    entities_.push_back(190); 
    // Right circular cylindrical surface.
    entities_.push_back(192); 
    /// Color definition.
    entities_.push_back(314); 
    /// Group without back pointers assoc. (form 7).
    entities_.push_back(402); 
    /// Vertex List.
    entities_.push_back(502); 
    /// Edge List.
    entities_.push_back(504); 
    /// Loop.
    entities_.push_back(508); 
    /// Face.
    entities_.push_back(510); 
//     entities_.push_back(514); // Shell.
  }

  bool validEntity(int ent);

private:
  std::vector<int> entities_;
};


/// Storage of all data contained in an IGES directory entity
struct IGESdirentry
{
public:
    int entity_type_number;
    int param_data_start;
    int structure;          // Typically 0
    int line_font_pattern;  // Typically 1 (solid)
    int level;              // Typically 0
    int view;               // Typically 0
    int trans_matrix;       // Typically 0
    int label_display;      // Typically 0
    std::string status;          // 8 digits, 
    double line_weight;        // Typically 0
    int color;              // Typically 0 (no color specified)
    int line_count;
    int form;               // 0 for free-form nurbs
    std::string entity_label;
    int entity_number;      // Typically 0

};

/// The converter between the IGES fileformat, the file format
/// used in GoTools, a file format used for SISL. The disp file format
/// is outdated.
class GO_API IGESconverter
{
public:
  /// Default constructor
    IGESconverter();
    /// This constructor calls the appropriate read and write members
    IGESconverter(std::istream& is, FileFormat from_type,
		  std::ostream& os, FileFormat to_type);
    /// Destructure
    ~IGESconverter();

    /// Read a g2-file
    void readgo(std::istream& is);
    // Expecting sisl input stream to start with the number of objects
    // (srfs or crvs).
    /// Read a number of sisl surfaces
    void readsislsrfs(std::istream& is);
    /// Read a number of sisl curves
    void readsislcrvs(std::istream& is);
    void readdisp(std::istream& is);
    /// Read an IGES file
    void readIGES(std::istream& is);

    /// Write the content of this converter to a g2-file
    void writego(std::ostream& os);
    void writedisp(std::ostream& os);
    /// Write the content of this converter to an IGES-file
    void writeIGES(std::ostream& os);

    /// Add one more geometry object to this converter
    void addGeom(shared_ptr<Go::GeomObject> sp);

    /// Get all geometry objects stored in this converter
    const std::vector<shared_ptr<Go::GeomObject> >& getGoGeom()
	{ return geom_; }

    /// Get the colour of entity number i
    const std::vector<double>& getColour(int i)
	{ return colour_[i]; }

    /// Get group information given in the IGES file (specified by form number)
    const std::vector<ftGroupGeom>& getGroup()
	{ return group_; }

    /// Fetch the geometry tolerance given in an IGES file
    double minResolution() const
        { return header_.min_resolution; }

    /// Number of geometric entities
    int num_geom() { return (int)geom_.size(); }
    /// Number of groups of entities
    int num_group() { return (int)group_.size(); }
    /// Gives dangerous but necessary direct access to header info
    IGESheader& header() { return header_; }

private:
    // Data members
    bool filled_with_data_;
    std::vector<shared_ptr<Go::GeomObject> > geom_;
    // For geom_[i]: colour_[i] refers to an RGB-triplet (empty if no
    // colour spec).  Values are in the range 0.0 to 100.0 (percents
    // of standard scale).
    std::vector<std::vector<double> > colour_;
    std::vector<ftGroupGeom> group_;
    IGESheader header_;
    std::vector<IGESdirentry> direntries_;
    EntityList supp_ent_; // The supported iges entities.

    // Local storage of pointers to entities
    std::vector<shared_ptr<Go::GeomObject> > local_geom_;
    std::vector<int> local_colour_;
    std::vector<int> geom_id_;
    std::vector<int> geom_used_;
    std::vector<Go::Point> plane_normal_;
    std::map<int, int> pnumber_to_plane_normal_index_;
    std::vector<std::vector<shared_ptr<Go::CurveOnSurface> > >
      local_loop_;
    std::vector<int> loop_id_;
    std::vector< std::vector<int> > local_comp_curve_;
    std::vector<int> comp_curve_id_;
    // We store information about colour entities (314).
    std::vector<std::vector<double> > colour_objects_;
    std::vector<std::string> colour_name_; // A colour object may also
					   // be given a name.
    std::vector<int> colour_id_; // local_colour_ may have (negated)
				 // index referring to colour_id_.
    // ... and transformation matrices.
    std::map< int, Go::CoordinateSystem<3> > coordsystems_;

#if 0
    // Topological vertices.
    std::vector<std::vector<shared_ptr<Vertex> > > vertices_;
    // Topological edges.
    std::vector<std::vector<shared_ptr<ftEdgeBase> > > edges_;
    // Topological loops.
    std::vector<std::vector<shared_ptr<Loop> > > top_loops_;
    // Topological faces.
    std::vector<std::vector<shared_ptr<ftFaceBase> > > faces_;
#endif

    // std::vector<std::vector<shared_ptr<Go::Point> > > local_points;
    std::vector<int> Pnumber_;
    int num_lines_[5];

    typedef const char* ccp;

    ccp start_of_P_section_;

    // Utility members

    bool readSingleIGESLine(std::istream& is, char line_terminated[81],
			    int& line_number, IGESSection& sect);
    void writeSingleIGESLine(std::ostream& os, const char line_terminated[73],
			     int line_number, IGESSection sect);
    /// If whereami is within the P section, it gives the current line
    /// number.
    int whichPLine(ccp whereami);
    /// Skips the next occurence of the given delimiter. If the first
    /// non-whitespace character is not pd, the function prints a
    /// warning.
    void skipDelimiter(ccp& whereami, char pd);
    /// Returns true if next non-whitespace character encountered is pd,
    /// false if it is not. If next character is neither pd nor rd,
    /// the functions prints a warning and return false.
    bool checkDelimiter(ccp& whereami,
			char wanted_delimiter,
			char alternative_delimiter);
    /// Skips the trailing ";" or "0,0;" part of a P-section entry.
    /// The arguments at the end may be significant in some settings,
    /// used for pointers to various properties, but we just skip them. 
    void skipOptionalTrailingArguments(ccp& whereami, char pd, char rd,
				       int max_to_skip = 2);
    std::string readIGESstring(ccp& start, char pd, char rd = ';');
    std::string writeIGESstring(const std::string& instring);
    double readIGESdouble(ccp& start, char pd, char rd);
    std::string writeIGESdouble(double d);
    int readIGESint(ccp& start, char pd, char rd);
    std::string writeIGESint(int i);
    // Given an IGES colour number (int in the range [0, 8]), return
    // the RGB-value.  In case of no colour specification, return an
    // empty vector.
    std::vector<double> transformFromIGEScolour(int i);
    // If colour_[i] is present, outine adds IGES-entities of type 314
    // and sets negated number of entity. Otherwise 0 is set.
    int transformToIGESColour();
    std::vector<std::vector<double> >
      uniqueColours(std::vector<std::vector<double> >& colours);

    void readIGESheader(std::string g);   // Reads from a string, into header_
    void writeIGESheader(std::string& g); // Writes to a string, from header_
    void writeIGESdirectory(std::string& g,
			    const std::vector<IGESdirentry>& dirent);
    void writeIGESparsect(std::string& g, std::vector<IGESdirentry>& dirent);

    void writeIGEScolour(const std::vector<double>& colour, std::string& g,
			 std::vector<IGESdirentry>& dirent, int& Pcurr);

    IGESdirentry readIGESdirentry(const char* start);
    shared_ptr<Go::SplineSurface>
      readIGESsurface(const char* start, int num_lines);
    void writeIGESsurface(Go::SplineSurface* surf, int colour, std::string& g,
                          std::vector<IGESdirentry>& dirent, int& Pcurr,
			  int dependency = 0);
    shared_ptr<Go::SplineCurve>
      readIGEScurve(const char* start, int num_lines, int direntry_index);
//     shared_ptr<Go::SplineCurve>
    shared_ptr<Go::BoundedCurve>
      readIGESline(const char* start, int num_lines, int direntry_index);
    shared_ptr<Go::PointCloud3D> readIGESpointCloud(const char* start,
						int num_lines);
    shared_ptr<Go::PointCloud3D>
      readIGESdirection(const char* start, int num_lines);
    shared_ptr<Go::ElementaryCurve>
      readIGESconicArc(const char* start, int num_lines, int form);
    shared_ptr<Go::Plane>
      readIGESplaneSurface(const char* start, int num_lines, int form);
    shared_ptr<Go::Cylinder>
      readIGESrightCircularCylindricalSurface(const char* start,
					      int num_lines, int form);
    shared_ptr<Go::SplineCurve>
    readIGEScircularsegment(const char* start,
			    int num_lines, int direntry_index);
    shared_ptr<Go::Plane> readIGESplane(const char* start,
					       int num_lines,
					       int form);
    void writeIGEScurve(Go::SplineCurve* curve, int colour, std::string& g,
                        std::vector<IGESdirentry>& dirent, int& Pcurr,
			int dependency = 0);
    // Composite curve may include entities of type 'point' and
    // type 'connect point'. Currently not supported.
    void readIGEScompositeCurve(const char* start, int num_lines,
				int direntry_index, std::vector<int>& crv_vec);

    // Currently handling (entity 106) form 12.
    shared_ptr<Go::SplineCurve>
      readIGESlinearPath(const char* start, int num_lines, int form);

    shared_ptr<Go::BoundedSurface>
        readIGESboundedSurf(const char* start, int num_lines);
    void writeIGESboundedSurf(Go::BoundedSurface *surf, int colour,
			      std::string& g, 
                              std::vector<IGESdirentry>& dirent, int& Pcurr);
    void readIGESboundary(const char* start, int num_lines,
                          std::vector<shared_ptr<Go::CurveOnSurface> >&
			  crv_vec);
    void writeIGESboundary(shared_ptr<Go::CurveLoop>  loop,
                           std::string& g, std::vector<IGESdirentry>& dirent,
                           int Psurf, int& Pcurr,
			   int dependency = 0);
    // Helper function
    shared_ptr<Go::CurveOnSurface>
    makeCurveOnSurface(int pc_ind,
		       int sc_ind,
		       shared_ptr<Go::ParamSurface> surf,
		       bool pref_param, int ccm);
    void readIGEScurveOnSurf(const char* start, int num_lines,
			     std::vector<shared_ptr<Go::CurveOnSurface> >& crv_vec);
    // void writeIGEScurveOnSurf(Go::CurveOnSurface *curve, std::string& g,
    //                          std::vector<IGESdirentry>& dirent, int& Pcurr);
    shared_ptr<Go::BoundedSurface> readIGEStrimmedSurf(const char* start,
						     int num_lines);
//     shared_ptr<SplineSurface> readIGEStrimmedSurf(const char* start,
// 						     int num_lines);

    shared_ptr<Go::SplineSurface>
        readIGESruledSurface(const char* start, int num_lines, int form);

    shared_ptr<Go::SplineSurface>
        readIGESsurfOfRevolution(const char* start, int num_lines);

    shared_ptr<Go::SplineSurface>
	readIGEStabulatedCylinder(const char* start, int num_lines);


    //    void writeIGEStrimmedSurf(Go::BoundedSurface *surf, std::string& g,
    //                        std::vector<IGESdirentry>& dirent, int& Pcurr);
    void readIGEScolour(const char* start, int num_lines,
			std::vector<double>& colour, std::string& cname);
    // void writeIGEScolour(std::string& g, std::vector<IGESdirentry>& dirent,
    //		 int colour_id, int& Pcurr);

    shared_ptr< Go::CoordinateSystem<3> >
    readIGEStransformation(const char* start, int num_lines);

    ftGroupGeom readIGESgroupAssembly(const char* start, int num_lines);

    shared_ptr<Go::SplineSurface> readdispsurface(std::istream& is);
    void writedispsurface(std::ostream& os, Go::SplineSurface* surf);
    void writedispcurve(std::ostream& os, Go::SplineCurve* crv);
    void writedispboundedSurf(std::ostream& os, Go::BoundedSurface *bdsf);

    void SkipSislComments(std::istream& is);


//     double
//     computeLoopGap(const vector<shared_ptr<Go::CurveOnSurface> >& curves,
// 		   bool pref_par);


};

} // namespace Go

#endif // This is what is 'ended': #ifndef IGESLIB_H
