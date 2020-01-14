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

#ifndef _COMPOSITEMODELFILEHANDLER_H
#define _COMPOSITEMODELFILEHANDLER_H


#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/GeomObject.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/Body.h"
//#include "GoTools/trivariatemodel/VolumeModel.h"

#include <string>
#include <vector>
#include <map>
#include "pugixml.hpp"


namespace Go
{


// Writing to and reading from the g22 file format.
// Supports topology information (as opposed to the g2 format for geometries only).
  class ParamSurface;

class CompositeModelFileHandler
{

public:
    CompositeModelFileHandler()
      : MAJOR_VERSION_(0), MINOR_VERSION_(1), indent_("  "), fix_geom_(false)
        {}

    ~CompositeModelFileHandler();

    void writeStart(std::ostream& os);

    void writeEnd(std::ostream& os);

    void writeHeader(const std::string& file_content_info,
                     std::ostream& os);

    void writeGeomObj(const std::vector<shared_ptr<GeomObject> >& geom_obj,
                      const std::vector<int>& obj_id,
                      std::ostream& os);

    void writeSurfModels(const std::vector<shared_ptr<Go::SurfaceModel> >& surf_models,
                         std::ostream& os);

    void writeSurfModel(Go::SurfaceModel& surf_model,
                        std::ostream& os, int surfmodel_id=-1,
			bool write_faces=true);

    void writeBody(shared_ptr<Body>& body,
			   std::ostream& os, int body_id=-1);

    std::vector<shared_ptr<GeomObject> > readGeomObj(const char* filein);

    std::vector<shared_ptr<ParamSurface> > readSurface(const char* filein);

    SurfaceModel readSurfModel(const char* g22_filein, int id=-1);

    std::vector<shared_ptr<SurfaceModel> > readSurfModels(const char* g22_filein);

    shared_ptr<SurfaceModel> readShell(const char* filein, int id=-1);

    shared_ptr<Body> readBody(const char* filein, int id=-1);

    /// Set flag for geometry repair. Currently only repair of surface
    /// trimming curves are performed. Default is false
    void setGeomFix(bool fix_geom)
    {
      fix_geom_ = fix_geom;
    }
    
protected:

    struct sfcvinfo
    {
      sfcvinfo(int ccm, int constdir, double constpar, int bd, int orientation)
      {
	ccm_ = ccm;
	constdir_ = constdir;
	constpar_ = constpar;
	bd_ = bd;
	orientation_ = orientation;
	infoset_ = true;
      }

      sfcvinfo()
      {
	infoset_ = false;
      }
      
      bool infoset_;
      int ccm_, constdir_, bd_, orientation_;
      double constpar_;
    };

    const int MAJOR_VERSION_;
    const int MINOR_VERSION_;
    const std::string indent_;

    // We store all id's used and issue a warning if we encounter one that is already used.
    // @@sbr201601 This should be handled in a more robust manner. Possibly fetch list of used id's prior to
    // writing to file, allowing the user to replace id's? I guess this matters the most when mixing models
    // in the same file.
//    std::vector<int> id_used_;

    // All GeomObject's are stored in a map to allow us to set a unique id.
    // We start the indexing at 0 and increment by 1 as the map grows.
    std::map<shared_ptr<GeomObject>, int> geom_objects_;

    std::map<shared_ptr<SurfaceModel>, int> shells_;
    std::map<shared_ptr<ftFaceBase>, int> faces_;
    std::map<shared_ptr<Loop>, int> loops_;
    std::map<shared_ptr<ftEdgeBase>, int> edges_;
    std::map<shared_ptr<Vertex>, int> vertices_;

    // Correspondence between entities and id for reading
    std::map<int, shared_ptr<ftSurface> > faces2_;
    std::map<int, shared_ptr<GeomObject> > geom_objects2_;
    std::map<int, shared_ptr<ftEdgeBase> > edges2_;

    // Whether or not fixes of the geometry is to be performed
    bool fix_geom_;
    
    // Creator for our supported GeomObject's.
    virtual shared_ptr<GeomObject> 
      createGeomObject(const ObjectHeader& obj_header) const;

    void writeFaces(std::ostream& os);
    void readFaces(const char* filein);

    // Remove all geometric and topological data.
    void clear();

    shared_ptr<SurfaceModel> getSurfModel(const pugi::xml_node& shell_node);

    shared_ptr<ParamSurface> findSurface(int edge_id, const pugi::xml_node& parent);

};


} // namespace Go


#endif // _COMPOSITEMODELFILEHANDLER_H

