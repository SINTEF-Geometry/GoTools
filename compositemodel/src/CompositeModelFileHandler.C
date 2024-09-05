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

#include "GoTools/compositemodel/CompositeModelFileHandler.h"
#include "GoTools/utils/Logger.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/Plane.h"
#include "GoTools/geometry/SurfaceOfLinearExtrusion.h"
#include "GoTools/geometry/Line.h"
#include "GoTools/geometry/Cone.h"
#include "GoTools/geometry/Circle.h"
#include "GoTools/geometry/Ellipse.h"
#include "GoTools/geometry/Sphere.h"
#include "GoTools/geometry/Torus.h"
#include "GoTools/geometry/SurfaceOfRevolution.h"
#include "GoTools/creators/OffsetSurface.h"

#include <sstream>
#include <fstream>

using std::vector;
using std::cout;
using std::endl;

namespace Go
{


//===========================================================================
CompositeModelFileHandler::~CompositeModelFileHandler()
{

}


//===========================================================================
void CompositeModelFileHandler::writeStart(std::ostream& os)
//===========================================================================
{
    os << "<?xml version='1.0' encoding='UTF-8'?>\n";
    os << "<g22>\n"; // Out parent node.
}

//===========================================================================
void CompositeModelFileHandler::writeEnd(std::ostream& os)
//===========================================================================
{
    os << "\n</g22>" << std::endl;
    clear();   // Make ready for a new model
}


//===========================================================================
void CompositeModelFileHandler::writeHeader(const std::string& file_content_info,
                                     std::ostream& os)
//===========================================================================
{
    os << "\n" << indent_ << "<HEADER>\n";
    os << indent_ << indent_  << "<VERSION>1.0</VERSION>\n";
    os << indent_ << indent_ << "<DESCRIPTION>\n";
    os << indent_ << indent_ << indent_ << file_content_info << "\n";
    os << indent_ << indent_ << "</DESCRIPTION>\n";
    os << indent_ << "</HEADER>\n";

    // @@sbr201601 It makes sense to require all ID's to be unique.
    // But this conflicts with ID's stored in topological objects handled in a
    // write function separate from the geometries.
    LOG_WARN("Consider clearing list of used ID's!");
}


//===========================================================================
void CompositeModelFileHandler::writeGeomObj(const vector<shared_ptr<GeomObject> >& geom_obj,
                                      const vector<int>& obj_id,
                                      std::ostream& os)
//===========================================================================
{
    // We put al the geom_obj into a map.
    for (size_t ki = 0; ki < geom_obj.size(); ++ki)
    {
        int id = (geom_obj.size() == obj_id.size()) ? obj_id[ki] : 
	  (int)geom_objects_.size();
        geom_objects_.insert(std::pair<shared_ptr<GeomObject>, int>(geom_obj[ki], id));

        os << "\n" << indent_ << "<GeomObject ID=\"" << id << "\">\n";
        geom_obj[ki]->writeStandardHeader(os);
        geom_obj[ki]->write(os);
        os << indent_ << "</GeomObject>\n";
    }
}


//===========================================================================
void CompositeModelFileHandler::writeBody(shared_ptr<Body>& body,
					  std::ostream& os, 
					  int body_id)
//===========================================================================
{
    // We start from the top, writing the body (list of shells).
  int nmb_shells = body->nmbOfShells();
  os << "\n" << indent_ << "<Body ID=\"" << body_id << "\">\n";
  if (body->getMaterial() >= 0)
    os << indent_ << indent_ << "<Material>" << body->getMaterial() <<  "</Material>\n";    
  os << indent_ << indent_ << "<Shells>" << nmb_shells;
  for (int ki = 0; ki < nmb_shells; ++ki)
    {
      shared_ptr<SurfaceModel> shell = body->getShell(ki);
      int shell_id = -1;
      auto iter = shells_.find(shell);
      if (iter == shells_.end())
        {
	  shell_id = (int)shells_.size();
	  shells_.insert(std::make_pair(shell, shell_id));
        }
      else
        {
	  shell_id = iter->second;
        }
      os << " " << shell_id;
     }
    os << "</Shells>\n";
    os << indent_ << " </Body>\n";
    
    // Write shells (SurfaceModel)
    for (auto iter = shells_.begin(); iter != shells_.end(); ++iter)
    {
        int shell_id = iter->second;
	writeSurfModel(*(iter->first), os, shell_id, false);
    }

    // Write associated faces
    writeFaces(os);
    
}


//===========================================================================
void CompositeModelFileHandler::writeSurfModels(const std::vector<shared_ptr<Go::SurfaceModel> >& surf_models,
                                                std::ostream& os)
//===========================================================================
{
    for (size_t ki = 0; ki < surf_models.size(); ++ki)
    {
        const bool write_faces = false;
        writeSurfModel(*surf_models[ki], os, ki, write_faces);
    }

    writeFaces(os);
}


//===========================================================================
void CompositeModelFileHandler::writeSurfModel(Go::SurfaceModel& surf_model,
					       std::ostream& os, 
					       int surf_model_id,
					       bool write_faces)
//===========================================================================
{
    // @@sbr201601 I think we should run through all entities and make sure that they have all
    // been assigned an id.

    // We start from the top, writing the surface model (list of faces constituting the shell).
    int nmb_faces = surf_model.nmbEntities();
    os << "\n" << indent_ << "<Shell ID=\"" << surf_model_id << "\">\n";
    os << indent_ << indent_ << "<Closed>" << surf_model.isClosed() <<  "</Closed>\n";
    os << indent_ << indent_ << "<Faces>" << nmb_faces;
    for (int ki = 0; ki < nmb_faces; ++ki)
    {
        shared_ptr<ftFaceBase> face = surf_model.getFace(ki);
        int face_id = -1;
        auto iter = faces_.find(face);
        if (iter == faces_.end())
        {
	  face_id = (int)faces_.size();
            faces_.insert(std::make_pair(face, face_id));
        }
        else
        {
            face_id = iter->second;
        }
        os << " " << face_id;
    }
    os << " </Faces>\n";
    // Open/closed.
    // <# Faces> <Face id 0> <Face id 1> ...
    // Do we need to add orientation flag to the faces? Would assume that this is not necessary
    // for files with 1 SurfaceModel only as a ftSurface will then occur 1 time only.
    // Loop id for faces?
    // We write the tolerances: approxtol, gap, neighbour, kink, bend.
    os << indent_ << indent_ << "<Approxtol>" << surf_model.getApproximationTol() << "</Approxtol>\n";
    tpTolerances top_tol = surf_model.getTolerances();
    os << indent_ << indent_ << "<Gap>" << top_tol.gap << "</Gap>\n";
    os << indent_ << indent_ << "<Neighbour>" << top_tol.neighbour << "</Neighbour>\n";
    os << indent_ << indent_ << "<Kink>" << top_tol.kink << "</Kink>\n";
    os << indent_ << indent_ << "<Bend>" << top_tol.bend << "</Bend>\n";

    os << indent_ << "</Shell>\n";

    // We move on to writing the faces.
    if (write_faces)
      writeFaces(os);
}

//===========================================================================
void CompositeModelFileHandler::writeFaces(std::ostream& os)
//===========================================================================
{

    for (auto iter = faces_.begin(); iter != faces_.end(); ++iter)
    {
        int face_id = iter->second;
        os << "\n" << indent_ << "<Face ID=\"" << face_id << "\">\n";

        shared_ptr<ftSurface> face = dynamic_pointer_cast<ftSurface>(iter->first);
        assert(face.get() != NULL);
        shared_ptr<ParamSurface> surf = face->surface();
        if (surf->instanceType() == Class_BoundedSurface)
        {
            shared_ptr<BoundedSurface> bd_sf = dynamic_pointer_cast<BoundedSurface>(surf);
            surf = bd_sf->underlyingSurface();
        }
        auto surf_iter = geom_objects_.find(surf);
        int surf_id = -1;
	int parent_id = -1;
        if (surf_iter == geom_objects_.end())
        {
            // We add the surface to the geom obj map.
	  surf_id = (int)geom_objects_.size();
            geom_objects_.insert(std::make_pair(surf, surf_id));
        }
        else
        {
            surf_id = surf_iter->second;
        }
	// if (surf->instanceType() == Class_SplineSurface)
	//   {
        //     shared_ptr<SplineSurface> spline_sf = dynamic_pointer_cast<SplineSurface>(surf);
	//     if (spline_sf->isElementarySurface())
	//       {
	if (surf->hasParentSurface())
	  {
	    shared_ptr<ParamSurface> parent = surf->getParentSurface();
	    //shared_ptr<ParamSurface> parent = spline_sf->getElementarySurface();
		auto parent_iter = geom_objects_.find(parent);
		if (parent_iter == geom_objects_.end())
		  {
		    // We add the curve to the geom obj map.
		    parent_id = (int)geom_objects_.size();
		    geom_objects_.insert(std::make_pair(parent, parent_id));
		  }
		else
		  {
		    parent_id = parent_iter->second;
		  }
	      }
	//	  }

        // Pointer to underlying surface.
        os << indent_ << indent_ << "<Surface>" << surf_id << "</Surface>\n";

	if (parent_id >= 0)
	  {
	    // Pointer to original surface
	    os << indent_ << indent_ << "<Parentsurface>" << parent_id << "</Parentsurface>\n";
	  }

        // Pointers to loops.
        int nmb_loops = face->nmbBoundaryLoops();
        os << indent_ << indent_ << "<Loops>" << nmb_loops;
        // Since the geometric surface does not contain an index we must locate the index in the geometry map.
        for (int kj = 0; kj < nmb_loops; ++kj)
        {
            shared_ptr<Loop> bd_loop = face->getBoundaryLoop(kj);
            auto iter2 = loops_.find(bd_loop);
            int loop_id = -1;
            if (iter2 == loops_.end())
            {
	      loop_id = (int)loops_.size();
              loops_.insert(std::make_pair(bd_loop, loop_id));
            }
            else
            {
                loop_id = iter2->second;
            }
            os << " " << loop_id;
        }
        os << "</Loops>\n";

	if (face->hasTwin())
	  {
	    ftSurface *twin = face->twin();
	    int twin_id = -1;
            auto iter_face = faces_.begin();
            while (iter_face != faces_.end()) {
	      if (iter_face->first.get() == twin) {
		break;
	      }
	      ++iter_face;
            }
            if (iter_face == faces_.end())
	      {
                LOG_WARN("Failed finding twin face!");
	      }
            else
	      {
                twin_id = iter_face->second;                
	      }
	    os << indent_ << indent_ << "<Twin>" << twin_id <<  "</Twin>\n"; 
	  }

	// Check for boundary conditions
	if (face->hasBoundaryConditions())
	  {
	    int bd_cond_type, bd_cond;
	    face->getBoundaryConditions(bd_cond_type, bd_cond);

	    os << indent_ << indent_ << "<Boundarycondition>" 
	       << bd_cond_type << " " << bd_cond << "</Boundarycondition>\n"; 
	  }
        os << indent_ << "</Face>\n";
    }

    // We then write the boundary loops.
    for (auto iter = loops_.begin(); iter != loops_.end(); ++iter)
    {
        int loop_id = iter->second;
        os << "\n" << indent_ << "<Loop ID=\"" << loop_id << "\">\n";
        // We must write the number of edges as well as their indices.
        int num_edges = (int)(iter->first->size());
        os << indent_ << indent_ << "<Edges>" << num_edges;
        for (int kj = 0; kj < num_edges; ++kj)
        {
            shared_ptr<ftEdgeBase> edge = iter->first->getEdge(kj);

            auto iter_edge = edges_.find(edge);
            int edge_id;
            if (iter_edge == edges_.end())
            {
	      edge_id = (int)edges_.size();
                edges_.insert(std::make_pair(edge, edge_id));
            }
            else
            {
                edge_id = iter_edge->second;
            }

            os << " " << edge_id;
        }
        os << "</Edges>\n";
	os << indent_ << indent_ << "<Spacetol>" << iter->first->getTol();
	os << "</Spacetol>\n";
        os << indent_ << "</Loop>\n";
    }

    // We write the edges.
#ifndef NDEBUG
    int num_missing_twin = 0;
    std::ofstream fileout_no_twin("tmp/edges_no_twin.g2");
#endif
    for (auto iter = edges_.begin(); iter != edges_.end(); ++iter)
    {
        int edge_id = iter->second;
        os << "\n" << indent_ << "<Edge ID=\"" << edge_id << "\">\n";

        int twin_id = -1;
        ftEdgeBase* twin = iter->first->twin();
        if (twin != NULL) {
            auto iter_edge = edges_.begin();
            while (iter_edge != edges_.end()) {
                if (iter_edge->first.get() == twin) {
                    break;
                }
                ++iter_edge;
            }
            if (iter_edge == edges_.end())
            {
                LOG_WARN("Failed finding twin edge!");
            }
            else
            {
                twin_id = iter_edge->second;                
            }
//            std::cout << "DEBUG INFO: edge_id: " << edge_id << ", twin_id: " << twin_id << std::endl;
        } else {
#ifndef NDEBUG
            ++num_missing_twin;
            // We fetch the 3d curve and write to file.
            shared_ptr<ftEdgeBase> edge = iter->first;
            ftEdge* ft_edge = edge->geomEdge();
            shared_ptr<ParamCurve> geom_cv = ft_edge->geomCurve();
            if (geom_cv->instanceType() == Class_CurveOnSurface) {
                shared_ptr<CurveOnSurface> cv_on_sf = dynamic_pointer_cast<CurveOnSurface>(geom_cv);
                shared_ptr<ParamCurve> space_cv = cv_on_sf->spaceCurve();
                if (space_cv.get() != NULL) {
                    space_cv->writeStandardHeader(fileout_no_twin);
                    space_cv->write(fileout_no_twin);
                } else {
                    LOG_DEBUG("The ft_edge is missing the space curve!");
                }
            } else {
                LOG_WARN("Not a CurveOnSurface, did not see that one coming!");
            }
#endif
//            std::cout << "Edge id " << edge_id << " has no twin. Is surf model not a closed shell?" << std::endl;
        }

        // We need a reference to the geometry curve as well as the boundary points.
        // @@sbr201602 Do we need a reference to the (possibly NULL) twin?
        shared_ptr<Vertex> vertex_start, vertex_end;
        shared_ptr<ftEdge> edge = dynamic_pointer_cast<ftEdge>(iter->first);

        if (edge->face() == NULL)
        {
            LOG_DEBUG("Edge with id " + std::to_string(edge_id) + " is missing its face pointer!");
        }

        bool orientation_ok = edge->orientationOK();
        if (!orientation_ok) {
            LOG_INFO("Orientation is not OK!");
        }
        int par_curve_id = -1;
        int space_curve_id = -1;
        int surf_id = -1;
	int parentcv_id = -1;
        int vertex_start_id = -1;
        int vertex_end_id = -1;
        bool prefer_param = false;
	bool sfcv = true;
	int ccm = 0;
	int constdir = 0;
	double constpar = 0.0;
	int at_bd = -1;
	bool same_orientation = true;
        if (edge.get() != NULL)
        {
            vertex_start = edge->getVertex(true);
            auto iter_start = vertices_.find(vertex_start);
            if (iter_start == vertices_.end())
            {
	      vertex_start_id = (int)vertices_.size();
                vertices_.insert(std::make_pair(vertex_start, vertex_start_id));
            }
            else
            {
                vertex_start_id = iter_start->second;
            }

            vertex_end = edge->getVertex(false);
            auto iter_end = vertices_.find(vertex_end);
            if (iter_end == vertices_.end())
            {
	      vertex_end_id = (int)vertices_.size();
                vertices_.insert(std::make_pair(vertex_end, vertex_end_id));
            }
            else
            {
                vertex_end_id = iter_end->second;
            }

            if (vertex_start_id == vertex_end_id)
            {
                LOG_INFO("Id of start and end vertex are identical, check that this is correct!");
            }
            shared_ptr<ParamCurve> geom_cv = edge->geomCurve();
            shared_ptr<ParamCurve> par_cv, space_cv;
            shared_ptr<ParamSurface> sf;
            if (geom_cv->instanceType() == Class_CurveOnSurface)
            {
//                std::cout << "Encountered a CurveOnSurface!" << std::endl;
                shared_ptr<CurveOnSurface> cv_on_sf = dynamic_pointer_cast<CurveOnSurface>(geom_cv);
                par_cv = cv_on_sf->parameterCurve();
                space_cv = cv_on_sf->spaceCurve();
                prefer_param = cv_on_sf->parPref();
		ccm = cv_on_sf->getCurveTypeInfo();
		cv_on_sf->getConstantCurveInfo(constdir, constpar, at_bd, same_orientation);
                sf = cv_on_sf->underlyingSurface();
            }
            else
            {
                space_cv = geom_cv;
		sfcv = false;
            }

            if (par_cv.get() != NULL)
            {
                auto iter2 = geom_objects_.find(par_cv);
                if (iter2 == geom_objects_.end())
                {
                    // We add the curve to the geom obj map.
		  par_curve_id = (int)geom_objects_.size();
                    geom_objects_.insert(std::make_pair(par_cv, par_curve_id));
                }
                else
                {
                    par_curve_id = iter2->second;
                }
            }

            if (sf.get() != nullptr)
            {
                if (sf->instanceType() == Class_BoundedSurface)
                {
                    shared_ptr<BoundedSurface> bd_sf = dynamic_pointer_cast<BoundedSurface>(sf);
                    sf = bd_sf->underlyingSurface();
                }
                auto surf_iter = geom_objects_.find(sf);
                surf_id = surf_iter->second;
            }

            if (space_cv.get() != NULL)
            {
                auto iter2 = geom_objects_.find(space_cv);
                if (iter2 == geom_objects_.end())
                {
                    // We add the curve to the geom obj map.
		  space_curve_id = (int)geom_objects_.size();
                    geom_objects_.insert(std::make_pair(space_cv, space_curve_id));
                }
                else
                {
                    space_curve_id = iter2->second;
                }
		if (space_cv->instanceType() == Class_SplineCurve)
		  {
		    shared_ptr<SplineCurve> spline_cv = 
		      dynamic_pointer_cast<SplineCurve>(space_cv);
		    if (spline_cv->isElementaryCurve())
		      {
			shared_ptr<ParamCurve> parentcv = spline_cv->getElementaryCurve();
			auto parentcv_iter = geom_objects_.find(parentcv);
			if (parentcv_iter == geom_objects_.end())
			  {
			    // We add the curve to the geom obj map.
			    parentcv_id = (int)geom_objects_.size();
			    geom_objects_.insert(std::make_pair(parentcv, parentcv_id));
			  }
			else
			  {
			    parentcv_id = parentcv_iter->second;
			  }
		      }
		  }
	    }	
	}

        // For CurveOnSurface objects we also have a reference to the surface. If a parameter curve is
        // given this surface is needed for the object to be well defined.
        os << indent_ << indent_ << "<Surface>" << surf_id << "</Surface>\n";

        // We need to handle both parameter curve and geometry curve.
        os << indent_ << indent_ << "<Parametercurve>" << par_curve_id << "</Parametercurve>\n";
        os << indent_ << indent_ << "<Spacecurve>" << space_curve_id << "</Spacecurve>\n";

        int pref_par = (prefer_param) ? 1 : 0;
        os << indent_ << indent_ << "<Preferparam>" << pref_par << "</Preferparam>\n";
	if (sfcv)
	  {
	    // Surface curve information related to constant parameter direction etc
	    os << indent_ << indent_ << "<Surfacecurveinfo>" << ccm << " " << constdir;
	    os << " " << constpar << " " << at_bd << " " << (int)same_orientation << "</Surfacecurveinfo>\n";
	  }
	
	if (parentcv_id >= 0)
	  os << indent_ << indent_ << "<Parentcurve>" << parentcv_id << "</Parentcurve>\n";

        // If reversed then both param_cv + space_cv + vertices are reversed.
        int reversed = (edge->isReversed()) ? 1 : 0;
        os << indent_ << indent_ << "<Reversed>" << reversed << "</Reversed>\n";

        os << indent_ << indent_ << "<Startvertex>" << vertex_start_id << "</Startvertex>\n";
        os << indent_ << indent_ << "<Endvertex>" << vertex_end_id << "</Endvertex>\n";
        os << indent_ << indent_ << "<Twin>" << twin_id << "</Twin>\n";
        // We need reference to the geometry curve as well as bd vertices.
        os << indent_ << "</Edge>\n";
    }
#ifndef NDEBUG
    LOG_DEBUG("Write: Number of edges without a twin: " + std::to_string(num_missing_twin));

#if 0
    std::string log_message_str("writeSurfaceModel: # edges: " + std::to_string(edges_.size()) +
                                ", # edges without a twin: " + std::to_string(num_without_twin) +
                                ", # deg circles without a twin: " + std::to_string(num_without_twin_deg_circle));
    if (num_without_twin > 0)
    {
        BOOST_LOG_TRIVIAL(warning) << log_message_str;
    }
    else
    {
        BOOST_LOG_TRIVIAL(debug) << log_message_str;
    }
#endif

#endif

    // We write the vertices.
    std::streamsize prev = os.precision(15);
    for (auto iter = vertices_.begin(); iter != vertices_.end(); ++iter)
    {
        Point vertex_pt = iter->first->getVertexPoint();
        const int vertex_id = iter->second;
        assert(vertex_pt.dimension() == 3);
        os << "\n" << indent_ << "<Node ID=\"" << vertex_id << "\">\n";
        os << indent_ << indent_ << "<Loc>" << vertex_pt[0] << " " << vertex_pt[1] << " " << vertex_pt[2] << "</Loc>\n";
        os << indent_ << "</Node>\n";
    }
    os.precision(prev);   // Reset precision to it's previous value

    // And finally we write the geometric objects.
    for (auto iter = geom_objects_.begin(); iter != geom_objects_.end(); ++iter)
    {
        int geom_id = iter->second;
        os << "\n" << indent_ << "<GeomObject ID=\"" << geom_id << "\">\n";
        iter->first->writeStandardHeader(os);
        iter->first->write(os);
        os << indent_ << "</GeomObject>\n";
    }

}


//===========================================================================
vector<shared_ptr<GeomObject> > CompositeModelFileHandler::readGeomObj(const char* filein)
//===========================================================================
{
    vector<shared_ptr<GeomObject> > geom_objs;
    // We search after 'GeomObject' nodes.
    pugi::xml_document xml_doc; 
    pugi::xml_parse_result result = xml_doc.load_file(filein);
#ifndef NDEBUG
    LOG_DEBUG("Load result fetchGeomObj: " + std::string(result.description()) + ".");
#endif

    Go::ObjectHeader obj_header;
    pugi::xml_node parent = xml_doc.first_child();
    for (pugi::xml_node node = parent.child("GeomObject"); node; node = node.next_sibling("GeomObject"))
    {
        // We fetch the string with the data.
        std::string geom_obj_string = node.child_value();
//        std::cout << geom_obj_string << std::endl;
        std::istringstream geom_obj_stream(geom_obj_string);
        obj_header.read(geom_obj_stream);
//        std::cout << "classType: " << obj_header.classType() << std::endl;
        if (obj_header.classType() == Class_BoundedSurface)
        {
            shared_ptr<BoundedSurface> bd_sf(new BoundedSurface());
            bd_sf->read(geom_obj_stream);
            geom_objs.push_back(bd_sf);
        }
        else
        {
            shared_ptr<GeomObject> geom_obj = createGeomObject(obj_header);
            geom_obj->read(geom_obj_stream);
            geom_objs.push_back(geom_obj);
//            std::cout << "readGeomObj(): Not yet supporting objects of type: " << obj_header.classType() << std::endl;
        }
    }

    return geom_objs;
}


//===========================================================================
vector<shared_ptr<ParamSurface> > 
CompositeModelFileHandler::readSurface(const char* filein)
//===========================================================================
{
  readFaces(filein);
  
  vector<shared_ptr<ParamSurface> > surfs;
  for (std::map<int, shared_ptr<ftSurface> >::iterator it=faces2_.begin();
       it != faces2_.end(); ++it)
    surfs.push_back(it->second->surface());

  return surfs;
}


//===========================================================================
SurfaceModel CompositeModelFileHandler::readSurfModel(const char* filein, int id)
//===========================================================================
{
    pugi::xml_document xml_doc; 
    pugi::xml_parse_result result = xml_doc.load_file(filein);
#ifndef NDEBUG
    LOG_DEBUG("Load result fetchGeomObj: " + std::string(result.description()) + ".");
#endif

    // If not previously done, read all faces and store them
    if (faces2_.size() == 0)
      readFaces(filein);

    // We start by fetching all the tolerances.
    pugi::xml_node parent = xml_doc.first_child();

    double gap_val, approxtol_val, neighbour_val, kink_val, bend_val;
    for (pugi::xml_node node = parent.child("Shell"); node; node = node.next_sibling("Shell"))
    {
      int shell_id = node.attribute("ID").as_int();
      if (id >= 0 && id != shell_id)
	continue;

      shared_ptr<SurfaceModel> surf_model = getSurfModel(node);
      return *surf_model;
    }

    vector<shared_ptr<ftSurface> > vec;
    SurfaceModel dummy(vec, 1.0e-4); // Requested surface model not found
    return dummy;
}


//===========================================================================
vector<shared_ptr<SurfaceModel> > CompositeModelFileHandler::readSurfModels(const char* g22_filein)
//===========================================================================
{
    vector<shared_ptr<SurfaceModel> > surf_models;

    pugi::xml_document xml_doc; 
    pugi::xml_parse_result result = xml_doc.load_file(g22_filein);
#ifndef NDEBUG
    LOG_DEBUG("Load result fetchGeomObj: " + std::string(result.description()) + ".");
#endif

    // If not previously done, read all faces and store them
    if (faces2_.size() == 0)
      readFaces(g22_filein);

    // We start by fetching all the tolerances.
    pugi::xml_node parent = xml_doc.first_child();

    for (pugi::xml_node node = parent.child("Shell"); node; node = node.next_sibling("Shell"))
    {
      shared_ptr<SurfaceModel> surf_model = getSurfModel(node);
      surf_models.push_back(surf_model);
    }

    return surf_models;
}


//===========================================================================
  shared_ptr<SurfaceModel> 
  CompositeModelFileHandler::readShell(const char* filein,
					   int id)
//===========================================================================
{
  shared_ptr<SurfaceModel> surf_model;
  pugi::xml_document xml_doc; 
  pugi::xml_parse_result result = xml_doc.load_file(filein);
#ifndef NDEBUG
  LOG_DEBUG("Load result fetchGeomObj: " + std::string(result.description()) + ".");
#endif

  // If not previously done, read all faces and store them
  if (faces2_.size() == 0)
    readFaces(filein);

  // We start by fetching all the tolerances.
  pugi::xml_node parent = xml_doc.first_child();

  double gap_val, approxtol_val, neighbour_val, kink_val, bend_val;
  for (pugi::xml_node node = parent.child("Shell"); node; node = node.next_sibling("Shell"))
    {
      int shell_id = node.attribute("ID").as_int();
      if (id >= 0 && id != shell_id)
	continue;

      pugi::xml_node gap_node = node.child("Gap");
      const std::string gap_string = gap_node.child_value();
      std::istringstream gap_ss(gap_string);
      gap_ss >> gap_val;

      pugi::xml_node approxtol_node = node.child("Approxtol");
      const std::string approxtol_string = approxtol_node.child_value();
      std::istringstream approxtol_ss(approxtol_string);
      approxtol_ss >> approxtol_val;

      pugi::xml_node neighbour_node = node.child("Neighbour");
      const std::string neighbour_string = neighbour_node.child_value();
      std::istringstream neighbour_ss(neighbour_string);
      neighbour_ss >> neighbour_val;

      pugi::xml_node kink_node = node.child("Kink");
      const std::string kink_string = kink_node.child_value();
      std::istringstream kink_ss(kink_string);
      kink_ss >> kink_val;

      pugi::xml_node bend_node = node.child("Bend");
      const std::string bend_string = bend_node.child_value();
      std::istringstream bend_ss(bend_string);
      bend_ss >> bend_val;

      // Assemble faces belonging to the current shell
      pugi::xml_node face_nodes = node.child("Faces");
      const std::string face_id_string = face_nodes.child_value();
      std::istringstream ss(face_id_string);
      int num_faces;
      ss >> num_faces;
      vector<int> face_id(num_faces);
      
      vector<shared_ptr<ftSurface> > shell_faces(num_faces);
      for (int ki = 0; ki < num_faces; ++ki)
        {
	  ss >> face_id[ki];
	  shell_faces[ki] = faces2_.find(face_id[ki])->second;
	}

      bool adjacency_set = true;
      surf_model =  
	shared_ptr<SurfaceModel>(new SurfaceModel(approxtol_val,
						  gap_val,
						  neighbour_val,
						  kink_val,
						  bend_val,
						  shell_faces,
						  adjacency_set));

      break;
    }

    return surf_model;
}

//===========================================================================
  shared_ptr<Body> 
  CompositeModelFileHandler::readBody(const char* filein, int id)
//===========================================================================
{
  shared_ptr<Body> body;
  pugi::xml_document xml_doc; 
  pugi::xml_parse_result result = xml_doc.load_file(filein);
#ifndef NDEBUG
  LOG_DEBUG("Load result fetchGeomObj: " + std::string(result.description()) + ".");
#endif

  // If not previously done, read all faces and store them
  if (faces2_.size() == 0)
    readFaces(filein);

  pugi::xml_node parent = xml_doc.first_child();
  for (pugi::xml_node node = parent.child("Body"); node; node = node.next_sibling("Shell"))
    {
      int body_id = node.attribute("ID").as_int();
      if (id >= 0 && id != body_id)
	continue;

      // Read body properties 
      pugi::xml_node material_node = node.child("Material");
      int material_val = -1;
      if (material_node)
	{
	  const std::string material_string = material_node.child_value();
	  std::istringstream material_ss(material_string);
	  material_ss >> material_val;
	}

      // Read all shells
      pugi::xml_node shell_nodes = node.child("Shells");
      const std::string shell_id_string = shell_nodes.child_value();
      std::istringstream ss(shell_id_string);
      int num_shells;
      ss >> num_shells;
      
      vector<shared_ptr<SurfaceModel> > shells(num_shells);
      for (int ki = 0; ki < num_shells; ++ki)
        {
	  int shell_id;
	  ss >> shell_id;
	  shells[ki] = readShell(filein, shell_id);
	}

      // Create Body
      body = shared_ptr<Body>(new Body(shells, material_val));

      // Set body pointers in all associated faces
      int nmb1 = body->nmbOfShells();
      for (int ki=0; ki<nmb1; ++ki)
	{
	  shared_ptr<SurfaceModel> curr_shell = body->getShell(ki);
	  int nmb2 = curr_shell->nmbEntities();
	  for (int kj=0; kj<nmb2; ++kj)
	    curr_shell->getFace(kj)->setBody(body.get());
	}
      break;
    }

  return body;
}

//===========================================================================
void CompositeModelFileHandler::readFaces(const char* filein)
//===========================================================================
{
  double space_eps = 1.0e-4; // Default. To be set from file
  pugi::xml_document xml_doc; 
  pugi::xml_parse_result result = xml_doc.load_file(filein);

  // We start by fetching all geom_objects_.
  ObjectHeader obj_header;
  pugi::xml_node parent = xml_doc.first_child();
  for (pugi::xml_node node = parent.child("GeomObject"); node; node = node.next_sibling("GeomObject"))
    {
      int geom_obj_id = node.attribute("ID").as_int();
      //        std::cout << "node ID: " << geom_obj_id << std::endl;
      // We fetch the string with the data.
      std::string geom_obj_string = node.child_value();
      //        std::cout << geom_obj_string << std::endl;
      std::istringstream geom_obj_stream(geom_obj_string);

      obj_header.read(geom_obj_stream);
      //        std::cout << "classType: " << obj_header.classType() << std::endl;
      shared_ptr<GeomObject> geom_obj = createGeomObject(obj_header);
      if (geom_obj.get() == NULL)
        {
	  LOG_DEBUG("readFaces(): Not yet supporting objects of type: " + std::to_string(obj_header.classType()));
	  continue;
        }

      geom_obj->read(geom_obj_stream);

      geom_objects2_.insert(std::make_pair(geom_obj_id, geom_obj));
    }
  
  // We then move on to the topological entities.
  // As opposed to when we write we want the index to be our key in the map.

  std::map<int, shared_ptr<Vertex> > vertices;

  std::map<int, shared_ptr<ftEdgeBase> > edges;
  // We store pairs of param_cv & space_cv. The index is the same as that in edges.
  std::map<int, std::pair<shared_ptr<ParamCurve>, shared_ptr<ParamCurve> > > edge_curves;
  std::map<int, bool> edge_curves_pref_par;
  std::map<int, shared_ptr<sfcvinfo> > edge_curves_info;

  std::map<int, shared_ptr<Loop> > loops;
  std::map<int, vector<std::pair<shared_ptr<ParamCurve>, shared_ptr<ParamCurve> > > > loop_curves;
  std::map<int, vector<bool> > loop_curves_pref_par;
  std::map<int, vector<shared_ptr<sfcvinfo> > > loop_curves_cvinfo;

  int num_nodes = 0;
  for (pugi::xml_node node = parent.child("Node"); node; node = node.next_sibling("Node"))
    {
      //        std::cout << "node.attribute(): " << node.attribute("ID").value() << std::endl;
      const int node_id = node.attribute("ID").as_int();

      pugi::xml_node loc_array = node.child("Loc");
      Point loc_pt(3);
      std::string loc_string = loc_array.child_value();
      std::istringstream loc_stream(loc_string);
      loc_pt.read(loc_stream);

      shared_ptr<Vertex> vertex(new Vertex(loc_pt));
      vertices.insert(std::make_pair(node_id, vertex));
      ++num_nodes;
    }

  int num_edges = 0;
#ifndef NDEBUG
  int num_missing_twin = 0;
#endif
  vector<std::pair<int, int> > twin_ids;
  for (pugi::xml_node node = parent.child("Edge"); node; node = node.next_sibling("Edge"))
    {
      const int edge_id = node.attribute("ID").as_int();
      //        std::cout << "edge_id: " << edge_id << std::endl;

      // For edges with a paramater curve we need the surface to define the CurveOnSurface object.
      pugi::xml_node sf_node = node.child("Surface");
      const std::string sf_id_string = sf_node.child_value();
      std::istringstream ss3(sf_id_string);
      int surf_id;
      ss3 >> surf_id;
      auto iter3 = geom_objects2_.find(surf_id);
      shared_ptr<ParamSurface> sf;
      if (iter3 != geom_objects2_.end()) // Surface exists.
        {
	  sf = dynamic_pointer_cast<ParamSurface>(iter3->second);
	  assert(sf.get() != NULL);
        }

      pugi::xml_node par_cv_node = node.child("Parametercurve");
      const std::string par_cv_id_string = par_cv_node.child_value();
      std::istringstream ss2(par_cv_id_string);
      int par_cv_id;
      ss2 >> par_cv_id;
      auto iter2 = geom_objects2_.find(par_cv_id);
      shared_ptr<ParamCurve> par_cv;
      if (iter2 != geom_objects2_.end()) // Parameter curve exists.
        {
	  par_cv = dynamic_pointer_cast<ParamCurve>(iter2->second);
	  assert(par_cv.get() != NULL);
        }

      pugi::xml_node space_cv_node = node.child("Spacecurve");
      const std::string space_cv_id_string = space_cv_node.child_value();
      std::istringstream ss(space_cv_id_string);
      int space_cv_id;
      ss >> space_cv_id;
      auto iter = geom_objects2_.find(space_cv_id);
      shared_ptr<ParamCurve> space_cv;
      if (iter != geom_objects2_.end()) // Space curve exists.
      {
          space_cv = dynamic_pointer_cast<ParamCurve>(iter->second);
          assert(space_cv.get() != NULL);
      }

      pugi::xml_node prefpar_node = node.child("Preferparam");
      const std::string prefpar_string = prefpar_node.child_value();
      std::istringstream ss5(prefpar_string);
      int prefpar_val;
      ss5 >> prefpar_val;
      const bool pref_par = (prefpar_val == 1);

      pugi::xml_node parentcv = node.child("Parentcurve");
      if (parentcv && space_cv->instanceType() == Class_SplineCurve)
	{
	  const std::string parentcv_id_string = parentcv.child_value();
	  std::istringstream ss(parentcv_id_string);
	  int parentcv_id = -1;
	  ss >> parentcv_id;

	  shared_ptr<ParamCurve> parent_cv =
	    dynamic_pointer_cast<ParamCurve>(geom_objects2_.find(parentcv_id)->second); 
	  shared_ptr<ElementaryCurve> parent_cv2 =
	    dynamic_pointer_cast<ElementaryCurve>(parent_cv);
	  if (parent_cv2.get())
	    {
	      shared_ptr<SplineCurve> spline_cv = dynamic_pointer_cast<SplineCurve>(space_cv);
	      spline_cv->setElementaryCurve(parent_cv2);
	    }
	}

      int ccm_val, constdir_val, bd_val, orientation_val;
      double constpar_val;
      bool cvinfo = false;
      pugi::xml_node cvinfo_node = node.child("Surfacecurveinfo");
      if (cvinfo_node)
	{
	  const std::string cvinfo_string = cvinfo_node.child_value();
	  std::istringstream cvinfo_ss(cvinfo_string);
	  cvinfo_ss >> ccm_val;
	  cvinfo_ss >> constdir_val;
	  cvinfo_ss >> constpar_val;
	  cvinfo_ss >> bd_val;
	  cvinfo_ss >> orientation_val;
	  cvinfo = true;
	}

      pugi::xml_node reversed_node = node.child("Reversed");
      const std::string reversed_string = reversed_node.child_value();
      std::istringstream ss6(reversed_string);
      int reversed_val;
      ss6 >> reversed_val;
      const bool reversed = (reversed_val == 1);

      // Expecting the edge to be bounded by two vertices.
      pugi::xml_node v1_node = node.child("Startvertex");
      const std::string v1_id_string = v1_node.child_value();
      std::istringstream ss4(v1_id_string);
      int v1_id;
      ss4 >> v1_id;
      shared_ptr<Vertex> v1 = vertices.find(v1_id)->second;
      assert(v1.get() != NULL);

      pugi::xml_node v2_node = node.child("Endvertex");
      const std::string v2_id_string = v2_node.child_value();
      std::istringstream ss7(v2_id_string);
      int v2_id;
      ss7 >> v2_id;
      shared_ptr<Vertex> v2 = vertices.find(v2_id)->second;
      assert(v2.get() != NULL);

      pugi::xml_node twin_node = node.child("Twin");
      const std::string twin_id_string = twin_node.child_value();
      std::istringstream twin_ss(twin_id_string);
      int twin_id;
      twin_ss >> twin_id;
      //        std::cout << "Read: edge_id: " << edge_id << ", twin_id: " << twin_id << std::endl;
#ifndef NDEBUG
      if (twin_id == -1)
        {
	  ++num_missing_twin;
        }
#endif

      // If the file does not contain surface info for the edge we need to search for it.
      shared_ptr<CurveOnSurface> cv_on_sf;
      if ((par_cv.get() != nullptr))
      {
          if (sf.get() == nullptr)
          {
              sf = findSurface(edge_id, parent);
          }

          // If the surface is missing we at least need the space curve to define the geometry.
          if (sf.get() == nullptr)
          {
              assert(space_cv.get() != nullptr);
          }
          else
          {
              cv_on_sf = shared_ptr<CurveOnSurface>(new CurveOnSurface(sf,
                                                                       par_cv,
                                                                       space_cv,
                                                                       pref_par));
          }
      }

      // @@sbr202001 We need the sf geometry when creating the ftEdge (required by setVertices()).
      // The sf geometry is part of a Face, together with a Loop that contains this edge.

      shared_ptr<ParamCurve> cv = (cv_on_sf.get() != nullptr) ? cv_on_sf : space_cv;
      if (cv.get() == nullptr)
      {
          LOG_WARN("Missing space curve, par curve misses the surface!");
      }
      assert(cv.get() != nullptr);

      shared_ptr<ftEdge> edge(new ftEdge(cv, v1, v2, reversed));
      bool orientation_ok = edge->orientationOK();
      if (!orientation_ok) {
	LOG_INFO("Orientation is not OK!");
      }

      twin_ids.push_back(std::make_pair(edge_id, twin_id));

      edges2_.insert(std::make_pair(edge_id, edge));
      edge_curves.insert(std::make_pair(edge_id, std::make_pair(par_cv, space_cv)));
      edge_curves_pref_par.insert(std::make_pair(edge_id, pref_par));
      edge_curves_info.insert(std::make_pair(edge_id, cvinfo ? 
					     shared_ptr<sfcvinfo>(new sfcvinfo(ccm_val, 
									       constdir_val, 
									       constpar_val,
									       bd_val, 
									       orientation_val)) :
					     shared_ptr<sfcvinfo>(new sfcvinfo())));
      ++num_edges;
    }

#ifndef NDEBUG
  LOG_DEBUG("Read: Number of edges without a twin: " + std::to_string(num_missing_twin));
#endif

  // We connect the twins.
  for (auto iter = twin_ids.begin(); iter != twin_ids.end(); ++iter)
    {
      int id1 = iter->first;
      int id2 = iter->second;
      if ((id1 != -1) && (id2 != -1))
        {
	  // Connecting the twin edges.
	  auto iter1 = edges2_.find(id1);
	  auto iter2 = edges2_.find(id2);
	  int status = -1;
	  if (iter1->second->twin() == NULL)
            {
	      iter1->second->connectTwin(iter2->second.get(), status);
            }
	  // I do not think we need to do the opposite as both connections are made in the first connect call.
	  //            iter2->second->connectTwin(iter1->second.get(), status);
        }
    }

  int num_loops = 0;
  for (pugi::xml_node node = parent.child("Loop"); node; node = node.next_sibling("Loop"))
    {
      const int node_id = node.attribute("ID").as_int();

      pugi::xml_node edge_nodes = node.child("Edges");
      const std::string edges_id_string = edge_nodes.child_value();
      std::istringstream ss(edges_id_string);
      int num_edges;
      ss >> num_edges;
      vector<int> edge_id(num_edges);
      vector<shared_ptr<ftEdgeBase> > loop_edges(num_edges);
      vector<std::pair<shared_ptr<ParamCurve>, shared_ptr<ParamCurve> > > curves(num_edges);
      vector<bool> par_pref(num_edges);
      vector<shared_ptr<sfcvinfo> > curves_info(num_edges);
      for (int ki = 0; ki < num_edges; ++ki)
        {
	  ss >> edge_id[ki];
	  loop_edges[ki] = edges2_.find(edge_id[ki])->second;

	  curves[ki] = edge_curves.find(edge_id[ki])->second;
	  par_pref[ki] = edge_curves_pref_par.find(edge_id[ki])->second;
	  curves_info[ki] = edge_curves_info.find(edge_id[ki])->second;
        }

      pugi::xml_node eps_node = node.child("Spacetol");
      if (eps_node)
	{
	  const std::string eps_string = eps_node.child_value();
	  std::istringstream eps_ss(eps_string);
	  eps_ss >> space_eps;
	}
      shared_ptr<Loop> loop(new Loop(loop_edges, space_eps));
      loops.insert(std::make_pair(node_id, loop));
      loop_curves.insert(std::make_pair(node_id, curves));
      loop_curves_pref_par.insert(std::make_pair(node_id, par_pref));
      loop_curves_cvinfo.insert(std::make_pair(node_id, curves_info));
      ++num_loops;
    }

  int num_faces = 0;
  vector<std::pair<int, int> > facetwin_ids;
  for (pugi::xml_node node = parent.child("Face"); node; node = node.next_sibling("Face"))
    {
      const int node_id = node.attribute("ID").as_int();

      pugi::xml_node surface = node.child("Surface");
      const std::string surface_id_string = surface.child_value();
      std::istringstream ss(surface_id_string);
      int surface_id;
      ss >> surface_id;

      pugi::xml_node parentsurf = node.child("Parentsurface");
      int parentsurf_id = -1;
      if (parentsurf)
	{
	  const std::string parentsurf_id_string = parentsurf.child_value();
	  std::istringstream ss(parentsurf_id_string);
	  ss >> parentsurf_id;
	}
     pugi::xml_node loop_nodes = node.child("Loops");
      const std::string loops_id_string = loop_nodes.child_value();
      std::istringstream ss2(loops_id_string);
      int num_loops;
      ss2 >> num_loops;
      vector<int> loop_id(num_loops);
      for (int ki = 0; ki < num_loops; ++ki)
        {
	  ss2 >> loop_id[ki];
        }

      assert(loop_id.size() > 0);

      shared_ptr<ParamSurface> sf = dynamic_pointer_cast<ParamSurface>(geom_objects2_.find(surface_id)->second);
      assert(sf.get() != NULL);
      if (parentsurf_id >= 0)
	{
	  shared_ptr<ParamSurface> parent_sf =
	    dynamic_pointer_cast<ParamSurface>(geom_objects2_.find(parentsurf_id)->second); 
	  shared_ptr<ElementarySurface> parent_sf2 =
	    dynamic_pointer_cast<ElementarySurface>(parent_sf);
	  // if (sf->instanceType() == Class_SplineSurface && parent_sf2.get())
	  //   {
	  //     shared_ptr<SplineSurface> spline_sf = dynamic_pointer_cast<SplineSurface>(sf);
	  //     spline_sf->setElementarySurface(parent_sf2);
	  //   }
	  // Fetch pointer to underlying spline surface, if any
	  SplineSurface* spline_sf = sf->getSplineSurface();
	  if (spline_sf != NULL && parent_sf2.get())
	    spline_sf->setElementarySurface(parent_sf2);
	}

      vector<vector<shared_ptr<CurveOnSurface> > > all_loop_cvs(loop_id.size());
      vector<shared_ptr<Loop> > all_loops(loop_id.size());
      double max_loop_eps = -1.0;
      for (size_t ki = 0; ki < loop_id.size(); ++ki)
        {
	  all_loops[ki] = loops.find(loop_id[ki])->second;
	  assert(all_loops[ki].get() != NULL);

          double loop_eps = all_loops[ki]->getTol();
          // To ensure that all loops are valid we must pick the largest loop tol.
          max_loop_eps = (ki == 0) ? loop_eps : std::max(max_loop_eps, loop_eps);

	  // We fetch the corresponding edge_curves, both par and space cvs.
	  vector<std::pair<shared_ptr<ParamCurve>, shared_ptr<ParamCurve> > >  curves =
	    loop_curves.find(loop_id[ki])->second;
	  vector<shared_ptr<CurveOnSurface> > loop_cvs(curves.size());
	  vector<bool> prefer_param = loop_curves_pref_par.find(loop_id[ki])->second;
	  vector<shared_ptr<sfcvinfo> > cv_info = loop_curves_cvinfo.find(loop_id[ki])->second;
	  size_t kd = 1;
	  for (size_t kj = 0; kj < curves.size(); kj=kd)
            {
	      if (cv_info[kj]->infoset_)
		loop_cvs[kj] = 
		  shared_ptr<CurveOnSurface>(new CurveOnSurface(sf,
								curves[kj].first,
								curves[kj].second,
								prefer_param[kj],
								cv_info[kj]->ccm_,
								cv_info[kj]->constdir_,
								cv_info[kj]->constpar_,
								cv_info[kj]->bd_,
								(bool)cv_info[kj]->orientation_));
	      else
		loop_cvs[kj] = 
		  shared_ptr<CurveOnSurface>(new CurveOnSurface(sf,
								curves[kj].first,
								curves[kj].second,
								prefer_param[kj]));
	      all_loop_cvs[ki].push_back(loop_cvs[kj]);

	      // Check for identity of curves
	      for (kd=kj+1; kd<curves.size(); ++kd)
		{
		  if (curves[kj].first.get() &&
		      curves[kj].first.get() != curves[kd].first.get())
		    break;
		  if (curves[kj].second.get() &&
		      curves[kj].second.get() != curves[kd].second.get())
		  break;   // Assumes parallelity of different curve identificators.
		}

	      for (size_t ka=kj+1; ka<kd; ++ka)
		loop_cvs[ka] = loop_cvs[kj];
            }

	  // Update curve pointer in edge
	  for (size_t kj=0; kj<loop_cvs.size(); ++kj)
	    all_loops[ki]->getEdge(kj)->geomEdge()->setGeomCurve(loop_cvs[kj]);
        }

      const bool fix_trim_cvs = fix_geom_;
      shared_ptr<BoundedSurface> bd_sf(new BoundedSurface(sf, all_loop_cvs, max_loop_eps, fix_trim_cvs));

      //        shared_ptr<ftSurface> face(new ftSurface(bd_sf, loop, node_id));
      shared_ptr<ftSurface> face(new ftSurface(bd_sf, node_id));
      // This call sets up pointers from all edges to the face.
      face->addBoundaryLoops(all_loops);

      // Check for twins
      int twin_id = -1;
      pugi::xml_node twin_node = node.child("Twin");
      if (twin_node)
	{
	  const std::string twin_string = twin_node.child_value();
	  std::istringstream twin_ss(twin_string);
	  twin_ss >> twin_id;
	}
      
      faces2_.insert(std::make_pair(node_id, face));
      facetwin_ids.push_back(std::make_pair(node_id, twin_id));
      ++num_faces;

      // Check for boundary conditions
      pugi::xml_node bd_node = node.child("Boundarycondition");
      if (bd_node)
	{
	  const std::string bd_string = bd_node.child_value();
	  std::istringstream bd_ss(bd_string);
	  int bd_cond_type, bd_cond;
	  bd_ss >> bd_cond_type;
	  bd_ss >> bd_cond;
	  face->setBoundaryConditions(bd_cond_type, bd_cond);
	}
      
    }

#ifndef NDEBUG
  LOG_DEBUG("\nnum_faces: " + std::to_string(num_faces) + ", num_loops: " + std::to_string(num_loops) +
    ", num_edges: " + std::to_string(num_edges) + ", num_nodes: " + std::to_string(num_nodes));

  // We run through all edges and see if any is missing a face.
  LOG_DEBUG("Checking edge face existence for " + std::to_string(edges2_.size()) + " edges.");
  for (auto iter = edges2_.begin(); iter != edges2_.end(); ++iter)
    {
      if (iter->second->face() == NULL)
        {
	  int edge_id = iter->first;
	  LOG_DEBUG("Edge with id " + std::to_string(edge_id) + " is missing it's face!");
        }
    }
  LOG_DEBUG("Done checking edge face existence.");
#endif

  // We connect the face twins.
  for (auto iter = facetwin_ids.begin(); iter != facetwin_ids.end(); ++iter)
    {
      int id1 = iter->first;
      int id2 = iter->second;
      if ((id1 != -1) && (id2 != -1))
        {
	  // Connecting the twin edges.
	  auto iter1 = faces2_.find(id1);
	  auto iter2 = faces2_.find(id2);
	  
	  if (iter1->second->twin() == NULL)
            {
	      iter1->second->connectTwin(iter2->second.get(), space_eps);
            }
        }
    }
  
}

  //===========================================================================
  shared_ptr<GeomObject> CompositeModelFileHandler::createGeomObject(const ObjectHeader& obj_header) const
  //===========================================================================
  {
    shared_ptr<GeomObject> geom_obj;
    // We do not support CurveOnSurface and BoundedSurface as those entities must be expressed
    // using their components.
    if (obj_header.classType() == Class_SplineCurve)
    {
        geom_obj = shared_ptr<SplineCurve>(new SplineCurve());
    }
    else if (obj_header.classType() == Class_Line)
    {
        geom_obj = shared_ptr<Line>(new Line());
    }
    else if (obj_header.classType() == Class_Circle)
    {
        geom_obj = shared_ptr<Circle>(new Circle());
    }
    else if (obj_header.classType() == Class_Ellipse)
    {
        geom_obj = shared_ptr<Ellipse>(new Ellipse());
    }
    else if (obj_header.classType() == Class_SplineSurface)
    {
        geom_obj = shared_ptr<SplineSurface>(new SplineSurface());
    }
    else if (obj_header.classType() == Class_Cylinder)
    {
        geom_obj = shared_ptr<Cylinder>(new Cylinder());
    }
    else if (obj_header.classType() == Class_Plane)
    {
        geom_obj = shared_ptr<Plane>(new Plane());
    }
    else if (obj_header.classType() == Class_Cone)
    {
        geom_obj = shared_ptr<Cone>(new Cone());
    }
    else if (obj_header.classType() == Class_Sphere)
    {
        geom_obj = shared_ptr<Sphere>(new Sphere());
    }
    else if (obj_header.classType() == Class_Torus)
    {
        geom_obj = shared_ptr<Torus>(new Torus());
    }
    else if (obj_header.classType() == Class_SurfaceOfLinearExtrusion)
    {
        geom_obj = shared_ptr<SurfaceOfLinearExtrusion>(new SurfaceOfLinearExtrusion());
    }
    else if (obj_header.classType() == Class_SurfaceOfRevolution)
    {
        geom_obj = shared_ptr<SurfaceOfRevolution>(new SurfaceOfRevolution());
    }
    else if (obj_header.classType() == Class_OffsetSurface)
    {
        geom_obj = shared_ptr<OffsetSurface>(new OffsetSurface());
    }
    else
    {
        LOG_WARN("createGeomObject(): Not yet supporting objects of type: " + std::to_string(obj_header.classType()));
    }

    return geom_obj;
}


//===========================================================================
void CompositeModelFileHandler::clear()
//===========================================================================
{
    geom_objects_.clear();

    shells_.clear();
    faces_.clear();
    loops_.clear();
    edges_.clear();
    vertices_.clear();

    edges2_.clear();
    faces2_.clear();
}


//===========================================================================
shared_ptr<SurfaceModel> CompositeModelFileHandler::getSurfModel(const pugi::xml_node& shell_node)
//===========================================================================
{
    double gap_val, approxtol_val, neighbour_val, kink_val, bend_val;

    pugi::xml_node gap_node = shell_node.child("Gap");
    const std::string gap_string = gap_node.child_value();
    std::istringstream gap_ss(gap_string);
    gap_ss >> gap_val;

    pugi::xml_node approxtol_node = shell_node.child("Approxtol");
    const std::string approxtol_string = approxtol_node.child_value();
    std::istringstream approxtol_ss(approxtol_string);
    approxtol_ss >> approxtol_val;

    pugi::xml_node neighbour_node = shell_node.child("Neighbour");
    const std::string neighbour_string = neighbour_node.child_value();
    std::istringstream neighbour_ss(neighbour_string);
    neighbour_ss >> neighbour_val;

    pugi::xml_node kink_node = shell_node.child("Kink");
    const std::string kink_string = kink_node.child_value();
    std::istringstream kink_ss(kink_string);
    kink_ss >> kink_val;

    pugi::xml_node bend_node = shell_node.child("Bend");
    const std::string bend_string = bend_node.child_value();
    std::istringstream bend_ss(bend_string);
    bend_ss >> bend_val;

    // Assemble faces belonging to the current shell
    pugi::xml_node face_nodes = shell_node.child("Faces");
    const std::string face_id_string = face_nodes.child_value();
    std::istringstream ss(face_id_string);
    int num_faces;
    ss >> num_faces;
      
    vector<shared_ptr<ftSurface> > shell_faces(num_faces);
    for (int ki = 0; ki < num_faces; ++ki)
    {
        int face_id;
        ss >> face_id;
        shell_faces[ki] = faces2_.find(face_id)->second;
    }

    bool adjacency_set = true;
    shared_ptr<SurfaceModel> surf_model(new SurfaceModel(approxtol_val,
                                                         gap_val,
                                                         neighbour_val,
                                                         kink_val,
                                                         bend_val,
                                                         shell_faces,
                                                         adjacency_set));

    return surf_model;
}


//===========================================================================
shared_ptr<ParamSurface> CompositeModelFileHandler::findSurface(int edge_id, const pugi::xml_node& parent)
//===========================================================================
{
    // We expect the edge_id to be used by 1 loop only (i.e. 1 surface only).

    // The sf geometry is part of a Face, together with a Loop that contains this edge.
    // We first search for the edge_id among the loops.
    vector<int> loop_id; // There should be 1 loop exactly. Other values means the input (or logic) is wrong.
    for (pugi::xml_node node = parent.child("Loop"); node; node = node.next_sibling("Loop"))
    {
        const int node_id = node.attribute("ID").as_int();

        pugi::xml_node edge_nodes = node.child("Edges");
        const std::string edges_id_string = edge_nodes.child_value();
        std::istringstream ss(edges_id_string);
        int num_edges;
        ss >> num_edges;
        vector<int> edge_id_vec(num_edges);
        // vector<shared_ptr<ftEdgeBase> > loop_edges(num_edges);
        // vector<std::pair<shared_ptr<ParamCurve>, shared_ptr<ParamCurve> > > curves(num_edges);
        // vector<bool> par_pref(num_edges);
        // vector<shared_ptr<sfcvinfo> > curves_info(num_edges);
        for (int ki = 0; ki < num_edges; ++ki)
        {
            ss >> edge_id_vec[ki];
            // loop_edges[ki] = edges2_.find(edge_id[ki])->second;

            // curves[ki] = edge_curves.find(edge_id[ki])->second;
            // par_pref[ki] = edge_curves_pref_par.find(edge_id[ki])->second;
            // curves_info[ki] = edge_curves_info.find(edge_id[ki])->second;
            if (edge_id_vec[ki] == edge_id)
            {
                loop_id.push_back(node_id);
            }
        }
    }

    assert(loop_id.size() == 1);

    // We then search for the Face with reference to this loop.
    int num_faces = 0;
    // vector<std::pair<int, int> > facetwin_ids;
    vector<shared_ptr<ParamSurface> > sfs;
    for (pugi::xml_node node = parent.child("Face"); node; node = node.next_sibling("Face"))
    {
        //        const int node_id = node.attribute("ID").as_int();

        pugi::xml_node loop_nodes = node.child("Loops");
        const std::string loops_id_string = loop_nodes.child_value();
        std::istringstream ss2(loops_id_string);
        int num_loops;
        ss2 >> num_loops;
        vector<int> loop_id_vec(num_loops);
        for (int ki = 0; ki < num_loops; ++ki)
        {
            ss2 >> loop_id_vec[ki];
            if (loop_id_vec[ki] == loop_id[0])
            {
                pugi::xml_node surface = node.child("Surface");
                const std::string surface_id_string = surface.child_value();
                std::istringstream ss(surface_id_string);
                int surface_id;
                ss >> surface_id;
                shared_ptr<ParamSurface> sf = dynamic_pointer_cast<ParamSurface>(geom_objects2_.find(surface_id)->second);
                sfs.push_back(sf);
            }
        }
    }

    assert(sfs.size() == 1);

    return sfs[0];
}


} // namespace Go

