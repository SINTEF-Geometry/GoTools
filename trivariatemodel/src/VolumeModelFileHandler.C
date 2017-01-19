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

#include "GoTools/trivariatemodel/VolumeModelFileHandler.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/CylinderVolume.h"
#include "GoTools/trivariate/SphereVolume.h"
#include "GoTools/trivariate/ConeVolume.h"
#include "GoTools/trivariate/Parallelepiped.h"
#include "GoTools/trivariate/CurveOnVolume.h"
#include "GoTools/trivariate/SurfaceOnVolume.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/Plane.h"
#include "GoTools/geometry/SurfaceOfLinearExtrusion.h"
#include "GoTools/geometry/Line.h"
#include "GoTools/geometry/Cone.h"
#include "GoTools/geometry/Circle.h"
#include "GoTools/geometry/Sphere.h"
#include "GoTools/geometry/Torus.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"

#include "pugixml.hpp"

#include <sstream>
#include <fstream>

using std::vector;
using namespace Go;

//===========================================================================
VolumeModelFileHandler::VolumeModelFileHandler()
  :CompositeModelFileHandler()
//===========================================================================
{
  GoTools::init();
}

//===========================================================================
VolumeModelFileHandler::~VolumeModelFileHandler()
//===========================================================================
{

}

//===========================================================================
void VolumeModelFileHandler::writeVolumeModel(VolumeModel& vol_model,
					      std::ostream& os)
//===========================================================================
{
  int volume_model_id = -1;  // Top level entity
  os << "\n" << indent_ << "<VolumeModel ID=\"" << volume_model_id << "\">\n";
  int nmb_vols = vol_model.nmbEntities();
  os << indent_ << indent_ << "<Volumes>" << nmb_vols;
  for (int ki = 0; ki < nmb_vols; ++ki)
    {
      shared_ptr<ftVolume> vol = vol_model.getBody(ki);
      int vol_id = -1;
      auto iter = volumes_.find(vol);
      if (iter == volumes_.end())
        {
	  vol_id = (int)volumes_.size();
	  volumes_.insert(std::make_pair(vol, vol_id));
        }
      else
        {
	  vol_id = iter->second;
        }
      os << " " << vol_id;
    }
  os << "</Volumes>\n";
  os << indent_ << indent_ << "<Approxtol>" << vol_model.getApproximationTol() << "</Approxtol>\n";
  os << indent_ << " </VolumeModel>\n";
  
  // Write volumes (ftVolume)
  for (auto iter = volumes_.begin(); iter != volumes_.end(); ++iter)
    {
      int vol_id = iter->second;
      writeVolume(iter->first, os, vol_id, false);
    }
  
  // Write associated faces
  writeFaces(os);

}

//===========================================================================
void VolumeModelFileHandler::writeVolume(const shared_ptr<ftVolume>& body,
					 std::ostream& os, 
					 int body_id, bool faces)
//===========================================================================
{
  // Remember geometry volume
  int volume_id = -1;
  shared_ptr<ParamVolume> vol = body->getVolume();
  auto iter_vol = geom_objects_.find(vol);
  if (iter_vol == geom_objects_.end())
    {
      volume_id = (int)geom_objects_.size();
      geom_objects_.insert(std::make_pair(vol , volume_id));
    }
  else
    {
      volume_id = iter_vol->second;
    }
 
  // Write body (list of shells).
  int nmb_shells = body->nmbOfShells();
  vector<shared_ptr<SurfaceModel> > sfmod;
  os << "\n" << indent_ << "<Volume ID=\"" << body_id << "\">\n";
  os << indent_ << indent_ << "<Geovolume>" << volume_id << "</Geovolume>\n";
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
      sfmod.push_back(shell);
      os << " " << shell_id;
    }
  os << "</Shells>\n";
  os << indent_ << " </Volume>\n";
    
  // Write shells (SurfaceModel)
  for (int ki=0; ki<(int)sfmod.size(); ++ki)
    //  for (auto iter = shells_.begin(); iter != shells_.end(); ++iter)
    {
      auto iter = shells_.find(sfmod[ki]);
      int shell_id = iter->second;
      writeSurfModel(*(iter->first), os, shell_id, false);
    }

  // Write associated faces
  if (faces)
    writeFaces(os);
    
}

//===========================================================================
shared_ptr<VolumeModel>
VolumeModelFileHandler::readVolumeModel(const char* filein)
//===========================================================================
{
  shared_ptr<VolumeModel> model;
  shared_ptr<ftVolume> body;
  pugi::xml_document xml_doc; 
  pugi::xml_parse_result result = xml_doc.load_file(filein);
#ifndef NDEBUG
  std::cout << "Load result fetchGeomObj: " << result.description() << "." << std::endl;
#endif
    // If not previously done, read all faces and store them
    if (faces2_.size() == 0)
      readFaces(filein);

    // We start by fetching all the tolerances.
    pugi::xml_node parent = xml_doc.first_child();

    double gap_val, approxtol_val, neighbour_val, kink_val, bend_val;
    for (pugi::xml_node node = parent.child("VolumeModel"); node; node = node.next_sibling("VolumeModel"))
    {
      int shell_id = node.attribute("ID").as_int();
      pugi::xml_node approxtol_node = node.child("Approxtol");
      const std::string approxtol_string = approxtol_node.child_value();
      std::istringstream approxtol_ss(approxtol_string);
      approxtol_ss >> approxtol_val;

     // Assemble volumes belonging to the current model
      pugi::xml_node vol_nodes = node.child("Volumes");
      const std::string vol_id_string = vol_nodes.child_value();
      std::istringstream ss(vol_id_string);
      int num_vols;
      ss >> num_vols;

      vector<shared_ptr<ftVolume> > volumes(num_vols);
      for (int ki = 0; ki < num_vols; ++ki)
        {
	  int vol_id;
	  ss >> vol_id;
	  volumes[ki] = readVolume(filein, vol_id);
	}

      // Create volume model
      tpTolerances toptol = volumes[0]->getTolerances();
      model = shared_ptr<VolumeModel>(new VolumeModel(volumes, toptol.gap, 
						      toptol.neighbour, 
						      toptol.kink, toptol.bend,
						      true));
      break;
    }
    return model;
}

//===========================================================================
shared_ptr<ftVolume> 
VolumeModelFileHandler::readVolume(const char* filein, int id)
//===========================================================================
{
  shared_ptr<ftVolume> body;
  pugi::xml_document xml_doc; 
  pugi::xml_parse_result result = xml_doc.load_file(filein);
#ifndef NDEBUG
  std::cout << "Load result fetchGeomObj: " << result.description() << "." << std::endl;
#endif

  // If not previously done, read all faces and store them
  if (faces2_.size() == 0)
    readFaces(filein);

  pugi::xml_node parent = xml_doc.first_child();
  for (pugi::xml_node node = parent.child("Volume"); node; node = node.next_sibling("Volume"))
    {
      int body_id = node.attribute("ID").as_int();
      if (id >= 0 && id != body_id)
	continue;

      // Read body properties 
      // Geometry volume
      pugi::xml_node geo_vol_node = node.child("Geovolume");
      int geo_vol_id = -1;
      const std::string geo_vol_string = geo_vol_node.child_value();
      std::istringstream geo_vol_ss(geo_vol_string);
      geo_vol_ss >> geo_vol_id;
      auto iter2 = geom_objects2_.find(geo_vol_id);
      shared_ptr<ParamVolume> vol;
      if (iter2 != geom_objects2_.end()) // Parameter curve exists.
        {
	  vol = dynamic_pointer_cast<ParamVolume>(iter2->second);
	  assert(vol.get() != NULL);
        }

      // Material
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
      body = shared_ptr<ftVolume>(new ftVolume(vol, shells));
      body->setMaterial(material_val);

      // Set body pointers in all associated faces
      int nmb1 = body->nmbOfShells();
      for (int ki=0; ki<nmb1; ++ki)
	{
	  shared_ptr<SurfaceModel> curr_shell = body->getShell(ki);
	  int nmb2 = curr_shell->nmbEntities();
	  for (int kj=0; kj<nmb2; ++kj)
	    curr_shell->getFace(kj)->setBody(body.get());
	}

      // Set volume pointers in all SurfaceOnVolume and CurveOnVolume entities
      addVolumePointers(body);

      break;
    }

  return body;
}

//===========================================================================
void VolumeModelFileHandler::addVolumePointers(shared_ptr<ftVolume> body)
//===========================================================================
{
  shared_ptr<ParamVolume> vol = body->getVolume();
  int nmb_shells = body->nmbOfShells();
  for (int ki = 0; ki < nmb_shells; ++ki)
    {
      shared_ptr<SurfaceModel> shell = body->getShell(ki);
      int nmb_faces = shell->nmbEntities();
      for (int kj=0; kj<nmb_faces; ++kj)
	{
	  shared_ptr<ParamSurface> surf = shell->getSurface(kj);
	  shared_ptr<SurfaceOnVolume> vol_sf = 
	    dynamic_pointer_cast<SurfaceOnVolume>(surf);
	  if (!vol_sf.get())
	    {
	      shared_ptr<BoundedSurface> bd_sf = 
		dynamic_pointer_cast<BoundedSurface>(surf);
	      if (bd_sf.get())
		vol_sf = dynamic_pointer_cast<SurfaceOnVolume>(bd_sf->underlyingSurface());
	    }
	  if (vol_sf.get())
	    vol_sf->setVolume(vol);

	  shared_ptr<ftSurface> face = shell->getFace(kj);
	  vector<shared_ptr<ftEdge> > edgs = face->getAllEdges();
	  for (size_t kr=0; kr<edgs.size(); ++kr)
	    {
	      shared_ptr<ParamCurve> curve = edgs[kr]->geomCurve();
	      shared_ptr<CurveOnVolume> vol_cv = 
		dynamic_pointer_cast<CurveOnVolume>(curve);
	      if (!vol_cv.get())
		{
		  shared_ptr<CurveOnSurface> sf_cv = 
		    dynamic_pointer_cast<CurveOnSurface>(curve);
		  if (sf_cv.get())
		    vol_cv = dynamic_pointer_cast<CurveOnVolume>(sf_cv->spaceCurve());
		}
	      if (vol_cv.get())
		vol_cv->setUnderlyingVolume(vol);
	      
	    }
	}
    }
}

//===========================================================================
shared_ptr<GeomObject> VolumeModelFileHandler::createGeomObject(const ObjectHeader& obj_header) const
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
  else if (obj_header.classType() == Class_SplineVolume)
    {
      geom_obj = shared_ptr<SplineVolume>(new SplineVolume());
    }
  else if (obj_header.classType() == Class_CylinderVolume)
    {
      geom_obj = shared_ptr<CylinderVolume>(new CylinderVolume());
    }
  else if (obj_header.classType() == Class_ConeVolume)
    {
      geom_obj = shared_ptr<ConeVolume>(new ConeVolume());
    }
  else if (obj_header.classType() == Class_SphereVolume)
    {
      geom_obj = shared_ptr<SphereVolume>(new SphereVolume());
    }
  else if (obj_header.classType() == Class_Parallelepiped)
    {
      geom_obj = shared_ptr<Parallelepiped>(new Parallelepiped());
    }
  else if (obj_header.classType() == Class_CurveOnVolume)
    {
      geom_obj = shared_ptr<CurveOnVolume>(new CurveOnVolume());
    }
  else if (obj_header.classType() == Class_SurfaceOnVolume)
    {
      geom_obj = shared_ptr<SurfaceOnVolume>(new SurfaceOnVolume());
    }
  else
    {
      std::cout << "createGeomObject(): Not yet supporting objects of type: " << obj_header.classType() << std::endl;
    }

  return geom_obj;
}





