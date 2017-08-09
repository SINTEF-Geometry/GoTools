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

#include "GoTools/igeslib/IGESconverter.h"
#include "GoTools/viewlib/gvData.h"
#include "GoTools/viewlib/gvObserver.h"
#include "GoTools/viewlib/gvGroupPropertySheet.h"
#include <QObject>
#include <QString>

#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/viewlib/gvCurvePaintable.h"
#include "GoTools/viewlib/gvRectangularSurfacePaintable.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/tesselator/NoopTesselator.h"
#include "GoTools/viewlib/gvNoopPaintable.h"

#include "GoTools/viewlib/DataHandler.h"


using namespace Go;
using std::vector;
// using std::shared_ptr;


//===========================================================================
void gvData::readIges(std::istream& is)
//===========================================================================
{
    //    clear();
    IGESconverter converter;
    try {
        converter.readIGES(is);
    } catch (...) {
        MESSAGE("Something went wrong while reading iges file!");
    }
    vector<shared_ptr<GeomObject> > new_objects = converter.getGoGeom();
    int nmb_new_objs = 0;
    for (size_t i = 0; i < new_objects.size(); ++i) {
        GeomObject* obj = new_objects[i].get();
        if (obj == 0)
	{
            MESSAGE("Missing object (NULL pointer)! Continuing anyway.");
        }

	ALWAYS_ERROR_IF(obj->dimension() != 2 && obj->dimension() != 3, 
			"Dimension must be 2 or 3.");

	shared_ptr<gvColor> gv_col;
	vector<double> iges_col = converter.getColour((int)i);
	if (iges_col.size() != 0)
	    gv_col = shared_ptr<gvColor>(new gvColor((float)iges_col[0]/100.0f,
						     (float)iges_col[1]/100.0f,
						     (float)iges_col[2]/100.0f,
						     1.0f));
	++nmb_new_objs;
	objects_.push_back(new_objects[i]);
	object_colors_.push_back(gv_col);
    }

    recreateDataStructure(nmb_new_objs);
}

//===========================================================================
void gvData::readSislSrfs(std::istream& is)
//===========================================================================
{
    //    clear();
    IGESconverter converter;
    converter.readsislsrfs(is);
    vector<shared_ptr<GeomObject> > new_objects = converter.getGoGeom();
    int nmb_new_objs = (int)new_objects.size();
    objects_.insert(objects_.end(), new_objects.begin(), new_objects.end());

    for (size_t i = 0; i < new_objects.size(); ++i) {
        GeomObject* obj = new_objects[i].get();
        ALWAYS_ERROR_IF(obj->dimension() != 2 && obj->dimension() != 3, 
                        "Dimension must be 2 or 3.");

        shared_ptr<gvColor> gv_col;
        vector<double> iges_col = converter.getColour((int)i);
        if (iges_col.size() != 0)
            gv_col = shared_ptr<gvColor>(new gvColor((float)iges_col[0]/100.0f,
                                                     (float)iges_col[1]/100.0f,
                                                     (float)iges_col[2]/100.0f,
                                                     1.0f));
        object_colors_.push_back(gv_col);
    }
    recreateDataStructure(nmb_new_objs);
}

//===========================================================================
void gvData::readSislCrvs(std::istream& is)
//===========================================================================
{
    //    clear();
    IGESconverter converter;
    converter.readsislcrvs(is);
    vector<shared_ptr<GeomObject> > new_objects = converter.getGoGeom();
    int nmb_new_objs = (int)new_objects.size();
    objects_.insert(objects_.end(), new_objects.begin(), new_objects.end());

    for (size_t i = 0; i < new_objects.size(); ++i) {
        GeomObject* obj = new_objects[i].get();
        ALWAYS_ERROR_IF(obj->dimension() != 2 && obj->dimension() != 3, 
                        "Dimension must be 2 or 3.");

        shared_ptr<gvColor> gv_col;
        vector<double> iges_col = converter.getColour((int)i);
        if (iges_col.size() != 0)
            gv_col = shared_ptr<gvColor>(new gvColor((float)iges_col[0]/100.0f,
                                                     (float)iges_col[1]/100.0f,
                                                     (float)iges_col[2]/100.0f,
                                                     1.0f));
        object_colors_.push_back(gv_col);
    }
    recreateDataStructure(nmb_new_objs);
}

//===========================================================================
void gvData::readGo(std::istream& is)
//===========================================================================
{
    //    clear();

    ObjectHeader header;
    int nmb_new_objs = 0;
    while (!is.eof()) {
        try {
            header.read(is);
        } catch (...) {
            MESSAGE("Failed reading header, using the geometry read so far.");
            break;
        }
        //Read(is, header);
        shared_ptr<GeomObject> obj(Factory::createObject(header.classType()));
        try {
	    if (obj->instanceType() == Class_BoundedSurface)
	    {
		shared_ptr<BoundedSurface> bd_sf = dynamic_pointer_cast<BoundedSurface>(obj);
		bool fix_trim_cvs = false;
		bd_sf->read(is, fix_trim_cvs);
	    }
	    else
	    {
		obj->read(is);
	    }
            // //Read(is, *obj);
            // ALWAYS_ERROR_IF(obj->dimension() != 3 && obj->dimension() != 2, 
            //                 "Dimension must be 2 or 3.");
        } catch (...) {
            MESSAGE("Failed reading (Go) object!");
            obj = shared_ptr<GeomObject>();
        }
        objects_.push_back(obj);
        ++nmb_new_objs;
        shared_ptr<gvColor> gv_col(new gvColor(0.3f, 0.3f, 1.0f, 1.0f)); // Blue.
        if (header.auxdataSize() == 4) {
            gv_col = shared_ptr<gvColor>(new gvColor((float)header.auxdata(0)/255.0f,
                                                     (float)header.auxdata(1)/255.0f,
                                                     (float)header.auxdata(2)/255.0f,
                                                     (float)header.auxdata(3)/255.0f));
        }
        object_colors_.push_back(gv_col);
        Utils::eatwhite(is);
        //SkipComments(is);
    }
    try {
        recreateDataStructure(nmb_new_objs);
    } catch (...) {
        MESSAGE("Failed recreating data structure!");
    }
}


//===========================================================================
void gvData::readGo(const std::vector<shared_ptr<Go::GeomObject> >& new_objects,
                    std::vector<shared_ptr<gvColor> >& new_colors)
//===========================================================================
{
    int nmb_new_objs = (int)new_objects.size();
    if (new_objects.size() != new_colors.size()) {
        new_colors.resize(new_objects.size());
    }

    for (size_t i = 0; i < new_objects.size(); ++i) {
        ALWAYS_ERROR_IF(new_objects[i]->dimension() != 3 &&
                        new_objects[i]->dimension() != 2, "Dimension must be 2 or 3.");
    }
    objects_.insert(objects_.end(), new_objects.begin(), new_objects.end());
    object_colors_.insert(object_colors_.end(), new_colors.begin(), new_colors.end());
    recreateDataStructure(nmb_new_objs);
}


//===========================================================================
void gvData::writeSelectedGo(std::ostream& os)
//===========================================================================
{
    os.precision(15);
    for (int i = 0; i < numObjects(); ++i) {
        if (getSelectedStateObject(i)) {
            vector<int> go_col;
            if (object_colors_[i].get() != 0) {
                float *rgba = object_colors_[i]->rgba; // 4 elements
                for (int j = 0; j < 4; ++j)
                    go_col.push_back((int)rgba[j]*255);
            }
            ObjectHeader header(objects_[i]->instanceType(), MAJOR_VERSION, MINOR_VERSION,
                                  go_col);
            header.write(os);
// 	    objects_[i]->writeStandardHeader(os);
            os << '\n';
            objects_[i]->write(os);
            os << '\n';
        }
    }
    os << std::flush;
}

//===========================================================================
void gvData::writeSelectedIges(std::ostream& os)
//===========================================================================
{
    os.precision(15);
    IGESconverter conv;
    for (int i = 0; i < numObjects(); ++i) {
        if (getSelectedStateObject(i)) {
          conv.addGeom(objects_[i]);
        }
    }

    conv.writeIGES(os);
}


// //===========================================================================
// void gvData::groupSelected()
// //===========================================================================
// {
//     vector<int> selected;
//     for (int i = 0; i < objects_.size(); ++i)
// 	if (getSelectedStateObject(i))
// 	    selected.push_back(i);

//     if (selected.size() > 0) {
// 	gvGroupPropertySheet group_sheet(selected);
// 	group_sheet.createSheet(0, 0);

// 	// If user clicks OK, group is created and pushed back.
// 	QObject::connect(&group_sheet, SIGNAL(value_changed(vector<int>&, QString)),
// 			 this, SLOT(add_group(vector<int>&, QString)));
//     }
// }

//===========================================================================
void gvData::extractSelectedObjects(std::vector< shared_ptr<Go::GeomObject> >& sel_obj)
//===========================================================================
{
    sel_obj.clear();
    for (int i = 0; i < numObjects(); ++i) {
        if (getSelectedStateObject(i)) {
           if (objects_[i].get())
              sel_obj.push_back(objects_[i]);
        }
    }
}



//===========================================================================
void gvData::recreateDataStructure(int nmb_new_objs)
//===========================================================================
{
    //    cout << "gvData::recreateDataStructure... " << flush;

    gvColor blue(0.3, 0.3, 1.0, 1.0);
    gvColor light_blue(0.7, 0.7, 1.0, 1.0);
    gvColor red(1.0, 0.3, 0.3, 1.0);
    gvColor light_red(1.0, 0.7, 0.7, 1.0);

//     painter_.removeAllPaintables();
//     paintables_.clear();
//     tesselators_.clear();
//     property_sheets_.clear();
    int nmb_objs = (int)objects_.size();
    // We start by inserting default values into vectors.
    // Some elements may already have well defined tesselators.
    int i;

    // We make sure that vectors have the right size. Size of vectors should correspond.
    int nmb_old_objs = nmb_objs - nmb_new_objs;
    if (int(tesselators_.size()) != nmb_old_objs) {
        int nmb_to_remove = (int)tesselators_.size() - nmb_old_objs;
        tesselators_.resize(nmb_old_objs);
        paintables_.resize(nmb_old_objs);
        property_sheets_.resize(nmb_old_objs);
        for (i = 0; i < nmb_to_remove; ++i) {
            painter_.removeLastPaintable();
        }
    }

    for (i = nmb_old_objs; i < nmb_old_objs + nmb_new_objs; ++i) {

        //	cout << "Iteration " << i << endl;
        if (objects_[i].get()==NULL)
	{
            MESSAGE("Object missing!");
//            continue;
        }
        gvColor col;
        // @@ Hack to get red curves by default, other things are blue...
//  	if (dynamic_cast<ParamCurve*>(objects_[i].get())) {
//  	    col = red;
//  	}
        if (object_colors_[i].get() != 0) {
            col = *object_colors_[i];
        } else {
            col = blue;
// 	    col = gvColor(0.3, 0.3, 1.0, 1.0);
        }
        datahandler_->clear();
        try {
            datahandler_->create(objects_[i], col, i);
	    if (datahandler_->tesselator())
	    {
		datahandler_->tesselator()->tesselate();
	    }
	    tesselators_.push_back(datahandler_->tesselator());
	    paintables_.push_back(datahandler_->paintable());
	    property_sheets_.push_back(datahandler_->propertySheet());
	    painter_.addPaintable(datahandler_->paintable());
	    // We wait until the object has been created before we require a 2D or 3D object, allowing
	    // the create() routine to lift or project for any other dimension.
	    ALWAYS_ERROR_IF((objects_[i].get() != NULL) && (objects_[i]->dimension() != 3) && (objects_[i]->dimension() != 2),
			    "Dimension must be 2 or 3.");
        } catch (...) {
            MESSAGE("Failed creating some object (index: " << i << ")!");
            tesselators_.push_back(shared_ptr<Tesselator>(new NoopTesselator()));
            paintables_.push_back(shared_ptr<gvPaintable>(new gvNoopPaintable(col, i)));
            property_sheets_.push_back(shared_ptr<gvPropertySheet>());
            painter_.addPaintable(shared_ptr<gvPaintable>());
        }
    }
    try {
        computeBox();
        updateObservers();
    } catch (...) {
        MESSAGE("Failed ...");
    }
    //    cout << " finished" << endl;
}


//===========================================================================
void gvData::computeBox()
//===========================================================================
{
   int i;
   const double unbounded_limit = 1.0e06;//08; // Consider the object unbounded if box diagonal is larger than value.
   for (i = 0; i < numObjects(); ++i) {
      if (object(i).get())
      {
          try {
              box_ = object(i)->boundingBox();
              double box_diag = (box_.low()).dist(box_.high());
              // if ((box_diag < 1.0e08) && (box_diag > 1.0e06)) {
              //     // std::cout << "Box diag larger than 1.0e06! Instance type: " << object(i)->instanceType() << std::endl;
              //     if (object(i)->instanceType() == Class_BoundedSurface) {
              //         shared_ptr<BoundedSurface> bd_sf = dynamic_pointer_cast<BoundedSurface>(object(i));
              //         // std::cout << "Instance type of under_sf: " << bd_sf->underlyingSurface()->instanceType() << std::endl;
              //     }
              // }
              bool bounded = (box_diag < unbounded_limit);
              if (box_.valid() && bounded) {
                  break;
              } else {
                  ;//MESSAGE("Box not valid or unbounded!");
              }
          } catch (...) {
              MESSAGE("Failed fetching bounding box.");
              continue;
          }
      } 
   }
   
   for (; i < numObjects(); ++i) {
      if (object(i).get())
      {
          try {
              BoundingBox next_box = object(i)->boundingBox();
              double box_diag = (next_box.low()).dist(next_box.high());
              // if ((box_diag < 1.0e08) && (box_diag > 1.0e06)) {
              //     std::cout << "Box diag larger than 1.0e06 and smaller than 1.0e08! Instance type: " <<
              //         object(i)->instanceType() << std::endl;
              //     if (object(i)->instanceType() == Class_BoundedSurface) {
              //         shared_ptr<BoundedSurface> bd_sf = dynamic_pointer_cast<BoundedSurface>(object(i));
              //         std::cout << "Instance type of under_sf: " << bd_sf->underlyingSurface()->instanceType() << std::endl;
              //     }
              // }
              bool bounded = (box_diag < unbounded_limit);
              if (box_.valid() && bounded) {
                  box_.addUnionWith(next_box);
              } else {
                  ;//MESSAGE("Box not valid or unbounded!");
              }
          } catch (...) {
              MESSAGE("Failed adding union with next box!");
          }
      }
   }

   if (!box_.valid()) {
       MESSAGE("Final box not valid!");
   }//  else {
   //     double box_diag = (box_.low()).dist(box_.high());
   //     std::cout << "box_diag: " << box_diag << std::endl;
   // }
}

//===========================================================================
Go::BoundingBox gvData::boundingBox(const std::vector<int>& objs) const
//===========================================================================
{
   Go::BoundingBox box;
   const double unbounded_limit = 1.0e06;//08; // Consider the object unbounded if box diagonal is larger than value.
   size_t i;
   for (i = 0; i < objs.size(); ++i)  {
       if (objs.size()>0 && numObjects()>objs[i] && object(objs[i]).get()) {
           try {
               box=object(objs[i])->boundingBox();
               double box_diag = (box.low()).dist(box.high());
               if (box_diag < unbounded_limit) {
                   break;
               }
           } catch (...) {
               MESSAGE("Failed fetching bounding box.");
               continue;
           }
       }
   }
   
   for (; i < objs.size(); ++i) {
       if (objs.size()>0 && numObjects()>objs[i] && object(objs[i]).get()) {
           try {
               BoundingBox next_box = object(objs[i])->boundingBox();
               double next_box_diag = (next_box.low()).dist(next_box.high());
               if (next_box_diag < unbounded_limit) {
                   double box_diag = (box.low()).dist(box.high());
                   box.addUnionWith(next_box);
                   box_diag = (box.low()).dist(box.high());
               }
           } catch (...) {
               MESSAGE("Failed fetching bounding box.");
               continue;
           }
       }
   }

   return box;
}


//===========================================================================
void gvData::updateObservers()
//===========================================================================
{
    if (!updates_enabled_) return;
    std::set<gvObserver*>::iterator it;
    for (it = observers_.begin(); it != observers_.end(); ++it) {
	if (*it)
	{
	    (*it)->observedChanged();
	}
    }
}

//===========================================================================
void gvData::setSelectedStateObject(int id, bool state)
//===========================================================================
{
   if (id < int(objects_.size()) && objects_[id].get())
   {
      bool oldstate = paintables_[id]->selected();
      if (state != oldstate) {
         paintables_[id]->setSelected(state);
         updateObservers();
      }
   }
}


//===========================================================================
void gvData::setVisibleStateObject(int id, bool state)
//===========================================================================
{
   if (id < int(objects_.size()) && objects_[id].get())
   {
      bool oldstate = paintables_[id]->visible();
      if (state != oldstate) {
         paintables_[id]->setVisible(state);
         updateObservers();
      }
   }
}

//===========================================================================
bool gvData::getSelectedStateObject(int id)
//===========================================================================
{
   if (paintables_[id].get())
    return paintables_[id]->selected();
   else
      return false;
}


//===========================================================================
bool gvData::getVisibleStateObject(int id)
//===========================================================================
{
   if (paintables_[id].get())
      return paintables_[id]->visible();
   else
      return false;
}

//===========================================================================
void gvData::clear()
//===========================================================================
{
    objects_.clear();
    object_colors_.clear();
    // box_ is unchanged
    object_groups_.clear();

    property_sheets_.clear();
    tesselators_.clear();
    paintables_.clear();

    painter_.removeAllPaintables();
}

//===========================================================================
void gvData::deleteObj(int obj)
//===========================================================================
{
   if (obj == int(objects_.size())-1)
      clearLast();
   else
   {
      objects_[obj].reset();
      object_colors_[obj].reset();
      property_sheets_[obj].reset();
      tesselators_[obj].reset();
      paintables_[obj].reset();
      painter_.removePaintable(obj);
   }
   updateObservers();
}

//===========================================================================
void gvData::clearLast()
//===========================================================================
{
    objects_.erase(objects_.end() - 1);
    object_colors_.erase(object_colors_.end() - 1);
    // box_ is unchanged
    // @@sbr Hmm, maybe we should make sure that object is not part of a group?

    property_sheets_.erase(property_sheets_.end() - 1);
    tesselators_.erase(tesselators_.end() - 1);
    paintables_.erase(paintables_.end() - 1);

    painter_.removeLastPaintable();
    if (objects_.size()>0 && objects_[objects_.size()-1].get()==NULL)
      clearLast();
}
