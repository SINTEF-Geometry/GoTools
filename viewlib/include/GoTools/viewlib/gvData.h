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

#ifndef _GVDATA_H
#define _GVDATA_H

// Standard includes
#include <iostream>
#include <vector>
#include <set>
#include <utility>

// Other includes
#include "GoTools/geometry/GeomObject.h"
#include "GoTools/viewlib/DataHandler.h"
#include "GoTools/tesselator/Tesselator.h"
#include "GoTools/viewlib/gvPaintable.h"
#include "GoTools/viewlib/gvPropertySheet.h"
#include "GoTools/viewlib/gvPainter.h"
#include "GoTools/viewlib/gvGroup.h"

#include "GoTools/utils/config.h"

class gvObserver;

/** gvData reads and stores the geometric data.
 */

class gvData
{
public:
    gvData(std::auto_ptr<DataHandler> datahandler)
	: updates_enabled_(true),
	  datahandler_(datahandler)
    {}

    virtual ~gvData()
    {}


    //    void read3ds(std::string filename);
    void readIges(std::istream& is);
    void readSislSrfs(std::istream& is);
    void readSislCrvs(std::istream& is);
    void readGo(std::istream& is);
    // If size of vectors do not match, new_colors is not used.
    void readGo(const std::vector<shared_ptr<Go::GeomObject> >& new_objects,
		std::vector<shared_ptr<gvColor> >& new_colors);
//      void readSislformat(std::istream& is);

    void writeSelectedGo(std::ostream& os);

    void writeSelectedIges(std::ostream& os);

    void setTexture(int index, shared_ptr<gvTexture> tex)
    { painter_.setTexture(index, tex); }

    shared_ptr<gvTexture> getTexture(int index)
    { return painter_.getTexture(index); }

    int numObjects() const
    { return (int)objects_.size(); }
    shared_ptr<Go::GeomObject> object(int i)
    { 
	if (i>=0 && i < int(objects_.size())) {
	    return objects_[i]; 
	}
	return shared_ptr<Go::GeomObject>();
    }

    shared_ptr<const Go::GeomObject> object(int i) const
    { return objects_[i]; }
    shared_ptr<gvColor> color(int i)
    { return object_colors_[i]; }
    shared_ptr<gvPropertySheet> propertySheet(int i)
    { return property_sheets_[i]; }
    shared_ptr<Go::Tesselator> tesselator(int index)
    { return tesselators_[index]; }
    shared_ptr<gvPaintable> paintable(int index)
    { return paintables_[index]; }

/*     void groupSelected(); */
    void addGroup(gvGroup group)
    { object_groups_.push_back(group); }
    int nmbGroups() const
    { return (int)object_groups_.size(); }
    std::vector<gvGroup> objectGroups() const
    { return object_groups_; }
    gvGroup getGroup(int index) const
    { return object_groups_[index]; }
    void clearGroup()
    { object_groups_.clear(); }

    gvPainter& painter()
    { return painter_; }

    void registerObserver(gvObserver* obs)
    { observers_.insert(obs); }

    const Go::BoundingBox& boundingBox()
    { computeBox(); return box_; }

    Go::BoundingBox boundingBox(const std::vector<int> &objs) const;

    void extractSelectedObjects(std::vector<shared_ptr<Go::GeomObject> >& sel_objs);


    void setSelectedStateObject(int id, bool state);
    void setVisibleStateObject(int id, bool state);
    bool getSelectedStateObject(int id);
    bool getVisibleStateObject(int id);


    void disableUpdates()
    { updates_enabled_ = false; }
    void enableUpdates()
    { updates_enabled_ = true; }
    void updateObservers();

    void clear();

    void deleteObj(int obj);
    /// Removes pointers for last object only.
    void clearLast();

    /// Not always any need to retesselate all objects.
    void recreateDataStructure(int nmb_new_objs);

protected:
    std::vector< shared_ptr<Go::GeomObject> > objects_;
    std::vector< shared_ptr<gvColor> > object_colors_;
    gvPainter painter_;
private:
    Go::BoundingBox box_;

    std::vector<gvGroup> object_groups_; // As given by objects_. Two
					 // groups need not be
					 // disjoint.

    std::set< gvObserver* > observers_;

    std::vector< shared_ptr<Go::Tesselator> > tesselators_;
    std::vector< shared_ptr<gvPaintable> > paintables_;
    std::vector< shared_ptr<gvPropertySheet> > property_sheets_;

    bool updates_enabled_;
    std::auto_ptr<DataHandler> datahandler_;

    void computeBox();
};


#endif // _GVDATA_H

