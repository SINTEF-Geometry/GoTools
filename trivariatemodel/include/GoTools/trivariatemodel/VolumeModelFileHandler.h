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

#ifndef _VOLUMEMODELFILEHANDLER_H
#define _VOLUMEFILEHANDLER_H


#include "GoTools/compositemodel/CompositeModelFileHandler.h"
#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/trivariatemodel/VolumeModel.h"

namespace Go
{


// Writing to and reading from the g22 file format.
// Supports topology information (as opposed to the g2 format for geometries only).
  class VolumeModelFileHandler : public CompositeModelFileHandler
  {
    
  public:
    VolumeModelFileHandler();
 
    ~VolumeModelFileHandler();
    
    void writeVolume(const shared_ptr<ftVolume>& body,
		     std::ostream& os, int body_id=-1, 
		     bool faces=true);


    shared_ptr<ftVolume> readVolume(const char* filein, int id=-1);

    void writeVolumeModel(VolumeModel& vol_model,
			  std::ostream& os);

    shared_ptr<VolumeModel> readVolumeModel(const char* filein);

  private:
    std::map<shared_ptr<ftVolume>, int> volumes_;
 
    // Creator for our supported GeomObject's.
    virtual shared_ptr<GeomObject> 
      createGeomObject(const ObjectHeader& obj_header) const;
    
    void addVolumePointers(shared_ptr<ftVolume> body);
  };
  
  
} // namespace Go


#endif // _VOLUMEMODELFILEHANDLER_H
