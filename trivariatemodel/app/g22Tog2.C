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

#include <iostream>
#include <fstream>
#include <string>

#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/trivariatemodel/VolumeModelFileHandler.h"

using namespace Go;
using std::vector;

// Convert a set of volumes or surface models to a g2 file. Currently handling trimmed sfs only. Not including any
// curves that are not used by the trimmed sfs. Volumes not yet supported.

int main( int argc, char* argv[] )
{
  if (argc != 3) {
    std::cout << "Input parameters : Input file(.g22), output file(.g2)" << std::endl;
    exit(-1);
  }

  std::string file1_name(argv[1]);
  std::ifstream file1(file1_name);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  std::ofstream out_file(argv[2]);

  VolumeModelFileHandler filehandler;
  std::vector<shared_ptr<ftVolume> > volumes = filehandler.readVolumes(file1_name.c_str());
  if (volumes.size() > 0) {
      std::cout << "Not yet supporting volumes: volumes.size(): " << volumes.size() << std::endl;
  }

  std::vector<shared_ptr<SurfaceModel> > sf_mods = filehandler.readSurfModels(file1_name.c_str());
  std::cout << "sf_mods.size(): " << sf_mods.size() << std::endl;
  for (auto surf_model : sf_mods)
  {
      std::vector<shared_ptr<ftSurface> > faces = surf_model->allFaces();
      std::cout << "Number of sfs: " << faces.size() << std::endl;
      for (auto face : faces)
      {
          shared_ptr<ParamSurface> sf = face->surface();
          sf->writeStandardHeader(out_file);
          sf->write(out_file);
      }
  }

}
