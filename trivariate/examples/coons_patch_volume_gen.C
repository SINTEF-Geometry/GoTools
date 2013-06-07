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

#include <fstream>
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/errormacros.h"
#include "GoTools/trivariate/CoonsPatchVolumeGen.h"

using namespace Go;
using namespace std;

//===========================================================================
//                                                                           
// File: coons_patch_volume_gen.C                                            
//                                                                           
//                                                                           
/// Description:
///
/// This program demonstrates the use of the static function 'createCoonsPatch'
/// in namespace 'CoonsPatchVolumeGen'.
/// The function can create a new 'SplineVolume' representing the coons patch of
/// six SplineSurfaces, the six faces of the volume.
///
/// This program reads a non-rational 'SplineVolume' object from file and extract
/// the six isosurfaces defined by the min and max parameter values in u, v and w
/// direction. The surfaces must be non-rational, lie in the same space, and have
/// conciding edge curves where needed. 
/// The B-spline bases of the surfaces might be reversed, swapped, have order
/// raised, knots inserted and parameter interval rescaled to [0,1].
/// The input filename is hardcoded to 'vol2.g2'. (Located in 'trivariate/data/')
///
/// Output is a file in Go-format. The file name is hard-coded to
/// 'coons_patch_volume.g2'. The program 'goview' can't display volumes, but you
/// can use the programs 'makeShield' or 'getBoundarySfs' to extract the boundary
/// faces. They write a new file which can be used by 'goview'.
/// Both programs have inputfilename and outputfilename as arguments, but
/// 'makeShield' has a third optional parameter. If this third parameter is set
/// to 0 (zero), and the file has more then one volume, the volumes will be
/// displayed in different colours.
/// 
//===========================================================================

int main(int argc, char* argv[] )
{
    cout << "\nRunning program " << argv[0]
	 << ".\nFile 'trivariate/data/vol2.g2' used as input file." << endl;

  // Open input volume file.
  ifstream is("vol2.g2");
  ALWAYS_ERROR_IF(is.bad(), "Bad or no volume input filename");

  // Read volume from file.
  ObjectHeader head;
  SplineVolume vol;
  is >> head >> vol;
  is.close();
  cout << "Old volume. Bounding box = " << vol.boundingBox() << endl;
  cout << "Volume is rational?    " << boolalpha << vol.rational() << endl;
  cout << "Volume is left handed? " << boolalpha << vol.isLeftHanded() << endl;

  // Extract the six boundary surfaces.
  vector<SplineSurface*> faces;
  faces.push_back(vol.constParamSurface(vol.startparam(0), 0));
  faces.push_back(vol.constParamSurface(vol.endparam(0), 0));
  faces.push_back(vol.constParamSurface(vol.startparam(1), 1));
  faces.push_back(vol.constParamSurface(vol.endparam(1), 1));
  faces.push_back(vol.constParamSurface(vol.startparam(2), 2));
  faces.push_back(vol.constParamSurface(vol.endparam(2), 2));

  // Create a new SplineVolume. Vector 'faces' are not changed. A copy is used.
  SplineVolume* volnew = Go::CoonsPatchVolumeGen::createCoonsPatch
    (faces[0], faces[1], faces[2], faces[3], faces[4], faces[5]);
  cout << "\nNew volume. Bounding box = " << volnew->boundingBox() << endl;
  cout << "Volume is rational?    " << boolalpha << volnew->rational() << endl;
  cout << "Volume is left handed? " << boolalpha << volnew->isLeftHanded() << endl;

  // Open output volume file.
  ofstream os("coons_patch_volume.g2");

  // Write the new volume to file.
  volnew->writeStandardHeader(os);
  volnew->write(os);
  os.close();

  // cout << "\nRun: makeShield coons_patch_volume.g2 coons_patch_volume_sf.g2\n"
  //      << "and open the file 'coons_patch_volume_sf.g2' in 'goview' to look at"
  //      << " the result.\n"
  //      << endl;

  // Clean up
  for (size_t i = 0; i < faces.size(); ++i) {
    delete faces[i];
  }
  delete volnew;

  return 0;
}

