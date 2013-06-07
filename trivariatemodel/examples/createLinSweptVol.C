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
#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/trivariatemodel/VolumeModelCreator.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"

//===========================================================================
///                                                                           
/// Description:
///  
/// Create a block structured volume model from a face set
/// describing a boundary represented solid that may be created by
/// sweeping a planar face set along a stright line
///
/// Input/Output:
///
/// Input to the geometry construction is the file data/two_holes_block.g2
/// The output volumes are stored in the file data/two_holes_block_out.g2
///                                                                           
//===========================================================================

using namespace Go;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::vector;

int main(int argc, char* argv[] )
{
  // Prepare for input and output files
  std::string infile("data/two_holes_block.g2");
  std::string outfile("data/two_holes_block_out.g2");

  ifstream input(infile.c_str());

  ofstream of(outfile.c_str());

  // Define tolerances
  // The neighbour tolerance is used in topology build if more than
  // one surface are given. Two surfaces lying closer than this tolerance
  // are assumed to be neighbours 
  double neighbour = 0.01;
  // Two surface lying closer than the neighbour tolerance, but more distant
  // than the gap tolerance are assumed to be neighbours, but the surface
  // set is not found to be C0 continuous within the given tolerance (gap)
  double gap = 0.001;
  // Two neighbouring surfaces where the angle between some corresponding
  // surface normals are larger than the bend tolerances are found to
  // create an intential corner. Angular tolerances are given in radians.
  double bend = 0.01;
  // If the angle between corresponding surface normals are larger than the
  // kink tolerance but smaller than the bend tolerance the surface set is 
  // found to be intentially, but not truely G1. If all angles are less than
  // the kink tolerance, the surface set is seen as G1
  double kink = 0.001;
  // Tolerance intended for approximations
  double approxtol = 0.001;

  // Create a factory class to read/create composite models, most often
  // one or more surfaces where the adjacency relationship between the
  // surfaces are known
  std::cout << "Reading input data" << std::endl;
  CompositeModelFactory factory(approxtol, gap, neighbour, kink, bend);

  // Read data from file. At this stage, it is not known whether the
  // model consists of curves or surfaces
  shared_ptr<CompositeModel> model(factory.createFromG2(input));

  // A surface model inherits composite model, check if we have
  // a surface model
  shared_ptr<SurfaceModel> sfmodel = 
    dynamic_pointer_cast<SurfaceModel,CompositeModel>(model);

if (sfmodel.get())
  {
    // Check if the model can be created by linear sweep. In that case
    // extract the faces at the top or bottom of the model and represent
    // the face set as a block structured model where each block is
    // a NURBS surface. Finally, perform the sweep to create volume
    // blocks
    shared_ptr<VolumeModel> volmod;
    bool linswept = VolumeModelCreator::linearSweptModel(sfmodel, volmod);
    std::cout << "Linswept: " << linswept << std::endl;

    if (linswept)
      {
	// The model can be created by linear sweep and a result is
	// obtained. Write the result to file. Note that the file contains
	// volumes. To visualize them use trivariate/app/makeShields to
	// extract the boundary surfaces of each volume
        int nmb_vols = volmod->nmbEntities();
	for (int kr=0; kr<nmb_vols; ++kr)
	  {
	    shared_ptr<ParamVolume> curr_vol2 = volmod->getVolume(kr);
	    curr_vol2->writeStandardHeader(of);
	    curr_vol2->write(of);
	  }
    }
  }

}

