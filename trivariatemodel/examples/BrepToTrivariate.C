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
#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/RegularizeFaceSet.h"

using namespace Go;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::vector;

//===========================================================================
//                                                                           
/// Description:
///  
/// The idea of this program is to read a Brep model in g2-format and create
/// trivariate spline model.
/// Note that the functionality works only for some classes of Brep models.
///
/// Input/Output:
///
/// argument: mode=(1/2/3). Strategy for remodelling. For most models
/// mode=1 is the most reliable
/// Input is the file data/plate_with_holes.g2
/// The output volumes are written to data/plate_with_holes_vol2.g2
///                                                                           
//===========================================================================

int main(int argc, char* argv[] )
{
  if (argc != 2)
    {
      cout << "Usage: " << "Block structuring mode (1,2,3)" << endl;
      exit(-1);
    }

  // Prepare for input and output files
  std::string infile("data/plate_with_holes.g2");
  std::string outfile("data/plate_with_holes_vol2.g2");
  std::ifstream input(infile.c_str());

  // Block structuring mode
  int split_mode = atoi(argv[1]);
  if (split_mode < 1 || split_mode > 3)
    split_mode = 1;  // Default

  int degree = 3;  // Create cubic splines

  // Define tolerances
  // The tolerances must be set according to the properties of the model.
  // The neighbour tolerance must be smaller than the smallest entity in the
  // model, but larger than the largest gap.
  // The gap tolerance must be smaller than the neighbour tolerance
  // The neighbour tolerance is used in topology build if more than
  // one surface are given. Two surfaces lying closer than this tolerance
  // are assumed to be neighbours 
  double neighbour = 0.001;
  // Two surface lying closer than the neighbour tolerance, but more distant
  // than the gap tolerance are assumed to be neighbours, but the surface
  // set is not found to be C0 continuous within the given tolerance (gap)
  double gap = 0.0001;
  // Two neighbouring surfaces where the angle between some corresponding
  // surface normals are larger than the bend tolerances are found to
  // create an intential corner. Angular tolerances are given in radians.
  double bend = 0.01;
  // If the angle between corresponding surface normals are larger than the
  // kink tolerance but smaller than the bend tolerance the surface set is 
  // found to be intentially, but not truely G1. If all angles are less than
  // the kink tolerance, the surface set is seen as G1
  double kink = 0.001;
  // Tolerance intended for approximations. Not used in this context
  double approxtol = 0.0001;

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
      // Check input
      if (sfmodel->nmbBoundaries() > 0)
	{
	  std::cout << "Not a brep solid. Consider increasing the neighbour tolerance" << std::endl;
	  exit(-1);
	}
    
      // Create trimmed volume
      shared_ptr<ftVolume> ftvol = 
	shared_ptr<ftVolume>(new ftVolume(sfmodel));

      // Volume model for the volume blocks to be created
      shared_ptr<VolumeModel> volmod;

      // Check if the trimmed volume can be replaced by one non-trimmed one
      bool reg = ftvol->isRegularized();
      bool pattern_split = false; //true;
      if (!reg)
	{
	  // More than one block is required. Perform block structuring
	  vector<SurfaceModel*> modified_adjacent;
	  vector<shared_ptr<ftVolume> > reg_vols = 
	    ftvol->replaceWithRegVolumes(degree, modified_adjacent,
					 false, split_mode, pattern_split);

	  if (reg_vols.size() > 0)
	    // Block structuring performed
	    volmod = shared_ptr<VolumeModel>(new VolumeModel(reg_vols, gap, neighbour,
							     kink, 10.0*kink));
	  else
	    {
	      // No success. Return initial underlying volume
	      vector<shared_ptr<ftVolume> > reg_vols(1);
	      reg_vols[0] = ftvol;
	      volmod = shared_ptr<VolumeModel>(new VolumeModel(reg_vols, gap, neighbour,
							       kink, 10.0*kink));
	    }    
	}
      else
	{
	  // One block
	  vector<shared_ptr<ftVolume> > reg_vols(1);
	  reg_vols[0] = ftvol;
	  volmod = shared_ptr<VolumeModel>(new VolumeModel(reg_vols, gap, neighbour,
							   kink, 10.0*kink));
	}
      std::cout << "Number of volumes: " << volmod->nmbEntities() << std::endl;

      // Replace trimmed volumes with non-trimmed spline volumes
      int nmb_vols = volmod->nmbEntities();
      for (int kr=0; kr<nmb_vols; ++kr)
	{
	  shared_ptr<ftVolume> curr_vol = volmod->getBody(kr);
	  bool reg = curr_vol->isRegularized();
	  if (reg)
	    curr_vol->untrimRegular(degree);
	}

      // Ensure C0 continuity between blocks
      volmod->makeCommonSplineSpaces();
      volmod->averageCorrespondingCoefs();

      // Write to output file
      std::ofstream of(outfile.c_str());
      for (int kr=0; kr<nmb_vols; ++kr)
	{
	  shared_ptr<ParamVolume> curr_vol2 = volmod->getVolume(kr);
	  curr_vol2->writeStandardHeader(of);
	  curr_vol2->write(of);
	}
    }
}


