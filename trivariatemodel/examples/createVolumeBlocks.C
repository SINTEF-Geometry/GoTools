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

#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/trivariate/ParamVolume.h"
#include "GoTools/compositemodel/RegularizeFaceSet.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/geometry/SplineSurface.h"
#include <fstream>

using namespace Go;
using std::cout;

//===========================================================================
///                                                                           
/// Description:
///  
/// Create a block structured volume model from a face set
/// consisting of possibly trimmed surfaces with arbitrary topology
/// (no corner-to-corner conditions). The face set represents a boundary
/// represented solid and is described in a g2-file.
///
/// Input/Output:
///
/// Input to the geometry construction is the file data/hole_block.g2
/// Current surfaces are written to g2-files as we go along and the final
/// version is stored in the file data/hole_block_out.g2
///                                                                           
//===========================================================================

int main( int argc, char* argv[] )
{
  // Prepare for input and output files
  std::string infile("data/hole_block.g2");
  std::string outfile1("data/trimmed_volumes_sfs.g2");
  std::string outfile2("data/output_volumes.g2");
  std::string outfile3("data/output_model.g2");
  
  // Prepare for input data
  std::ifstream input(infile.c_str());

  // Prepare for output data
  std::ofstream of1(outfile1.c_str());
  std::ofstream of2(outfile2.c_str());
  std::ofstream of3(outfile3.c_str());
  
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

  int degree = 3;  // Make cubic volumes
  vector<SurfaceModel*> modified_adjacent;  // Dummy vector

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
    // Create a trimmed volume
    shared_ptr<ftVolume> ftvol = 
      shared_ptr<ftVolume>(new ftVolume(sfmodel));

    // Check if the trimmed volume has 6 boundary surfaces with
    // 4 boundary curves each
    int nmb;
    int ki;
    shared_ptr<VolumeModel> volmod;
    bool reg = ftvol->isRegularized();
    if (!reg)
      {
	// The current trimmed volume cannot be represented/approximated 
	// by one NURBS volume. Split the current volume into a set of
	// trimmed blocks where each can be represented by a NURBS volume
	// In the first attempt only one approach for extending the
	// wireframe of the boundary surfaces to the wireframe of the
	// volumemodel is used. Normally, this is enough as the
	// alternative approach is invoked automatically after a first
	// split of the trimmed volume. The parameter false indicates
	// that the second approach is not to be used initially
	vector<shared_ptr<ftVolume> > reg_vols = 
	  ftvol->replaceWithRegVolumes(degree, modified_adjacent, false);

	// Assemble the volume blocks into a volume model
      if (reg_vols.size() > 0)
	// A number of blocks is found
	volmod = shared_ptr<VolumeModel>(new VolumeModel(reg_vols, gap, 
							 neighbour,
							 kink, 10.0*kink));
      else
	{
	  // No block structuring has been performed. Either, there is only
	  // one block after some merging of surfaces, the block 
	  // structuring failed or there is a special configuration where
	  // that needs the second approach for extending the wireframe
	  // and where no initial split is performed. Try again using
	  // the second approach.
	  vector<shared_ptr<ftVolume> > reg_vols2 = 
	    ftvol->replaceWithRegVolumes(degree, modified_adjacent, true);
	  if (reg_vols2.size() > 0)
	    // A number of blocks is found
	    volmod = shared_ptr<VolumeModel>(new VolumeModel(reg_vols2, gap, 
							     neighbour,
							     kink, 10.0*kink));
	  else
	    {
	      vector<shared_ptr<ftVolume> > reg_vols2(1);
	      reg_vols2[0] = ftvol;
	      volmod = shared_ptr<VolumeModel>(new VolumeModel(reg_vols2, 
							       gap, neighbour,
							       kink, 
							       10.0*kink));
	    }    
	}
      }
    else
      {
	// The initial trimmed volume can be replaced by a NURBS block
	vector<shared_ptr<ftVolume> > reg_vols(1);
	reg_vols[0] = ftvol;
	volmod = shared_ptr<VolumeModel>(new VolumeModel(reg_vols, gap, 
							 neighbour,
							 kink, 10.0*kink));
      }


  std::cout << "Number of volumes: " << volmod->nmbEntities() << std::endl;
	  
  int nmb_vols = volmod->nmbEntities();
  for (int kr=0; kr<nmb_vols; ++kr)
    {
      // For each trimmed volume, check if it is iso-trimmed, i.e. 
      // is trimmed along isoparametric surfaces, or regular, i.e. can
      // be replaced by a NURBS block
      shared_ptr<ftVolume> curr_vol = volmod->getBody(kr);
      bool bd_trim = curr_vol->isBoundaryTrimmed();
      bool iso_trim = curr_vol->isIsoTrimmed();
      bool reg = curr_vol->isRegularized();

      std::cout << "Volume nr " << kr << ": " << bd_trim;
      std::cout << " " << iso_trim << " " << reg << std::endl;

      shared_ptr<SurfaceModel> mod = curr_vol->getOuterShell();
      nmb = mod->nmbEntities();
      for (ki=0; ki<nmb; ++ki)
	{
	  shared_ptr<ParamSurface> sf = mod->getSurface(ki);
	  sf->writeStandardHeader(of1);
	  sf->write(of1);
	}

      if (reg)
	{
	  // The trimmed volume is regular. Replace with a NURBS block
	  curr_vol->untrimRegular(degree);

	  // Fetch the NURBS block
	  shared_ptr<ParamVolume> curr_vol2 = volmod->getVolume(kr);

	  // Write the blocks to a file. Note that the file will contain
	  // volumes. To visualize it in the viewer, the application program
	  // makeShields in the module trivariate is used. The parameters
	  // to this program is the infile, the outfile and possibly an
	  // integer indicating if different colours is to be used for the
	  // boundary surfaces if the different volumes. 0 is a good choice
	  // for this integer
	  curr_vol2->writeStandardHeader(of2);
	  curr_vol2->write(of2);
	}
    }

  // Ensure that the blocks in the volume model have corresponding spline
  // spaces along common boundary surfaces
  volmod->makeCommonSplineSpaces();
  volmod->averageCorrespondingCoefs();

  nmb_vols = volmod->nmbEntities();
  for (int kr=0; kr<nmb_vols; ++kr)
    {
      // Write the final volume blocks to a file
      shared_ptr<ParamVolume> curr_vol2 = volmod->getVolume(kr);
      curr_vol2->writeStandardHeader(of3);
      curr_vol2->write(of3);
    }

  }
}

