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

#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/trivariate/SweepVolumeCreator.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include <fstream>

using namespace Go;
using std::cout;


//===========================================================================
//                                                                           
/// Description:
///  
/// The idea of this program is to sweep a set surfaces to create a multi patch 
/// volume model.
/// The surface set will be created by the example programs in compositemodel:
/// createSplitDisc and createBlockStructuredDisc
///
/// Input/Output:
///
/// Input is the file data/split_disc4.g2
/// The output volumes are written to data/swept_vol.g2
///                                                                           
//===========================================================================


int main( int argc, char* argv[] )
{
  // Prepare for input and output files
  std::string infile("data/split_disc4.g2");
  std::string outfile("data/swept_vol.g2");
  std::ifstream input(infile.c_str());

  // Read input data
  // Define tolerances
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
  // Tolerance intended for approximations
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
  
    // Make sure that neighbouring surfaces have the same spline space
    // and coincidence of corresponding coefficients
    // First check if all surfaces are non-trimmed spline surfaces, i.e.
    // the previous opeation succeeded
    bool isOK = sfmodel->allSplines();
    if (!isOK)
      {
	std::cout << "Not all surfaces are splines. Stopping computation" << std::endl;
	exit(-1);
      }
    
    // Create linear swept volumes
    std::cout << "Create swept volumes" << std::endl;
    vector<shared_ptr<ftVolume> > blocks;  // Storage for volumes including
                                           // topology information
    double sweep_length_fac = 5.0;
    // The centre and z-axis should correspond to the corresponding parameters
    // in compositemodel/examples/createSplitDisc
    Point centre(0.0, 0.0, 0.0);  
    Point z_axis(0.0, 0.0, 1.0);
    int nmb = sfmodel->nmbEntities();
    for (int ki=0; ki<nmb; ++ki)
      {
	// Fetch current spline surface
	shared_ptr<ParamSurface> sf = sfmodel->getSurface(ki);
	shared_ptr<SplineSurface> surf = 
	  dynamic_pointer_cast<SplineSurface,ParamSurface>(sf);

	if (surf.get())
	  {
	    // Create swept volume
	    shared_ptr<SplineCurve> cv(new SplineCurve(centre, 
						       centre+sweep_length_fac*z_axis));
	    shared_ptr<ParamVolume> vol = 
	      shared_ptr<ParamVolume>(SweepVolumeCreator::linearSweptVolume(*surf, 
									    *cv, 
									    centre));
	    // Add data structures for volume topology
	    shared_ptr<ftVolume> ftvol = 
	      shared_ptr<ftVolume>(new ftVolume(vol, gap, kink));
	    blocks.push_back(ftvol);
	  }
      }
    
    // Create a volume model. The topology build performs coincidence testing
    // between possible adjacent volumes. This may be time consuming
    std::cout << "Volume topology build"  << std::endl;
    shared_ptr<VolumeModel> volmodel = 
      shared_ptr<VolumeModel>(new VolumeModel(blocks, gap, neighbour, 
					      kink, 10.0*kink));


    // Fetch the number of volumes
    int nmb_vol = volmodel->nmbEntities();

    // Write the volumes to the output file
    // With the current file format, the topology information is lost
    std::ofstream of(outfile.c_str());
    for (int ki=0; ki<nmb_vol; ++ki)
      {
	shared_ptr<ParamVolume> vol = volmodel->getVolume(ki);
	vol->writeStandardHeader(of);
	vol->write(of);
      }
  }
}

  
		      
