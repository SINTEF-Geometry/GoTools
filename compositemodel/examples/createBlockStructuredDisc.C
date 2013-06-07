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
/// Create a block structured set of spline surfaces from a face set
/// consisting of possibly trimmed surfaces with arbitrary topology
/// (no corner-to-corner conditions)
///
/// Input/Output:
///
/// Input to the geometry construction is the file data/split_disc.g2
/// Current surfaces are written to g2-files as we go along and the final
/// version is stored in the file data/split_disc4.g2
///                                                                           
//===========================================================================

int main( int argc, char* argv[] )
{
  // Prepare for input and output files
  std::string infile("data/split_disc.g2");
  std::string outfile1("data/split_disc2.g2");
  std::string outfile2("data/split_disc3.g2");
  std::string outfile3("data/split_disc4.g2");
  
  // Prepare for input data
  std::ifstream input(infile.c_str());

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
    // Create class for splitting a given face set into a face set where
    // all associated surfaces has 4 boundaries
    std::cout << "Replace current surfaces with 4-sided surfaces"  << std::endl;
    RegularizeFaceSet reg(sfmodel);

    // Perform splitting
    // Note that this functionality is under developement and unexpected results
    // may occur. 
    shared_ptr<SurfaceModel> model2 = reg.getRegularModel();
    int nmb = model2->nmbEntities();
    std::ofstream of1(outfile1.c_str());
    for (int ki=0; ki<nmb; ++ki)
      {
	shared_ptr<ParamSurface> sf = model2->getSurface(ki);
	sf->writeStandardHeader(of1);
	sf->write(of1);
      }

    // Replace by spline surfaces
    // Boundary curves are approximated by spline curves, or fetched from
    // an already created adjacent spline surface. Spline surfaces are 
    // initially created as Coons patches interpolating the boundary curves
    // then updated with respect a point set fetched from the initial 
    // trimmed surface
    // The gap tolerance is used in the curve and surface approximations
    std::cout << "Replace trimmed surfaces by spline surfaces"  << std::endl;
    model2->replaceRegularSurfaces();

    // Write current surface set to a file
    std::ofstream of2(outfile2.c_str());
    nmb = model2->nmbEntities();
    for (int ki=0; ki<nmb; ++ki)
      {
	shared_ptr<ParamSurface> sf = model2->getSurface(ki);
	sf->writeStandardHeader(of2);
	sf->write(of2);
      }

    // Make sure that neighbouring surfaces have the same spline space
    // and coincidence of corresponding coefficients
    // First check if all surfaces are non-trimmed spline surfaces, i.e.
    // the previous opeation succeeded
    bool isOK = model2->allSplines();
    if (!isOK)
      {
	std::cout << "Not all surfaces are splines. Stopping computation" << std::endl;
	exit(-1);
      }
    
    // Ensure common spline spaces and corresponding coefficients
    std::cout << "Ensure common spline space"  << std::endl;
    model2->makeCommonSplineSpaces();

    // Write result to file
    std::ofstream of3(outfile3.c_str());
    for (int ki=0; ki<nmb; ++ki)
      {
	shared_ptr<ParamSurface> sf = model2->getSurface(ki);
	sf->writeStandardHeader(of3);
	sf->write(of3);
      }

  }
}
