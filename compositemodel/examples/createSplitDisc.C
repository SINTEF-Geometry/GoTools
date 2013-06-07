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

#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SweepSurfaceCreator.h"
#include <fstream>

using namespace Go;
using std::cout;



//===========================================================================
//                                                                           
/// Description:
///  
/// This programs creates a face set representing a disc as two trimmed
/// surfaces: The trimmed disc and a rectangular surface lying inside this
/// disc. This is a starting point for the creation of a disc as a multi patch
/// model with spline surfaces and no degeneracies.
/// The construction uses planar, rectangular surfaces and a truncated sylinder,
/// but the operations performed using these surfaces, do not depend on that
/// level of regularity.
///
/// Input/Output:
///
/// Input to the geometry construction is hardcoded
/// Current surfaces written to g2-files as we go along. The final surface set
/// is written to the file data/split_disc.g2
//                                                                           
//===========================================================================
int main( int argc, char* argv[] )
{
  // Prepare for output files
  std::string outfile1("data/rect_sf.g2");
  std::string outfile2("data/cylinder_sf.g2");
  std::string outfile3("data/disc.g2");
  std::string outfile4("data/rect_tube.g2");
  std::string outfile5("data/split_disc.g2");
  
  // First a rectangular spline surface is created as the base surface
  // for a disc represented as a bounded surface
  // Define data
  std::cout << "Creating underlying surface" << std::endl;
  int dim = 3;
  int numcoefs = 2;
  int order = 2;
  double knots[4] = {0.0, 0.0, 1.0, 1.0};
  double coefs[12] = {-2.0, -2.0, 0.0,
		      2.0, -2.0, 0.0,
		      -2.0, 2.0, 0.0,
		      2.0, 2.0, 0.0};

  // Create spline rectangle
  shared_ptr<SplineSurface> rect(new SplineSurface(numcoefs, numcoefs,
						   order, order,
						   knots, knots, coefs,
						   dim));

  std::ofstream of1(outfile1.c_str());
  rect->writeStandardHeader(of1);
  rect->write(of1);

  // Create cylinder
  std::cout << "Creating cylinder " << std::endl;
  double radius = 1.5;
  Point centre(0.0, 0.0, 0.0);
  Point z_axis(0.0, 0.0, 1.0);
  Point x_axis(1.0, 0.0, 0.0);

  shared_ptr<Cylinder> cylinder(new Cylinder(radius, centre, z_axis, x_axis));

  // Restrict the cylinder height. The 1. parameter directions is
  // running from 0 to 2*pi, while the 2. is infinite
  cylinder->setParameterBounds(0.0, -1, 2*M_PI, 1);

  // Represent the cylinder as a NURBS surface
  shared_ptr<SplineSurface> cyl(cylinder->createSplineSurface());

  std::ofstream of2(outfile2.c_str());
  cyl->writeStandardHeader(of2);
  cyl->write(of2);

  // Represent each surface as one surface set to perform a Boolean operation
  // between these two surface set
  // First define topology tolerances
  double neighbour = 0.01;
  double gap = 0.0001;
  double bend = 0.5;
  double kink = 0.01;
  // Approximation tolerance. Not used in this context
  double approxtol = 0.01;

  vector<shared_ptr<ParamSurface> > set1(1);
  set1[0] = rect;
  shared_ptr<SurfaceModel> model1(new SurfaceModel(approxtol, gap, neighbour,
						   kink, bend, set1));
  vector<shared_ptr<ParamSurface> > set2(1);
  set2[0] = cyl;
  shared_ptr<SurfaceModel> model2(new SurfaceModel(approxtol, gap, neighbour,
						   kink, bend, set2));
  
  // Split both surface models with respect to the intersection curves
  // between them, and fetch the circular disc
  cout << "Intersecting base rectangle with cylinder" << std::endl;
  vector<shared_ptr<SurfaceModel> > sub_models = 
    model1->splitSurfaceModels(model2);
  if (sub_models.size() != 4)
    {
      std::cout << "Unexpected number of pieces" << std::endl;
      exit(-1);
    }

  // By construction, the number of models after split should be 4 and
  // the first one is the part of the first model lying inside the
  // second surface set
  shared_ptr<SurfaceModel> disc = sub_models[0];

  std::ofstream of3(outfile3.c_str());
  disc->getSurface(0)->writeStandardHeader(of3);
  disc->getSurface(0)->write(of3);

  // Create a quadratic tube intersecting the disc in the inner in order
  // to split the disc into two parts
  // Start by creating the lower boundary curve in each of the surface
  // defining the tube
  std::cout << "Creating rectangular cube as swept surfaces" << std::endl;
  vector<Point> pnts(4);
  pnts[0] = Point(-0.5, -0.5, -1);
  pnts[1] = Point(0.5, -0.5, -1);
  pnts[2] = Point(0.5, 0.5, -1);
  pnts[3] = Point(-0.5, 0.5, -1);
  
  // Each curve is a linear spline curve on the parameter interval
  // [0,1] which interpolates two tube corners
  vector<shared_ptr<SplineCurve> > cvs(4);
  int ki, kj;
  for (ki=0; ki<4; ++ki)
    {
      kj = (ki+1)%4;
      cvs[ki] = 
	shared_ptr<SplineCurve>(new SplineCurve(pnts[ki], pnts[kj]));
    }

  // Perform a linear sweep of each curve in the direction of the
  // cylinder axis to create surfaces
  // Create swept surfaces
  SweepSurfaceCreator sweep;
  vector<shared_ptr<ParamSurface> > swept_sfs(4);
  for (ki=0; ki<4; ++ki)
    {
      shared_ptr<SplineCurve> axis(new SplineCurve(pnts[ki],
						   pnts[ki]+2.0*z_axis));
      swept_sfs[ki] = 
	shared_ptr<SplineSurface>(sweep.linearSweptSurface(*cvs[ki],
							   *axis, 
							   pnts[ki]));
    }

  std::ofstream of4(outfile4.c_str());
  for (ki=0; ki<4; ++ki)
    {
      swept_sfs[ki]->writeStandardHeader(of4);
      swept_sfs[ki]->write(of4);
    }

  // Create surface model containing the swept surfaces
  shared_ptr<SurfaceModel> model3(new SurfaceModel(approxtol, gap, neighbour,
						   kink, bend, swept_sfs));

  // Split the dist with respect to the tube model and the other way around
  std::cout << "Split disc by intersecting with rectangular tube" << std::endl;
  vector<shared_ptr<SurfaceModel> > sub_models2 = 
    disc->splitSurfaceModels(model3);
  if (sub_models2.size() != 4)
    {
      std::cout << "Unexpected number of pieces" << std::endl;
      exit(-1);
    }
 
  // Now we want to keept both pieces of the disc surface. They are stored 
  // the first two surface models. Join both parts into one model
  sub_models2[0]->append(sub_models2[1]);
  
  std::ofstream of5(outfile5.c_str());
  int nmb = sub_models2[0]->nmbEntities();
  for (ki=0; ki<nmb; ++ki)
    {
      shared_ptr<ParamSurface> sf = sub_models2[0]->getSurface(ki);
      sf->writeStandardHeader(of5);
      sf->write(of5);
    }
}

