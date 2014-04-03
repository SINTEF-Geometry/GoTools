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
#include <fstream>

using namespace Go;
using std::cout;

int main( int argc, char* argv[] )
{
  // Prepare for output files
  std::ofstream of("curved_plate.g2");
  std::ofstream of0("tmp_plate.g2");
  
  // First a rectangular spline surface is created as the base surface
  // for a disc represented as a bounded surface
  // Define data
  std::cout << "Creating underlying surface" << std::endl;
  int dim = 3;
  int numcoefs1 = 4;
  int order1 = 4;
  int numcoefs2 = 2;
  int order2 = 2;
  double knots1[8] = {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};
  double knots2[4] = {0.0, 0.0, 1.0, 1.0};
  double coefs[24] = {-4.5, -3, 0,
		      -2, 2, 0,
		      2, 2, 0,
		      4.5, -3, 0,
		      -6, 0, 0,
		      -2, 6, 0,
		      2, 6, 0,
		      6, 0, 0};

  // Create spline rectangle
  shared_ptr<SplineSurface> surf(new SplineSurface(numcoefs1, numcoefs2,
						   order1, order2,
						   knots1, knots2, coefs,
						   dim));
  surf->writeStandardHeader(of0);
  surf->write(of0);

  // Create cylinders
  std::cout << "Creating cylinder " << std::endl;
  double radius = 0.9;
  Point centre1(-4.2, 0.1, 0.0);
  Point centre2(-1.7, 2, 0.0);
  Point centre3(1.7, 2, 0.0);
  Point centre4(4.2, 0.1, 0.0);
  Point z_axis(0.0, 0.0, 1.0);
  Point x_axis(1.0, 0.0, 0.0);

  shared_ptr<Cylinder> cyl1(new Cylinder(radius, centre1, z_axis, x_axis));
  cyl1->setParameterBounds(0.0, -1, 2*M_PI, 1);
  shared_ptr<SplineSurface> cyl1_2(cyl1->createSplineSurface());

  shared_ptr<Cylinder> cyl2(new Cylinder(radius, centre2, z_axis, x_axis));
  cyl2->setParameterBounds(0.0, -1, 2*M_PI, 1);
  shared_ptr<SplineSurface> cyl2_2(cyl2->createSplineSurface());

  shared_ptr<Cylinder> cyl3(new Cylinder(radius, centre3, z_axis, x_axis));
  cyl3->setParameterBounds(0.0, -1, 2*M_PI, 1);
  shared_ptr<SplineSurface> cyl3_2(cyl3->createSplineSurface());

  shared_ptr<Cylinder> cyl4(new Cylinder(radius, centre4, z_axis, x_axis));
  cyl4->setParameterBounds(0.0, -1, 2*M_PI, 1);
  shared_ptr<SplineSurface> cyl4_2(cyl4->createSplineSurface());

  cyl1_2->writeStandardHeader(of0);
  cyl1_2->write(of0);
  cyl2_2->writeStandardHeader(of0);
  cyl2_2->write(of0);
  cyl3_2->writeStandardHeader(of0);
  cyl3_2->write(of0);
  cyl4_2->writeStandardHeader(of0);
  cyl4_2->write(of0);

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
  set1[0] = surf;
  shared_ptr<SurfaceModel> model1(new SurfaceModel(approxtol, gap, neighbour,
						   kink, bend, set1));
  vector<shared_ptr<ParamSurface> > set2(1);
  set2[0] = cyl1_2;
  shared_ptr<SurfaceModel> model2(new SurfaceModel(approxtol, gap, neighbour,
						   kink, bend, set2));
  
  // Split both surface models with respect to the intersection curves
  // between them, and fetch the outer part
  cout << "Intersecting base rectangle with cylinder" << std::endl;
  vector<shared_ptr<SurfaceModel> > sub_models = 
    model1->splitSurfaceModels(model2);

  model1 = sub_models[1];
  set2[0] = cyl2_2;
  model2 = shared_ptr<SurfaceModel>(new SurfaceModel(approxtol, gap, neighbour,
						     kink, bend, set2));

  sub_models = model1->splitSurfaceModels(model2);

  model1 = sub_models[1];
  set2[0] = cyl3_2;
  model2 = shared_ptr<SurfaceModel>(new SurfaceModel(approxtol, gap, neighbour,
						     kink, bend, set2));

  sub_models = model1->splitSurfaceModels(model2);

  model1 = sub_models[1];
  set2[0] = cyl4_2;
  model2 = shared_ptr<SurfaceModel>(new SurfaceModel(approxtol, gap, neighbour,
						     kink, bend, set2));
  sub_models = model1->splitSurfaceModels(model2);

  shared_ptr<SurfaceModel> plate = sub_models[1];

  int nmb = plate->nmbEntities();
  for (int ki=0; ki<nmb; ++ki)
    {
      shared_ptr<ParamSurface> sf = plate->getSurface(ki);
      sf->writeStandardHeader(of);
      sf->write(of);
    }
}

