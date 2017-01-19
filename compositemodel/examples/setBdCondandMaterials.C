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


#include "GoTools/compositemodel/Body.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/compositemodel/CompositeModelFileHandler.h"
#include "GoTools/geometry/Plane.h"
#include "GoTools/geometry/Line.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/BoundedSurface.h"
#include <fstream>

using std::vector;
using namespace Go;

//===========================================================================
///                                                                           
/// Description:
///  
/// Creation of a toy Brep solid containing spline surfaces and trimmed
/// planes. Some faces are equipped with boundary condition tags (one for
/// type of boundary condition and one for the actual condition) and the
/// solid with a material tag.
/// The solid is written to a g22 file and read in again.
///
///                                                                           
//===========================================================================


int main( int argc, char* argv[] )
{
  double eps = 1.0e-6;

  vector<shared_ptr<ParamSurface> > face_sfs(6);

  // Create planes
  Point pos1(1.0, 0.0, 3.0);
  Point pos2(0.0, 0.0, 0.0);
  Point normal(0.0, 0.0, 1.0);
  shared_ptr<ParamSurface> plane1(new Plane(pos1, normal));
  shared_ptr<ParamSurface> plane2(new Plane(pos2, normal, true));  // Normal pointing downwards
  
  // Create trimming curves for the first plane in geometry space
  Point dir1(1.0, 0.0, 0.0);
  Point pos3(5.0, 0.0, 3.0);
  Point pos4(5.0, 3.0, 3.0);
  Point pos5(1.0, 3.0, 3.0);
  shared_ptr<ParamCurve> line1(new Line(pos1, dir1, 4.0));
  shared_ptr<ParamCurve> line2(new Line(pos5, dir1, 4.0, true));
  shared_ptr<ParamCurve> spline1(new SplineCurve(pos3, pos4));
  shared_ptr<ParamCurve> spline2(new SplineCurve(pos5, pos1));
  
  // Create CurveOnSurface entities for the first plane
  vector<shared_ptr<CurveOnSurface> > loop_cvs1(4);
  loop_cvs1[0] = shared_ptr<CurveOnSurface>(new CurveOnSurface(plane1, line1,
							       false));
  loop_cvs1[1] = shared_ptr<CurveOnSurface>(new CurveOnSurface(plane1, spline1,
							       false));
  loop_cvs1[2] = shared_ptr<CurveOnSurface>(new CurveOnSurface(plane1, line2,
							       false));
  loop_cvs1[3] = shared_ptr<CurveOnSurface>(new CurveOnSurface(plane1, spline2,
							       false));

  // Create trimmed plane
  face_sfs[0] = shared_ptr<ParamSurface>(new BoundedSurface(plane1, loop_cvs1,
							    eps));

  // Create trimming curves for the first plane in geometry space
  Point pos6(0.0, 4.0, 0.0);
  Point pos7(4.0, 4.0, 0.0);
  Point pos8(4.0, 0.0, 0.0);
  shared_ptr<ParamCurve> spline3(new SplineCurve(pos2, pos6));
  shared_ptr<ParamCurve> spline4(new SplineCurve(pos6, pos7));
  shared_ptr<ParamCurve> spline5(new SplineCurve(pos7, pos8));
  shared_ptr<ParamCurve> spline6(new SplineCurve(pos8, pos2));
  
  // Create CurveOnSurface entities for the first plane
  vector<shared_ptr<CurveOnSurface> > loop_cvs2(4);
  loop_cvs2[0] = shared_ptr<CurveOnSurface>(new CurveOnSurface(plane2, spline3,
							       false));
  loop_cvs2[1] = shared_ptr<CurveOnSurface>(new CurveOnSurface(plane2, spline4,
							       false));
  loop_cvs2[2] = shared_ptr<CurveOnSurface>(new CurveOnSurface(plane2, spline5,
							       false));
  loop_cvs2[3] = shared_ptr<CurveOnSurface>(new CurveOnSurface(plane2, spline6,
							       false));

  // Create trimmed plane
  face_sfs[1] = shared_ptr<ParamSurface>(new BoundedSurface(plane2, loop_cvs2,
							    eps));

  // Create non-trimmed linear spline surfaces for the remaining sides
  double knots[4] = {0.0, 0.0, 1.0, 1.0};
  double coefs1[12] = {1.0, 0.0, 3.0, 1.0, 3.0, 3.0,
		       0.0, 0.0, 0.0, 0.0, 4.0, 0.0};
  double coefs2[12] = {1.0, 0.0, 3.0,  0.0, 0.0, 0.0, 
		       5.0, 0.0, 3.0, 4.0, 0.0, 0.0};
  double coefs3[12] = {4.0, 4.0, 0.0, 0.0, 4.0, 0.0,
		       5.0, 3.0, 3.0, 1.0, 3.0, 3.0};
  double coefs4[12] = {4.0, 0.0, 0.0, 4.0, 4.0, 0.0,
		       5.0, 0.0, 3.0, 5.0, 3.0, 3.0};

  face_sfs[2] = shared_ptr<ParamSurface>(new SplineSurface(2, 2, 2, 2, 
							   &knots[0],
							   &knots[0], 
							   &coefs1[0], 3));
  face_sfs[3] = shared_ptr<ParamSurface>(new SplineSurface(2, 2, 2, 2, 
							   &knots[0],
							   &knots[0], 
							   &coefs2[0], 3));
  face_sfs[4] = shared_ptr<ParamSurface>(new SplineSurface(2, 2, 2, 2, 
							   &knots[0],
							   &knots[0], 
							   &coefs3[0], 3));
  face_sfs[5] = shared_ptr<ParamSurface>(new SplineSurface(2, 2, 2, 2, 
							   &knots[0],
							   &knots[0], 
							   &coefs4[0], 3));

  // Test surface configuration
  std::ofstream of("solid_sfs.g2");
  for (int ki=0; ki<6; ++ki)
    {
      face_sfs[ki]->writeStandardHeader(of);
      face_sfs[ki]->write(of);
    }

  // Create faces and attach boundary condition flags for some surfaces
  vector<shared_ptr<ftSurface> > faces(6);
  for (int ki=0; ki<6; ++ki)
    {
      faces[ki] = shared_ptr<ftSurface>(new ftSurface(face_sfs[ki], ki));
      if (ki % 2 == 0)
	faces[ki]->setBoundaryConditions(ki, ki+1);
    }

  // Create solid. First create the outer shell
  // Topology tolerances
  double gap = 1.0e-6;
  double neighbour = 1.0e-3;
  double kink = 0.01;
  double bend = 0.05;
  shared_ptr<SurfaceModel> shell(new SurfaceModel(gap, gap, neighbour, kink,
						  bend, faces));

  // Create solid and set material flag to 1
  shared_ptr<Body> solid(new Body(shell, 1));
  
  // Write solid to file
  std::ofstream out_file("solid.g22");
  CompositeModelFileHandler filehandler;
  filehandler.writeStart(out_file);
  filehandler.writeHeader("Solid with material flag and boundary condition flags", out_file);
  filehandler.writeBody(solid, out_file);
  filehandler.writeEnd(out_file);

  // Read the file back
  CompositeModelFileHandler filehandler2;
  shared_ptr<Body> solid2 = filehandler2.readBody("solid.g22");

  // Fetch information and check that it is the same geometry
  if (solid2->hasMaterialInfo())
    {
      int material = solid2->getMaterial();
      std::cout << "Material flag: " << material << std::endl;
    }
  else
    std::cout << "The solid contains no material information" << std::endl;

  std::ofstream of2("solid_sfs2.g2");
  shared_ptr<SurfaceModel> shell2 = solid2->getOuterShell();
  int nmb_faces = shell2->nmbEntities();
  std::cout << "Number of faces: " << nmb_faces << std::endl;
  
  // Check that it really is a closed shell
  int nmb_bd = shell2->nmbBoundaries();
  std::cout << "Number of boundaries: " << nmb_bd << std::endl;

  for (int ki=0; ki<nmb_faces; ++ki)
    {
      shared_ptr<ParamSurface> surf = shell2->getSurface(ki);
      shared_ptr<ftSurface> face = shell2->getFace(ki);
      if (face->hasBoundaryConditions())
	{
	  int bd_cond_type, bd_cond;
	  face->getBoundaryConditions(bd_cond_type, bd_cond);
	  std::cout << "Face number " << ki << " has boundary condition: ";
	  std::cout << bd_cond_type << ", " << bd_cond << std::endl;
	}
      else 
	std::cout << "Face number " << ki << " has no boundary condition" << std::endl;
	
      surf->writeStandardHeader(of2);
      surf->write(of2);
    }
  
}

