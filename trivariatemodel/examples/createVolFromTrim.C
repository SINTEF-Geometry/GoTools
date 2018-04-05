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
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/ParamVolume.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/trivariatemodel/ftVolumeTools.h"
#include "GoTools/trivariatemodel/VolumeModelFileHandler.h"
#include <fstream>

using namespace Go;
using std::vector;

//===========================================================================
//
/// Description: 
///
///              
///  
/// 
///
/// Input/Output: 
///
///               
/// 
/// Note:       
///
///             
///
//   
//===========================================================================

int main( int argc, char* argv[] )
{
  // Create spline volume as a linear block
  double knots1[4] = {4.0, 4.0, 7.0, 7.0};
  double knots2[4] = {4.0, 4.0, 7.0, 7.0};
  double knots3[4] = {4.0, 4.0, 7.0, 7.0};
  double coefs[24] = {0.0, 0.0, 0.0, 3.0, 0.0, 0.0,
		     0.0, 3.0, 0.0, 3.0, 3.0, 0.0,
		     0.0, 0.0, 3.0, 3.0, 0.0, 3.0,
		     0.0, 3.0, 3.0, 3.0, 3.0, 3.0};
  shared_ptr<SplineVolume> vol1(new SplineVolume(2, 2, 2, 2, 2, 2,
						 &knots1[0], &knots2[0],
						 &knots3[0], &coefs[0],
						 3));

  // Degree elevation
  vol1->raiseOrder(2, 2, 2);

  // Knot insertion in all parameter direction
  vector<double> internal_knots(2);
  internal_knots[0] = 5.0;
  internal_knots[1] = 6.0;
  vol1->insertKnot(0, internal_knots);
  vol1->insertKnot(1, internal_knots);
  vol1->insertKnot(2, internal_knots);

  // Create splitting surface 
  double knots_u[9] = {1.5707963267949, 1.5707963267949, 1.5707963267949, 1.5707963267949,  3.14159265358979, 4.71238898038469, 4.71238898038469, 4.71238898038469, 4.71238898038469};
  double knots_v[7] = {-1, -1, -1, 1.5, 4, 4, 4};
  double cfs[60] = {4.2, 5.5, -1, 
		    0.14142135623731, 3.88908729652601, -0.707106781186547, 
		    0.2, 1.5, -1, 
		    0.14142135623731, -1.76776695296637, -0.707106781186547,
		    4.2, -2.5, -1, 
		    4.2, 5.5, 1.5, 
		    0.14142135623731, 3.88908729652601, 1.06066017177982, 
		    0.2, 1.5, 1.5, 
		    0.44142135623731, -1.76776695296637, 1.06066017177982,
		    4.2, -2.5, 1.5, 
		    4.2, 5.5, 3, 
		    0.14142135623731, 3.88908729652601, 2.82842712474619, 
		    0.2, 1.5, 3, 
		    0.14142135623731, -1.76776695296637, 2.82842712474619,
		    4.2, -2.5, 3, 
		    4.2, 5.5, 4, 
		    0.14142135623731, 3.88908729652601, 4,
		    0.2, 1.5, 4, 
		    0.14142135623731, -1.76776695296637, 4,
		    4.2, -2.5, 4};
  shared_ptr<SplineSurface> sf1(new SplineSurface(5, 4, 4, 3,
						 &knots_u[0], &knots_v[0],
						 &cfs[0], 3));
  
  std::ofstream of1("tmp1.g2");
  vol1->writeStandardHeader(of1);
  vol1->write(of1);
  std::ofstream of2("tmp2.g2");
  sf1->writeStandardHeader(of2);
  sf1->write(of2);
  
  // Prepare for splitting
  double gap = 1.0e-4;
  double kink = 0.01;
  shared_ptr<ParamVolume> vol2 = vol1;
  shared_ptr<ftVolume> topvol(new ftVolume(vol2, gap, kink));
  shared_ptr<ftSurface> face(new ftSurface(sf1, 0));

  // Split volume
  double eps = 1.0e-4;
  vector<shared_ptr<ftVolume> > topvols2 = 
    ftVolumeTools::splitVolumes(topvol, face, eps);
  
  std::ofstream of3("tmp3.g2");
  for (size_t kj=0; kj<topvols2.size(); ++kj)
    {
      shared_ptr<SurfaceModel> mod = topvols2[kj]->getOuterShell();
      int nmb = mod->nmbEntities();
      int ki;
      for (ki=0; ki<nmb; ++ki)
	{
	  shared_ptr<ftSurface> face = mod->getFace(ki);
	  int bd_stat = ftVolumeTools::boundaryStatus(topvols2[kj].get(), 
						      face, eps);
	  std::cout << "Bd status " << kj << ", " << ki << ": " << bd_stat << std::endl;
	  shared_ptr<ParamSurface> sf = mod->getSurface(ki);
	  sf->writeStandardHeader(of3);
	  sf->write(of3);
	}
    }

  // Select the first volume and pass through all elements and check if
  // they intersect the non-boundary trimming surface
  shared_ptr<ftVolume> curr_vol = topvols2[0];
  
  std::ofstream out_file("volmodel2.g22");
  VolumeModelFileHandler filehandler;
  filehandler.writeStart(out_file);
  filehandler.writeHeader("Test ftVolume", out_file);
  filehandler.writeVolume(curr_vol, out_file);
  filehandler.writeEnd(out_file);

  VolumeModelFileHandler filehandler2;
  shared_ptr<ftVolume> curr_vol2 = filehandler2.readVolume("volmodel2.g22");

  // Number of elements in underlying volume
  shared_ptr<SplineVolume> curr_under = 
    dynamic_pointer_cast<SplineVolume>(curr_vol2->getVolume());
  int nmb_elem = curr_under->numElem();
  std::cout << "No of elements: " << nmb_elem << std::endl;

  int degree = 3;
  std::ofstream of5("tmp5.g2");
  std::ofstream of6("tmp6.g2");
  for (int ki=0; ki<nmb_elem; ++ki)
    {
      int elem_stat = curr_vol2->ElementBoundaryStatus(ki);
      std::cout << "Boundary status, element " << ki+1 << ": " << elem_stat << std::endl;

      if (elem_stat == 1)
	{
	  // Element intersects trimming surface
	  // Split element with trimming shell
	  vector<shared_ptr<ftVolume> > sub_elem;
	  vector<int> is_inside; // Equal 1 if the sub element is inside the trimmed volume
	  curr_vol2->splitElementByTrimSfs(ki, eps, sub_elem, is_inside);

	  std::ofstream of4("tmp4.g2");
	  for (size_t kj=0; kj<sub_elem.size(); ++kj)
	    {
	      shared_ptr<SurfaceModel> mod = sub_elem[kj]->getOuterShell();
	      int nmb = mod->nmbEntities();
	      for (int kr=0; kr<nmb; ++kr)
		{
		  shared_ptr<ParamSurface> sf = mod->getSurface(kr);
		  sf->writeStandardHeader(of4);
		  sf->write(of4);
		}
	    }
	  int stop_break = 1;

 	  // Check if the remaining element is hexagonal
	  for (size_t kj=0; kj<sub_elem.size(); ++kj)
	    {
	      bool regular = sub_elem[kj]->isRegularized();
	      std::cout << "Sub element nr " << kj+1 << ": " << regular << std::endl;
	      if (regular)
		{
		  // Create non-trimmed parameter element
		  int bd_cond[6][2];
		  shared_ptr<ParamVolume> reg_vol = 
		    sub_elem[kj]->getRegParVol(degree, bd_cond);
		  if (reg_vol.get())
		    {
		      reg_vol->writeStandardHeader(of5);
		      reg_vol->write(of5);
		    }

		  // Create non-trimmed element
		  sub_elem[kj]->untrimRegular(degree);
		  shared_ptr<ParamVolume> tmp_vol = sub_elem[kj]->getVolume();
		  tmp_vol->writeStandardHeader(of6);
		  tmp_vol->write(of6);
		}
	      else
		{
		  std::cout << "Number of surfaces: " << sub_elem[kj]->getOuterShell()->nmbEntities() << std::endl;
		}
	    }
	}
    }
}
