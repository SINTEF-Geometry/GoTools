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

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <memory> // For std::unique_ptr
#include <iomanip> // For std::setprecision

// GoTools headers
#include "GoTools/lrsplines3D/LRSplineVolume.h" 
#include "GoTools/trivariate/SplineVolume.h"  // For defining the initial B-spline volume
#include "GoTools/lrsplines3D/Mesh3D.h"        // The LR mesh
#include "GoTools/lrsplines3D/Direction3D.h" // For XDIR, YDIR, ZDIR enum 
#include "GoTools/geometry/ObjectHeader.h"   // For reading the input spline volume
#include "GoTools/geometry/PointCloud.h"     // For visualization purposes
#include "GoTools/geometry/LineCloud.h"      // For visualization purposes


// Use Go namespace for GoTools types
using namespace Go;
using std::vector;

//===========================================================================
//                                                                           
/// Description:
/// The program reads a tensor-product spline volume and represents this 
/// volume as an LR B-spline volume. This volume is refined in all three
/// parameter directions in the sequence:  
///  
/// Input to the geometry construction is read from data/tpvol.g2
/// Current volumes  are written to g2-files as we go
/// along. The final volume is written to the file
/// data/lrvol_fin.g2. The outer boundary surfaces of the volumes can
/// be visualized in gotools/viewlib/app/goview_vol_and_lr.
/// Also the element midpoints and boundary curves can be visualized in
/// goview_vol_and_lr
//                                                                           
//===========================================================================

// Extract information from a current LR spline volume in order to
// visualize the structure of the parameter domain
void extractElementMidAndBoundary(shared_ptr<LRSplineVolume> vol,
				  vector<double>& mid, vector<double>& bd)
{
  int num_el = vol->numElements();
  mid.resize(3*num_el);   // Mid parameters of each element (xmid, ymid, zmid)
  bd.resize(72*num_el);   // The element corners orginized as corner curves
  
  int ki=0, kj=0;
  Point pos, pos2, dir1, dir2, dir3;
  for (auto it=vol->elementsBegin(); it != vol->elementsEnd(); ++it)
    {
      const Element3D* elem = it->second.get();
      double umin = elem->umin();
      double umax = elem->umax();
      double vmin = elem->vmin();
      double vmax = elem->vmax();
      double wmin = elem->wmin();
      double wmax = elem->wmax();
      mid[ki++] = 0.5*(umin+umax);
      mid[ki++] = 0.5*(vmin+vmax);
      mid[ki++] = 0.5*(wmin+wmax);

      pos = Point(umin, vmin, wmin);
      dir1 = Point(umax-umin, 0.0, 0.0);
      dir2 = Point(0.0, vmax-vmin, 0.0);
      dir3 = Point(0.0, 0.0, wmax-wmin);

      pos2 = pos+dir1;
      for (int ka=0; ka<3; ++ka)
	bd[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd[kj++] = pos2[ka];

      pos2 = pos+dir2;
      for (int ka=0; ka<3; ++ka)
	bd[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd[kj++] = pos2[ka];

      pos2 = pos+dir3;
      for (int ka=0; ka<3; ++ka)
	bd[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd[kj++] = pos2[ka];

      pos += dir1;
      pos2 = pos+dir2;
      for (int ka=0; ka<3; ++ka)
	bd[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd[kj++] = pos2[ka];

      pos2 = pos+dir3;
      for (int ka=0; ka<3; ++ka)
	bd[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd[kj++] = pos2[ka];

      pos += dir2;
      pos2 = pos-dir1;
      for (int ka=0; ka<3; ++ka)
	bd[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd[kj++] = pos2[ka];

      pos2 = pos+dir3;
      for (int ka=0; ka<3; ++ka)
	bd[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd[kj++] = pos2[ka];

      pos += dir3;
      pos2 = pos-dir1;
      for (int ka=0; ka<3; ++ka)
	bd[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd[kj++] = pos2[ka];

      pos2 = pos-dir2;
      for (int ka=0; ka<3; ++ka)
	bd[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd[kj++] = pos2[ka];

      pos -= dir1;
      pos2 = pos-dir2;
      for (int ka=0; ka<3; ++ka)
	bd[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd[kj++] = pos2[ka];

      pos2 = pos-dir3;
      for (int ka=0; ka<3; ++ka)
	bd[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd[kj++] = pos2[ka];

      pos -= dir2;
      pos2 = pos+dir1;
      for (int ka=0; ka<3; ++ka)
	bd[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd[kj++] = pos2[ka];
     }
}

int main(int argc, char *argv[])
{
  std::cout << "--- refine_lrvol: LR-spline Volume Refinement Example ---" << std::endl;

  // --- Read B-spline volume from file
    
  std::string infile("../../gotools-data/lrsplines3D/examples/data/tpvol.g2"); // Input file name
  std::cout << "\nAttempting to read SplineVolume from file: " << infile << std::endl;

  shared_ptr<SplineVolume> spl_vol(new SplineVolume());

  // Specific file reading and volume creation.
  // In GoTools, you would typically use a factory or a specific reader method:
  std::ifstream input_stream(infile.c_str());
  if (!input_stream.is_open()) {
    std::cerr << "Error: Could not open file " << infile << std::endl;
    return 1;
  }
    
  // Read header specifying the type of geometry entity
  // The function throws if the entity header is invalid
  ObjectHeader header;
  try {
    header.read(input_stream);
  }
  catch (...)
    {
      std::cerr << "Exiting" << std::endl;
      exit(-1);
    }
  
  try {
    spl_vol->read(input_stream); // Assumed static read method
  } catch (const std::exception& e) {
    std::cerr << "Error reading SplineVolume: " << e.what() << std::endl;
    return 1;
  }

  // Prepare for output files
  std::string outfile1("data/lrvol_in.g2");
  std::string outfile1_2("data/lrvol_in_mid_elem.g2");
  std::string outfile1_3("data/lrvol_in_bd_elem.g2");
  std::string outfile2("data/lrvol_ref1.g2");
  std::string outfile2_2("data/lrvol_ref1_mid_elem.g2");
  std::string outfile2_3("data/lrvol_ref1_bd_elem.g2");
  std::string outfile3("data/lrvol_ref2.g2");
  std::string outfile3_2("data/lrvol_ref2_mid_elem.g2");
  std::string outfile3_3("data/lrvol_ref2_bd_elem.g2");
  std::string outfile4("data/lrvol_fin.g2");
  std::string outfile4_2("data/lrvol_fin_mid_elem.g2");
  std::string outfile4_3("data/lrvol_fin_bd_elem.g2");


  // Represent the initial spline volume as an LR spline volume
  // The second parameter indicates when two consequtive knots are
  // regarded as identical
  std::cout << "\n--- Converting SplineVolume to LRSplineVolume  ---" << std::endl;
  shared_ptr<LRSplineVolume> lr_vol(new LRSplineVolume(spl_vol.get(),
							       1.0e-6));

  std::ofstream of1(outfile1.c_str());
  lr_vol->writeStandardHeader(of1);
  lr_vol->write(of1);

  // Fetch mid point and boundary curves of elements for visualization of
  // the parameter domain
  vector<double> mid1, bd1;
  extractElementMidAndBoundary(lr_vol, mid1, bd1);
  
  // Make point cloud from mid parameters
  PointCloud3D elem_cloud1(mid1.begin(), mid1.size()/3);

  std::ofstream of1_2(outfile1_2.c_str());
  elem_cloud1.writeStandardHeader(of1_2);
  elem_cloud1.write(of1_2);

  // Make line cloud from element boundaries
  std::ofstream of1_3(outfile1_3.c_str());
  LineCloud elem_lines1(bd1.begin(), bd1.size()/6);
  elem_lines1.writeStandardHeader(of1_3);
  elem_lines1.write(of1_3);

  // Fetch parameter domain of volume
  const Array<double,6> dom = lr_vol->parameterSpan();
  // The sequence is: umin, umax, vmin, vmax, wmin, wmax

  // Fetch LR mesh
  const Mesh3D mesh = lr_vol->mesh();

  // Fetch information about knot vectors in the three parameter directions.
  // Knot multiplicity (including multiplicity 0, i.e. not an active knot)
  // is not included
  // Note that the current volume is really a tensor product volume so all
  // knots traverse the entire parameter domain.
  vector<double> knots1 = mesh.allKnots(XDIR); // 1. parameter direction
  vector<double> knots2 = mesh.allKnots(YDIR); // 2. parameter direction
  vector<double> knots3 = mesh.allKnots(ZDIR); // 3. parameter direction
  std::cout << "Knots in 1. parameter direction: ";
  for (size_t ki=0; ki<knots1.size(); ++ki)
    std::cout << knots1[ki] << " ";
  std::cout << std::endl;
  std::cout << "Knots in 2. parameter direction: ";
  for (size_t ki=0; ki<knots2.size(); ++ki)
    std::cout << knots2[ki] << " ";
  std::cout << std::endl;
  std::cout << "Knots in 3. parameter direction: ";
  for (size_t ki=0; ki<knots3.size(); ++ki)
    std::cout << knots3[ki] << " ";
  std::cout << std::endl;
  
  bool absolute_refinement = true; // true = Do not increment multiplicity 
                                   // during multiple knot insertions,
                                   // false = increment multiplicity
  
  // Refine in second parameter direction
  Direction3D dir_y = YDIR;  // The constant direction of the new knot line segment
  double par_y = 0.5;      // Position of new knot in the corresponding knot vector
  double start_x = 0.0;    // The start value of the new segment in 1. parameter direction
  double end_x = 1.0;      // The end value of the new segment in 1. parameter direction
  // The new knot line segment must cover at least one knot interval and must
  // traverse the support of at least one B-splines. The first condition is met as
  // it starts and ends at distinct knots. Due to multiple knots in the start
  // and end of the parameter domain of the volume, the segment traverses one
  // B-spline
  double start_z = 0.0;    // The start value of the new segment in 3. parameter direction
  double end_z = 1.0;      // The end value of the new segment in 3. parameter direction
  // The new knot line segment traverses the entire parameter domain in the
  // 3. parameter direction
  int mult = 1;            // Knot multiplicity of the new knot is one

  int num_Bsplines = lr_vol->numBasisFunctions();
  int num_element = lr_vol->numElements();
  std::cout << "Number of B-splines prior to refinement: " << num_Bsplines << std::endl;
  std::cout << "Number of elements prior to refinement: " << num_element << std::endl;

  // Permform refinement
  LRSplineVolume::Refinement3D ref_y;
  ref_y.setVal(par_y, start_z, end_z, start_x, end_x, dir_y, mult);
  // Note the sequence of the parameter information. 
  lr_vol->refine(ref_y, absolute_refinement);

  num_Bsplines = lr_vol->numBasisFunctions();
  num_element = lr_vol->numElements();
  std::cout << "Number of B-splines after refinement in 2. parameter directon: " << num_Bsplines << std::endl;
  std::cout << "Number of elements after refinement in 2. parameter directon: " << num_element << std::endl;

  std::ofstream of2(outfile2.c_str());
  lr_vol->writeStandardHeader(of2);
  lr_vol->write(of2);

  // Fetch mid point and boundary curves of elements for visualization of
  // the parameter domain
  vector<double> mid2, bd2;
  extractElementMidAndBoundary(lr_vol, mid2, bd2);
  
  // Make point cloud from mid parameters
  PointCloud3D elem_cloud2(mid2.begin(), mid2.size()/3);

  std::ofstream of2_2(outfile2_2.c_str());
  elem_cloud2.writeStandardHeader(of2_2);
  elem_cloud2.write(of2_2);

  // Make line cloud from element boundaries
  std::ofstream of2_3(outfile2_3.c_str());
  LineCloud elem_lines2(bd2.begin(), bd2.size()/6);
  elem_lines2.writeStandardHeader(of2_3);
  elem_lines2.write(of2_3);

  // Refine in first parameter direction
  Direction3D dir_x = XDIR;  // The constant direction of the new knot line segment
  double par_x1 = 0.25;      // Position of new knot in the corresponding knot vector
  double par_x2 = 0.5;      // Position of new knot in the corresponding knot vector
  double par_x3 = 0.75;      // Position of new knot in the corresponding knot vector
  double start_y = 0.0;    // The start value of the new segment in 2. parameter direction
  double end_y = 0.5;      // The end value of the new segment in 2. parameter direction
  // The new knot line segments end at the previously inserted knot line segment in
  // the second parameter direction
  start_z = 0.0;    // The start value of the new segment in 3. parameter direction
  end_z = 1.0;      // The end value of the new segment in 3. parameter direction
 
  // Permform refinement
  vector<LRSplineVolume::Refinement3D> ref_x(3);
  ref_x[0].setVal(par_x1, start_y, end_y, start_z, end_z, dir_x, mult);
  ref_x[1].setVal(par_x2, start_y, end_y, start_z, end_z, dir_x, mult);
  ref_x[2].setVal(par_x3, start_y, end_y, start_z, end_z, dir_x, mult);
  // Note the sequence of the parameter information. 
  lr_vol->refine(ref_x, absolute_refinement);

  num_Bsplines = lr_vol->numBasisFunctions();
  num_element = lr_vol->numElements();
  std::cout << "Number of B-splines after refinement in 1. parameter directon: " << num_Bsplines << std::endl;
  std::cout << "Number of elements after refinement in 1. parameter directon: " << num_element << std::endl;

  std::ofstream of3(outfile3.c_str());
  lr_vol->writeStandardHeader(of3);
  lr_vol->write(of3);

  // Fetch mid point and boundary curves of elements for visualization of
  // the parameter domain
  vector<double> mid3, bd3;
  extractElementMidAndBoundary(lr_vol, mid3, bd3);
  
  // Make point cloud from mid parameters
  PointCloud3D elem_cloud3(mid3.begin(), mid3.size()/3);

  std::ofstream of3_2(outfile3_2.c_str());
  elem_cloud3.writeStandardHeader(of3_2);
  elem_cloud3.write(of3_2);

  // Make line cloud from element boundaries
  std::ofstream of3_3(outfile3_3.c_str());
  LineCloud elem_lines3(bd3.begin(), bd3.size()/6);
  elem_lines3.writeStandardHeader(of3_3);
  elem_lines3.write(of3_3);

  // Refine in third parameter direction
  Direction3D dir_z = ZDIR;  // The constant direction of the new knot line segment
  double par_z = 0.25;        // Position of new knot in the corresponding knot vector
  start_x = 0.25;    // The start value of the new segment in 1. parameter direction 
  end_x = 1.0;      // The end value of the new segment in 1. parameter direction
  // The start and end of the new segment corresponds to previously inserted
  // knot line segments
  start_y = 0.0;    // The start value of the new segment in 2. parameter direction
  end_y = 0.5;      // The end value of the new segment in 2. parameter direction

 
  // Permform refinement
  LRSplineVolume::Refinement3D ref_z;
  ref_z.setVal(par_z, start_x, end_x, start_y, end_y, dir_z, mult);
  // Note the sequence of the parameter information. 
  lr_vol->refine(ref_z, absolute_refinement);

  num_Bsplines = lr_vol->numBasisFunctions();
  num_element = lr_vol->numElements();
  std::cout << "Number of B-splines after refinement in 3. parameter directon: " << num_Bsplines << std::endl;
  std::cout << "Number of elements after refinement in 3. parameter directon: " << num_element << std::endl;

  std::ofstream of4(outfile4.c_str());
  lr_vol->writeStandardHeader(of4);
  lr_vol->write(of4);

  // Fetch mid point and boundary curves of elements for visualization of
  // the parameter domain
  vector<double> mid4, bd4;
  extractElementMidAndBoundary(lr_vol, mid4, bd4);
  
  // Make point cloud from mid parameters
  PointCloud3D elem_cloud4(mid4.begin(), mid4.size()/3);

  std::ofstream of4_2(outfile4_2.c_str());
  elem_cloud4.writeStandardHeader(of4_2);
  elem_cloud4.write(of4_2);

  // Make line cloud from element boundaries
  std::ofstream of4_3(outfile4_3.c_str());
  LineCloud elem_lines4(bd4.begin(), bd4.size()/6);
  elem_lines4.writeStandardHeader(of4_3);
  elem_lines4.write(of4_3);

  std::cout << "Program finished" << std::endl;
  return 0;
}

