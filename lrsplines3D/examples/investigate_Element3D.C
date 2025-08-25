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
#include <map>
#include <memory> // For std::unique_ptr and shared_ptr
#include <string> // For std::string
#include <fstream> // For file operations

// GoTools headers
#include "GoTools/lrsplines3D/LRSplineVolume.h" //
#include "GoTools/lrsplines3D/Element3D.h"     //
#include "GoTools/lrsplines3D/LRBSpline3D.h"   // Basis functions for Element3D's support
#include "GoTools/utils/Point.h"             // For geometric points, e.g., for 'contains' method
#include "GoTools/lrsplines3D/Direction3D.h" // Used in Element3D::split
#include "GoTools/geometry/ObjectHeader.h"


using namespace Go; // To use Go:: classes directly
using std::vector;

//===========================================================================
//                                                                           
/// Description:
/// Read LR B-spline volume from file.
/// Iterate through all elements in a surface and demonstrate available 
/// enquiries excluding those connected to scattered data approximation.
///  
/// Input to the example is a volume constructed in the example
/// program refine_lrvol.
/// For each element, the Bezier coefs of the corresponding patch are computed
/// and stored in the file data/Bezier_coefs.g2. The corresponding patches
/// are represented as spline surfaces and written to data/Bezier_patches.g2.
//                                                                           
//===========================================================================

int main(int argc, char *argv[]) {
    std::cout << "Example program for LRSplineVolume and Element3D functionality." << std::endl;

    // --- Read LR B-spline volume from file
    
    std::string infile("data/lrvol_ref2.g2"); // Input file name
    std::cout << "\nAttempting to read LRSplineVolume from file: " << infile << std::endl;

    std::unique_ptr<LRSplineVolume> lrSplineVolume(new LRSplineVolume());

    // Specific file reading and volume creation.
    // In GoTools, you can use a factory or a specific reader method:
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
        lrSplineVolume->read(input_stream);
    } catch (const std::exception& e) {
        std::cerr << "Error reading LRSplineVolume: " << e.what() << std::endl;
        return 1;
    }

    // Fetch parameter domain of volume
    const Array<double,6> dom = lrSplineVolume->parameterSpan();
    // The sequence is: umin, umax, vmin, vmax, wmin, wmax
  // Mid point of parameter domain
    double umid = 0.5*(dom[0] + dom[1]);
    double vmid = 0.5*(dom[2] + dom[3]);
    double wmid = 0.5*(dom[4] + dom[5]);
      
    std::cout << "\n--- Iterating through Elements in LRSplineVolume" << std::endl;

    int element_count = 0;
    // Iteration loop:
    if (lrSplineVolume) {
      for (LRSplineVolume::ElementMap::const_iterator el_it = lrSplineVolume->elementsBegin();
	   el_it != lrSplineVolume->elementsEnd(); ++el_it)
	{
	  Go::Element3D* element = el_it->second.get(); // Get the Element3D pointer
	  element_count++;
    
	  std::cout << "\n--- Element " << element_count << " ---" << std::endl;
    
	  // Enquire parameter domain of the current element
	  double umin = element->umin();
	  double umax = element->umax();
	  double vmin = element->vmin();
	  double vmax = element->vmax();
	  double wmin = element->wmin();
	  double wmax = element->wmax();
	  std::cout << "  U-Domain: [" << umin << ", " << umax << "]" << std::endl;
	  std::cout << "  V-Domain: [" << vmin << ", " << vmax << "]" << std::endl;
	  std::cout << "  W-Domain: [" << wmin << ", " << wmax << "]" << std::endl;
    
	  // Volume of element domain
	  double volume = element->volume();
	  std::cout << "  Element Volume: " << volume << std::endl;
    
	  // Number of B-splines having this element in its support
	  int num_Bspline = element->nmbBasisFunctions();
	  std::cout << "  Number of LR B-splines in support: " << num_Bspline << std::endl;
    
	  // Fetch all B-splines having this element in its support
	  const std::vector<LRBSpline3D*>& Bsplines = element->getSupport();
	  if (!Bsplines.empty()) {
	    std::cout << "    (e.g., First LR B-spline pointer in support: " << Bsplines[0] << ")" << std::endl;
	    // Access a specific supporting B-spline
	    LRBSpline3D* curr_bspline = element->supportFunction(0);
	    std::cout << "    (Accessing first support function via supportFunction(0): " << curr_bspline << ")" << std::endl;
	  }
    
	  // Check if a parameter point is contained in the element domain
	  bool is_inside = element->contains(umid, vmid, wmid);
	  std::cout << "  Contains mid-point (" << umid << ", " << vmid << ", " << wmid << ")? "
		    << (is_inside ? "Yes" : "No") << std::endl;
    
	  // Information related to approximation data (from Approx3DData struct)
	  std::cout << "  --- Approximation Data (if applicable) ---" << std::endl;
	  if (element->hasDataPoints()) { //
	    std::cout << "    Has associated data points: Yes" << std::endl;
	    std::cout << "    Number of scattered data points: " << element->nmbDataPoints() << std::endl; //
	    // You could also fetch and iterate through getDataPoints()
	  } else {
	    std::cout << "    Has associated data points: No" << std::endl;
	  }
    
	  bool info = element->hasAccuracyInfo();
	  std::cout << "    Has accuracy info: " << info << std::endl;
	  if (info)
	    {
	      double avg_err, max_err;
	      int nmb_out_tol;
	      element->getAccuracyInfo(avg_err, max_err, nmb_out_tol); //
	      std::cout << "    Has accuracy info: Yes" << std::endl;
	      std::cout << "    Average Error: " << avg_err << std::endl;
	      std::cout << "    Maximum Error: " << max_err << std::endl;
	      std::cout << "    Number of points outside tolerance: " << nmb_out_tol << std::endl;
            } 
    

        }
    } else {
      std::cout << "  LRSplineVolume object is not valid. Cannot iterate." << std::endl;
    }
    std::cout << "*/" << std::endl;

    return 0;
}


