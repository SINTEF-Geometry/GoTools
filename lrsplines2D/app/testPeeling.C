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

#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LinDepUtils.h"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>

using namespace std;
using namespace Go;

namespace { // begin anonymous namespace

  // -----------------------------------------------------------------
  // Construct the LR spline example from the paper that results in
  // linear dependence.
  LRSplineSurface construct_lrspline_from_paper( ) {
   const int deg_u     = 2;
   const int deg_v     = 2;
   const int coefs_u   = 9;
   const int coefs_v   = 8;
   const int dimension = 1;
   const double k_u[12] =  {0, 0, 0, 1, 2, 3, 6, 8, 9, 10, 10, 10};
   const double k_v[11] =  {0, 0, 0, 1, 2, 4, 6, 7, 8, 8, 8};
   const vector<double> knots_u(k_u,k_u+12);
   const vector<double> knots_v(k_v,k_v+11);
   const vector<double> coefs(coefs_u * coefs_v * dimension, 0);
   const LRSplineSurface::Refinement2D rfs[5] = {
	   {5,2,7,XFIXED,1},
	   {7,2,6,XFIXED,1},
	   {5,1,7,YFIXED,1},
	   {3,3,9,YFIXED,1},
	   {4,1,5,XFIXED,1}
   };
   const vector<LRSplineSurface::Refinement2D> refs(rfs,rfs+5);
   LRSplineSurface lrs(deg_u, deg_v, coefs_u, coefs_v, dimension, knots_u.begin(), knots_v.begin(), coefs.begin());
   lrs.refine(refs);
   return lrs;
  }

  // -----------------------------------------------------------------
  // Construct a "hierarchical" LR spline.
  LRSplineSurface construct_hierarchical_lrspline( ) {
   const int deg_u     = 2;
   const int deg_v     = 1;
   const int coefs_u   = 7;
   const int coefs_v   = 5;
   const int dimension = 1;
   const double k_u[10] = {0, 0, 0, 4, 8, 12, 16, 20, 20, 20};
   const double k_v[7] = {0, 0, 4, 8, 12, 16, 16};
   const vector<double> knots_u(k_u,k_u+10);
   const vector<double> knots_v(k_v,k_v+7);
   const vector<double> coefs(coefs_u * coefs_v * dimension, 0);

   const LRSplineSurface::Refinement2D rf_1[5] = { 
	   { 6,4,12,XFIXED,2},
	   {10,4,12,XFIXED,2},
	   {14,4,12,XFIXED,2},
	   { 6,4,16,YFIXED,1},
	   {10,4,16,YFIXED,1}
   };
   const LRSplineSurface::Refinement2D rf_2[6] = { 
	   { 7,6,10,XFIXED,2},
	   { 9,6,10,XFIXED,2},
	   {11,6,10,XFIXED,2},
	   {13,6,10,XFIXED,2},
	   { 7,6,14,YFIXED,1},
	   { 9,6,14,YFIXED,1}
   };
   vector< vector<LRSplineSurface::Refinement2D> > refs;
   refs.push_back(vector<LRSplineSurface::Refinement2D>(rf_1,rf_1+5));
   refs.push_back(vector<LRSplineSurface::Refinement2D>(rf_2,rf_2+6));
   
   LRSplineSurface lrs(deg_u, deg_v, coefs_u, coefs_v, dimension, knots_u.begin(), knots_v.begin(), coefs.begin());
   for (auto it=refs.begin(); it!=refs.end(); ++it)
     lrs.refine(*it);
   return lrs;
  }

  // -----------------------------------------------------------------
  // Reads an LR spline from a file.
  LRSplineSurface read_lrspline_from_file( std::string filename ) {
    cout << "Reading LR spline from file: " << filename << " ..." << endl;
    LRSplineSurface lrs;
    ifstream infile (filename);
    if (!infile)
      throw runtime_error("Could not open file: " + filename + ".");
    lrs.read(infile);
    if (!infile)
      throw runtime_error("Could not read file: " + filename + ".");
    infile.close();
    if (!infile)
      throw runtime_error("Could not close file: " + filename + ".");
    cout << " ... Done!" << endl;
    return lrs;
  }

} // end anonymous namespace

// -----------------------------------------------------------------
// Constructs LR spline, tests its peelability, and exits.
// -----------------------------------------------------------------
int main () {

  // Flush output buffer.
  cout << endl;

  // Construct LR spline.
  LRSplineSurface lrs = construct_lrspline_from_paper();
  //LRSplineSurface lrs = construct_hierarchical_lrspline();
  
  // Print sizes
  cout << "numElements: " << lrs.numElements() << endl;
  cout << "numBasisFunctions: " << lrs.numBasisFunctions() << endl;

  // Check for overloading.
  vector<LRBSpline2D*> funs = LinDepUtils::unpeelableBasisFunctions(lrs);
  cout << "Number of unpeelable LR B-splines: " << funs.size() << endl;
  cout << "LR spline is peelable: " << ((funs.size()==0) ? "Yes" : "No") << "!" << endl;
  //bool lrs_is_ovl = isPeelable(lrs);
  //cout << "LR spline is peelable Y/N: " << lrs_is_ovl << endl;
  
  std::ofstream of("non_peelable_Bsplines.g2");
  lrs.writeStandardHeader(of);
  lrs.write(of);
  of << "===========================================================" << std::endl << std::endl << std::endl;
  for (size_t ki=0; ki<funs.size(); ++ki)
    {
      funs[ki]->write(of);
      of << std::endl;
    }

  // That's it.
  return 0;

}
