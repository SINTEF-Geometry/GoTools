//===========================================================================
//
// File: testPeeling.cpp
//
// Created: Wed Oct 31 12:00:00 2012
//
// Author: Peter Nortoft <penn@win.dtu.dk>
//
// Revision: $Id: $
//
// Description: A simple test of the peeling methods in LinDepUtils.
//
//===========================================================================

#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LinDepUtils.h"
#include <iostream>
#include <fstream>
#include <stdexcept>

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
   const vector<double> knots_u = {0, 0, 0, 1, 2, 3, 6, 8, 9, 10, 10, 10};
   const vector<double> knots_v = {0, 0, 0, 1, 2, 4, 6, 7, 8, 8, 8};
   const vector<double> coefs(coefs_u * coefs_v * dimension, 0);
   const vector<LRSplineSurface::Refinement2D> refs = {
           {5,2,7,XFIXED,1},
           {7,2,6,XFIXED,1},
           {5,1,7,YFIXED,1},
           {3,3,9,YFIXED,1},
           {4,1,5,XFIXED,1}
         };
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
   const vector<double> knots_u = {0, 0, 0, 4, 8, 12, 16, 20, 20, 20};
   const vector<double> knots_v = {0, 0, 4, 8, 12, 16, 16};
   const vector<double> coefs(coefs_u * coefs_v * dimension, 0);
   const vector< vector<LRSplineSurface::Refinement2D> > refs = { 
         {
           { 6,4,12,XFIXED,2},
           {10,4,12,XFIXED,2},
           {14,4,12,XFIXED,2},
           { 6,4,16,YFIXED,1},
           {10,4,16,YFIXED,1}
         }
         ,
         {
           { 7,6,10,XFIXED,2},
           { 9,6,10,XFIXED,2},
           {11,6,10,XFIXED,2},
           {13,6,10,XFIXED,2},
           { 7,6,14,YFIXED,1},
           { 9,6,14,YFIXED,1}
         }
   };
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
