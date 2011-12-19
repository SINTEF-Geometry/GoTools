//===========================================================================
//
// File : findCommonFaces.C
//
// Created: Fri Jan  2 10:39:13 2009
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: findCommonFaces.C,v 1.3 2009-06-12 08:45:17 vsk Exp $
//
// Description: Read a family of volumes, and store pairs of boundary surfaces that are the same
//
//===========================================================================

/*
  The boundary surfaces that coincide are stored as numbers in  a file, one line for each
  relation, on the following format:
  <vol1> <surf1> <vol2> <surf2> <swapParDir?> <reverse_u> <reverse_v>
  where
  - <vol1> is the index of the volume of the first coinciding boundary surface. The index
    is given by counting the volumes as they are stored in the input file, starting at 0.
  - <surf1> is the boundary surface index on <vol1> for the first coinciding boundary surface,
    a number from 0 to 5, corresponding to u_min, u_max, v_min, v_max, w_min, w_max
  - <vol2> is the index of the volume of the second coinciding boundary surface
  - <surf2> is the boundary surface index on <vol2> for the second coinciding boundary surface
  - <swapParDir?> is 0 or 1.
    - 0 if the u- and v-parameter directions on <surf1> coincide with u and v respectively on <surf2>
    - 1 if the u- and v-parameter directions on <surf1> coincide with v and u respectively on <surf2>
  - <reverse_u> is 0 or 1
    - 0 if the start- and end-points in the u-parameter direction on <surf1> coincide with start and end
      respectively on <surf2> (either along its u- or v-parameter direction, depending on <swapParDir?>
    - 1 if the start- and end-points in the u-parameter direction on <surf1> coincide with end and start
      respectively on <surf2> (either along its u- or v-parameter direction, depending on <swapParDir?>
  - <reverse_v> is 0 or 1, describing start and end in v-parameter direction on <surf_2>, see <reverse_u>
*/


#include <fstream>
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/errormacros.h"
#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/geometry/Utils.h"

using namespace Go;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::vector;

int main(int argc, char* argv[] )
{
  if (argc != 3)
    {
      cout << "Usage: " << argv[0] << " volumesinfile relationsoutfile" << endl;
      exit(-1);
    }

  // Open input volumes file
  ifstream is(argv[1]);
  ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

  // Open outfile
  ofstream os(argv[2]);
  ALWAYS_ERROR_IF(os.bad(), "Bad or no output filename");


  // Read volumes from files
  vector<vector<shared_ptr<SplineSurface> > > faces;
  while (!is.eof())
    {

      // Read volume from file
      ObjectHeader head;
      SplineVolume vol;
      is >> head >> vol;

      
      faces.push_back(vol.getBoundarySurfaces());

      eatwhite(is);
    }

  // Run through all choices for first face
  double tol = 1.0e-4;
  for (size_t i = 0; i < faces.size()-1 ; ++i)
    for (size_t j = 0; j < faces[i].size(); ++j)
      {

	shared_ptr<SplineSurface> surf1 = faces[i][j];
	double start1_u = surf1 -> startparam_u();
	double end1_u = surf1 -> endparam_u();
	double start1_v = surf1 -> startparam_v();
	double end1_v = surf1 -> endparam_v();
	double range1_u = end1_u - start1_u;
	double range1_v = end1_v - start1_v;
	int order1_u = surf1 -> order_u();
	int order1_v = surf1 -> order_v();
	int numCoefs1_u = surf1 -> numCoefs_u();
	int numCoefs1_v = surf1 -> numCoefs_v();
	vector<double>::const_iterator knots1_u = surf1 -> basis_u().begin();
	vector<double>::const_iterator knots1_v = surf1 -> basis_v().begin();
	int kdim = surf1->dimension();
	vector<double>::const_iterator coefs1 = surf1 -> ctrl_begin();
	if (surf1->rational())
	  ++kdim;

	// Run through all choices for second face
	for (size_t k = i+1; k < faces.size() ; ++k)
	  for (size_t l = 0; l < faces[k].size(); ++l)
	    {

	      shared_ptr<SplineSurface> surf2 = faces[k][l];

	      // Test dimesnsion of space of surface and rationality
	      if (surf1->dimension() != surf2->dimension() || surf1->rational() != surf2->rational()) continue;

	      // We now run through the eight ways surf2 can be mirrored or rotated, and test if it coincides with surf1
	      // First run through two choices of parameter order for second face: Either fixed or u <-> v
	      for (int par_order = 0; par_order < 2; ++par_order)
		{

		  // Test if order and number of knots coincide
		  if (par_order == 0 && (order1_u != surf2 -> order_u() ||
					 order1_v != surf2 -> order_v() ||
					 numCoefs1_u != surf2 -> numCoefs_u() ||
					 numCoefs1_v != surf2 -> numCoefs_v()))
		    continue;
		  if (par_order == 1 && (order1_u != surf2 -> order_v() ||
				     order1_v != surf2 -> order_u() ||
				     numCoefs1_u != surf2 -> numCoefs_v() ||
				     numCoefs1_v != surf2 -> numCoefs_u()))
		    continue;

		  // Swap u <-> v if necessary
		  if (par_order == 1) surf2 -> swapParameterDirection();

		  vector<double>::const_iterator knots2_u = surf2 -> basis_u().begin();
		  vector<double>::const_iterator knots2_v = surf2 -> basis_v().begin();
		  vector<double>::const_iterator coefs2 = surf2 -> ctrl_begin();

		  // Run through the two ways to run along first parameter direction: From start to end or end to start
		  for (int reverse_u = 0; reverse_u < 2; ++reverse_u)
		    {
		      double start2_u;
		      double end2_u;
		      if (reverse_u == 0)
			{
			  start2_u = surf2 -> startparam_u();
			  end2_u = surf2 -> endparam_u();
			}
		      else
			{
			  start2_u = surf2 -> endparam_u();
			  end2_u = surf2 -> startparam_u();
			}
		      double range2_u = end2_u - start2_u;

		      // Test if knot vectors in first direction are equal up to linear transformation
		      int m;
		      for (m = 0; m < order1_u + numCoefs1_u; ++m)
			if ((reverse_u == 0 && fabs((knots1_u[m] - start1_u) * range2_u - (knots2_u[m] - start2_u) * range1_u) > tol) ||
			    (reverse_u == 1 && fabs((knots1_u[m] - start1_u) * range2_u - (knots2_u[order1_u+numCoefs1_u-1-m] - start2_u) * range1_u) > tol))
			  break;
		      if (m < order1_u + numCoefs1_u) continue;

		      // Run through the two ways to run along second parameter direction: From start to end or end to start
		      int reverse_v;
		      for (reverse_v = 0; reverse_v < 2; ++reverse_v)
			{
			  double start2_v;
			  double end2_v;
			  if (reverse_v == 0)
			    {
			      start2_v = surf2 -> startparam_v();
			      end2_v = surf2 -> endparam_v();
			    }
			  else
			    {
			      start2_v = surf2 -> endparam_v();
			      end2_v = surf2 -> startparam_v();
			    }
			  double range2_v = end2_v - start2_v;

			  // Test if knot vectors in second direction are equal up to linear transformation
			  if (reverse_u == 0)
			    {
			      for (m = 0; m < order1_v + numCoefs1_v; ++m)
				if ((reverse_v == 0 && fabs((knots1_v[m] - start1_v) * range2_v - (knots2_v[m] - start2_v) * range1_v) > tol) ||
				    (reverse_v == 1 && fabs((knots1_v[m] - start1_v) * range2_v - (knots2_v[order1_v+numCoefs1_v-1-m] - start2_v) * range1_v) > tol))
				  break;
			      if (m < order1_v + numCoefs1_v) 
				continue;
			    }

			  // At this point, the Bspline basis of surf1 and surf2 (as it is rotated) are equal. Remains to check coefficients
			  int c_start = 0;
			  int step_u = kdim;
			  int step_v = numCoefs1_u*kdim;
			  if (reverse_u == 1)
			    {
			      step_u *= -1;
			      c_start += kdim*(numCoefs1_u-1);
			    }
			  if (reverse_v == 1)
			    {
			      step_v *= -1;
			      c_start += kdim*numCoefs1_u*(numCoefs1_v-1);
			    }

			  if (surf1->rational())
			    {
			      int n, d;
			      for (m = 0; m < numCoefs1_u; ++m)
				{
				  for (n = 0; n < numCoefs1_v; ++n)
				    {
				      for (d = 0; d < kdim-1; ++d)
					if (fabs(coefs1[n*numCoefs1_u*kdim+m*kdim+d]*coefs2[c_start+n*step_v+m*step_u+kdim-1]-coefs2[c_start+n*step_v+m*step_u+d]*coefs1[n*numCoefs1_u*kdim+m*kdim+kdim-1]) > tol)
					  break;
				      if (d < kdim-1) break;
				    }
				  if (n < numCoefs1_v)
				    break;
				}
			      if (m < numCoefs1_u)
				continue;
			    }
			  else
			    {
			      int n, d;
			      for (m = 0; m < numCoefs1_u; ++m)
				{
				  for (n = 0; n < numCoefs1_v; ++n)
				    {
				      for (d = 0; d < kdim; ++d)
					if (fabs(coefs1[n*numCoefs1_u*kdim+m*kdim+d]-coefs2[c_start+n*step_v+m*step_u+d]) > tol)
					  break;
				      if (d < kdim) 
					break;
				    }
				  if (n < numCoefs1_v)
				    break;
				}
			      if (m < numCoefs1_u)
				continue;
			    }

			  // surf1 and surf2 coincide. Store result.
			  os << i << " " << j << " " << k << " " << l << " " << par_order << " " << reverse_u << " " << reverse_v << endl;
			  break;  // Match found

			}  //  for (int reverse_v = 0; reverse_v < 2; ++reverse_v)

		      if (reverse_v < 2)
			break;  // Match already found
		    }  //  for (int reverse_u = 0; reverse_u < 2; ++reverse_u)

		  // Swap u <-> v back, to leave surf2 unchanged
		  if (par_order == 1) surf2 -> swapParameterDirection();

		}  //  for (int par_order = 0; par_order < 2; ++par_order)

	    }  //  for (size_t k = i+1; k < faces.size() ; ++k)     for (size_t l = 0; l < faces[k].size(); ++l)

      }  //  for (size_t i = 0; i < faces.size()-1 ; ++i)    for (size_t j = 0; j < faces[i].size(); ++j)

}
