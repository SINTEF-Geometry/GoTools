//===========================================================================
//
// File : hemisphereShell.C
//
// Created: Fri Dec 19 11:03:51 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: hemisphereShell.C,v 1.1 2009-01-09 11:02:10 kfp Exp $
//
// Description: Program to make a hemisphere with stiffener
//
//===========================================================================


/*

  For description of solid, see page 4168 in 'Isogeometric analysis: CAD,
  finite elements, NURBS, exact geometry and mesh refinement',
  Comput. Methods Appl. Mech. Engrg. 194 (2005)

  The program takes three arguments.
  - arg. 1 is the input filename for data describing the geometry.
  - arg. 2 is the output filename for the volumes describing the hemisphere and stiffeners
  - arg. 3 is the output filename for the coinciding boundary surfaces of the volumes

  About the file formats:

  Arg 1:
  Data input file format er numbers in following order
  - 0 if angles are radians, 1 if degrees
  - Angle between baseline and line from center of shell to intersection of stiffener and inner side of shell.
  - Angle between vertical rotation axis and line from shell to top of solid
  - Inner radius of shell
  - Shell thickness
  - Stiffener width
  - Stiffener height

  Arg 2:
  The volumes are stored as four solids (see below) in g2-format

  Arg 3 - see also the program findCommonFaces.C
  The boundary surfaces that coincide are stored as numbers in  a file, one line for each
  relation, on the following format:
  <vol1> <surf1> <vol2> <surf2> <swapParDir?> <reverse_u> <reverse_v>
  where
  - <vol1> is the index of the volume of the first coinciding boundary surface. The index
    is 0, 1, 2 or 3, for hemisphere (0) and then inner (1), center (2) and outer (3) part
    of the stiffener (see below for the four solid parts)
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


  The program produces four solids, the hemisphere and three parts of the stiffener separated by vertical cylindrical faces:
  Solid 0: The hemisphere (without the stiffener)
  Solid 1: The inner stiffener part, where the top face lies inside the hemisphere.
  Solid 2: The center stiffener part, where the top face is the same as the bottom face of the hemisphere
  Solid 3: The outer stiffener part, where the top face lies outside the hemisphere.



 */




#include <iostream>
#include <fstream>
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SweepSurfaceCreator.h"
#include "GoTools/trivariate/SweepVolumeCreator.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/utils/Point.h"

using std::endl;
using std::ifstream;
using std::ofstream;
using std::vector;
using namespace Go;




int main(int argc, char* argv[] )
{

  ALWAYS_ERROR_IF(argc != 4, "Usage: " << argv[0]
		  << " datainfile volumeoutfile endfacerelationsoutfile" << endl);

  // Open input dat a file
  ifstream is(argv[1]);
  ALWAYS_ERROR_IF(is.bad(), "Bad or no data input filename");

  // Open output volumes file
  ofstream osvols(argv[2]);
  ALWAYS_ERROR_IF(osvols.bad(), "Bad volumes output filename");

  // Open output relations file
  ofstream osrels(argv[3]);
  ALWAYS_ERROR_IF(osrels.bad(), "Bad face relations output filename");

  // Read data
  int is_degrees;
  double bot_angle_in;
  double top_angle;
  double shell_rad_in;
  double shell_thick;
  double stiff_width;
  double stiff_height;
  is >> is_degrees >> bot_angle_in >> top_angle >> shell_rad_in >> shell_thick >> stiff_width >> stiff_height;

  if (is_degrees)
    {
      bot_angle_in *= M_PI / (double)(180.0);
      top_angle *= M_PI / (double)(180.0);
    }

  // Some other data
  double shell_rad_out = shell_rad_in + shell_thick;
  double bot_angle_out = asin( shell_rad_in * sin(bot_angle_in) / shell_rad_out);

  vector<double> quart_circle_knots(6);
  quart_circle_knots[0] = quart_circle_knots[1] = quart_circle_knots[2] = 0.0;
  quart_circle_knots[3] = quart_circle_knots[4] = quart_circle_knots[5] = M_PI * (double)0.5;

  vector<double> quart_circle_rcoefs(12);
  quart_circle_rcoefs[1] = quart_circle_rcoefs[2] = quart_circle_rcoefs[5]
    = quart_circle_rcoefs[8] = quart_circle_rcoefs[9] = 0.0;
  quart_circle_rcoefs[0] = quart_circle_rcoefs[3]
    = quart_circle_rcoefs[10] = quart_circle_rcoefs[11] = 1.0;
  quart_circle_rcoefs[4] = quart_circle_rcoefs[6] = quart_circle_rcoefs[7] = sqrt((double)0.5);

  SplineCurve quart_circle(3, 3, quart_circle_knots.begin(), quart_circle_rcoefs.begin(), 3, true);

  // Calculate start and end parameters for circle segments
  double diag_bot_in = bot_angle_in - M_PI * double(0.25);
  double diag_bot_out = bot_angle_out - M_PI * double(0.25);
  double diag_top = M_PI * double(0.25) - top_angle;

  double param_bot_in = M_PI * (double)(0.25) * (1.0 
						 + (sqrt ((double) 2.0) + 1.0)
						 * tan (diag_bot_in / 2.0));
  double param_bot_out = M_PI * (double)(0.25) * (1.0 
						  + (sqrt ((double) 2.0) + 1.0)
						  * tan (diag_bot_out / 2.0));
  double param_top = M_PI * (double)(0.25) * (1.0 
					      + (sqrt ((double) 2.0) + 1.0)
					      * tan (diag_top / 2.0));

  shared_ptr<SplineCurve> outer_circle(quart_circle.subCurve(param_bot_out, param_top));
  shared_ptr<SplineCurve> inner_circle(quart_circle.subCurve(param_bot_in, param_top));

  vector<double>::iterator rcoefs_circle = inner_circle->rcoefs_begin();
  for (int i = 0; i < inner_circle->numCoefs(); ++i)
    {
      rcoefs_circle[i<<2 | 0] *= shell_rad_in;
      rcoefs_circle[i<<2 | 2] *= shell_rad_in;
    }

  rcoefs_circle = outer_circle->rcoefs_begin();
  for (int i = 0; i < outer_circle->numCoefs(); ++i)
    {
      rcoefs_circle[i<<2 | 0] *= shell_rad_out;
      rcoefs_circle[i<<2 | 2] *= shell_rad_out;
    }

  vector<shared_ptr<SplineCurve> > shell_curves;
  shell_curves.push_back(inner_circle);
  shell_curves.push_back(outer_circle);

  vector<double> loftParams(2);
  loftParams[0] = 0.0;
  loftParams[1] = 1.0;
  outer_circle->basis().rescale(inner_circle->startparam(), inner_circle->endparam());

  SplineSurface* shell_section = CoonsPatchGen::loftSurface(shell_curves.begin(), loftParams.begin(), 2);

  // Description of stiffener
  double stiff_top_z = shell_rad_in * sin(bot_angle_in);
  double stiff_pos2 = shell_rad_in * cos(bot_angle_in);
  double stiff_pos3 = shell_rad_out * cos(bot_angle_out);
  double stiff_smallwidth = (stiff_width - stiff_pos3 + stiff_pos2) * 0.5;
  double stiff_pos1 = stiff_pos2 - stiff_smallwidth;
  double stiff_pos4 = stiff_pos3 + stiff_smallwidth;

  Point stiff_top1(stiff_pos1, 0.0, stiff_top_z);
  Point stiff_top2(stiff_pos2, 0.0, stiff_top_z);
  Point stiff_top3(stiff_pos3, 0.0, stiff_top_z);
  Point stiff_top4(stiff_pos4, 0.0, stiff_top_z);
  Point stiff_bot(stiff_pos1, 0.0, stiff_top_z - stiff_height);

  SplineCurve* stiff_hline_1 = new SplineCurve(stiff_top1, stiff_top2);
  SplineCurve* stiff_hline_2 = new SplineCurve(stiff_top2, stiff_top3);
  SplineCurve* stiff_hline_3 = new SplineCurve(stiff_top3, stiff_top4);
  SplineCurve* stiff_vline = new SplineCurve(stiff_top1, stiff_bot);

  SplineSurface* stiff_section1
    = SweepSurfaceCreator::linearSweptSurface(*stiff_hline_1,
					      *stiff_vline,
					      stiff_top1);
  SplineSurface* stiff_section2
    = SweepSurfaceCreator::linearSweptSurface(*stiff_hline_2,
					      *stiff_vline,
					      stiff_top1);
  SplineSurface* stiff_section3
    = SweepSurfaceCreator::linearSweptSurface(*stiff_hline_3,
					      *stiff_vline,
					      stiff_top1);


  Point origo(0.0, 0.0, 0.0);
  Point z_axis(0.0, 0.0, 1.0);

  SplineVolume* shell
    = SweepVolumeCreator::rotationalSweptVolume(*shell_section,
						2 * M_PI,
						origo,
						z_axis);
  SplineVolume* stiff1
    = SweepVolumeCreator::rotationalSweptVolume(*stiff_section1,
						2 * M_PI,
						origo,
						z_axis);
  SplineVolume* stiff2
    = SweepVolumeCreator::rotationalSweptVolume(*stiff_section2,
						2 * M_PI,
						origo,
						z_axis);
  SplineVolume* stiff3
    = SweepVolumeCreator::rotationalSweptVolume(*stiff_section3,
						2 * M_PI,
						origo,
						z_axis);

  shell->writeStandardHeader(osvols);
  shell->write(osvols);
  stiff1->writeStandardHeader(osvols);
  stiff1->write(osvols);
  stiff2->writeStandardHeader(osvols);
  stiff2->write(osvols);
  stiff3->writeStandardHeader(osvols);
  stiff3->write(osvols);

  // Remark: At this moment, the code producing the coinciding faces releations is only found in a separate program findCommonFaces.
  // We therefore store the relations here hardcoded.
  // However, findCommonFaces has been run on the volumes, and did produce the same result

  osrels << "0 2 2 4 0 0 0" << endl;   // Hemisphere and center stiffener part
  osrels << "1 3 2 2 0 0 0" << endl;   // Inner and center stiffener part
  osrels << "2 3 3 2 0 0 0" << endl;   // Center and outer stiffener part

}
