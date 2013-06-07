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
