//===========================================================================
//
// File : stiffedPanelShipHull.C
//
// Created: Mon Jan 19 12:17:36 2009
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: stiffedPanelShipHull.C,v 1.2 2009-01-23 12:36:22 kfp Exp $
//
// Description:
//
//===========================================================================


/*

  For description, 'see stiffed panel in a ship hull' in the Expectations.doc document,
  ICADA project.

  The program takes three arguments.
  - arg. 1 is the input filename for data describing the geometry.
  - arg. 2 is the output filename for the volumes describing the volumes

  About the file formats:

  Arg 1:
  Data input file format er numbers in following order (see document)
  - Frame spacing (mm)
  - Stiffener spacing (mm)
  - Plate thickness (mm)
  - Stiffener height (mm)
  - Web thickness (mm)
  - Flange width (mm)
  - Flange thickness (mm)
  - Number of stiffeners (integer)

  Arg 2:
  The volumes are stored as a set of solids in g2-format

 */




#include <iostream>
#include <fstream>
#include <memory>
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/SweepSurfaceCreator.h"
#include "GoTools/trivariate/SweepVolumeCreator.h"
#include "GoTools/geometry/ObjectHeader.h"
// #include "GoTools/geometry/GeometryTools.h"
#include "GoTools/utils/Point.h"

using std::cin;
using std::cout;
using std::vector;
using std::ostream;
using std::ofstream;
using std::ifstream;
using std::endl;
using namespace Go;


vector<shared_ptr<SplineVolume> > split(vector<shared_ptr<SplineVolume> >& vols, vector<double>& splits, int pardir, unsigned int mask4 = 15)
{
  vector<shared_ptr<SplineVolume > > result;
  for (vector<shared_ptr<SplineVolume> >::iterator it = vols.begin(); it != vols.end(); ++it)
    {
      vector<shared_ptr<SplineVolume > > split_vols = (*it)->split(splits, pardir);
      int pos = 0;
      for (vector<shared_ptr<SplineVolume> >::iterator it2 = split_vols.begin(); it2 != split_vols.end(); ++it2, ++pos)
	if (mask4 & 1<<(pos&3))
	  result.push_back(*it2);
    }
  return result;
}

void store(const vector<shared_ptr<SplineVolume> >& vols, ostream& os)
{
  for (vector<shared_ptr<SplineVolume> >::const_iterator it = vols.begin(); it != vols.end(); ++it)
    {
      (*it)->writeStandardHeader(os);
      (*it)->write(os);
    }
}



int main(int argc, char* argv[] )
{

  ALWAYS_ERROR_IF(argc != 3, "Usage: " << argv[0]
		  << " datainfile volumeoutfile" << endl);

  // Open input dat a file
  ifstream is(argv[1]);
  ALWAYS_ERROR_IF(is.bad(), "Bad or no data input filename");

  // Open output volumes file
  ofstream os(argv[2]);
  ALWAYS_ERROR_IF(os.bad(), "Bad volumes output filename");

  // Read data
  double frame_spacing;
  double stiffener_spacing;
  double plate_thickness;
  double stiffener_height;
  double web_thickness;
  double flange_width;
  double flange_thickness;
  int numb_stiffeners;

  is >> frame_spacing >> stiffener_spacing >> plate_thickness >> stiffener_height >> web_thickness >> flange_width >> flange_thickness >> numb_stiffeners;

  double stiffener_length = 3.0 * frame_spacing;
  double frame_length = (double) (numb_stiffeners+1) * stiffener_spacing;

  // Make the bottom plane as one solid
  Point plate_corner_0(0.0, 0.0, 0.0);
  SplineCurve plate_edge_x(plate_corner_0, Point(frame_length, 0.0, 0.0));
  SplineCurve plate_edge_y(plate_corner_0, Point(0.0, stiffener_length, 0.0));
  SplineCurve plate_edge_z(plate_corner_0, Point(0.0, 0.0, plate_thickness));
  SplineSurface* plateBottom
    = SweepSurfaceCreator::linearSweptSurface(plate_edge_x,
					      plate_edge_y,
					      plate_corner_0);
  shared_ptr<SplineVolume> entire_plate_vol(SweepVolumeCreator::linearSweptVolume(*plateBottom,
									      plate_edge_z,
									      plate_corner_0));
  vector<shared_ptr<SplineVolume> > entire_plate;
  entire_plate.push_back(entire_plate_vol);

  // Make webs and flanges pointing in y-direction as entire solids
  Point web_corner_0(-web_thickness/2.0, 0.0, plate_thickness);
  SplineCurve web_edge_x(web_corner_0, web_corner_0 + Point(web_thickness, 0.0, 0.0));
  SplineCurve web_edge_y(web_corner_0, web_corner_0 + Point(0.0, stiffener_length, 0.0));
  SplineCurve web_edge_z(web_corner_0, web_corner_0 + Point(0.0, 0.0, stiffener_height - flange_thickness));
  SplineSurface* web_bottom_ydir
    = SweepSurfaceCreator::linearSweptSurface(web_edge_x,
					      web_edge_y,
					      web_corner_0);
  Point sweep_web_pt = web_corner_0;

  Point flange_corner_0(-flange_width/2.0, 0.0, plate_thickness + stiffener_height - flange_thickness);
  SplineCurve flange_edge_x(flange_corner_0, flange_corner_0 + Point(flange_width, 0.0, 0.0));
  SplineCurve flange_edge_y(flange_corner_0, flange_corner_0 + Point(0.0, stiffener_length, 0.0));
  SplineCurve flange_edge_z(flange_corner_0, flange_corner_0 + Point(0.0, 0.0, flange_thickness));
  SplineSurface* flange_bottom_ydir
    = SweepSurfaceCreator::linearSweptSurface(flange_edge_x,
					      flange_edge_y,
					      flange_corner_0);
  Point sweep_flange_pt = flange_corner_0;

  vector<shared_ptr<SplineVolume> > entire_webs_ydir;
  vector<shared_ptr<SplineVolume> > entire_flanges_ydir;
  Point stiffener_step(stiffener_spacing, 0.0, 0.0);

  for (int i = 0; i < numb_stiffeners; ++i)
    {
      sweep_web_pt -= stiffener_step;
      sweep_flange_pt -= stiffener_step;

      shared_ptr<SplineVolume> current_web(SweepVolumeCreator::linearSweptVolume(*web_bottom_ydir,
										 web_edge_z,
										 sweep_web_pt));
      entire_webs_ydir.push_back(current_web);
      shared_ptr<SplineVolume> current_flange(SweepVolumeCreator::linearSweptVolume(*flange_bottom_ydir,
										    flange_edge_z,
										    sweep_flange_pt));
      entire_flanges_ydir.push_back(current_flange);
    }

  // Make webs and flanges pointing in x-direction as entire solids
  web_corner_0 = Point(0.0, -web_thickness/2.0, plate_thickness);
  web_edge_x = SplineCurve(web_corner_0, web_corner_0 + Point(frame_length, 0.0, 0.0));
  web_edge_y = SplineCurve(web_corner_0, web_corner_0 + Point(0.0, web_thickness, 0.0));
  web_edge_z = SplineCurve(web_corner_0, web_corner_0 + Point(0.0, 0.0, stiffener_height - flange_thickness));
  SplineSurface* web_bottom_xdir
    = SweepSurfaceCreator::linearSweptSurface(web_edge_x,
					      web_edge_y,
					      web_corner_0);
  sweep_web_pt = Point(0.0, (frame_spacing-web_thickness)/2.0, plate_thickness);

  flange_corner_0 = Point(0.0, -flange_width/2.0, plate_thickness + stiffener_height - flange_thickness);
  flange_edge_x = SplineCurve(flange_corner_0, flange_corner_0 + Point(frame_length, 0.0, 0.0));
  flange_edge_y = SplineCurve(flange_corner_0, flange_corner_0 + Point(0.0, flange_width, 0.0));
  flange_edge_z = SplineCurve(flange_corner_0, flange_corner_0 + Point(0.0, 0.0, flange_thickness));
  SplineSurface* flange_bottom_xdir
    = SweepSurfaceCreator::linearSweptSurface(flange_edge_x,
					      flange_edge_y,
					      flange_corner_0);
  sweep_flange_pt = Point(0.0, (frame_spacing-flange_width)/2.0, plate_thickness + stiffener_height - flange_thickness);
 
  vector<shared_ptr<SplineVolume> > entire_webs_xdir;
  vector<shared_ptr<SplineVolume> > entire_flanges_xdir;
  stiffener_step = Point(0.0, frame_spacing, 0.0);

  for (int i = 0; i < 3; ++i)
    {
      sweep_web_pt -= stiffener_step;
      sweep_flange_pt -= stiffener_step;

      shared_ptr<SplineVolume> current_web(SweepVolumeCreator::linearSweptVolume(*web_bottom_xdir,
										 web_edge_z,
										 sweep_web_pt));
      entire_webs_xdir.push_back(current_web);
      shared_ptr<SplineVolume> current_flange(SweepVolumeCreator::linearSweptVolume(*flange_bottom_xdir,
										    flange_edge_z,
										    sweep_flange_pt));
      entire_flanges_xdir.push_back(current_flange);
    }


  // Split vectors
  vector<double> flange_split(2);
  vector<double> global_split_ydir;
  vector<double> global_split_xdir;

  flange_split[0] = (flange_width-web_thickness)/2.0;
  flange_split[1] = (flange_width+web_thickness)/2.0;

  double split_position = (frame_spacing-flange_width)/2.0;
  for (int i = 0; i < 3; ++i, split_position += frame_spacing)
    {
      global_split_ydir.push_back(split_position);
      global_split_ydir.push_back(split_position + flange_split[0]);
      global_split_ydir.push_back(split_position + flange_split[1]);
      global_split_ydir.push_back(split_position + flange_width);
    }

  split_position = stiffener_spacing - flange_width/2.0;
  for (int i = 0; i < numb_stiffeners; ++i, split_position += stiffener_spacing)
    {
      global_split_xdir.push_back(split_position);
      global_split_xdir.push_back(split_position + flange_split[0]);
      global_split_xdir.push_back(split_position + flange_split[1]);
      global_split_xdir.push_back(split_position + flange_width);
    }

  // Split volumes
  vector<shared_ptr<SplineVolume> >
    plate_parts_ydir, plates,
    webs_ydir, webs_xdir,
    flange_parts_ydir, flange_parts_xdir, flanges_ydir, flanges_xdir;

  char answ;

  cout << "Split plate (y/n)? ";
  cin >> answ;
  if (answ == 'y')
    {
      cout << "Planes splitting" << endl;
      plate_parts_ydir = split(entire_plate, global_split_xdir, 0);
      plates = split(plate_parts_ydir, global_split_ydir, 1);
      store(plates, os);
    }
  else
    store(entire_plate, os);

  cout << "Split webs y-dir (y/n)? ";
  cin >> answ;
  if (answ == 'y')
    {
      cout << "Webs y" << endl;
      webs_ydir = split(entire_webs_ydir, global_split_ydir, 1);
      store(webs_ydir, os);
    }
  else
    store(entire_webs_ydir, os);

  cout << "Split webs x-dir (y/n)? ";
  cin >> answ;
  if (answ == 'y')
    {
      cout << "Webs x" << endl;
      webs_xdir = split(entire_webs_xdir, global_split_xdir, 0, 11);
      store(webs_xdir, os);
    }
  else
    store(entire_webs_xdir, os);

  cout << "Split flanges y-dir (y/n)? ";
  cin >> answ;
  if (answ == 'y')
    {
      cout << "Flanges y" << endl;
      flange_parts_ydir = split(entire_flanges_ydir, flange_split, 0);
      flanges_ydir = split(flange_parts_ydir, global_split_ydir, 1);
      store(flanges_ydir, os);
    }
  else
    store(entire_flanges_ydir, os);

  cout << "Split flanges x-dir (y/n)? ";
  cin >> answ;
  if (answ == 'y')
    {
      cout << "Flanges x" << endl;
      flange_parts_xdir = split(entire_flanges_xdir, flange_split, 1);
      flanges_xdir = split(flange_parts_xdir, global_split_xdir, 0, 1);
      store(flanges_xdir, os);
    }
  else
    store(entire_flanges_xdir, os);

}
