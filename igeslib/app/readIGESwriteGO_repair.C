//===========================================================================
//                                                                           
// File: readIGESwriteGO_repair.C
//                                                                           
// Created: July 2009
//                                                                           
// Author: Sverre Briseid
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/igeslib/IGESconverter.h"
#include <fstream>
#include <stdlib.h>  // For atof()
#include <memory>
#include "GoTools/creators/CreatorsUtils.h"
#include "GoTools/geometry/SplineDebugUtils.h"
#include "GoTools/geometry/BoundedUtils.h"


using namespace Go;
using std::vector;


int main( int argc, char* argv[] )
{

  std::ifstream infile(argv[1]);
  if (infile.bad()) {
    std::cout << "Infile not found or file corrupt" << std::endl;
    return -1;
  }
  if (argc != 4) {
    std::cout << "Expecting 3 arguments (infile outfile fix_currupt_geometry)."
	      << std::endl;
    return -1;
  }
  IGESconverter conv1;
  try {
      conv1.readIGES(infile);
  } catch (...) {
      std::cout << "Failed reading input IGES file, exiting." << std::endl;
      return -1;  
  }

  std::ofstream outfile(argv[2]);

  std::ofstream outfile_fixed("tmp/fixed_bd_sfs.g2");
  std::ofstream outfile_failures("tmp/failed_bd_sfs.g2");
  std::ofstream outfile_valid("tmp/valid_bd_sfs.g2");
  std::ofstream outfile_trim_space_cvs("tmp/trim_space_cvs.g2");
  std::ofstream outfile_spline_sfs("tmp/spline_sfs.g2");

  int fix_corrupt_geometry(atoi(argv[3]));
  if (fix_corrupt_geometry) {
      // We should run through all objects and fix those that are not
      // valid. Such as loops which are not closed within tolerances,
      // and trimmed surfaces without valid loops.
      MESSAGE("Trying to fix corrupt/invalid geometry!");

      vector<shared_ptr<GeomObject> > geom = conv1.getGoGeom();
      int nmb_bd_sfs = 0;
      int nmb_spline_sfs = 0;
      int nmb_spline_cvs = 0;
      int nmb_cvs_on_sfs = 0;
      int nmb_pt_clouds = 0;
      int nmb_valid_sfs = 0;
      int nmb_fixed_sfs = 0;
      int nmb_failures = 0;
      // int nmb_possibly_fixed_sfs = 0;
      int nmb_invalid_sfs = 0;

//       double max_epsgeo = 1e-02; // Loop tol should not be set larger.

#ifdef SBR_DBG
      std::cout << "Total number of objects: " << geom.size() << std::endl;
#endif
      for (size_t ki = 0; ki < geom.size(); ++ki) {
	  if (geom[ki].get() == NULL)
	      MESSAGE("Object missing!");
	  else if (geom[ki]->instanceType() == Class_BoundedSurface) {
	      ++nmb_bd_sfs;
	      shared_ptr<BoundedSurface> bd_sf =
		  dynamic_pointer_cast<BoundedSurface, GeomObject>(geom[ki]);

	      std::cout << "Object number: " << ki << std::endl;
	      MESSAGE("Object number: " << ki << " (type: " <<
		      bd_sf->underlyingSurface()->instanceType() << ")");

	      int init_state = 0;
	      bool sf_ok = bd_sf->isValid(init_state);

	      if (sf_ok)
	      {
		  ASSERT(init_state > 0);
		  MESSAGE("State: Input surface ok, nothing to be done!");
		  bd_sf->writeStandardHeader(outfile_valid);
		  bd_sf->write(outfile_valid);
		  ++nmb_valid_sfs;
		  continue;
	      }
	      else
	      {
		  BoundedUtils::fixInvalidBoundedSurface(bd_sf);
		  int state = 0;
		  sf_ok = bd_sf->isValid(state);
		  if (!sf_ok)
		      MESSAGE("Failed fixing bd_sf!");
	      }
#if 0
	      std::ofstream outfile_curr_bd_sf("tmp/curr_bd_sf.g2");
	      shared_ptr<ParamSurface> under_sf = bd_sf->underlyingSurface();
	      under_sf->writeStandardHeader(outfile_spline_sfs);
	      under_sf->write(outfile_spline_sfs);
	      SplineDebugUtils::writeTrimmedInfo(*bd_sf, outfile_curr_bd_sf, 0.0);
	      double debug_val = 0.0;
#endif
	  } else if (geom[ki]->instanceType() == Class_SplineSurface) {
	      ++nmb_spline_sfs;
	      geom[ki]->writeStandardHeader(outfile_spline_sfs);
	      geom[ki]->write(outfile_spline_sfs);
	  } else if (geom[ki]->instanceType() == Class_SplineCurve)
	      ++nmb_spline_cvs;
	  else if (geom[ki]->instanceType() == Class_CurveOnSurface)
	      ++nmb_cvs_on_sfs;
	  else if (geom[ki]->instanceType() == Class_PointCloud)
	      ++nmb_pt_clouds;
      }

      std::cout << "# valid sfs: " << nmb_valid_sfs << ", # invalid sfs: " <<
	  nmb_invalid_sfs << ", # fixed sfs: " << nmb_fixed_sfs <<
	  ", # failures: " << nmb_failures << std::endl;

      int nmb_other_objs = (int)geom.size() - nmb_spline_sfs - nmb_bd_sfs -
	  nmb_spline_cvs - nmb_cvs_on_sfs - nmb_pt_clouds;
      std::cout << "# spline_sfs: " << nmb_spline_sfs << ", # bd_sfs: " <<
	  nmb_bd_sfs << ", nmb_spline_cvs: " << nmb_spline_cvs <<
	  ", nmb_cvs_on_sfs: " << nmb_cvs_on_sfs << ", # pt_clouds: " <<
	  nmb_pt_clouds << ", # other_objs: " << nmb_other_objs << std::endl;
  }

  conv1.writego(outfile);

}
