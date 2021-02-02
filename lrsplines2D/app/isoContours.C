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

#include <fstream>
#include <iostream>
#include <chrono>

#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/lrsplines2D/LRTraceIsocontours.h"


using namespace std;
using namespace Go;

namespace {
  vector<double> contour_vals(const shared_ptr<ParamSurface>& surf, 
			      int num_contours);
}// end anonymous namespace

void print_help_text()
{
  std::cout << "Purpose: Compute contour curves corresponding to an LR B-spline surface \n";
  std::cout << "Mandatory parameters: input surface (.g2), output curves (.g2), number of planes \n";
  std::cout << "Optional input parameters: \n";
  std::cout << "-tol : tolerance for tracing of curves \n";
  std::cout << "-curves2D <filename>: Output curves in the parameter domain of the surface. Default: no 2D curves \n";
  std::cout << "-min <value>: Lower iso value \n";
  std::cout << "-max <value>: Upper iso value \n";
  std::cout << "-knots <integer>: Maximum number of missing knots. Default=100 \n";
  std::cout << "-h or --help : Write this text\n";
}

int fetchIntParameter(int argc, char *argv[], int ki, int& parameter, 
		      int& nmb_par, vector<bool>& par_read)
{
  if (ki == argc-1)
    {
      std::cout << "ERROR: Missing input" << std::endl;
      print_help_text();
      return -1;
    }
  parameter = atoi(argv[ki+1]);
  par_read[ki-1] = par_read[ki] = true;
  nmb_par -= 2;
  return 0;
}

int fetchDoubleParameter(int argc, char *argv[], int ki, double& parameter, 
		      int& nmb_par, vector<bool>& par_read)
{
  if (ki == argc-1)
    {
      std::cout << "ERROR: Missing input" << std::endl;
      print_help_text();
      return -1;
    }
  parameter = atof(argv[ki+1]);
  par_read[ki-1] = par_read[ki] = true;
  nmb_par -= 2;
  return 0;
}

int fetchCharParameter(int argc, char *argv[], int ki, char*& parameter, 
		      int& nmb_par, vector<bool>& par_read)
{
  if (ki == argc-1)
    {
      std::cout << "ERROR: Missing input" << std::endl;
      print_help_text();
      return -1;
    }
  parameter = argv[ki+1];
  par_read[ki-1] = par_read[ki] = true;
  nmb_par -= 2;
  return 0;
}

int main(int argc, char* argv[])
{
  char *surffile = 0;       // Surface output file
  char *cvfile = 0;       // Curve output file
  char *cvfile2 = 0;       // 2D curve output file
  int nmb_iso = 0;
  double tol = 0.1;
  int threshold_missing = 100;
  double min = 1.0;
  double max = -1.0;
  


   int ki, kj;
  vector<bool> par_read(argc-1, false);

  // Read optional parameters
  int nmb_par = argc-1;
  for (ki=1; ki<argc; ++ki)
    {
      string arg(argv[ki]);
      if (arg == "-h" || arg == "--help")
	{
	  print_help_text();
	  exit(0);
	}
      else if (arg == "-tol")
	{
	  int stat = fetchDoubleParameter(argc, argv, ki, tol, 
					  nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-min")
	{
	  int stat = fetchDoubleParameter(argc, argv, ki, min, 
					  nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-max")
	{
	  int stat = fetchDoubleParameter(argc, argv, ki, max, 
					  nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if(arg == "-knots")
 	{
	  int stat = fetchIntParameter(argc, argv, ki, threshold_missing, 
			      nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
     else if(arg == "-curves2D")
 	{
	  int stat = fetchCharParameter(argc, argv, ki, cvfile2, 
					nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
   }
  
  // Read remaining parameters
  if (nmb_par != 3)
    {
      std::cout << "ERROR: Number of parameters is not correct" << std::endl;
      print_help_text();
      return 1;
    }

  for (ki=1; ki<argc; ++ki)
    {
      if (par_read[ki-1])
	continue;
      if (nmb_par == 3)
	{
	  surffile = argv[ki];
	  nmb_par--;
	}
      else if (nmb_par == 2)
	{
	  cvfile = argv[ki];
	  nmb_par--;
	}
      else if (nmb_par == 1)
	{
	  nmb_iso = atoi(argv[ki]);
	  nmb_par--;
	}
    }
  
  // Create the default factory
  GoTools::init();
  Registrator<LRSplineSurface> r293;

  // Read surface
  ifstream is(surffile);
  ObjectHeader header;
  header.read(is);
  shared_ptr<GeomObject> geom_obj(Factory::createObject(header.classType()));
  shared_ptr<ParamSurface> surf =
    dynamic_pointer_cast<ParamSurface,GeomObject>(geom_obj);
  if (!surf.get())
    {
      std::cout << "Object one is not a surface" << std::endl;
      exit(1);
    }
  surf->read(is);
  is.close();

  // Computing isocontours
  vector<double> isovals;
  if (min > max)
    isovals = contour_vals(surf, nmb_iso); // 70
  else
    {
      isovals.resize(nmb_iso);
      double del = (nmb_iso == 1) ? 0 : (max - min)/(double)(nmb_iso-1);
      for (int ka=0; ka<nmb_iso; ++ka)
	isovals[ka] = min + (double)ka*del;
    }

  //const vector<double> isovals {1689.7636324695331};  // this isocontour causes topology problems with surface: data/256_lr_1d.g2
  auto t1 = chrono::high_resolution_clock::now();
  const vector<CurveVec> curves = LRTraceIsocontours(surf,
  						     isovals,
						     threshold_missing,
  						     tol);

  // const auto ssurf = lrsurf.asSplineSurface();
  // const vector<CurveVec> curves = SSurfTraceIsocontours(*ssurf,
  // 							isovals,
  // 							1e-5 * span,
  // 							include_3D,
  // 							use_sisl_marching);

  auto t2 = chrono::high_resolution_clock::now();
  cout << "Curves found in " << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << " milliseconds." << endl;
  
  ofstream os(cvfile);
  //cout << "Number of curves found: " << endl;
  for (size_t i = 0; i != curves.size(); ++i) {
    //cout << "At height " << isovals[i] << ": " << curves[i].size() << " curves." << endl;
    for (auto cv : curves[i]) {
      if (cv.second.get())
	{
	  cv.second->writeStandardHeader(os);
	  cv.second->write(os);
	}
    }
  }
  os.close();

  if (cvfile2 != 0)
    {
      ofstream os2(cvfile2);
      for (size_t i = 0; i != curves.size(); ++i) {
	//cout << "At height " << isovals[i] << ": " << curves[i].size() << " curves." << endl;
	for (auto cv : curves[i]) {
	  if (cv.second.get())
	    {
	      cv.first->writeStandardHeader(os2);
	      cv.first->write(os2);
	    }
	}
      }
      os2.close();
    }
  
  return 0;
}

namespace {


// =============================================================================
  vector<double> contour_vals(const shared_ptr<ParamSurface>& surf, 
			      int num_contours)
// =============================================================================
{
  BoundingBox bbox = surf->boundingBox();

  const double minval = bbox.low()[0];
  const double maxval = bbox.high()[0];
  
  vector<double> result(num_contours, 0);
  for (int i = 0; i != num_contours; ++i)
    result[i] = minval + ((maxval - minval)/(num_contours-1)) * i;

  return result;
}


} // end anonymous namespace
