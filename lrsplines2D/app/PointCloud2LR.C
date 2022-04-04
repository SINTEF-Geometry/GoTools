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

#include "GoTools/utils/config.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/FileUtils.h"
#include "GoTools/utils/Array.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSurfApprox.h"
#include "GoTools/lrsplines2D/LRApproxApp.h"
#include "GoTools/utils/timeutils.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <sys/stat.h>
#include <boost/timer.hpp>
#include <time.h>

//#define DEBUG
//#define DEBUG_EL
//#define DEBUG2


using namespace Go;
using std::vector;
using std::string;

void print_help_text()
{
  std::cout << "Purpose: Approximate a point cloud by an LR B-spline surface. \n";
  std::cout << "Mandatory parameters: input point cloud (.txt, .xyz or .g2), output surface (.g2), tolerance, number of iterations. \n";
  std::cout << "An adaptive approximation procedure is applied which for the";
  std::cout << " specified number of iterations: \n";
  std::cout << " - Approximates the points with a surface in the current spline space \n";
  std::cout << " - Computes the approximation accuracy \n";
  std::cout << " - Refines the surfaces in areas where the tolerance is not met \n";
  std::cout << "The approximation is completed when the tolerance is met or the";
  std::cout << " specified number of iterations is exceeded.";
  std::cout << "The number of iterations is recommended to lie in the interval [4:7]. \n";
  std::cout << "The points are expected to be given as x, y, z and be parameterized on x and y, but 3D parameterized points are accepted.\n";
  std::cout << "Then the points are given as u, v, x, y, z. \n";
  std::cout << "Optional input parameters: \n";
  std::cout << "-par <0/1> : True (1) if parameterized points are given. Default value is 0.\n";
  std::cout << "-dist <filename (.txt)> : Write distance field to file (x, y, z, distance) \n";
  std::cout << "-info <filename> : Write accuracy information to file \n";
  std::cout << "-smooth <weight> : Overrule default smoothing weight (1.0e-10) in least squares approximation \n";
  std::cout << "-mba <0/1/n/-1>: 0 = use only least squares approximation \n";
  std::cout << "                 1 = use only multilevel B-spline approximation (MBA) \n";
  std::cout << "                 n = start with least squares, turn to MBA after n iterations \n";
  std::cout << "                -1 = initiate computation using MBA \n";
  std::cout << "Default setting is start with least squares, turn to MBA for the last iterations \n";
  std::cout << "-degree <polynomial degree> : 2 or 3 recommended \n";
  std::cout << "-nmb_coef <initial value> : Initial number of coefficients in each parameter direction \n";
  std::cout << "-distributecf <0/1> : Modify initial number of coefficients according to relative size of parameter domain \n";
  std::cout << "-outlier <0/1>: Flag for removal of outliers (0=false, 1=true). Default false \n";
  std::cout << "-minsize <size> : Minimum element size, all directions \n";
  std::cout << "-reltol <0/1>: Apply relative tolerance flag. Default false \n";
  std::cout << "-tolfac1: Factor for modification of tolerance, positive heights. Default 0.0 \n";
  std::cout << "-tolfac2: Factor for modification of tolerance, negative heights. Default 0.0 \n";
  std::cout << "-signpoints: Filename significant points, same formats as point cloud \n";
  std::cout << "-signtol: Tolerance for significant points, default 0.5*tolerance \n";
  std::cout << "-signpost: Flag for post processing significant points outside tolerance (0=false, 1=true). Default false \n";
  std::cout << "-verbose <0/1>: Write accuracy information at each iteration level. Default = 0 \n";
  std::cout << "-tolfile: File specifying domains with specific tolerances, global tolerance apply outside domains. PointCloud2LR -tolfile for file format \n";
  std::cout << "-toldoc: Documentation on file format for tolerance domains. \n";
  std::cout << "-featuredoc: Show feature documentation \n";
  std::cout << "-refstrat: Print parameter information to define refinement strategy. \n";
  std::cout << "-AIC <filename>: Compute AIC at each iteration level and output to the given file. Time consuming. Default: No AIC computations \n";
  std::cout << "-h or --help : Write this text\n";
}

void print_refine_info()
{
  std::cout << "Define refinement strategy from command line with the following parameters: \n";
  std::cout << "-refcat <code>: F = full span, Ml/Mu/Mc = minimum span, S = structured mesh, \n";
  std::cout << "R = restricted mesh, RL = restricted mesh with element extension. Default = F. \n";
  std::cout << "-alterdir <code>: B = refine in both parameter directions, A = refine in alternating directions. Default is A. \n";
  std::cout << "-threshold <code>: no = none, td = with respect to distance, tn = with respect to number. Default is tn. \n";
  std::cout << "tdk = distance and number for category R and S. Default = no \n";
  std::cout << "-swapcat <size>: Swap strategy when approximation efficience drops below size. Default = -100 \n";
  std::cout << "-refcat2 <code>: After swap, as -refcat \n";
  std::cout << "-threshold2 <code>: After swap, as -threshold \n";
}

void print_feature_info()
{
  std::cout << "-feature: <resolution> : Command line parameter to write feature information to file according to given grid resolution \n \n";
  std::cout << "Compute grid based feature information for every iteration level. \n";
  std::cout << "The information is stored in files called cellinfox where the number x represent the iteration level." << std::endl;
  std::cout << "The surface is parameterized on x and y and the surface value represents intensity/height" << std::endl;
  std::cout << "The output files are text based and must be translated to tiff with a different application. \n";
  std::cout << "Each line represent one grid point and the coloumns store the following information: \n";
  std::cout << "0: Average slope in cell (9 samples) \n";
  std::cout << "1: Average value of surface in cell (9 samples)  \n";
  std::cout << "2: Maximum difference of surface values in cell (9 samples)  \n";
  std::cout << "3: Average distance between surface and points for each cell \n";
  std::cout << "4: Maximum distance between surface and points in cell \n";
  std::cout << "5: Average intensity/height value of points in cell \n";
  std::cout << "6: Maximum difference of intensity values in cell \n";
  std::cout << "7: Standard deviation of distances between point cloud and surface in cell \n";
  std::cout << "8: Standard deviation of intensity values in cell \n";
  std::cout << "9: Average distance between surface and points in cell divided by maximum distance \n";
  std::cout << "10: Maximum difference between signed distances between points and surface in cell \n";
  std::cout << "11: Average distance between points with higher intensity than the surface and surface in cell \n";
  std::cout << "12: Average distance between points with lower intensity than the surface and surface in cell \n";
  std::cout << "13: Number of point with lower intensity than the surface where the intensity difference is larger than threshold divided by the number of points in the cell \n";
  std::cout << "14: Number of point with higher intensity than the surface where the intensity difference is larger than threshold divided by the number of points in the cell \n";
  std::cout << "15: Number of surface elements in cell \n";
  std::cout << "16: Average laplacian in cell (9 samples) \n";
  std::cout << "The entries are scaled to represent a number in the range [0,10]. \n";
}

void print_tol_file_format()
{
  std::cout << "File specifying domains/boxes with different tolerances. \n";
  std::cout << "If not otherwise stated, the global tolerance will applies outside the boxes \n";
  std::cout << "Numbers are given as floats and separated by space \n";
  std::cout << "Line1: Number of boxes (integer), whether or not the tolerance is specified outside the boxes (0/1) \n";
  std::cout << "Following lines: \n";
  std::cout << "xmin ymin xmax ymax tolerance \n";
  std::cout << "Positive numbers for tolerance means absolute value, \n";
  std::cout << "negative numbers mean multiplication factor to standard deviation of points \n";
  std::cout << "Last line (if given): Tolerance, positive or negative float as for previous lines. Value overrules global tolerance. \n";
  std::cout << "Ensure non-overlapping boxes. No test applied. \n";
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

int main(int argc, char *argv[])
{
  char* input_type = 0;    // Type of point file
  char *pointfile = 0;     // Input point file
  char *surffile = 0;       // Surface output file
  char *infofile = 0;      // Accuracy information output file
  char *tolfile = 0;       // File specifying varying tolerances
  int del = 3;             // Number of entries for each point
  double AEPSGE = 0.5;     // Requested accuracy
  int max_iter = 6;        // Maximum number of iterations in adaptive alogrithm
  char *field_out = 0;     // Distance field output file
  double smoothwg = 1.0e-9; 
  int initmba = 0; //1;  // Initiate surface using the mba method
  int mba = 0;      // Use least squares approximation
  int tomba = std::min(5, max_iter-1);    // Turn to the mba method at 
  // iteration level 5 or in the last iteration
  int degree = 2;
  int outlierflag = 0;
  int reltol = 0;
  double tolfac1 = 0.0, tolfac2 = 0.0;
  char *signpointfile = 0;  // Input significant points
  double signtol = -1.0;  // Tolerance for significant points
  int signpost = 0;  // Flag for post procession of significant points
  double minsize = -1.0;
  int feature_out = -1;
  int verbose = 0;
  int AIC = 0;
  char *AIC_file = 0;
  
  int initncoef = 10;
  int distribute_ncoef = 0;

  int refcat1=1, refcat2=0, threshold1=2, threshold2=-1, alter=1;
  double swap = -100.0;

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
      else if (arg == "-toldoc")
	{
	  print_tol_file_format();
	  exit(0);
	}
      else if (arg == "-featuredoc")
	{
	  print_feature_info();
	  exit(0);
	}
      else if (arg == "-refstrat")
	{
	  print_refine_info();
	  exit(0);
	}
      else if (arg == "-par")
	{
	  int tmp;
	  int stat = fetchIntParameter(argc, argv, ki, tmp, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	  del = (tmp == 1) ? 5 : 3;
	}
      else if (arg == "-dist")
	{
	  int stat = fetchCharParameter(argc, argv, ki, field_out, 
					nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-info")
	{
	  int stat = fetchCharParameter(argc, argv, ki, infofile, 
					nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-smooth")
	{
	  int stat = fetchDoubleParameter(argc, argv, ki, smoothwg, 
					  nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-mba")
	{
	  int mm;
	  int stat = fetchIntParameter(argc, argv, ki, mm, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	  if (mm == 0)
	    tomba = 100;
	  else if (mm == 1)
	    mba = 1;
	  else if (mm < 0)
	    initmba = 1;
	  else
	    tomba = mm;
	}
      else if (arg == "-degree")
	{
	  int stat = fetchIntParameter(argc, argv, ki, degree, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-nmb_coef")
	{
	  int stat = fetchIntParameter(argc, argv, ki, initncoef, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-distributecf")
	{
	  int stat = fetchIntParameter(argc, argv, ki, distribute_ncoef, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
       else if (arg == "-outlier")
	{
	  int stat = fetchIntParameter(argc, argv, ki, outlierflag, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-minsize")
	{
	  int stat = fetchDoubleParameter(argc, argv, ki, minsize, 
					  nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-reltol")
	{
	  int stat = fetchIntParameter(argc, argv, ki, reltol, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-tolfac1")
	{
	  int stat = fetchDoubleParameter(argc, argv, ki, tolfac1, 
					  nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-tolfac2")
	{
	  int stat = fetchDoubleParameter(argc, argv, ki, tolfac2, 
					  nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-signpoints")
	{
	  int stat = fetchCharParameter(argc, argv, ki, signpointfile, 
					nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-signtol")
	{
	  int stat = fetchDoubleParameter(argc, argv, ki, signtol, 
					  nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-signpost")
	{
	  int stat = fetchIntParameter(argc, argv, ki, signpost, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-verbose")
	{
	  int stat = fetchIntParameter(argc, argv, ki, verbose, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-tolfile")
	{
	  int stat = fetchCharParameter(argc, argv, ki, tolfile, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-feature")
	{
	  int stat = fetchIntParameter(argc, argv, ki, feature_out, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-refcat")
	{
	  char* cat0;
	  int stat = fetchCharParameter(argc, argv, ki, cat0, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;

	  string cat(cat0);
	  if (cat == "F")
	    refcat1 = 1;
	  else if (cat == "Ml")
	    refcat1 = 2;
	  else if (cat == "Mu")
	    refcat1 = 3;
	  else if (cat == "Mc")
	    refcat1 = 4;
	  else if (cat == "S")
	    refcat1 = 5;
	  else if (cat == "R")
	    refcat1 = 6;
	  else if (cat == "RL")
	    refcat1 = 7;
	  if (refcat2 < 1 || refcat2 > 7)
	    refcat2 = refcat1;
	}
      else if (arg == "-refcat2")
	{
	  char* cat0;
	  int stat = fetchCharParameter(argc, argv, ki, cat0, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;

	  string cat(cat0);
	  if (cat == "F")
	    refcat2 = 1;
	  else if (cat == "Ml")
	    refcat2 = 2;
	  else if (cat == "Mu")
	    refcat2 = 3;
	  else if (cat == "Mc")
	    refcat2 = 4;
	  else if (cat == "S")
	    refcat2 = 5;
	  else if (cat == "R")
	    refcat2 = 6;
	  else if (cat == "RL")
	    refcat2 = 7;
	}
      else if (arg == "-threshold")
	{
	  char* thresh0;
	  int stat = fetchCharParameter(argc, argv, ki, thresh0, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;

	  string thresh(thresh0);
	  if (thresh == "no")
	    threshold1 = 0;
	  else if (thresh == "td")
	    threshold1 = 1;
	  else if (thresh == "tn")
	    threshold1 = (refcat1 >= 6) ? 3 : 2;
	  else if (thresh == "tdk")
	    threshold1 = 4;
	  if (threshold2 < 0)
	    threshold2 = threshold1;
	}
      else if (arg == "-threshold2")
	{
	  char* thresh0;
	  int stat = fetchCharParameter(argc, argv, ki, thresh0, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;

	  string thresh(thresh0);
	  if (thresh == "no")
	    threshold2 = 0;
	  else if (thresh == "td")
	    threshold2 = 1;
	  else if (thresh == "tn")
	    threshold2 = (refcat1 >= 6) ? 3 : 2;
	  else if (thresh == "tdk")
	    threshold2 = 4;
	}
      else if (arg == "-alterdir")
	{
	  char *altdir;
	  int stat = fetchCharParameter(argc, argv, ki, altdir, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	  alter = (altdir[0] == 'A') ? 1 : 0;
	}
      else if (arg == "-swap")
	{
	  int stat = fetchDoubleParameter(argc, argv, ki, swap, 
					  nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-AIC")
	{
	  AIC = 1;
	  int stat = fetchCharParameter(argc, argv, ki, AIC_file, 
					nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
    }

  // Read remaining parameters
  if (nmb_par != 4)
    {
      std::cout << "ERROR: Number of parameters is not correct" << std::endl;
      print_help_text();
      return 1;
    }

  for (ki=1; ki<argc; ++ki)
    {
      if (par_read[ki-1])
	continue;
      if (nmb_par == 4)
	{
	  pointfile = argv[ki];
	  nmb_par--;
	}
      else if (nmb_par == 3)
	{
	  surffile = argv[ki];
	  nmb_par--;
	}
      else if (nmb_par == 2)
	{
	  AEPSGE = atof(argv[ki]);
	  nmb_par--;
	}
      else
	{
	  max_iter = atoi(argv[ki]);
	}
    }
  if (signtol < 0)
    signtol = 0.5*AEPSGE;  // Default value

  // Read point cloud
  vector<double> data;
  vector<double> extent(2*del);   // Limits for points in all coordinates
  // Possible types of input files
  char keys[6][8] = {"g2", "txt", "TXT", "xyz", "XYZ", "dat"};
  int ptstype = FileUtils::fileType(pointfile, keys, 6);
  if (ptstype < 0)
    {
      std::cout << "ERROR: File type not recognized" << std::endl;
      return 1;
    }

  int nmb_pts = 0;
  std::ifstream pointsin(pointfile);
  if (ptstype == 0)
    {
     // g2 file
      ObjectHeader header;
      PointCloud3D points;
      try {
	header.read(pointsin);
	points.read(pointsin);
      }
      catch (...)
	{
	  std::cout << "ERROR: Not a valid point file" << std::endl;
	  return -1;
	}
      BoundingBox box = points.boundingBox();
      Point low = box.low();
      Point high = box.high();
      nmb_pts = points.numPoints();
      data.insert(data.end(), points.rawData(), points.rawData()+3*nmb_pts);
      for (int ki=0; ki<3; ++ki)
	{
	  extent[2*ki] = low[ki];
	  extent[2*ki+1] = high[ki];
	}
      std::cout << "Domain: [" << extent[0] << ", " << extent[1] << "] x [" << extent[2];
      std::cout << ", " << extent[3] << "]" << std::endl;
      std::cout << "Range: " << extent[4] << " - " << extent[5] << std::endl;
    }
  else
    FileUtils::readTxtPointFile(pointsin, del, data, nmb_pts, extent);

  int nmb_sign = 0;
  vector<double> sign_data;
  vector<double> sign_extent(2*del);
  if (signpointfile != 0)
    {
      int sgntype = FileUtils::fileType(signpointfile, keys, 6);
      if (sgntype < 0)
	{
	  std::cout << "ERROR: File type not recognized" << std::endl;
	  return 1;
	}

      std::ifstream signpointsin(signpointfile);
      if (sgntype == 0)
	{
	  // g2 file
	  ObjectHeader header;
	  PointCloud3D points;
	  try {
	    header.read(signpointsin);
	    points.read(signpointsin);
	  }
	  catch (...)
	    {
	      std::cout << "ERROR: Not a valid point file" << std::endl;
	      return -1;
	    }
	  BoundingBox box = points.boundingBox();
	  Point low = box.low();
	  Point high = box.high();
	  nmb_sign = points.numPoints();
	  sign_data.insert(sign_data.end(), points.rawData(), points.rawData()+3*nmb_sign);
	  for (int ki=0; ki<3; ++ki)
	    {
	      sign_extent[2*ki] = low[ki];
	      sign_extent[2*ki+1] = high[ki];
	    }
	}
      else
	FileUtils::readTxtPointFile(signpointsin, del, sign_data,
				    nmb_sign, sign_extent);
    }

  if (nmb_sign > 0)
    {
      // Modify data extent
      for (int ka=0; ka<2*del; ka+=2)
	{
	  extent[ka] = std::min(extent[ka], sign_extent[ka]);
	  extent[ka+1] = std::max(extent[ka+1], sign_extent[ka+1]);
	}
    }
  
  bool use_stdd = false;
  vector<LRSurfApprox::TolBox> tolerances;
  if (tolfile != 0)
    {
      std::ifstream tolin(tolfile);
      int nmb_box, last;
      double umin, umax, vmin, vmax, tol;
      tolin >> nmb_box >> last;
      tolerances.resize(nmb_box);
      for (int ka=0; ka<nmb_box; ++ka)
	{
	  tolin >> umin >> umax >> vmin >> vmax >> tol;
	  tolerances[ka].setVal(std::max(umin,extent[0]), std::min(umax,extent[1]),
				std::max(vmin,extent[2]), std::min(vmax,extent[3]), tol);
	  if (tol < 0.0)
	    use_stdd = true;
	}
      if (last > 0)
	{
	  tolin >> AEPSGE;
	  if (AEPSGE < 0.0)
	    use_stdd = true;
	}
    }
  if (AEPSGE < 0.0)
    use_stdd = true;
	  
#ifdef DEBUG2
  double maxtol = 0.0;
  double mintol = 1.0e8;
  double minheight = 1.0e8;
  double maxheight = -minheight;
  for (ki=0; ki<nmb_pts; ++ki)
    {
      double height = data[del*ki+del-1];
      double tol = (height < 0.0) ? AEPSGE - tolfac2*height :
	AEPSGE + tolfac1*height;
      minheight = std::min(height, minheight);
      maxheight = std::max(height, maxheight);
      mintol = std::min(tol, mintol);
      maxtol = std::max(tol, maxtol);
    }
  std::cout << "Elevation: [" << minheight << ", " << maxheight << "]" << std::endl;
    std::cout << "Tolerances: [" << mintol << ", " << maxtol << "]" << std::endl;
#endif

  // Move point cloud to origo
  Point mid(0.5*(extent[2*(del-3)] + extent[2*(del-3)+1]),
	    0.5*(extent[2*(del-2)] + extent[2*(del-2)+1]), 0.0);
  for (ki=0; ki<nmb_pts; ++ki)
    for (kj=del-3; kj<del-1; ++kj)
      {
	data[del*ki+kj] -= mid[kj-del+3];
      }
  for (ki=0; ki<nmb_sign; ++ki)
    for (kj=del-3; kj<del-1; ++kj)
      {
	sign_data[del*ki+kj] -= mid[kj-del+3];
      }

  for (size_t kj=0; kj<tolerances.size(); ++kj)
    tolerances[kj].translateBox(-mid[0], -mid[1]);

  if (/*true)/*/use_stdd)
    {
      double avheight = 0.0;
      for (ki=0; ki<nmb_pts; ++ki)
	{
	  double height = data[del*ki+del-1];
	  avheight += (height/(double)nmb_pts);
	}
      double stdd = 0.0;
     for (ki=0; ki<nmb_pts; ++ki)
	{
	  double height = data[del*ki+del-1];
	  stdd += (pow(avheight-height,2)/(double)nmb_pts);
	}
     stdd = sqrt(stdd);
     if (use_stdd)
       {
	 for (size_t kj=0; kj<tolerances.size(); ++kj)
	   {
	     if (tolerances[kj].tol < 0.0)
	       tolerances[kj].setTol(fabs(tolerances[kj].tol)*stdd);
	   }
	 if (AEPSGE < 0.0)
	   AEPSGE = fabs(AEPSGE)*stdd;
       }
     std::cout << "Standard deviation: " << stdd << std::endl;
    }
     

#ifdef DEBUG
  // Write translated surface and points to g2 format
  vector<double> data2;
  data2.reserve(nmb_pts*3);
  for (ki=0, kj=0; ki<nmb_pts; ++ki, kj+=del)
    data2.insert(data2.end(), data.begin()+kj, data.begin()+kj+3);
  PointCloud3D cloud(data2.begin(), nmb_pts);

  std::ofstream of1("translated_sf.g2");
  std::ofstream of2("translated_points.g2");
  cloud.writeStandardHeader(of2);
  cloud.write(of2);
#endif
  

  time_t start = time(NULL);



 boost::timer t;
  double duration;

  t.restart();

  // if (del > 3)
  //   {
  //     initmba = 0;
  //     mba = 0;
  //   }
  int order = degree + 1; 
  int nmb_coef = std::max(order, initncoef); //std::max(order, 14);
  int nm = nmb_coef*nmb_coef;
  double dom = (extent[1]-extent[0])*(extent[3]-extent[2]);
  double c1 = std::pow((double)nm/dom, 1.0/2.0);
  int nc[2];
  for (kj=0; kj<2; ++kj)
    {
      double tmp = c1*(extent[2*kj+1]-extent[2*kj]);
      nc[kj] = (int)tmp;
      //if (tmp - (double)nc[kj] > 0.5)
      ++nc[kj];
      nc[kj] = std::max(nc[kj], order);
    }
  double mba_coef = 0.0;
  if (initmba)
    mba_coef = 0.5*(extent[2*(del-1)] + extent[2*(del-1)+1]);
  shared_ptr<LRSurfApprox> approx;
  if (distribute_ncoef)
    approx = shared_ptr<LRSurfApprox>(new LRSurfApprox(nc[0], order, nc[1], order, data, del-2, 
						       AEPSGE, initmba ? true : false, mba_coef,
						       true, true));
  else
    approx = shared_ptr<LRSurfApprox>(new LRSurfApprox(nmb_coef, order, nmb_coef, order, data, del-2, 
						       AEPSGE, initmba ? true : false, mba_coef,
						       true, true));
  approx->setSmoothingWeight(smoothwg);
  approx->setSmoothBoundary(true);
  if (mba)
    approx->setUseMBA(true);
  else
    {
      // if (initmba)
      // 	approx->setInitMBA(initmba, 0.5*(low[2]+high[2]));
      //if (del == 3)
	approx->setSwitchToMBA(tomba);
	approx->setMakeGhostPoints(false /*true*/);
    }
  if (outlierflag > 0)
    approx->setOutlierFlag(true);
  if (minsize > 0.0)
    approx->setMinimumElementSize(minsize, minsize);
  if (reltol > 0)
    {
      double tol1 = extent[4]<0.0 ? AEPSGE - tolfac2*extent[4] : AEPSGE + tolfac1*extent[4];
      double tol2 = extent[5]<0.0 ? AEPSGE - tolfac2*extent[5] : AEPSGE + tolfac1*extent[5];
      std::cout << "Variable tolerance: " << tol1 << " - " << tol2 << std::endl;
      approx->setVarTol(tolfac1, tolfac2);
    }

  if (sign_data.size() > 0)
    {
      approx->addSignificantPoints(sign_data, signtol);
      double sign_fac = 200.0; //100.0; //20.0; //1.0; //10.0;
      approx->setSignificantFactor(sign_fac);
    }

  if (tolerances.size() > 0)
    approx->setVarTolBox(tolerances);
  
  approx->setVerbose(verbose);

  // Feature output
  if (feature_out > 0)
    approx->setFeatureOut(feature_out);

  // Refine strategy
  approx->setRefinementStrategy(refcat1, alter, threshold1, swap, refcat2, threshold2);

  if (del == 3)
    {
      double zrange = extent[5] - extent[4];
      double zfac = std::max(AEPSGE, 0.005*zrange);
      approx->addLowerConstraint(extent[4] - zfac);
      approx->addUpperConstraint(extent[5] + zfac);
      approx->setLocalConstraint(zfac);
    }

  if (AIC)
    approx->setAIC(true);
  
#ifdef DEBUG
  std::cout << "Range: " << extent[1]-extent[0] << ", " << extent[3]-extent[2];
  std::cout << ", " << extent[5]-extent[4] << std::endl;
  std::cout << "Heights: " << extent[4] << " " << extent[5] << std::endl;
#endif
  double maxdist, avdist, avdist_total; // will be set below
  int nmb_out_eps;        // will be set below
  double maxout, avout;
  shared_ptr<LRSplineSurface> surf;
  try {
  surf = approx->getApproxSurf(maxdist, avdist_total,
			      avdist, nmb_out_eps, 
			      max_iter);
  }
  catch (...)
    {
      std::cout << "ERROR: Surface approximation failed" << std::endl;
      return 1;
    }

  approx->fetchOutsideTolInfo(maxout, avout);


  duration = t.elapsed();
  std::cout << "Duration: " << duration << std::endl;
  double min = floor(duration/60);
  double sec = duration - 60*min;
  std::cout << min << "m" << sec << "s" << std::endl;
  time_t end = time(NULL);
  std::cout<<"Execution Time: "<< (double)(end-start)<<" Seconds"<<std::endl;



  double maxdist_sign, avdist_sign;
  int nmb_outside_sign;
  if (sign_data.size() > 0)
    {
      // Fetch information about significant point approximation
      approx->getSignificantPointInfo(maxdist_sign, avdist_sign,
				     nmb_outside_sign);
      if (maxdist_sign > signtol && signpost > 0)
	{
	  // Post process approximation to reach tolerance in significant
	  // points
	  double fac1 = 0.0, fac2 = 0.0;
	  LRApproxApp::updateSurfWithSignificantPts(surf, AEPSGE, signtol,
						    fac1, fac2,
						    maxdist, avdist_total,
						    avdist, nmb_out_eps, 
						    maxdist_sign, avdist_sign,
						    nmb_outside_sign);
	}
    }

  int nmb_outlier = approx->getNmbOutliers();
  if (nmb_outlier > 0)
    {
      vector<double> outliers, regular;
      int nmb_outliers = 0, nmb_regular = 0;
      approx->getClassifiedPts(outliers, nmb_outliers, regular, nmb_regular);
      for (ki=0; ki<nmb_outliers; ++ki)
	{
	  for (kj=del-3; kj<del-1; ++kj)
	    outliers[del*ki+kj] += mid[kj-del+3];
	}
      for (ki=0; ki<nmb_regular; ++ki)
	{
	  for (kj=del-3; kj<del-1; ++kj)
	    regular[del*ki+kj] += mid[kj-del+3];
	}
      PointCloud3D outlier_cloud(outliers.begin(), nmb_outlier);
      PointCloud3D regular_cloud(regular.begin(), nmb_regular);
      std::ofstream ofo("outlier_pts.g2");
      ofo.precision(15);
      outlier_cloud.writeStandardHeader(ofo);
      outlier_cloud.write(ofo);
      std::ofstream ofr("regular_pts.g2");
      ofr.precision(15);
      regular_cloud.writeStandardHeader(ofo);
      regular_cloud.write(ofr);
    }
      
  if (AIC)
    {
      // Fetch Mahalanobis distance and log likelyhood
      vector<double> AICinfo;
      vector<int> ncoef;
      approx->getAICInfo(AICinfo, ncoef);
      std::ofstream ofAIC(AIC_file);
      ofAIC.precision(10);
      ofAIC << "For each iteration level the following is reported: " << std::endl;
      ofAIC << "level, number of coefficients, Mahalanobix distance with students t-distribution, \n";
      ofAIC << "AIC student t-distribution, Mahalanobix distance normal distribution, AIC normal distribution" << std::endl;
      for (size_t ki=0; ki<ncoef.size(); ++ki)
	{
	  ofAIC << ki << "\t" << ncoef[ki];
	  for (size_t kj=0; kj<4; ++kj)
	    ofAIC << "\t" << AICinfo[4*ki+kj];
	  ofAIC << std::endl;
	}
    }
 
  std::ofstream sfout(surffile);     // Surface output stream

  if (infofile)
    {
      std::ofstream infoout(infofile);   // Accuracy information output stream

      infoout << "Total number of points: " << nmb_pts << std::endl;
      infoout << "Number of elements: " << surf->numElements() << std::endl;
      infoout << "Number of coefficients: " << surf->numBasisFunctions() << std::endl;
      infoout << "Maximum distance: " << maxdist << std::endl;
      infoout << "Average distance: " << avdist_total << std::endl;
      infoout << "Average distance for points outside of the tolerance: " << avdist << std::endl;
      infoout << "Number of points outside the tolerance: " << nmb_out_eps << std::endl;
      infoout << "Maximum distance exceeding tolerance (dist-tol): " << maxout << std::endl;
      infoout << "Average distance exceeding tolerance (dist-tol): " << avout << std::endl;
      infoout << "Number of classified outliers: " << nmb_outlier << std::endl;
      if (sign_data.size() > 0)
	{
	  infoout << "INFO: Maximum distance in significant points: " << maxdist_sign << std::endl;
	  infoout << "INFO: Average distance in significant points: " << avdist_sign << std::endl;
	  infoout << "INFO: Number of significant points outside tolerance: " << nmb_outside_sign << std::endl;
	}
    } 
  else
    {
      std::cout << "INFO: Total number of points: " << nmb_pts << std::endl;
      std::cout << "INFO: Number of elements: " << surf->numElements() << std::endl;
      std::cout << "INFO: Number of coefficients: " << surf->numBasisFunctions() << std::endl;
      std::cout << "INFO: Maximum distance: " << maxdist << std::endl;
      std::cout << "INFO: Average distance: " << avdist_total << std::endl;
      std::cout << "INFO: Average distance for points outside of the tolerance: " << avdist << std::endl;
      std::cout << "INFO: Number of points outside the tolerance: " << nmb_out_eps << std::endl;
      std::cout << "INFO: Maximum distance exceeding tolerance (dist-tol): " << maxout << std::endl;
      std::cout << "INFO: Average distance exceeding tolerance (dist-tol): " << avout << std::endl;
      std::cout << "INFO: Number of classified outliers: " << nmb_outlier << std::endl;
      if (sign_data.size() > 0)
	{
	  std::cout << "INFO: Maximum distance in significant points: " << maxdist_sign << std::endl;
	  std::cout << "INFO: Average distance in significant points: " << avdist_sign << std::endl;
	  std::cout << "INFO: Number of significant points outside tolerance: " << nmb_outside_sign << std::endl;
	}
    } 


  if (surf.get())
    {
#ifdef DEBUG
      std::ofstream of1("translated_sf_3d.g2");
      if (surf->dimension() == 3)
	{
	  surf->writeStandardHeader(of1);
	  surf->write(of1);
	}
      else
	{
	  shared_ptr<LRSplineSurface> surf2(surf->clone());
	  surf2->to3D();
	  surf2->writeStandardHeader(of1);
	  surf2->write(of1);
	}
#endif
	  
      // Translate
      if (surf->dimension() == 3)
	{
	  surf->translate(mid);
	}
      else
	{
	  // Update parameter domain
	  double umin = surf->paramMin(XFIXED);
	  double umax = surf->paramMax(XFIXED);
	  double vmin = surf->paramMin(YFIXED);
	  double vmax = surf->paramMax(YFIXED);

	  surf->setParameterDomain(umin + mid[0], umax + mid[0],
				   vmin + mid[1], vmax + mid[1]);
	}

      surf->writeStandardHeader(sfout);
      surf->write(sfout);

      if (field_out)
	{
	  // Fetch data points with distance information
	  vector<double> pnts_dist;
	  pnts_dist.reserve(4*nmb_pts);
	  try {
	   LRSplineSurface::ElementMap::const_iterator elem = surf->elementsBegin();
	   LRSplineSurface::ElementMap::const_iterator last = surf->elementsEnd();
	  for (; elem != last; ++elem)
	    {
	      if (!elem->second->hasDataPoints())
		continue;
	      vector<double>& points = elem->second->getDataPoints();
	      pnts_dist.insert(pnts_dist.end(), points.begin(), points.end());
	    }
	  }
	  catch (...)
	    {
	      std::cout << "ERROR: Extraction of distance information failed" << std::endl;
	      return 1;
	    }

	  // Translate to initial domain
	  for (size_t kj=0; kj<pnts_dist.size(); kj+=4)
	    {
	      pnts_dist[kj] += mid[0];
	      pnts_dist[kj+1] += mid[1];
	    }

	  // Write to file
	  std::ofstream field_info(field_out);
	  (void)field_info.precision(15);
	  for (size_t kj=0; kj<pnts_dist.size(); kj+=4)
	    {
	      for (ki=0; ki<4; ++ki)
		field_info << pnts_dist[kj+ki] << " ";
	      field_info << std::endl;
	    }
	}

#ifdef DEBUG
      if (surf->dimension() == 1)
	{
	  std::ofstream of2("surf_3D.g2");
	  surf->to3D();
	  surf->writeStandardHeader(of2);
	  surf->write(of2);
	}
#endif
#ifdef DEBUG_EL
      std::ofstream ofel("bd_el.g2");
      LineCloud lines3 = surf->getElementBds();
      lines3.writeStandardHeader(ofel);
      lines3.write(ofel);
 #endif
    }
  return 0;
}

