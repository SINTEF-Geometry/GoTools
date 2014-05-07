#include "GoTools/utils/config.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/utils/Array.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;
using std::vector;

int main(int argc, char *argv[])
{
  if (argc != 3) {
    std::cout << "Usage: input data, output data " << std::endl;
    return -1;
  }
  int ki, kj;

  std::ifstream datain(argv[1]);
  char* outfile(argv[2]);

  char location[40];
  vector<double> xyz;
  vector<double> rain;
  double val;
  int nmb_loc = 0;
  char curr[40];

  datain >> location;
  while (!datain.eof())
    {
      nmb_loc++;
      for (ki=0; ki<3; ++ki)
	{
	  datain >> val;
	  xyz.push_back(val);
	}
      while (true)
	{
	  if (datain.eof())
	    break;
	  datain >> curr;
	  if (isdigit(curr[0]))
	    {
	      val = atof(curr);
	      rain.push_back(val);
	    }
	  else if (curr[0] == '-')
	    {
	      val = atof(curr);
	      rain.push_back(val);
	    }
	  else 
	    {
	      break;
	    }
	}
    }

  // Write to files
  char xyz_file[80];
  strcpy(xyz_file, outfile);
  strncat(xyz_file, ".xyz", 4);

  std::ofstream ofxyz(xyz_file);
  (void)ofxyz.precision(15);
  for (ki=0; ki<nmb_loc; ++ki)
    ofxyz << xyz[3*ki] << " " << xyz[3*ki+1] << " " << xyz[3*ki+2] << std::endl;

  int nmb_rain = (int)rain.size()/nmb_loc;
  for (kj=0; kj<nmb_rain; ++kj)
    {
      char rain_file[80];
      strcpy(rain_file, outfile);
      char tmp[2];
      sprintf(tmp, "%d", kj+1);
      strncat(rain_file, tmp, 2);
      strncat(rain_file, ".txt", 4);
      std::ofstream ofrain(rain_file);
      (void)ofrain.precision(15);
      for (ki=0; ki<nmb_loc; ++ki)
	{
	  if (rain[ki*nmb_rain+kj] >= 0)
	    {
	      ofrain << xyz[3*ki] << " " << xyz[3*ki+1] << " ";
	      ofrain << rain[ki*nmb_rain+kj] << std::endl;
	    }
	}
    }
      
}

  
	  
