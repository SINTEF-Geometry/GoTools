#include "GoTools/utils/config.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/utils/Array.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSurfApprox.h"
#include <iostream>
#include <fstream>
#include <string.h>

#define DEBUG

using namespace Go;
using std::vector;
using std::pair;

int main(int argc, char *argv[])
{
  if (argc != 4)
    {
      std::cout << "Input parameters: input file (.csv), input result file, output file " << std::endl;
      return -1;
    }

 // Define parameters
  std::ifstream in1(argv[1]);
  std::ifstream in2(argv[2]);
  std::ofstream out(argv[3]);
  double aeps = 0.05;

 // Read rain data
  char location[40];
  vector<double> xyz;
  vector<double> rain;
  double val;
  int nmb_loc = 0;
  char curr[40];

  int ki, kj, kr;
  in1 >> location;
  while (!in1.eof())
    {
      nmb_loc++;
      for (ki=0; ki<3; ++ki)
	{
	  in1 >> val;
	  xyz.push_back(val);
	}
      while (true)
	{
	  if (in1.eof())
	    break;
	  in1 >> curr;
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
  int nmb_rain = (int)rain.size()/nmb_loc;

  // Read result
  int nmb_loc2 = 0;
  double tmp;
  vector<double> result;
  while (!in2.eof())
    {
      for (ki=0; ki<nmb_rain; ++ki)
	{
	  in2 >> tmp;
	  result.push_back(tmp);
	}
      nmb_loc2++;
    }
  nmb_loc2 = std::min(nmb_loc2, nmb_loc);

  // Compare
  int first = nmb_loc - nmb_loc2;

  double maxdiff = 0.0;
  double avdiff = 0.0;
  double avdiffsgn = 0.0;
  double diff;
  int nmb = 0;
  int nmbaeps = 0;
  for (ki=first, kj=0; ki<nmb_loc; ++ki, ++kj)
    {
      for (kr=0; kr<nmb_rain; ++kr)
	{
	  if (rain[ki*nmb_rain+kr] >= 0)
	    {
	      diff = result[kj*nmb_rain+kr] - rain[ki*nmb_rain+kr];
	      maxdiff = std::max(maxdiff, fabs(diff));
	      avdiff += fabs(diff);
	      avdiffsgn += diff;
	      out << diff << " ";
	      nmb++;
	      if (fabs(diff) > aeps)
		nmbaeps++;
	    }
	  else
	    out << "null ";
	}
      out << std::endl;
    }
  avdiff /= (double)nmb;
  avdiffsgn /= (double)nmb;
  out << "Max diff: " << maxdiff << std::endl;
  out << "Average diff: " << avdiff << std::endl;
  out << "Signed average diff: " << avdiffsgn << std::endl;
  out << "Number of diffs > 0.05: " << nmbaeps << std::endl;

  std::cout << "Max diff: " << maxdiff << std::endl;
  std::cout << "Average diff: " << avdiff << std::endl;
  std::cout << "Signed average diff: " << avdiffsgn << std::endl;
  std::cout << "Number of diffs > 0.05: " << nmbaeps << std::endl;
}



	  
  
  
