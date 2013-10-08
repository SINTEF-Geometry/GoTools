#include "GoTools/geometry/PointCloud.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;
using std::vector;

int main(int argc, char *argv[])
{
  if (argc != 3) {
    std::cout << "Usage: infile outfile.g2 " << std::endl;
    return -1;
  }

  std::ifstream filein(argv[1]);
  std::ofstream fileout(argv[2]); 
  int in, im;
  double x1, y1;
  double del;
  double no_data;
  filein >> in;
  filein >> im;
  filein >> x1;
  filein >> y1;
  filein >> del;
  filein >> no_data;

  int ki, kj;
  double x, y, z;
  vector<double> data;
  for (kj=0, y=y1; kj<im; ++kj, y+=del)
    for (ki=0, x=x1; ki<in; ++ki, x+=del)
      {
	filein >> z;
	if (z == no_data)
	  continue;
	data.push_back(x);
	data.push_back(y);
	data.push_back(z);
      }
  
  PointCloud3D points(data.begin(), data.size()/3);
  points.writeStandardHeader(fileout);
  points.write(fileout);
}

	
