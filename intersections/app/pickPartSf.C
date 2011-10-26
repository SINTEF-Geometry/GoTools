#include <fstream>
#include <iomanip>
#include <memory>
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineSurface.h"

using namespace std;
using namespace Go;

int main(int argc, char** argv)
{
 
  if (argc != 3) {
    cout << "Usage: pickPartSf infile outfile"
	 << endl;
    return 0;
  }


  ObjectHeader header;

  // Read the first curve from file
  ifstream input1(argv[1]);
  if (input1.bad()) {
    cerr << "File #1 error (no file or corrupt file specified)."
	 << std::endl;
    return 1;
  }
  header.read(input1);
  std::shared_ptr<SplineSurface> surf1(new SplineSurface());
  surf1->read(input1);
  input1.close();

  double mima[4];
  std::cout << "Give parameter domain, umin umax vmin vmax :" << std::endl;
  std::cin >> mima[0];
  std::cin >> mima[1];
  std::cin >> mima[2];
  std::cin >> mima[3];

  std::shared_ptr<SplineSurface> surf2(surf1->subSurface(mima[0], mima[2],
							   mima[1], mima[3]));
  ofstream output(argv[2]);
  surf2->writeStandardHeader(output);
  surf2->write(output);
}


 
