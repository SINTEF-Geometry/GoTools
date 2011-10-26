#include "GoTools/parametrization/PrTriangulation_OP.h"
#include "GoTools/parametrization/PrPlanarGraph_OP.h"
#include "GoTools/parametrization/PrParametrizeBdy.h"
#include "GoTools/parametrization/PrPrmShpPres.h"
#include <fstream>
#include <memory>
using std::shared_ptr;
using std::cout;
using std::endl;

int main()
{
  cout << "Here is an explicit triangulation:\n";
  cout << "\n";
  cout << "   2-------3 \n";
  cout << "   |\\      | \n";
  cout << "   |\\\\     | \n";
  cout << "   | |\\    | \n";
  cout << "   | | \\   | \n";
  cout << "   | |  \\  | \n";
  cout << "   | 4-  \\ | \n";
  cout << "   |/  \\--\\| \n";
  cout << "   0-------1 \n";

  int numpnts=5, numtrs=4;
  vector<double> points(3*numpnts);
  vector<int> triangles(3*numtrs);
  points[0] = 0.0; points[1] = 0.0; points[2] = 0.0; 
  points[3] = 1.0; points[4] = 0.0; points[5] = 0.0; 
  points[6] = 0.0; points[7] = 1.0; points[8] = 0.0; 
  points[9] = 1.0; points[10] = 1.0; points[11] = 1.0; 
  points[12] = 0.25; points[13] = 0.25; points[14] = 0.5; 
  triangles[0] = 0; triangles[1] = 1; triangles[2] = 4;
  triangles[3] = 4; triangles[4] = 2; triangles[5] = 0;
  triangles[6] = 1; triangles[7] = 3; triangles[8] = 2;
  triangles[9] = 1; triangles[10] = 2; triangles[11] = 4;

  
  shared_ptr<PrTriangulation_OP>
      pr_triang(new PrTriangulation_OP(&points[0],
				       numpnts,
				       &triangles[0],
				       numtrs));

  cout << "The data of pr_triang is \n";
  pr_triang->print(cout);
  cout << "The information concerning pr_triang is \n";
  pr_triang->printInfo(cout);

  cout << "The nodes of pr_triang are \n";
  pr_triang->printXYZNodes(cout);
  cout << "The edges of pr_triang are \n";
  pr_triang->printXYZEdges(cout);
  cout << "The triangles of pr_triang are \n";
  pr_triang->printXYZFaces(cout);



  //----------------------- Added testing code from Atgeirr ----------------


  PrParametrizeBdy bdy;
  bdy.attach(pr_triang);
  bdy.setParamKind(PrCENTRIPETAL);

  cout << "Parametrizing boundary..." << endl;
  bdy.parametrize();

  if (pr_triang->getNumNodes() - pr_triang->findNumBdyNodes() > 0) {
      PrPrmShpPres interior;
      interior.attach(pr_triang);
      interior.setStartVectorKind(PrBARYCENTRE);
      interior.setBiCGTolerance(1.0e-6);
      cout << "Parametrizing interior..." << endl;
      interior.parametrize();
  }

  cout << "Saving parametrization to file..." << endl;
  std::ofstream pout("edges_out");
  pr_triang->printUVEdges(pout);

  cout << "Quitting..." << endl;



  return 0;
}
