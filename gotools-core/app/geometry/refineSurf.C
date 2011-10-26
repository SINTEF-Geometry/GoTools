#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <fstream>
#include <vector>

using namespace Go;
using namespace std;

int main(int argc, char** argv)
{
  if (argc < 3) {
      cerr << "Usage: " << argv[0]
	   << " inputfile outputfile [max_coefs_u max_coefs_v]" << endl;
      return 1;
  }

  ifstream in(argv[1]);
  ofstream out(argv[2]);

  if (!in || !out) {
    cout << "Bad file(s) or filename(s)." << endl;
    return 1;
  }

  ObjectHeader oh;
  SplineSurface sf;

  in >> oh >> sf;


  int m = sf.numCoefs_v() - sf.order_v() + 1;
  int n = sf.numCoefs_u() - sf.order_u() + 1;
  if (argc >= 5) {
      // Note the weird order (v then u)
      m = min(atoi(argv[4])-sf.numCoefs_v(), m);
      n = min(atoi(argv[3])-sf.numCoefs_u(), n);
  }
  int i;
  vector<double> newknots_v;
  vector<double> newknots_u;
  for (i = 0; i < m; ++i) {
    vector<double>::const_iterator it = sf.basis_v().begin();
    double newknot = 0.5*it[sf.order_v()+i-1] + 0.5*it[sf.order_v()+i];
    newknots_v.push_back(newknot);
  }
  for (i = 0; i < n; ++i) {
    vector<double>::const_iterator it = sf.basis_u().begin();
    double newknot = 0.5*it[sf.order_u()+i-1] + 0.5*it[sf.order_u()+i];
    newknots_u.push_back(newknot);
  }

  sf.insertKnot_v(newknots_v);
  sf.insertKnot_u(newknots_u);

  out << oh << sf;
  return 0;
}


