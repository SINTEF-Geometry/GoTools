#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/intersections/Identity.h"
#include "GoTools/igeslib/IGESconverter.h"
#include <fstream>

using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using namespace Go;
using std::shared_ptr;
using std::dynamic_pointer_cast;

int main(int argc, char** argv)
{
    if (argc != 4) {
	cout << "Usage: test_SfSfIdentity FileSf1 FileSf2 "
	     <<" aepsge"
	     << endl;
	return 0;
    }
    ObjectHeader header;

    // Read the first surface from file
    std::ifstream input1(argv[1]);
    if (input1.bad()) {
	cerr << "File #1 error (no file or corrupt file specified)."
	     << std::endl;
	return 1;
    }
    IGESconverter conv1;
    conv1.readgo(input1);
    vector<shared_ptr<GeomObject> > geomobj1 = conv1.getGoGeom();
    std::shared_ptr<ParamSurface> surf1;
    if (geomobj1[0]->instanceType() == Class_SplineSurface || 
	geomobj1[0]->instanceType() == Class_BoundedSurface)
	surf1 = dynamic_pointer_cast<ParamSurface, GeomObject>(geomobj1[0]);

    
    // Read the second surface from file
    ifstream input2(argv[2]);
    if (input2.bad()) {
	cerr << "File #2 error (no file or corrupt file specified)."
	     << std::endl;
	return 1;
    }

    IGESconverter conv2;
    conv2.readgo(input2);
    vector<shared_ptr<GeomObject> > geomobj2 = conv2.getGoGeom();
    std::shared_ptr<ParamSurface> surf2;
    if (geomobj2[0]->instanceType() == Class_SplineSurface || 
	geomobj2[0]->instanceType() == Class_BoundedSurface)
	surf2 = dynamic_pointer_cast<ParamSurface, GeomObject>(geomobj2[0]);

    double aepsge;
    aepsge = atof(argv[3]);

    Identity ident;
    int coincident = ident.identicalSfs(surf1, surf2, aepsge);

    std::cout << "Identity: " << coincident << std::endl;
}
