//===========================================================================
//                                                                           
// File: make_grid.C                                                         
//                                                                           
// Created: Thu Feb  1 14:58:27 2007                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: make_grid.C,v 1.1 2007-03-08 10:45:18 afr Exp $
//                                                                           
//===========================================================================



#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/RectGrid.h"

using namespace std;
using namespace Go;

int main(int argc, char** argv)
{
    if (argc < 3) {
	cerr << "Usage: " << argv[0] << " u_res v_res" << endl;
	return 1;
    }
    int ures = atoi(argv[1]);
    int vres = atoi(argv[2]);

    ObjectHeader head;
    cin >> head;
    ASSERT(head.classType() == SplineSurface::classType());
    SplineSurface sf;
    cin >> sf;
    vector<double> points;
    vector<double> param_u;
    vector<double> param_v;
    sf.gridEvaluator(ures, vres, points, param_u, param_v);
    RectGrid grid(ures, vres, sf.dimension(), &points[0]);
    grid.writeStandardHeader(cout);
    grid.write(cout);
}
