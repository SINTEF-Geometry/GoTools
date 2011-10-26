//===========================================================================
//                                                                           
// File: test_normalcone.C                                                   
//                                                                           
// Created: Fri Oct 22 13:18:43 2004                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: test_normalcone.C,v 1.2 2005-04-07 11:26:39 jbt Exp $
//                                                                           
//===========================================================================

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"

using namespace std;
using namespace Go;

int main()
{
    ObjectHeader head;
    SplineSurface sf;
    cin >> head;
    ASSERT(head.classType() == SplineSurface::classType());
    cin >> sf;

    DirectionCone dcone = sf.normalCone();
    dcone.write(cout);
    cout << endl;
    dcone = sf.normalCone(SplineSurface::SMCornersFirst);
    dcone.write(cout);
    cout << endl;
}











