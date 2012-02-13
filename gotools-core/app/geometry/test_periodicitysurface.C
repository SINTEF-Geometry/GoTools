//===========================================================================
//                                                                           
// File: test_periodicitysurface.C                                           
//                                                                           
// Created: Tue Mar  1 13:20:07 2005                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: test_periodicitysurface.C,v 1.1 2005-03-04 14:46:06 afr Exp $
//                                                                           
//===========================================================================


#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/GeometryTools.h"

using namespace Go;
using namespace std;

int main()
{
     ObjectHeader h;
     SplineSurface s;
     cin >> h >> s;

     int peru = GeometryTools::analyzePeriodicityDerivs(s, 0, s.basis_u().order()-2);
     cout << "U-periodicity is (derivative-based): " << peru << endl;
     int perku = GeometryTools::analyzePeriodicity(s, 0);
     cout << "U-periodicity is (knot and control point-based): " << perku << endl;
     int perv = GeometryTools::analyzePeriodicityDerivs(s, 1, s.basis_v().order()-2);
     cout << "V-periodicity is (derivative-based): " << perv << endl;
     int perkv = GeometryTools::analyzePeriodicity(s, 1);
     cout << "V-periodicity is (knot and control point-based): " << perkv << endl;
}
