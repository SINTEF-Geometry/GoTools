//===========================================================================
//                                                                           
// File: test_periodicity.C                                                  
//                                                                           
// Created: Thu Feb 17 09:37:08 2005                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: test_periodicity.C,v 1.2 2005-02-21 12:48:09 afr Exp $
//                                                                           
//===========================================================================


#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/GeometryTools.h"

using namespace Go;
using namespace std;

int main()
{
     ObjectHeader h;
     SplineCurve c;
     cin >> h >> c;

     int per = GeometryTools::analyzePeriodicityDerivs(c, c.basis().order()-2);
     cout << "Periodicity is (derivative-based): " << per << endl;
     int perk = GeometryTools::analyzePeriodicity(c);
     cout << "Periodicity is (knot and control point-based): " << perk << endl;
}
