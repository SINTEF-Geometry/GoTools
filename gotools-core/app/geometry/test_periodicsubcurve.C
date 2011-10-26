//===========================================================================
//                                                                           
// File: test_periodicity.C                                                  
//                                                                           
// Created: Thu Feb 17 09:37:08 2005                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: test_periodicsubcurve.C,v 1.1 2005-03-04 14:46:06 afr Exp $
//                                                                           
//===========================================================================


#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <memory>

using namespace Go;
using namespace std;

int main()
{
     ObjectHeader h;
     SplineCurve c;
     cin >> h >> c;

     double a = c.startparam();
     double b = c.endparam();
     double tmin = a + 0.7*(b-a);
     double tmax = tmin + 0.6*(b-a);
     std::shared_ptr<SplineCurve> sc(c.subCurve(tmin, tmax));
     cout << h << (*sc);
}
