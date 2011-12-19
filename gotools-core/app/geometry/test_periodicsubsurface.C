//===========================================================================
//                                                                           
// File: test_periodicity.C                                                  
//                                                                           
// Created: Thu Feb 17 09:37:08 2005                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: test_periodicsubsurface.C,v 1.1 2005-03-04 14:46:06 afr Exp $
//                                                                           
//===========================================================================


#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"

using namespace Go;
using namespace std;

int main()
{
     ObjectHeader h;
     SplineSurface s;
     cin >> h >> s;

     double a = s.startparam_u();
     double b = s.endparam_u();
     double tmin = a + 0.7*(b-a);
     double tmax = tmin + 0.6*(b-a);
     shared_ptr<SplineSurface>
	 ss(s.subSurface(tmin, s.startparam_v()+0.1, tmax, s.endparam_v()-0.1));
     cout << h << (*ss);
}
