//===========================================================================
//                                                                           
// File: surf2cloud.C                                                        
//                                                                           
// Created: Wed May 15 10:23:18 2002                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: surf2cloud.C,v 1.2 2006-04-19 09:27:33 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/ObjectHeader.h"

using namespace std;
using namespace Go;

int main()
{
    ObjectHeader header;
    SplineSurface surf;
    cin >> header >> surf;

    PointCloud<3> cloud(surf.coefs_begin(),
		     surf.numCoefs_u()*surf.numCoefs_v());
    cloud.writeStandardHeader(cout);
    cout << cloud;
}
