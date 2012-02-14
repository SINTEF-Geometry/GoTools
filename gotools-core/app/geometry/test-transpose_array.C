//===========================================================================
//                                                                           
// File: test-transpose_array.C                                              
//                                                                           
// Created: Tue Jul  2 14:49:42 2002                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: test-transpose_array.C,v 1.2 2003-05-08 15:37:15 afr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/geometry/SplineUtils.h"
#include <iostream>

using namespace std;
using namespace Go;

int main()
{
    double v[] = { 0, 1, 2, 3, 4, 5 };
    SplineUtils::transpose_array(1, 2, 3, v);
    cout << "We want: 0 3 1 4 2 5" << endl;
    cout << v[0] << ' ' << v[1]<< ' ' << v[2]<< ' ' << v[3]<< ' ' << v[4]<< ' ' << v[5] << endl;
}
