//===========================================================================
//                                                                           
// File: test_rotatedbox.C                                                   
//                                                                           
// Created: Wed Jan 12 16:03:11 2005                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: test_rotatedbox.C,v 1.1 2005-01-12 15:22:10 afr Exp $
//                                                                           
//===========================================================================

#include "GoTools/utils/RotatedBox.h"

using namespace std;

int main()
{
    double sq2h = std::sqrt(2.0)/2.0;
    double pt[] = {0,0,  1, 1.1,   2.2, 2 };
    Go::Point axis[1];
    axis[0].resize(2);
    axis[1].resize(2);
    axis[0][0] = sq2h;
    axis[0][1] = sq2h;
    Go::RotatedBox rb(pt, 2, 3, 1, axis);
    Go::Point p(3);
    p[0] = 1.5;
    p[1] = 0.6;
    cout << "RotBox contains: " << rb.containsPoint(p) << endl;

    Go::CompositeBox cb(pt, 2, 3, 1);
    cout << "CompBox contains: " << cb.containsPoint(p) << endl;
}
