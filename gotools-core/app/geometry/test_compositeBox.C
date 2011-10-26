//===========================================================================
//                                                                           
// File: test_compositeBox.C                                                 
//                                                                           
// Created: Wed Dec  8 13:22:37 2004                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: test_compositeBox.C,v 1.1 2004-12-08 13:36:44 afr Exp $
//                                                                           
//===========================================================================

#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"

using namespace std;
using namespace Go;


int main()
{
    ObjectHeader head;
    cin >> head;
    if (head.classType() == SplineCurve::classType()) {
	SplineCurve curve;
	cin >> curve;
	CompositeBox b = curve.compositeBox();
	b.write(cout);
	BoundingBox b2 = curve.boundingBox();
	b2.write(cout);
    } else if (head.classType() == SplineSurface::classType()) {
	SplineSurface surf;
	cin >> surf;
	CompositeBox b = surf.compositeBox();
	b.write(cout);
	BoundingBox b2 = surf.boundingBox();
	b2.write(cout);
    } else {
	MESSAGE("Unknown object.");
    }
}
