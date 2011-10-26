//==========================================================================
//                                                                          
// File: test_BaryCoordSystem.C                                              
//                                                                          
// Created: Mon Dec  2 18:28:17 2002                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: test_BaryCoordSystem.C,v 1.6 2006-02-06 13:26:50 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include "GoTools/utils/BaryCoordSystem.h"
#include <iostream>


using namespace std;
using namespace Go;


int main()
{
    // Defining coord system in 2D
    Vector2D corn_tri[3];
    corn_tri[0] = Vector2D(0.0, 0.0);
    corn_tri[1] = Vector2D(1.0, 0.0);
    corn_tri[2] = Vector2D(0.5, 1.0);
    BaryCoordSystem2D bc2d(corn_tri);
    cout << "bc2d:\n" << bc2d << endl;

    // Defining coord system in 3D
    Vector3D corn_tetra[4];
    corn_tetra[0] = Vector3D(0.0, 0.0, 0.0);
    corn_tetra[1] = Vector3D(3.0, 0.0, 0.0);
    corn_tetra[2] = Vector3D(0.0, 3.0, 0.0);
    corn_tetra[3] = Vector3D(0.0, 0.0, 3.0);
    BaryCoordSystem3D bc3d(corn_tetra);
    cout << "bc3d:\n" << bc3d << endl;

    // Testing 2D
    cout << "*** Testing BaryCoordSystem2D ***" << endl;
    Vector2D cart2d(0.17, 0.29);
    Vector3D bary2d(0.3, 0.3, 0.4);

    Vector3D b2 = bc2d.cartToBary(cart2d);
    Vector2D c2 = bc2d.baryToCart(b2);

    Vector2D c3 = bc2d.baryToCart(bary2d);
    Vector3D b3 = bc2d.cartToBary(c3);

    cout << cart2d << endl << "->" << endl << b2 << endl
	 << "->" << endl << c2 << endl;
    cout << "---" << endl;
    cout << bary2d << endl << "->" << endl << c3 << endl
	 << "->" << endl << b3 << endl;

    // Testing 3D
    cout << "*** Testing BaryCoordSystem3D ***" << endl;
    Vector3D cart3d(0.17, 0.29, 0.08);
    Vector4D bary3d(0.1, 0.3, 0.4, 0.2);

    Vector4D bb2 = bc3d.cartToBary(cart3d);
    Vector3D cc2 = bc3d.baryToCart(bb2);

    Vector3D cc3 = bc3d.baryToCart(bary3d);
    Vector4D bb3 = bc3d.cartToBary(cc3);

    cout << cart3d << endl << "->" << endl << bb2 << endl
	 << "->" << endl << cc2 << endl;
    cout << "---" << endl;
    cout << bary3d << endl << "->" << endl << cc3 << endl
	 << "->" << endl << bb3 << endl;

    return 0;
}
