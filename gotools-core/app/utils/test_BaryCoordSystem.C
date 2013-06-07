/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

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
