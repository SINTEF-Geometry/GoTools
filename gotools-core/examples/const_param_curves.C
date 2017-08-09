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

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <fstream>

using namespace std;
using namespace Go;

//===========================================================================
//                                                                           
// File: const_param_curves.C                                                  
//                                                                           
/// Description:
///
/// This program demonstrates the use of the function
/// 'SplineSurface::constParamCurve(double parameter, bool pardir_is_u).
/// The declaration of the function is in 'SplineSurface.h'.
///
/// The function generates and returns a SplineCurve that represents a parameter
/// curve on the surface with a constant parameter in one parameter direction 
/// and all the parameters in the other direction.
///
/// It reads the spline surface file 'cpc_surface.g2' and generates two curves.
/// One curve with constant u parameter and one with constant v parameter.
/// The parameter values are the midpoint of the parameter range.
///
/// Output:
///
/// The program will write the output files 'const_u_paramcurve.g2' and 
/// 'const_v_paramcurve.g2'. The surface and the curves can be displayed with the
/// program 'goview'.
///
//===========================================================================


int main(int argc, char** argv)
{
    // Read the surface from a file in Go-format.
    string surf_filename("data/surface.g2");
    ifstream sfile(surf_filename.c_str());
    if (!sfile) {
	cerr << "\nFile error. Could not open file: " << surf_filename.c_str() << endl;
	return 1;
    }
    ObjectHeader head;
    SplineSurface surf;
    sfile >> head;
    if (!(head.classType() == SplineSurface::classType())) {
	THROW("Object type is NOT SplineSurface.");
    }
    sfile >> surf;
    sfile.close();

    cout << "\nProgram '" << argv[0] << "' using input file '" << surf_filename.c_str()
	 << endl;
    cout << "Parameter u from " << surf.startparam_u() << " to "
	 << surf.endparam_u() << endl;
    cout << "Parameter v from " << surf.startparam_v() << " to "
	 << surf.endparam_v() << endl;

    // Generate and return a SplineCurve that represents a constant u-parameter 
    // curve on the surface. The running parameter direction is v. 
    double upar = 0.5*(surf.startparam_u() +   // Fixed parameter
		       surf.endparam_u());
    bool pardir_is_u = false;    // The fixed parameter is in the first
    // parameter direction, thus the running parameter direction is in
    // the second parameter direction and consequently not u.
    SplineCurve* const_u_curve = surf.constParamCurve(upar, pardir_is_u);   

    // Write curve to file
    ofstream ufout("const_u_paramcurve.g2");
    //const_u_curve->writeStandardHeader(ufout);
    ufout << "100 1 0 4 0 255 0  255" << endl; // write header. Green curve
    ufout << *const_u_curve;    // write spline curve data.
    ufout.close();

    cout << "Curve with constant u-parameter " << upar 
	 << " written to file : 'const_u_paramcurve.g2'" << endl;
    cout << "Start parameter = " << const_u_curve->startparam() << "\t "
	 << "End parameter = " << const_u_curve->endparam() << endl;

    //It is the user's reponsibility to delete it when it is no longer needed.
    delete const_u_curve;


    // Generate and return a SplineCurve that represents a constant v-parameter 
    // curve on the surface. The running direction is u.
    double vpar = 0.5*(surf.startparam_v() + surf.endparam_v());
    pardir_is_u = true;    //  The fixed parameter is in the second
    // parameter direction, thus the running parameter direction is in
    // the first parameter direction and consequently u.
    SplineCurve* const_v_curve = surf.constParamCurve(vpar, pardir_is_u);   

    // Write curve to file
    ofstream vfout("const_v_paramcurve.g2");
    vfout << "100 1 0 4 255 0 0 255" << endl;  // write header. Red curve
    vfout << *const_v_curve;    // write spline curve data.
    vfout.close();

    cout << "Curve with constant v-parameter " << vpar
	 << " written to file : 'const_v_paramcurve.g2'" << endl;
    cout << "Start parameter = " << const_v_curve->startparam() << "\t "
	 << "End parameter = " << const_v_curve->endparam() << endl;

    //It is the user's reponsibility to delete it when it is no longer needed.
    delete const_v_curve;
}
