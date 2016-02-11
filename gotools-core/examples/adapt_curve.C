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

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/creators/LiftCurve.h"
#include "GoTools/creators/AdaptCurve.h"
#include <fstream>

using std::cout;
using std::endl;
using std::string;
using std::cerr;
using std::ifstream;
using std::ofstream;
using namespace Go;

//===========================================================================
//                                                                           
// File: adapt_curve.C                                                  
//                                                                           
/// Description:
///  
/// This program demonstrates the use of the class AdaptCurve.
/// The class can generate a B-spline curve that approximates
/// an evaluator based curve to satisfy a given accuracy.
///
/// The program reads a 2D parameter curve and a surface from file, and
/// lifts the curve onto the surface to make an evaluator based curve.
/// This evaluator based curve are then used as input to the AdaptCurve
/// constructor.
///
/// Input/Output:
///
/// Input from the command line is the accuracy. (The maximum allowed distance
/// from the generated curve to the input curve.)
/// The file names of the input curve and surface are hardcoded to
/// "uproj_param_curve.g2" and "surface.g2".
///
/// Output is a file in Go-format for plotting the curve.
/// The file name is hard-coded to "adapt_lift_curve.g2". Together with
/// 'uproj_space_curve.g2 it can be viewed by program 'goview'.
/// If debug data are written in AdaptCurve, look at 'crv.out.g2' as well.
//   
//===========================================================================

int main(int argc, char** argv)
{
    double epsge;
    string inp_curve_filename("data/uproj_param_curve.g2");
    string inp_surf_filename("data/surface.g2");   

    // Get geometric tolerance from the argument list or use a default value.
    if (argc != 2)
    {
	epsge = 0.5;
	cout << "\nNo tolerance from command line. Using default value" << endl;
    } else {
	epsge = atof(argv[1]);
    }
    cout << "\nRunning program " << argv[0] << " with tolerance= " << epsge
	 << "\nParameter curve file name= " << inp_curve_filename.c_str()
	 << ".  Surface file name= " << inp_surf_filename.c_str() << endl;

    // Make an EvalCurve by lifting a parameter curve onto the spline surface.
    // Read parameter curve file
    ifstream cfile(inp_curve_filename.c_str());
    if (!cfile) {
	cerr << "\nFile error. Could not open file: " << inp_curve_filename.c_str()
	     << endl;
	return 1;
    }
    shared_ptr<SplineCurve> spline_curve(new SplineCurve);
    ObjectHeader header;
    cfile >> header;
    if (!(header.classType() == SplineCurve::classType())) {
	THROW("Object type is NOT SplineCurve.");
    }
    cfile >> (*spline_curve);
    cfile.close();

    Point pnt2d(2);
    spline_curve->point(pnt2d, spline_curve->startparam()); 
    cout << "\nSplineCurve:  Dimension= " << spline_curve->dimension()
	 << "\nStart.  Param= " << spline_curve->startparam() << "  Point= "
	 << pnt2d << endl;
    spline_curve->point(pnt2d, spline_curve->endparam());    
    cout << "End.  Param= " << spline_curve->endparam() << "  Point= "
	 << pnt2d << endl;
    cout << "Bounding box =   " << spline_curve->boundingBox() << endl;


    // Read surface file
    ifstream sfile(inp_surf_filename.c_str());
    if (!sfile) {
	cerr << "\nFile error. Could not open file: " << inp_surf_filename.c_str()
	     << endl;
	return 1;
    }
    shared_ptr<SplineSurface> surf(new SplineSurface);
    sfile >> header;
    if (!(header.classType() == SplineSurface::classType())) {
	THROW("Object type is NOT SplineSurface.");
    }
    sfile >> (*surf);
    sfile.close();

    cout << "\nSurface: Dimension= " << surf->dimension() << endl;
    cout << "Parameter u from "
	 << surf->startparam_u() << " to " << surf->endparam_u() << endl;
    cout << "Parameter v from "
	 << surf->startparam_v() << " to " << surf->endparam_v() << endl;
    cout << "Bounding box    = "; 
    cout << surf->boundingBox() << endl;

    // Create the evaluator based curve by lifting the parameter curve onto the
    // spline surface. 
    shared_ptr<ParamCurve> pcurve = spline_curve;
    shared_ptr<ParamSurface> psurf = surf;
    LiftCurve* eval_curve(new LiftCurve(pcurve, psurf, epsge));
    cout << "\nEvalCurve: Dimension= " << eval_curve->dim()
	 << "\nStart  Param= " <<  eval_curve->start()
	 << "  Point= " << eval_curve->eval(eval_curve->start()) << endl;
    cout << "End=   " << eval_curve->end() << 
	"  Point= " << eval_curve->eval(eval_curve->end()) << endl;

 
    // Create the AdaptCurve object. The simplest constructor takes only the
    // curve to adapt to and the geometric tolerance as arguments.
    // Other constructors can take arguments to specify an initial spline curve.
    AdaptCurve adapt_curve(eval_curve, epsge);
    // Perform the approximation
    adapt_curve.approximate(10); // Maximum ten iterations
    // Get the generated SplineCurve
    double maxdist;   // Maximum distance between the evaluator based curve and
                      // the generated B-spline curve.
    double avdist;    // Average distance between the evaluator based curve and
                      // the generated B-spline curve.
    shared_ptr<SplineCurve> adapt_spline_curve =
	adapt_curve.getAdaptCurve(maxdist, avdist, 10); // Max. ten iterations
    Point pnt(adapt_spline_curve->dimension());
    adapt_spline_curve->point(pnt, adapt_spline_curve->startparam()); 
    cout << "\nAdaptCurve:  Dimension= " << adapt_spline_curve->dimension()
	 << "\nStart.  Param= " << adapt_spline_curve->startparam()
	 << "  Point= " << pnt << endl;
    adapt_spline_curve->point(pnt, adapt_spline_curve->endparam());    
    cout << "End.  Param= " << adapt_spline_curve->endparam()
	 << "  Point= " << pnt << endl;
    cout << "Bounding box =   " << adapt_spline_curve->boundingBox() << endl; 
    cout << "\nMaximum distance between input and output curves= "
	 << maxdist << endl;
    cout << "Average distance between  input and output curves= "
	 << avdist << endl;

    // Write spline curve to file.
    ofstream fout("adapt_lift_curve.g2");
    fout << "\n100 1 0 4 255 0 0  255" << endl; // Red colour.
    adapt_spline_curve->write(fout);
    fout.close();
    // cout << "\nOpen the file 'adapt_lift_curve.g2' and 'surface.g2'"
    //      << " in 'goview' to look at the results." << endl;
    // cout << "If they are present, you may also look at the files : "
    //      << "'uproj_space_curve.g2' and  a debug file 'crv_out.g2'\n" << endl;
    delete eval_curve;

    return 0;
}
