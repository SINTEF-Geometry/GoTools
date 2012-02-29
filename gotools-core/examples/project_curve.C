#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/creators/CurveCreators.h"
#include <fstream>

using std::string;
using std::ifstream;
using std::cerr;
using std::endl;
using std::ofstream;
using namespace Go;

//===========================================================================
//                                                                           
// File: project_curve.C                                                  
///                                                                           
/// Description:
///  
/// This program demonstrates the use of the function
/// void CurveCreators::projectCurve(shared_ptr<ParamCurve>& space_cv,
///                                  shared_ptr<ParamSurface>& surf,
///   			            double epsge,
///			            shared_ptr<SplineCurve>& proj_cv,
///			            shared_ptr<SplineCurve>& par_cv)
///
/// The function generates a cubic spline curve(order four) which lies on a given
/// SplineSurface and is the projection of a given space curve (in 3D space) onto
/// that surface, within a given tolerance. The given space curve should be close
/// to the surface.
/// 
/// Input/Output:
///
/// Input to this program from the command line is the accuracy. (The maximum
/// allowed distance from the generated curve to the surface.)
/// The file names of the given spline curve and the surface are hardcoded 
/// to 'const_v_paramcurve.g2' and 'surface.g2'.
/// Output is one file with a 3D space curve 'proj_space_curve.g2' and a file
/// with a 2D parameter curve 'proj_param_curve.g2'.
/// The main purpose of the function is to get the parameter curve.
///
/// If the wanted tolerance can not be achieved, the message "Knot interval too
/// small" is displayed, but the curves are generated.
///   
//===========================================================================


int main(int argc, char** argv)
{
    // Get geometric tolerance from the argument list.    
    double epsge = atof(argv[1]);

    // Read curve file
    //string inp_curve_filename("approj_curve.g2"); 
    //string inp_curve_filename("const_u_paramcurve.g2"); 
    string inp_curve_filename("const_v_paramcurve.g2"); 
    ifstream cfile(inp_curve_filename.c_str());
    if (!cfile) {
	cerr << "\nFile error. Could not open file: " << inp_curve_filename.c_str() << endl;
	return 1;
    }
    shared_ptr<SplineCurve> curve(new SplineCurve);
    ObjectHeader header;
    cfile >> header;
    if (!header.classType() == SplineCurve::classType()) {
	THROW("Object type is NOT SplineCurve.");
    }
    cfile >> (*curve);
    cfile.close();
    shared_ptr<ParamCurve> pcurve = curve;

    // Read surface file
    string inp_surf_filename("surface.g2");   
    ifstream sfile(inp_surf_filename.c_str());
    if (!sfile) {
	cerr << "\nFile error. Could not open file: " << inp_surf_filename.c_str() << endl;
	return 1;
    }
    shared_ptr<SplineSurface> surf(new SplineSurface);
    sfile >> header;
    if (!header.classType() == SplineSurface::classType()) {
	THROW("Object type is NOT SplineSurface.");
    }
    sfile >> (*surf);
    sfile.close();
    shared_ptr<ParamSurface> psurf = surf;

    // Create projection curves
    shared_ptr<SplineCurve> proj_cv;
    shared_ptr<SplineCurve> par_cv;
    CurveCreators::projectCurve(pcurve, psurf, epsge, proj_cv, par_cv);

    // Write 3D space curve to file.
    ofstream scvout("proj_space_curve.g2");
    //proj_cv->writeStandardHeader(scvout);
    scvout << "100 1 0 4 0 255 0  255" << endl; // write header. Green curve
    scvout << *proj_cv;    // write space spline curve data.
    scvout.close();

    // Write 2D parameter curve to file.
    ofstream pcvout("proj_param_curve.g2");
    par_cv->writeStandardHeader(pcvout); // write header
    pcvout << *par_cv;    // write parameter spline curve data.
    pcvout.close();

    return 0;
}
