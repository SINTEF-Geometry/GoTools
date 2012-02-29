#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineCurve.h"
#include <fstream>

using namespace std;
using namespace Go;


//===========================================================================
//                                                                           
// File: append_curve.C                                                   
//                                                                           
/// Description:
///
/// This program demonstrates the use of the function
/// SplineCurve::append_curve(void appendCurve(ParamCurve* other_curve,
///			     int continuity, 
///			     double& dist, 
///			     bool repar=true);
/// declared in SplineCurve.h.
/// The function joins two SplineCurves by appending the start of the second
/// curve to the end of the first curve. The two curves must be of the same type.
/// NOTE! The second curve will also be changed.
///
/// Input/Output
/// From the command line : Infile1  Infile2  continuity  Outfile
/// where Infile1 and  Infile2 are the name of files with the two input curves.
/// The files 'data/spline_ellipse_segm_org.g2' and
/// 'data/interpol_curve1_free.curve.g2' can be used as input files.
/// continuity is the level of continuity we demand at the transition between the
/// two curves (can be from -1 to order()-1, but the higher the value the more
/// will the curves have to be locally modified.
/// continuity = -1 inserts a straight line between the end point of the first
/// curve and the start point of the second curve.
/// continuity = 0 moves the end point of the first curve and the start point of
/// the second curve to the midpoint between them.
/// continuity = 1 is the default value when calling the short hand version
/// void appendCurve(ParamCurve* cv, bool repar=true);
/// Outfile is the name of the file where the new curve will be written.
//
//===========================================================================

int main(int argc, char** argv)
{
    if (argc != 5) {
	cout << "\nUsage: " << argv[0]
	     << " Infile1 Infile2 continuity Outfile\n" << endl;
	exit(-1);
    }
    cout << "\nRunning program " << argv[0]
	 << "\nInfile1    = " << argv[1]
	 << "\nInfile2    = " << argv[2]
	 << "\ncontinuity = " << argv[3]
	 << "\nOutfile    = " << argv[4]
	 <<  '\n' << endl;
    // Read first curve file
    ifstream infile(argv[1]);
    if (!infile) {
	cerr << "\nFile error. Could not open file: " << argv[1] << endl;
	return 1;
    }
    SplineCurve curve, other_curve;
    ObjectHeader header;
    infile >> header >> curve;
    if (!header.classType() == SplineCurve::classType()) {
	THROW("Object type is NOT SplineCurve.");
    }
    infile.close();

    // Read second curve file
    ifstream infile2(argv[2]);
    if (!infile2) {
	cerr << "\nFile error. Could not open file: " << argv[2] << endl;
	return 1;
    }
    infile2 >> header >> other_curve;
    if (!header.classType() == SplineCurve::classType()) {
	THROW("Object type is NOT SplineCurve.");
    }
    infile2.close();

    cout << "Curve orders : " << curve.order() << " and " << other_curve.order()
	 << endl;
    // Append the start of the second curve to the end of the first curve.
    double dist = 0; // The estimated maximum distorsion after 'smoothing' of
                     // joined curve to achieve the desired continuity.
    int continuity = atoi(argv[3]);  // Continuity. From -1 to order()-1.
    bool repar = true;     // The reparametrizatin of the second curve will also
                           // be scaled as a function of position of control
                           // points close to the transition.
    curve.appendCurve(&other_curve, continuity, dist, repar);
    cout << "Curve orders : " << curve.order() << " and " << other_curve.order()
 	 << " after appending curve." << endl;
    cout << "Estimated difference between original and smooth curve: "
	 << dist << endl;

    // Write the the new curve to output file. 
    ofstream outfile(argv[4]);
    outfile << header << curve;	

    return 0;

}










