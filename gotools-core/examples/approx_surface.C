//===========================================================================
//                                                                           
// File: approx_surface.C                                                       
//                                                                           
// Description:
//
// This program demonstrates the use of the class ApproxSurf.
// The class can generate a new tensor product B-spline surface with four boundary
// curves that approximates a set of parametrized points for a given accuracy, or
// modify an old surface by a set of parametrized points.
//
// Input to this program from the command line is the accuracy (the maximum
// allowed distance from one of the points to the surface). 
//
// The input points and their parameter values, the boundary curves and other
// arguments to ApproxSurf's constructor are made by this program.
// All the points will not necessarily be within the wanted distance from the surface.
//
// Output is a file in Go-format for plot of the points and the surfaces.
// The file name is "approx_surface.g2".
//
//===========================================================================

#include "GoTools/creators/ApproxSurf.h"
#include "GoTools/creators/ApproxCurve.h"
#include "GoTools/utils/Point.h"
#include <vector>
#include <fstream>

using std::cout;
using std::endl;
using std::vector;
using std::ofstream;
using std::endl;
using namespace Go;
using std::shared_ptr;

// Help function. Make a vector of n equidistant values from start to end.
void fill(double start, double end, int n, vector<double>& values)
{
    values.resize(n);
    double step = (end-start)/(n-1);
    for(int i=0; i<n; ++i) {
	values[i] = start + i*step;
    }
}


int main(int argc, char** argv)
{
    if (argc != 2) {
      cout << "\nUsage: " << argv[0]
	   << " tolerance\n" << endl;
      exit(-1);
    }

    // Get geometric tolerance from the argument list.    
    double aepsge = atof(argv[1]);
    const int dim   = 3;  // Space dimension.
    cout << "\nRunning program " << argv[0] << " with tolerance = " << aepsge
	 << endl;

    // Define the four boundary curves of the surface.
    // The endpoints of the curves must be connected so that they form a loop.
    std::vector<std::shared_ptr<SplineCurve> > surf_boundary_curves(4);
    double start_param = 0.0;
    double end_param   = 1.0;
    int num_crv_pts = 3;   // 
    std::vector<double> crv_param(num_crv_pts);  // Parameter values
    fill(start_param, end_param, num_crv_pts, crv_param); 

    // Surface corner points. 
    const int num_bnd_crvs = 4;
    double bnd_corner_x[] = {1, 2, 2, 1};
    double bnd_corner_y[] = {3, 3, 4, 4};
    double bnd_corner_z[] = {5, 6, 5, 6};
    
    std::vector<double> crv_points(dim*num_crv_pts);  // Curve points
 
    for (int i=0; i<num_bnd_crvs; ++i) {
	vector<double> x(num_crv_pts), y(num_crv_pts), z(num_crv_pts);
	fill(bnd_corner_x[i], bnd_corner_x[(i+1)%4], num_crv_pts, x);  
	fill(bnd_corner_y[i], bnd_corner_y[(i+1)%4], num_crv_pts, y);  
	fill(bnd_corner_z[i], bnd_corner_z[(i+1)%4], num_crv_pts, z);
	for (int j=0, p=-1; j<num_crv_pts; ++j) {
	    crv_points[++p] = x[j];
	    crv_points[++p] = y[j];
	    crv_points[++p] = z[j];
	}

	double maxdist;   // Maximum distance between the generated curve and
	// the data points.
	double avdist;    // Average distance between the generated curve and
	// the data points.
	const int max_iter = 5; // Maximum number of iterations to use.
	ApproxCurve approx_curve(crv_points, crv_param, dim, aepsge);
	surf_boundary_curves[i] =    // Get the surface boundary spline curve.  
	    approx_curve.getApproxCurve(maxdist, avdist, max_iter);    
    }

    // Define approximation points.
    int numpt = 1;  // Number of input points.
    // Points. Consecutively in "xyzxyzxyz-fashion".
    // Parameters. Consecutively in "uvuvuvuv-fashion".
    vector<double> points(numpt * dim);
    vector<double> parvals(numpt * 2);
    for (int i = 0; i < numpt; ++i) {
	points[i*dim + 0] = 1.5;  // x
	points[i*dim + 1] = 3.5;  // y
	points[i*dim + 2] = 5.5;  // z
	parvals[2*i]   = 0.5;  // u
	parvals[2*i+1] = 0.5;  // v
    }

    // Parametric domain for the surface.
    double param_domain[] = {0,1,0,1};   // u_min, u_max, v_min, v_max
    int constdir = 0;  // Reparameterize the points in both directions.
    bool repar = true;
 
    // Approximate a surface through the points. Constructor 1.
    ApproxSurf approx_surf1(surf_boundary_curves, points, parvals, param_domain,
    			    dim, aepsge, constdir, repar);

    // Get the spline surface.
    double maxdist, avdist;
    int nmb_out_eps;
    const int max_iter = 5;  // Maximum number of iterations to use.
    const int keep_init = 0; // For each iteration step, the surface is 
    // updated to approximate the data points, the accuracy is checked, 
    // the spline space of the surface is refined, and the modified 
    // surface is used as input for the new iteration step. If keep_init
    // is set to a number larger than 0, for the first keep_init iteration
    // steps, a refined version the initial surface will be used as input
    // to parameter iteration and approximation.
    shared_ptr<SplineSurface> spline_surf1 = 
	approx_surf1.getApproxSurf(maxdist, avdist, nmb_out_eps, max_iter,
				   keep_init);
    cout << "\nSurface 1\n";
    cout << "Maximum distance between surface 1 and the data points= "
	 << maxdist << endl;
    cout << "Average distance between surface 1 and the data points= "
	 << avdist << endl;
    cout << "Number of data points with distance greater than tolerance= " 
	 << nmb_out_eps << endl;
    cout << "Parameter u from " << spline_surf1->startparam_u() << " to "
	 << spline_surf1->endparam_u() << endl;
    cout << "Parameter v from " << spline_surf1->startparam_v() << " to "
	 << spline_surf1->endparam_v() << endl;
    cout << "Number of control points:  u= " << spline_surf1->numCoefs_u()
	 << "   v= " <<  spline_surf1->numCoefs_v()<< endl;

    // Evaluate the position of the surface at some parameter values u and v.
    // These points will get new z-values, and will be used to create a
    // modified surface.
    int num_extrapoints = 4;
    vector<Point> extrapoints(num_extrapoints);
    vector<double> u_params(num_extrapoints), v_params(num_extrapoints); 
    u_params[0] = 0.25;  v_params[0] = 0.25; 
    u_params[1] = 0.75;  v_params[1] = 0.25; 
    u_params[2] = 0.25;  v_params[2] = 0.75; 
    u_params[3] = 0.75;  v_params[3] = 0.75;
    for (int i=0; i<num_extrapoints; ++i) {
	Point pnt;
	spline_surf1->point(pnt, u_params[i], v_params[i]);
	cout << "Parameters : u=" << u_params[i] << "  v=" <<  v_params[i]
	     << "\t Position= " << pnt << endl;
	extrapoints[i] = pnt;
    }
 
    // Write surface to file
    ofstream fout("approx_surface.g2");
    // Class_SplineSurf=200 MAJOR_VERSION=1 MINOR_VERSION=0 auxillary_data=0
    spline_surf1->writeStandardHeader(fout); // write header.
    fout << *spline_surf1;    // write spline spline_surf1 data.

    // Write input points to file
    // Class_PointCloud=400 MAJOR_VERSION=1 MINOR_VERSION=0 auxillary data=4
    // The four auxillary data values defines the colour (r g b alpha)
    fout << "400 1 0 4 255 0 0 255" << endl; // Header.
    fout << numpt << endl;
    for (int i = 0; i < numpt; ++i) {
	int ip = i*dim;
	Go::Point inp_point(points[ip], points[ip+1], points[ip+2]);
	fout << inp_point << endl;  // write input point coordinates.
    }


    // -------------------------------------------------------------------------

    // Modify a surface by adding points. Constructor 2.
    // Points.     Consecutively in "xyzxyzxyz-fashion".
    // Parameters. Consecutively in "uvuvuvuv-fashion".
    // Define modifying points and their parameters.
    points.resize(num_extrapoints * dim);
    parvals.resize(num_extrapoints * 2);
    for (int i = 0; i < num_extrapoints; ++i) {
	points[i*dim + 0] = extrapoints[i][0];         // x
	points[i*dim + 1] = extrapoints[i][1];         // y
	points[i*dim + 2] = extrapoints[i][2] + 0.25;  // z
	parvals[2*i]   = u_params[i];  // u
	parvals[2*i+1] = v_params[i];  // v
    }

    constdir = 0; // Reparameterize the points in both their u and v parameters.
                  // Only v:constdir=1, Only u:constdir=2.
    bool approx_orig = false;   // Indicates if an approximation of the given
    // input surface should be included in the functional. Defining this
    // parameter as true, may stabilize the computations if the data points
    // are badly distributed. It will, on the other hand, lead to a poorer
    // approximation.
    bool close_belt  = true;    // True if only coeffiecients close to the 
                                // sampling points should be modified.
    int nmb_stabil = 0;         // The last nmb_stabil data points are
    // included in the approximation, but not in the error computation. Thus,
    // a larger error for some of these points will not alter the approximation.
    repar = true;               // Indicates if a reparameterization of the
    // data points is to be performed at each step
    ApproxSurf approx_surf2(spline_surf1, points, parvals, dim, aepsge,
			    constdir, approx_orig, close_belt, nmb_stabil, repar);
 
   // Get the new spline surface.
    shared_ptr<SplineSurface> spline_surf2 =
	approx_surf2.getApproxSurf(maxdist, avdist, nmb_out_eps, max_iter,
				   keep_init);
    cout << "\nSurface 2\n";
    cout << "Maximum distance between surface 2 and the data points= "
	 << maxdist << endl;
    cout << "Average distance between surface 2 and the data points= "
	 << avdist << endl;
    cout << "Number of data points with distance greater than tolerance= " 
	 << nmb_out_eps << endl;
    cout << "Parameter u from " << spline_surf2->startparam_u() << " to "
	 << spline_surf2->endparam_u() << endl;
    cout << "Parameter v from " << spline_surf2->startparam_v() << " to "
	 << spline_surf2->endparam_v() << endl;
    cout << "Number of control points:  u= " << spline_surf2->numCoefs_u()
	 << "   v= " <<  spline_surf2->numCoefs_v()<< endl;

    // Write surface 2 to file
    fout << "200 1 0  4 50 100 0 255" << std::endl;
    fout << *spline_surf2;    // write spline surface data.

    // Write extra input points to file
    fout << "400 1 0 4 0 255 0 255" << endl; // Header.
    fout << num_extrapoints << endl;
    for (int i = 0; i < num_extrapoints; ++i) {
	int ip = i*dim;
	Go::Point inp_point(points[ip], points[ip+1], points[ip+2]);
	fout << inp_point << endl;  // write input point coordinates.
    }
    fout.close();

    // cout << "Open the file 'approx_surface.g2' in 'goview' to look at the results"
    //      << endl;
}
