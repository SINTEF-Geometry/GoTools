//===========================================================================
//                                                                           
// File: coons_patch_gen.C                                                       
//                                                                           
// Description:
//
// This program demonstrates the use of some of the functions in namespace
// CoonsPatchGen.
// The functions can be used to create a Coons Patch or a Gordon Surface.
// The functions returns a SplineSurface pointer to the created surface.
//
// The first example creates a Coons patch defined by four connected boundary
// curves. They must all be of type 'SplineCurve', nonrational, and form a loop.
// The end point of one curve must coincide with the start point of the next
// curve.
// Output is a file in Go-format for plot of the curve points and the surface.
// The file name is "coons_patch_surface.g2".
//
// The second example creates a Gordon surface defined by a curve mesh.
// We demand a regular grid in the parameter domain, i.e. we require a u- (v-)
// curve to cross v- (u-) curves at same parameter value in the v- (u-) curve.
// Here we extract the mesh curves from the surface of the first example,
// makes a Gordon surface from the curves and writes the Gordon surface and
// curves to a file named 'gordon_surface.g2'
//
//===========================================================================

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/creators/ApproxCurve.h"
#include "GoTools/utils/Point.h"
#include <vector>
#include <fstream>
#include <set>
#include <iterator>

using std::cout;
using std::vector;
using std::set;
using std::ofstream;
using std::endl;
using std::ostream_iterator;
using namespace Go;

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
    const int dim = 3;         // Space dimension.
    const double tol = 1.e-4;  // Geometric tolerance. 
    cout << "\nRunning program " << argv[0] << endl;

    // Example 1. Create a Coons patch defined by four boundary curves.
    // Define the four boundary curves of the surface.
    // The endpoints of the curves must be connected so that they form a loop.
    std::vector<shared_ptr<ParamCurve> > surf_boundary_curves(4);
    double start_param = 0.0;
    double end_param   = 1.0;
    int num_crv_pts = 3;   // Three points per boundary curve.
    std::vector<double> crv_param(num_crv_pts);  // Parameter values
    fill(start_param, end_param, num_crv_pts, crv_param); 

    // Surface corner points. 
    const int num_bnd_crvs = 4;
    double bnd_corner_x[] = {1, 2, 2, 1};
    double bnd_corner_y[] = {3, 3, 4, 4};
    double bnd_corner_z[] = {5, 6, 5, 6};
    
    vector<double> crv_points(dim*num_crv_pts);  // Curve points
    set<Point> bnd_pnts;
    for (int i=0; i<num_bnd_crvs; ++i) {
	// Create points on a boundary curve
	vector<double> x(num_crv_pts), y(num_crv_pts), z(num_crv_pts);
	fill(bnd_corner_x[i], bnd_corner_x[(i+1)%4], num_crv_pts, x);  
	fill(bnd_corner_y[i], bnd_corner_y[(i+1)%4], num_crv_pts, y);  
	fill(bnd_corner_z[i], bnd_corner_z[(i+1)%4], num_crv_pts, z);
	z[1]+=0.1;
	for (int j=0, p=-1; j<num_crv_pts; ++j) {
	    crv_points[++p] = x[j];
	    crv_points[++p] = y[j];
	    crv_points[++p] = z[j];
	    Point go_point(x[j], y[j], z[j]);
	    bnd_pnts.insert(go_point);
	}

	// Create a boundary curve
	double maxdist, avdist;  // Maximum and average distance between the
	                         // generated curve and the data points.
	const int max_iter = 5;  // Maximum number of iterations to use.
	ApproxCurve approx_curve(crv_points, crv_param, dim, tol);
	surf_boundary_curves[i] =    // Get the surface boundary spline curve.  
	    approx_curve.getApproxCurve(maxdist, avdist, max_iter);    
    }

    // Create the boundary curve loop
    CurveLoop bnd_curve_loop(surf_boundary_curves, tol);

    // Create the surface
    SplineSurface* 
	coons_patch_surface(CoonsPatchGen::createCoonsPatch(bnd_curve_loop));
    cout << "\nCoons patch surface\n";
    cout << "Parameter u from " << coons_patch_surface->startparam_u() << " to "
	 << coons_patch_surface->endparam_u() << endl;
    cout << "Parameter v from " << coons_patch_surface->startparam_v() << " to "
	 << coons_patch_surface->endparam_v() << endl;
    cout << "Number of control points: u= " << coons_patch_surface->numCoefs_u()
	 << "   v= " <<  coons_patch_surface->numCoefs_v()<< endl;
    cout << "Bounding box = " << coons_patch_surface->boundingBox() << endl;
 
    // Write the surface to file
    ofstream fout("coons_patch_surface.g2");
    coons_patch_surface->writeStandardHeader(fout); // write header.
    fout << *coons_patch_surface;    // write coons_patch_surface data.

    // Write the input points to file
    // Class_PointCloud=400 MAJOR_VERSION=1 MINOR_VERSION=0 auxillary data=4
    // The four auxillary data values defines the colour (r g b alpha)
    fout << "400 1 0 4 255 0 0 255" << endl; // Header.
    fout << bnd_pnts.size() << endl;
    copy(bnd_pnts.begin(), bnd_pnts.end(), ostream_iterator<Point>(fout, "\n"));  
    fout.close();


    // Example 2. Create a Gordon surface defined by a curve mesh.
    // Extract mesh curves from the Coons patch surface by calling it's
    // 'constParamCurve' function.
    vector<shared_ptr<SplineCurve> > mesh_curves;
    int nmb_u_crvs = 4;
    int nmb_v_crvs = 3;
    vector<double> params(nmb_u_crvs + nmb_v_crvs);
    // Constant v parameter values for curves parametrized in the u-direction.
    params[0]=0.0; params[1]=0.333; params[2]=0.667; params[3]=1.0;
    // Constant u parameter values for curves parametrized in the v-direction.
    params[4]=0.0; params[5]=0.500; params[6]=1.0;

    // u-curves at fixed v-values
    for (int i=0; i<nmb_u_crvs; ++i) {
	mesh_curves.push_back(shared_ptr<SplineCurve>
			      (coons_patch_surface->
			       constParamCurve(params[i], true)));
    }
    // v-curves at fixed u-values 
    for (int i=nmb_u_crvs; i<nmb_u_crvs+nmb_v_crvs; ++i) {
	mesh_curves.push_back(shared_ptr<SplineCurve>
			      (coons_patch_surface->
			       constParamCurve(params[i], false)));
    }

    // Copy the mesh curves. They might be changed in createGordonSurface
    vector<shared_ptr<SplineCurve> > org_mesh_curves = mesh_curves;

    // Create Gordon surface
    SplineSurface* gordon_surface =
	CoonsPatchGen::createGordonSurface(mesh_curves, params, nmb_u_crvs,
					   true);
    cout << "\nGordon surface\n";
    cout << "Parameter u from " << gordon_surface->startparam_u() << " to "
	 << gordon_surface->endparam_u() << endl;
    cout << "Parameter v from " << gordon_surface->startparam_v() << " to "
	 << gordon_surface->endparam_v() << endl;
    cout << "Number of control points:  u= " << gordon_surface->numCoefs_u()
	 << "   v= " <<  gordon_surface->numCoefs_v()<< endl;
    cout << "Bounding box = " << gordon_surface->boundingBox() << endl;

    // Write Gordon surface and mesh curves to file.
    ofstream fout2("gordon_surface.g2");
    gordon_surface->writeStandardHeader(fout2); // write header.
    gordon_surface->write(fout2);  // write surface data.
     for (int i = 0; i < (int)org_mesh_curves.size(); ++i) // 
	{
	org_mesh_curves[i]->writeStandardHeader(fout2);
	org_mesh_curves[i]->write(fout2);
    }
    fout2.close();

    delete coons_patch_surface;
    delete gordon_surface;

    // cout << "\nOpen the files 'coons_patch_surface.g2' and 'gordon_surface.g2'"
    //      << " in 'goview' \nto look at the results\n"
    //      << endl;
}
