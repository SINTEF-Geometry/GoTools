//===========================================================================
//                                                                           
// File: test_cpu_speed.C                                                    
//                                                                           
// Created: Wed Oct 12 15:59:08 2005                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: test_cpu_speed.C,v 1.2 2006-03-03 15:50:53 jbt Exp $
//                                                                           
// Description: Application to test the speed of a cpu-based implementation
//              of bezier patch self-intersection scheme.
//              Part of the routine will be used to create input to the gpu
//              test (the creation of the transformation matrix).
//              Currently assuming cubic in both directions.
//                                                                           
//===========================================================================


#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/GeometryTools.h"
#include "sislP.h"
#include "GoTools/geometry/SplineDebugUtils.h"

#include <vector>
#include <fstream>
#include <iomanip>


using namespace Go;
using std::vector;
using std::max;
using std::ifstream;
using std::shared_ptr;


// Create transformation matrix from input bezier curve to spline
// space with the inserted knots.
// Assuming that isert_knots all lie inside end parameters, sorted.

void transformationMatrix(const SISLCurve* cv,
			  const vector<double>& insert_knots,
			  vector<double>& transf_mat);
// 			  vector<int>& first_ind, vector<int>& last_ind);


// Create the cv refined coefficients by multiplying current
// coefficients with the refinement matrix.

vector<double> refinedCoefficients(const SISLCurve* cv,
				   vector<double> transf_mat); //,
// 				   vector<int> first_ind,
// 				   vector<int> last_ind);


int main(int argc, char** argv)
{

    if (argc != 3) {
	std::cout << "Usage: bezier_patch (or spline_sf, "
		  << "if you insist) nmb_divisions" << std::endl;
	return -1;
    }

    ObjectHeader header;

    // Read the curve from file.
    ifstream filein(argv[1]);
    if (filein.bad()) {
	std::cerr << "File #1 error (no file or corrupt file specified)."
		  << std::endl;
	return -1;
    }
    int nmb_div(atoi(argv[2])); // The bezier patch is split into
			        // (2^nmb_div)*(2^nmb_div) patches.

    header.read(filein);
    shared_ptr<SplineSurface> sf(new SplineSurface());
    sf->read(filein);
    filein.close();

    shared_ptr<SplineSurface> total_normal_sf(sf->normalSurface());

#ifdef INTERSECTIONS_DEBUG
    objToFile(total_normal_sf.get(), "total_normal_sf.g2");
#endif //INTERSECTIONS_DEBUG

    // If input surface is not bezier we split.
    vector<SplineSurface> patches;
    splitSurfaceIntoPatches(*sf, patches);

    // For each Bezier segment we create the transformation matrix
    // (and do some other stuff).
    for (size_t ki = 0; ki < patches.size(); ++ki) {
	// Assuming we're operating on the unit square.
	double tmin = 0.0;
	double tmax = 1.0;
	if ((patches[ki].startparam_u() != tmin) ||
	    (patches[ki].startparam_v() != tmin) ||
	    (patches[ki].endparam_u() != tmax) ||
	    (patches[ki].endparam_v() != tmax)) {
	    patches[ki].setParameterDomain(tmin, tmax, tmin, tmax);
	}
	int dim = patches[ki].dimension();
	bool rat = patches[ki].rational();
	int rdim = dim + rat;
	if (patches[ki].order_u()
	    != patches[ki].order_v()) { // We'll focus on equal order,
					// so ...
	    int raise_u = max(0, patches[ki].order_v()
			      - patches[ki].order_u());
	    int raise_v = max(0, patches[ki].order_u()
			      - patches[ki].order_v());
	    Point test_pt1 = patches[ki].ParamSurface::point(0.6132, 0.12994);
	    patches[ki].raiseOrder(raise_u, raise_v);
	    Point test_pt2 = patches[ki].ParamSurface::point(0.6132, 0.12994);
	    double dist = test_pt1.dist(test_pt2);
	    if (dist > 0.0) {
		std::cout << "Distance between points pre and "
			  << "post degree raising: "
			  << dist << std::endl;
	    }
	}
// 	const int order = 4;
// 	ASSERT(ik1 == ik2 == order);
	shared_ptr<SplineSurface> normal_sf(patches[ki].normalSurface());

	// A normal on input surface should equal the corr point on
	// the normal surface:
	Point pt1;
	patches[ki].normal(pt1, 0.6132, 0.12994);
	Point pt2 = normal_sf->ParamSurface::point(0.6132, 0.12994);
	pt2.normalize();
	double dist = pt1.dist(pt2);
	if (dist > 0.0) {
	    std::cout << "Distance between normal on surface and "
		      << "point on normal-surface: "
		      << dist << std::endl;
	}

	int ik1 = normal_sf->order_u();
	int ik2 = normal_sf->order_v();
// 	int in1 = normal_sf->numCoefs_u();
	int in2 = normal_sf->numCoefs_v();
	
	bool normal_sf_rat = (normal_sf->rational());

	// Using sisl routines we transfer the surface to a
	// SISLSurface.
	shared_ptr<SISLSurf> sisl_sf(GoSurf2SISL(*normal_sf));

	// The sisl_sf represented as a cv parametrized in the
	// v-direction (with dimension dim*in1).
	// Evaluated in a parameter results in an iso-curve
	// parametrized in the u-direction.
	double* coefs = (normal_sf_rat) ? sisl_sf->rcoef : sisl_sf->ecoef;
	shared_ptr<SISLCurve> u_cvs(newCurve(in2, ik2, sisl_sf->et2, coefs,
					     1, sisl_sf->in1*rdim, 1));

	// We're assuming uniform knot distribution (could of course
	// alter this, for instance based on rough curvature
	// analysis).
	vector<double> insert_knots; // Inserting 4 knots for each
				     // split parameter.
	int nmb_segments = int(std::pow((double)2, (double)nmb_div));
	double tstep = (tmax - tmin)/nmb_segments;
	// Not inserting tmin & tmax as they already exist.
	int order = ik1;
	for (int kj = 1; kj < nmb_segments; ++kj) {
	    double tpar = tmin + kj*tstep;
	    insert_knots.insert(insert_knots.end(), order, tpar);
	}

	// We then create the transformation matrix.
	// Hmm, as we're possibly operating on the normal surface we
	// cat not restrict ourselves to order 4.
	vector<double> transf_mat; // The compressed transformation
				   // matrix.
// 	vector<int> first_ind, last_ind; // Index of first and last
// 					 // non-zero basis function.
	transformationMatrix(u_cvs.get(), insert_knots, transf_mat);

	// We next verify that the created transformation matrix does
	// indeed yield the wanted spline space.
	vector<double> refined_coefs
	    = refinedCoefficients(u_cvs.get(), transf_mat);

	// The coefficients should equal those created when inserting
	// the knots into the cv.
	SISLCurve* refined_cv = NULL;
	int kstat = 0;
	s1018(u_cvs.get(), &insert_knots[0], (int)insert_knots.size(),
	      &refined_cv, &kstat);
	if (kstat < 0) {
	    std::cout << "Failed refining curve, exiting!" << std::endl;
	    return -1;
	}

	// The number of coefficients should equal size of
	// refined_coefs.
	bool cv_rat = (u_cvs->ikind == 2 || u_cvs->ikind == 4);
	int ref_cv_rdim = refined_cv->idim + (cv_rat);
	if (int(refined_coefs.size()) != (refined_cv->in)*ref_cv_rdim) {
	    std::cout << "Size mismatch between refined coefficients and "
		      << "refined curve!" << std::endl;
	    if (refined_cv) freeCurve(refined_cv);
	    return -1;
	}

	// We then run through the refined_cv calculating the
	// greateset difference.
	double max_dist = -1.0;
	int max_ind = -1;
	coefs = (cv_rat) ? refined_cv->rcoef : refined_cv->ecoef;
	for (int kj = 0; kj < int(refined_coefs.size()); ++kj) {
	    double dist = fabs(refined_coefs[kj] - coefs[kj]);
	    if (dist > max_dist) {
		max_dist = dist;
		max_ind = (kj*rdim)/rdim;
	    }
	}
	std::cout << "Max dist (achieved in pt # "
		  << max_ind << "): " << max_dist << std::endl;

	// We then create the refined spline surface by inserting
	// knots and then extract the Bezier patches to be tested for
	// intersection with a suitable zero-box.
	shared_ptr<SplineSurface> refined_normal_sf(normal_sf->clone());
	refined_normal_sf->insertKnot_u(insert_knots);
	refined_normal_sf->insertKnot_v(insert_knots);

	vector<SplineSurface> refined_normal_patches;
	splitSurfaceIntoPatches(*refined_normal_sf, refined_normal_patches);
	
	BoundingBox zero_box(3);
	double geo_deg_tol = 1e-17;
	zero_box.setFromPoints(Point(-geo_deg_tol, -geo_deg_tol, -geo_deg_tol),
			       Point(geo_deg_tol, geo_deg_tol, geo_deg_tol));
	int nmb_overlaps = 0;
	for (size_t kj = 0; kj < refined_normal_patches.size(); ++kj) {
	    BoundingBox bd_box = refined_normal_patches[kj].boundingBox();
	    if (bd_box.overlaps(zero_box)) {
		++nmb_overlaps;
	    }
	}

	std::cout << "Number of overlapping bd_boxes: "
		  << nmb_overlaps << std::endl;
    }

    return 1;
}


// Routine based on s1018 ('Insert a given set of knots into the
// description of a B-spline curve).

void transformationMatrix(const SISLCurve* cv,
			  const vector<double>& insert_knots,
			  vector<double>& transf_mat)
// 			  vector<int>& first_ind, vector<int>& last_ind)
{

    int kpl, kfi, kla; // To position elements in trans.-matrix.
    int kstat = 0;

//     int kdim = cv->idim;
    int kk = cv->ik;
    int kn = cv->in;
    int ref_in = kn + (int)insert_knots.size(); // The number of vertices
					        // in refined space.
    double* st = cv->et;
    vector<double> new_knots(st, st + kk + kn);
    new_knots.insert(new_knots.begin() + kk,
		     insert_knots.begin(), insert_knots.end());
    transf_mat.resize(ref_in*kk, 0.0); // @@sbr Currently not using
				       // the compression.
    vector<double> sp(kk, 0.0);

    // Updating the coefficientvector to the new curve.
    int ki, kmy;
    for (ki = 0, kmy = 0; ki < ref_in; ki++) {
	// Here we compute a new line with line number ki of the knot
	// insertion matrix.
	
	while (st[kmy + 1] <= new_knots[ki])
	    kmy++;
	s1701 (ki, kmy, kk, kn, &kpl, &kfi, &kla, &new_knots[0], st,
	       &sp[0], &transf_mat[ki*kn], &kstat);
	if (kstat) {
	    throw;
	}
	
// 	for (int kj = kfi; kj < kla + 1; ++kj) {
// 	    transf_mat.push_back(salfa[ki*kk+kj+kpl]);
// 	}
// 	first_ind.push_back(kfi);
// 	last_ind.push_back(kla);

	// Currently not using the compression.
	// @@sbr Remove if compression is to be utilized!
	if (kpl != 0) {
	    for (int kj = kfi; kj < kla + 1; ++kj) {
		transf_mat[ki*kk+kj] = transf_mat[ki*kk+kj+kpl];
		transf_mat[ki*kk+kj+kpl] = 0.0;
	    }
	}

// 	// Compute the kdim vertices with the same "index".
// 	for (int kj = 0; kj < kdim; kj++, s1++) {
// 	    for (*s1 = 0, kj1 = kfi, kj2 = kfi + kpl;
// 		 kj1 <= kla; kj1++, kj2++) {
// 		ki2 = kj1 * kdim + kj;
// 		*s1 += salfa[kj2] * coef[ki2];
// 	    }
// 	}

    }

    return;
}


vector<double> refinedCoefficients(const SISLCurve* cv,
				   vector<double> transf_mat) //,
// 				   vector<int> first_ind, vector<int> last_ind)
{
    int dim = cv->idim;
    bool rat = (cv->ikind ==2 || cv->ikind == 4);
    int rdim = dim + rat;
//     int kk = cv->ik;
    int kn = cv->in;
    // @@sbr As the gpu is based around 4-tuples we'll possibly let
    // the matrices have the same size.
    int nmb_matrices =  (int)transf_mat.size()/kn; //first_ind.size();
    int new_kn = nmb_matrices; //kn + 1;
    vector<double> refined_coefs(rdim*new_kn, 0.0);
    double* coefs = (rat) ? cv->rcoef : cv->ecoef;
//     int mat_curr_id = 0;
    for (int ki = 0; ki < new_kn; ++ki) {
// 	for (int kj = first_ind[ki]; kj < last_ind[ki] + 1; ++kj) {
	for (int kj = 0; kj < kn; ++kj) {
	    for (int kk = 0; kk < rdim; ++kk) {
// 		refined_coefs[ki*rdim+kk]
// 		    += coefs[kj*rdim+kk]
// 		    * transf_mat[mat_curr_id+kj-first_ind[ki]];
		refined_coefs[ki*rdim+kk]
		    += coefs[kj*rdim+kk] * transf_mat[ki*kn+kj];
	    }
	}
// 	mat_curr_id += last_ind[ki] + 1 - first_ind[ki];
    }

    return refined_coefs;
}
