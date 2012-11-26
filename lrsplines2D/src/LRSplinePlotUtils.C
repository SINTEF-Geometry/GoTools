//===========================================================================
//                                                                           
// File: LRSplinePlotUtils.C                                                 
//                                                                           
// Created: Fri Nov 16 14:23:20 2012                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/lrsplines2D/LRSplinePlotUtils.h"
#include "GoTools/lrsplines2D/Mesh2DIterator.h"

namespace Go
{

    void writePostscriptMesh(Go::LRSplineSurface& lr_spline_sf, std::ostream &out)
    {
	const bool close = true;
	const bool colorDiag = false;
	const Mesh2D& mesh = lr_spline_sf.mesh();
	const int num_diff_knots_u = mesh.numDistinctKnots(XFIXED);
	const int num_diff_knots_v = mesh.numDistinctKnots(YFIXED);
	// std::vector<double> knot_u(num_diff_knots_u), knot_v(num_diff_knots_v);
	const double* const knot_u = mesh.knotsBegin(XFIXED);
	const double* const knot_v = mesh.knotsBegin(YFIXED);
//    lr_spline_sf.getGlobalUniqueKnotVector(knot_u, knot_v);
	double min_span_u = knot_u[1] - knot_u[0];
	double min_span_v = knot_v[1] - knot_v[0];
	for(size_t i=1; i<num_diff_knots_u-1; i++)
	    min_span_u = (min_span_u<knot_u[i+1]-knot_u[i]) ? min_span_u : knot_u[i+1]-knot_u[i];
	for(size_t i=1; i<num_diff_knots_v-1; i++)
	    min_span_v = (min_span_v<knot_v[i+1]-knot_v[i]) ? min_span_v : knot_v[i+1]-knot_v[i];

	// get date
	time_t t = time(0);
	tm* lt = localtime(&t);
	char date[11];
	sprintf(date, "%02d/%02d/%04d", lt->tm_mday, lt->tm_mon + 1, lt->tm_year+1900);

	// get bounding box
	double umin = lr_spline_sf.paramMin(XFIXED);
	double umax = lr_spline_sf.paramMax(XFIXED);
	double vmin = lr_spline_sf.paramMin(YFIXED);
	double vmax = lr_spline_sf.paramMax(YFIXED);
	double dx = umax - umin;
	double dy = vmax - vmin;
	double scale = 1.0;//(dx>dy) ? 1000.0/dx : 1000.0/dy;
	// set the duplicate-knot-line (dkl) display width
	double dkl_range = (min_span_u>min_span_v) ? min_span_v*scale/6.0 : min_span_u*scale/6.0; 
	double xmin = (umin - dx/100.0)*scale;
	double ymin = (vmin - dy/100.0)*scale;
	double xmax = (umax   + dx/100.0)*scale + dkl_range;
	double ymax = (vmax   + dy/100.0)*scale + dkl_range;

	// print eps header
	out << "%!PS-Adobe-3.0 EPSF-3.0\n";
	out << "%%Creator: LRSplinePlotUtils.C object\n";
	out << "%%Title: LRSplineSurface parameter domain\n";
	out << "%%CreationDate: " << date << std::endl;
	out << "%%Origin: 0 0\n";
	out << "%%BoundingBox: " << xmin << " " << ymin << " " << xmax << " " << ymax << std::endl;

	out << "0 setgray\n";
	double linewidth = std::min(dx/100.0, dy/100.0);
	out << linewidth << " setlinewidth\n";
	Mesh2DIterator mesh_beg = mesh.begin();
	Mesh2DIterator mesh_end = mesh.end();
	Mesh2DIterator mesh_it = mesh_beg;
	while (mesh_it != mesh_end)
	{
	    const std::array<int, 4> elem_corners = *mesh_it; // ll_u, ll_v, ur_u, ur_v.
	    // Currently we do not care about multiplicity.
	    double dm = 0.0;// (mesh_iter[i]->multiplicity_==1) ? 0 : dkl_range/(mesh_iter[i]->multiplicity_-1);
	    double m = 1.0; // The multiplicity,
//	    int mult = mesh_it->multiplicity_;
	    // First we create the lines in the u-dir.
	    // We also do not care about double lines (for neighbour elements) ...
	    // out << mesh_iter[i]->start_*scale << " " << mesh_iter[i]->const_par_*scale + dm*m << " moveto\n";
	    // if(mesh_iter[i]->stop_ == umax)
	    //     out << mesh_iter[i]->stop_*scale+dkl_range << " " << mesh_iter[i]->const_par_*scale + dm*m << " lineto\n";
	    // else
	    //     out << mesh_iter[i]->stop_*scale << " " << mesh_iter[i]->const_par_*scale + dm*m << " lineto\n";
	    // umin
	    out << "newpath\n";
	    out << knot_u[elem_corners[0]]*scale << " " << knot_v[elem_corners[1]]*scale + dm*m << " moveto\n";
	    out << knot_u[elem_corners[2]]*scale << " " << knot_v[elem_corners[1]]*scale + dm*m << " lineto\n";
	    out << "stroke\n";

	    // umax
	    out << "newpath\n";
	    out << knot_u[elem_corners[0]]*scale << " " << knot_v[elem_corners[3]]*scale + dm*m << " moveto\n";
	    out << knot_u[elem_corners[2]]*scale << " " << knot_v[elem_corners[3]]*scale + dm*m << " lineto\n";
	    out << "stroke\n";

	    // vmin
	    out << "newpath\n";
	    out << knot_u[elem_corners[0]]*scale << " " << knot_v[elem_corners[1]]*scale + dm*m << " moveto\n";
	    out << knot_u[elem_corners[0]]*scale << " " << knot_v[elem_corners[3]]*scale + dm*m << " lineto\n";
	    out << "stroke\n";

	    // vmax
	    out << "newpath\n";
	    out << knot_u[elem_corners[2]]*scale << " " << knot_v[elem_corners[1]]*scale + dm*m << " moveto\n";
	    out << knot_u[elem_corners[2]]*scale << " " << knot_v[elem_corners[3]]*scale + dm*m << " lineto\n";
	    out << "stroke\n";

	    ++mesh_it;
	}

	if(close)
	    out << "%%EOF\n";
    }

} // end of namespace Go.
