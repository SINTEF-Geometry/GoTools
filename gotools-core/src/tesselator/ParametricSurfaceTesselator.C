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

#include "GoTools/tesselator/ParametricSurfaceTesselator.h"
#include <fstream>
// #include <qmessagebox.h>
#include <memory>
#include "GoTools/creators/CurveCreators.h"
#include "GoTools/creators/CreatorsUtils.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/tesselator/spline2mesh.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/Plane.h"

//#define VIEWLIB_DEBUG

using std::vector;

namespace Go
{


//===========================================================================
ParametricSurfaceTesselator::~ParametricSurfaceTesselator()
//===========================================================================
{
}


//===========================================================================
void ParametricSurfaceTesselator::changeRes(int n, int m)
//===========================================================================
{
    if ((m != m_) || (n != n_)) {
	m_ = m;
	n_ = n;
	tesselate();
    }
}


//===========================================================================
void ParametricSurfaceTesselator::tesselate()
//===========================================================================
{
    vector<shared_ptr<ParamCurve> > par_cv;
    shared_ptr<SplineSurface> spline_sf;
    shared_ptr<BoundedSurface> bd_sf;
    int dim = surf_.dimension();

    // @@sbr201506 The tolerance should be given as input to make_trimmed_mesh.
    double tol2d = 1.0e-8;//12;//4; // Tolerance used to check if a surface
    // is trimmed along iso parametric curves

    int ki;
    double umin, umax, vmin, vmax;
    bool rectangular_domain = false;

    // Check if the domain is rectangular. In that case, a simpler
    // tesselation may be applied
    if (surf_.instanceType() == Class_BoundedSurface) {
        bd_sf = shared_ptr<BoundedSurface> (
            dynamic_cast<BoundedSurface*> (surf_.clone()));
        RectDomain domain = surf_.containingDomain(); // Equals real domain
        umin = domain.umin();
        umax = domain.umax();
        vmin = domain.vmin();
        vmax = domain.vmax();
    }
    else {
        // All other surfaces are fine, except for unbounded
        // ElementarySurfaces
        const ElementarySurface* elemsf
            = dynamic_cast<const ElementarySurface*>(&surf_);
        bool is_unbounded_elementary = (elemsf && !elemsf->isBounded());
        if (!is_unbounded_elementary) {
            RectDomain domain = surf_.containingDomain(); // Equals real domain
            umin = domain.umin();
            umax = domain.umax();
            vmin = domain.vmin();
            vmax = domain.vmax();
            rectangular_domain = true;

        }
    }

    if (bd_sf.get() && bd_sf->isIsoTrimmed(tol2d)) {
        // Get surrounding domain
        RectDomain domain = bd_sf->containingDomain();

        // Get smallest surrounding surface
        shared_ptr<ParamSurface> base_sf = bd_sf->underlyingSurface();
        while (base_sf->instanceType() == Class_BoundedSurface)
            base_sf = dynamic_pointer_cast<BoundedSurface, ParamSurface>(
                    base_sf)->underlyingSurface();
        RectDomain dom2 = base_sf->containingDomain(); // To avoid
                                                // problems due to numerics
        umin = std::max(domain.umin(), dom2.umin());
        umax = std::min(domain.umax(), dom2.umax());
        vmin = std::max(domain.vmin(), dom2.vmin());
        vmax = std::min(domain.vmax(), dom2.vmax());
        rectangular_domain = true;
    }

    if (rectangular_domain) {
        Point pt(dim);
        mesh_->resize(n_ * m_, 2 * (n_ - 1) * (m_ - 1));
        int iu, iv, idx;
        for (iu = 0; iu < n_; ++iu) {
            for (iv = 0; iv < m_; ++iv) {
                double ru = double(iu) / double(n_ - 1);
                double rv = double(iv) / double(m_ - 1);
                double u = umin * (1.0 - ru) + ru * umax;
                double v = vmin * (1.0 - rv) + rv * vmax;
                surf_.point(pt, u, v);
                //	    std::cout << pt << std::endl;
		int j;
		for (j=0; j<dim; ++j)
		  mesh_->vertexArray()[(iv*n_ + iu)*3+j] = pt[j];
		for (; j<3; ++j)
		  mesh_->vertexArray()[(iv*n_ + iu)*3 + j] = 0.0;
                mesh_->paramArray()[(iv * n_ + iu) * 2] = u;
                mesh_->paramArray()[(iv * n_ + iu) * 2 + 1] = v;
                mesh_->boundaryArray()[iv * n_ + iu] = (iv == 0 || iv == m_ - 1
                        || iu == 0 || iu == n_ - 1) ? 1 : 0;
                if (mesh_->useNormals()) {
                    surf_.normal(pt, u, v);
                    mesh_->normalArray()[(iv * n_ + iu) * 3] = pt[0];
                    mesh_->normalArray()[(iv * n_ + iu) * 3 + 1] = pt[1];
                    mesh_->normalArray()[(iv * n_ + iu) * 3 + 2] = pt[2];
                }
                if (mesh_->useTexCoords()) {
                    mesh_->texcoordArray()[(iv * n_ + iu) * 2] = ru;
                    mesh_->texcoordArray()[(iv * n_ + iu) * 2 + 1] = rv;
                }
            }
        }
        // This is really a rectangular mesh. It remains to
        // create the triangle indicies
        for (iv = 0, idx = 0; iv < m_ - 1; ++iv) {
            for (iu = 0; iu < n_ - 1; ++iu) {
                mesh_->triangleIndexArray()[idx++] = iv * n_ + iu;
                mesh_->triangleIndexArray()[idx++] = iv * n_ + iu + 1;
                mesh_->triangleIndexArray()[idx++] = (iv + 1) * n_ + iu + 1;

                mesh_->triangleIndexArray()[idx++] = iv * n_ + iu;
                mesh_->triangleIndexArray()[idx++] = (iv + 1) * n_ + iu + 1;
                mesh_->triangleIndexArray()[idx++] = (iv + 1) * n_ + iu;
            }
        }
    }
    else if (bd_sf.get()) {
        // We must first extract the boundary domain.
        shared_ptr<ParamSurface> under_sf = bd_sf->underlyingSurface();
        //shared_ptr<SplineSurface> spline_sf;
        if (under_sf->instanceType() >= Class_Plane
                && under_sf->instanceType() <= Class_Torus) {
            shared_ptr<ElementarySurface> elem 
                = dynamic_pointer_cast<ElementarySurface>(under_sf);
            //spline_sf = shared_ptr<SplineSurface> (elem->geometrySurface());
            RectDomain domain = surf_.containingDomain();
            double umin = domain.umin();
            double umax = domain.umax();
            double vmin = domain.vmin();
            double vmax = domain.vmax();
            //double udel = umax - umin;
            //double vdel = vmax - vmin;
            //RectDomain domain2 = spline_sf->containingDomain();
            //umin = std::max(umin - 0.1 * udel, domain2.umin());
            //umax = std::min(umax + 0.1 * udel, domain2.umax());
            //vmin = std::max(vmin - 0.1 * vdel, domain2.vmin());
            //vmax = std::min(vmax + 0.1 * vdel, domain2.vmax());
            under_sf = shared_ptr<ParamSurface> ((under_sf->subSurfaces(umin,
                    vmin, umax, vmax))[0]);
        }

        vector<CurveLoop> bd_loops = bd_sf->absolutelyAllBoundaryLoops();
        for (int crv = 0; crv < int(bd_loops.size()); crv++) {
            for (ki = 0; ki < bd_loops[crv].size(); ++ki) {
                shared_ptr<CurveOnSurface> cv_on_sf(dynamic_pointer_cast<
                        CurveOnSurface, ParamCurve> (bd_loops[crv][ki]));
                if (cv_on_sf.get() == 0) {
                    THROW("Missing curve on surface, needed for tesselation!");
                }
                double eps = bd_loops[0].getSpaceEpsilon();
                cv_on_sf->ensureParCrvExistence(eps);
                shared_ptr<ParamCurve> pcv = cv_on_sf->parameterCurve();
		if (pcv.get() == NULL) {
                    THROW("Missing parameter curve, needed for tesselation!");
                }
                shared_ptr<SplineCurve> spline_cv(pcv->geometryCurve());
                if (ki == 0) {
                    // We do not want to alter sf...
                    par_cv.push_back(spline_cv);
                }
                else {
                    double dummy_dist;
                    par_cv[crv]->appendCurve(spline_cv->clone(), 0, dummy_dist,
                            false);
                }
            }
        }

        // We then tesselate the object.
        vector<Vector3D> trimmed_vert; // 3D vertices.
        vector<Vector2D> trimmed_par; // Corresponding 2D vertices, includes the regular (m_+1)x(n_+1)-grid.
        vector<int > trimmed_bd; // 1 == at_boundary, 0 == !at_boundary.
        vector<Vector3D> trimmed_norm; // Corresponding normal.
//        vector<Vector3D> trimmed_col;
        vector<int> trimmed_mesh; // Index of triangles sent to OpenGL. Refers to trimmed_vert (and trimmed_par).
        vector<Vector3D> trim_curve; // Not used on the outside.
        vector<Vector3D> trim_curve_p; // Not used on the outside.
        //int n = (m_ + n_)/2;
        //vector< Vector3D > extra_v;
        double bd_res_ratio = 1.0;
        {
            make_trimmed_mesh(under_sf, par_cv, trimmed_vert, trimmed_par,
                    trimmed_bd, trimmed_norm, trimmed_mesh, trim_curve,
                    trim_curve_p, n_, m_, bd_res_ratio);
        }

        // Finally we must update values in mesh_.
        int nmb_vert = (int)trimmed_vert.size();
        int nmb_triangles = (int)trimmed_mesh.size() / 3;
        mesh_->resize(nmb_vert, nmb_triangles);
        for (ki = 0; ki < int(trimmed_vert.size()); ++ki) {
            mesh_->vertexArray()[ki * 3] = trimmed_vert[ki][0];
            mesh_->vertexArray()[ki * 3 + 1] = trimmed_vert[ki][1];
            mesh_->vertexArray()[ki * 3 + 2] = trimmed_vert[ki][2];
            mesh_->paramArray()[ki * 2] = trimmed_par[ki][0];
            mesh_->paramArray()[ki * 2 + 1] = trimmed_par[ki][1];
            mesh_->boundaryArray()[ki] = trimmed_bd[ki];
            if (mesh_->useNormals()) {
                mesh_->normalArray()[ki * 3] = trimmed_norm[ki][0];
                mesh_->normalArray()[ki * 3 + 1] = trimmed_norm[ki][1];
                mesh_->normalArray()[ki * 3 + 2] = trimmed_norm[ki][2];
            }
            if (mesh_->useTexCoords()) {
                //mesh_->texcoordArray()[ki * 2] = trimmed_par[ki][0];
                //mesh_->texcoordArray()[ki * 2 + 1] = trimmed_par[ki][1];
                double s = (trimmed_par[ki][0] - umin) / (umax - umin);
                double t = (trimmed_par[ki][1] - vmin) / (vmax - vmin);
                mesh_->texcoordArray()[ki * 2] = s;
                mesh_->texcoordArray()[ki * 2 + 1] = t;
            }
        }

#ifdef VIEWLIB_DEBUG
	{
	    std::ofstream debug("tmp/debug.g2");
	    vector<double> pts;
	    vector<double> par_z;
	    for (ki = 0; ki < int(trim_curve.size()); ++ki)
	    {
		Point sf_pt = surf_.point(trim_curve_p[ki][0], trim_curve_p[ki][1]);
		pts.insert(pts.end(), sf_pt.begin(), sf_pt.end());
		par_z.push_back(trim_curve_p[ki][0]);
		par_z.push_back(trim_curve_p[ki][1]);
		par_z.push_back(0.0);
	    }

	    vector<double> vert_par;
	    for (ki = 0; ki < int(trimmed_par.size()); ++ki)
	    {
		vert_par.push_back(trimmed_par[ki][0]);
		vert_par.push_back(trimmed_par[ki][1]);
		vert_par.push_back(0.0);
	    }

	    vector<double> triang_nodes_par;
	    for (ki = 0; ki < int(trimmed_mesh.size()); ++ki)
	    {
		int ind = trimmed_mesh[ki];
		triang_nodes_par.push_back(trimmed_par[ind][0]);
		triang_nodes_par.push_back(trimmed_par[ind][1]);
		triang_nodes_par.push_back(0.0);
	    }

	    int nmb_points = pts.size()/3;
	    PointCloud<3> pt_cloud(pts.begin(), nmb_points);
	    pt_cloud.writeStandardHeader(debug);
	    pt_cloud.write(debug);
	    PointCloud<3> pt_cloud2(par_z.begin(), nmb_points);
	    pt_cloud2.writeStandardHeader(debug);
	    pt_cloud2.write(debug);
	    PointCloud<3> pt_cloud3(vert_par.begin(), vert_par.size()/3);
	    pt_cloud3.writeStandardHeader(debug);
	    pt_cloud3.write(debug);
	    PointCloud<3> pt_cloud4(triang_nodes_par.begin(), triang_nodes_par.size()/3);
	    pt_cloud4.writeStandardHeader(debug);
	    pt_cloud4.write(debug);
	    double debug_val = 0.0;
	}
#endif // VIEWLIB_DEBUG

	if (trimmed_mesh.size() > 0)
	{
	    copy(trimmed_mesh.begin(), trimmed_mesh.end(),
		 mesh_->triangleIndexArray());
	}
	else
	{
	    ;//MESSAGE("No trimmed mesh in output.");
	}
    }
    else
    {
        MESSAGE("Unexpected surface type, returning.");
        // 	    QMessageBox::warning( this, "Tesselating surface:
        // 				  ", "Unexpected surface type,
        // 				  returning.",
        // 				  QMessageBox::Ok,
        // 				  QMessageBox::NoButton);
        return;
    }

}


} // namespace Go

