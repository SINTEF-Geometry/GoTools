//===========================================================================
//                                                                           
// File: ParametricSurfaceTesselator.C                                          
//                                                                           
// Created:
//                                                                           
// Author:
//                                                                           
// Revision:
//                                                                           
// Description:
//                                                                           
//===========================================================================


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


using std::vector;
using std::shared_ptr;
using std::dynamic_pointer_cast;
using std::dynamic_pointer_cast;


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
    m_ = m;
    n_ = n;
    tesselate();
}


//===========================================================================
void ParametricSurfaceTesselator::tesselate()
//===========================================================================
{
    vector<shared_ptr<SplineCurve> > par_cv;
    shared_ptr<SplineSurface> spline_sf;
    shared_ptr<BoundedSurface> bd_sf;

    double tol2d = 1.0e-4; // Tolerance used to check if a surface
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
            base_sf = std::dynamic_pointer_cast<BoundedSurface, ParamSurface>(
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
        Point pt(3);
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
                mesh_->vertexArray()[(iv * n_ + iu) * 3] = pt[0];
                mesh_->vertexArray()[(iv * n_ + iu) * 3 + 1] = pt[1];
                mesh_->vertexArray()[(iv * n_ + iu) * 3 + 2] = pt[2];
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
        shared_ptr<SplineSurface> spline_sf;
        if (under_sf->instanceType() == Class_SplineSurface)
            spline_sf = dynamic_pointer_cast<SplineSurface> (under_sf);
        else if (under_sf->instanceType() >= Class_Plane
                && under_sf->instanceType() <= Class_Torus) {
            shared_ptr<ElementarySurface> elem = dynamic_pointer_cast<
                    ElementarySurface> (under_sf);
            spline_sf = shared_ptr<SplineSurface> (elem->geometrySurface());
            RectDomain domain = surf_.containingDomain();
            double umin = domain.umin();
            double umax = domain.umax();
            double vmin = domain.vmin();
            double vmax = domain.vmax();
            double udel = umax - umin;
            double vdel = vmax - vmin;
            RectDomain domain2 = spline_sf->containingDomain();
            umin = std::max(umin - 0.1 * udel, domain2.umin());
            umax = std::min(umax + 0.1 * udel, domain2.umax());
            vmin = std::max(vmin - 0.1 * vdel, domain2.vmin());
            vmax = std::min(vmax + 0.1 * vdel, domain2.vmax());
            spline_sf = shared_ptr<SplineSurface> (spline_sf->subSurface(umin,
                    vmin, umax, vmax));
        }
        ASSERT(spline_sf.get() != 0);

        vector<CurveLoop> bd_loops = bd_sf->absolutelyAllBoundaryLoops();
        for (int crv = 0; crv < int(bd_loops.size()); crv++) {
            for (ki = 0; ki < bd_loops[crv].size(); ++ki) {
                shared_ptr<CurveOnSurface> cv_on_sf(dynamic_pointer_cast<
                        CurveOnSurface, ParamCurve> (bd_loops[crv][ki]));
                ASSERT(cv_on_sf.get() != 0);
                // 	      shared_ptr<SplineCurve> spline_cv =
                // 		dynamic_pointer_cast<SplineCurve, ParamCurve>
                // 		(cv_on_sf->parameterCurve());
                shared_ptr<ParamCurve> pcv = cv_on_sf->parameterCurve();
                // 	    ASSERT(spline_cv.get() != 0);
                shared_ptr<SplineCurve> spline_cv =
                        (pcv.get() != 0) ? shared_ptr<SplineCurve> (
                                pcv->geometryCurve())
                                : shared_ptr<SplineCurve> ();
                if (spline_cv.get() == 0) {
                    // The parameter curve is not given - we must project
                    shared_ptr<ParamCurve> sc = cv_on_sf->spaceCurve();
                    // 			=
                    // 			dynamic_pointer_cast<SplineCurve,
                    // 			ParamCurve>
                    if (sc.get() == NULL)
                        THROW("Missing data needed for tesselating surface.");
                    shared_ptr<ParamSurface> sf = cv_on_sf->underlyingSurface();
                    shared_ptr<Point> pt;
                    double eps = bd_loops[0].getSpaceEpsilon();
                    spline_cv.reset(CurveCreators::projectSpaceCurve(sc, sf,
                            pt, pt, eps));
                    if (spline_cv.get() == 0) {
                        THROW("Error: Failed to project space curve.");
                    }
                }
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
        vector<Vector3D> trimmed_vert;
        vector<Vector2D> trimmed_par;
        vector<int > trimmed_bd;
        vector<Vector3D> trimmed_norm;
        vector<Vector3D> trimmed_col;
        vector<int> trimmed_mesh;
        vector<Vector3D> trim_curve;
        vector<Vector3D> trim_curve_p;
        //int n = (m_ + n_)/2;
        //vector< Vector3D > extra_v;
        double bd_res_ratio = 1.0;
        {
            make_trimmed_mesh(spline_sf, par_cv, trimmed_vert, trimmed_par,
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

        copy(trimmed_mesh.begin(), trimmed_mesh.end(),
                mesh_->triangleIndexArray());
    }
    else {
        MESSAGE("Unexpected surface type, returning.");
        // 	    QMessageBox::warning( this, "Tesselating surface:
        // 				  ", "Unexpected surface type,
        // 				  returning.",
        // 				  QMessageBox::Ok,
        // 				  QMessageBox::NoButton);
        return;
    }

#ifdef VIEWLIB_DEBUG
    std::ofstream debug("data/debug.g2");
    vector<double> pts;
    for (ki = 0; ki < int(trim_curve.size()); ++ki)
    {
        Point sf_pt = surf_.point(trim_curve_p[ki][0], trim_curve_p[ki][1]);
        pts.insert(pts.end(), sf_pt.begin(), sf_pt.end());
    }

    int nmb_points = pts.size()/3;
    PointCloud<3> pt_cloud(pts.begin(), nmb_points);
    pt_cloud.writeStandardHeader(debug);
    pt_cloud.write(debug);
#endif // VIEWLIB_DEBUG

//     int m = mesh_->numStrips() + 1;
//     int n = mesh_->numVertices()/m;
//     /// @@@ We can only tesselate properly rectangular-domain surfaces.
//     Go::RectDomain dom = surf_.containingDomain();
//     Go::Point pt(3);
//     for (int iu = 0; iu < n; ++iu) {
// 	for (int iv = 0; iv < m; ++iv) {
// 	    double ru = double(iu)/double(n-1);
// 	    double rv = double(iv)/double(m-1);
// 	    double u = dom.umin()*(1.0-ru) + ru*dom.umax();
// 	    double v = dom.vmin()*(1.0-rv) + rv*dom.vmax();
// 	    surf_.point(pt, u, v);
// 	    //	    std::cout << pt << std::endl;
// 	    mesh_->vertexArray()[(iv*n + iu)*3] = pt[0];
// 	    mesh_->vertexArray()[(iv*n + iu)*3 + 1] = pt[1];
// 	    mesh_->vertexArray()[(iv*n + iu)*3 + 2] = pt[2];
// 	    if (mesh_->useNormals()) {
// 		surf_.normal(pt, u, v);
// 		mesh_->normalArray()[(iv*n + iu)*3] = pt[0];
// 		mesh_->normalArray()[(iv*n + iu)*3 + 1] = pt[1];
// 		mesh_->normalArray()[(iv*n + iu)*3 + 2] = pt[2];
// 	    }
// 	    if (mesh_->useTexCoords()) {
// 		mesh_->texcoordArray()[(iv*n + iu)*2] = ru;
// 		mesh_->texcoordArray()[(iv*n + iu)*2+1] = rv;
// 	    }
// 	}
//     }
}


} // namespace Go

