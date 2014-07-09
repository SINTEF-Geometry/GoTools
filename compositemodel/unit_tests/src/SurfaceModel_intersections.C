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
//#define DEBUG

#include "GoTools/compositemodel/CellDivision.h"
#include "GoTools/utils/Point.h"
#include "GoTools/utils/Array.h"
#include "GoTools/compositemodel/ftEdgeBase.h"
#include "GoTools/compositemodel/ftCurve.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeCurve.h"
#include "GoTools/compositemodel/Body.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "sislP.h"
#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/creators/CurveCreators.h"
#include "GoTools/compositemodel/IntResultsSfModel.h"
#include "GoTools/topology/FaceAdjacency.h"
#include "GoTools/topology/FaceConnectivityUtils.h"
#include <fstream>


using std::vector;
using std::make_pair;

namespace Go
{

// ----------------- Helper functions and classes ---------------------

namespace // anon namespace
{



class ftEdgePointInfo
{
public:
    Point spacept_;
    Point parampt_;
    int face_id_;
    ftEdgeBase* edge_;
    bool closed_;
    ftEdgePointInfo(Point spacept, Point parampt, 
		    int face_id, ftEdgeBase* edge, bool cl = false)
	: spacept_(spacept), parampt_(parampt), face_id_(face_id), edge_(edge),
	closed_(cl)
	{}
};



//===========================================================================
vector<ftEdgePointInfo>
MakeEdgePointInfo(vector<ftCurveSegment>::iterator beg,
		  vector<ftCurveSegment>::iterator end,
		  double tolerance)
//===========================================================================
{
    int n = (int)(end - beg);
    vector<ftEdgePointInfo> epinfo;
    for (int i = 0; i < n; ++i) {
	const ftCurveSegment& s = *(beg+i);
	for (int j = 0; j < 2; ++j) {
	    // We assume that the first face and parameter curve is the correct
	    // one (i.e. that we have an intersection curve and not a
	    // kink curve for instance)
	    const int side = 0;
	    Point pt = (j == 0) ? s.startPoint() : s.endPoint();
	    Point par;
	    s.paramcurvePoint(side, j == 0 ?
			      s.startOfSegment()
			      : s.endOfSegment(),
			      par);
	    ftFaceBase* face = s.face(side);
	    // Now we have the parameters and the face

	    FaceConnectivityUtils<ftEdgeBase,ftFaceBase> connectivity;
	    vector<ftEdgeBase*> bedges = connectivity.edgesBoundingFace(face);
	    int bes = (int)bedges.size();
	    ALWAYS_ERROR_IF(bes == 0,
			    "Face had no edges in topology table.");

	    // Do a closest point on those edges
	    // @@ We do not handle the corner case very well
	    ftEdgeBase* edge = 0;
	    for (int k = 0; k < bes; ++k) {
		double clo_t, clo_dist;
		Point clo_pt;
		bedges[k]->closestPoint(pt, clo_t, clo_pt, clo_dist);
		if (clo_dist < tolerance) {
		    edge = bedges[k];
		}
	    }

	    // Now, if we're not on an edge, seal it off by setting the
	    // neighbour pointer to 0.
	    bool cl;
	    if (edge == 0) {
		cl = true;
	    } else {
		if (edge->twin()) {
		    cl = false;
		} else {
		    cl = true;
		}
	    }
	    epinfo.push_back(ftEdgePointInfo(pt, par, face->getId(),
					     edge, cl));
	}
    }
    return epinfo;
}

//===========================================================================
inline bool IsEven(int i)
//===========================================================================
{
    return (i % 2) == 0;
}


} // anon namespace



//===========================================================================
shared_ptr<IntResultsModel> SurfaceModel::intersect_plane(const ftPlane& plane)
//===========================================================================
{
  shared_ptr<IntResultsSfModel> intersections = 
    shared_ptr<IntResultsSfModel>(new IntResultsSfModel(this,
							plane)); // Empty storage for output
  ftCurve int_curves = intersect(plane);

  intersections->addIntCvs(int_curves);
  
  return intersections;
}


//===========================================================================
ftCurve SurfaceModel::intersect(const ftPlane& plane)
//===========================================================================
{
    // First, we make a list of cells that overlap the plane
    // Then, run intersection on each of the surfaces touching those cells,
    // if the bounding boxes overlap.

    ftCurve intcurve(CURVE_INTERSECTION);

    if (!celldiv_.get())
	initializeCelldiv();
    if (!plane.intersectsBox(celldiv_ -> big_box())) return intcurve;

    highest_face_checked_ = 0;
    for (int i = 0; i < celldiv_ -> numCells(); ++i) {
      ftCell aCell = celldiv_ -> getCell(i);
	if (plane.intersectsBox(aCell.box())) {
	    for (int j = 0; j < aCell.num_faces(); ++j) {
		ftSurface* face = aCell.face(j) -> asFtSurface();
		int id = face->getId();
		if (!face_checked_[id]) {
		    face_checked_[id] = true;
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
		    highest_face_checked_ = max(highest_face_checked_, id);
#else
		    highest_face_checked_ = 
		      std::max(highest_face_checked_, id);
#endif
		    if (plane.intersectsBox(face->boundingBox()))
			intcurve += localIntersect(plane, face);
		}
	    }
	}
    }
    std::fill(face_checked_.begin(),
	      face_checked_.begin() + highest_face_checked_ + 1, false);
    if (limit_box_.valid())
	intcurve.chopOff(limit_box_);
    intcurve.orientSegments(toptol_.neighbour);
    intcurve.joinSegments(toptol_.gap, toptol_.neighbour, toptol_.kink, toptol_.bend);
    return intcurve;
}







//===========================================================================
ftCurve SurfaceModel::localIntersect(const ftPlane& plane,
				     ftSurface* sf)
//===========================================================================
{
    int i, j, k1, k2;

    // Get all the intersection segments from this surface
    vector<ftCurveSegment> segments = intersect(plane, sf);
    int num_curves = (int)segments.size();

    ftCurve intcurve(CURVE_INTERSECTION);
    if (num_curves == 0) return intcurve;

    // Let each of those separate segments be the start of a new curve,
    // compute initial endpoints for those curves
    vector<ftCurve> curves;
    vector<ftEdgePointInfo> epinfo = MakeEdgePointInfo(segments.begin(),
						       segments.end(),
						       toptol_.neighbour);
    for (i = 0; i < num_curves; ++i) {
	curves.push_back(ftCurve(CURVE_INTERSECTION));
	curves[i].appendSegment(segments[i]);
    }

    // Scan through every edgepointinfo, and replace it by its connection
    // in the next face (connecting the curves), if any.

    vector<ftCurveSegment> new_segments;
    vector<ftEdgePointInfo> new_info;
    for (i = 0; i < 2*num_curves; ++i) {
	while (!epinfo[i].closed_) {
	    // We will trace this curve until it either runs into
	    // a face we've already treated or until it ends.
	    // Afterwards, we may have to join some curves.
	    ftEdgeBase* adjacent_edge = epinfo[i].edge_->twin();
	    ftSurface* adjacent_face = adjacent_edge -> face() -> asFtSurface();
	    int adjacent_face_id = adjacent_face->getId();
	    if (face_checked_[adjacent_face_id])
		break; // Out of the while-loop
	    // We're tracing the curve into adjacent_face.
	    // First we mark it:
	    face_checked_[adjacent_face_id] = true;
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
	    highest_face_checked_ = max(highest_face_checked_,
					adjacent_face_id);
#else
	    highest_face_checked_ = std::max(highest_face_checked_,
					adjacent_face_id);
#endif
	    // Then we get segments from it and analyze their endpoints
	    new_segments = intersect(plane, adjacent_face);
	    new_info = MakeEdgePointInfo(new_segments.begin(),
					 new_segments.end(),
					 toptol_.neighbour);
	    // Scan the new_info vector to identify the
	    // segment matching our previous curve
	    int nnseg = (int)new_segments.size();
	    int matching_newinfo = -1;
	    int j;
	    for (j = 0; j < 2*nnseg; ++j) {
		double d = epinfo[i].spacept_.dist(new_info[j].spacept_);
		if (d < toptol_.neighbour) {
		    matching_newinfo = j;
		    break; // Out of this local for-loop
		}
	    }
	    // Now matching_newinfo is either -1 (oops) or the index of a
	    // new info matching epinfo[i]
	    if (matching_newinfo == -1) {
		MESSAGE("Could not find a matching intcurve.");
		break; // Out of while-loop
	    }
	    // We have a match.
	    int matching_segment = matching_newinfo / 2;
	    if (IsEven(i)) {
		// We have to insert the new segment in front of the
		// matching curve. Reverse segment if necessary
		if(IsEven(matching_newinfo)) {
		    new_segments[matching_segment].reverse();
		}
		curves[i/2].prependSegment(new_segments[matching_segment]);
	    } else {
		// We append the new segment. Reverse if necessary
		if (!IsEven(matching_newinfo)) {
		    new_segments[matching_segment].reverse();
		}
		curves[i/2].appendSegment(new_segments[matching_segment]);
	    }
	    // Update the edgepoint info of the curve we have lengthened
	    int idx = matching_newinfo;
	    if (IsEven(matching_newinfo))
		++idx;
	    else
		--idx;
	    epinfo[i] = new_info[idx];
	    // Add any extra segments as their own curves
	    num_curves += nnseg - 1;
	    for (j = 0; j < nnseg; ++j) {
		if (j != matching_segment) {
		    curves.push_back(ftCurve(CURVE_INTERSECTION));
		    curves.back().appendSegment(new_segments[j]);
		    epinfo.push_back(new_info[2*j]);
		    epinfo.push_back(new_info[2*j + 1]);
		}
	    }
	    // We are done with this face
	} // End of: while(!epinfo[i].closed_)
    } // End of: for (int i = 0; i < 2*num_curves; ++i)

    // We now have all the intersection curves connected to this
    // local face, we have checked all the faces involved and
    // done all the curves connecting to them and so on.
    // We may have some cases where the curves we have should be
    // connected.

    // We have to connect any curves that should connect
    num_curves = (int)curves.size();
    Point pts[2][2];
    vector<bool> joined(num_curves, false);
    for (i = 0; i < num_curves - 1; ++i) {
	for (j = i + 1; j < num_curves; ++j) {
	    // Let's compare the endpoints of curves[i] and curves[j]
	    int seg = 0;
	    curves[i].point(curves[i].startOfSegment(seg), seg, pts[0][0]);
	    seg = curves[i].numSegments() - 1;
	    curves[i].point(curves[i].endOfSegment(seg), seg, pts[0][1]);
	    seg = 0;
	    curves[j].point(curves[j].startOfSegment(seg), seg, pts[1][0]);
	    seg = curves[j].numSegments() - 1;
	    curves[j].point(curves[j].endOfSegment(seg), seg, pts[1][1]);
	    for (k1 = 0; k1 < 2; ++k1) {
		for (k2 = 0; k2 < 2; ++k2) {
		    if (pts[0][k1].dist(pts[1][k2]) < toptol_.neighbour) {
			joined[i] = true;
			if (k2 == 0) {
			    // Prepend
			    if (k1 == 0)
				curves[i].reverse();
			    curves[j].prependCurve(curves[i]);
			} else {
			    // Append
			    if (k1 == 1)
				curves[i].reverse();
			    curves[j].appendCurve(curves[i]);
			}
		    }
		}
	    }
	    if (joined[i])
		break; // Out of the j-loop to start on the next i-value
	}
    }

    // We should now build intcurve from curves[]
    for (i = 0; i < (int)curves.size(); ++i)
	if (!joined[i])
	    intcurve += curves[i];
    return intcurve;
}




//===========================================================================
vector<ftCurveSegment> SurfaceModel::intersect(const ftPlane& plane,
					       ftSurface* sf)
//===========================================================================
{
    // Convert the surface to a SISLSurf in order to use SISL functions
    // on it. The "false" argument dictates that the SISLSurf will only
    // copy pointers to arrays, not the arrays themselves.
    SplineSurface* splinesf
	= dynamic_cast<SplineSurface*>(sf->surface().get());
    const CurveBoundedDomain* bdomain = 0;
    if (splinesf == 0) {
	BoundedSurface* bs
	    = dynamic_cast<BoundedSurface*>(sf->surface().get());
	if (bs != 0) {
	    splinesf
		= dynamic_cast<SplineSurface*>(bs->underlyingSurface().get());
	    bdomain = &(bs->parameterDomain());
	} else {
	    THROW("You cannot intersect this class of surface: "
		  << sf->surface()->instanceType());
	}
    }
    ASSERT(splinesf != 0);

    SISLSurf* sislsf = GoSurf2SISL(*splinesf, false);
    int dim = 3;
    double epsco = 1e-15; // Not used
    double epsge = 1e-6;
    int numintpt;  // number of single intersection points
    double* pointpar = 0; // array containing the parameter values of single intersect. pt.
    int numintcr; // number of intersection curves
    SISLIntcurve** intcurves = 0;
    int stat;
    Point pnt = plane.point();
    Point nrm = plane.normal();
    // Find the topology of the intersection
    s1851(sislsf, pnt.begin(), nrm.begin(), dim, epsco, epsge,
	  &numintpt, &pointpar, &numintcr, &intcurves, &stat);
    MESSAGE_IF(stat!=0, "s1851 returned code: " << stat);
    // pointpar is not used any further
    free(pointpar);
    double maxstep = 0.0;
    int makecurv = 2;     // Make both geometric and parametric curves
    int graphic = 0;      // Do not draw the curve
    epsge = toptol_.gap;
    vector<ftCurveSegment> segs;
    //    std::ofstream debug("data/debug.g2");
    for (int i = 0; i < numintcr; ++i) {
	// March out the intersection curves
	intcurves[i]->ipoint = 1;
	s1314(sislsf, pnt.begin(), nrm.begin(), dim, epsco, epsge,
	      maxstep, intcurves[i], makecurv, graphic, &stat);
	SISLCurve* sc = intcurves[i]->pgeom;
	if (sc == 0) {
	    MESSAGE("s1314 returned code: " << stat << ", returning.");
	    freeIntcrvlist(intcurves, numintcr);
	    freeSurf(sislsf);
	    return segs;
	}
	double* t = sc->et;
	double* c = (sc->ikind==2 || sc->ikind==4)? sc->rcoef : sc->ecoef;
	// Convert the geometric curve to Go format
	SplineCurve* gcv = new SplineCurve(sc->in, sc->ik,
					   t, c, 3,
					   (sc->ikind==2 || sc->ikind==4));
//  	// Debug
//  	gcv->writeStandardHeader(debug);
//  	gcv->write(debug);
//  	// End debug
	
	sc = intcurves[i]->ppar1;
	t = sc->et;
	c = (sc->ikind==2 || sc->ikind==4)? sc->rcoef : sc->ecoef;
	SplineCurve* pcv = new SplineCurve(sc->in, sc->ik,
					   t, c, 2,
					   (sc->ikind==2 || sc->ikind==4));

	vector<SplineCurve*> final_param_curves;
	vector<SplineCurve*> final_space_curves;
	if (bdomain != 0) {
	    // the surface was trimmed.  We must check for intersections with
	    // trimming curves
	    vector<double> params_start_end;
	    bdomain->findPcurveInsideSegments(*pcv, 
					      toptol_.gap, //@ other tolerance here???
					      params_start_end);
	    int num_segments = (int)params_start_end.size() / 2;
	    //cout << "Num segments found: " << num_segments << endl;

	    for (int j = 0; j < num_segments; ++j) {
		SplineCurve* pcv_sub = pcv->subCurve(params_start_end[2 * j],
						     params_start_end[2 * j + 1]);
		final_param_curves.push_back(pcv_sub->clone());
		final_space_curves.push_back(0); //@ change this? (not necessary)
		delete pcv_sub;
	    }
	    // deleting curves that will not be directly used later
	    delete(gcv);
	    delete(pcv);
	} else {
	    final_param_curves.push_back(pcv);
	    final_space_curves.push_back(gcv);
	}

	// pushing back segments
	for (size_t j = 0; j < final_space_curves.size(); ++j) {
	    ftCurveSegment seg(CURVE_INTERSECTION, 
			       JOINT_DISC, 
			       sf,  // underlying surface
			       0,   // second underlying surface (null)
			       shared_ptr<ParamCurve>(final_param_curves[j]),
			       shared_ptr<ParamCurve>(), // second parameter curve (null)
			       shared_ptr<ParamCurve>(final_space_curves[j]));
	    segs.push_back(seg);
	}
    }
    freeIntcrvlist(intcurves, numintcr);
    freeSurf(sislsf);

    return segs;
}


//===========================================================================
void SurfaceModel::booleanIntersect(const ftPlane& plane)
//===========================================================================
{
  // First, we make a list of cells containing the faces
  // Then, check if the cell overlaps the plane, trim the faces within
  // the cell, otherwise check if the faces lies on the positive side of
  // the plane

  double eps = std::min(toptol_.gap, 1.0e-4);

  int nmb_faces = (int)faces_.size();
  for (int ki=0; ki<nmb_faces; ++ki)
    {
	shared_ptr<ParamSurface> sf = getSurface(ki);
	vector<shared_ptr<BoundedSurface> > split_sf;

      // Box test
	if (plane.intersectsBox(sf->boundingBox()))
	  {
	    // Possibility for intersection. Trim
	    split_sf = BoundedUtils::trimWithPlane(sf, plane.point(),
						   plane.normal(), eps);
	  }
	  
	FaceAdjacency<ftEdgeBase,ftFaceBase> adjacency(toptol_);
	if (split_sf.size() == 0)
	  {
	    // No intersection possible or found. Fetch a point in 
	    // the surface and check which side of the plane the
	    // surface lies
	    // First find a point inside the surface
	    RectDomain dom = sf->containingDomain();
	    double upar = 0.5*(dom.umin() + dom.umax());
	    double vpar = 0.5*(dom.vmin() + dom.vmax());
	    Point par = sf->closestInDomain(upar, vpar);
	    
	    Point mid = sf->point(par[0], par[1]);
	    double scpr = plane.normal()*(mid - plane.point());
	    if (scpr >= 0)
	      {
		// The surface lies at the positive side of the plane
		// Remove the surface from the surface set
		adjacency.releaseFaceAdjacency(faces_[ki]);
		faces_.erase(faces_.begin()+ki);
		nmb_faces--;
		ki--;
	      }
	  }
	else
	  {
	    // The surface intersects the plane
	    // Remove the initialsurface from the surface set
	    adjacency.releaseFaceAdjacency(faces_[ki]);
	    faces_.erase(faces_.begin()+ki);
	    nmb_faces--;
	    ki--;

	    // Add new surfaces
	    for (size_t kj=0; kj<split_sf.size(); ++kj)
	      {
		shared_ptr<ftSurface> newSurf;
		newSurf.reset(new ftSurface(split_sf[kj], (int)faces_.size()));
		append(newSurf);
	      }
	  }
    }
}

//===========================================================================
shared_ptr<SurfaceModel> SurfaceModel::trimWithPlane(const ftPlane& plane)
//===========================================================================
{
  // First, we make a list of cells containing the faces
  // Then, check if the cell overlaps the plane, trim the faces within
  // the cell, otherwise check if the faces lies on the positive side of
  // the plane

  vector<shared_ptr<ParamSurface> > inside;
  double eps = std::min(toptol_.gap, 1.0e-4);

  for (size_t ki=0; ki<faces_.size(); ++ki)
    {
	shared_ptr<ParamSurface> sf = getSurface((int)ki);
	vector<shared_ptr<BoundedSurface> > split_sf;

      // Box test
	if (plane.intersectsBox(sf->boundingBox()))
	  {
	    // Possibility for intersection. Trim
	    split_sf = BoundedUtils::trimWithPlane(sf, plane.point(),
						   plane.normal(), eps);
	  }
	  
	if (split_sf.size() == 0)
	  {
	    // No intersection possible or found. Fetch a point in 
	    // the surface and check which side of the plane the
	    // surface lies
	    // First find a point inside the surface
	    RectDomain dom = sf->containingDomain();
	    double upar = 0.5*(dom.umin() + dom.umax());
	    double vpar = 0.5*(dom.vmin() + dom.vmax());
	    Point par = sf->closestInDomain(upar, vpar);
	    
	    Point mid = sf->point(par[0], par[1]);
	    double scpr = plane.normal()*(mid - plane.point());
	    if (scpr < 0.0)
	      {
		// The surface lies at the negative side of the plane
		inside.push_back(sf);
	      }
	  }
	else
	  inside.insert(inside.end(), split_sf.begin(), split_sf.end());
    }
  
  shared_ptr<SurfaceModel> trimmed_model =
    shared_ptr<SurfaceModel>(new SurfaceModel(approxtol_, toptol_.gap,
					      toptol_.neighbour, toptol_.kink,
					      toptol_.bend, inside));

  return trimmed_model;
}


//===========================================================================
// Split surface model by intersection with a different surface model.
// The result is returned in a number of surface models in the following order
// The part of this surface model being inside the other model
// The part of this surface model being outside the other model
// The part of the other surface model being inside this model
// The part of the other surface model being outside this model
// The returned surface models need to to be connected.
// The input surface models are expected to be connected in order to
// get a consistent result
  vector<shared_ptr<SurfaceModel> > 
  SurfaceModel::splitSurfaceModels(shared_ptr<SurfaceModel>& model2)

//===========================================================================
{
  vector<shared_ptr<ParamSurface> > inside1, outside1, inside2, outside2;
  //double eps = std::min(toptol_.gap, 1.0e-4);
  double eps = toptol_.gap;

  // Prepare for storage of intersection curves and bounded surfaces
  int nmb1 = nmbEntities();
  int nmb2 = model2->nmbEntities();
  vector<vector<shared_ptr<CurveOnSurface> > > all_int_cvs1(nmb1);
  vector<vector<shared_ptr<CurveOnSurface> > > all_int_cvs2(nmb2);
  vector<shared_ptr<BoundedSurface> > bd_sfs1(nmb1);
  vector<shared_ptr<BoundedSurface> > bd_sfs2(nmb2);

  // Perform all intersections and store results
  int ki, kj;
  for (ki=0; ki<nmb1; ++ki)
    {
      shared_ptr<ParamSurface> surf1 = faces_[ki]->surface();
      BoundingBox box1 = surf1->boundingBox();
      for (kj=0; kj<nmb2; ++kj)
	{
	  shared_ptr<ParamSurface> surf2 = model2->getSurface(kj);
	  BoundingBox box2 = surf2->boundingBox();

#ifdef DEBUG
	  std::ofstream out("curr_sf_int.g2");
	  surf1->writeStandardHeader(out);
	  surf1->write(out);
	  surf2->writeStandardHeader(out);
	  surf2->write(out);
#endif

	  if (box1.overlaps(box2, eps))
	    {
	      shared_ptr<BoundedSurface> bd1, bd2;
	      vector<shared_ptr<CurveOnSurface> > int_cv1, int_cv2;
	      BoundedUtils::getSurfaceIntersections(surf1, surf2, eps,
						    int_cv1, bd1,
						    int_cv2, bd2);
	      bd_sfs1[ki] = bd1;
	      bd_sfs2[kj] = bd2;
	      if (int_cv1.size() > 0)
		{
		  all_int_cvs1[ki].insert(all_int_cvs1[ki].end(), 
					 int_cv1.begin(), int_cv1.end());
		  all_int_cvs2[kj].insert(all_int_cvs2[kj].end(), 
					 int_cv2.begin(), int_cv2.end());
		}
	    }
	}
    }

#ifdef DEBUG
  std::ofstream of0("intcurves.g2");
  for (ki=0; ki<nmb1; ++ki)
    {
      for (size_t km=0; km<all_int_cvs1[ki].size(); ++km)
	{
	  shared_ptr<ParamCurve> tmpcv = all_int_cvs1[ki][km]->spaceCurve();
	  tmpcv->writeStandardHeader(of0);
	  tmpcv->write(of0);
	}
    }
  for (ki=0; ki<nmb2; ++ki)
    {
      for (size_t km=0; km<all_int_cvs2[ki].size(); ++km)
	{
	  shared_ptr<ParamCurve> tmpcv = all_int_cvs2[ki][km]->spaceCurve();
	  tmpcv->writeStandardHeader(of0);
	  tmpcv->write(of0);
	}
    }
  std::ofstream of01("parcurves.g2");
  for (ki=0; ki<nmb1; ++ki)
    {
      for (size_t km=0; km<all_int_cvs1[ki].size(); ++km)
	{
	  shared_ptr<ParamCurve> tmpcv = all_int_cvs1[ki][km]->parameterCurve();
	  tmpcv->writeStandardHeader(of01);
	  tmpcv->write(of01);
	}
    }
  for (ki=0; ki<nmb2; ++ki)
    {
      for (size_t km=0; km<all_int_cvs2[ki].size(); ++km)
	{
	  shared_ptr<ParamCurve> tmpcv = all_int_cvs2[ki][km]->parameterCurve();
	  tmpcv->writeStandardHeader(of01);
	  tmpcv->write(of01);
	}
    }
#endif

  // Make trimmed surfaces and sort trimmed an non-trimmed surface according
  // to whether they are inside or outside the other surface model
  // First this surface model
  for (ki=0; ki<nmb1; ki++)
    {
      if (all_int_cvs1[ki].size() == 0)
	{
	  // The surface is not involved in any intersections. Check if
	  // it lies inside or outside the other surface model
	  // Fetch a point in the surface
	  shared_ptr<ParamSurface> surf = faces_[ki]->surface();
	  double u, v;
	  Point pnt = surf->getInternalPoint(u,v);

#ifdef DEBUG
	  int state;
	  shared_ptr<BoundedSurface> bdsf = 
	    dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);
	  if (bdsf.get())
	    {
	      bdsf->analyzeLoops();
	      bool valid = bdsf->isValid(state);
	      if (!valid)
		std::cout << "Surface not valid: " << state << std::endl;
	    }
	  std::ofstream of1("curr1.g2");
	  surf->writeStandardHeader(of1);
	  surf->write(of1);
#endif

	  double pt_dist;
 	  bool inside = model2->isInside(pnt, pt_dist);
	  if (inside)
	    {
	      shared_ptr<ParamSurface> tmp_surf = shared_ptr<ParamSurface>(surf->clone());
	      inside1.push_back(tmp_surf);
	      if (fabs(pt_dist) < toptol_.gap)
		outside1.push_back(tmp_surf);
	    }
			      
	  else
      	    outside1.push_back(shared_ptr<ParamSurface>(surf->clone()));
	}
      else
	{
	  // Make bounded surfaces
	  vector<shared_ptr<BoundedSurface> > trim_sfs;
	  try {
	    trim_sfs = 
	      BoundedUtils::splitWithTrimSegments(bd_sfs1[ki], all_int_cvs1[ki],
						  eps);
	  }
	  catch(...)
	    {
	      std::cout << "Trimmed surfaces missing" << std::endl;
	    }
	  for (size_t kr=0; kr<trim_sfs.size(); ++kr)
	    {
#ifdef DEBUG
	      int state;
	      trim_sfs[kr]->analyzeLoops();
	      bool valid = trim_sfs[kr]->isValid(state);
	      if (!valid)
		std::cout << "Surface not valid: " << state << std::endl;
#endif

	  // Check if the trimmed surface lies inside or outside the 
	  // other surface model.
	      double u, v;
	      Point pnt =  trim_sfs[kr]->getInternalPoint(u,v);

#ifdef DEBUG
	      std::ofstream of1("curr1.g2");
	      trim_sfs[kr]->writeStandardHeader(of1);
	      trim_sfs[kr]->write(of1);
#endif

	      double pt_dist;
	      bool inside = model2->isInside(pnt, pt_dist);
	      if (inside)
		{
		  inside1.push_back(trim_sfs[kr]);
		  if (fabs(pt_dist) < toptol_.gap)
		    outside1.push_back(shared_ptr<ParamSurface>(trim_sfs[kr]->clone()));
		}
			      
	      else
		outside1.push_back(trim_sfs[kr]);
	    }
	}
    }
 
  // The other surface model
  for (ki=0; ki<nmb2; ki++)
    {
      if (all_int_cvs2[ki].size() == 0)
	{
	  // The surface is not involved in any intersections. Check if
	  // it lies inside or outside the other surface model
	  // Fetch a point in the surface
	  shared_ptr<ParamSurface> surf = model2->getSurface(ki);
	  double u, v;
	  Point pnt = surf->getInternalPoint(u,v);

#ifdef DEBUG
	  int state;
	  shared_ptr<BoundedSurface> bdsf = 
	    dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);
	  if (bdsf.get())
	    {
	      bdsf->analyzeLoops();
	      bool valid = bdsf->isValid(state);
	      if (!valid)
		std::cout << "Surface not valid: " << state << std::endl;
	    }

	  std::ofstream of1("curr2.g2");
	  surf->writeStandardHeader(of1);
	  surf->write(of1);
#endif

	  double pt_dist;
	  bool inside = isInside(pnt, pt_dist);
	  if (inside)
	    {
	      shared_ptr<ParamSurface> tmp_surf = shared_ptr<ParamSurface>(surf->clone());
	      inside2.push_back(tmp_surf);
	      if (fabs(pt_dist) < toptol_.gap)
		outside2.push_back(tmp_surf);
	    }
	  else
      	    outside2.push_back(shared_ptr<ParamSurface>(surf->clone()));
	}
      else
	{
	  // Make bounded surfaces
	  vector<shared_ptr<BoundedSurface> > trim_sfs;
	  try {
	    trim_sfs = 
	      BoundedUtils::splitWithTrimSegments(bd_sfs2[ki], all_int_cvs2[ki],
						  eps);
	  }
	  catch(...)
	    {
	      std::cout << "Trimmed surfaces missing" << std::endl;
	    }
	  for (size_t kr=0; kr<trim_sfs.size(); ++kr)
	    {
#ifdef DEBUG
	      int state;
	      trim_sfs[kr]->analyzeLoops();
	      bool valid = trim_sfs[kr]->isValid(state);
	      if (!valid)
		std::cout << "Surface not valid: " << state << std::endl;
#endif

	  // Check if the trimmed surface lies inside or outside the 
	  // other surface model.
	      double u, v;
	      Point pnt =  trim_sfs[kr]->getInternalPoint(u,v);

#ifdef DEBUG
	      std::ofstream of1("curr2.g2");
	      trim_sfs[kr]->writeStandardHeader(of1);
	      trim_sfs[kr]->write(of1);
#endif

	      double pt_dist;
	      bool inside = isInside(pnt, pt_dist);
	      if (inside)
		{
		  inside2.push_back(trim_sfs[kr]);
		  if (fabs(pt_dist) < toptol_.gap)
		    outside2.push_back(shared_ptr<ParamSurface>(trim_sfs[kr]->clone()));
		}			      
	      else
		outside2.push_back(trim_sfs[kr]);
	    }
	}
    }
 
  vector<shared_ptr<SurfaceModel> > split_models(4);
  if (inside1.size() > 0)
    split_models[0] = 
      shared_ptr<SurfaceModel>(new SurfaceModel(approxtol_, toptol_.gap,
						toptol_.neighbour,
						toptol_.kink, toptol_.bend,
						inside1));
  if (outside1.size() > 0)
    split_models[1] = 
      shared_ptr<SurfaceModel>(new SurfaceModel(approxtol_, toptol_.gap,
						toptol_.neighbour,
						toptol_.kink, toptol_.bend,
						outside1));
  if (inside2.size() > 0)
    split_models[2] = 
      shared_ptr<SurfaceModel>(new SurfaceModel(approxtol_, toptol_.gap,
						toptol_.neighbour,
						toptol_.kink, toptol_.bend,
						inside2));
  if (outside2.size() > 0)
    split_models[3] = 
      shared_ptr<SurfaceModel>(new SurfaceModel(approxtol_, toptol_.gap,
						toptol_.neighbour,
						toptol_.kink, toptol_.bend,
						outside2));
  return split_models;
}

//===========================================================================
// Check if a given point lies on the positive or negative side of this
// SurfaceModel with respect to the model normal.
// NB! If the surface set is open, this function requires a consistent
// normal behaviour for all surfaces
bool SurfaceModel::isInside(const Point& pnt, double& dist) 
//===========================================================================
    {
      // Check if this surface set belongs to a solid. In that case check
      // if the point lies inside this solid
      if (faces_.size() == 0)
	return false;
      
      dist = -1;  // Initiate to no information
      ftSurface *curr = faces_[0]->asFtSurface();
      Point pnt1 = pnt;
      Point clo_pnt;
      int idx;
      double par[2];
      
      closestPoint(pnt1, clo_pnt, idx, par, dist);

      if (dist < toptol_.gap)
	return true;  // On boundary
      else if (curr && curr->hasBody())
	{
	  return curr->getBody()->isInside(pnt);
	}
      else
	{
	  // Current simple solution
	  // Check surface normal
	  Point normal = faces_[idx]->normal(par[0], par[1]);
	  Point vec = pnt - clo_pnt;
	  if (normal*vec < 0.0 || vec.length() < toptol_.gap)
	    return true;
	  else 
	    return false;
	}

	
    }

//===========================================================================
  // Intersection with a line. Expected output is points, probably one point. Curves 
  // can occur in special configurations
     shared_ptr<IntResultsModel> SurfaceModel::intersect(const ftLine& line)
//===========================================================================
{
  shared_ptr<IntResultsSfModel> intersections = 
    shared_ptr<IntResultsSfModel>(new IntResultsSfModel(this,
							line)); // Empty storage for output
  ftCurve int_curves;
  vector<ftPoint> int_points;

  intersect(line, int_curves, int_points);

  intersections->addIntPts(int_points);
  intersections->addIntCvs(int_curves);
  
  return intersections;
}



// Intersection with a line. Expected output is points, probably one point. Curves 
// can occur in special configurations
void 
SurfaceModel::intersect(const ftLine& line, // New class, consist of one point and one direction
			// represented by Go::Point. Just storage, not much content (yet)
			ftCurve& int_curves, // Intersection curves, one curve may
			// cross several of the surfaces in the model, but each curve is connected
			// and simple.
			std::vector<ftPoint>& int_points)  // Found intersection points
//===========================================================================
{
  // First, we make a list of cells that intersect the line
  // Then, run intersection on each of the surfaces touching those cells,
  // if the bounding boxes overlap.

  vector<ftPoint> result;
  vector<ftCurveSegment> line_segments;

  initializeCelldiv();
  if (!line.intersectsBox(celldiv_ -> big_box())) 
      return;

  highest_face_checked_ = 0;
  int i, j;
  for (i = 0; i < celldiv_ -> numCells(); ++i) {
    ftCell aCell = celldiv_ -> getCell(i);
    if (line.intersectsBox(aCell.box()))
      {
      for (j = 0; j < aCell.num_faces(); ++j) {
	ftSurface* face = aCell.face(j) -> asFtSurface();
	int id = face->getId();
	if (!face_checked_[id]) {
	  face_checked_[id] = true;
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
	  highest_face_checked_ = max(highest_face_checked_, id);
#else
	  highest_face_checked_ = 
	    std::max(highest_face_checked_, id);
#endif
	  if (line.intersectsBox(face->boundingBox()))
	    localIntersect(line, face, result, line_segments);
	}
      }
    }
  }
  std::fill(face_checked_.begin(),
	    face_checked_.begin() + highest_face_checked_ + 1, false);
  
    // We have to connect any curves that should connect
  int num_curves = (int)line_segments.size();
  vector<ftCurve> curves;
  vector<ftEdgePointInfo> epinfo = MakeEdgePointInfo(line_segments.begin(),
						     line_segments.end(),
						     toptol_.neighbour);
  for (i = 0; i < num_curves; ++i) {
      curves.push_back(ftCurve(CURVE_INTERSECTION));
      curves[i].appendSegment(line_segments[i]);
  }

  Point pts[2][2];
  int k1, k2;
  vector<bool> joined(num_curves, false);
  for (i = 0; i < num_curves - 1; ++i) {
      for (j = i + 1; j < num_curves; ++j) {
	  // Let's compare the endpoints of curves[i] and curves[j]
	  int seg = 0;
	  curves[i].point(curves[i].startOfSegment(seg), seg, pts[0][0]);
	  seg = curves[i].numSegments() - 1;
	  curves[i].point(curves[i].endOfSegment(seg), seg, pts[0][1]);
	  seg = 0;
	  curves[j].point(curves[j].startOfSegment(seg), seg, pts[1][0]);
	  seg = curves[j].numSegments() - 1;
	  curves[j].point(curves[j].endOfSegment(seg), seg, pts[1][1]);
	  for (k1 = 0; k1 < 2; ++k1) {
	      for (k2 = 0; k2 < 2; ++k2) {
		  if (pts[0][k1].dist(pts[1][k2]) < toptol_.neighbour) {
		      joined[i] = true;
		      if (k2 == 0) {
			  // Prepend
			  if (k1 == 0)
			      curves[i].reverse();
			  curves[j].prependCurve(curves[i]);
		      } else {
			  // Append
			  if (k1 == 1)
			      curves[i].reverse();
			  curves[j].appendCurve(curves[i]);
		      }
		  }
	      }
	    }
	    if (joined[i])
		break; // Out of the j-loop to start on the next i-value
	}
    }


    // We should now build intcurve from curves[]
    ftCurve intcurve(CURVE_INTERSECTION);
    for (i = 0; i < (int)curves.size(); ++i)
	if (!joined[i])
	    intcurve += curves[i];

    intcurve.joinSegments(toptol_.gap, toptol_.neighbour, toptol_.kink, toptol_.bend);

    int_curves = intcurve;
    int_points = result;
      
}

//===========================================================================
vector<ftPoint> SurfaceModel::intersect(const ftLine& line, 
					vector<bool>& represent_segment) 
//===========================================================================
{
  // First, we make a list of cells that intersect the line
  // Then, run intersection on each of the surfaces touching those cells,
  // if the bounding boxes overlap.

  vector<ftPoint> result;
  vector<ftCurveSegment> line_segments;

  initializeCelldiv();
  if (!line.intersectsBox(celldiv_ -> big_box())) return result;

  highest_face_checked_ = 0;
  for (int i = 0; i < celldiv_ -> numCells(); ++i) {
    ftCell aCell = celldiv_ -> getCell(i);
    if (line.intersectsBox(aCell.box()))
      {
      for (int j = 0; j < aCell.num_faces(); ++j) {
	ftSurface* face = aCell.face(j) -> asFtSurface();
	int id = face->getId();
	if (!face_checked_[id]) {
	  face_checked_[id] = true;
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
	  highest_face_checked_ = max(highest_face_checked_, id);
#else
	  highest_face_checked_ = 
	    std::max(highest_face_checked_, id);
#endif
	  if (line.intersectsBox(face->boundingBox()))
	    localIntersect(line, face, result, line_segments);
	}
      }
    }
  }
  std::fill(face_checked_.begin(),
	    face_checked_.begin() + highest_face_checked_ + 1, false);
  
  size_t kr;
  for (kr=0; kr<result.size(); kr++)
      represent_segment.push_back(false);
  for (kr=0; kr<line_segments.size(); kr++)
  {
      Point pt_2D;
      double start = line_segments[kr].startOfSegment();
      double end = line_segments[kr].endOfSegment();
      line_segments[kr].paramcurvePoint(0, 0.5*(start+end), pt_2D);
      
      Point pt_3D;
      line_segments[kr].point(0.5*(start+end), pt_3D);
      result.push_back(ftPoint(pt_3D, line_segments[kr].face(0)->asFtSurface(), pt_2D[0], pt_2D[1]));

      represent_segment.push_back(true);
  }
      
  return result;
}




//===========================================================================
vector<pair<ftPoint, double> > 
SurfaceModel::intersect(shared_ptr<SplineCurve> crv,
			vector<bool>& represent_segment) 
//===========================================================================
{
  // First, we make a list of cells that intersect the line
  // Then, run intersection on each of the surfaces touching those cells,
  // if the bounding boxes overlap.

  vector<pair<ftPoint, double> > result;
  vector<ftCurveSegment> crv_segments;
  vector<pair<double, double> > segment_bound;

  initializeCelldiv();
  BoundingBox cv_box = crv->boundingBox();

  // Check if an intersection is possible
  if (!cv_box.overlaps(celldiv_->big_box(), toptol_.gap)) 
    return result;

  highest_face_checked_ = 0;
  for (int i = 0; i < celldiv_ -> numCells(); ++i) {
    ftCell aCell = celldiv_ -> getCell(i);
    if (cv_box.overlaps(aCell.box(), toptol_.gap))
      {
      for (int j = 0; j < aCell.num_faces(); ++j) {
	ftSurface* face = aCell.face(j) -> asFtSurface();
	int id = face->getId();
	if (!face_checked_[id]) {
	  face_checked_[id] = true;
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
	  highest_face_checked_ = max(highest_face_checked_, id);
#else
	  highest_face_checked_ = 
	    std::max(highest_face_checked_, id);
#endif
	  if (cv_box.overlaps(face->boundingBox(), toptol_.gap))
	    localIntersect(crv, face, result, crv_segments, segment_bound);
	}
      }
    }
  }
  std::fill(face_checked_.begin(),
	    face_checked_.begin() + highest_face_checked_ + 1, false);
  
  size_t kr;
  for (kr=0; kr<result.size(); kr++)
      represent_segment.push_back(false);
  for (kr=0; kr<crv_segments.size(); kr++)
  {
    // Evaluate crv
    double tpar = 0.5*(segment_bound[kr].first + segment_bound[kr].second);
    Point pos = crv->ParamCurve::point(tpar);
    Point seed;
    double start = crv_segments[kr].startOfSegment();
    double end = crv_segments[kr].endOfSegment();
    crv_segments[kr].paramcurvePoint(0, 0.5*(start+end), seed);
    ftSurface* face = crv_segments[kr].face(0)->asFtSurface();

    // Closest point in face
    double upar, vpar, dist;
    Point face_pos;
    face->surface()->closestPoint(pos, upar, vpar, face_pos, dist,
				  toptol_.gap, NULL, seed.begin());
      
    result.push_back(make_pair(ftPoint(face_pos, face, upar, vpar), tpar));

    represent_segment.push_back(true);
  }
      
  return result;
}

//===========================================================================
bool SurfaceModel::hit(const Point& point, const Point& dir, ftPoint& result) 
//===========================================================================
{
  // Fetch the closest point to the given input point of the intersections
  // between this surface model and the specified line, if any

  // First, we make a list of cells that intersect the line
  // Then, run intersection on each of the surfaces touching those cells,
  // if the bounding boxes overlap.

  bool hit = false;
  vector<ftPoint> current;
  vector<ftCurveSegment> line_segments;
  initializeCelldiv();
  BoundingBox box = celldiv_->big_box();
  ftLine line(dir, point);  // Represent beam as line
  if (!line.intersectsBox(box)) 
    return false;

  Point mid = 0.5*(box.low() + box.high());  // Midpoint in the box
  double rad = mid.dist(box.low());          // Radius in surronding sphere
  double min_dist = point.dist(mid) + rad;   // A long distance

  highest_face_checked_ = 0;
  for (int i = 0; i < celldiv_ -> numCells(); ++i) 
    {
      ftCell aCell = celldiv_ -> getCell(i);
      BoundingBox cell_box = aCell.box();
      if (line.intersectsBox(cell_box))
	{
	  Point cell_mid = 0.5*(cell_box.low()+cell_box.high());
	  double cell_dist = point.dist(cell_mid) - cell_mid.dist(cell_box.low());
	  if (cell_dist > min_dist)
	    continue;  // No minimum distance can be found

	  for (int j = 0; j < aCell.num_faces(); ++j) {
	    ftSurface* face = aCell.face(j) -> asFtSurface();
	    int id = face->getId();
	    if (!face_checked_[id]) 
	      {
		face_checked_[id] = true;
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
		highest_face_checked_ = max(highest_face_checked_, id);
#else
		highest_face_checked_ = 
		  std::max(highest_face_checked_, id);
#endif
		BoundingBox face_box = face->boundingBox();
		if (line.intersectsBox(face_box))
		  {
		    Point face_mid = 0.5*(face_box.low()+face_box.high());
		    double face_dist = point.dist(face_mid) - face_mid.dist(face_box.low());
		    if (face_dist > min_dist)
		      continue;  // No minimum distance can be found
		    localIntersect(line, face, current, line_segments);

		    // Find closest intersction and update smallest distance
		    size_t kd;
		    for (kd=0; kd<current.size(); ++kd)
		      {
			Point pos = current[kd].position();

			// Make sure that the point is on the correct side
			// of the point on line
			if (dir*(pos - point) < -toptol_.gap)
			  continue;

			hit = true;
			double dist = point.dist(pos);
			if (dist < min_dist)
			  {
			    result = current[kd];
			    min_dist = dist;
			  }
		      }
		    for (kd=0; kd<line_segments.size(); ++kd)
		      {
			hit = true;
			Point pos = line_segments[kd].startPoint();
			double dist = point.dist(pos);
			if (dist < min_dist)
			  {
			    Point param; 
			    line_segments[kd].paramcurvePoint(0, line_segments[kd].startOfSegment(), 
							      param);
			    result = ftPoint(pos, current[kd].face()->asFtSurface(), 
					     param[0], param[1]);
			    min_dist = dist;
			  }
			pos = line_segments[kd].endPoint();
			dist = point.dist(pos);
			if (dist < min_dist)
			  {
			    Point param; 
			    line_segments[kd].paramcurvePoint(0, line_segments[kd].endOfSegment(), 
							      param);
			    result = ftPoint(pos, current[kd].face()->asFtSurface(), 
					     param[0], param[1]);
			    min_dist = dist;
			  }
		      }
		  }
	      }
	  }
	}
    }
      
  return hit;
}




//===========================================================================
void SurfaceModel::localIntersect(const ftLine& line,
				  ftSurface* sf,
				  vector<ftPoint>& result,
				  vector<ftCurveSegment>& line_segments) const
//===========================================================================
{
  // Convert the surface to a SISLSurf in order to use SISL functions
  // on it. The "false" argument dictates that the SISLSurf will only    
  // copy pointers to arrays, not the arrays themselves.
    shared_ptr<ParamSurface> psurf = sf->surface();
    shared_ptr<SplineSurface> surf;
    if (psurf->instanceType() == Class_SplineSurface)
	surf = dynamic_pointer_cast<SplineSurface, ParamSurface>(psurf);
    else
    {
	shared_ptr<BoundedSurface> bsurf = dynamic_pointer_cast<BoundedSurface, ParamSurface>(psurf);
	if (bsurf.get())
	    surf = dynamic_pointer_cast<SplineSurface, ParamSurface>(bsurf->underlyingSurface());
    }
  SplineSurface* splinesf
    = dynamic_cast<SplineSurface*>(sf->surface().get());
  const CurveBoundedDomain* bdomain = 0;
  if (splinesf == 0) {
    BoundedSurface* bs
      = dynamic_cast<BoundedSurface*>(sf->surface().get());
    if (bs != 0) {
      splinesf
	= dynamic_cast<SplineSurface*>(bs->underlyingSurface().get());
      bdomain = &(bs->parameterDomain());
    } else {
      THROW("You cannot intersect this class of surface: "
	    << sf->surface()->instanceType());
    }
  }
  ASSERT(splinesf != 0);

  SISLSurf* sislsf = GoSurf2SISL(*splinesf, false);
  int dim = 3;
  double epsco = 1e-15; // Not used
  double epsge = 1e-6;
  int numintpt;  // number of single intersection points
  double* pointpar = 0; // array containing the parameter values of single intersect. pt.
  int numintcr; // number of intersection curves
  SISLIntcurve** intcurves = 0;
  int stat;
  Point pnt = line.point();
  Point dir = line.direction();
  // Find the intersection points
  s1856(sislsf, pnt.begin(), dir.begin(), dim, epsco, epsge,
	&numintpt, &pointpar, &numintcr, &intcurves, &stat);
  MESSAGE_IF(stat!=0, "s1856 returned code: " << stat);

  int i;
  for (i = 0; i < numintpt; ++i)
    {
      double u = pointpar [i<<1];
      double v = pointpar [i<<1 | 1];
      Point pt = sf -> point(u, v);

      bool in_domain = true;
      if (bdomain != 0)
	{
	  // Check if the point is inside the trimmed surface
	  Array<double,2> tmp_pt(u,v);
	  in_domain = bdomain->isInDomain(tmp_pt, epsge);
	  
	}
      if (in_domain)
	result.push_back(ftPoint(pt, sf, u, v));
    }

  for (i=0; i<numintcr; i++)
  {
      // Evaluate endpoints of line segment and make geometry curve
      int npt = intcurves[i]->ipoint;
      Point pt1 = sf->point(intcurves[i]->epar1[0],intcurves[i]->epar1[1]);
      Point pt2 = sf->point(intcurves[i]->epar1[2*(npt-1)],intcurves[i]->epar1[2*(npt-1)+1]);
      SplineCurve *gcv = new SplineCurve(pt1, pt2);

      // Project the curve into the parameter space of the surface
      shared_ptr<Point> pt1_2D = shared_ptr<Point>(new Point(intcurves[i]->epar1[0],intcurves[i]->epar1[1]));
      shared_ptr<Point> pt2_2D = shared_ptr<Point>(new Point(intcurves[i]->epar1[2*(npt-1)],intcurves[i]->epar1[2*(npt-1)+1]));
      shared_ptr<ParamCurve> gcv2 = shared_ptr<ParamCurve>(gcv->clone());
      double tol = approxtol_;
      SplineCurve *pcv = CurveCreators::projectSpaceCurve(gcv2, psurf, 
							  pt1_2D, pt2_2D, tol);
      
 	vector<SplineCurve*> final_param_curves;
	vector<SplineCurve*> final_space_curves;
	if (bdomain != 0) {
	    // the surface was trimmed.  We must check for intersections with
	    // trimming curves
	    vector<double> params_start_end;
	    bdomain->findPcurveInsideSegments(*pcv, 
					      toptol_.gap, //@ other tolerance here???
					      params_start_end);
	    int num_segments = (int)params_start_end.size() / 2;
	    //cout << "Num segments found: " << num_segments << endl;

	    for (int j = 0; j < num_segments; ++j) {
		SplineCurve* pcv_sub = pcv->subCurve(params_start_end[2 * j],
						     params_start_end[2 * j + 1]);
		final_param_curves.push_back(pcv_sub->clone());
		final_space_curves.push_back(0); //@ change this? (not necessary)
		delete pcv_sub;
	    }
	    // deleting curves that will not be directly used later
	    delete(gcv);
	    delete(pcv);
	} else {
	    final_param_curves.push_back(pcv);
	    final_space_curves.push_back(gcv);
	}

	// pushing back segments
	for (size_t j = 0; j < final_space_curves.size(); ++j) {
	    ftCurveSegment seg(CURVE_INTERSECTION, 
			       JOINT_DISC, 
			       sf,  // underlying surface
			       0,   // second underlying surface (null)
			       shared_ptr<ParamCurve>(final_param_curves[j]),
			       shared_ptr<ParamCurve>(), // second parameter curve (null)
			       shared_ptr<ParamCurve>(final_space_curves[j]));
	    line_segments.push_back(seg);
	}
    }

  free(pointpar);
  freeIntcrvlist(intcurves, numintcr);
  freeSurf(sislsf);
}

//===========================================================================
void SurfaceModel::localIntersect(shared_ptr<SplineCurve> crv,
				  ftSurface* sf,
				  vector<pair<ftPoint, double> >& result,
				  vector<ftCurveSegment>& crv_segments,
				  vector<pair<double,double> >& crv_bound) const
//===========================================================================
{
  // Convert the surface to a SISLSurf in order to use SISL functions
  // on it. The "false" argument dictates that the SISLSurf will only    
  // copy pointers to arrays, not the arrays themselves.
    shared_ptr<ParamSurface> psurf = sf->surface();
    shared_ptr<SplineSurface> surf;
    if (psurf->instanceType() == Class_SplineSurface)
	surf = dynamic_pointer_cast<SplineSurface, ParamSurface>(psurf);
    else
    {
	shared_ptr<BoundedSurface> bsurf = dynamic_pointer_cast<BoundedSurface, ParamSurface>(psurf);
	if (bsurf.get())
	    surf = dynamic_pointer_cast<SplineSurface, ParamSurface>(bsurf->underlyingSurface());
    }
  SplineSurface* splinesf
    = dynamic_cast<SplineSurface*>(sf->surface().get());
  const CurveBoundedDomain* bdomain = 0;
  if (splinesf == 0) {
    BoundedSurface* bs
      = dynamic_cast<BoundedSurface*>(sf->surface().get());
    if (bs != 0) {
      splinesf
	= dynamic_cast<SplineSurface*>(bs->underlyingSurface().get());
      bdomain = &(bs->parameterDomain());
    } 
      //else {
//       THROW("You cannot intersect this class of surface: "
// 	    << sf->surface()->instanceType());
//     }
  }
  if (splinesf == 0)
    splinesf = psurf->asSplineSurface();

  ASSERT(splinesf != 0);

  SISLSurf* sislsf = GoSurf2SISL(*splinesf, false);
  SISLCurve* sislcv = Curve2SISL(*crv, false);
  // int dim = 3;
  double epsco = 1e-15; // Not used
  double epsge = 1e-6;
  int numintpt;  // number of single intersection points
  double* pointpar1 = 0; // par. values of single intersect. pt. for surface
  double* pointpar2 = 0; // par. values of single intersect. pt. for curve
  int numintcr; // number of intersection curves
  SISLIntcurve** intcurves = 0;
  int stat;
  // Find the intersection points
  s1858(sislsf, sislcv, epsco, epsge, &numintpt, 
	&pointpar1, &pointpar2, &numintcr, &intcurves, &stat);
  MESSAGE_IF(stat!=0, "s1856 returned code: " << stat);

  int i;
  for (i = 0; i < numintpt; ++i)
    {
      double u = pointpar1 [i<<1];
      double v = pointpar1 [i<<1 | 1];
      Point pt = sf -> point(u, v);

      bool in_domain = true;
      if (bdomain != 0)
	{
	  // Check if the point is inside the trimmed surface
	  Array<double,2> tmp_pt(u,v);
	  try {
	    in_domain = bdomain->isInDomain(tmp_pt, epsge);
	  }
	  catch (...)
	    {
#ifdef DEBUG
	      std::ofstream of("domain_sf.g2");
	      psurf->writeStandardHeader(of);
	      psurf->write(of);
	      crv->writeStandardHeader(of);
	      crv->write(of);
#endif
	      THROW("Error in domain check");
	    }
	  
	}
      if (in_domain)
	result.push_back(make_pair(ftPoint(pt, sf, u, v), pointpar2[i]));
    }

  for (i=0; i<numintcr; i++)
  {
      // Fetch curve piece
      int npt = intcurves[i]->ipoint;
      double t1 = intcurves[i]->epar2[0];
      double t2 = intcurves[i]->epar2[npt-1];
      shared_ptr<ParamCurve> gcv = 
	shared_ptr<ParamCurve>(crv->subCurve(std::min(t1,t2),
					     std::max(t1,t2)));
      // Project the curve into the parameter space of the surface
      shared_ptr<Point> pt1_2D = shared_ptr<Point>(new Point(intcurves[i]->epar1[0],intcurves[i]->epar1[1]));
      shared_ptr<Point> pt2_2D = shared_ptr<Point>(new Point(intcurves[i]->epar1[2*(npt-1)],intcurves[i]->epar1[2*(npt-1)+1]));
      if (t1 > t2)
// 	{
// 	  gcv->reverseParameterDirection();
 	  std::swap(pt1_2D, pt2_2D);
// 	}

      double tol = approxtol_;
      SplineCurve *pcv = CurveCreators::projectSpaceCurve(gcv, psurf, 
							  pt1_2D, pt2_2D, tol);
      if (pcv == NULL)
	{
#ifdef DEBUG
	  std::ofstream of("sf_cv.g2");
	  psurf->writeStandardHeader(of);
	  psurf->write(of);
	  gcv->writeStandardHeader(of);
	  gcv->write(of);
	  std::cout << "localIntersect, no pcurve computed" << std::endl;
#endif
	  continue;  // Curve outside of trimmed surface
	}
      
      vector<SplineCurve*> final_param_curves;
      vector<SplineCurve*> final_space_curves;
      if (bdomain != 0) {
	// the surface was trimmed.  We must check for intersections with
	// trimming curves
	vector<double> params_start_end;
	bdomain->findPcurveInsideSegments(*pcv, 
					  toptol_.gap, //@ other tolerance here???
					  params_start_end);
	int num_segments = (int)params_start_end.size() / 2;
	//cout << "Num segments found: " << num_segments << endl;
	
	for (int j = 0; j < num_segments; ++j) {
	  crv_bound.push_back(make_pair(params_start_end[2 * j],
					params_start_end[2 * j + 1]));
	  shared_ptr<ParamCurve> pcv_sub = 
	    shared_ptr<ParamCurve>(pcv->subCurve(params_start_end[2 * j],
						 params_start_end[2 * j + 1]));
	  shared_ptr<ParamCurve> gcv_sub = 
	    shared_ptr<ParamCurve>(crv->subCurve(params_start_end[2 * j],
						 params_start_end[2 * j + 1]));	  
	  ftCurveSegment seg(CURVE_INTERSECTION, 
			     JOINT_DISC, 
			     sf,  // underlying surface
			     0,   // second underlying surface (null)
			     pcv_sub,
			     shared_ptr<ParamCurve>(), // second parameter curve (null)
			     gcv_sub); 
	  crv_segments.push_back(seg);
	}
      } else {
	crv_bound.push_back(make_pair(gcv->startparam(), gcv->endparam()));
	ftCurveSegment seg(CURVE_INTERSECTION, 
			   JOINT_DISC, 
			   sf,  // underlying surface
			   0,   // second underlying surface (null)
			   shared_ptr<ParamCurve>(pcv),
			   shared_ptr<ParamCurve>(), // second parameter curve (null)
			   gcv);
	crv_segments.push_back(seg);
      }
    }

  free(pointpar1);
  free(pointpar2);
  freeIntcrvlist(intcurves, numintcr);
  freeSurf(sislsf);
  freeCurve(sislcv);
}

//===========================================================================
void 
SurfaceModel::getOverlappingEdges(double tol,
				  vector<pair<shared_ptr<ftEdgeBase>, shared_ptr<ftEdgeBase> > >& edges)
//===========================================================================
{
      // Check if there are any faces. @jbt
      if (faces_.empty()) {
	  MESSAGE("No faces - return empty CellDivision object.");
	  celldiv_ = shared_ptr<CellDivision>();
	  return;
      }

    // First look for overlapping cells
    if (!celldiv_.get())
	initializeCelldiv();

    int num_cell = celldiv_->numCells();
    int ki, kj, kr, kh;
    highest_face_checked_ = 0;
    vector<ftSurface*> checked1;
    vector<ftSurface*> checked2;
    for (ki=0; ki<num_cell; ++ki)
    {
	ftCell cell1 = celldiv_ -> getCell(ki);
	int num_face1 = cell1.num_faces();
	if (num_face1 == 0)
	    continue;
	for (kj=ki; kj<num_cell; ++kj)
	{
	    // Check cell overlap
	    ftCell cell2 = celldiv_ -> getCell(kj);
	    int num_face2 = cell2.num_faces();
	    if (num_face2 == 0)
		continue;

	    if (!cell1.box().overlaps(cell2.box(), tol))
		continue;

	    // Find edge overlaps within these cells. First look for faces
	    for (kr=0; kr<num_face1; ++kr)
	    {
		ftSurface *face1 = cell1.face(kr);
		int id = face1->getId();
		if (face_checked_[id])
		    continue;

		face_checked_[id] = true;
		highest_face_checked_ = 
		    std::max(highest_face_checked_, id);
		int first = (ki == kj) ? kr+1 : 0;
		for (kh=first; kh<num_face2; ++kh)
		{
		    ftSurface *face2 = cell2.face(kh);
		    if (face1 == face2)
			continue;

		    // Check face overlap
		    if (!cell1.faceBox(kr).overlaps(cell2.faceBox(kh), tol))
			continue;

		    size_t j1;
		    for (j1=0; j1<checked1.size(); ++j1)
			if ((checked1[j1] == face1 && checked2[j1] == face2) ||
			    (checked1[j1] == face2 && checked2[j1] == face1))
			    break;
		    if (j1 < checked1.size())
			continue;

		    checked1.push_back(face1);
		    checked2.push_back(face2);

		    // Check edge overlap
		    vector<shared_ptr<ftEdgeBase> > edges1 = 
			face1->createInitialEdges();
		    vector<shared_ptr<ftEdgeBase> > edges2 = 
			face2->createInitialEdges();

		    size_t i1, i2;
		    for (i1=0; i1<edges1.size(); ++i1)
		    {
			BoundingBox edgebox1 = edges1[i1]->geomEdge()->geomCurve()->boundingBox();
			for (i2=0; i2<edges2.size(); ++i2)
			{
			    if (edges1[i1]->twin() && edges1[i1]->twin() == edges2[i2].get())
				continue;

			    BoundingBox edgebox2 = edges2[i2]->geomEdge()->geomCurve()->boundingBox();
			    if (edgebox1.overlaps(edgebox2, tol))
			    {
				// A candidate is found
				edges.push_back(make_pair(edges1[i1],edges2[i2]));
			    }
			}
		    }
		}
	    }
	}
    }

}

//===========================================================================
void 
SurfaceModel::getOverlappingFaces(double tol,
				  vector<pair<ftSurface*, ftSurface*> >& faces)
//===========================================================================
{
    // First look for overlapping cells
    if (!celldiv_.get())
	initializeCelldiv();

    int num_cell = celldiv_->numCells();
    int ki, kj, kr, kh;
    highest_face_checked_ = 0;
    vector<ftSurface*> checked1;
    vector<ftSurface*> checked2;
    for (ki=0; ki<num_cell; ++ki)
    {
	ftCell cell1 = celldiv_ -> getCell(ki);
	int num_face1 = cell1.num_faces();
	if (num_face1 == 0)
	    continue;
	for (kj=ki; kj<num_cell; ++kj)
	{
	    // Check cell overlap
	    ftCell cell2 = celldiv_ -> getCell(kj);
	    int num_face2 = cell2.num_faces();
	    if (num_face2 == 0)
		continue;

	    if (!cell1.box().overlaps(cell2.box(), tol))
		continue;

	    // Find face overlaps within these cells.
	    for (kr=0; kr<num_face1; ++kr)
	    {
		ftSurface *face1 = cell1.face(kr);
		int id = face1->getId();
		if (face_checked_[id])
		    continue;

		face_checked_[id] = true;
		highest_face_checked_ = 
		    std::max(highest_face_checked_, id);
		int first = (ki == kj) ? kr+1 : 0;
		for (kh=first; kh<num_face2; ++kh)
		{
		    ftSurface *face2 = cell2.face(kh);
		    if (face1 == face2)
			continue;

		    // Check face overlap
		    if (!cell1.faceBox(kr).overlaps(cell2.faceBox(kh), tol))
			continue;

		    size_t j1;
		    for (j1=0; j1<checked1.size(); ++j1)
			if ((checked1[j1] == face1 && checked2[j1] == face2) ||
			    (checked1[j1] == face2 && checked2[j1] == face1))
			    break;
		    if (j1 < checked1.size())
			continue;

		    checked1.push_back(face1);
		    checked2.push_back(face2);

		    faces.push_back(make_pair(face1, face2));
		}
	    }
	}
    }
}

//===========================================================================
void 
SurfaceModel::extremalPoint(Point& dir, Point& ext_pnt, int& idx, 
			    double ext_par[]) 
//===========================================================================
{
  // Division of surface set into cells to speed up the compuations
  if (!celldiv_.get())
    initializeCelldiv();

  highest_face_checked_ = 0;  // Keep track on the faces checked already
  idx = -1;                   // No candidate found so far
  bool possible;
  for (int ki = 0; ki < celldiv_ -> numCells(); ++ki) 
    {
      ftCell aCell = celldiv_ -> getCell(ki);

      if (idx >= 0)
	{
	  // Check if the current cell can provide a more extreme point
	  possible = boxExtreme(aCell.box(), dir, ext_pnt);
	}
      else
	possible = true;

      if (possible)
	{
	  for (int kj=0; kj<aCell.num_faces(); ++kj)
	    {
	      ftSurface* face = aCell.face(kj)->asFtSurface();
	      int id = face->getId();

	      if (!face_checked_[id])
		{
		  face_checked_[id] = true;
		  highest_face_checked_ = 
		    std::max(highest_face_checked_, id);

		  if (idx<0 || boxExtreme(face->boundingBox(), dir, ext_pnt))
		    {
		      localExtreme(face, dir, ext_pnt, idx, ext_par);
		    }
		}
	    }
	}
    }
}

//===========================================================================
void
SurfaceModel::localExtreme(ftSurface *face, Point& dir, 
			   Point& ext_pnt, int& ext_id,
			   double ext_par[]) 
//===========================================================================
{
  int id = face->getId();
  double tol2d = 1.0e-4;

  // Convert the surface to a SISLSurf in order to use SISL functions
  // on it. The "false" argument dictates that the SISLSurf will only    
  // copy pointers to arrays, not the arrays themselves.
    shared_ptr<ParamSurface> psurf = face->surface();
    shared_ptr<SplineSurface> surf;
    shared_ptr<BoundedSurface> bsurf;
    const CurveBoundedDomain* bddomain = 0;
    if (psurf->instanceType() == Class_SplineSurface)
	surf = dynamic_pointer_cast<SplineSurface, ParamSurface>(psurf);
    else
    {
	bsurf = dynamic_pointer_cast<BoundedSurface, ParamSurface>(psurf);
	if (bsurf.get())
	  {
	    while (bsurf->underlyingSurface()->instanceType() == 
		   Class_BoundedSurface)
	      bsurf = dynamic_pointer_cast<BoundedSurface, ParamSurface>(bsurf->underlyingSurface());	    
	    RectDomain domain = bsurf->containingDomain();
	    surf = dynamic_pointer_cast<SplineSurface, ParamSurface>(bsurf->underlyingSurface());
	    ASSERT(surf.get() != 0);
	    if (bsurf->isIsoTrimmed(tol2d))
	      {
		RectDomain dom2 = surf->containingDomain();  // To avoid problems due to numerics
		double umin = std::max(domain.umin(), dom2.umin());
		double umax = std::min(domain.umax(), dom2.umax());
		double vmin = std::max(domain.vmin(), dom2.vmin());
		double vmax = std::min(domain.vmax(), dom2.vmax());
    
		vector<shared_ptr<ParamSurface> > sfs = surf->subSurfaces(umin, vmin, umax, vmax);
		surf = dynamic_pointer_cast<SplineSurface, ParamSurface>(sfs[0]);
	      }
	    else
	      bddomain = &(bsurf->parameterDomain());
	  }

	// Make the smallest possible underlying surface
    }
  ASSERT(surf.get() != 0);

  SISLSurf* sislsf = GoSurf2SISL(*(surf.get()), false);
  double epsge = 1.0e-6;
  int numintpt;  // number of single extremal points
  double* pointpar = 0; // array containing the parameter values of single extremal. pt.
  int numintcr; // number of extremal curves
  SISLIntcurve** intcurves = 0;
  int stat = 0;

  s1921(sislsf, dir.begin(), dir.dimension(), 0.0, epsge, 
	&numintpt, &pointpar, &numintcr, &intcurves, &stat);
  MESSAGE_IF(stat!=0, "s1921 returned code: " << stat); 

  // Check if any of the found extremal points are better than the
  // current most extreme point
  vector<Point> curr_pnt;
  vector<double> curr_par;
  int ki;
  for (ki=0; ki<numintpt; ++ki)
    {
      // Evaluate surface
      Point pos = surf->ParamSurface::point(pointpar[2*ki],pointpar[2*ki+1]);
      if (ext_id < 0 || pos*dir > ext_pnt*dir)
	{
	  curr_pnt.push_back(pos);
	  curr_par.insert(curr_par.end(), pointpar+2*ki, pointpar+2*(ki+1));
	}
    }

  for (ki=0; ki<numintcr; ++ki)
    {
      Point pos = surf->ParamSurface::point(intcurves[ki]->epar1[0],
					    intcurves[ki]->epar1[1]);
      if (ext_id < 0 || pos*dir > ext_pnt*dir)
	{
	  curr_pnt.push_back(pos);
	  curr_par.insert(curr_par.end(), intcurves[ki]->epar1, 
			  intcurves[ki]->epar1+2);
	}

      double *pp = intcurves[ki]->epar1 + 2*(intcurves[ki]->ipoint-1);
      pos = surf->ParamSurface::point(pp[0], pp[1]);
      if (ext_id < 0 || pos*dir > ext_pnt*dir)
	{
	  curr_pnt.push_back(pos);
	  curr_par.insert(curr_par.end(), pp, pp+2);
	}
    }

  if (curr_pnt.size() == 0)
    {
      // No better extremal point is found. 
      return;
    }

  if (bddomain != 0)
    {
      // The surface was trimmed
      // Remove extremal points which are outside the domain
      for (ki=0; ki<(int)curr_pnt.size(); ++ki)
	{
	  Vector2D param(curr_par[2*ki], curr_par[2*ki+1]);
	  if (!bddomain->isInDomain(param, epsge))
	    {
	      curr_pnt.erase(curr_pnt.begin()+ki);
	      curr_par.erase(curr_par.begin()+2*ki, curr_par.begin()+2*(ki+1));
	    }
	}
    }

  if (curr_pnt.size() > 0)
    {
      //  Fetch the best extremal point
      ext_id = id;
      ext_pnt = curr_pnt[0];
      ext_par[0] = curr_par[0];
      ext_par[1] = curr_par[1];
      for (ki=1; ki<(int)curr_pnt.size(); ++ki)
	{
	  if (curr_pnt[ki]*dir > ext_pnt*dir)
	    {
	      ext_pnt = curr_pnt[ki];
	      ext_par[0] = curr_par[2*ki];
	      ext_par[1] = curr_par[2*ki+1];
 	    }
	}
    }  

  else if (bddomain != 0)
    {
      // It can be an extremal point inside the face that is better than
      // the previous one.
      // First get the extreme points on the boundary
      vector<CurveLoop> bd_loops = bsurf->allBoundaryLoops();
      for (ki=0; ki<(int)bd_loops.size(); ++ki)
	{
	  vector<shared_ptr<ParamCurve> > curr_loop(bd_loops[ki].begin(),
						    bd_loops[ki].end());
	  shared_ptr<CompositeCurve> comp_cv = 
	    shared_ptr<CompositeCurve>(new CompositeCurve(toptol_.gap,
							  toptol_.neighbour,
							  toptol_.kink,
							  toptol_.bend,
							  curr_loop));
				       
	  int idx;
	  Point bd_ext;
	  double bd_par;
	  comp_cv->extremalPoint(dir, bd_ext, idx, &bd_par);
	  if (ext_id < 0 || bd_ext*dir > ext_pnt*dir)
	    {
	      ext_id = id;
	      ext_pnt = bd_ext;
	      Point param = bsurf->getSurfaceParameter(ki, idx, bd_par);
	      ext_par[0] = param[0];
	      ext_par[1] = param[1];
	    }
	}
	      

      // Triangulate trimmed surface
      int n, m;
      double density = 1.0;
      int min_nmb = 4, max_nmb = 50;
      setResolutionFromDensity(psurf, density, min_nmb, max_nmb, n, m);

      shared_ptr<GeneralMesh> mesh;
      tesselateOneSrf(psurf, mesh, n, m);

      // Get the most extreme triangulation nodes
      double *nodes = mesh->vertexArray();
      // int nmb_nodes = mesh->numVertices();
      int num_triang = mesh->numTriangles();
      double *par_nodes = mesh->paramArray();
      unsigned int *triang_idx = mesh->triangleIndexArray();
      for (ki=0; ki<num_triang; ++ki)
	{
	  // Due to the structure of the tesselation, the points must
	  // be handled more than once
	  for (int kj=0; kj<3; ++kj)
	    {
	      Point node_ext(nodes+triang_idx[ki+kj], 
			     nodes+triang_idx[ki+kj]+3, false);
	      if (ext_id < 0 || node_ext*dir > ext_pnt*dir)
		{
		  ext_id = id;
		  ext_pnt = node_ext;
		  ext_par[0] = par_nodes[2*ki];
		  ext_par[1] = par_nodes[2*ki+1];
		}
	    }
	}

      // Use this value as a start point for an extreme point iteration

    }
  else
    {
      // No better extremal point is found
      return;
    }
}

// //===========================================================================
// void
// SurfaceModel::iterateExtreme(const ftSurface *face, const Point& dir, 
// 			     const Point& ext_pnt, int& ext_it,
// 			     double ext_par[]) const)
// //===========================================================================
// {
// }

} // namespace Go
