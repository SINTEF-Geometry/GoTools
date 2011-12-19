//===========================================================================
//
// File : SplitModelUtils.C
//
// Created: Mars 2010
//
// Author: Vibeke Skytt
//
// Revision: 
//
// Description:
//
//===========================================================================


#include "GoTools/compositemodel/SplitModelUtils.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/ParamSurface.h"

#include <fstream>

using std::vector;
using std::make_pair;


namespace Go
{
//===========================================================================
  void SplitModelUtils::splitInFreeCorners(shared_ptr<SurfaceModel> sfmodel,
					   const Point& pnt, const Point& axis)
//===========================================================================
  {

    // Simplify boundary loops to avoid false free corners
    double max_dist;
    bool simplified;
    simplified = sfmodel->simplifyTrimLoops(max_dist);
    
    // Fetch vertices belonging only to one face
    vector<shared_ptr<Vertex> > vertices;
    sfmodel->getAllVertices(vertices);

    vector<pair<Point,ftSurface*> > vertices2;
    size_t ki;
    for (ki=0; ki<vertices.size(); ++ki)
      {
	vector<pair<ftSurface*, Point> > faces = vertices[ki]->getFaces();
	if (faces.size() == 1)
	  vertices2.push_back(make_pair(vertices[ki]->getVertexPoint(),
					faces[0].first));
      }

    // Split faces corresponding to the identified vertices
    for (ki=0; ki<vertices2.size(); ++ki)
      {
	// Define plane
	Point mid = vertices2[ki].first;
	Point norm = (mid - pnt).cross(axis);
	norm.normalize();

	// Split corresponding surface
	shared_ptr<ParamSurface> srf = vertices2[ki].second->surface();

// 	std::ofstream out_file0("A2_0.g2");
// 	srf->writeStandardHeader(out_file0);
// 	srf->write(out_file0);

	vector<shared_ptr<BoundedSurface> > srf_pieces = 
	  BoundedUtils::splitWithPlane(srf, mid, norm, sfmodel->getTolerances().gap);

	// If a split is performed, replace the modified surface
	if (srf_pieces.size() > 1)
	  {
	    shared_ptr<ftSurface> face = 
	      sfmodel->fetchAsSharedPtr(vertices2[ki].second);
	    sfmodel->removeFace(face);
	    for (size_t kr=0; kr<srf_pieces.size(); ++kr)
	      {
		shared_ptr<ParamSurface> srf2 = srf_pieces[kr];

// 		srf2->writeStandardHeader(out_file0);
// 		srf2->write(out_file0);

		shared_ptr<ftSurface> face2 = 
		  shared_ptr<ftSurface>(new ftSurface(srf2, -1));
		sfmodel->append(face2);
	      }
	  }
      }

  }

//===========================================================================
  void SplitModelUtils::splitInNonCorners(shared_ptr<SurfaceModel> sfmodel,
					  const Point& pnt, const Point& axis)
//===========================================================================
  {
    bool modified = true;
    while (modified)
      {
	modified = false;

	int nmb_sfs = sfmodel->nmbEntities();
	for (int ki=0; ki<nmb_sfs; ++ki)
	  {
	    shared_ptr<ftSurface> face = sfmodel->getFace(ki);
	    shared_ptr<ParamSurface> srf = face->surface();
	    vector<shared_ptr<Vertex> > vertices = 
	      face->getNonCornerVertices(sfmodel->getTolerances().kink);

	    size_t kj, kr;
	    for (kj=0; kj<vertices.size(); ++kj)
	      {
		// Do not split at boundary vertices
		bool boundary = vertices[kj]->isBoundaryVertex();
		if (boundary)
		  continue;

		// Define plane
		Point mid1 = vertices[kj]->getVertexPoint();
		Point norm = (mid1 - pnt).cross(axis);
		norm.normalize();
		for (kr=kj+1; kr<vertices.size(); ++kr)
		  {
		    // Check if the two vertices approximately
		    // defines the same plane
		    Point mid2 = vertices[kr]->getVertexPoint();
		    double dist1 = mid1.dist(mid2);
		    double dist2 = fabs((mid1 - mid2)*norm);
		    // Check if the two vertices lie at the same side of
		    // the axis
		    Point tmp1 = mid1 - ((mid1 - pnt)*axis)*axis;
		    Point tmp2 = mid2 - ((mid2 - pnt)*axis)*axis;
		    double side = tmp1*tmp2;

		    if (side > 0.0 &&
			dist1 > sfmodel->getTolerances().neighbour &&
			dist2 <= sfmodel->getTolerances().neighbour)
		      break;
		  }

// 		std::ofstream out_file0("A2_0.g2");
// 		srf->writeStandardHeader(out_file0);
// 		srf->write(out_file0);
		if (kr < vertices.size())
		  {
		    // Corresponding vertices are found
		    // Fetch parameter values corresponding to the
		    // vertices
		    Point par1 = vertices[kj]->getFacePar(face.get());
		    Point par2 = vertices[kr]->getFacePar(face.get());
	
		    vector<shared_ptr<BoundedSurface> > srf_pieces = 
		      BoundedUtils::splitBetweenParams(srf, par1, par2, 
						       sfmodel->getTolerances().gap);

		    // If a split is performed, replace the modified surface
		    if (srf_pieces.size() > 1)
		      {
			sfmodel->removeFace(face);
			for (kr=0; kr<srf_pieces.size(); ++kr)
			  {
			    shared_ptr<ParamSurface> srf2 = srf_pieces[kr];

// 			    srf2->writeStandardHeader(out_file0);
// 			    srf2->write(out_file0);

			    shared_ptr<ftSurface> face2 = 
			      shared_ptr<ftSurface>(new ftSurface(srf2, -1));
			    sfmodel->append(face2);
			  }

			modified = true;
		      }
		  }
		else
		  {
		    vector<shared_ptr<BoundedSurface> > srf_pieces = 
		      BoundedUtils::splitWithPlane(srf, mid1, norm, 
						   sfmodel->getTolerances().gap);

		    // If a split is performed, replace the modified surface
		    if (srf_pieces.size() > 1)
		      {
			sfmodel->removeFace(face);
			for (kr=0; kr<srf_pieces.size(); ++kr)
			  {
			    shared_ptr<ParamSurface> srf2 = srf_pieces[kr];

// 			    srf2->writeStandardHeader(out_file0);
// 			    srf2->write(out_file0);

			    shared_ptr<ftSurface> face2 = 
			      shared_ptr<ftSurface>(new ftSurface(srf2, -1));
			    sfmodel->append(face2);
			  }

			modified = true;
		      }
		  }

	    if (modified)
	      break;
	      }

	    if (kj < vertices.size())
	      break;
      }
      }
  }

//===========================================================================
  void SplitModelUtils::splitInOuterVertices(shared_ptr<SurfaceModel> sfmodel,
					     shared_ptr<ftSurface> face,
					     const Point& pnt, const Point& axis)
//===========================================================================
  {
    double eps = sfmodel->getTolerances().gap;
    int nmb_sfs = sfmodel->nmbEntities();
    int idx = sfmodel->getIndex(face);
    if (idx < 0 || idx >= nmb_sfs)
      return;   // Face not in face set

    int nmb_loop = face->nmbBoundaryLoops();
    if (nmb_loop <= 1)
      return;   // No inner trimming

    // Get all vertices in outer loop
    vector<shared_ptr<Vertex> > vertices = 
      face->getBoundaryLoop(0)->getVertices();

    // Get all vertices in inner loops
    vector<shared_ptr<Vertex> > inner_vx;
    int ki;
    for (ki=1; ki<nmb_loop; ++ki)
      {
	vector<shared_ptr<Vertex> > tmp_vx =
	  face->getBoundaryLoop(ki)->getVertices();
	inner_vx.insert(inner_vx.end(), tmp_vx.begin(), tmp_vx.end());
      }

    // Get bounded surface corresponding to the face, and 
    // get the underlying surface
    shared_ptr<ParamSurface> surf = face->surface();
    shared_ptr<BoundedSurface> bd_srf;
    if (surf->instanceType() == Class_BoundedSurface)
      bd_srf = dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);
    else if (surf->instanceType() == Class_SplineSurface)
      try {
	shared_ptr<SplineSurface> spline_sf =
	  dynamic_pointer_cast<SplineSurface, ParamSurface>(surf);
	bd_srf = 
	  shared_ptr<BoundedSurface>(BoundedUtils::convertToBoundedSurface(*spline_sf, eps));
      } catch (...) {
	THROW("Something went wrong, returning.");
      }
    else
      THROW("Surface type not supported.");

    shared_ptr<SplineSurface> underlying_sf;
    bool found = bd_srf->hasUnderlyingSpline(underlying_sf);
    if (!found)
      return; 

//     std::ofstream out_file0("A2_trimcrvs.g2");

    // For all outer vertices, define the splitting curve
    size_t kj, kr, kv;
    vector<shared_ptr<CurveOnSurface> > segments;
    for (kj=0; kj<vertices.size(); ++kj)
      {
	// Define plane
	Point mid1 = vertices[kj]->getVertexPoint();
	Point norm = (mid1 - pnt).cross(axis);
	norm.normalize();

	// Check if an inner vertex approximately defines the same plane
	for (kr=0; kr<inner_vx.size(); ++kr)
	  {
	    Point mid2 = inner_vx[kr]->getVertexPoint();
	    double dist1 = mid1.dist(mid2);
	    double dist2 = fabs((mid1 - mid2)*norm);
	    // Check if the two vertices lie at the same side of
	    // the axis
	    Point tmp1 = mid1 - ((mid1 - pnt)*axis)*axis;
	    Point tmp2 = mid2 - ((mid2 - pnt)*axis)*axis;
	    double side = tmp1*tmp2;

	    if (side > 0.0 &&
		dist1 > sfmodel->getTolerances().neighbour &&
		dist2 <= sfmodel->getTolerances().neighbour)
	      break;
	  }

	if (kr < inner_vx.size())
	  {
	    // Corresponding vertices are found
	    // Fetch parameter values corresponding to the
	    // vertices
	    Point par1 = vertices[kj]->getFacePar(face.get());
	    Point par2 = inner_vx[kr]->getFacePar(face.get());

	    // Make parameter curve between parameter values
	    shared_ptr<ParamCurve> pcrv = 
	      shared_ptr<ParamCurve>(new SplineCurve(par1, par2));

	    // Make trimming curve
	    shared_ptr<CurveOnSurface> trimcrv = 
	      shared_ptr<CurveOnSurface>(new CurveOnSurface(underlying_sf, 
							    pcrv, true));
	    bool updated;
	    updated = trimcrv->ensureSpaceCrvExistence(eps);
	    segments.push_back(trimcrv);
	  }
	else
	  {
	    shared_ptr<BoundedSurface> bd_sf2;
	    vector<shared_ptr<CurveOnSurface> > curr_segments =
	      BoundedUtils::getPlaneIntersections(surf, mid1, norm, eps,
						  bd_sf2);

	    // Make sure that the segments lie on the correct side
	    // of the axis
	    for (kv=0; kv<curr_segments.size(); ++kv)
	      {
		Point tmp1 = mid1 - ((mid1 - pnt)*axis)*axis;
		Point tmp2;
		curr_segments[kv]->point(tmp2, curr_segments[kv]->startparam());
		Point tmp3 = tmp2 - ((tmp2 - pnt)*axis)*axis;
		if (tmp1*tmp2 < 0.0)
		  continue;

		segments.push_back(curr_segments[kv]);
	      }
	  }
      }

//     for (kr=0; kr<segments.size(); ++kr)
//       {
// 	segments[kr]->spaceCurve()->writeStandardHeader(out_file0);
// 	segments[kr]->spaceCurve()->write(out_file0);
//       }
	
    // Perform splitting
    vector<shared_ptr<BoundedSurface> > srf_pieces = 
      BoundedUtils::splitWithTrimSegments(bd_srf, segments, eps);

    // If a split is performed, replace the modified surface
    if (srf_pieces.size() > 1)
      {
	sfmodel->removeFace(face);
	for (kr=0; kr<srf_pieces.size(); ++kr)
	  {
	    shared_ptr<ParamSurface> srf2 = srf_pieces[kr];

	    shared_ptr<ftSurface> face2 = 
	      shared_ptr<ftSurface>(new ftSurface(srf2, -1));
	    sfmodel->append(face2);
	  }
      }
  }


}

