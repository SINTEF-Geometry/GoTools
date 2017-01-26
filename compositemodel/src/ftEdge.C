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

#include "GoTools/compositemodel/ftEdge.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/Line.h"
#include "GoTools/geometry/Circle.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/compositemodel/Vertex.h"
#include "GoTools/compositemodel/EdgeVertex.h"


using std::vector;

namespace Go
{



//===========================================================================
ftEdge::ftEdge(ftFaceBase* face,
	       shared_ptr<ParamCurve> cv, 
	       double tmin,
	       double tmax,
	       int entry_id)
//===========================================================================
  : ftEdgeBase(), face_(face), geom_curve_(cv),
      low_param_(tmin), high_param_(tmax),
      entry_id_(entry_id), is_reversed_(false)
{
    ALWAYS_ERROR_IF(tmin > tmax,
		"TMin must be not be greater than TMax");

    Point v1 = cv->point(tmin);
    Point v2 = cv->point(tmax);

    v1_ = shared_ptr<Vertex>(new Vertex(v1, this));
    v2_ = shared_ptr<Vertex>(new Vertex(v2, this));
}

//===========================================================================
ftEdge::ftEdge(shared_ptr<ParamCurve> cv, 
	       double tmin,
	       double tmax,
	       int entry_id)
//===========================================================================
    : ftEdgeBase(), face_(0), geom_curve_(cv),
      low_param_(tmin), high_param_(tmax),
      entry_id_(entry_id), is_reversed_(false)
{
    ALWAYS_ERROR_IF(tmin > tmax,
		"TMin must be not be greater than TMax");

    Point v1 = cv->point(tmin);
    Point v2 = cv->point(tmax);

    v1_ = shared_ptr<Vertex>(new Vertex(v1, this));
    v2_ = shared_ptr<Vertex>(new Vertex(v2, this));
}

//===========================================================================
ftEdge::ftEdge(ftFaceBase* face,
	       shared_ptr<ParamCurve> cv, 
	       shared_ptr<Vertex> v1,
	       shared_ptr<Vertex> v2,
               bool is_reversed,
	       int entry_id)
//===========================================================================
    : ftEdgeBase(), face_(face), geom_curve_(cv),
      entry_id_(entry_id), is_reversed_(is_reversed)
{
    setVertices(v1, v2);
}

//===========================================================================
ftEdge::ftEdge(shared_ptr<ParamCurve> cv, 
	       shared_ptr<Vertex> v1,
	       shared_ptr<Vertex> v2,
               bool is_reversed,
	       int entry_id)
//===========================================================================
    : ftEdgeBase(), face_(0), geom_curve_(cv),
      entry_id_(entry_id), is_reversed_(is_reversed)
{
    setVertices(v1, v2);
}

//===========================================================================
ftEdge::ftEdge(ftFaceBase* face,
	       shared_ptr<ParamCurve> cv, 
	       double tmin,
	       shared_ptr<Vertex> v1,
	       double tmax,
	       shared_ptr<Vertex> v2,
	       int entry_id)
//===========================================================================
    : ftEdgeBase(), face_(face), geom_curve_(cv),
      low_param_(tmin), high_param_(tmax), v1_(v1),
      v2_(v2), entry_id_(entry_id), /*is_turned_(false),*/ is_reversed_(false)
{
    ALWAYS_ERROR_IF(tmin > tmax,
		"TMin must be not be greater than TMax");

    v1_->addEdge(this);
    v2_->addEdge(this);
}

//===========================================================================
ftEdge::~ftEdge()
//===========================================================================
{
  //std::cout << this << std::endl;
}


//===========================================================================
void ftEdge::setVertices(shared_ptr<Vertex> v1,
			 shared_ptr<Vertex> v2)
//===========================================================================
{
    // If the curve is closed, the order of v1 and v2 is significant and
    // determines the direction of traversal of the edge from v1 to v2.
    // tmin or tmax may be shifted by 2pi in order to reflect this. If
    // the curve is open, v1/v2 will be set as v_start/v_end according
    // to the corresponding parameter values.

    Point close1, close2;
    double t1, t2, td1, td2;

#ifndef NDEBUG
    {
        Point start_debug = v1->getVertexPoint();
        Point end_debug = v2->getVertexPoint();
        double dist_debug = start_debug.dist(end_debug);
        double val_debug = dist_debug;
    }
#endif

    geom_curve_->closestPoint(v1->getVertexPoint(), t1, close1, td1);
    geom_curve_->closestPoint(v2->getVertexPoint(), t2, close2, td2);

    double startpar = geom_curve_->startparam();
    double endpar = geom_curve_->endparam();

    // If the curve is closed, i.e. periodic, we do certain things.
    const double pareps = 1.0e-8;
    const double geoeps = 1.0e-6;
    const bool geom_cv_closed = geom_curve_->isClosed();
    const bool edge_cv_closed = (fabs(t1-t2) < pareps);
    if (geom_cv_closed) {
        // For the special case of a circle we must check if the seam should be moved.
        // @@sbr201701 We must also consider cases when a circle sector crosses the seam.
        if (edge_cv_closed) {
            if (geom_curve_->instanceType() == Class_Circle) {
                MESSAGE("We must move the seam of the circle!");
                shared_ptr<Circle> circle_cv = dynamic_pointer_cast<Circle>(geom_curve_);
                // We move the seam by rotating the curve.
#if 0
                Point new_start_pt = circle_cv->ParamCurve::point(circle_cv->startparam());
                GeometryTools::rotatePoint(circle_cv->getNormal(), t1, new_start_pt);
#else
                Point new_start_pt = v1->getVertexPoint();
#endif                
                Point x_axis = new_start_pt - circle_cv->getCentre();
                x_axis.normalize();
                shared_ptr<Circle> rot_circle(new Circle(circle_cv->getRadius(), circle_cv->getCentre(),
                                                         circle_cv->getNormal(), x_axis));
                geom_curve_ = rot_circle;
                // We verify the rotation ...
                geom_curve_->closestPoint(v1->getVertexPoint(), t1, close1, td1);
                geom_curve_->closestPoint(v2->getVertexPoint(), t2, close2, td2);
                // t1 = startpar;
                // t2 = endpar;
                double sum_t_params = t1 + t2; // These should add up to 0.0.
                if (sum_t_params > 0.1) {
                    std::cout << "DEBUG: sum_t_params: " << sum_t_params << std::endl;
                }
            }
        }

	// First snap the endpoints if necessary
	Point startpt, endpt;
	geom_curve_->point(startpt, startpar);
	geom_curve_->point(endpt, endpar);
	if ((startpt - close1).length() < geoeps) {
	    t1 = startpar;
	}
	if ((endpt - close2).length() < geoeps) {
	    t2 = endpar;
	}

	// In the periodic case, we assume that the order of the
	// vertices determines the orientation. Thus - if we cross a
	// seam - the value of t1 found above may be greater than
	// t2. If this is the case, we must subtract 2pi from t1. But:
	// Not if t2=0, which means that t2 actually is 2pi (!).
	if (!is_reversed_ && t1 > t2) {
	    if (t2 == startpar) {
		t2 = endpar;
	    }
	    if (t1 == endpar) {
		t1 = startpar;
	    }
	}
	if (is_reversed_ && t1 < t2) {
	    if (t1 == startpar) {
		t1 = endpar;
	    }
            if (t2 == endpar) {
		t2 = startpar;
	    }
	}
    }

    if (geom_curve_->instanceType() == Class_Line) {
        bool is_bounded = (dynamic_pointer_cast<Line>(geom_curve_))->isBounded();
        if (!is_bounded) { // It should not matter if curve is not bounded, only used for snapping to end parameters.
            ;//MESSAGE("The line is not bounded, did not expect that!");
        }
    }
    
    if (fabs(t2 - t1) < pareps) {
        // @@sbr201701 If the loop is closed this means that the edge is a "self-loop". We must either split
        // the edge into two separate edges or move the seem of the closed curve.
	MESSAGE("t1 ~ t2: Edge is degenerate. edge_cv_closed: " << edge_cv_closed << ". Continuing...");
    }

    // Set the vertices according to parameter order. Snap parameters
    // to endpoints if necessary.
    if (t1 < t2)
    {	
	if (fabs(t1 - startpar) < pareps)
	    t1 = startpar;
	if (fabs(t2 - endpar) < pareps)
	    t2 = endpar;
	low_param_ = t1;
	high_param_ = t2;
	v1_ = v1;
	v2_ = v2;
    }
    else
    {
	if (fabs(t2 - startpar) < pareps)
	    t2 = startpar;
	if (fabs(t1 - endpar) < pareps)
	    t1 = endpar;
	low_param_ = t2;
	high_param_ = t1;
//        MESSAGE("Swapping input vertices to match geometry curve!");
	v1_ = v2;
	v2_ = v1;
        bool debug_is_ok = orientationOK();
        if (!debug_is_ok) {
            std::cout << "Orientation is not ok!" << std::endl;
        }
    }
    v1_->addEdge(this);
    v2_->addEdge(this);
}


//===========================================================================
shared_ptr<Vertex> ftEdge::getVertex(bool at_start)
//===========================================================================
{
    if (at_start) {
	return (isReversed() ? v2_ : v1_);
    }
    else {
	return (isReversed() ? v1_ : v2_);
    }
}


// //===========================================================================
// void ftEdge::turnOrientation()
// //===========================================================================
// {
// //    std::cout << "Turning edge : " << this << " prev: " << prev_;
// //    std::cout << " next: " << next_ << std::endl;
//   is_turned_ = (is_turned_) ? false : true;
// }

// //===========================================================================
// void ftEdge::setOrientation()
// //===========================================================================
// {
//   if (is_turned_)
//     {
//       bool switchparam = true; // @@VSK@@ Since the face is not turned 
//                                 // physically and the surface of a bounded
//                                 // surface does not know that it should be
//                                 // turned, the order of the parameter in 
//                                 // the parameter curve must not be turned.
//                                 // This solution is too complicated to
//                                 // be good, but I do not know any better
//                                 // at the moment (010627).
//       // @@VSK@@ A face is evaluate either after a closest point or after
//       // the evaluation of a parameter curve. In the first case, the
//       // parameter of the closest point is not turned. Thus, it must not
//       // be in the second case either.
// #if ((_MSC_VER > 0) && (_MSC_VER < 1300))
//       ParamCurve *turned_curve = dynamic_cast<ParamCurve*>(geom_curve_->clone());
//       ALWAYS_ERROR_IF(turned_curve==NULL, "Unsuccessfull cast");
// #else
//       ParamCurve *turned_curve = geom_curve_->clone();
// #endif
//       turned_curve->reverseParameterDirection(switchparam);
//       geom_curve_ = shared_ptr<ParamCurve>(turned_curve);
// //          cv1 = geom_curve_->geometryCurve();
// //          cv2 = cv1.get();
// //         std::cout << "turned: " << geom_curve_.get() << ", " << cv2 << std::endl;
//       ftEdgeBase *edum = next_;
//       next_ = prev_;
//       prev_ = edum;

//       shared_ptr<Vertex> tmp_vertex = v1_;
//       v1_ = v2_;
//       v2_ = tmp_vertex;
//     }
//   is_turned_ = false;
// }

// //===========================================================================
// bool ftEdge::isTurned()
// //===========================================================================
// {
//   return is_turned_;
// }

//===========================================================================
void ftEdge::setReversed(bool is_reversed)
//===========================================================================
{
    is_reversed_ = is_reversed;
}

//===========================================================================
bool ftEdge::isReversed()
//===========================================================================
{
    return is_reversed_;
}

//===========================================================================
void ftEdge::reverseGeomCurve()
//===========================================================================
{
    shared_ptr<ParamCurve> newcrv(geom_curve_->subCurve(low_param_,
							high_param_));
    newcrv->reverseParameterDirection();

    geom_curve_ = newcrv;
    v1_.swap(v2_);
    is_reversed_ = !is_reversed_;

    return;
}

//===========================================================================
ftFaceBase* ftEdge::face()
//===========================================================================
{
    return face_;
}


//===========================================================================
BoundingBox ftEdge::boundingBox()
//===========================================================================
{
    return geom_curve_->boundingBox();
}

//===========================================================================
void ftEdge::closestPoint(const Point& pt,
			  double& clo_t,
			  Point& clo_pt,
			  double& clo_dist,
			  double const *seed) const
//===========================================================================
{
    geom_curve_
        ->closestPoint(pt, low_param_, high_param_, clo_t,
		       clo_pt, clo_dist, seed);

    // We may experience that clo_t is outside legal t-values.
    if (clo_t < low_param_) {
	clo_t = low_param_;
	clo_pt = geom_curve_->point(low_param_);
	clo_dist = pt.dist(clo_pt);
    } else if (clo_t > high_param_) {
	clo_t = high_param_;
	clo_pt = geom_curve_->point(high_param_);
	clo_dist = pt.dist(clo_pt);
    }
}

//===========================================================================
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
    ftEdgeBase* ftEdge::split(double t)
#else
    ftEdge* ftEdge::split(double t)
#endif
//===========================================================================
{
#ifdef DEBUG
  if (all_edges_ && !all_edges_->hasEdge(this))
    std::cout << "Split1. Radial edge missing" << std::endl;
#endif

  // If t is close, but not equal, to existing knot, we make it equal.
    double knot_diff_tol = 1e-05;
    shared_ptr<SplineCurve> spline_cv =
      dynamic_pointer_cast<SplineCurve, ParamCurve>(geom_curve_);
    if (spline_cv.get() != 0)
      spline_cv->basis().knotIntervalFuzzy(t, knot_diff_tol);

//     ALWAYS_ERROR_IF(twin(),
// 		"Cannot split edge, already has twin!");
//     ALWAYS_ERROR_IF(!prev() || !next(), 
// 		"Cannot split edge, not fully connected");
    ALWAYS_ERROR_IF(t <= low_param_ || t >= high_param_,
		"Split parameter not in interior of edge range");

    // Save initial vertices
    shared_ptr<Vertex> v1 = v1_;
    shared_ptr<Vertex> v2 = v2_;

    Point split_pt = geom_curve_->point(t);
    shared_ptr<Vertex> split_vx = shared_ptr<Vertex>(new Vertex(split_pt));
    split_vx->addEdge(this);

    //shared_ptr<Vertex> tmp_vx = is_turned_ ? v1_ : v2_;
    shared_ptr<Vertex> tmp_vx = v2_;
    if (tmp_vx.get() != v1_.get())
      tmp_vx->removeEdge(this);  // Don't remove edge for a closed one-edge loop

    ftEdge* newedge;
    newedge = new ftEdge(face_, geom_curve_, t, split_vx, high_param_, tmp_vx);
    high_param_ = t;
//     if (is_turned_)
// 	v1_ = split_vx;
//     else
    v2_ = split_vx;
    newedge->connectAfter(this);
//     if (is_turned_)
//       newedge->turnOrientation();

    // Split radial edge and all associated half edges
    if (twin_)
      {
	int status = 1;
	ftEdge* e2 = twin_->geomEdge();
	if (e2)
	  {
	    e2->ftEdgeBase::disconnectTwin();
	    shared_ptr<ftEdge> e3 = e2->splitAtVertex(split_vx);
	    shared_ptr<Vertex> v3, v4;
	    e2->getVertices(v3, v4);
	    if (v3.get() == v1_.get() || v4.get() == v1_.get())
	      {
		e2->ftEdgeBase::connectTwin(this, status);
		e3->ftEdgeBase::connectTwin(newedge, status);
	      }
	    else
	      {
		e3->ftEdgeBase::connectTwin(this, status);
		e2->ftEdgeBase::connectTwin(newedge, status);
	      }
	  }
      }

    if (all_edges_.get())
      {
	all_edges_->addEdge(newedge);
	newedge->setEdgeVertex(all_edges_);
	all_edges_->splitAtVertex(v1, v2, split_vx);
      }
		

#ifdef DEBUG
  if (all_edges_ && !all_edges_->hasEdge(this))
    std::cout << "Split2. Radial edge missing" << std::endl;
  if (newedge->all_edges_ && !newedge->all_edges_->hasEdge(newedge))
    std::cout << "Split3. Radial edge missing" << std::endl;
#endif

    return newedge;
}

//===========================================================================
shared_ptr<ftEdge> ftEdge::split2(double t)
//===========================================================================
{
  ftEdge *e1 = split(t);
  shared_ptr<ftEdge> e2 =  shared_ptr<ftEdge>(e1);
  e2->face()->updateBoundaryLoops(e2);
  return e2;
}

//===========================================================================
shared_ptr<ftEdge> ftEdge::splitAtVertex(shared_ptr<Vertex> vx)
//===========================================================================
  {
    shared_ptr<ftEdge> dummy;
    if (vx.get() == v1_.get() || vx.get() == v2_.get())
      return dummy;   // Cannot split

    // Find split parameter
    double par, dist;
    Point pt;
    closestPoint(vx->getVertexPoint(), par, pt, dist);

    ALWAYS_ERROR_IF(par <= low_param_ || par >= high_param_,
		"Split parameter not in interior of edge range");

    // If the split parameter is close, but not equal, to existing knot, 
    // we make it equal.
    double knot_diff_tol = 1e-05;
    shared_ptr<SplineCurve> spline_cv =
      dynamic_pointer_cast<SplineCurve, ParamCurve>(geom_curve_);
    if (spline_cv.get() != 0)
      spline_cv->basis().knotIntervalFuzzy(par, knot_diff_tol);

    vx->addEdge(this);

    shared_ptr<Vertex> tmp_vx = v2_;
    if (tmp_vx.get() != v1_.get())
      tmp_vx->removeEdge(this);  // Don't remove edge for a closed one-edge loop

    ftEdge* e1 = new ftEdge(face_, geom_curve_, par, vx, high_param_, tmp_vx);
    high_param_ = par;

    v2_ = vx;
    e1->connectAfter(this);

    // This split function is expected to be called in connection with 
    // splitting of radial edges and does not concernt about all_edges_

    shared_ptr<ftEdge> e2 =  shared_ptr<ftEdge>(e1);
    face_->updateBoundaryLoops(e2);
    return e2;
  }

//===========================================================================
void ftEdge::connectAfter(ftEdgeBase* edge)
//===========================================================================
{
    ftEdgeBase::connectAfter(edge);

    ftEdge *tmp_edge = edge->geomEdge();
    if (tmp_edge)
    {
	// We must use getVertex() - not v1_ or v2_ - in order to get
	// correct vertex. The difference is the orientation.
	shared_ptr<Vertex> vx = getVertex(true);
	shared_ptr<Vertex> tmp_vx = tmp_edge->getVertex(false);
	if (vx.get() != tmp_vx.get())
	    joinVertex(vx, tmp_vx);
    }
}

//===========================================================================
void ftEdge::closeLoop(ftEdgeBase* last)
//===========================================================================
{
    ftEdgeBase::closeLoop(last);

    ftEdge *tmp_edge = last->geomEdge();
    if (tmp_edge)
    {
	// We must use getVertex() - not v1_ or v2_ - in order to get
	// correct vertex. The difference is the orientation.
	shared_ptr<Vertex> vx = getVertex(true);
	shared_ptr<Vertex> tmp_vx = tmp_edge->getVertex(false);
	if (vx.get() != tmp_vx.get())
	    joinVertex(vx, tmp_vx);
    }
}

//===========================================================================
void ftEdge::disconnectThis()
//===========================================================================
{
    ftEdgeBase::disconnectThis();

    v1_->removeEdge(this);
    v2_->removeEdge(this);

    if (all_edges_)
      all_edges_->removeEdge(this);
    all_edges_.reset();

    v1_ = shared_ptr<Vertex>(new Vertex(this, true));
    v2_ = shared_ptr<Vertex>(new Vertex(this, false));
}

//===========================================================================
void ftEdge::connectTwin(ftEdgeBase* newtwin, int& status)
//===========================================================================
{
  bool twin1 = (twin_ != NULL);
  bool twin2 = (newtwin->twin() != NULL);
  bool reorganize = true;

#ifdef DEBUG
  // TEST
  // if (all_edges_ && !all_edges_->checkTwins())
  //   std::cout << "Connect1. Radial edge inconsistency" << std::endl;
  MESSAGE("EdgeVertex::checkTwins() removed!");
  if (all_edges_ && !all_edges_->hasEdge(this))
    std::cout << "Connect1. Radial edge missing" << std::endl;
#endif     

  if (twin1 && twin_ == newtwin && twin2 && this == newtwin->twin())
    {
      // Already connected
      reorganize = false;
      if (newtwin->geomEdge())
	newtwin->geomEdge()->addEdgeMultiplicityInstance(this);
    }
  else if (twin1 || twin2)
    {
      // One of the edges has already a twin
      // Use an EdgeVertex instance to store the information
      if (newtwin->geomEdge())
	newtwin->geomEdge()->addEdgeMultiplicityInstance(this);
    }
  else
    {
      ftEdgeBase::connectTwin(newtwin, status);
      if (all_edges_.get())
	all_edges_->addEdge(newtwin->geomEdge());
      if (newtwin->geomEdge()->hasEdgeMultiplicity())
	newtwin->geomEdge()->getEdgeMultiplicityInstance()->addEdge(this);
      if (all_edges_.get())
	newtwin->geomEdge()->joinEdgeVertex(all_edges_);
      if (newtwin->geomEdge()->hasEdgeMultiplicity())
	joinEdgeVertex(newtwin->geomEdge()->getEdgeMultiplicityInstance());
    }

  joinVertices(newtwin);
  if (reorganize && (twin1 || twin2))
    {
      // The pairing of half edges may not be optimal. Reorganize
      if (all_edges_.get())
	all_edges_->reOrganize();
    }
#ifdef DEBUG
  // TEST
  // if (all_edges_ && !all_edges_->checkTwins())
  //   std::cout << "Connect2. Radial edge inconsistency" << std::endl;
  MESSAGE("EdgeVertex::checkTwins() removed!");
  if (all_edges_ && !all_edges_->hasEdge(this))
    std::cout << "Connect2. Radial edge missing" << std::endl;
#endif
}

//===========================================================================
void ftEdge::joinVertices(ftEdgeBase* newtwin)
//===========================================================================
{
    ftEdge *tmp_twin = newtwin->geomEdge();

    if (tmp_twin)
    {
	if (v1_->getDist(tmp_twin->getVertex(true)) < 
	    v1_->getDist(tmp_twin->getVertex(false)))
	{
	  if (v1_.get() == tmp_twin->v1_.get() &&
	      v2_.get() == tmp_twin->v2_.get())
	    {
                v1_->reOrganize();
                v2_->reOrganize();
	    }
	  else
	    {
	      shared_ptr<Vertex> tmp_vx1 = tmp_twin->getVertex(true);
	      joinVertex(v1_, tmp_vx1);
	      // 	    tmp_twin->getVertex(true)->joinVertex(v1_);
	      // 	    v1_ = tmp_twin->getVertex(true);
	      shared_ptr<Vertex> tmp_vx2 = tmp_twin->getVertex(false);
	      joinVertex(v2_, tmp_vx2);
	      // 	    tmp_twin->getVertex(false)->joinVertex(v2_);
	      // 	    v2_ = tmp_twin->getVertex(false);
	    }
	}
	else
	{
	  if (v1_.get() == tmp_twin->v2_.get() &&
	      v2_.get() == tmp_twin->v1_.get())
	    {
                v1_->reOrganize();
                v2_->reOrganize();
	    }
	  else
	    {
	      shared_ptr<Vertex> tmp_vx1 = tmp_twin->getVertex(false);
	      joinVertex(v1_, tmp_vx1);
	      // 	    tmp_twin->getVertex(false)->joinVertex(v1_);
	      // 	    v1_ = tmp_twin->getVertex(false);
	      shared_ptr<Vertex> tmp_vx2 = tmp_twin->getVertex(true);
	      joinVertex(v2_, tmp_vx2);
	      // 	    tmp_twin->getVertex(true)->joinVertex(v2_);
	      // 	    v2_ = tmp_twin->getVertex(true);
	    }
	}
    }
}

//===========================================================================
void ftEdge::joinVertex(shared_ptr<Vertex> this_vertex,
			shared_ptr<Vertex> other_vertex) 
//===========================================================================
{
    other_vertex->joinVertex(this_vertex);
    vector<ftEdge*> edges = this_vertex->allEdges();
    for (size_t ki=0; ki<edges.size(); ++ki)
	edges[ki]->replaceVertex(this_vertex, other_vertex);
}

//===========================================================================
void ftEdge::replaceVertex(shared_ptr<Vertex>& this_vertex, 
			   shared_ptr<Vertex>& other_vertex)
//===========================================================================
{
    if (v1_.get() == this_vertex.get())
	v1_ = other_vertex;
    else if (v2_.get() == this_vertex.get())
	v2_ = other_vertex;
}

//===========================================================================
void ftEdge::disconnectTwin()
//===========================================================================
{
  if (!twin_)
    return;   // Nothing to do

#ifdef DEBUG
  // TEST
  // if (all_edges_ && !all_edges_->checkTwins())
  //   std::cout << "Disconnect1. Radial edge inconsistency" << std::endl;
  MESSAGE("EdgeVertex::checkTwins() removed!");
  if (all_edges_ && !all_edges_->hasEdge(this))
    std::cout << "Disonnect1. Radial edge missing" << std::endl;
#endif

  ftEdge *prev = dynamic_cast<ftEdge*>(prev_);
  ftEdge *next = dynamic_cast<ftEdge*>(next_);
  shared_ptr<Vertex> prev_start = prev->getVertex(true);
  shared_ptr<Vertex> next_start = next->getVertex(true);
  shared_ptr<Vertex> prev_end = prev->getVertex(false);
  shared_ptr<Vertex> next_end = next->getVertex(false);
  bool prev_v1 = (v1_->hasEdge(prev));
  bool next_v2 = (v2_->hasEdge(next));
  bool at_start1 = (prev_v1 && prev_start == v1_) || 
    ((!prev_v1) && prev_start == v2_) ;
  bool at_start2 = (next_v2 && next_start == v2_) ||
    ((!next_v2) && next_start == v1_);

  if (!twin_)
    return;  // Should not happen

  ftEdge *twin = twin_->geomEdge();

    ftEdgeBase::disconnectTwin();
//     v1_->disconnectTwin(this);
//     v2_->disconnectTwin(this);
    v1_->removeEdge(this);
    v2_->removeEdge(this);

    v1_ = shared_ptr<Vertex>(new Vertex(this, true));
    v2_ = shared_ptr<Vertex>(new Vertex(this, false));

    shared_ptr<Vertex> tmp_vx1 = prev_v1 ? prev->getVertex(at_start1) : 
      next->getVertex(!at_start2);
    joinVertex(v1_, tmp_vx1);
    shared_ptr<Vertex> tmp_vx2 = next_v2 ? next->getVertex(at_start2) : 
      prev->getVertex(!at_start1);
    joinVertex(v2_, tmp_vx2);

    // Remove from radial edge
    if (all_edges_.get())
      {
	all_edges_->disconnectTwin(this, twin);
      }
    
#ifdef DEBUG
  // TEST
  // if (all_edges_ && !all_edges_->checkTwins())
  //   std::cout << "Disconnect2. Radial edge inconsistency" << std::endl;
  MESSAGE("EdgeVertex::checkTwins() removed!");
  if (all_edges_ && !all_edges_->hasEdge(this))
    std::cout << "Disonnect2. Radial edge missing" << std::endl;
#endif
}

//===========================================================================
Point ftEdge::point(double t) const
//===========================================================================
{
    return geom_curve_->point(t);
}

//===========================================================================
Point ftEdge::tangent(double t) const
//===========================================================================
{
  std::vector<Point> point = geom_curve_->point(t, 1);
  if (is_reversed_)
      point[1] *= -1.0;

  return point[1];
}

//===========================================================================
void ftEdge::point(double t, int der, std::vector<Point>& derivs) const
//===========================================================================
{
    derivs = geom_curve_->point(t, der);
    if (is_reversed_)
    {
	for (int ki=1; ki<=der; ++ki)
	    derivs[ki] *= -1.0;
    }
}

//===========================================================================
Point ftEdge::normal(double t) const
//===========================================================================
{
  Point pt2 = faceParameter(t);

  return face_->normal(pt2[0], pt2[1]);
}

//===========================================================================
Point ftEdge::normal(double t, Point& face_par_pt, double* face_seed) const
//===========================================================================
{
  Point pt2 = faceParameter(t, face_seed);
  face_par_pt = pt2;

  return face_->normal(pt2[0], pt2[1]);
}

//===========================================================================
ftEdge* ftEdge::geomEdge()
//===========================================================================
{
  return this;
}

//===========================================================================
double ftEdge::estimatedCurveLength()
//===========================================================================
{
    // Computes four (hardcoded) points and the cordlengths.
    // This gives a low estimate.
//    std::cout << "Edge turned? " << isTurned() << ". Face turned? " << Face()->getOrientation() << std::endl;
//    std::cout << "Instance type: " << geom_curve_->instanceType() << std::endl;
    const int numpts = 4;
    Point pprev = point(tMin());
    Point pnext;
    double length = 0;
    for (int i = 1; i < numpts; ++i) {
	double fac = double(i)/double(numpts-1);
	pnext = point((1.0-fac)*tMin() + fac*tMax());
	length += pnext.dist(pprev);
	pprev = pnext;
    }
    return length;
}

//===========================================================================
  double ftEdge::estimatedCurveLength(double min_par, double max_par)
//===========================================================================
{
    // Computes four (hardcoded) points and the cordlengths.
    // This gives a low estimate.
//    std::cout << "Edge turned? " << isTurned() << ". Face turned? " << Face()->getOrientation() << std::endl;
//    std::cout << "Instance type: " << geom_curve_->instanceType() << std::endl;
    const int numpts = 4;
    Point pprev = point(min_par);
    Point pnext;
    double length = 0;
    for (int i = 1; i < numpts; ++i) {
	double fac = double(i)/double(numpts-1);
	pnext = point((1.0-fac)*min_par + fac*max_par);
	length += pnext.dist(pprev);
	pprev = pnext;
    }
    return length;
}

//===========================================================================
Point ftEdge::faceParameter(double t, double* seed) const
//===========================================================================
{
  // Check if the underlying curve is a curve on surface. In that case
  // effective search methods may exist
  shared_ptr<CurveOnSurface> sf_cv = 
    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(geom_curve_);
  if (sf_cv.get())
    return sf_cv->faceParameter(t);
  else
    {
      // Find the point on the curve
      Point pt = point(t);
      Point clo_pt;
      double clo_u, clo_v, clo_dist;
      // Find the closest point on the surface
      // std::cout << srf.get() << std::endl;
      face_->surface()->closestBoundaryPoint(pt, clo_u, clo_v, clo_pt, clo_dist,
					     1e-10, NULL, seed);

      return Point(clo_u, clo_v);
    }
}

//===========================================================================
bool ftEdge::orientationOK() const
//===========================================================================
{
    // @jbt: This function is probably redundant, since evaluating
    // point() is by definition equivalent to evaluating
    // geom_curve_->point(). It is probably a remainder from before
    // is_reversed_ was implemented (?). We simply return 'true' and
    // hope the best...

    // The low_param_ is always related to v1_, even if reversed_ is true.
    Point low = geom_curve_->point(low_param_);
    Point high = geom_curve_->point(high_param_);
    double d1 = low.dist(v1_->getVertexPoint()) + high.dist(v2_->getVertexPoint());
    double d2 = low.dist(v2_->getVertexPoint()) + high.dist(v1_->getVertexPoint());
    bool isOK = (d1 <= d2); // For a closed curve with snapped end params the test is inconclusive.
    if (!isOK) {
        MESSAGE("orientationOK(): Not OK!");
    }

    return isOK;


//     // Evaluate edge
//     Point pe1, pe2;
//     pe1 = point(tMin());
//     pe2 = point(tMax());

//     // Evaluate curve
//     Point pc1, pc2;
//     pc1 = geom_curve_->point(tMin());
//     pc2 = geom_curve_->point(tMax());

//     double d1 = pe1.dist(pc1) + pe2.dist(pc2);
//     double d2 = pe1.dist(pc2) + pe2.dist(pc1);
//     bool isOK = ((d1 < d2 && !is_reversed_) || d1 >= d2 && is_reversed_);

//     return isOK;
}


//===========================================================================
void ftEdge::setGeomCurve(shared_ptr<ParamCurve> geom_curve)
//===========================================================================
{
    geom_curve_ = geom_curve;
}


//===========================================================================
void ftEdge::updateGeomCurve(double tol)
//===========================================================================
{

    Point par1 = faceParameter(tMin());
    Point par2 = faceParameter(tMax());

    face()->clearInitialEdges();
    vector<shared_ptr<ftEdgeBase> > tmp_edges = face()->createInitialEdges(tol);
    for (size_t ki=0; ki<tmp_edges.size(); ki++)
    {
	// Identify the current edge
	Point par3 = tmp_edges[ki]->geomEdge()->faceParameter(tMin());
	Point par4 = tmp_edges[ki]->geomEdge()->faceParameter(tMax());

	if (par1.dist(par3) < tol && par2.dist(par4) < tol)
	{
	    geom_curve_ = tmp_edges[ki]->geomEdge()->geomCurve();
	    break;
	}
	else if (par1.dist(par4) < tol && par2.dist(par3) < tol)
	{
	    geom_curve_ = tmp_edges[ki]->geomEdge()->geomCurve();
	    //turnOrientation();  !!!!!! VSK NB!!!!!!!!
	    break;
	}
    }
}

//===========================================================================
bool ftEdge::updateEdgeInfo(double tol)
//===========================================================================
{
  // Fetch geometry curve as curve on surface
  shared_ptr<CurveOnSurface> sfcv = 
    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(geom_curve_);
  if (sfcv.get())
    return sfcv->updateCurves(tol);
  else
    return false;
}
//===========================================================================
int ftEdge::getCurveIndex() const
//===========================================================================
{
  // Fetch geometry curve as curve on surface
  shared_ptr<CurveOnSurface> sfcv = 
    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(geom_curve_);
  if (sfcv.get())
    {
      double tol = 1.0e-8;
      bool same_orient;
      return sfcv->whichBoundary(tol, same_orient);
    }
  else
    return -1;
}

//===========================================================================
double ftEdge::parAtVertex(const Vertex* vx) const
//===========================================================================
{
    if (v1_.get() == vx)
	return low_param_;
    else if (v2_.get() == vx)
	return high_param_;
    else
	return -MAXDOUBLE;  //Nonsense

//     if (v1_.get() == vx && !is_turned_)
// 	return low_param_;
//     else if (v1_.get() == vx)
// 	return high_param_;
//     else if (v2_.get() == vx && !is_turned_)
// 	return high_param_;
//     else if (v2_.get() == vx)
// 	return low_param_;
//     else
// 	return -MAXDOUBLE;  //Nonsense
}

//===========================================================================
  void ftEdge::addEdgeMultiplicityInstance(ftEdge *other)
//===========================================================================
{
  // Check input
  if (all_edges_ && other->all_edges_ && 
      all_edges_.get() == other->all_edges_.get())
    return;  // Nothing to do

  // Collect edges
  vector<ftEdge*> edges;
  edges.push_back(this);
  if (twin_ && twin_ != other)
    edges.push_back(twin_->geomEdge());
  edges.push_back(other);
  if (other->twin_ && other->twin_ != this)
    edges.push_back(other->twin_->geomEdge());

  // Check if an EdgeVertex instance exists
  size_t kj;
  shared_ptr<EdgeVertex> all_edges;
  for (kj=0; kj<edges.size(); ++kj)
    {
      if (edges[kj]->hasEdgeMultiplicity())
	{
	  all_edges = edges[kj]->getEdgeMultiplicityInstance();
#ifdef DEBUG
	  // if (!all_edges->checkTwins())
	  //   std::cout << "Add edge multiplicity 1(" << kj <<") inconsistency" << std::endl;
  MESSAGE("EdgeVertex::checkTwins() removed!");
#endif
	  for (size_t ki=0; ki<edges.size(); ++ki)
	    all_edges->addEdge(edges[ki]);

	  // If two instances exist, let all_edges point to the one
	  // belonging to the other edge
	  if (all_edges.get() != all_edges_.get())
	    break;
	}
    }
	  
  if (!all_edges.get())
    all_edges_ = shared_ptr<EdgeVertex>(new EdgeVertex(edges));
  else
    {
 #ifdef DEBUG
     // if (!all_edges->checkTwins())
     // 	std::cout << "Add edge multiplicity 2 inconsistency" << std::endl;
  MESSAGE("EdgeVertex::checkTwins() removed!");
#endif
      if (all_edges_.get())
	{
	  // Make sure that all edges will have pointers to the
	  // merged edge vertex instance and that the edge vertex 
	  // instance has pointers to all edges
	  vector<ftEdge*> this_edges = all_edges_->allEdges();
	  for (size_t kr=0; kr<this_edges.size(); ++kr)
	    {
	      this_edges[kr]->setEdgeVertex(all_edges);
	      all_edges->addEdge(this_edges[kr]);
	    }
	  
#ifdef DEBUG
	  // if (!all_edges->checkTwins())
	  //   std::cout << "Add edge multiplicity 3 inconsistency" << std::endl;
  MESSAGE("EdgeVertex::checkTwins() removed!");
 	  for (size_t kr=0; kr<this_edges.size(); ++kr)
	    if (!all_edges->hasEdge(this_edges[kr]))
	      std::cout << "Add edge multiplicity. Missing edge " << kr << std::endl;
#endif

	  all_edges_ = all_edges;
	}
      else
	all_edges_ = all_edges;
    }

  for (kj=0; kj<edges.size(); ++kj)
    {
      if (edges[kj] == this)
	continue;

      edges[kj]->joinEdgeVertex(all_edges_);
    }
#ifdef DEBUG
  // TEST
  // if (all_edges_ && !all_edges_->checkTwins())
  //   std::cout << "Add radial edge. Radial edge inconsistency" << std::endl;
  MESSAGE("EdgeVertex::checkTwins() removed!");
  if (all_edges_ && !all_edges_->hasEdge(this))
    std::cout << "Add radial edge. Radial edge missing" << std::endl;
  for (size_t kr=0; kr<edges.size(); ++kr)
    if (edges[kr]->all_edges_ && !edges[kr]->all_edges_->hasEdge(edges[kr]))
      std::cout << "Add radial edge. Missing edge nr: " << kr << std::endl;
#endif
}

//===========================================================================
  void ftEdge::joinEdgeVertex(shared_ptr<EdgeVertex> radial_edge)
//===========================================================================
  {
#ifdef DEBUG
  // TEST
  // if (radial_edge && !radial_edge->checkTwins())
  //   std::cout << "Join radial edge 1. Radial edge inconsistency" << std::endl;
  MESSAGE("EdgeVertex::checkTwins() removed!");
  if (all_edges_ && !all_edges_->hasEdge(this))
    std::cout << "Join radial edge. Radial edge missing" << std::endl;
#endif

    if (!all_edges_.get())
      all_edges_ = radial_edge;
    else
      {
#ifdef DEBUG
	// TEST
	// if (all_edges_ && !all_edges_->checkTwins())
	//   std::cout << "Join radial edge 2. Radial edge inconsistency" << std::endl;
  MESSAGE("EdgeVertex::checkTwins() removed!");
#endif
	radial_edge->addEdgeVertex(all_edges_.get());
	all_edges_ = radial_edge;

	// Ensure that all involved edges point to the same
	// radial edge instance
	vector<ftEdge*> all_edges = all_edges_->allEdges();

	for (size_t ki=0; ki<all_edges.size(); ++ki)
	  all_edges[ki]->setEdgeVertex(radial_edge);
      }
#ifdef DEBUG
  // TEST
  // if (all_edges_ && !all_edges_->checkTwins())
  //   std::cout << "Join radial edge 3. Radial edge inconsistency" << std::endl;
  MESSAGE("EdgeVertex::checkTwins() removed!");
  if (all_edges_ && !all_edges_->hasEdge(this))
    std::cout << "Join radial edge3. Radial edge missing" << std::endl;
#endif
  }


//===========================================================================
void ftEdge::removeEdgeVertex()
//===========================================================================
  {
    shared_ptr<EdgeVertex> dummy;
    all_edges_ = dummy;
  }

//---------------------------------------------------------------------------
  vector<ftSurface*> ftEdge::getAdjacentFaces() const
//---------------------------------------------------------------------------
  {
    vector<ftSurface*> faces;
    if (face_)
      {
	ftSurface *curr = face_->asFtSurface();
	if (curr)
	  faces.push_back(curr);
	if (twin_)
	  {
	    curr = 0;
	    if (twin_->geomEdge() && twin_->geomEdge()->face_ &&
		twin_->geomEdge()->face_ != curr)
	      curr = twin_->geomEdge()->face_->asFtSurface();
	    if (curr)
	      faces.push_back(curr);
	  }
      }

    return faces;
  }

//---------------------------------------------------------------------------
  vector<ftSurface*> ftEdge::getAllAdjacentFaces() const
//---------------------------------------------------------------------------
  {
    if (all_edges_.get())
      return all_edges_->getAdjacentFaces();
    else
      return getAdjacentFaces();
  }

//---------------------------------------------------------------------------
bool ftEdge::checkEdgeTopology()
//---------------------------------------------------------------------------
  {
    bool isOK = true;
    Point pos1 = geom_curve_->point(low_param_);
    Point pos2 = geom_curve_->point(high_param_);
    double dist1 = pos1.dist(v1_->getVertexPoint());
    double dist2 = pos2.dist(v1_->getVertexPoint());
    double dist3 = pos1.dist(v2_->getVertexPoint());
    double dist4 = pos2.dist(v2_->getVertexPoint());
    if (dist1 > 0.01 && dist2 > 0.01)
      {
      std::cout << "Vertex - point inconsistence, edge = " << this;
      std::cout << ", vertex = " << v1_ << std::endl;
      isOK = false;
      }
     if (dist3 > 0.01 && dist4 > 0.01)
      {
      std::cout << "Vertex - point inconsistence, edge = " << this;
      std::cout << ", vertex = " << v2_ << std::endl;
      isOK = false;
      }

     int highval = 1000;
     if (entry_id_ >= highval || entry_id_ < -1)
       {
	 std::cout << "Edge entry: " << this << "(" << entry_id_ << ")" << std::endl;
	 isOK = false;
       }
     if (face_->getId() >= highval || face_->getId() < -1)
       {
	 std::cout << "Face entry: " << this << ", " << face_;
	 std::cout << "(" << face_->getId() << ")" << std::endl;
	 isOK = false;
       }
      return isOK;
  }

} // namespace Go
