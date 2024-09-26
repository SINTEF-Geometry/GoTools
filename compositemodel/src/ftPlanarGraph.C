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

#include <algorithm>
#include "GoTools/compositemodel/ftPlanarGraph.h"
#include <fstream>
#include "GoTools/compositemodel/ftSurfaceSetPoint.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/geometry/LineCloud.h"

using std::make_pair;

namespace Go
{

const double knot_tol = 1e-18;

// Return the intersection point between the two line segments. If lines are parallell the returned
// point has dimension 0.
Vector2D intersect2DLines(Vector2D from1, Vector2D dir1, Vector2D from2, Vector2D dir2);
    
// Ordering on ftSearchNode, given by higher y-value.
// If equal y-value, sort by higher x-value.
bool nodeSortBool(const ftSearchNode& node1, const ftSearchNode& node2)
{
    if (fabs(node1.node()[1] - node2.node()[1]) < knot_tol)
	if (fabs(node1.node()[0] - node2.node()[0]) < knot_tol)
	    return true;
	else
	    return (node1.node()[0] < node2.node()[0]);
    else
	return (node1.node()[1] < node2.node()[1]);
}

// Used for locating included node.
class nodeRuntimeBool
{
public :
    nodeRuntimeBool(Vector2D node) {
	node_ = node;
    }

    bool operator()(ftSearchNode& point) {
	return ((point.node() - node_).length() == 0);
    }
private:
//     PointIter node_;
    Vector2D node_;
};

// Ordering of edges (with common start point), given by smaller angle with x-axis.
bool edgeSortBool(const ftGraphEdge& edge1, const ftGraphEdge& edge2) {

    Vector2D first_leg = edge1.endPoint() - edge1.startPoint();
    Vector2D second_leg = edge2.endPoint() - edge2.startPoint();

    double cross = first_leg[0]*second_leg[1] - second_leg[0]*first_leg[1];
    if (cross == 0)
	return (first_leg[0] < second_leg[0]); // Vectors are parallell.
    else
	return (cross < 0);
}


// The class decides on which side of a (directed) edge a point is.
// Used when searching ordered_segments_[i].
// False if point is to the right. Use upper_bound (not lower_bound)!
class edgeRuntimeBool
{
public:
    edgeRuntimeBool(const Vector2D& point)
	: point_(point)
    {};
    //	bool operator()(ftGraphEdge& edge1, ftGraphEdge& edge2) {
    bool operator()(const ftGraphEdge& edge1) {
	if (point_.dist(edge1.endPoint()) == 0)
	    return 1; // We then consider the edge as to the right of point.
	else {
	    Vector2D first_leg = point_ - edge1.startPoint();
	    Vector2D second_leg = edge1.endPoint() - edge1.startPoint();
	    double cross =
		first_leg[0]*second_leg[1] - second_leg[0]*first_leg[1];
	    return (cross < 0);
	}
    }
    // We use edge2 as dummy var (whith lower_bound). Also used for sorting.
    bool operator()(const ftGraphEdge& edge1, const ftGraphEdge& edge2) {
	    return (this->operator()(edge2));
    }
private:
    Vector2D point_;
};



//===========================================================================
ftGraphEdge::ftGraphEdge()
//===========================================================================
{
}

//===========================================================================
ftGraphEdge::ftGraphEdge(ftSurfaceSetPoint*& lower, ftSurfaceSetPoint*& upper)
//===========================================================================
    :lower_(lower->ftSamplePoint::getPar()), upper_(upper->ftSamplePoint::getPar())
{
    // We check whether faces of nodes correspond. Push_back existing faces of edge.
    for (int i = 0; i < lower->nmbFaces(); ++i)
	for (int j = 0; j < upper->nmbFaces(); ++j)
	    if (lower->face(i) == upper->face(j)) {
		faces_.push_back((lower->face(i)));
		face_params_.push_back(make_pair(lower->parValue(i),
						 upper->parValue(j)));
	    }

    ALWAYS_ERROR_IF(faces_.size() == 0,
		"The two nodes did not correspond to an edge.");
}

//===========================================================================
ftGraphEdge::~ftGraphEdge()
//===========================================================================
{
}


//===========================================================================
double ftGraphEdge::startparam() const
//===========================================================================
{
    return lower_[1];
}

//===========================================================================
double ftGraphEdge::endparam() const
//===========================================================================
{
    return upper_[1];
}

//===========================================================================
Vector2D ftGraphEdge::startPoint() const
//===========================================================================
{
    return lower_;
}

//===========================================================================
Vector2D ftGraphEdge::endPoint() const
//===========================================================================
{
    return upper_;
}

//===========================================================================
Vector2D ftGraphEdge::point(double v_par)
//===========================================================================
{
    if ((v_par < lower_[1]) || (v_par > upper_[1])) {
        MESSAGE(std::setprecision(17) << "v_par = " << v_par << ", lower_[1] = " << lower_[1] <<
                ", upper_[1] = " << upper_[1]);
        THROW("Illegal parameter, must lie between endparameters");
    }
    
    double height = upper_[1] - lower_[1];
    double s = (upper_[1] - v_par) / height; // s + (1 - s) = 1.

    Vector2D return_vector = s*lower_ + (1-s)*upper_;

    return return_vector;
}

//===========================================================================
Vector2D ftGraphEdge::point(double v_par, shared_ptr<ftFaceBase> face)
//===========================================================================
{
    // We first check whether the edge belongs to face.
    size_t ki;
    for (ki = 0; ki < faces_.size(); ++ki)
	if (face == faces_[ki])
	    break;
    ALWAYS_ERROR_IF(ki == faces_.size(),
		"Edge was not a member of input face.");
    Vector2D local_lower_pt = face_params_[ki].first;
    Vector2D local_upper_pt = face_params_[ki].second;

    ALWAYS_ERROR_IF((v_par < lower_[1]) || (v_par > upper_[1]),
		"Illegal parameter, must lie between endparameters.");

    double height = upper_[1] - lower_[1];
    double s = (upper_[1] - v_par) / height; // s + (1 - s) = 1.

    Vector2D return_vector = s*local_lower_pt + (1-s)*local_upper_pt;

    return return_vector;
}

//===========================================================================
const vector<shared_ptr<ftFaceBase> >& ftGraphEdge::getFaces() const
//===========================================================================
{
    return faces_;
}

//===========================================================================
Vector2D ftGraphEdge::lower() const
//===========================================================================
{
    return lower_;
}

//===========================================================================
Vector2D ftGraphEdge::upper() const
//===========================================================================
{
    return upper_;
}

//===========================================================================
ftSearchNode::ftSearchNode(ftSurfaceSetPoint* node)
//===========================================================================
    : node_(node->ftSamplePoint::getPar())
{
}

//===========================================================================
ftSearchNode::~ftSearchNode()
//===========================================================================
{
}

//===========================================================================
Vector2D ftSearchNode::node() const
//===========================================================================
{
    return node_;
}

// //===========================================================================
// const vector<ftFaceBase*>& ftSearchNode::faces() const
// //===========================================================================
// {
//     return faces_;
// }

// //===========================================================================
// double ftSearchNode::u() const
// //===========================================================================
// {
//     return node_[0];
// }

// //===========================================================================
// double ftSearchNode::v() const
// //===========================================================================
// {
//     return node_[1];
// }

//===========================================================================
void ftSearchNode::setOrderedSegments(const vector<ftGraphEdge>& edges)
//===========================================================================
{
    ordered_segments_ = edges;
}


//===========================================================================
const vector<ftGraphEdge>& ftSearchNode::getOrderedSegments() const
//===========================================================================
{
    return ordered_segments_;
}


//===========================================================================
ftPlanarGraph::ftPlanarGraph()
//===========================================================================
{
}

//===========================================================================
ftPlanarGraph::ftPlanarGraph(vector<ftSamplePoint*>& nodes) //ftPointSet& nodes)
//===========================================================================
{
    setGraph(nodes);
}

//===========================================================================
ftPlanarGraph::~ftPlanarGraph()
//===========================================================================
{
}

//===========================================================================
void ftPlanarGraph::setGraph(vector<ftSamplePoint*>& nodes) //ftPointSet& points)
//===========================================================================
{
    nodes_.clear();

    ftSurfaceSetPoint* node;
//     vetor<ftSurfaceSetPoint*> graph_nodes;
    for (size_t ki = 0; ki < nodes.size(); ++ki) {
// 	if (!(nodes[i]->isOnSubSurfaceBoundary()))
// 	    continue; // Graph to be defined by boundary points only.
	ftSamplePoint* dummy_node = nodes[ki];
	node = dynamic_cast<ftSurfaceSetPoint*>(dummy_node);
	ALWAYS_ERROR_IF(node == 0,
		    "Node was not of type ftSurfaceSetPoint.");
// 	graph_nodes.push_back(node);
 	nodes_.push_back(ftSearchNode(node));
    }

//     sort(graph_nodes.begin(), graph_nodes.end(), nodeSortBool);
    sort(nodes_.begin(), nodes_.end(), nodeSortBool);

    // For all nodes in nodes_, create edges ascending through node.
    createOrderedSegments(nodes);
}

//===========================================================================
shared_ptr<ftFaceBase> ftPlanarGraph::locateInGraph(double u, double v) const
//===========================================================================
{
    Vector2D pt(u, v);
    ftGraphEdge left, right;
    double y_lower, y_upper;
    shared_ptr<ftFaceBase> found_face;
    findBoundingTrapezoid(pt, left, right, y_lower, y_upper, found_face);

    return found_face;
}

//===========================================================================
void ftPlanarGraph::getLocalParameters(double& u, double& v,
				       shared_ptr<ftFaceBase>& face) const
//===========================================================================
{
    Vector2D pt(u, v);
    ftGraphEdge left_edge, right_edge;
    double v_lower, v_upper;
    findBoundingTrapezoid(pt, left_edge, right_edge, v_lower, v_upper, face);

    // We find the horizontal line intersecting left end right edges, on which the
    // point lies. We then describe point as a convex combination of the intersection
    // points. Finally we transform calculations to the original parameter domain.
    Vector2D left_pt = left_edge.point(pt[1]);
    Vector2D right_pt = right_edge.point(pt[1]);
    double width = right_pt[0] - left_pt[0];
    // t+(1-t)=1 => t*left_pt+(1-t)*right_pt=pt
    // We make sure that t-value is inside [0,1]
    double t = std::min(1.0, std::max(0.0, (right_pt[0] - pt[0]) / width));

    // Param pts in orig domain.
    Vector2D new_left_pt = left_edge.point(pt[1], face);
    Vector2D new_right_pt = right_edge.point(pt[1], face);
    Vector2D return_vector = t*new_left_pt + (1-t)*new_right_pt;

    // @@sbr201704 If there is a degenerate edge in any of the end points we must carefully select the
    // suitable parameter value along the degenerate edge. We should make sure that a straight line in
    // the global domain maps to a straight line in the local domain.
    // @@sbr201704 This call should be a separate function! Extract when it is complete.
    ftSurface* ft_sf = face->asFtSurface();
    shared_ptr<ParamSurface> sf = ft_sf->surface();
    bool deg_b, deg_r, deg_t, deg_l;
    const double deg_tol = 1e-12;
    sf->isDegenerate(deg_b, deg_r, deg_t, deg_l, deg_tol);

    bool degenerate = (deg_b || deg_t || deg_t || deg_l);
    if (degenerate) {
        // For the degenerate case we compute the convex combination of the corner points, i.e. we find
        // the barycentric coordinates. By looking at the contribution from the non-degenerate corner
        // points we get a well defined parameter value in that direction. The parameter value in the
        // other direction is given by the location on the line from the degenerate point which extends
        // through the non-degenerate edge. We currently only support cases with 1 degenerate edge.

        Vector2D left_edge_start_pt = left_edge.point(left_edge.startparam(), face);
        Vector2D left_edge_end_pt = left_edge.point(left_edge.endparam(), face);
        bool left_edge_start_deg = false;
        bool left_edge_end_deg = false;
            
        Vector2D right_edge_start_pt = right_edge.point(right_edge.startparam(), face);
        Vector2D right_edge_end_pt = right_edge.point(right_edge.endparam(), face);
        bool right_edge_start_deg = false;
        bool right_edge_end_deg = false;
        
        bool left_pt_deg = false;
        bool right_pt_deg = false;
        
        RectDomain dom = sf->containingDomain();
        double knot_tol = 1e-08;
        int num_deg = 0;
        if (deg_b) {
            ++num_deg;

            if (left_edge_start_pt[1] - dom.vmin() < knot_tol)
                left_edge_start_deg = true;
            if (left_edge_end_pt[1] - dom.vmin() < knot_tol)
                left_edge_end_deg = true;
            if (right_edge_start_pt[1] - dom.vmin() < knot_tol)
                right_edge_start_deg = true;
            if (right_edge_end_pt[1] - dom.vmin() < knot_tol)
                right_edge_end_deg = true;

            if (new_left_pt[1] - dom.vmin() < knot_tol)
                left_pt_deg = true;
            if (new_right_pt[1] - dom.vmin() < knot_tol)
                right_pt_deg = true;
        }
        if (deg_r) {
            ++num_deg;

            if (dom.umax() - left_edge_start_pt[0] < knot_tol)
                left_edge_start_deg = true;
            if (dom.umax() - left_edge_end_pt[0] < knot_tol)
                left_edge_end_deg = true;
            if (dom.umax() - right_edge_start_pt[0] < knot_tol)
                right_edge_start_deg = true;
            if (dom.umax() - right_edge_end_pt[0] < knot_tol)
                right_edge_end_deg = true;

            if (dom.umax() - new_left_pt[0] < knot_tol)
                left_pt_deg = true;
            if (dom.umax() - new_right_pt[0] < knot_tol)
                right_pt_deg = true;
        }
        if (deg_t) {
            ++num_deg;

            if (dom.vmax() - left_edge_start_pt[0] < knot_tol)
                left_edge_start_deg = true;
            if (dom.vmax() - left_edge_end_pt[0] < knot_tol)
                left_edge_end_deg = true;
            if (dom.vmax() - right_edge_start_pt[0] < knot_tol)
                right_edge_start_deg = true;
            if (dom.vmax() - right_edge_end_pt[0] < knot_tol)
                right_edge_end_deg = true;

            if (dom.vmax()- new_left_pt[1]< knot_tol)
                left_pt_deg = true;
            if (dom.vmax() - new_right_pt[1] < knot_tol)
                right_pt_deg = true;
        }
        if (deg_l) {
            ++num_deg;

            if (left_edge_start_pt[0] - dom.umin() < knot_tol)
                left_edge_start_deg = true;
            if (left_edge_end_pt[0] - dom.umin() < knot_tol)
                left_edge_end_deg = true;
            if (right_edge_start_pt[0] - dom.umin() < knot_tol)
                right_edge_start_deg = true;
            if (right_edge_end_pt[0] - dom.umin() < knot_tol)
                right_edge_end_deg = true;
            
            if (new_left_pt[0] - dom.umin() < knot_tol)
                left_pt_deg = true;
            if (new_right_pt[0] - dom.umin() < knot_tol)
                right_pt_deg = true;
        }
        if (num_deg > 1) {
            MESSAGE("More than 1 degenerate edge, not supported, expect failure!");
        }
        
//        if (left_pt_deg || right_pt_deg) {
        if (left_edge_start_deg || left_edge_end_deg || right_edge_start_deg || right_edge_end_deg) {
            // std::cout << "left_edge_start_deg: " << left_edge_start_deg << ", left_edge_end_deg: " << left_edge_end_deg <<
            //     ", right_edge_start_deg: " << right_edge_start_deg << ", right_edge_end_deg: " << right_edge_end_deg <<
            // std::endl;
            // std::cout << "left_pt_deg: " << left_pt_deg << ", right_pt_deg: " << right_pt_deg << std::endl;
            // std::cout << "u: " << u << ", v: " << v << std::endl;
            // std::cout << "new_left_pt: " << new_left_pt << ", new_right_pt: " << new_right_pt << std::endl;
            // std::cout << "return_vector: " << return_vector << std::endl;
            
            // If either new_left_pt or new_right_pt are on a degenerate edge the case requires special
            // handling.

            // We must check if the left_pt or right_pt is along a degenerate edge.
            Vector2D pt1, pt2, pt3;
            Vector2D local_pt1, local_pt2, local_pt3;
            double vmin, vmax;
            // The deg point (in global parameter domain) as well as point on the oppsite side along min/max edge.
            Point deg_pt, opp_min, opp_max;
            if (left_edge_start_deg || left_edge_end_deg) {
                vmin = right_edge.startparam();
                vmax = right_edge.endparam();
                pt1 = right_edge.lower();
                pt2 = right_edge.upper();
                pt3 = (left_edge_start_deg) ? left_edge.lower() : left_edge.upper();
                // Then the params in the local surface.
                local_pt1 = right_edge.point(right_edge.startparam(), face);
                local_pt2 = right_edge.point(right_edge.endparam(), face);
                local_pt3 = (left_edge_start_deg) ? left_edge.point(left_edge.startparam(), face) :
                    left_edge.point(left_edge.endparam(), face);
            } else {
                vmin = left_edge.startparam();
                vmax = left_edge.endparam();
                pt1 = left_edge.lower();
                pt2 = left_edge.upper();
                pt3 = (right_edge_start_deg) ? right_edge.lower() : right_edge.upper();
                // Then the params in the local surface.
                local_pt1 = left_edge.point(left_edge.startparam(), face);
                local_pt2 = left_edge.point(left_edge.endparam(), face);
                local_pt3 = (right_edge_start_deg) ? right_edge.point(right_edge.startparam(), face) :
                    right_edge.point(right_edge.endparam(), face);
            }

#if 0
            // Using Cramer's rule we find the barycentric coordinates.
            double area0 = areaTriangle(pt1, pt2, pt3);
            double area1 = areaTriangle(pt2, pt3, pt);
            double area2 = areaTriangle(pt3, pt1, pt);
            double area3 = areaTriangle(pt1, pt2, pt);
            // The parameter in the degenerate direction is given as a convex combination of
            // bar1 & bar2.
            double bar1 = area1/area0;
            double bar2 = area2/area0;
            double bar3 = area3/area0;
            // std::cout << "bar1: " << bar1 << ", bar2: " << bar2 << ", bar3: " << bar3 << std::endl;

            // std::cout << "old u: " << u << ", old v: " << v << std::endl;            
            // We can expect the surface normal of the global and local surf to coincide
            // (approximately). But the parameter directions may be rotated.
            MESSAGE("Relating to full parameter domain, we should relate to the sub-domain given by the trapezoid!");
            if (deg_b || deg_t) {
                // @@sbr201704 Working for this specific case ... Fix!
                //u = (bar1 + bar2 < knot_tol) ? vmin : (vmin*bar1 + vmax*bar2)/(bar1 + bar2);
                u = (bar1 + bar2 < knot_tol) ? dom.umin() : (dom.umin()*bar1 + dom.umax()*bar2)/(bar1 + bar2);
                v = (deg_t) ? dom.vmin() + bar3*(dom.vmax() - dom.vmin()) : dom.vmax() - bar3*(dom.vmax() - dom.vmin());
            } else {
                // @@sbr201704 Working for this specific case ... Fix!
                // v = (bar1 + bar2 < knot_tol) ? vmin : (vmin*bar1 + vmax*bar2)/(bar1 + bar2);
                v = (bar1 + bar2 < knot_tol) ? dom.vmin() : (dom.vmin()*bar1 + dom.vmax()*bar2)/(bar1 + bar2);
                u = (deg_r) ? dom.umin() + bar3*(dom.umax() - dom.umin()) : dom.umax() - bar3*(dom.umax() - dom.umin());
            }
            //MESSAGE("pt: " << pt);
            //MESSAGE("u: " << u << ", v: " << v);
#endif
            
#if 0
            // A second version, slightly improved ...
            if (deg_b || deg_t) {
                // @@sbr201704 Working for this specific case ... Fix!
                //u = (bar1 + bar2 < knot_tol) ? vmin : (vmin*bar1 + vmax*bar2)/(bar1 + bar2);
                u = (bar1 + bar2 < knot_tol) ? dom.umin() : (local_pt1[0]*bar1 + local_pt2[0]*bar2)/(bar1 + bar2);
                v = (deg_t) ? dom.vmin() + bar3*(dom.vmax() - dom.vmin()) : dom.vmax() - bar3*(dom.vmax() - dom.vmin());
            } else {
                // @@sbr201704 Working for this specific case ... Fix!
                // v = (bar1 + bar2 < knot_tol) ? vmin : (vmin*bar1 + vmax*bar2)/(bar1 + bar2);
                v = (bar1 + bar2 < knot_tol) ? dom.vmin() : (local_pt2[1]*bar1 + local_pt2[1]*bar2)/(bar1 + bar2);
                u = (deg_r) ? dom.umin() + bar3*(dom.umax() - dom.umin()) : dom.umax() - bar3*(dom.umax() - dom.umin());
            }
            //MESSAGE("Second try: u: " << u << ", v: " << v);
#endif

            // We find the intersection between the line going from the degenerate point and intersecting
            // the non-degenerate edge. We then use the convex combination of the local parameter values
            // along the edge to find the corresponding local parameter point, giving us the parameter
            // value in the degenerate direction. Finally the parameter in the opposite direction is
            // given by expressing pt as the convex combination of the degenerate global parameter point
            // and the intersection point.
            Vector2D pt_line = pt - pt3;
            Vector2D non_deg_line = pt2 - pt1;
            Vector2D int_pt = intersect2DLines(pt3, pt_line, pt1, non_deg_line);
            //MESSAGE("int_pt: " << int_pt);
            // We express the int_pt as a convex combination of pt1 & pt2 (it should line between the points).
            // t+(1-t)=1 => t*pt1+(1-t)*pt2=int_pt
            const double width0 = pt2[0] - pt1[0];
            const double width1 = pt2[1] - pt1[1];
            const int ind = (fabs(pt2[0] - pt1[0]) > fabs(pt2[1] - pt1[1])) ? 0 : 1;
            const double width =  pt2[ind] - pt1[ind];
            // We allow the t-value to be outside [0,1].
            const double t = (pt2[ind] - int_pt[ind])/width;
            Vector2D int_pt2 = t*pt1 + (1.0-t)*pt2;
            const double conv_dist = int_pt.dist(int_pt2);
            Vector2D local_int_pt = t*local_pt1 + (1.0-t)*local_pt2;
            //MESSAGE("DEBUG: conv_dist: " << conv_dist);
            if (deg_b || deg_t)
            {
                u = local_int_pt[0];
                v = local_pt3[1] + (local_int_pt[1] - local_pt3[1])*(pt.dist(pt3))/(int_pt.dist(pt3));
            }
            else
            {
                MESSAGE("Not yet tested!");
                v = local_int_pt[1];
                u = local_pt3[0] + (local_int_pt[0] - local_pt3[0])*(pt[0] - pt3[0])/(int_pt[0] - pt3[0]);
            }
            //MESSAGE("Third try: u: " << u << ", v: " << v);

#ifndef NDEBUG
            {
                std::ofstream fileout_debug("tmp/pts.g2");
                Point local_pt = face->point(u, v);
                vector<double> pts(local_pt.begin(), local_pt.end());
                PointCloud3D pt_cl(pts.begin(), pts.size()/3);
                pt_cl.writeStandardHeader(fileout_debug);
                pt_cl.write(fileout_debug);
            }
#endif
            
            return;
            
        } // else { // If the deg param pt is not included we should not run into trouble.
            // std::cout << "Neither left_pt or right_pt is degenerate. pt: " << pt
            //           << ", return_vector: " << return_vector << std::endl;
            //std::cout << "Surface degenerate, but not end points of edges! Did not expect this." << std::endl;
//        }      
    }
    
// #if 1 // Debugging
//     double knot_tol = 1e-04; // @@sbr201704 Rather large value ... Required for current case. Fix calling tolerance!
    // We expect the points to share v-value (i.e. lie on a horizontal line). No, not in local domain.
    // if (fabs(new_left_pt[1] - new_right_pt[1]) > knot_tol) {
    //     MESSAGE("ftPlanarGraph::getLocalParameters(): Method seems to have failed!");
    //     std::cout << "new_left_pt: " << new_left_pt << ", new_right_pt: " << new_right_pt << std::endl;
    // }
// #endif
    
    u = return_vector[0];
    v = return_vector[1];
}

//===========================================================================
void ftPlanarGraph::createOrderedSegments(vector<ftSamplePoint*>& nodes) //ftPointSet& points)
//===========================================================================
{
    vector<ftGraphEdge> out_edges; // Edges starting in a given node.
    vector<ftGraphEdge>::iterator edge_iter;

    vector<ftGraphEdge> edges;
    vector<ftGraphEdge> all_edges = createAllGraphEdges(nodes);
//     // debugging
//     std::ofstream fileout("data/output/gnuplot/edge_graph.dat");
//     for (i = 0; i < all_edges.size(); ++i) {
// 	fileout << all_edges[i].startPoint() << " " << all_edges[i].endPoint() <<
// 	    std::endl;
//     }
//     // end of debugging

    for (size_t ki = 0; ki < nodes_.size(); ++ki) {
	Vector2D from = nodes_[ki].node();
	out_edges.clear(); // We clean up.
	int nmb_in_edges = 0; // For each node: = # neighbours - # up_edges - # horizontal edges.
	ftSurfaceSetPoint* found_pt;
 	vector<PointIter> neighbours =
	    getNeighbours(nodes_[ki].node(), nodes, found_pt);

	// We extract the neighbours which have higher y-value (out edges).
 	for (size_t kj = 0; kj < neighbours.size(); ++kj)
// 	    if (!(neighbours[kj]->isOnSubSurfaceBoundary()))
// 		continue;
	    if (nodes_[ki].node()[1] < neighbours[kj]->getPar()[1]) {
		Vector2D to = neighbours[kj]->getPar();
		out_edges.push_back(findEdge(all_edges, from, to));
	    } else if (nodes_[ki].node()[1] > neighbours[kj]->getPar()[1]) {
		++nmb_in_edges;
 	    } //else {
// 		GO_ERROR("Unexpected incident!", UnknownError());
//  	    }

	// We must sort the out_edges (which all start in the same point),
	sort(out_edges.begin(), out_edges.end(), edgeSortBool);

	if (ki != 0) {
	    edges = nodes_[ki-1].getOrderedSegments();
	    ASSERT(edges.size() != 0);
	    edgeRuntimeBool ang_bool(from); // Angular bool (tests: angle > 0).
	    edge_iter = find_if(edges.begin(), edges.end(), ang_bool);
	    if (nmb_in_edges > 0)
		// @@sbr We're in trouble if parametrization has lead to overlapping edges...
		ALWAYS_ERROR_IF((edge_iter == edges.end()) ||
			    ((edges.end() - edge_iter) < nmb_in_edges),
			    "This should never happen!");
            int pos = edge_iter - edges.begin();
	    edges.erase(edge_iter, edge_iter + nmb_in_edges);
	    edges.insert(edges.begin() + pos, out_edges.begin(), out_edges.end());
	} else
	    edges = out_edges;

	if (edges.size() == 0) {
// 	    MESSAGE("Node contained no out-edges, remove it from vector?");
// 	    nodes_.erase(nodes_.begin() + i);
// 	    --i;
	} else {
	    nodes_[ki].setOrderedSegments(edges);
	}
    }

    // Almost done, however if equal y-value we remove the prior node.
    for (size_t ki = nodes_.size() - 1; ki > 0; --ki)
	if (nodes_[ki-1].node()[1] == nodes_[ki].node()[1]) {
// 	    nodes_[ki-1].setOrderedSegments(nodes_[ki].getOrderedSegments());
	    nodes_.erase(nodes_.begin() + (ki - 1));
	}

}

//===========================================================================
vector<ftGraphEdge> ftPlanarGraph::createAllGraphEdges(vector<ftSamplePoint*>& nodes) //ftPointSet& points)
//===========================================================================
{
//     GraphIter neighbour;
    vector<ftGraphEdge> return_edges;
    GraphIter iter = nodes_.begin();
    while (iter != nodes_.end()) {
	ftSurfaceSetPoint* found_pt;
	ftSurfaceSetPoint* neighbour;
	vector<PointIter> neighbours =
	    getNeighbours(iter->node(), nodes, found_pt);
	for (size_t kj = 0; kj < neighbours.size(); ++kj) {
// 	    if ((neighbours[kj]->isOnSubSurfaceBoundary()) &&
	    if (iter->node()[1] < neighbours[kj]->getPar()[1]) {
// 		nodeRuntimeBool node_bool(neighbours[kj]);
// 		neighbour = find_if(nodes_.begin(), nodes_.end(), node_bool);
		neighbour = dynamic_cast<ftSurfaceSetPoint*>(neighbours[kj]);
		ALWAYS_ERROR_IF(neighbour == 0,
				"Found neighbour was not of type ftSurfaceSetPoint!");

		return_edges.push_back(ftGraphEdge(found_pt, neighbour));
	    }
	}
	++iter;
    }

    return return_edges;
}

//===========================================================================
vector<PointIter> ftPlanarGraph::getNeighbours(const Vector2D& node,
					       vector<ftSamplePoint*>& nodes, //ftPointSet& points,
					       ftSurfaceSetPoint*& found_pt)
//===========================================================================
{
    size_t ki;
    for (ki = 0; ki < nodes.size(); ++ki) {
	// Maybe we should allow numerical error?
	if ((node - nodes[ki]->getPar()).length() < knot_tol) {
	    found_pt = dynamic_cast<ftSurfaceSetPoint*>(nodes[ki]);
	    ALWAYS_ERROR_IF(found_pt == 0,
			    "Member of nodes was not of type ftSurfaceSetPoint.");

	    break;
	}
    }
    if (ki == nodes.size()) {
	THROW("Point was not to be found. Should never happen.");
    }

    // We find all nodes which are connected to node.
    // I.e. for all bd nodes we continue until we find the first member
    vector<PointIter> neighbours = found_pt->getNeighbours();
    vector<PointIter> graph_neighbours;
    for (ki = 0; ki < neighbours.size(); ++ki) {
	if (neighbours[ki]->isOnSubSurfaceBoundary()) {
	    bool forward = (neighbours[ki]->getFirstNeighbour() != found_pt);
	    PointIter curr_iter = neighbours[ki];
	    PointIter graph_neighbour = NULL;
	    while (graph_neighbour == NULL) {
		size_t kj;
		for (kj = 0; kj < nodes.size(); ++kj) {
		    if (curr_iter == nodes[kj]) {
			graph_neighbour = nodes[kj];
			break;
		    }
		}
		if (kj == nodes.size()) {
		    vector<PointIter> local_neighbours = curr_iter->getNeighbours();
		    curr_iter = forward ? local_neighbours[0] :
			local_neighbours[local_neighbours.size() - 1];
		    ASSERT(curr_iter != found_pt);
		}
	    }
	    graph_neighbours.push_back(graph_neighbour);
	}
    }

    ASSERT(graph_neighbours.size() > 0);
    return graph_neighbours;
}

//===========================================================================
ftGraphEdge ftPlanarGraph::findEdge(const vector<ftGraphEdge>& edges,
				    Vector2D from, Vector2D to)
//===========================================================================
{
    for (size_t ki = 0; ki < edges.size(); ++ki) {
	if (((edges[ki].startPoint() - from).length() == 0)
	    && ((edges[ki].endPoint() - to).length() == 0))
	    return edges[ki];
    }

    THROW("Edge not to be found!");
}

// //===========================================================================
// void ftPlanarGraph::sortOutEdges(vector<ftGraphEdge>& edges) const
// //===========================================================================
// {
//     sort(edges.begin(), edges.end(), edgeSortBool);
// }

//===========================================================================
void ftPlanarGraph::findBoundingTrapezoid(Vector2D& pt,
					  ftGraphEdge& left, ftGraphEdge& right,
					  double& y_lower, double& y_upper,
					  shared_ptr<ftFaceBase>& face) const
//===========================================================================
{
    if (nodes_.size() == 0)
    {
        THROW("The bounding trapezoid was not set, this function should not have been called!");
    }

    // We start by moving parameter v inside upper and lower edges, if close.
    double vmin = nodes_[0].node()[1]; // nodes_ are ordered by increasing v-value.
    double vmax = nodes_[nodes_.size()-1].node()[1];
    if ((pt[1] < vmin) || (pt[1] > vmax)) {
        if (fabs(pt[1] - vmin) < knot_tol)
            pt = Vector2D(pt[0], vmin);
        else if (fabs(pt[1] - vmax) < knot_tol)
            pt = Vector2D(pt[0], vmax);
        else {
            THROW("Searching for point outside graph.");
        }
    }

    Vector3D dummy_pt;
    ftSurfaceSetPoint dummy_par(dummy_pt, -1);
    dummy_par.setPar(pt);
    ftSearchNode dummy(&dummy_par);
    vector<ftSearchNode>::const_iterator v_strip = // We find containing strip.
	upper_bound(nodes_.begin(), nodes_.end(), dummy, nodeSortBool);

    if (v_strip == nodes_.end()) {
	THROW("This should never happen!");
    }

    // If found strip is an upper bound, we choose previous strip.
    if (v_strip != nodes_.begin()) {
	--v_strip;
    }

    y_lower = v_strip->node()[1];
    y_upper = v_strip[1].node()[1]; // Set return value y_lower.

    edgeRuntimeBool edge_bool(pt);
    ftGraphEdge dummy_edge; // Not to be used, but required for formal reasons.
    vector<ftGraphEdge> edges = v_strip->getOrderedSegments();
    if (edges.size() == 0) {
	THROW("Should never happen!");
    }
    vector<ftGraphEdge>::const_iterator right_edge = // ++left_edge == right_edge.
	upper_bound(edges.begin(), edges.end(), dummy_edge, edge_bool);

    // To handle points on bounding edges, we may have to iterate.
    if (right_edge == edges.end()) {
	--right_edge;
    } else if (right_edge == edges.begin()) {
	++right_edge;
    }

    left = right_edge[-1]; // Setting return left.
    right = right_edge[0]; // Setting return right.
    vector<shared_ptr<ftFaceBase> > left_faces = left.getFaces();
    vector<shared_ptr<ftFaceBase> > right_faces = right.getFaces();
    vector<shared_ptr<ftFaceBase> > common_faces;
    for (size_t ki = 0; ki < left_faces.size(); ++ki) {
	for (size_t kj = 0; kj < right_faces.size(); ++kj) {
	    if (left_faces[ki] == right_faces[kj]) {
		common_faces.push_back(left_faces[ki]);
	    }
	}
    }

    if (common_faces.size() == 1) {
	face = common_faces[0];
    }else if (common_faces.size() == 2) {
	// It appears the two edges belong to two identical patches. We use the
	// right_edge iterator and reductio ad absurdum to decide which one it is.
	// Assuming a face is the one, and iterating towards the bnd, we search for
	// an edge contradicting the assumption. Or else we made the right guess!
	vector<ftGraphEdge>::const_iterator next_edge = right_edge + 1;
	vector<shared_ptr<ftFaceBase> > next_faces = next_edge->getFaces();
	while ((next_faces.size() != 1) && (
	       ((common_faces[0] == next_faces[0]) &&
		(common_faces[1] == next_faces[1])) ||
	       ((common_faces[1] == next_faces[0]) &&
		(common_faces[0] == next_faces[1])))) {
	    ++next_edge;
	    next_faces = next_edge->getFaces();
	    ALWAYS_ERROR_IF(next_faces.size() == 0, "This should never happen!");
	}

	// We are either at the boundary, or have met a new face; we summarize:
	int nmb_steps = (int)(next_edge - right_edge);
	shared_ptr<ftFaceBase> outer_face, inner_face;

	if ((common_faces[0] == next_faces[0]) ||
	    ((next_faces.size() > 1) && (common_faces[0] == next_faces[1]))) {
	    outer_face = common_faces[0];
	    inner_face = common_faces[1];
	} else {
	    outer_face = common_faces[1];
	    inner_face = common_faces[0];
	}

	if (nmb_steps % 2) // True if odd number of steps.
	    face = inner_face;
	else
	    face = outer_face;
    } else
	THROW("This should never happen!");

#ifndef NDEBUG
    {
        std::ofstream fileout_debug("tmp/graph_edges.g2");

        Vector2D left_lower = left.point(left.startparam(), face);
        Vector2D left_upper = left.point(left.endparam(), face);
        Point left_low = face->point(left_lower[0], left_lower[1]);
        Point left_high = face->point(left_upper[0], left_upper[1]);
        vector<double> left_pts(left_low.begin(), left_low.end());
        left_pts.insert(left_pts.end(), left_high.begin(), left_high.end());
        LineCloud line_cl(left_pts.begin(), left_pts.size()/6);
        line_cl.writeStandardHeader(fileout_debug);
        line_cl.write(fileout_debug);

        Vector2D right_lower = right.point(right.startparam(), face);
        Vector2D right_upper = right.point(right.endparam(), face);
        Point right_low = face->point(right_lower[0], right_lower[1]);
        Point right_high = face->point(right_upper[0], right_upper[1]);
        vector<double> right_pts(right_low.begin(), right_low.end());
        right_pts.insert(right_pts.end(), right_high.begin(), right_high.end());
        LineCloud line_cl2(right_pts.begin(), right_pts.size()/6);
        line_cl2.writeStandardHeader(fileout_debug);
        line_cl2.write(fileout_debug);
    }
#endif

}


//===========================================================================
double areaTriangle(const Vector2D& corner1, const Vector2D& corner2, const Vector2D& corner3)
//===========================================================================
{
    // We use Heron's formula to compute the area.
    double a = corner1.dist(corner2);
    double b = corner2.dist(corner3);
    double c = corner3.dist(corner1);
    double s = 0.5*(a + b + c); // The half-sum.
    double prod = s*(s - a)*(s - b)*(s - c);
    double tol = 1.0e-08; // To avoid nan.
    double area = (prod < tol) ? 0.0 : sqrt(prod);

    return area;
}


Vector2D intersect2DLines(Vector2D from1, Vector2D dir1, Vector2D from2, Vector2D dir2)
{
    Vector2D int_pt;
    const double ang = dir1.angle(dir2);
    const double ang_tol = 1.0e-04;
    if (ang < ang_tol) // The lines are parallell.
    {
        return int_pt;
    }
    
    // from1 + t1*dir1 = from2 + t2*dir2

    // from1[0] + t1*dir1[0] = from2[0] + t2*dir2[0]
    // from1[1] + t1*dir1[1] = from2[1] + t2*dir2[1]

    // t1 = ((from2[0] + t2*dir2[0] - from1[0])/dir1[0])
    // t2*dir1[1]*dir2[0]/dir1[0] - t2*dir2[1] = from2[1] - from1[1] + ((from1[0] - from2[0])/dir1[0])*dir1[1]
    if (fabs(dir1[0]) > fabs(dir2[0]))
    {
        double t2 = (from2[1] - from1[1] + ((from1[0] - from2[0])/dir1[0])*dir1[1])/((dir1[1]*dir2[0]/dir1[0]) - dir2[1]);
        int_pt = from2 + t2*dir2;
    }
    else
    {
        double t1 = (from1[1] - from2[1] - ((from1[0] + from2[0])/dir2[0])*dir2[1])/((dir2[1]*dir1[0]/dir2[0]) - dir1[1]);
        int_pt = from1 + t1*dir1;
    }

    return int_pt;
}

}
