//===========================================================================
//                                                                           
// File: ftPlanarGraph.C                                                     
//                                                                           
// Created: Mon Jan 28 12:21:37 2002                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@math.sintef.no>
//                                                                           
// Revision: $Id: ftPlanarGraph.C,v 1.19 2005-06-09 07:14:33 oan Exp $
//                                                                           
// Description: Implementation file for the classes ftSearchNode, ftGraphEdge
//              & ftPlanarGraph.
//                                                                           
//===========================================================================

#include <algorithm>
#include "GoTools/compositemodel/ftPlanarGraph.h"
#include <fstream>
#include "GoTools/compositemodel/ftSurfaceSetPoint.h"

using std::make_pair;
using namespace Go;

const double knot_tol = 1e-18;

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

    // We find the horizontal line interseting left end right edges, on which the
    // point lies. We then describe point as a convex combination of the intersection
    // points. Finally we transfrom calculations to the original parameter domain.
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

#if 1 // Debugging
    double knot_tol = 1e-04; // Rather large value ... Required for current case. Fix calling tolerance!
    // We expect the points to share v-value (i.e. lie on a horizontal line).
    if (fabs(new_left_pt[1] - new_right_pt[1]) > knot_tol) {
        MESSAGE("ftPlanarGraph::getLocalParameters(): Method seems to have failed!");
        std::cout << "new_left_pt: " << new_left_pt << ", new_right_pt: " << new_right_pt << std::endl;
    }
#endif
    
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
		int pos = edge_iter- edges.begin();
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
}

