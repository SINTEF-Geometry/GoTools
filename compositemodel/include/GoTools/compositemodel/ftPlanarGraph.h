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

#ifndef _FTPLANARGRAPH_H
#define _FTPLANARGRAPH_H

#include <utility>
#include <map>
#include <set>
#include "GoTools/compositemodel/ftFaceBase.h"
#include "GoTools/compositemodel/ftSurfaceSetPoint.h"


namespace Go
{


  class ftSearchNode;
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
  typedef std::vector<ftSearchNode>::iterator GraphIter;
#else
  typedef std::vector<ftSearchNode>::iterator GraphIter;
#endif


  /// Edge is defined by two end points. The face(s) to which it belongs is also known.
  /// y-value of startpoint() is less than that of endpoint().
  class ftGraphEdge
  {

  public:

    /// Constructor
    ftGraphEdge();
    //// Constructor
    ftGraphEdge(ftSurfaceSetPoint*& lower, ftSurfaceSetPoint*& upper);

    /// Destructor
    ~ftGraphEdge();

    double startparam() const;

    double endparam() const;

    Vector2D startPoint() const;

    Vector2D endPoint() const;
      
    Vector2D point(double v_par);

    /// Return parameter point in domain of face, given by its's v_par on graph edge.
    Vector2D point(double v_par, shared_ptr<ftFaceBase> face);

    /// Return the face(s) to which the edge belongs.
    const std::vector<shared_ptr<ftFaceBase> >& getFaces() const;

    // Return the lower node of the edge in the graph.
    Vector2D lower() const;

    // Return the upper node of the edge in the graph.
    Vector2D upper() const;

  private:

    // Value y of lower_ is less than value y of upper_.
    Vector2D lower_; // Parameter value in the graph.
    Vector2D upper_; // Parameter value in the graph.

    std::vector<shared_ptr<ftFaceBase> > faces_; // One or two faces (graph is planar).

    // Parameter values of (lower_, upper_) in faces_. Used for computing local parameters.
    std::vector<std::pair<Vector2D, Vector2D> > face_params_;

  };


  /// Node inside a graph. A ftSearchNode contains all edges in graph crossing
  /// horisontal y0-line (point = (x0, y0)), ordered from left to right, enabling
  /// fast searching of points inside graph. We also store ftFaceBase's of point.
  class ftSearchNode
  {

  public:

    /// Constructor
    ftSearchNode(ftSurfaceSetPoint* node);

    /// Destructor
    ~ftSearchNode();

    /// As we allow to set param value of point, function nor return val are const.
    Vector2D node() const;
    //     ftSurfaceSetPoint* node();

    void setOrderedSegments(const std::vector<ftGraphEdge>& edges);

    const std::vector<ftGraphEdge>& getOrderedSegments() const;


  private :

    Vector2D node_;
    //     ftSurfaceSetPoint* node_;
    std::vector<ftGraphEdge> ordered_segments_; // Vector is sorted => binary_search.

  };


  /// The chosen structure of graph is, in essense, taken from
  /// "Computational Geometry, 2.2.2.1, F. P. Preparata & M. I. Shamos, 1985",
  /// and is well suited for fast searching of points inside graph.
  class ftPlanarGraph
  {

  public:

    /// Constructor
    ftPlanarGraph();
    /// Warning! points are really of type ftSurfaceSetPoint! We cast...
    /// We expect that every pair of points differ in their parameter values.
    ftPlanarGraph(std::vector<ftSamplePoint*>& nodes); //ftPointSet& nodes);

    /// Destructor
    ~ftPlanarGraph();

    /// Warning! points are really of type ftSurfaceSetPoint! We cast...
    /// We expect that every pair of points differ in their parameter values.
    /// We only use those points which are on a subsurface-boundary.
    void setGraph(std::vector<ftSamplePoint*>& nodes); //ftPointSet& points);

    /// Given input of parameter points, return the corresponding ftFaceBase.
    shared_ptr<ftFaceBase> locateInGraph(double u, double v) const;

    /// Return local parameter values in corresponding patch. We locate a bounding
    /// trapezoid inside a patch, and then make a qualified guess from corner values.
    void getLocalParameters(double& u, double& v, shared_ptr<ftFaceBase>& face) const;


  private:

    // All nodes of the graph are sorted after increasing y-value.
    // @@ Should we also sort after increasing x-value? y-value may be equal...
    std::vector<ftSearchNode> nodes_;

    // Private member functions.

    // nodes_ has been constructed, except for edges. All topological information
    // is included in points.
    // Using sorting criterion given by edgeSortCriterion, sort edges for each strip.
    void createOrderedSegments(std::vector<ftSamplePoint*>& nodes); //ftPointSet& points);

    // We create the edges, given y-value is higher. Used to set values in nodes_[i].
    // We do not create horizontal edges.
    // @@ ?? If equal, we demand higher x-value.
    std::vector<ftGraphEdge> createAllGraphEdges(std::vector<ftSamplePoint*>& nodes); //tPointSet& points);

    // Expecting parameter points in points are unique.
    std::vector<PointIter> getNeighbours(const Vector2D& node,
					 std::vector<ftSamplePoint*>& nodes, //ftPointSet& points,
					 ftSurfaceSetPoint*& found_pt);

    // Utility function for finding an edge in edges, given it's unique endpoints.
    ftGraphEdge findEdge(const std::vector<ftGraphEdge>& edges,
			 Vector2D from, Vector2D to);

    //     // We sort the out edges, emerging from the same node, clockwise.
    //     /// Preprocessing.
    //     void sortOutEdges(vector<ftGraphEdge>& edges) const;

    // As we tolerate the input pt to lie outside domain (within num_tol), the value of pt may change.
    void findBoundingTrapezoid(Vector2D& pt,
			       ftGraphEdge& left, ftGraphEdge& right,
			       double& y_lower, double& y_upper,
			       shared_ptr<ftFaceBase>& face) const;


  };

    double areaTriangle(const Vector2D& corner1, const Vector2D& corner2, const Vector2D& corner3);

} // namespace Go


#endif // _FTPLANARGRAPH_H

