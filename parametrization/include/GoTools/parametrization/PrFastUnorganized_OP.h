/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1997 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#ifndef PRFASTUNORGANIZED_OP_H
#define PRFASTUNORGANIZED_OP_H

#include "GoTools/parametrization/PrOrganizedPoints.h"
#include "GoTools/parametrization/PrCellStructure.h"
#include "GoTools/utils/Array.h"
#include "GoTools/utils/errormacros.h"
using Go::Vector2D;

/*<PrFastUnorganized_OP-syntax: */

/** - PrFastUnorganized_OP represents a set of points in three dimensions.
 * All we know about them is their boundary in sequence.
 * No other topology is given.
 * This kind of data can be parametrized by PrParametrizeBdy and PrParametrizeInt.
 * The major difference to the class PrUnorganized_OP
 * is that this class stores the neighbours for each vertex.
 * Note that this requires the call of "initNeighbours" first! 
 */
class PrFastUnorganized_OP : public PrOrganizedPoints
{
private:
  PrCellStructure cellstruct_;
  vector<Vector2D> uv_;
  int nInt_; // no of interior points,
             // so xyz_(0),...,xyz_(nInt_-1) are the interior points
             // and xyz_(nInt_),...,xyz_(xyz_.size()-1) are the boundary points
             // in an anticlockwise order.
  double radius2_; // square of radius of epsilon ball.
  int knearest_; // number of points in neighbourhood (alternative to radius)
  int use_k_; // = 1 use k nearest, = 0 use fixed radius

  // additional function and memory to store neighbours explicitly
  void findNeighbours(int i, vector<int>& neighbours) const;
  vector< vector<int> > nbrs;

public:
  /// Constructor
  PrFastUnorganized_OP(int num_cells = 10);
  /// Constructor
  PrFastUnorganized_OP(int n, int n_int, double* xyz_points, int num_cells = 10);
  /// Destructor
  ~PrFastUnorganized_OP() {};

  /// Compute (and internally store) information about the neighbours of each 
  /// point.  Neighbours are defined either by the \em k nearest points or by
  /// all points within a specified radius.  Use the member functions useK() and
  /// useRadius() to specify this.
  void initNeighbours() ;

  //             Derived from base class
  /// Return the number of nodes in the graph.
  virtual int       getNumNodes() const  {return cellstruct_.getNumNodes(); }
  /// Return the i-th node in the graph if the nodes are three-dimensional.
  virtual Vector3D get3dNode(int i) const {return cellstruct_.get3dNode(i); }
  /// This inherited member function does not apply to this derived class.
  /// As of now, it will throw if you try to run it.
  virtual void        set3dNode(int i, const Vector3D& p)
    {
	THROW("PrFastUnorganized_OP::set3dNode() doesn't fit in, not implemented.");

    }

  /// Return the indices of the neighbours of the i-th node.
  /// (here there is no ordering: it's not a planar graph.)
  virtual void       getNeighbours(int i, vector<int>& neighbours) const;
  /// Is 'i' a boundary node or interior node?
  virtual bool    isBoundary(int i) const {return (i >= nInt_ ? 1 : 0);}

  virtual double getU(int i) const {return uv_[i].x(); }
  virtual double getV(int i) const {return uv_[i].y(); }
  virtual void setU(int i, double u) {uv_[i].x() = u; }
  virtual void setV(int i, double v) {uv_[i].y() = v; }

  /// get the specified radius for defining neighbours by distance
  double getRadius()
  {
#ifdef __BORLANDC__
    return std::sqrt(radius2_);
#else
    return sqrt(radius2_);
#endif
  }
  /// specify the radius for defining neighbourhood by distance
  void setRadius(double radius) {radius2_ = radius * radius; }
  /// get the specified number for defining neighbourhoods by a fixed
  /// number of neighbours
  int getKNearest() { return knearest_; }
  /// specify a number for defining neighbourhoods by a fixed number
  /// of neighbours.
  void setKNearest(int knearest) {knearest_ = knearest; }
  /// Use a fixed number of neighbours to define the neighbourhood of a 
  /// point.
  void useK() {use_k_ = 1; }
  /// Use an euclidean distance to define the neighbourhood of a point.
  void useRadius() {use_k_ = 0; }

  /// Write contents to stream
  void print(ostream& os);
  //void scan(istream& is);
  // print and scan raw data (no uv's)
  //void printRawData(ostream& os);

  /// Read in points from a stream and regenerate internal structures (no 
  /// parameterization of the points).
  /// \param is The stream containing the points.  The format of the stream
  ///           must be:
  ///           <ul>
  ///           <li> first: the total number of points </li>
  ///           <li> then: the number of interior points (not boundary points) </li>
  ///           <li> then: all internal point coordinates listed as xyz xyz, etc. </li>
  ///           <li> finally: the boundary point coordinates </li>
  ///           </ul>
  /// \param num_cells number of cells in each spatial direction when dividing up
  ///                  the space containing the points.
  /// \param noise magnitude of noise added to the internal points.  Zero by default. 
  void scanRawData(istream& is, int num_cells = 10, double noise = 0.0);
};

/*>PrFastUnorganized_OP-syntax: */

/*Class:PrFastUnorganized_OP

Name:              PrFastUnorganized_OP
Syntax:	           @PrFastUnorganized_OP-syntax
Keywords:
Description:       This class represents a set of poits in three dimensions.
                   All we know about them is their boundary in sequence.
                   No other topology is given.
                   This kind of data can be parametrized by PrParametrize.

		   The major difference to the class PrUnorganized_OP
		   is that this class stores the neighbours for each
		   vertex.
		   Note that this requires the call of "initNeighbours"
		   first!

Member functions:
                   "getNumNodes()" --\\
                   Return the number of nodes in the graph.
                 
                   "get3dNode(int i)" --\\
                   Return the i-th node in the graph if the nodes are 
                   three-dimensional

                   "getNeighbours(int i, PrListInt& neighbours)" --\\
                   Return the indices of the neighbours of the i-th node
                   (here there is no ordering: it's not a planar graph).

                   "isBoundary(int i)" --\\
                   Is i a boundary node or interior node?

Constructors:
Files:
Example:


See also:
Developed by:      SINTEF Applied Mathematics, Oslo, Norway
Author:	           Kai Hormann and Michael Floater, SINTEF
Date:              Feb. 2001
*/

#endif // PRFASTUNORGANIZED_OP_H
