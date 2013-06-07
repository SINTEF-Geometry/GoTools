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

#ifndef PRRECTANGULARGRID_OP_H
#define PRRECTANGULARGRID_OP_H

#include "GoTools/parametrization/PrExplicitConnectivity.h"
#include "GoTools/utils/Array.h"
using Go::Vector3D;
using Go::Vector2D;

/*<PrRectangularGrid_OP-syntax: */

/** PrRectangularGrid_OP - This class represents a (topologically) rectagular grid
 * of points in \f$R^3\f$ with or without parameter points in \f$R^2\f$.
 * It implements the virtual functions in PrOrganizedPoints specific to a grid.
 * This kind of grid can be parametrized by PrParametrize. 
 */
class PrRectangularGrid_OP : public PrExplicitConnectivity
{
private:
  //  Regular grid of size m x n.
  vector<Vector3D> grid_;
  vector<Vector2D> uvgrid_;
  int ncols_; // Number of columns
  int nrows_; // Number of rows
  // The grid is m * n with nodes
  // (P_{i,j})_{i=0,...m-1, j=0,...,n-1}. The nodes are stored
  // in the array grid_ in the order:
  // P_{00},P_{01},...,P{0(m-1)},P_{10},...,P_{(m-1)(n-1)}.

public:
  /// Default constructor
  PrRectangularGrid_OP() : ncols_(0), nrows_(0) {}

  /// Constructor. Construct the PrRectangularGrid_OP from an array of nodes.
  PrRectangularGrid_OP(int numcols, int numrows, const double *grid);

  /// Constructor. Construct the PrRectangularGrid_OP from an array of nodes.
  PrRectangularGrid_OP(int numcols, int numrows,
		       const double *grid, const double *uvgrid);

  /// Default destructor
  virtual ~PrRectangularGrid_OP();

  /// compute the index of the point located at column 'col' and row 'row'.
  int  gridToGraph(int col, int row) const
	{ return row*ncols_ + col; }
  /// Get the column and row position of the point indexed 'index'.
  void graphToGrid(int index, int& col, int& row) const
	{ col = index%ncols_; row = (index-col)/ncols_; }

  /** @name Derived from base class */
  //@{
    virtual int       getNumNodes() const {return (int)grid_.size(); }
  virtual Vector3D get3dNode(int i) const {return grid_[i]; }
  virtual void        set3dNode(int i, const Vector3D& p) { grid_[i] = p; }
  virtual void      getNeighbours(int i, vector<int>& neighbours) const;
  virtual bool    isBoundary(int i) const;

  virtual double getU(int i) const {return uvgrid_[i].x(); }
  virtual double getV(int i) const {return uvgrid_[i].y(); }
  virtual void setU(int i, double u) {uvgrid_[i].x() = u; }
  virtual void setV(int i, double v) {uvgrid_[i].y() = v; }
  //@}

  /** @name Other routines */
  //@{
  /// Set the dimensions of the grid.
  void setDim(int numcols, int numrows);
  /// Reset all the xyz points from an array
  void setXYZVertices(const double *grid);
  /// Reset all the uv points from an array
  void setUVVertices(const double *grid);
  /// Get the dimensions of the grid
  void getDim(int& numcols, int& numrows) const
  { numcols = ncols_; numrows = nrows_; }

  /// Get U value for gridpoint (i,j)    
  double getU(int i, int j) const {return uvgrid_[i + j*ncols_].x(); }
  /// Get V value for gridpoint (i,j)    
  double getV(int i, int j) const {return uvgrid_[i + j*ncols_].y(); }
    
  /// Get the indices of the four corner nodes in an anticlockwise
  /// direction, starting with the bottom left hand corner.
  void getCorners(int& c1, int& c2, int& c3, int& c4) const;

  /// Locate the triangle in the parameter plane that contains the point (u,v).
  /// Report back the indexes of the "lower left" node of this triangle, as well
  /// as the point's barycentric coordinates in this triangle
  /// \param u u-coordinate of the specified point
  /// \param v v-coordinate of the specified point 
  /// \retval i if triangle found, this will be the row index of the node in the 
  ///           "lower left" corner of the triangle.
  /// \retval j if triangle found, this will be the column index of the node in
  ///           the "lower left" corner of the triangle.
  /// \retval right if the triangle was found, this variable reports wether the
  ///               triangle consist of the "right" or "left" part of a quad 
  ///               in the grid.  (Where the right part is defined as (i, j), (i+1, j)
  ///               (i+1, j+1), whereas the left part is defined as (i, j), (i, j+1)
  ///               (i+1, j+1).
  /// \retval tau0 if triangle found, this parameter is set to the first barycentric
  ///              coordinate of the point (u,v) in the triangle.
  /// \retval tau1 if triangle found, this parameter is set to the second barycentric
  ///              coordinate of the point (u,v) in the triangle.
  /// \retval tau2 if triangle found, this parameter is set to the third barycentric
  ///              coordinate of the point (u,v) in the triangle.
  /// \return 'true' if the triangle was located, 'false' otherwise
  bool getUVTriangle(double u, double v, int& i, int&j, bool& right,
                     double& tau0, double& tau1, double& tau2) const;

  /// This function behaves as the other getUVTriangle() function, with the difference
  /// that the user can specify two additional parameters, 'ilast' and 'jlast' as a 
  /// suggestion for where to start the search for a triangle.
  bool getUVTriangle(double& u, double& v,
                     int ilast, int jlast,
                     int& i, int&j, bool& right,
                     double& tau0, double& tau1, double& tau2) const;


  /// Check if the parameter point (u,v) can be found in one of the two
  /// triangles specified by the nodes at indexes (ii, jj), (ii+1, jj),
  /// (ii, jj+1) and (ii+1, jj+1).  In that case, return its barycentric
  /// coordinates with respect to that triangle
  /// \param u u-coordinate of the specified point
  /// \param v v-coordinate of the specified point
  /// \param ii row-index of the "lower left" node of the rectangle specifying
  ///           the two triangles where we want to look for the point
  /// \param jj column-index of the "lower left" node of the rectangle
  ///           specifying the two triangles where we want to look for the point
  /// \retval right if the point was found to be in the triangle specified
  ///               by the nodes (ii,jj), (ii+1, jj) and (ii+1, jj+1), then the 
  ///               value of 'right' will be set to 'true'.  It the point was
  ///               found to be in the triangle specified by the nodes (ii, jj),
  ///               (ii, jj+1) and (ii+1, jj+1), then the value of 'right' will
  ///               be set to 'false'.  Otherwise, the value of this parameter
  ///               will remain unchanged
  /// \retval tau0 first barycentric coordinate of the point in the located triangle
  /// \retval tau1 second barycentric coordinate of the point in the located triangle
  /// \retval tau2 third barycentric coordinate of the point in the located triangle
  /// \return 'true' if the parameter point was found to be within one of the two
  ///                triangles.
  bool getTriangle(double u, double v,
		   int ii, int jj, bool& right,
		   double& tau0, double& tau1, double& tau2) const;
  //@}

  //print and scan routines
  /// print object to stream
  void print(std::ostream& os);
  /// read object from stream
  void scan(std::istream& is);
  // print and scan raw data (no uv's)
  /// print raw data to stream (no parameterization involved)
  void printRawData(std::ostream& os);
  /// scan raw data from stream (no parameterization involved)
  void scanRawData(std::istream& is);

    double* uvData()
    {
	return uvgrid_[0].begin();
    }
    double* xyzData()
    {
	return grid_[0].begin();
    }
    const double* uvData() const
    {
	return uvgrid_[0].begin();
    }
    const double* xyzData() const
    {
	return grid_[0].begin();
    }
};


/*>PrRectangularGrid_OP-syntax: */

/*Class:PrRectangularGrid_OP

Name:              PrRectangularGrid_OP
Syntax:	           @PrRectangularGrid_OP-syntax
Keywords:
Description:       This class represents a (topologically) rectagular grid
                   of points in R^3 with or without parameter points in R^2.
                   It implements the virtual functions in PrOrganizedPoints
                   specific to a grid.
                   This kind of grid can be parametrized by PrParametrize.
Member functions:
                   "getNumNodes()" --\\
                   Return the number of nodes in the graph.
                 
                   "get3dNode(int i)" --\\
                   Return the i-th node in the graph if the nodes are 
                   three-dimensional

                   "getNeighbours(int i, PrListInt& neighbours)" --\\
                   Return the indices of the neighbours of the i-th node in:
                   1. any anticlockwise order if i is an interior node
                   2. the unique anticlockwise order if i is a boundary node.

                   "isBoundary(int i)" --\\
                   Is i a boundary node or interior node?


Constructors:
Files:
Example:

| #include "GoTools/parametrization/PrRectangularGrid_OP.h"
| #include "GoTools/parametrization/PrPlanarGraph_OP.h"
| 
| main()
| {
|   s_o << "Here is a rectangular grid:\n";
|   s_o << "\n";
|   s_o << "   01-------11-------21 \n";
|   s_o << "   |         |        | \n";
|   s_o << "   |         |        | \n";
|   s_o << "   |         |        | \n";
|   s_o << "   |         |        | \n";
|   s_o << "   |         |        | \n";
|   s_o << "   00-------10-------20 \n";
| 
|   double *grid;
|   int m=2, n=1;
|   grid = new double[3*(m+1)*(n+1)];
|   grid[0] = 0.0; grid[1] = 0.0; grid[2] = 0.0; 
|   grid[3] = 1.0; grid[4] = 0.0; grid[5] = 5.0; 
|   grid[6] = 2.0; grid[7] = 0.0; grid[8] = 10.0; 
|   grid[9] = 0.0; grid[10] = 1.0; grid[11] = 0.0; 
|   grid[12] = 1.0; grid[13] = 1.0; grid[14] = 5.0; 
|   grid[15] = 2.0; grid[16] = 1.0; grid[17] = 10.0; 
| 
|   PrRectangularGrid_OP rect_grid(m,n,grid);
|   delete [] grid;
| 
|   s_o << "The data of rect_grid is \n";
|   rect_grid.print(s_o);
|   s_o << "The information concerning rect_grid is \n";
|   rect_grid.printInfo(s_o);
| 
|   s_o << "The nodes of rect_grid are \n";
|   rect_grid.printXYZNodes(s_o);
|   s_o << "The edges of rect_grid are \n";
|   rect_grid.printXYZEdges(s_o);
|   s_o << "The faces of rect_grid are \n";
|   rect_grid.printXYZFaces(s_o);
| 
|   // make a planar graph from the rectangular grid
|   s_o << "The numbering of the nodes will be like this:\n";
|   s_o << "\n";
|   s_o << "   4---------5--------6 \n";
|   s_o << "   |         |        | \n";
|   s_o << "   |         |        | \n";
|   s_o << "   |         |        | \n";
|   s_o << "   |         |        | \n";
|   s_o << "   |         |        | \n";
|   s_o << "   1---------2--------3 \n";
| 
|   PrPlanarGraph_OP pl_graph(rect_grid);
|   s_o << "The data of pl_graph is \n";
|   pl_graph.print(s_o);
| 
|   return 0;
| }

See also:
Developed by:      SINTEF Applied Mathematics, Oslo, Norway
Author:	           Michael Floater, SINTEF
Date:              Apr. 97
*/

#endif // PRRECTANGULARGRID_OP_H
