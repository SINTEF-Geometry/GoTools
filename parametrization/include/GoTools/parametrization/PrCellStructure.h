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

#ifndef PRCELLSTRUCTURE_H
#define PRCELLSTRUCTURE_H

#include "GoTools/utils/Array.h"
using Go::Vector3D;
#include <vector>
using std::vector;
using std::istream;
using std::ostream;

/*<PrCellStructure-syntax: */

/** PrCellStructure - Represents a set of points in three dimensions
 * and a cell structure: a 3D grid of cubical cells each with a local
 * list of those points which lie inside it.
 * This allows one to perform various queries, such as finding all points
 * in a ball, very efficiently.
 * Creating the cell structure requires only O(N) operations,
 * where \a N is the number of points.
 */
class PrCellStructure
{
private:
    vector<Vector3D> xyz_;
    //VecSimplest(CgPoint3d) xyz_;
    int n_; // no points,
    double min_[3]; // min value in each coordinate
    int ncells_[3]; // number of cells in each coordinate
    double cell_size_; // dimension of cell, each cell is a cube
    vector<vector<int> > ind_;
    //VecSimplest(PrIntList) ind_;
    int max_no_cells_;

    void makeCellStructure();

public:
    /// Default constructor
    PrCellStructure() {}
    /// Constructor
    /// \param n total number of points
    /// \param xyz_points pointer to an array of points stored xyz-wise (The 
    ///                   points will be copied to internal data structure.
    /// \param max_no_cells indicate maximum number of cells along the 
    ///                     X, Y and Z axis.
    PrCellStructure(int n, double* xyz_points, int max_no_cells = 10);

    /// Destructor
    virtual ~PrCellStructure() {};

    /// Reset the PrCellStructure to a set of new points, deleting old content.
    /// \param n number of points
    /// \param xyz_points pointer to the array of stored points.  (The points
    ///                   will be copied to internal data structure).
    void      attach(int n, const double* xyz_points);

    /// Set the maximum number of cells in X, Y and Z direction to 'max_no_cells'.
    /// This must be set before the cells are generated (ie. before the 'attach()'
    /// command.
    void      setNumCells(int max_no_cells) {max_no_cells_ = max_no_cells; }

    /// Get the previously set limit for the maximum number of cells in the 
    /// X, Y and Z direction.
    int       getNumCells() const {return max_no_cells_; }

    /// Get the number of points (nodes) stored in the cell structure.
    int       getNumNodes() const {return (int)xyz_.size(); }

    /// Get a specific point (node) in the PrCellStructure by its index.
    Vector3D get3dNode(int i) const {return xyz_[i]; }

    /// Change the coordinates of a specific point in the PrCellStructure.
    /// \param i index of node to change
    /// \param p the new coordinates
    void       set3dNode(int i, const Vector3D& p) 
    {xyz_[i].x() = p.x(); xyz_[i].y() = p.y(); xyz_[i].z() = p.z();}

    /// Return all points within the ball of radius radius around
    /// the point p. Don't include p itself if notP = 1.
    /// \param p the center of the ball
    /// \param radius the radius of the ball
    /// \retval neighbours will be filled with the indexes of the 
    ///                    neighbour points.
    /// \param notP if != 0, 'p' itself will not be included in
    ///             the list of returned points.
    void getBall(const Vector3D& p, double radius,
		 vector<int>& neighbours, int notP = 0) const;

    /// Return k nearest points to the point p.
    /// Don't include p itself if notP = 1.
    /// \param p the point for which we seek the k nearest neighbours
    /// \param k the number of neighbours we seek
    /// \retval neighbours will be filled with the indexes of the 
    ///                    neighbour points.
    /// \param notP if != 0, 'p' itself will not be included in the 
    ///             list of returned points.
    void getKNearest(const Vector3D& p, int k,
		     vector<int>& neighbours, int notP = 0) const;

    /// Get the index of the cell positioned at 'i', 'j', 'k'.
    int getI(int i, int j, int k) const
    {return i + ncells_[0] * (j + ncells_[1] * k); }

    /// Get the 'i', 'j', 'k'-position of the cell indexed 'ii'.
    void getIJK(int ii, int& i, int& j, int& k) const;

    /// Locate the cell which contains the coordinate 'xyz'
    /// \param xyz we want to find the cell containing these coordinates
    /// \retval i i-position of the located cell
    /// \retval j j-position of the located cell
    /// \retval k k-position of the located cell
    void whichCell(const Vector3D& xyz, int& i, int& j, int& k) const;

    //print and scan routines
    /// write points to stream
    void print(ostream& os);

    /// read points from stream and regenerate cell structure.
    void scan(istream& is);
};

/*>PrCellStructure-syntax: */

/*Class:PrCellStructure

Name:              PrCellStructure
Syntax:	           @PrCellStructure-syntax
Keywords:
Description:       This class represents a set of points in three dimensions
                   and a cell structure: a 3D grid of cubical cells
                   each with a local list of those points which lie inside it.
                   This allows one to perform various queries, such as
                   finding all points in a ball, very efficiently.
                   Creating the cell structure requires only O(N) operations,
                   where $N$ is the number of points.
Member functions:

Constructors:
Files:
Example:


See also:
Developed by:      SINTEF Applied Mathematics, Oslo, Norway
Author:	           Michael Floater, SINTEF
Date:              Sep. 99
*/

#endif // PRCELLSTRUCTURE_H
