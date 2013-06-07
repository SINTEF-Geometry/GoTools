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

#ifndef _CELLDIVISION_H
#define _CELLDIVISION_H

#include "GoTools/compositemodel/ftPoint.h"
#include "GoTools/utils/BoundingBox.h"
#include <vector>             // Standard library STL vector

namespace Go
{

class ftSurface;

// ------------- Some helper classes -------------------

/// Helper class handling one of the cells in class CellDivision.
class ftCell
{
protected:
    std::vector<ftSurface*> faces_;
    std::vector<BoundingBox> face_boxes_;
    BoundingBox box_;
public:
    ftCell() {}
    ~ftCell() {}
    void setBox(const BoundingBox& box)
        { box_ = box; }
    void addFace(ftSurface* f)
        { faces_.push_back(f); }
    void addFaceBox(BoundingBox& box)
        { face_boxes_.push_back(box); }
    const BoundingBox& box() const
        { return box_; }
    int num_faces() const
    { return (int)faces_.size(); }
    ftSurface* face(int i) const
        { return faces_[i]; }
    const BoundingBox& faceBox(int i) const
        { return face_boxes_[i];}

    void removeFace(ftSurface* face)
        {
            for (size_t ki=0; ki<faces_.size(); ki++)
                if (faces_[ki] == face)
                    faces_.erase(faces_.begin()+ki);
        }


};

/// Function object for sorting ftCells.
class ftCellInfo
{
public:
    int index_;
    double dist_;
    ftCellInfo() :  index_(-1), dist_(-1) {}
    ftCellInfo(int i, double d) :  index_(i), dist_(d) {}
    bool operator < (const ftCellInfo& ci) const { return dist_ < ci.dist_; }
};


/// Used internally in SurfaceModel. CellDivision divides the space occupied 
/// by a surface model into a number of equal boxes. 
/// The division is intended to improve the performance in closest
/// point and intersection computations by giving the means to perform an initial 
/// interception test before the actual computations are started
class CellDivision
{
 public:
    CellDivision()
        //: ncellsx_(0),ncellsy_(0), ncellsz_(0)
        : dim_(0)
        {
                ncells_.clear();
        }

    CellDivision(std::vector<ftSurface*> faces, int n1, int n2);
    CellDivision(std::vector<ftSurface*> faces, int n1, int n2, int n3);
    CellDivision(std::vector<ftSurface*> faces, std::vector<int> n);
    ~CellDivision();

    int numCells() const
        {
            return (int)cells_.size();
        }

    const ftCell& getCell(int idx) const
        {
            return cells_[idx];
        }

    void addFace(ftSurface* face);
    void removeFace(ftSurface* face);
    void constructCells();
    void constructCells(int n1, int n2);
    void constructCells(int n1, int n2, int n3);
    void constructCells(std::vector<int> n);

    BoundingBox big_box() const
      {
        return big_box_;
      }


 private:
    std::vector<ftSurface*> faces_;
    BoundingBox big_box_;
    //int ncellsx_, ncellsy_, ncellsz_;
    int dim_;
    std::vector<int> ncells_;
    Point cell_delta_;
    std::vector<ftCell> cells_;

    void makeBox();
    //void cellContaining(const Point& pt, int& i1, int& i2, int& i3);
    void cellContaining(const Point& pt, std::vector<int>& ii);

    std::vector<std::vector<int> > coords_;
    void makeCoords();
    void updateCoords(int curr_dim, std::vector<int>& curr_coord);
    
};

} // namespace Go

#endif // _CELLDIVISION_H

