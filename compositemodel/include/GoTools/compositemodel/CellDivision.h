//===========================================================================
//                                                                           
// File: CellDivision.h                                                    
//                                                                           
// Created: November 2007
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================

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

