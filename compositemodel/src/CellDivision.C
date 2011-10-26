//===========================================================================
//                                                                           
// File: CellDivision.C                                                   
//                                                                           
// Created: Nov 2007
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/compositemodel/CellDivision.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/utils/BoundingBox.h"

using namespace std;

namespace Go
{



//===========================================================================
CellDivision::CellDivision(vector<ftSurface*> faces, int n1, int n2, int n3)
//===========================================================================
    : ncellsx_(n1), ncellsy_(n2), ncellsz_(n3)
{
    faces_ = faces;
    makeBox();
    constructCells();
}

//===========================================================================
CellDivision::~CellDivision()
//===========================================================================
{
}

//===========================================================================
void CellDivision::addFace(ftSurface* face)
//===========================================================================
{
    size_t ki;
    for (ki=0; ki<faces_.size(); ki++)
	if (faces_[ki] == face)
	    break;

    if (ki < faces_.size())
	return;  // Face already in list

    faces_.push_back(face);
    
    BoundingBox box = face->boundingBox();
    if (!big_box_.containsBox(box))
    {
	big_box_.addUnionWith(box);
	constructCells();
    }
}

//===========================================================================
void CellDivision::removeFace(ftSurface* face)
//===========================================================================
{
    // Note that this function does not reset the size of the box and the cells
    size_t ki;
    for (ki=0; ki<faces_.size(); ki++)
	if (faces_[ki] == face)
	    faces_.erase(faces_.begin()+ki);

    BoundingBox box = face->boundingBox();
    for (ki=0; ki<cells_.size(); ki++)
    {
	if (box.overlaps(cells_[ki].box()))
	    cells_[ki].removeFace(face);
    }
}

//===========================================================================
void CellDivision::makeBox()
//===========================================================================
{
    // First, we make the big box:

    // We check to see if there are any faces. @jbt
    if (faces_.empty()) {
	MESSAGE("No faces - returning empty BoundingBox.");
	big_box_ = BoundingBox();
	return;
    }

    big_box_ = faces_[0]->boundingBox();
    for (int i=1; i<(int)faces_.size(); ++i) {
	BoundingBox curr_box = faces_[i]->boundingBox();
	big_box_.addUnionWith(curr_box);
    }

}

//===========================================================================
void CellDivision::constructCells()
//===========================================================================
{
    // @@@ NOTE
    //
    // This implementation depends upon using bounding boxes that are
    // aligned with the principal axis and box-shaped.

    // First we check if big_box_ is valid. @jbt
    if (!big_box_.valid()) {
	MESSAGE("Big box is invalid - returning empty cell vector.");
	cells_.clear();
	return;
    }

    int i, j, k;

    // We partition space into cells.
    cells_.resize(ncellsx_*ncellsy_*ncellsz_);
    Point low = big_box_.low();
    Point len = big_box_.high() - low;
    cell_delta_ = Point(len[0]/ncellsx_, len[1]/ncellsy_, len[2]/ncellsz_);
    const Point& delta = cell_delta_;

    double x = low[0];
    double y = low[1];
    double z = low[2];
    Point corner;

    for (i=0; i<ncellsx_; ++i) {
	y = low[1];
	for (j=0; j<ncellsy_; ++j) {
	    z = low[2];
	    for (k=0; k<ncellsz_; ++k) {
		corner.setValue(x, y, z);
		cells_[i + j*ncellsx_ + k*ncellsx_*ncellsy_]
		    .setBox(BoundingBox(corner, corner + delta));
		z += delta[2];
	    }
	    y += delta[1];
	}
	x += delta[0];
    }

    int x1, x2, y1, y2, z1, z2;
    Point lo, hi;
    // Now that we have the cells, we check every face, and add its pointer
    // to every cell overlapping its bounding box.
    for (size_t f=0; f<faces_.size(); ++f) {
	BoundingBox b = faces_[f]->boundingBox();
	lo = b.low();
	hi = b.high();
	cellContaining(lo, x1, y1, z1);
	cellContaining(hi, x2, y2, z2);
	for (i=x1; i<=x2; ++i)
	    for (j=y1; j<=y2; ++j)
		for (k=z1; k<=z2; ++k)
		{
		    cells_[i + j*ncellsx_ + k*ncellsx_*ncellsy_].addFace(faces_[f]);
		    cells_[i + j*ncellsx_ + k*ncellsx_*ncellsy_].addFaceBox(b);
		}
    }

#ifdef DEBUG_CELLS
    for (i=0; i<ncellsx_*ncellsy_*ncellsz_; ++i) {
	cout << "Cell " << i << " has "
	     << cells_[i].num_faces() << " faces ";
	for (j=0; j<cells_[i].num_faces(); ++j)
	    cout << cells_[i].face(j)->GetId() << " ";
	cout << endl;
    }
#endif

}



//===========================================================================
void CellDivision::constructCells(int n0, int n1, int n2)
//===========================================================================
// Like constructCells(), but resets number of cells
{
  ncellsx_ = n0;
  ncellsy_ = n1;
  ncellsz_ = n2;
  constructCells();
}



//===========================================================================
void CellDivision::cellContaining(const Point& pt, 
				  int& i1, int& i2, int& i3)
//===========================================================================
{
    Point low = big_box_.low();
    double d[3];
    int ii[3];
    int nc[3];
    nc[0] = ncellsx_;
    nc[1] = ncellsy_;
    nc[2] = ncellsz_;
    // @@ Hardcoded constant
    double eps = 1e-13;
    for (int i = 0; i < 3; ++i) {
	d[i] = (pt[i]-low[i])/cell_delta_[i];
	ii[i] = static_cast<int>(floor(d[i]));
	ii[i] = max(min(ii[i], nc[i]-1),0);
	MESSAGE_IF (d[i] < -eps,
		       "Cell index below 0");
	MESSAGE_IF (d[i] > nc[i] + eps,
		       "Cell index above ncells[xyz]");
    }
    i1 = ii[0];
    i2 = ii[1];
    i3 = ii[2];
}


} // namespace Go


