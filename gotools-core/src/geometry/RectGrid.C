//===========================================================================
//                                                                           
// File: RectGrid.C                                                          
//                                                                           
// Created: Mon Jan 17 14:38:40 2005                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: RectGrid.C,v 1.3 2007-03-08 10:43:29 afr Exp $
//                                                                           
//===========================================================================

#include "GoTools/geometry/RectGrid.h"
#include "GoTools/geometry/SplineUtils.h"

using namespace std;

namespace Go
{

//===========================================================================
RectGrid::~RectGrid()

//===========================================================================
{
}

//===========================================================================
void RectGrid::read (std::istream& is)
//===========================================================================
{
    // Canonical data
    is >> dim_;
    is >> numu_;
    is >> numv_;
    int n = numu_*numv_*dim_;
    points_.resize(n);
    for (int i = 0; i < n; ++i) {
	is >> points_[i];
    }
}

//===========================================================================
void RectGrid::write (std::ostream& os) const
//===========================================================================
{
    streamsize prev = os.precision(15);

    os << dim_ << ' ' << numu_ << ' ' << numv_ << '\n';
    int n = numu_*numv_;
    for (int i = 0; i < n; ++i) {
	os << points_[i*dim_];
	for (int d = 1; d < dim_; ++d) {
	    os << ' ' << points_[i*dim_ + d];
	}
	os << '\n';
    }
    os << endl;
    os.precision(prev);   // Reset precision to it's previous value
}

//===========================================================================
BoundingBox RectGrid::boundingBox() const
//===========================================================================
{
    BoundingBox box;
    box.setFromArray(&points_[0], &points_[0] + points_.size(), dim_);
    return box;
}

//===========================================================================
int RectGrid::dimension() const
//===========================================================================
{
    return dim_;
}

//===========================================================================
ClassType RectGrid::instanceType() const
//===========================================================================
{
    return classType();
}

//===========================================================================
void RectGrid::swapDirections() 
//===========================================================================
{
    SplineUtils::transpose_array(dim_, numv_, numu_, &points_[0]);
    swap(numu_, numv_);
}

} // namespace Go
