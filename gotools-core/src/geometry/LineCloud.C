//===========================================================================
//                                                                           
// File: LineCloud.C                                                       
//                                                                           
// Created: Tue Oct 29 09:37:45 2002                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: LineCloud.C,v 1.8 2005-09-19 09:24:42 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/geometry/LineCloud.h"


using namespace std;


namespace Go {


//===========================================================================
LineCloud::~LineCloud()
//===========================================================================
{
}


//===========================================================================
void LineCloud::read(std::istream& is)
//===========================================================================
{
    bool is_good = is.good();
    if (!is_good) {
	THROW("Invalid geometry file!");
    }
    int numl;
    is >> numl;
    ALWAYS_ERROR_IF(numl < 1, "Less than one line in cloud.");
    is_good = is.good();
    if (!is_good) {
	THROW("Invalid geometry file!");
    }
    points_.resize(numl*2);
    for (int i = 0; i < numl; ++i) {
	is >> points_[2*i][0];
	is >> points_[2*i][1];
	is >> points_[2*i][2];
	is >> points_[2*i+1][0];
	is >> points_[2*i+1][1];
	is >> points_[2*i+1][2];
    }

    is_good = is.good();
    if (!is_good) {
	THROW("Invalid geometry file!");
    }
}

//===========================================================================
void LineCloud::write(std::ostream& os) const
//===========================================================================
{
    streamsize prev = os.precision(15);

    int numl = (int)points_.size()/2;
    os << numl << '\n';
    for(int i = 0; i < numl; ++i) {
	os << points_[2*i][0] << ' ';
	os << points_[2*i][1] << ' ';
	os << points_[2*i][2] << "     ";
	os << points_[2*i+1][0] << ' ';
	os << points_[2*i+1][1] << ' ';
	os << points_[2*i+1][2] << '\n';
    }
    os << endl;
    os.precision(prev);   // Reset precision to it's previous value
}



//===========================================================================
BoundingBox LineCloud::boundingBox() const
//===========================================================================
{
    BoundingBox box;
    box.setFromArray(points_[0].begin(), points_.back().end(), 3);
    return box;
}

//===========================================================================
void LineCloud::setCloud(const double* points, int numlines)
//===========================================================================
{
    DEBUG_ERROR_IF(numlines < 1, "Less than one line in cloud.");
    points_.resize(numlines*2);
    std::copy(points, points + 6*numlines, points_[0].begin());
}


} // namespace Go









