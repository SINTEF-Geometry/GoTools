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









