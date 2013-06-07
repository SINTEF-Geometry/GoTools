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

#include "GoTools/utils/CompositeBox.h"

using namespace Go;
using std::streamsize;

namespace {
    void adjustPoints(Point& low, Point& high);

}; // End anonymous namespace

//===========================================================================
CompositeBox::~CompositeBox()
//===========================================================================
{
}

//===========================================================================
void CompositeBox::read(std::istream& is)
//===========================================================================
{
    inner_.read(is);
    edge_.read(is);
}


//===========================================================================
void CompositeBox::write(std::ostream& os) const
//===========================================================================
{
    streamsize prev = os.precision(15);
    inner_.write(os);
    edge_.write(os);
    os.precision(prev);   // Reset precision to it's previous value
}


/// The lower bound of the composite box.
Point CompositeBox::low(double toli, double tole) const
{
    Point l;
    if (obsolete_inner_)
    {
	l = edge_.low();
	for (int d = 0; d < dimension(); ++d) {
	    l[d] -= tole;
	}
    }
    else
    {
	l = inner_.low();
	for (int d = 0; d < dimension(); ++d) {
	    l[d] -= toli;
	    if (edge_.low()[d] - tole < l[d]) {
		l[d] = edge_.low()[d] - tole;
	    }
	}
    }
    return l;
}

/// The upper bound of the composite box.
Point CompositeBox::high(double toli, double tole) const
{
    Point h;
    if (obsolete_inner_)
    {
	h = edge_.high();
	for (int d = 0; d < dimension(); ++d) {
	    h[d] += tole;
	}
    }
    else
    {
	h = inner_.high();
	for (int d = 0; d < dimension(); ++d) {
	    h[d] += toli;
	    if (edge_.high()[d] + tole > h[d]) {
		h[d] = edge_.high()[d] + tole;
	    }
	}
    }
    return h;
}

//===========================================================================
bool CompositeBox::containsPoint(const Point& pt,
				 double toli,
				 double tole) const
//===========================================================================
{
    Point ll = low(toli, tole);
    Point hh = high(toli, tole);

    // Check that the size of the box is not less than zero
    adjustPoints(ll, hh);
    BoundingBox temp(ll, hh);
    return temp.containsPoint(pt);
}


//===========================================================================
bool CompositeBox::overlaps(const CompositeBox& box, 
			    double toli,
			    double tole) const
//===========================================================================
{
    double overlap;
    return getOverlap(box, overlap, toli, tole);
}

//===========================================================================
bool CompositeBox::getOverlap(const CompositeBox& box, 
			    double& overlap,
			    double toli,
			    double tole) const
//===========================================================================
{
    Point low1 = low(toli, tole);
    Point high1 = high(toli, tole);
    adjustPoints(low1, high1);
    BoundingBox temp1(low1, high1);
    Point low2 = box.low(toli, tole);
    Point high2 = box.high(toli, tole);
    adjustPoints(low2, high2);
    BoundingBox temp2(low2, high2);
    return temp1.getOverlap(temp2, overlap, 0.0);
}



//===========================================================================
bool CompositeBox::containsBox(const CompositeBox& box, 
			       double toli,
			       double tole) const
//===========================================================================
{
    Point low1 = low(toli, tole);
    Point high1 = high(toli, tole);
    adjustPoints(low1, high1);
    BoundingBox temp1(low1, high1);
    Point low2 = box.low(toli, tole);
    Point high2 = box.high(toli, tole);
    adjustPoints(low2, high2);
    BoundingBox temp2(low2, high2);
    return temp1.containsBox(temp2);
}

namespace {
//===========================================================================
void adjustPoints(Point& low, Point& high)
{
    // Check that the size of the box is not less than zero
    int ki;
    int dim = low.dimension();
    for (ki=0; ki<dim; ki++)
	if (high[ki] < low[ki])
	{
	    double dum = low[ki];
	    low.resetValue(ki, high[ki]);
	    high.resetValue(ki, dum);
	}
}

};
