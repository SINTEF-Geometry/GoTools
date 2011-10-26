//===========================================================================
//                                                                           
// File: CompositeBox.C                                                      
//                                                                           
// Created: Wed Dec  8 12:04:48 2004                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: CompositeBox.C,v 1.5 2005-02-23 07:48:41 vsk Exp $
//                                                                           
//===========================================================================


#include "GoTools/utils/CompositeBox.h"

using namespace Go;

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
    inner_.write(os);
    edge_.write(os);
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
