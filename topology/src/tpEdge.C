//===========================================================================
//                                                                           
// File: tpEdge.C                                                           
//                                                                           
// Created: Thu Jul 11 12:42:53 2002                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@math.sintef.no>
//                                                                           
// Revision: $Id: tpEdge.C,v 1.21 2009-01-30 09:53:14 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/topology/tpEdge.h"
#include "GoTools/topology/tpUtils.h"


//===========================================================================
tpEdge::tpEdge()
    : next_(0), prev_(0), twin_(0)
//=========================================================================== 
{
}

//===========================================================================
tpEdge::~tpEdge()
//===========================================================================
{
}

//===========================================================================
tpEdge* tpEdge::next()
//===========================================================================
{
    return next_;
}

//===========================================================================
tpEdge* tpEdge::prev()
//===========================================================================
{
    return prev_;
}

//===========================================================================
tpEdge* tpEdge::twin()
//===========================================================================
{
    return twin_;
}

//===========================================================================
void tpEdge::connectAfter(tpEdge* edge)
//===========================================================================
{
    if (edge == 0)
	return;
    next_ = edge->next_;
    prev_ = edge;
    edge->next_ = this;
    if (next_)
	next_->prev_ = this;
}

//===========================================================================
void tpEdge::closeLoop(tpEdge* last)
//===========================================================================
{
    ALWAYS_ERROR_IF(prev_ || last->next_, "Could not close loop");

    prev_ = last;
    last->next_ = this;
}

//===========================================================================
void tpEdge::disconnectThis()
//===========================================================================
{
    if (prev_)
	prev_->next_ = next_;
    if (next_)
	next_->prev_ = prev_;
    next_ = 0;
    prev_ = 0;
}

//===========================================================================
void tpEdge::connectTwin(tpEdge* newtwin, int& status)
//===========================================================================
{
    ALWAYS_ERROR_IF((twin_ || newtwin->twin_) &&
		twin_ != newtwin, // Allow existing twins to connect.
		    "Edge already has a twin.");

    if (!(twin_ && twin_ != newtwin && status == 2))
    {
	twin_ = newtwin;
	newtwin->twin_ = this;
	status = 3;
    }
}

//===========================================================================
void tpEdge::disconnectTwin()
//===========================================================================
{
    ALWAYS_ERROR_IF(twin_ == 0, "Edge has no twin.");
    twin_->twin_ = 0;
    twin_ = 0;
}

//===========================================================================
void tpEdge::adjacentEdges(bool at_start_of_edge,
			    std::vector<tpEdge*>& adjacent,
			    std::vector<bool>& at_start)
//===========================================================================
{
    tpUtils::adjacentEdges<tpEdge>(this, at_start_of_edge, adjacent, at_start);
}

//===========================================================================
tpJointType tpEdge::checkContinuity(tpEdge* nextedge, double neighbour,
					double gap, double bend, double kink) const
//===========================================================================
{
    return tpUtils::checkContinuity<tpEdge>(this, nextedge,
					     neighbour, gap, bend, kink);
}
