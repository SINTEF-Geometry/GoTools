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

#include "GoTools/topology/tpEdge.h"
#include "GoTools/topology/tpUtils.h"

using namespace Go;

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
