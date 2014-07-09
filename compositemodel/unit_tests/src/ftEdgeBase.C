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

#include "GoTools/compositemodel/ftEdgeBase.h"
#include "GoTools/topology/tpUtils.h"
#include "GoTools/utils/errormacros.h"
#include <fstream>

namespace Go
{



//===========================================================================
ftEdgeBase::ftEdgeBase()
    : next_(0), prev_(0), twin_(0)
//=========================================================================== 
{
}


//===========================================================================
ftEdgeBase::~ftEdgeBase()
//===========================================================================
{
}

//===========================================================================
ftEdgeBase* ftEdgeBase::next()
//===========================================================================
{
    return next_;
}

//===========================================================================
ftEdgeBase* ftEdgeBase::prev()
//===========================================================================
{
    return prev_;
}

//===========================================================================
ftEdgeBase* ftEdgeBase::twin()
//===========================================================================
{
    return twin_;
}

//===========================================================================
void ftEdgeBase::connectAfter(ftEdgeBase* edge)
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
void ftEdgeBase::closeLoop(ftEdgeBase* last)
//===========================================================================
{
    ALWAYS_ERROR_IF(prev_ || last->next_, "Could not close loop");
    prev_ = last;
    last->next_ = this;
}

//===========================================================================
void ftEdgeBase::disconnectThis()
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
void ftEdgeBase::connectTwin(ftEdgeBase* newtwin, int& status)
//===========================================================================
{
    /*ALWAYS_ERROR_IF((twin_ || newtwin->twin_) &&
		    twin_ != newtwin, // Allow existing twins to connect.		    
		    "Edge already has a twin.");*/
    if ((twin_ || newtwin->twin_) && twin_ != newtwin)
    {
// 	std::ofstream debug("twin.g2");
// 	ftFaceBase* f1 = geomEdge()->face();
// 	ftFaceBase* f2 = (twin_) ? twin_->geomEdge()->face() : 
// 	  newtwin->twin_->geomEdge()->face();
// 	ftFaceBase* f3 = newtwin->geomEdge()->face();
// 	ParamSurface *s1 = f1->surface().get();
// 	ParamSurface *s2 = f2->surface().get();
// 	ParamSurface *s3 = f3->surface().get();
// 	s1->writeStandardHeader(debug);
// 	s1->write(debug);
// 	s2->writeStandardHeader(debug);
// 	s2->write(debug);
// 	s3->writeStandardHeader(debug);
// 	s3->write(debug);
	MESSAGE("Edge already has a twin.");
    }
    if (!(twin_ && twin_ != newtwin && status == 2))
    {
	twin_ = newtwin;
	newtwin->twin_ = this;
	status = 3;
    }
}

//===========================================================================
void ftEdgeBase::disconnectTwin()
//===========================================================================
{
  //ALWAYS_ERROR_IF(twin_ == 0, "Edge has no twin.");
  if (twin_)
    twin_->twin_ = 0;
  twin_ = 0;
}

//===========================================================================
void ftEdgeBase::adjacentEdges(bool at_start_of_edge,
			       std::vector<ftEdgeBase*>& adjacent,
			       std::vector<bool>& at_start)
//===========================================================================
{
    tpUtils::adjacentEdges<ftEdgeBase>(this, at_start_of_edge, adjacent, at_start);
}

//===========================================================================
tpJointType ftEdgeBase::checkContinuity(ftEdgeBase* nextedge, double neighbour,
					double gap, double bend, double kink) const
//===========================================================================
{
    return tpUtils::checkContinuity<ftEdgeBase>(this, nextedge,
						neighbour, gap, bend, kink);
}

//===========================================================================
bool ftEdgeBase::orientationOK() const
{
    // No possibility to check at this level. Return OK.
    return true;
}

//===========================================================================
  bool ftEdgeBase::checkOverlap(ftEdgeBase *other, double tol, int nmbsample,
			      double& t1, double& t2, double& t3, 
				double& t4, bool& same_dir, bool no_snap) const
//===========================================================================
{
  // Check endpoints
  int ki;
  Point endpos[2][2];
  t1 = tMin();
  t2 = tMax();
  t3 = other->tMin();
  t4 = other->tMax();
  endpos[0][0] = point(t1);
  endpos[0][1] = point(t2);
  endpos[1][0] = other->point(t3);
  endpos[1][1] = other->point(t4);
  double fac = 1.0e-4;


  // Initial guess
  same_dir = true;

  // Check configuration
  double par, dist;
  Point pos;
  if (endpos[0][0].dist(endpos[1][0]) <= tol &&
      endpos[0][1].dist(endpos[1][1]) <= tol)
    {
      same_dir = true;
    }
  else if (endpos[0][0].dist(endpos[1][1]) <= tol &&
      endpos[0][1].dist(endpos[1][0]) <= tol)
    {
      same_dir = false;
    }
  else 
    {
      if (endpos[0][0].dist(endpos[1][0]) > tol &&
	  endpos[0][0].dist(endpos[1][1]) > tol)
	{
	  other->closestPoint(endpos[0][0], par, pos, dist);
	  if (dist <= tol)
	    {
	      if (endpos[1][1].dist(endpos[0][1]) <= tol)
		{
		  if (no_snap || par-t3 > fac*(t4 - t3))
		    t3 = par;
		  endpos[1][0] = pos;
		}
	      else
		{
		  if (no_snap || t4-par > fac*(t4-t3))
		    t4 = par;
		  endpos[1][1] = pos;
		  same_dir = false;
		}
	    }
	}
      if (endpos[0][1].dist(endpos[1][0]) > tol &&
	  endpos[0][1].dist(endpos[1][1]) > tol)
	{
	  other->closestPoint(endpos[0][1], par, pos, dist);
	  if (dist <= tol)
	    {
	      if (endpos[1][0].dist(endpos[0][0]) <= tol)
		{
		  if (no_snap || t4-par > fac*(t4-t3))
		    t4 = par;
		  endpos[1][1] = pos;
		}
	      else
		{
		  if (no_snap || par-t3 > fac*(t4-t3))
		    t3 = par;
		  endpos[1][0] = pos;
		  same_dir = false;
		}
	    }
	}
      if (endpos[1][0].dist(endpos[0][0]) > tol &&
	  endpos[1][0].dist(endpos[0][1]) > tol)
	{
	  closestPoint(endpos[1][0], par, pos, dist);
	  if (dist <= tol)
	    {
	      if (endpos[0][1].dist(endpos[1][1]) <= tol)
		{
		  if (no_snap || par-t1 > fac*(t2-t1))
		    t1 = par;
		  endpos[0][0] = pos;
		}
	      else
		{
		  if (no_snap || t2-par > fac*(t2-t1))
		    t2 = par;
		  endpos[0][1] = pos;
		  same_dir = false;
		}
	    }
	}
      if (endpos[1][1].dist(endpos[0][0]) > tol &&
	  endpos[1][1].dist(endpos[0][1]) > tol)
	{
	  closestPoint(endpos[1][1], par, pos, dist);
	  if (dist <= tol)
	    {
	      if (endpos[0][0].dist(endpos[1][0]) <= tol)
		{
		  if (no_snap || t2-par > fac*(t2-t1))
		    t2 = par;
		  endpos[0][1] = pos;
		}
	      else
		{
		  if (no_snap || par-t1 > fac*(t2-t1))
		    t1 = par;
		  endpos[0][0] = pos;
		  same_dir = false;
		}
	    }
	}

      if (t1 > t2 && t3 > t4)
	{
	  std::swap(t1, t2);
	  std::swap(endpos[0][0], endpos[0][1]);
	  std::swap(t3, t4);
	  std::swap(endpos[1][0], endpos[1][1]);
	}
      else if (t1 > t2)
	{
	  std::swap(t1, t2);
	  std::swap(endpos[0][0], endpos[0][1]);
	  same_dir = false;
	}
      else if (t3 > t4)
	{
	  std::swap(t3, t4);
	  std::swap(endpos[1][0], endpos[1][1]);
	  same_dir = false;
	}
    }

  // Check endpoint coincidence
  for (ki=0; ki<2; ki++)
    {
      if (same_dir && (endpos[0][0].dist(endpos[1][0]) > tol ||
		       endpos[0][1].dist(endpos[1][1]) > tol))
	return false;
      else if (!same_dir && (endpos[0][0].dist(endpos[1][1]) > tol ||
		       endpos[0][1].dist(endpos[1][0]) > tol))
	return false;
    }

  if (nmbsample <= 2)
    return true;  // Only endpoint check expected

  // Check internally
  double tint1 = (t2 - t1)/(double)(nmbsample-1);
  double tint2 = (t4 - t3)/(double)(nmbsample-1);
  int sgn = (same_dir) ? 1 : -1;
  double guess, par2;
  Point pos2;

  double tol2 = 2.0*tol;  // More loose checking in the inner
  guess = (same_dir) ? t3 : t4;
  for (ki=1, par=t1+tint1; ki<nmbsample-1; ki++, par+=tint1)
    {
      pos = point(par);
      guess += sgn*tint2;
      other->closestPoint(pos, par2, pos2, dist, &guess);
      if (dist > tol2)
	return false;  // No coincidence
    }

  return true;  // Coincidence found
}

} // namespace Go
