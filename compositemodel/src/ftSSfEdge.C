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

#include "GoTools/compositemodel/ftEdge.h"
#include "GoTools/compositemodel/ftSSfEdge.h"
//#include "GoTools/compositemodel/ftSurface.h"
//#include "GoTools/model_toolbox/ftSuperSurface.h"
#include "GoTools/geometry/SplineCurve.h"

using namespace Go;

//===========================================================================
ftSSfEdge::ftSSfEdge(ftFaceBase* face, ftEdge *edge, int entry_id)
//===========================================================================
    : face_(face), edg_(edge),
      entry_id_(entry_id)
{
    
}

//===========================================================================
ftSSfEdge::~ftSSfEdge()
//===========================================================================
{
  //std::cout << "Deletes edge" << this << std::endl;
}


// //===========================================================================
// void ftSSfEdge::turnOrientation()
// //===========================================================================
// {
//     edg_->turnOrientation();
// }

// //===========================================================================
// void ftSSfEdge::setOrientation()
//     //===========================================================================
// {
//     edg_->setOrientation();
// }

// //===========================================================================
// bool ftSSfEdge::isTurned()
//     //===========================================================================
// {
//     return edg_->isTurned();
// }

//===========================================================================
ftFaceBase* ftSSfEdge::face()
    //===========================================================================
{
    return face_;
}


//===========================================================================
BoundingBox ftSSfEdge::boundingBox()
    //===========================================================================
{
    return edg_->boundingBox();
}

//===========================================================================
void ftSSfEdge::closestPoint(const Point& pt,
			     double& clo_t,
			     Point& clo_pt,
			     double& clo_dist,
			     double const *seed) const
    //===========================================================================
{
    edg_->closestPoint(pt, clo_t, clo_pt, clo_dist, seed);
}


//===========================================================================
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
ftEdgeBase* ftSSfEdge::split(double t)
#else
ftSSfEdge* ftSSfEdge::split(double t)
#endif
    //===========================================================================
{
    ALWAYS_ERROR_IF(twin(), "Cannot split edge, already has twin!");
    ALWAYS_ERROR_IF(!prev() || !next(), 
		"Cannot split edge, not fully connected");

#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
    ftEdge *newedge = dynamic_cast<ftEdge*>(edg_->split(t));
    ALWAYS_ERROR_IF(newedge==NULL,
		    "Bad cast");

#else
    ftEdge *newedge = edg_->split(t);
#endif

    ftSSfEdge* newssfedge;
    newssfedge = new ftSSfEdge(face_, newedge);
    newssfedge->connectAfter(this);
    return newssfedge;
}


//===========================================================================
Point ftSSfEdge::point(double t) const
    //===========================================================================
{
    return edg_->point(t);
}


//===========================================================================
Point ftSSfEdge::tangent(double t) const
    //===========================================================================
{
    return edg_->tangent(t);
}


//===========================================================================
Point ftSSfEdge::normal(double t) const
    //===========================================================================
{
    return edg_->normal(t);
}

//===========================================================================
Point ftSSfEdge::normal(double t, Point& face_par_pt, double* face_seed) const
    //===========================================================================
{
    Point normal_pt = edg_->normal(t, face_par_pt, face_seed);

    return normal_pt;
}

//===========================================================================
ftEdge* ftSSfEdge::geomEdge()
    //===========================================================================
{
    return edg_; 
}

