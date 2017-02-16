//===========================================================================
//                                                                           
// File: ftSSfEdge.C                                                            
//                                                                           
// Created: 010802
//                                                                           
// Author: Vibeke Skytt, SINTEF
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================


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

