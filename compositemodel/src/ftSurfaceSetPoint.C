//===========================================================================
//                                                                           
// File: ftSurfaceSetPoint.C                                                 
//                                                                           
// Created: Mon Feb  4 00:27:47 2002                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@math.sintef.no>
//                                                                           
// Revision: $Id: ftSurfaceSetPoint.C,v 1.1 2009-01-23 13:34:38 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/compositemodel/ftSurfaceSetPoint.h"

using std::vector;

namespace Go
{


//===========================================================================
ftSurfaceSetPoint::ftSurfaceSetPoint(Vector3D xyz, int bnd)
    : ftSamplePoint(xyz, bnd)
//===========================================================================
{
}

//===========================================================================
ftSurfaceSetPoint::ftSurfaceSetPoint(Vector3D xyz, int bnd,
				     shared_ptr<ftFaceBase>& face, Vector2D par_pt)
    : ftSamplePoint(xyz, bnd)
//===========================================================================
{
    addPair(face, par_pt); // In theory we could avoid pushing back inner points.
}

//===========================================================================
ftSurfaceSetPoint::~ftSurfaceSetPoint()
//===========================================================================
{
}

//===========================================================================
bool ftSurfaceSetPoint::containsFace(ftFaceBase* face) const
//===========================================================================
{
    for (size_t ki=0; ki<par_pts_.size(); ++ki)
	if (par_pts_[ki].first.get() == face)
	    return true;

    return false;
}
//===========================================================================
void
ftSurfaceSetPoint::addPair(shared_ptr<ftFaceBase>& face, const Vector2D& par_pt)
//===========================================================================
{
    par_pts_.push_back(std::make_pair(face, par_pt));
}

//===========================================================================
int ftSurfaceSetPoint::nmbFaces()
//===========================================================================
{
    return (int)par_pts_.size();
}

//===========================================================================
shared_ptr<ftFaceBase> ftSurfaceSetPoint::face(int i)
//===========================================================================
{
    if (i > int(par_pts_.size()) - 1) {
        THROW("Not that many faces.");
    }
	return par_pts_[i].first;
}

//===========================================================================
Vector2D ftSurfaceSetPoint::parValue(int i)
//===========================================================================
{
    if (i > int(par_pts_.size()) - 1) {
        THROW("Not that many faces.");
    }
	return par_pts_[i].second;
}

//===========================================================================
Vector2D ftSurfaceSetPoint::getPar(ftFaceBase* face)
//===========================================================================
{
    size_t ki;
    for (ki = 0; ki < par_pts_.size(); ++ki)
        if ((par_pts_[ki].first).get() == face)
            break;

    if (ki == par_pts_.size()) {
        THROW("Point not member of face");
    }
	return par_pts_[ki].second;
}

//===========================================================================
void ftSurfaceSetPoint::addInfo(ftSurfaceSetPoint* other)
//
// Add face, parameter and connectivity info from the point other to this point
// NB! It is assumed that the position information is consistent
//===========================================================================
{
    size_t ki, kj;

    // Face and parameter information
    for (ki=0; ki<other->par_pts_.size(); ++ki)
    {
	// Check if the face exists already
	for (kj=0; kj<par_pts_.size(); ++kj)
	    if (par_pts_[kj].first.get() == other->par_pts_[ki].first.get())
		break;

	if (kj == par_pts_.size())
	    addPair(other->par_pts_[ki].first, other->par_pts_[ki].second);
    }

    // Connectivity information
    vector<PointIter> neighbours = other->getNeighbours();
    for (ki=0; ki<neighbours.size(); ++ki)
    {
	neighbours[ki]->addNeighbour(this);
	addNeighbour(neighbours[ki]);
    }
}

//===========================================================================
void ftSurfaceSetPoint::resetPosition(Vector3D pos, int bnd)
//
//===========================================================================
{
    xyz_ = pos;
    uv_ = Vector2D(0.0, 0.0);
    dist_ = -1.0;
    at_boundary_ = bnd;

    // Iterate to find parameter values
    Point pos2(pos.begin(), pos.end());
    Point close_pt;
    double par_u, par_v, dist;
    double eps = 1.0e-12;  // A small number
    for (size_t ki=0; ki<par_pts_.size(); ++ki)
    {
	par_pts_[ki].first->closestPoint(pos2, par_u, par_v, close_pt, dist, eps);
	Vector2D param(par_u, par_v);
	par_pts_[ki].second = param;
    }
}
	
//===========================================================================
    void ftSurfaceSetPoint::addFace(shared_ptr<ftFaceBase>& face)
//
//===========================================================================
{
    // Iterate to find parameter values
    Point pos2(xyz_.begin(), xyz_.end());
    Point close_pt;
    double par_u, par_v, dist;
    double eps = 1.0e-12;  // A small number

    face->closestPoint(pos2, par_u, par_v, close_pt, dist, eps);
    Vector2D param(par_u, par_v);
    addPair(face, param);
}


//===========================================================================
void ftSurfaceSetPoint::write2Dval(std::ostream& os) const
//===========================================================================
{
  par_pts_[0].second.write(os);
}

} // namespace Go
