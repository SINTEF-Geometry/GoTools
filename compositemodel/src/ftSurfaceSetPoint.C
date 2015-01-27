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
      if (neighbours[ki] == this)
	continue;
      neighbours[ki]->addNeighbour(this);
      addNeighbour(neighbours[ki]);
    }
}

//===========================================================================
void ftSurfaceSetPoint::resetPosition(Vector3D pos, int bnd)
//
//===========================================================================
{
  at_boundary_ = bnd;
  resetPosition(pos);
}

//===========================================================================
void ftSurfaceSetPoint::resetPosition(Vector3D pos)
//
//===========================================================================
{
    xyz_ = pos;
    uv_ = Vector2D(0.0, 0.0);
    dist_ = -1.0;

    // Iterate to find parameter values
    Point pos2(pos.begin(), pos.end());
    Point close_pt;
    double par_u, par_v, dist;
    double eps = 1.0e-12;  // A small number
    for (size_t ki=0; ki<par_pts_.size(); ++ki)
    {
      if (par_pts_[ki].first)
	{
	  par_pts_[ki].first->closestPoint(pos2, par_u, par_v, close_pt, dist, eps);
	  Vector2D param(par_u, par_v);
	  par_pts_[ki].second = param;
	}
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
