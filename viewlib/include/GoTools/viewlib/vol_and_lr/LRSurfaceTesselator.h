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

#ifndef LRSURFACETESSELATOR_H
#define LRSURFACETESSELATOR_H

#include "GoTools/tesselator/Tesselator.h"
#include "GoTools/tesselator/RegularMesh.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/RectDomain.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/tesselator/GenericTriMesh.h"
#include <memory>
#include "GoTools/utils/config.h"

namespace Go
{

/** LRSurfaceTesselator: create a mesh for an LR surface with a suitable
    triangulation. Trimmed surfaces are currently not supported.
*/
  class Element2D;
  class LRSplineSurface;

class GO_API LRSurfaceTesselator : public Tesselator
{
  
public:
  /// Constructor. Surface and mesh size are given. The mesh size relates to 
  /// the underlying surface in the case of bounded surfaces.
    LRSurfaceTesselator(const ParamSurface& surf)
	: surf_(surf), m_(20), n_(20)
    {
 	mesh_ = shared_ptr<GenericTriMesh>(new GenericTriMesh(0,0,true,true));
    }

    virtual ~LRSurfaceTesselator();

    virtual void tesselate();

    // virtual GeneralMesh* getMesh()
    // {
    // 	return mesh_.get();
    // }

    /// Fetch the resulting mesh
    shared_ptr<GenericTriMesh> getMesh()
    {
	return mesh_;
    }

    /// Change mesh size
    void changeRes(int n, int m);

    /// Fetch info about mesh size
    void getRes(int& n, int& m)
    {
	m = m_;
	n = n_;
    }

private:
    const ParamSurface& surf_;
    shared_ptr<GenericTriMesh> mesh_;
    int m_;
    int n_;

    void identifyNearTriangles(std::pair<Point,Point> gap_par,
			       std::vector<Point>& vx_par,
			       std::vector<unsigned int>& tri,
			       double del,
			       std::vector<size_t>& tri_near,
			       std::vector<BoundingBox>& bb_near);

    void identifyJointTriangles(const Point& joint,
			       std::vector<Point>& vx_par,
			       std::vector<unsigned int>& tri,
			       double del,
			       std::vector<size_t>& tri_near,
			       std::vector<BoundingBox>& bb_near);
    
    void  updateVertex(Point& par, const LRSplineSurface* lrsf,
		       RectDomain& dom, Point& vertex_par,
		       Point& vertex, Point& vertex_norm,
		       Point& vertex_tex, Element2D *elem=0);


    double edgeDist(Point curr, Point pa, Point pb, Point& pp);

    bool isInsideTri(Point par, Point pa, Point pb, Point pc);
    
    void getNextTriangle(int start_ix1, int start_ix2,
			 Point par1, Point par2,
			 std::vector<size_t>& tri_near,
			 std::vector<unsigned int>& tri,
			 std::vector<Point>& vx_par, double pdel,
			 int& tri_ix1, int& tri_ix2, Point& pp, int& vx_ix, 
			 int& next_ix1, int& next_ix2);

    int parConfiguration(const Point& param,
			 std::vector<Point>& vx_par,
			 std::vector<unsigned int>& tri,
			 double del,
			 std::vector<size_t>& tri_near,
			 int& ix1, int& ix2, 
			 int edge[], Point& pp);

    void sortAndConnect(std::vector<Point>& vx_par, std::vector<int>& vx_ix,
			std::vector<unsigned int*> triangle);

    void jointSplitAtEdge(const Point& joint,
			  std::vector<Point>& pvx,
			  std::vector<int>& vx_at,
			  double del, std::vector<Point>& par,
			  std::vector<int>& par_at,
			  std::vector<std::pair<int, int> >& parquart);

    void writeGapTri(std::ofstream& oftg, std::vector<Point>& vertex,
		     unsigned int ix1, unsigned int ix2,
		     unsigned int ix3);
};

} // namespace Go




#endif //  LRSURFACETESSELATOR_H

