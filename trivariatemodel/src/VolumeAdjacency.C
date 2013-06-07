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

#include "GoTools/trivariatemodel/VolumeAdjacency.h"
#include "GoTools/trivariate/SurfaceOnVolume.h"
#include "GoTools/compositemodel/ftFaceBase.h"
#include "GoTools/intersections/Identity.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/BoundedUtils.h"
#include <fstream>

//#define DEBUG_VOL

using std::vector;

using namespace Go;

//---------------------------------------------------------------------------
VolumeAdjacency::VolumeAdjacency(double gap, double neighbour)
  : gap_(gap), neighbour_(neighbour)
//---------------------------------------------------------------------------
{
}


//---------------------------------------------------------------------------
VolumeAdjacency::~VolumeAdjacency()
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void VolumeAdjacency::setAdjacency(vector<shared_ptr<Body> >& solids)
//---------------------------------------------------------------------------
{
  int num_solids = (int)solids.size();
  vector<BoundingBox> boxes;
  boxes.reserve(num_solids);
  for (int i = 0; i < num_solids; ++i)
    boxes.push_back(solids[i]->boundingBox());

  // For every combination of solids, do a box test to check if adjacency
  // is possible
  for (int i = 0; i < num_solids; ++i)
    for (int j = i+1; j < num_solids; ++j)
      if (boxes[i].overlaps(boxes[j], neighbour_))
	// Adjacency is possible.
	setAdjacency(solids[i], solids[j]);
}

//---------------------------------------------------------------------------
void VolumeAdjacency::setAdjacency(vector<shared_ptr<Body> >& solids, int new_solid_pos)
//---------------------------------------------------------------------------
{
  int num_solids = (int)solids.size();
  vector<BoundingBox> boxes;
  boxes.reserve(num_solids);
  for (int i = 0; i < num_solids; ++i)
    boxes.push_back(solids[i]->boundingBox());
  BoundingBox new_box = boxes[new_solid_pos];
  shared_ptr<Body> new_solid = solids[new_solid_pos];

  // For every combination of new solid with the other solids, do a box test to check if adjacency
  // is possible
  for (int i = 0; i < num_solids; ++i)
    if (i != new_solid_pos && new_box.overlaps(boxes[i], neighbour_))
      // Adjacency is possible.
      setAdjacency(solids[i], new_solid);
}

//---------------------------------------------------------------------------
void VolumeAdjacency::setAdjacency(shared_ptr<Body> solid1, shared_ptr<Body> solid2)
//---------------------------------------------------------------------------
{
  int kr, kh, idx1, idx2;

  // Check boundary faces
  int nmb1 = solid1->nmbOfShells();
  int nmb2 = solid2->nmbOfShells();
  for (kr=0; kr<nmb1; ++kr)
    {
      shared_ptr<SurfaceModel> shell1 = solid1->getShell(kr);
      int nmb_faces1 = shell1->nmbEntities();
      for (kh=0; kh<nmb2; ++kh)
	{
	  shared_ptr<SurfaceModel> shell2 = solid2->getShell(kh);
	  int nmb_faces2 = shell2->nmbEntities();

	  // For each combination of boundary faces, check adjacency
	  for (idx1=0; idx1<nmb_faces1; ++idx1)
	    {
	      shared_ptr<ftSurface> face1 = shell1->getFace(idx1);
	      BoundingBox box1 = face1->boundingBox();
	      for (idx2=0; idx2<nmb_faces2; ++idx2)
		{
		  shared_ptr<ftSurface> face2 = shell2->getFace(idx2);

		  // Check if the faces are connected already
		  if (face1->twin() && face2->twin() &&
		      face1->twin() == face2.get())
		    continue;

		  // Box testing to check if the faces may be adjacent
		  BoundingBox box2 = face2->boundingBox();
		  if (!box1.overlaps(box2, neighbour_))
		    continue;

		  // Adjacency analysis of the current faces
		  vector<shared_ptr<ftSurface> > new_faces1;
		  vector<shared_ptr<ftSurface> > new_faces2;
		  // int is_changed = faceAdjacency(face1, face2, new_faces1, 
		  // 				 new_faces2);
		  faceAdjacency(face1, face2, new_faces1, 
				new_faces2);

		  // Update involved shells
		  if (new_faces1.size() > 0)
		    {
		      shell1->removeFace(face1);
		      shell1->append(new_faces1);
		      nmb_faces1 = shell1->nmbEntities();
		    }
		  if (new_faces2.size() > 0)
		    {
		      shell2->removeFace(face2);
		      shell2->append(new_faces2);
		      nmb_faces2 = shell2->nmbEntities();
		    }
		  if (new_faces1.size() > 0 || new_faces2.size() > 0)
		    {
		      idx1--;
		      break;
		    }

		}
	    }		
	}
    }
}

//---------------------------------------------------------------------------
int VolumeAdjacency::faceAdjacency(shared_ptr<ftSurface> face1, 
			       shared_ptr<ftSurface> face2,
			       vector<shared_ptr<ftSurface> >& new_faces1,
			       vector<shared_ptr<ftSurface> >& new_faces2)
//---------------------------------------------------------------------------
{
  // @@@ VSK 0110
  // For the time being only identical and embedded surfaces are handled,
  // not overlapping surfaces

  Identity ident;
  shared_ptr<ParamSurface> srf1 = face1->surface();
  shared_ptr<ParamSurface> srf2 = face2->surface();

#ifdef DEBUG_VOL
  std::ofstream out_file("face_pair.g2");
  shared_ptr<SurfaceOnVolume> vol_sf1 =
    dynamic_pointer_cast<SurfaceOnVolume, ParamSurface>(srf1);
  shared_ptr<SurfaceOnVolume> vol_sf2 =
    dynamic_pointer_cast<SurfaceOnVolume, ParamSurface>(srf2);
  if (vol_sf1.get() && vol_sf2.get())
    {
      if (vol_sf1->spaceSurface().get())
	{
	  vol_sf1->spaceSurface()->writeStandardHeader(out_file);
	  vol_sf1->spaceSurface()->write(out_file);
	}
      if (vol_sf2->spaceSurface().get())
	{
	  vol_sf2->spaceSurface()->writeStandardHeader(out_file);
	  vol_sf2->spaceSurface()->write(out_file);
	}
    }
  else
    {
      srf1->writeStandardHeader(out_file);
      srf1->write(out_file);
      srf2->writeStandardHeader(out_file);
      srf2->write(out_file);
    }
#endif

  int res = ident.identicalSfs(srf1, srf2, neighbour_);
  if (res == 1)
    {
      // Coincidence
      // Set adjacency
      face1->connectTwin(face2.get(), neighbour_);
      return 1;
    }
  else if (res == 2)
    {
      // Face 1 is embedded in face 2
      // Trim face 2 according to the boundary curves of face1
      vector<shared_ptr<ParamSurface> > new_sfs;
      splitSurface(srf1, srf2, new_sfs);

      for (size_t kj=0; kj<new_sfs.size(); ++kj)
	{
	  shared_ptr<ftSurface> curr_face = 
	    shared_ptr<ftSurface>(new ftSurface(new_sfs[kj], -1));
	  curr_face->createInitialEdges(neighbour_);
	  new_faces2.push_back(curr_face);
	}

      face1->connectTwin(new_faces2[0].get(), neighbour_);
      return 2;
    }
  else if (res == 3)
    {
      // Face 2 is embedded in face 1
      // Trim face 1 according to the boundary curves of face1
      vector<shared_ptr<ParamSurface> > new_sfs;
      splitSurface(srf2, srf1, new_sfs);

      for (size_t kj=0; kj<new_sfs.size(); ++kj)
	{
	  shared_ptr<ftSurface> curr_face = 
	    shared_ptr<ftSurface>(new ftSurface(new_sfs[kj], -1));
	  curr_face->createInitialEdges(neighbour_);
	  new_faces1.push_back(curr_face);
	}

      face2->connectTwin(new_faces1[0].get(), neighbour_);
      return 3;
    }
  else
    return 0; // No adjacency
}

//---------------------------------------------------------------------------
void VolumeAdjacency::splitSurface(shared_ptr<ParamSurface> srf1,
				   shared_ptr<ParamSurface> srf2,
				   vector<shared_ptr<ParamSurface> >& new_sfs)
//---------------------------------------------------------------------------
{
  // srf1 is embedded in srf 2. Split srf2 accordingly

  shared_ptr<BoundedSurface> tmp_sf2 = 
    shared_ptr<BoundedSurface>(new BoundedSurface(srf2, gap_));

  // Fetch loop from face 1 and represent as surface curves with
  // regard to face 2
  vector<CurveLoop> loops1 = SurfaceTools::allBoundarySfLoops(srf1, gap_);

  vector<vector<shared_ptr<CurveOnSurface> > > bd_loops;
  double eps = 0.0;
  size_t kj;
  for (kj=0; kj<loops1.size(); ++kj)
    {
      eps = std::max(eps, loops1[kj].getSpaceEpsilon());
      vector<shared_ptr<CurveOnSurface> > new_bd;
      int nmb_bd = loops1[kj].size();
      int ki;
      for (ki=0; ki<nmb_bd; ++ki)
	{
	  shared_ptr<CurveOnSurface> bd_cv = 
	    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(loops1[kj][ki]);
	  shared_ptr<CurveOnSurface> curr_cv;
	  if (bd_cv)
	    curr_cv = 
	      shared_ptr<CurveOnSurface>(new CurveOnSurface(tmp_sf2->underlyingSurface(),
							    bd_cv->spaceCurve(),
							    false));
	  else
	    curr_cv = 
	      shared_ptr<CurveOnSurface>(new CurveOnSurface(tmp_sf2->underlyingSurface(),
							    loops1[kj][ki],
							    false));
	  curr_cv->ensureParCrvExistence(gap_);
	  new_bd.push_back(curr_cv);
	}
      bd_loops.push_back(new_bd);
    }

  // Make trimmed surface corresponding to srf1
  shared_ptr<BoundedSurface> new_sf2 = 
    shared_ptr<BoundedSurface>(new BoundedSurface(tmp_sf2->underlyingSurface(),
						  bd_loops, eps));

  // Make remaining surfaces
  vector<shared_ptr<BoundedSurface> > rest_sfs = 
    BoundedUtils::subtractSfPart(tmp_sf2, bd_loops[0], gap_);

  new_sfs.push_back(new_sf2);
  new_sfs.insert(new_sfs.end(), rest_sfs.begin(), rest_sfs.end());
}
