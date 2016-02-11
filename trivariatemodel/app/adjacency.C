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

#include <fstream>
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/geometry/Utils.h"

using namespace Go;
using std::cout;
using std::endl;
using std::ifstream;
using std::vector;
using std::cin;

int main(int argc, char* argv[] )
{
  if (argc != 3 && argc != 4)
      cout << "Usage: " << "infile2, average or not, (Insert knots)" << endl;

  ifstream is2(argv[1]);
  ALWAYS_ERROR_IF(is2.bad(), "Bad or no input filename");
  int average = atoi(argv[2]);
  int insert = 0;
  if (argc == 4)
    insert = atoi(argv[3]);

  vector<shared_ptr<ftVolume> > volumes;
  
  double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double bend = 0.05;

  int ki;
  while (!is2.eof())
    {
      // Read volume from file
      ObjectHeader head;
      is2 >> head;
      shared_ptr<SplineVolume> vol2;
      vol2 = shared_ptr<SplineVolume>(new SplineVolume());
      vol2->read(is2);

      shared_ptr<ParamVolume> pvol
          = dynamic_pointer_cast<ParamVolume, SplineVolume>(vol2);
      volumes.push_back(shared_ptr<ftVolume>(new ftVolume(pvol, gap, neighbour,
							  kink, bend)));
      // volumes.push_back(shared_ptr<ftVolume>(new ftVolume(vol2)));

      Utils::eatwhite(is2);
    }

  shared_ptr<VolumeModel> model = 
    shared_ptr<VolumeModel>(new VolumeModel(volumes, gap, neighbour, 
					    kink, bend));

  if (model && insert)
    {
      cout << "Number of volumes: " << model->nmbEntities() << std::endl;
      cout << "volume to refine: ";
      int idx;
      cin >> idx;

      cout << "Parameter direction: ";
      int dir;
      cin >> dir;

      cout << "Number of knots: ";
      int nmb;
      cin >> nmb;

      cout << "Knots: ";
      vector<double> knots(nmb);
      for (int ki=0; ki<nmb; ++ki)
	cin >> knots[ki];

      shared_ptr<ParamVolume> vol = model->getVolume(idx);
      shared_ptr<SplineVolume> v1 = 
	dynamic_pointer_cast<SplineVolume,ParamVolume>(vol);
      if (v1.get())
	{
	  v1->insertKnot(dir, knots);
	}
    }

  model->makeCommonSplineSpaces();
  if (average)
    model->averageCorrespondingCoefs();

  std::ofstream out_file2("adjacency.dat");
  int nmb = model->nmbEntities();
  for (ki=0; ki<nmb; ++ki)
    {
      shared_ptr<ftVolume> body = model->getBody(ki);
      shared_ptr<ParamVolume> vol1 = model->getVolume(ki);
      shared_ptr<SplineVolume> svol1 = 
	dynamic_pointer_cast<SplineVolume, ParamVolume>(vol1);
      vector<double>::const_iterator coef1 = (svol1->rational()) ?
	svol1->rcoefs_begin() : svol1->coefs_begin();
      vector<double>::const_iterator c1 = svol1->coefs_begin();
      int dim = svol1->dimension();
      int idim = svol1->dimension() + svol1->rational();
      
      vector<int> free_bd;
      bool found_free = body->getFreeBoundaryInfo(gap, free_bd);
      out_file2 << "Volume " << ki << ", result: " << found_free;
      out_file2 << ", number: " << free_bd.size() << std::endl;
      int kj;
      for (kj=0; kj < (int)free_bd.size(); ++kj)
	out_file2 << free_bd[kj] << "  ";
      out_file2 << std::endl << std::endl;

      for (kj=0; kj < (int)free_bd.size(); ++kj)
	{
	  vector<int> bd_enumeration;
	  found_free = body->getBoundaryCoefEnumeration(free_bd[kj],
							bd_enumeration);
	  out_file2 << "Boundary: " << kj << std::endl;
	  for (size_t kr=0; kr<bd_enumeration.size(); ++kr)
	    out_file2 << bd_enumeration[kr] << " ";
	  out_file2 << std::endl;
	}
      out_file2 << std::endl;

      vector<ftVolume*> neighbours;
      body->getAdjacentBodies(neighbours);
      for (kj=0; kj < (int)neighbours.size(); ++kj)
	{
	  int idx = model->getIndex(neighbours[kj]);
	  std::cout << "Found adjacency info for volume " ;
	  std::cout << ki << " and " << idx << std::endl;
	}

      for (kj=ki+1; kj<nmb; ++kj)
	{
	  shared_ptr<ftVolume> body2 = model->getBody(kj);
	  shared_ptr<ParamVolume> vol2 = model->getVolume(kj);
	  shared_ptr<SplineVolume> svol2 = 
	    dynamic_pointer_cast<SplineVolume, ParamVolume>(vol2);
	  vector<double>::const_iterator coef2 = (svol2->rational()) ?
	    svol2->rcoefs_begin() : svol2->coefs_begin();
	  vector<double>::const_iterator c2 = svol2->coefs_begin();
	  
	  int bd1, bd2;
	  int orientation;
	  bool seq;
	  bool found = body->getAdjacencyInfo(body2.get(), gap, bd1, bd2,
					    orientation, seq);
	  std::cout << "Found adjacency info for volume " ;
	  std::cout << ki << " and " << kj << ": " << found << std::endl;
	  if (found)
	    std::cout << bd1 << " " << bd2 << " " << orientation << " " << seq << std::endl;

	  out_file2 << "Found adjacency info for volume " ;
	  out_file2 << ki << " and " << kj << ": " << found << std::endl;
	  if (found)
	    out_file2 << bd1 << " " << bd2 << " " << orientation << " " << seq << std::endl;

	  vector<pair<int,int> > enumeration;
	  found = body->getCorrCoefEnumeration(body2.get(), gap, enumeration);
	  for (size_t kr=0; kr<enumeration.size(); ++kr)
	    {
	      out_file2 << enumeration[kr].first << " ";
	      out_file2 << enumeration[kr].second << "  ,  ";
	      for (int k2=0; k2<idim; ++k2)
		out_file2 << coef1[idim*enumeration[kr].first+k2] << "  ";
	      out_file2 << ",  ";
	      for (int k2=0; k2<idim; ++k2)
		out_file2 << coef2[idim*enumeration[kr].second+k2] << "  ";
	      double d2 = Utils::distance_squared(&coef1[idim*enumeration[kr].first],
					   &coef1[idim*enumeration[kr].first+idim],
					   &coef2[idim*enumeration[kr].second]);
	      double d3 = Utils::distance_squared(&c1[dim*enumeration[kr].first],
					   &c1[dim*enumeration[kr].first+dim],
					   &c2[dim*enumeration[kr].second]);
	      out_file2 << ", dist2 = " << d2 <<", " << d3 << std::endl;
	    }

	  VolumeAdjacencyInfo info = body->getCornerAdjacencyInfo(body2.get(), gap);
	  std::cout << "Corner found: " << ki << " " << kj << " " << info.corner_adjacency_ << std::endl;
	  if (info.corner_adjacency_)
	    {
	      std::cout << info.bd_idx_1_ << " " << info.bd_idx_2_ << " " << info.edg_idx_1_;
	      std::cout << " " << info.edg_idx_2_ << " " << info.same_orient_edge_ << std::endl;
	    }
	}
      out_file2 << std::endl;
    }

  out_file2 << std::endl;

  vector<shared_ptr<Vertex> > vx;
  model->getAllVertices(vx);
  for (ki=0; ki<(int)vx.size(); ++ki)
    {
      Point pos = vx[ki]->getVertexPoint();
      vector<Body*> bodies = vx[ki]->getBodies();
      for (size_t kj=0; kj<bodies.size(); ++kj)
	{
	  Point param;
	  int corner, coef_nmb;
	  ftVolume *curr = dynamic_cast<ftVolume*>(bodies[kj]);
	  // bool found = curr->getVertexEnumeration(vx[ki], param,
	  //					  corner, coef_nmb);
	  curr->getVertexEnumeration(vx[ki], param,
				     corner, coef_nmb);
	  int idx = model->getIndex(curr);
	  out_file2 << "Volume " << idx <<", vertex: " << pos;
	  out_file2 << ", param: " << param <<", corner: ";
	  out_file2 << corner << " " << coef_nmb << std::endl;
	}
    }
	  
  std::ofstream out_file("volumes.g2");
  for (ki=0; ki<nmb; ++ki)
    {
      shared_ptr<ParamVolume> vol = model->getVolume(ki);
      vol->writeStandardHeader(out_file);
      vol->write(out_file);
    }
  
}
