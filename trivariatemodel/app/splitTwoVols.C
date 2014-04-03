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
#include "GoTools/trivariatemodel/ftVolumeTools.h"
#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/compositemodel/EdgeVertex.h"

using namespace Go;
using std::cout;
using std::endl;
using std::ifstream;

int main(int argc, char* argv[] )
{
#define DEBUG_VOL1

  if (argc != 3)
      cout << "Usage: " << "<infile1> <infile2>" << endl;

  ifstream infile1(argv[1]);
  ALWAYS_ERROR_IF(infile1.bad(), "Bad or no input filename");

  ifstream infile2(argv[2]);
  ALWAYS_ERROR_IF(infile2.bad(), "Bad or no input filename");

 
  ObjectHeader header;
  header.read(infile1);

  shared_ptr<SplineVolume> vol1(new SplineVolume());
  vol1->read(infile1);

  header.read(infile2);

  shared_ptr<SplineVolume> vol2(new SplineVolume());
  vol2->read(infile2);

  double gap_eps = 1.0e-4;
  double kink_eps = 1.0e-2;

  int degree = 3;

  shared_ptr<ftVolume> ftvol1 = 
    shared_ptr<ftVolume>(new ftVolume(vol1, gap_eps, kink_eps));

  shared_ptr<ftVolume> ftvol2 = 
    shared_ptr<ftVolume>(new ftVolume(vol2, gap_eps, kink_eps));

  std::ofstream of1("Vol_1.g2");
  shared_ptr<SurfaceModel> model1 = ftvol1->getOuterShell();
  int nmb1 = model1->nmbEntities();
  int ki;
  for (ki=0; ki<nmb1; ++ki)
    {
      shared_ptr<ParamSurface> sf = model1->getSurface(ki);
      sf->writeStandardHeader(of1);
      sf->write(of1);
    }

  std::ofstream of2("Vol_2.g2");
  shared_ptr<SurfaceModel> model2 = ftvol2->getOuterShell();
  int nmb2 = model2->nmbEntities();
  for (ki=0; ki<nmb2; ++ki)
    {
      shared_ptr<ParamSurface> sf = model2->getSurface(ki);
      sf->writeStandardHeader(of2);
      sf->write(of2);
    }
  
  vector<int> dummy;
  vector<shared_ptr<ftVolume> > vols = 
    ftVolumeTools::splitVolumes(ftvol1, ftvol2, gap_eps, dummy);
  std::cout << "Number of volumes: " << vols.size() << std::endl;

  if (vols.size() > 0)
    {
      std::ofstream of2("TrimVol_1.g2");
      shared_ptr<SurfaceModel> mod = vols[0]->getOuterShell();
      int nmb = mod->nmbEntities();
      int ki;
      for (ki=0; ki<nmb; ++ki)
	{
	  shared_ptr<ParamSurface> sf = mod->getSurface(ki);
	  sf->writeStandardHeader(of2);
	  sf->write(of2);
	}
    }

  if (vols.size() > 1)
    {
      std::ofstream of2("TrimVol_2.g2");
      shared_ptr<SurfaceModel> mod = vols[1]->getOuterShell();
      int nmb = mod->nmbEntities();
      int ki;
      for (ki=0; ki<nmb; ++ki)
	{
	  shared_ptr<ParamSurface> sf = mod->getSurface(ki);
	  sf->writeStandardHeader(of2);
	  sf->write(of2);
	}
    }

  if (vols.size() > 2)
    {
      std::ofstream of2("TrimVol_3.g2");
      shared_ptr<SurfaceModel> mod = vols[2]->getOuterShell();
      int nmb = mod->nmbEntities();
      int ki;
      for (ki=0; ki<nmb; ++ki)
	{
	  shared_ptr<ParamSurface> sf = mod->getSurface(ki);
	  sf->writeStandardHeader(of2);
	  sf->write(of2);
	}
    }

  if (vols.size() > 3)
    {
      std::ofstream of2("TrimVol_4.g2");
      shared_ptr<SurfaceModel> mod = vols[3]->getOuterShell();
      int nmb = mod->nmbEntities();
      int ki;
      for (ki=0; ki<nmb; ++ki)
	{
	  shared_ptr<ParamSurface> sf = mod->getSurface(ki);
	  sf->writeStandardHeader(of2);
	  sf->write(of2);
	}
    }

  if (vols.size() > 4)
    {
      std::ofstream of2("TrimVol_5.g2");
      shared_ptr<SurfaceModel> mod = vols[4]->getOuterShell();
      int nmb = mod->nmbEntities();
      int ki;
      for (ki=0; ki<nmb; ++ki)
	{
	  shared_ptr<ParamSurface> sf = mod->getSurface(ki);
	  sf->writeStandardHeader(of2);
	  sf->write(of2);
	}
    }

  if (vols.size() > 5)
    {
      std::ofstream of2("TrimVol_6.g2");
      shared_ptr<SurfaceModel> mod = vols[5]->getOuterShell();
      int nmb = mod->nmbEntities();
      int ki;
      for (ki=0; ki<nmb; ++ki)
	{
	  shared_ptr<ParamSurface> sf = mod->getSurface(ki);
	  sf->writeStandardHeader(of2);
	  sf->write(of2);
	}
    }

  if (vols.size() > 6)
    {
      std::ofstream of2("TrimVol_7.g2");
      shared_ptr<SurfaceModel> mod = vols[6]->getOuterShell();
      int nmb = mod->nmbEntities();
      int ki;
      for (ki=0; ki<nmb; ++ki)
	{
	  shared_ptr<ParamSurface> sf = mod->getSurface(ki);
	  sf->writeStandardHeader(of2);
	  sf->write(of2);
	}
    }

  if (vols.size() > 7)
    {
      std::ofstream of2("TrimVol_8.g2");
      shared_ptr<SurfaceModel> mod = vols[7]->getOuterShell();
      int nmb = mod->nmbEntities();
      int ki;
      for (ki=0; ki<nmb; ++ki)
	{
	  shared_ptr<ParamSurface> sf = mod->getSurface(ki);
	  sf->writeStandardHeader(of2);
	  sf->write(of2);
	}
    }

  if (vols.size() > 8)
    {
      std::ofstream of2("TrimVol_9.g2");
      shared_ptr<SurfaceModel> mod = vols[8]->getOuterShell();
      int nmb = mod->nmbEntities();
      int ki;
      for (ki=0; ki<nmb; ++ki)
	{
	  shared_ptr<ParamSurface> sf = mod->getSurface(ki);
	  sf->writeStandardHeader(of2);
	  sf->write(of2);
	}
    }

  if (vols.size() > 9)
    {
      std::ofstream of2("TrimVol_10.g2");
      shared_ptr<SurfaceModel> mod = vols[9]->getOuterShell();
      int nmb = mod->nmbEntities();
      int ki;
      for (ki=0; ki<nmb; ++ki)
	{
	  shared_ptr<ParamSurface> sf = mod->getSurface(ki);
	  sf->writeStandardHeader(of2);
	  sf->write(of2);
	}
    }

  // Make volume model
  shared_ptr<VolumeModel> volmod =
    shared_ptr<VolumeModel>(new VolumeModel(vols, gap_eps, kink_eps));

  // // Remove inner piece
  volmod->removeSolid(vols[3]);

  std::ofstream of12("Ver1.g2");
  int n1 = volmod->nmbEntities();
  for (int kr=0; kr<n1; ++kr)
    {
      shared_ptr<SurfaceModel> mod = volmod->getBody(kr)->getOuterShell();
      int nmb = mod->nmbEntities();
      int ki;
      for (ki=0; ki<nmb; ++ki)
	{
	  shared_ptr<ParamSurface> sf = mod->getSurface(ki);
	  sf->writeStandardHeader(of12);
	  sf->write(of12);
	}
    }

  shared_ptr<SurfaceModel> outer_bd = volmod->getOuterBoundary(0);
  int nmb_bd = outer_bd->nmbEntities();
  std::ofstream of3("Vol_bd.g2");
  for (ki=0; ki<nmb_bd; ++ki)
    {
      shared_ptr<ParamSurface> bd_sf = outer_bd->getSurface(ki);
      bd_sf->writeStandardHeader(of3);
      bd_sf->write(of3);
    }

  std::ofstream ofradial("radialedges.txt");
  vector<shared_ptr<EdgeVertex> > radial1;
  volmod->getRadialEdges(radial1);
  ofradial << "Radial1: " << std::endl;
  for (size_t kv=0; kv<radial1.size(); ++kv)
    {
      vector<ftEdge*> edges = radial1[kv]->allEdges();
      ofradial << "radial edge " << kv << "  " << radial1[kv].get() << std::endl;
      ofradial << "edges: " << std::endl;
      for (size_t kh=0; kh<edges.size(); ++kh)
	ofradial << edges[kh] << "  ";
      ofradial << std::endl;
    }
      
  vector<shared_ptr<Vertex> > vx1;
  volmod->getAllVertices(vx1);
  ofradial << std::endl << "Vx1: " << std::endl;
  for (size_t kv=0; kv<vx1.size(); ++kv)
    {
      vector<ftEdge*> edges = vx1[kv]->allEdges();
      ofradial << "vertex " << kv << "  " << vx1[kv].get();
      ofradial << "  " << vx1[kv]->getVertexPoint() << std::endl;
      ofradial << "edges: " << std::endl;
      for (size_t kh=0; kh<edges.size(); ++kh)
	ofradial << edges[kh] << "  ";
      ofradial << std::endl;
    }

  // Regularize
  volmod->regularizeBdShells();

  std::ofstream of13("Ver2.g2");
  n1 = volmod->nmbEntities();
  for (int kr=0; kr<n1; ++kr)
    {
      shared_ptr<SurfaceModel> mod = volmod->getBody(kr)->getOuterShell();
      int nmb = mod->nmbEntities();
      int ki;
      for (ki=0; ki<nmb; ++ki)
	{
	  shared_ptr<ParamSurface> sf = mod->getSurface(ki);
	  sf->writeStandardHeader(of13);
	  sf->write(of13);
	}
    }

   vector<shared_ptr<EdgeVertex> > radial2;
  volmod->getRadialEdges(radial2);
  ofradial << std::endl << "Radial2: " << std::endl;
  for (size_t kv=0; kv<radial2.size(); ++kv)
    {
      vector<ftEdge*> edges = radial2[kv]->allEdges();
      ofradial << "radial edge " << kv << "  " << radial2[kv].get() << std::endl;
      ofradial << "edges: " << std::endl;
      for (size_t kh=0; kh<edges.size(); ++kh)
	ofradial << edges[kh] << "  ";
      ofradial << std::endl;
    }

  volmod->getAllVertices(vx1);
  ofradial << std::endl << "Vx2: " << std::endl;
  for (size_t kv=0; kv<vx1.size(); ++kv)
    {
      vector<ftEdge*> edges = vx1[kv]->allEdges();
      ofradial << "vertex " << kv << "  " << vx1[kv].get();
      ofradial << "edges: " << std::endl;
      for (size_t kh=0; kh<edges.size(); ++kh)
	ofradial << edges[kh] << "  ";
      ofradial << std::endl;
    }

    outer_bd = volmod->getOuterBoundary(0);
  nmb_bd = outer_bd->nmbEntities();
  std::ofstream of4("Vol_bd2.g2");
  for (ki=0; ki<nmb_bd; ++ki)
    {
      shared_ptr<ParamSurface> bd_sf = outer_bd->getSurface(ki);
      bd_sf->writeStandardHeader(of4);
      bd_sf->write(of4);
    }
  std::ofstream of5("Notreg.g2");
  for (size_t kr=0; kr<vols.size(); ++kr)
    {
      bool reg = vols[kr]->isRegularized();
      if (!reg)
	{
	  shared_ptr<SurfaceModel> mod = vols[kr]->getOuterShell();
	  int nmb = mod->nmbEntities();
	  int ki;
	  for (ki=0; ki<nmb; ++ki)
	    {
	      shared_ptr<ParamSurface> sf = mod->getSurface(ki);
	      sf->writeStandardHeader(of5);
	      sf->write(of5);
	    }
	}
    }

  volmod->replaceNonRegVolumes();

  vector<shared_ptr<EdgeVertex> > radial3;
  volmod->getRadialEdges(radial3);
  ofradial << "Radial3: " << std::endl;
  for (size_t kv=0; kv<radial3.size(); ++kv)
    {
      vector<ftEdge*> edges = radial3[kv]->allEdges();
      ofradial << "radial edge " << kv << "  " << radial3[kv].get() << std::endl;
      ofradial << "edges: " << std::endl;
      for (size_t kh=0; kh<edges.size(); ++kh)
	ofradial << edges[kh] << "  ";
      ofradial << std::endl;
    }

  volmod->getAllVertices(vx1);
  ofradial << std::endl << "Vx3: " << std::endl;
  for (size_t kv=0; kv<vx1.size(); ++kv)
    {
      vector<ftEdge*> edges = vx1[kv]->allEdges();
      ofradial << "vertex " << kv << "  " << vx1[kv].get() << std::endl;
      ofradial << "edges: " << std::endl;
      for (size_t kh=0; kh<edges.size(); ++kh)
	ofradial << edges[kh] << "  ";
      ofradial << std::endl;
    }
  std::cout << "Number of volumes: " << volmod->nmbEntities() << std::endl;
	  
  std::ofstream of6("output_volumes.g2");
  int nmb_vols = volmod->nmbEntities();
  for (int kr=0; kr<nmb_vols; ++kr)
    {
      shared_ptr<ftVolume> curr_vol = volmod->getBody(kr);
      bool bd_trim = curr_vol->isBoundaryTrimmed();
      bool iso_trim = curr_vol->isIsoTrimmed();
      bool reg = curr_vol->isRegularized();

      std::cout << "Volume nr " << kr << ": " << bd_trim;
      std::cout << " " << iso_trim << " " << reg << std::endl;

      std::ofstream of7("Curr_vol.g2");
      shared_ptr<SurfaceModel> mod = curr_vol->getOuterShell();
      int nmb = mod->nmbEntities();
      int ki;
      for (ki=0; ki<nmb; ++ki)
	{
	  shared_ptr<ParamSurface> sf = mod->getSurface(ki);
	  sf->writeStandardHeader(of7);
	  sf->write(of7);
	}

      if (reg)
	{
	  vector<ftVolume*> ng1;
	  curr_vol->getAdjacentBodies(ng1);
	  std::cout << "Number of neighbours before untrim: " << ng1.size() << std::endl;
	  curr_vol->untrimRegular(degree);
	  vector<ftVolume*> ng2;
	  curr_vol->getAdjacentBodies(ng2);
	  std::cout << "Number of neighbours after untrim: " << ng2.size() << std::endl;

	  std::ofstream of11("adj_vol.g2");
	  shared_ptr<SurfaceModel> mod = curr_vol->getOuterShell();
	  int nmb = mod->nmbEntities();
	  int ki;
	  for (ki=0; ki<nmb; ++ki)
	    {
	      shared_ptr<ParamSurface> sf = mod->getSurface(ki);
	      sf->writeStandardHeader(of11);
	      sf->write(of11);
	    }
	  for (size_t kf=0; kf<ng2.size(); ++kf)
	    {
	      if (!ng2[kf])
		continue;
	      mod = ng2[kf]->getOuterShell();
	      nmb = mod->nmbEntities();
	      for (ki=0; ki<nmb; ++ki)
		{
		  shared_ptr<ParamSurface> sf = mod->getSurface(ki);
		  sf->writeStandardHeader(of11);
		  sf->write(of11);
		}
	    }
	      
	    vector<shared_ptr<EdgeVertex> > radial4;
	    volmod->getRadialEdges(radial4);
	    ofradial << "Radial4: " << std::endl;
	    std::cout << "Radial4: " << std::endl;
	    for (size_t kv=0; kv<radial4.size(); ++kv)
	      {
		vector<ftEdge*> edges = radial4[kv]->allEdges();
		ofradial << "radial edge " << kv << "  " << radial4[kv].get() << std::endl;
		ofradial << "edges: " << std::endl;
		for (size_t kh=0; kh<edges.size(); ++kh)
		  ofradial << edges[kh] << "  ";
		ofradial << std::endl;
	      }

	    volmod->getAllVertices(vx1);
	    ofradial << std::endl << "Vx4: " << std::endl;
	    for (size_t kv=0; kv<vx1.size(); ++kv)
	      {
		vector<ftEdge*> edges = vx1[kv]->allEdges();
		ofradial << "vertex " << kv << "  " << vx1[kv].get() << std::endl;
		ofradial << "edges: " << std::endl;
		for (size_t kh=0; kh<edges.size(); ++kh)
		  ofradial << edges[kh] << "  ";
		ofradial << std::endl;
	      }
	}
      shared_ptr<ParamVolume> curr_vol2 = volmod->getVolume(kr);
      curr_vol2->writeStandardHeader(of6);
      curr_vol2->write(of6);
    }
    
  std::ofstream of8("output_volumes_2.g2");
  volmod->makeCommonSplineSpaces();
  volmod->averageCorrespondingCoefs();

  nmb_vols = volmod->nmbEntities();
  for (int kr=0; kr<nmb_vols; ++kr)
    {
      shared_ptr<ParamVolume> curr_vol2 = volmod->getVolume(kr);
      vector<ftVolume*> ng3;
      volmod->getBody(kr)->getAdjacentBodies(ng3);
      std::cout << "Vol nr" << kr << ", nmb neighbours: " << ng3.size() << std::endl;
      curr_vol2->writeStandardHeader(of8);
      curr_vol2->write(of8);
    }

}
