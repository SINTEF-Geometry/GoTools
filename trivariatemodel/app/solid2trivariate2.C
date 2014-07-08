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
#include "GoTools/compositemodel/CompositeModelFactory.h"

using namespace Go;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;

int main(int argc, char* argv[] )
{

  if (argc != 5)
    {
      cout << "Usage: " << "<infile1> <infile2> <outfile> <block structuring mode (1,2,3)>" << endl;
      exit(-1);
    }

  ifstream infile1(argv[1]);
  ALWAYS_ERROR_IF(infile1.bad(), "Bad or no input filename");

  ifstream infile2(argv[2]);
  ALWAYS_ERROR_IF(infile2.bad(), "Bad or no input filename");

 
  ofstream outfile(argv[3]);
  int split_mode = atoi(argv[4]);
  if (split_mode < 1 || split_mode > 3)
    split_mode = 1;  // Default

  // The tolerances must be set according to the properties of the model.
  // The neighbour tolerance must be smaller than the smallest entity in the
  // model, but larger than the largest gap.
  // The gap tolerance must be smaller than the neighbour tolerance
  double gap = 0.0001; //0.001;
  double neighbour = 0.001; //0.01;
  double kink = 0.01;
  double approxtol = 0.01;
  int degree = 3;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  CompositeModel *model = factory.createFromG2(infile1);

  shared_ptr<SurfaceModel> sfmodel = 
    shared_ptr<SurfaceModel>(dynamic_cast<SurfaceModel*>(model));
  if (!sfmodel.get())
    {
      std::cout << "No input model read" << std::endl;
      exit(-1);
    }
 
  if (sfmodel->nmbBoundaries() > 0)
    {
      std::cout << "Not a brep solid. Consider increasing the neighbour tolerance" << std::endl;
      exit(-1);
    }
      
  shared_ptr<ftVolume> ftvol = 
    shared_ptr<ftVolume>(new ftVolume(sfmodel));

  CompositeModel *model2 = factory.createFromG2(infile2);

  shared_ptr<SurfaceModel> sfmodel2 = 
    shared_ptr<SurfaceModel>(dynamic_cast<SurfaceModel*>(model2));
  if (!sfmodel2.get())
    {
      std::cout << "No splitting surface read" << std::endl;
      exit(-1);
    }

  vector<shared_ptr<ftVolume> > vols;
  vols.push_back(ftvol);
  int nmb = sfmodel2->nmbEntities();
  for (int ki=0; ki<nmb; ++ki)
    {
      shared_ptr<ftSurface> face = sfmodel2->getFace(ki);

      size_t nmb_vols = vols.size();
      for (size_t kj=0; kj<nmb_vols; )
	{
	  vector<int> dummy;
	  vector<shared_ptr<ftVolume> > vols2 = 
	    ftVolumeTools::splitVolumes(vols[kj], face, gap/*, dummy*/);
	  std::cout << "Number of volumes: " << vols2.size() << std::endl;

	  if (vols2.size() > 1)
	    {
	      vols.erase(vols.begin()+kj);
	      nmb_vols--;
	      vols.insert(vols.end(), vols2.begin(), vols2.end());
	    }
	  else
	    kj++;
	}
    }

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
    shared_ptr<VolumeModel>(new VolumeModel(vols, gap, kink));

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
  for (int ki=0; ki<nmb_bd; ++ki)
    {
      shared_ptr<ParamSurface> bd_sf = outer_bd->getSurface(ki);
      bd_sf->writeStandardHeader(of3);
      bd_sf->write(of3);
    }

  // Regularize
  //volmod->regularizeBdShells();
  volmod->replaceNonRegVolumes(degree, split_mode);

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


  outer_bd = volmod->getOuterBoundary(0);
  nmb_bd = outer_bd->nmbEntities();
  std::ofstream of4("Vol_bd2.g2");
  for (int ki=0; ki<nmb_bd; ++ki)
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

  //volmod->replaceNonRegVolumes();

  std::cout << "Number of volumes: " << volmod->nmbEntities() << std::endl;
	  
  std::ofstream of6("output_volumes.g2");
  int nmb_vols = volmod->nmbEntities();
  for (int kr=0; kr<nmb_vols; ++kr)
    {
      // if (kr == 19 || kr == 21 || kr == 22)
      // 	continue;

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
	      
	}
      shared_ptr<ParamVolume> curr_vol2 = volmod->getVolume(kr);
      curr_vol2->writeStandardHeader(of6);
      curr_vol2->write(of6);
    }
    
  volmod->makeCommonSplineSpaces();
  volmod->averageCorrespondingCoefs();

  nmb_vols = volmod->nmbEntities();
  for (int kr=0; kr<nmb_vols; ++kr)
    {
      shared_ptr<ParamVolume> curr_vol2 = volmod->getVolume(kr);
      vector<ftVolume*> ng3;
      volmod->getBody(kr)->getAdjacentBodies(ng3);
      std::cout << "Vol nr" << kr << ", nmb neighbours: " << ng3.size() << std::endl;
      curr_vol2->writeStandardHeader(outfile);
      curr_vol2->write(outfile);
    }

}
  
