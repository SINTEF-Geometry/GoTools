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
#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/RegularizeFaceSet.h"

using namespace Go;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::vector;

int main(int argc, char* argv[] )
{
  if (argc != 4)
    {
      cout << "Usage: " << "<infile> <outfile> <block structuring mode (1,2,3)>" << endl;
      exit(-1);
    }

  ifstream infile(argv[1]);
  ALWAYS_ERROR_IF(infile.bad(), "Bad or no input filename");

  ofstream outfile(argv[2]);
  int split_mode = atoi(argv[3]);
  if (split_mode < 1 || split_mode > 3)
    split_mode = 1;  // Default

  // The tolerances must be set according to the properties of the model.
  // The neighbour tolerance must be smaller than the smallest entity in the
  // model, but larger than the largest gap.
  // The gap tolerance must be smaller than the neighbour tolerance
  double gap = 0.0001; //0.001;
  double neighbour = 0.001; //0.01;
  double kink = 0.01;
  double approxtol = 0.001;
  int degree = 3;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  CompositeModel *model = factory.createFromG2(infile);

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
      
  bool isOK = sfmodel->checkShellTopology();
  std::cout << "Shell topology: " << isOK << std::endl;

  // RegularizeFaceSet regularize(sfmodel);
  // shared_ptr<SurfaceModel> sfmodel2 = regularize.getRegularModel();
  
  // std::ofstream of6_1("bd_split.g2");
  // int nmb = sfmodel2->nmbEntities();
  // int ki;
  // for (ki=0; ki<nmb; ++ki)
  //   {
  //     shared_ptr<ParamSurface> sf = sfmodel2->getSurface(ki);
  //     sf->writeStandardHeader(of6_1);
  //     sf->write(of6_1);
  //   }

  // shared_ptr<ftVolume> ftvol = 
  //   shared_ptr<ftVolume>(new ftVolume(sfmodel2));
  shared_ptr<ftVolume> ftvol = 
    shared_ptr<ftVolume>(new ftVolume(sfmodel));

  int nmb;
  int ki;
  shared_ptr<VolumeModel> volmod;
  bool reg = ftvol->isRegularized();
  bool pattern_split = false; //true;
  if (!reg)
    {
      vector<SurfaceModel*> modified_adjacent;
      vector<shared_ptr<ftVolume> > reg_vols = 
	ftvol->replaceWithRegVolumes(degree, modified_adjacent,
				     false, split_mode, pattern_split);

      // // Check each entity
      // nmb = (int)reg_vols.size();
      // int kn;
      // for (kn=0; kn<nmb; ++kn)
      // 	{
      // 	  bool reg2 = reg_vols[kn]->isRegularized();
      // 	  if (!reg2)
      // 	    {
      // 	      std::ofstream of6_2("notreg_vol.g2");
      // 	      shared_ptr<SurfaceModel> mod =  reg_vols[kn]->getOuterShell();
      // 	      int nmb_vol = mod->nmbEntities();
      // 	      for (ki=0; ki<nmb_vol; ++ki)
      // 		{
      // 		  shared_ptr<ParamSurface> sf = mod->getSurface(ki);
      // 		  sf->writeStandardHeader(of6_2);
      // 		  sf->write(of6_2);
      // 		}

      // 	      vector<shared_ptr<ftVolume> > reg_vols2 = 
      // 		reg_vols[kn]->replaceWithRegVolumes(true);
      // 	      if (reg_vols2.size() > 1)
      // 		{
      // 		  reg_vols.erase(reg_vols.begin()+kn);
      // 		  reg_vols.insert(reg_vols.end(), reg_vols2.begin(), 
      // 				  reg_vols2.end());
      // 		  nmb--;
      // 		  kn--;
      // 		}
      // 	    }
      // 	}
      if (reg_vols.size() > 0)
	volmod = shared_ptr<VolumeModel>(new VolumeModel(reg_vols, gap, neighbour,
							 kink, 10.0*kink));
      else
	{
	  vector<shared_ptr<ftVolume> > reg_vols(1);
	  reg_vols[0] = ftvol;
	  volmod = shared_ptr<VolumeModel>(new VolumeModel(reg_vols, gap, neighbour,
							   kink, 10.0*kink));
	}    
    }
  else
    {
      vector<shared_ptr<ftVolume> > reg_vols(1);
      reg_vols[0] = ftvol;
      volmod = shared_ptr<VolumeModel>(new VolumeModel(reg_vols, gap, neighbour,
						       kink, 10.0*kink));
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
      nmb = mod->nmbEntities();
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
	      
	  shared_ptr<ParamVolume> curr_vol2 = volmod->getVolume(kr);
	  curr_vol2->writeStandardHeader(of6);
	  curr_vol2->write(of6);
	}
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
