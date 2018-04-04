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
#include "GoTools/trivariatemodel/VolumeModelFileHandler.h"
#include "GoTools/trivariate/CurveOnVolume.h"
#include "GoTools/trivariate/SurfaceOnVolume.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"

using namespace Go;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::vector;

int main(int argc, char* argv[] )
{
  if (argc != 4 && argc != 5 && argc != 6 && argc != 7)
    {
      cout << "Usage: " << "<infile> <outfile> <block structuring mode (1,2,3)> (<file type, 2 for g22>) (<file type, 2 for output g22>) (<parameter block structuring (0/1)>" << endl;
      exit(-1);
    }

  std::string infile(argv[1]);

  std::ofstream  outfile(argv[2]);
  int split_mode = atoi(argv[3]);
  if (split_mode < 1 || split_mode > 3)
    split_mode = 1;  // Default

  int file_type = 1;
  if (argc > 4)
    file_type = atoi(argv[4]);

  int file_type_out = 1;
  if (argc > 5)
    file_type_out = atoi(argv[5]);

  int block_par = 0;
  if (argc == 7)
    block_par = atoi(argv[6]);

  double gap = 0.001; //0.001;
  double neighbour = 0.01; //0.01;
  double kink = 0.01;
  double approxtol = 0.001;
  int degree = 3;

  shared_ptr<ftVolume> ftvol;
  vector<shared_ptr<ftVolume> > tmpvols;
  if (file_type == 3)
    {
      VolumeModelFileHandler fileread;
      tmpvols = fileread.readVolumes(infile.c_str());
      ftvol = tmpvols[0];
      tpTolerances top = ftvol->getTolerances();
      gap = top.gap;
      neighbour = top.neighbour;
      kink = top.kink;
    }
  else if (file_type == 2)
    {
      VolumeModelFileHandler fileread;
      ftvol = fileread.readVolume(infile.c_str());
      if (!ftvol.get())
	{
	  shared_ptr<Body> body = fileread.readBody(infile.c_str());
	  ftvol = shared_ptr<ftVolume>(new ftVolume(body));
	}

      // Try to simplify model
     ftvol->getShell(0)->simplifyShell();

      tpTolerances top = ftvol->getTolerances();
      gap = top.gap;
      neighbour = top.neighbour;
      kink = top.kink;
    }
  else
    {
      GoTools::init();
      Registrator<SurfaceOnVolume> r211;
      Registrator<CurveOnVolume> r111;

      // The tolerances must be set according to the properties of the model.
      // The neighbour tolerance must be smaller than the smallest entity in the
      // model, but larger than the largest gap.
      // The gap tolerance must be smaller than the neighbour tolerance
      CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

      std::ifstream is(infile);
      CompositeModel *model = factory.createFromG2(is);
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

      // Try to simplify model
      sfmodel->simplifyShell();

      ftvol = shared_ptr<ftVolume>(new ftVolume(sfmodel));
    }

  int nmb;
  int ki;
  shared_ptr<VolumeModel> volmod;
  bool reg = ftvol->isRegularized(true);
  bool pattern_split = false; //true;
  if (!reg)
    {
      vector<SurfaceModel*> modified_adjacent;
      vector<shared_ptr<ftVolume> > reg_vols;
      try {
	reg_vols = 
	  ftvol->replaceWithRegVolumes(degree, modified_adjacent,
				       false, split_mode, pattern_split,
				       true);
      }
      catch (...)
	{
	  std::cout << "Failed createing trivariate block structured model" << std::endl;
	  return 1;
	}

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
  int nmb_vols0 = volmod->nmbEntities();
  for (int kr=0; kr<nmb_vols0; ++kr)
    {
      shared_ptr<ftVolume> curr_vol = volmod->getBody(kr);
      bool bd_trim = curr_vol->isBoundaryTrimmed();
      bool iso_trim = curr_vol->isIsoTrimmed();
      bool reg = curr_vol->isRegularized(true);

      vector<ftVolume*> ng0;
      curr_vol->getAdjacentBodies(ng0);

      std::cout << "Volume nr " << kr << ": " << bd_trim;
      std::cout << " " << iso_trim << " " << reg;
      std::cout << ", no of neighbours: " << ng0.size() << std::endl;
    }
	  
  // std::ofstream out_file("volmodel3.g22");
  // VolumeModelFileHandler filehandler;
  // filehandler.writeStart(out_file);
  // filehandler.writeHeader("Test VolumeModel", out_file);
  // filehandler.writeVolumeModel(*volmod, out_file);
  // filehandler.writeEnd(out_file);

  // VolumeModelFileHandler filehandler2;
  // shared_ptr<VolumeModel> volmod2 = filehandler2.readVolumeModel("volmodel3.g22");
  // std::cout << "Number of volumes: " << volmod2->nmbEntities() << std::endl;

  std::ofstream of6("output_volumes.g2");
  std::ofstream ofpar("output_par_volumes.g2");
  int nmb_vols = volmod->nmbEntities();
  for (int kr=0; kr<nmb_vols; ++kr)
    {
      shared_ptr<ftVolume> curr_vol = volmod->getBody(kr);
      bool bd_trim = curr_vol->isBoundaryTrimmed();
      bool iso_trim = curr_vol->isIsoTrimmed();
      bool reg = curr_vol->isRegularized(true);

      if (block_par && reg)
	{
	  int bd_cond[6][2];
	  shared_ptr<ParamVolume> reg_vol = 
	    curr_vol->getRegParVol(degree, bd_cond, true);
	  if (reg_vol.get())
	    {
	      reg_vol->writeStandardHeader(ofpar);
	      reg_vol->write(ofpar);
	    }
	}
      vector<ftVolume*> ng0;
      curr_vol->getAdjacentBodies(ng0);

      std::cout << "Volume nr " << kr << ": " << bd_trim;
      std::cout << " " << iso_trim << " " << reg;
      std::cout << ", no of neighbours: " << ng0.size() << std::endl;

      std::ofstream of7("Curr_vol.g2");
      shared_ptr<SurfaceModel> mod = curr_vol->getOuterShell();
      vector<shared_ptr<Vertex> > vxs;
      mod->getAllVertices(vxs);
      nmb = mod->nmbEntities();
      for (ki=0; ki<nmb; ++ki)
	{
	  shared_ptr<ParamSurface> sf = mod->getSurface(ki);
	  sf->writeStandardHeader(of7);
	  sf->write(of7);
	}
      of7 << std::endl << "400 1 0 4 255 0 0 255" << std::endl;
      of7 << vxs.size() << std::endl;
      for (ki=0; ki<(int)vxs.size(); ++ki)
	of7 << vxs[ki]->getVertexPoint() << std::endl;

      if (reg)
	{
	  vector<ftVolume*> ng1;
	  curr_vol->getAdjacentBodies(ng1);
	  std::cout << "Number of neighbours before untrim: " << ng1.size() << std::endl;
	  curr_vol->untrimRegular(degree, true);
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

  if (file_type_out == 1)
    {
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
  else
    {
      VolumeModelFileHandler filehandler_out;
      filehandler_out.writeStart(outfile);
      filehandler_out.writeHeader("Block structured volume model", outfile);
      filehandler_out.writeVolumeModel(*volmod, outfile);
      filehandler_out.writeEnd(outfile);
    }
}
