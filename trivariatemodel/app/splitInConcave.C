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

#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/ModifyFaceSet.h"
#include "GoTools/trivariatemodel/VolumeModelFileHandler.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"

using namespace Go;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::vector;

// Test splitting of trimmed volumes in cancave edges in the context
// of trivariate block structuring

int main(int argc, char* argv[] )
{
  if (argc != 3 && argc != 4) {
      cout << "Usage: " << "<infile> <outfile> (<file type, 2 for g22>) " << endl;
      exit(-1);
    }

  std::string infile(argv[1]);

  ofstream out_file(argv[2]);

  int file_type = 1;
  if (argc > 3)
    file_type = atoi(argv[3]);

  double gap = 0.001; //0.001;
  double neighbour = 0.01; //0.01;
  double kink = 0.01;
  double approxtol = 0.001;
  int degree = 3;

  shared_ptr<ftVolume> ftvol;
  if (file_type == 2)
    {
      VolumeModelFileHandler fileread;
      ftvol = fileread.readVolume(infile.c_str());
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


      ftvol = shared_ptr<ftVolume>(new ftVolume(sfmodel));
    }

  vector<shared_ptr<ftVolume> > reg_vols;
  try {
    reg_vols = 
      ftvol->splitConcaveVol(degree, true);
  }
  catch (...)
    {
      std::cout << "Failed splitting model" << std::endl;
      return 1;
    }

  if (reg_vols.size() > 0)
    {
      VolumeModelFileHandler filewrite;
      filewrite.writeStart(out_file);
      filewrite.writeHeader("Divided trimmed volumes", out_file);
      filewrite.writeVolumes(reg_vols, out_file);
      filewrite.writeEnd(out_file);
    }

  vector<shared_ptr<ftVolume> > blocks;
  bool failed = false;

  std::ofstream of5("par_sub.g2");
  std::ofstream of6("geom_sub.g2");
  for (size_t kr=0; kr<reg_vols.size(); ++kr)
    {
      // Block structuring
      vector<SurfaceModel*> modified_adjacent;
      bool pattern_split = false;
      int split_mode = 1;
      vector<shared_ptr<ftVolume> > blocks0;
      try {
	blocks0 = 
	  reg_vols[kr]->replaceWithRegVolumes(degree, modified_adjacent,
						false, split_mode, 
						pattern_split, true);
      }
      catch (...)
	{
	  failed = true;
	}
      if (blocks0.size() > 0)
	blocks.insert(blocks.end(), blocks0.begin(), blocks0.end());
      else
	blocks.push_back(reg_vols[kr]);
    }

  if (blocks.size() == 0)
    failed = true;

  for (size_t kr=0; kr<blocks.size(); ++kr)
    {
      bool regular = blocks[kr]->isRegularized(true);
      if (regular)
	{
	  // if (false)
	  //   {
	  // Create non-trimmed parameter element
	  int bd_cond[6][2];
	  shared_ptr<ParamVolume> reg_vol = 
	    blocks[kr]->getRegParVol(degree, bd_cond, true);
	  if (reg_vol.get())
	    {
	      reg_vol->writeStandardHeader(of5);
	      reg_vol->write(of5);
	    }
	  else
	    failed = true;
	  // }
	  // Create non-trimmed element

	  std::ofstream of7("elem_sub.g2");
	  shared_ptr<SurfaceModel> mod = blocks[kr]->getOuterShell();
	  int nmb = mod->nmbEntities();
	  for (int kc=0; kc<nmb; ++kc)
	    {
	      shared_ptr<ParamSurface> sf = mod->getSurface(kc);
	      sf->writeStandardHeader(of7);
	      sf->write(of7);
	    }

	  blocks[kr]->untrimRegular(degree, true);
	  shared_ptr<ParamVolume> tmp_vol = blocks[kr]->getVolume();
	  if (tmp_vol.get())
	    {
	      tmp_vol->writeStandardHeader(of6);
	      tmp_vol->write(of6);
	    }
	  else
	    failed = true;
	}
      else
	failed = true;
    }
}
