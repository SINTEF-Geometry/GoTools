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
#include <string.h>
#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/trivariatemodel/CreateTrimVolume.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/trivariatemodel/VolumeModelFileHandler.h"

#define DEBUG

using namespace Go;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::string;

int compare(const char *str1, char str2[][8], int nmb)
{
  for (int ki=0; ki<nmb; ++ki)
    if (strcmp(str1, str2[ki]) == 0)
      return ki;
  return -1;
}

int main(int argc, char* argv[] )
{
  if (argc != 3)
    {
      cout << "Usage: " << "<infile (g2/g22)> <outfile (g22)>" << endl;
      exit(-1);
    }

  char* infile(argv[1]);
  ofstream outfile(argv[2]);

  // Check file type
  // Find file extension
  char *loc;
  char *last = 0;
  loc = strchr(infile, '.');
  while (loc != NULL)
    {
      last = loc;
      loc = strchr(loc+1, '.');
    }
  if (last == NULL)
    {
      std::cout << "Missing file extension of input file" << std::endl;
      exit(1);
    }
  char *input_type = last+1;

  // Check type of input file
  char keys[2][8] = {"g2", "g22"};
  int type_in;
  try {
    type_in = compare(input_type, keys, 2);
  }
  catch (...)
    {
      std::cout << "Invalide file extension of input file" << std::endl;
      exit(1);
    }
  if (type_in < 0)
    {
      std::cout << "Not a valid file type" << std::endl;
      exit(1);
    }

  double gap, neighbour, kink;
  shared_ptr<SurfaceModel> sfmodel;
  vector<shared_ptr<SurfaceModel> > voids;
  int material_id = -1;
  if (type_in == 1)
    {
      VolumeModelFileHandler fileread;
      shared_ptr<Body> body = fileread.readBody(infile);
      // if (!body.get())
      // 	body = fileread.readVolume(infile.c_str());
      if (!body.get())
	{
	  sfmodel = fileread.readShell(infile);
	  if (!sfmodel.get())
	    exit(1);
	}
      else
	{
	  int nmb = body->nmbOfShells();
	  if (nmb == 1)
	    sfmodel = body->getOuterShell();
	  else
	    {
	      vector<shared_ptr<SurfaceModel> > all_shells = 
		body->getAllShells();		
	      sfmodel = all_shells[0];
	      voids.insert(voids.end(), all_shells.begin()+1,
			   all_shells.end());
	    }
	  material_id = body->getMaterial();
	}
      tpTolerances top = sfmodel->getTolerances();
      gap = top.gap;
      neighbour = top.neighbour;
      kink = top.kink;
    }
  else
    {
      // The tolerances must be set according to the properties of the model.
      // The neighbour tolerance must be smaller than the smallest entity in the
      // model, but larger than the largest gap.
      // The gap tolerance must be smaller than the neighbour tolerance
      double gap = 0.0001; //0.0001;
      double neighbour = 0.001; //0.001;
      double kink = 0.01;
      double approxtol = 0.001;

      CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

      std::ifstream is(infile);
      CompositeModel *model = factory.createFromG2(is);

      sfmodel = 
	shared_ptr<SurfaceModel>(dynamic_cast<SurfaceModel*>(model));
      if (!sfmodel.get())
	{
	  std::cout << "No input model read" << std::endl;
	  exit(1);
	}
 
      if (sfmodel->nmbBoundaries() > 0)
	{
	  std::cout << "Not a brep solid. Consider increasing the neighbour tolerance" << std::endl;
	  exit(1);
	}
      
#ifdef DEBUG
      bool isOK = sfmodel->checkShellTopology();
      std::cout << "Shell topology: " << isOK << std::endl;
#endif
    }
  CreateTrimVolume trim(sfmodel, material_id);
  if (voids.size() > 0)
    trim.addVoids(voids);

  shared_ptr<ftVolume> vol;
  try {
    vol = trim.fetchRotationalTrimVol();
  }
  catch (...)
    {
    }
  if (!vol.get())
    {
      try {
	vol = trim.fetchOneTrimVol();
      }
      catch (...)
	{
	}
      }

  if (!vol.get())
    {
      try {
	vol = shared_ptr<ftVolume>(new ftVolume(sfmodel));
      }
      catch (...)
	{
	}
      }

#ifdef DEBUG
  std::ofstream of1("under_vol_final.g2");
  vol->getVolume()->writeStandardHeader(of1);
  vol->getVolume()->write(of1);
  std::ofstream of2("vol_shell.g2");
  int nmb_shell = vol->nmbOfShells();
  for (int ka=0; ka<nmb_shell; ++ka)
    {
      shared_ptr<SurfaceModel> shell = vol->getShell(ka);
      int nmb_face = shell->nmbEntities();
      for (int kb=0; kb<nmb_face; ++kb)
	{
	  shared_ptr<ParamSurface> surf = shell->getSurface(kb);
	  surf->writeStandardHeader(of2);
	  surf->write(of2);
	}
	}
#endif

  VolumeModelFileHandler filehandler;
  filehandler.writeStart(outfile);
  filehandler.writeHeader("Trimmed volume", outfile);
  filehandler.writeVolume(vol, outfile);
  filehandler.writeEnd(outfile);
  int stop_break = 1;
}
