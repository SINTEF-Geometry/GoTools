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
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/ParamVolume.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/trivariatemodel/ftVolumeTools.h"
#include "GoTools/trivariatemodel/VolumeModelFileHandler.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include <fstream>

using namespace Go;
using std::vector;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;

//===========================================================================
//
/// Description: 
///
///              
///  
/// 
///
/// Input/Output: 
///
///               
/// 
/// Note:       
///
///             
///
//   
//===========================================================================

int main( int argc, char* argv[] )
{
  if (argc != 5 && argc != 6 && argc != 7)
    {
      cout << "Usage: " << "<infile> " << " <file type> " << "min nmb knot" << "distribute knots (0/1) " << "(start index)" << "(stop index)"<< endl;
      exit(-1);
    }
  std::string infile(argv[1]);
  int file_type = atoi(argv[2]);
  int min_nmb = atoi(argv[3]);
  int distribute = atoi(argv[4]);
  int start_ix = 0;
  if (argc >= 6)
    start_ix = atoi(argv[5]);
  int stop_ix = -1;
  if (argc == 7)
    stop_ix = atoi(argv[6]);

  shared_ptr<ftVolume> curr_vol;
  int degree = 3;
  if (file_type == 1)
    {
      // The tolerances must be set according to the properties of the model.
      // The neighbour tolerance must be smaller than the smallest entity in the
      // model, but larger than the largest gap.
      // The gap tolerance must be smaller than the neighbour tolerance
      double gap = 0.0001; 
      double neighbour = 0.001; 
      double kink = 0.01;
      double approxtol = 0.001;

      CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

      ifstream is(infile);
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

 
      // Select the first volume and pass through all elements and check if
      // they intersect the non-boundary trimming surface
      curr_vol = shared_ptr<ftVolume>(new ftVolume(sfmodel));

      std::ofstream out_file("volmodel.g22");
      VolumeModelFileHandler filehandler;
      filehandler.writeStart(out_file);
      filehandler.writeHeader("Test ftVolume", out_file);
      filehandler.writeVolume(curr_vol, out_file);
      filehandler.writeEnd(out_file);
    }
  else
    {
      VolumeModelFileHandler filehandler;
      curr_vol = filehandler.readVolume(infile.c_str());
    }

  // Add new knot lines if necessary
  shared_ptr<SplineVolume> curr_under = 
    dynamic_pointer_cast<SplineVolume>(curr_vol->getVolume());
  double parspan[3];
  double mean_span = 0.0;
  for (int dir=0; dir<3; ++dir)
    {
      parspan[dir] = curr_under->endparam(dir) - curr_under->startparam(dir);
      mean_span += parspan[dir];
    }
  mean_span /= 3.0;

  for (int dir=0; dir<3; ++dir)
    {
      int curr_min = 
	(distribute) ? (int)(parspan[dir]*min_nmb/mean_span) : min_nmb;
      curr_min = std::min(4*min_nmb, std::max(2, curr_min));
      int num_el = curr_under->numElem(dir);
      if (num_el < curr_min)
	{
	  double tstart =  curr_under->startparam(dir);
	  double tend =  curr_under->endparam(dir);
	  double del1 = (tend-tstart)/(double)curr_min;
	  vector<double> knots;
	  curr_under->basis(dir).knotsSimple(knots);
	  vector<double> newknots;
	  for (size_t kj=1; kj<knots.size(); ++kj)
	    {
	      double del2 = (knots[kj] - knots[kj-1]);
	      if (del2 < del1)
		continue;
	      int nmb = std::max(1, (int)(del2/del1));
	      double del3 = del2/(double)(nmb+1);
	      int kr;
	      double par;
	      for (kr=0, par=knots[kj-1]+del3; kr<nmb; ++kr, par+=del3)
		newknots.push_back(par);
	    }
	  
	  curr_under->insertKnot(dir, newknots);
	}
    }


  shared_ptr<SplineVolume> under = 
    dynamic_pointer_cast<SplineVolume>(curr_vol->getVolume());
  double gap = curr_vol->getTolerances().gap;
  
  int nmb_elem = under->numElem();
  std::cout << "No of elements: " << nmb_elem << std::endl;

  std::ofstream of5("tmp5.g2");
  std::ofstream of6("tmp6.g2");
  std::ofstream of7("tmp7.g2");
  if (stop_ix < 0)
    stop_ix = nmb_elem;
  for (int ki=start_ix; ki<stop_ix; ++ki)
    {
      std::ofstream of8("tmp8.g2");

      // Fetch parameter values surrounding the specified element
      double elem_par[6];
      under->getElementBdPar(ki, elem_par);
	  
      // Create an ftVolume entity corresponding to the element. 
      // First create underlying SplineVolume
      shared_ptr<ParamVolume> elem_vol(under->subVolume(elem_par[0], elem_par[2],
							elem_par[4], elem_par[1],
							elem_par[3], elem_par[5]));
      elem_vol->writeStandardHeader(of8);
      elem_vol->write(of8);


      int elem_stat = curr_vol->ElementBoundaryStatus(ki);
      std::cout << "Boundary status, element " << ki << ": " << elem_stat << std::endl;

      int nmb_par=0;
      int nmb_geom=0;
      int nmb_blocks=0;
      if (elem_stat == 1)
	{
	  // Element intersects trimming surface
	  // Split element with trimming shell
	  vector<shared_ptr<ftVolume> > sub_elem;
	  vector<int> is_inside; // Equal 1 if the sub element is inside
	  // the trimmed volume
	  curr_vol->splitElementByTrimSfs(ki, gap, sub_elem, is_inside);

	  std::ofstream of4("tmp4.g2");
	  for (size_t kj=0; kj<sub_elem.size(); ++kj)
	    {
	      shared_ptr<SurfaceModel> mod = sub_elem[kj]->getOuterShell();
	      int nmb = mod->nmbEntities();
	      for (int kr=0; kr<nmb; ++kr)
		{
		  shared_ptr<ParamSurface> sf = mod->getSurface(kr);
		  sf->writeStandardHeader(of4);
		  sf->write(of4);
		}
	    }
	  int stop_break = 1;

	  // if (sub_elem.size() < 2)
	  //   continue;

 	  // Check if the remaining element is hexagonal
	  for (size_t kj=0; kj<sub_elem.size(); ++kj)
	    {
	      if (!is_inside[kj])
		continue;

	      std::ofstream of9("tmp9.g2");
	      shared_ptr<SurfaceModel> mod = sub_elem[kj]->getOuterShell();
	      int nmb = mod->nmbEntities();
	      for (int kr=0; kr<nmb; ++kr)
		{
		  shared_ptr<ParamSurface> sf = mod->getSurface(kr);
		  sf->writeStandardHeader(of9);
		  sf->write(of9);
		}

	      if (sub_elem[kj]->getOuterShell()->nmbBoundaries() > 0)
		{
		  std::cout << "Open shell. Check" << std::endl;
		  vector<shared_ptr<ftEdge> > bd_edgs = 
		    sub_elem[kj]->getOuterShell()->getBoundaryEdges();
		  std::ofstream of_bd("tmp9_bd.g2");
		  for (size_t kr=0; kr<bd_edgs.size(); ++kr)
		    {
		      shared_ptr<ParamCurve> tmp_cv = bd_edgs[kr]->geomCurve();
		      shared_ptr<ParamCurve> tmp_cv2(tmp_cv->geometryCurve());
		      tmp_cv2->writeStandardHeader(of_bd);
		      tmp_cv2->write(of_bd);
		    }
		  continue;  
		}

	      bool regular = sub_elem[kj]->isRegularized(true);
	      //std::cout << "Sub element nr " << kj+1 << ": " << regular << std::endl;
	      if (regular)
		{
		  nmb_blocks++;
		  // if (false)
		  //   {
		  // Create non-trimmed parameter element
		  int bd_cond[6][2];
		  shared_ptr<ParamVolume> reg_vol = 
		    sub_elem[kj]->getRegParVol(degree, bd_cond, true);
		  if (reg_vol.get())
		    {
		      std::cout << "Boundary conditions: ";
		      for (int ka=0; ka<6; ++ka)
			std::cout << bd_cond[ka][0] << " " << bd_cond[ka][1] << " ";
		      std::cout << std::endl;

		      reg_vol->writeStandardHeader(of5);
		      reg_vol->write(of5);
		      nmb_par++;
		    }

		  // Create non-trimmed element
		  bool done = sub_elem[kj]->untrimRegular(degree, true);
		  if (done)
		    {
		      shared_ptr<ParamVolume> tmp_vol = sub_elem[kj]->getVolume();
		      tmp_vol->writeStandardHeader(of6);
		      tmp_vol->write(of6);
		      nmb_geom++;
		    }
		    // }
		}
	      else
		{
		  std::ofstream pre_block("pre_block.g22");
		  VolumeModelFileHandler filewrite0;
		  filewrite0.writeStart(pre_block);
		  filewrite0.writeHeader("Irregular volume", pre_block);
		  filewrite0.writeVolume(sub_elem[kj], pre_block);
		  filewrite0.writeEnd(pre_block);
		  
		  // Split in concave edges
		  vector<shared_ptr<ftVolume> > split_elem; 
		  bool do_split = true;
		  if (do_split)
		    {
		      split_elem = sub_elem[kj]->splitConcaveVol(degree, true);
		      std::cout << "split concave: " << split_elem.size() << std::endl;
		      std::ofstream pre_block1("pre_block1.g22");
		      VolumeModelFileHandler filewrite1;
		      filewrite1.writeStart(pre_block1);
		      filewrite1.writeHeader("Irregular volume after pre split", pre_block1);
		      filewrite1.writeVolumes(sub_elem, pre_block1);
		      filewrite1.writeEnd(pre_block1);
		    }
		  if (split_elem.size() == 0)
		    split_elem.push_back(sub_elem[kj]);

		  vector<shared_ptr<ftVolume> > blocks;
		  bool failed = false;

		  for (size_t kr=0; kr<split_elem.size(); ++kr)
		    {
		      std::ofstream pre_block2("pre_block2.g22");
		      VolumeModelFileHandler file2;
		      file2.writeStart(pre_block2);
		      file2.writeHeader("Pre block structureing", pre_block2);
		      file2.writeVolume(split_elem[kr], pre_block2);
		      file2.writeEnd(pre_block2);
		  
		      // Block structuring
		      vector<SurfaceModel*> modified_adjacent;
		      bool pattern_split = false;
		      int split_mode = 1;
		      vector<shared_ptr<ftVolume> > blocks0;
		      try {
			blocks0 = 
			  split_elem[kr]->replaceWithRegVolumes(degree, modified_adjacent,
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
			blocks.push_back(split_elem[kr]);
		    }

		  if (blocks.size() == 0)
		    failed = true;

		  std::ofstream post_block("post_block.g22");
		  VolumeModelFileHandler filewrite2;
		  filewrite2.writeStart(post_block);
		  filewrite2.writeHeader("Block structured volume", post_block);
		  filewrite2.writeVolumes(blocks, post_block);
		  filewrite2.writeEnd(post_block);

		  nmb_blocks += (int)blocks.size();
		  for (size_t kr=0; kr<blocks.size(); ++kr)
		    {
		      regular = blocks[kr]->isRegularized(true);
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
			      std::cout << "Boundary conditions: ";
			      for (int ka=0; ka<6; ++ka)
				std::cout << bd_cond[ka][0] << " " << bd_cond[ka][1] << " ";
			      std::cout << std::endl;

			      reg_vol->writeStandardHeader(of5);
			      reg_vol->write(of5);
			      nmb_par++;
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
			      nmb_geom++;
			    }
			  else
			    failed = true;
			}
		      else
			failed = true;
		    }
		  
		  if (failed)
		    {
		      std::ofstream of4_2("tmp4_2.g2");
		      shared_ptr<SurfaceModel> mod = sub_elem[kj]->getOuterShell();
		      int nmb = mod->nmbEntities();
		      for (int kr=0; kr<nmb; ++kr)
			{
			  shared_ptr<ParamSurface> sf = mod->getSurface(kr);
			  sf->writeStandardHeader(of4_2);
			  sf->write(of4_2);
			  sf->writeStandardHeader(of7);
			  sf->write(of7);
			}

		      std::ofstream of_mod("tmp_mod.g22");
		      VolumeModelFileHandler filewrite;
		      filewrite.writeStart(of_mod);
		      filewrite.writeHeader("Irregular volume", of_mod);
		      filewrite.writeVolume(sub_elem[kj], of_mod);
		      filewrite.writeEnd(of_mod);
		  

		      std::cout << "Number of surfaces: " << sub_elem[kj]->getOuterShell()->nmbEntities() << std::endl;
		    }
		}
	    }
	}
      else if (elem_stat == 2)
	{
	  elem_vol->writeStandardHeader(of6);
	  elem_vol->write(of6);
	  nmb_blocks++;
	  nmb_geom++;
	  nmb_par++;
	}
      if (nmb_par != nmb_blocks /*|| nmb_geom != nmb_blocks*/)
	std::cout << "Nmb blocks: " << nmb_blocks << ", nmb geom: " << nmb_geom << ", nmb par: " << nmb_par << std::endl;
    }
}
