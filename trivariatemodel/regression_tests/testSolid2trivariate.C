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

#define BOOST_TEST_MODULE trivariatemodel/solid2trivariate
#include <boost/test/included/unit_test.hpp>

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
using std::string;


struct Config {
public:
    Config()
    {
        datadir = "../compositemodel/data/"; // Relative to build/compositemodel

        infiles.push_back("b25.g2");
        numfaces.push_back(2);

        infiles.push_back("b26.g2");
        numfaces.push_back(4);

        gap = 0.001; // 0.001;
        neighbour = 0.01; // 0.01;
        kink = 0.01;
        approxtol = 0.01;

    }

public:
    string datadir;
    vector<string> infiles;
    vector<int> numfaces;
    double gap;
    double neighbour;
    double kink;
    double approxtol;
};


BOOST_FIXTURE_TEST_CASE(solid2trivariate, Config)
{
    int nfiles = infiles.size();
    for (int i = 0; i < nfiles; ++i) {

        string infilestr = datadir + infiles[i];

        ifstream infile(infilestr.c_str());
        BOOST_CHECK_MESSAGE(infile.good(), "Bad or no input filename");

        ofstream outfile("outfile.g2");

        CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

        CompositeModel *model = factory.createFromG2(infile);

        shared_ptr<SurfaceModel> sfmodel = 
            shared_ptr<SurfaceModel>(dynamic_cast<SurfaceModel*>(model));
        if (!sfmodel.get())
            exit(-1);

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
	int degree = 3;
        shared_ptr<VolumeModel> volmod;
        bool reg = ftvol->isRegularized();
        if (!reg)
        {
	  vector<SurfaceModel*> modified_adjacent;  // Dummy vector
             vector<shared_ptr<ftVolume> > reg_vols = 
	       ftvol->replaceWithRegVolumes(degree, modified_adjacent,
					    false);

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
}
