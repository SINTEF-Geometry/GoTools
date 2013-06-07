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

#include <iostream>
#include <fstream>
#include <cstdio>
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/igeslib/IGESconverter.h"

using std::cin;
using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using std::ofstream;
using std::vector;

using namespace Go;

int main(int argc, char* argv[])
{
    ifstream infile(argv[1]);
    ALWAYS_ERROR_IF(infile.bad(), "Wrong filename or corrupted file.");
    ALWAYS_ERROR_IF(argc != 5,
		    "Arguments: Infile Offset-dist Convertedfile Outfile");


    IGESconverter iges_converter;
    iges_converter.readdisp(infile);
    vector<shared_ptr<GeomObject> > geom_obj = iges_converter.getGoGeom();
    vector<shared_ptr<SplineSurface> > sfs;
    for (size_t ki = 0; ki < geom_obj.size(); ++ki)
	if (geom_obj[ki]->instanceType() == Class_SplineSurface)
	    sfs.push_back(dynamic_pointer_cast<SplineSurface>(geom_obj[ki]));

    ASSERT(sfs.size() == 1);

    ofstream outfile(argv[2]); // Write the input surface to output
			       // using Go-format.
    sfs[0]->writeStandardHeader(outfile);
    sfs[0]->write(outfile);


    // We convert the spline surface to SISL.
    SISLSurf* sisl_sf = GoSurf2SISL(*sfs[0]);

    int kstat = 0;          /* Local status flag */
    int index;

    SISLSurf *offset_sf = NULL; // The offset surface.

    double aoffset;
    double aepsge;
    double amax;
    int idim;

    /*
     * Fetch command line parameters
     * -----------------------------
     */

    aoffset = atof(argv[3]); // = 1.0; //atof(argv[2]);
    aepsge  = 1.0e-02; //atof(argv[3]);
    amax = 0; // atof(argv[4]);
    idim = sfs[0]->dimension(); // atoi(argv[5]);
    // fileout = argv[6];


    /*
     * Perform offset
     * --------------
     */
    s1365(sisl_sf, aoffset, aepsge, amax, idim, &offset_sf, &kstat);

    printf("s1365 : kstat=%d \n", kstat);
    printf("s1365 : in1=%d in2=%d ik1=%d ik2=%d \n", offset_sf->in1,
	   offset_sf->in2, offset_sf->ik1, offset_sf->ik2);
    printf("s1365 : u_startpar=%lf u_endpar=%lf \n",
	   offset_sf->et1[offset_sf->ik1 - 1],
	   offset_sf->et1[offset_sf->in1]);
    printf("s1365 : v_startpar=%lf v_endpar=%lf \n",
	   offset_sf->et2[offset_sf->ik2 - 1],
	   offset_sf->et2[offset_sf->in2]);
    printf("s1365 : (u_min, v_min)-corner=%lf %lf %lf \n",
	   offset_sf->ecoef[0], offset_sf->ecoef[1], offset_sf->ecoef[2]);
    index = (offset_sf->in1 - 1)*(offset_sf->idim);
    printf("s1365 : (u_max, v_min)-corner=%lf %lf %lf \n",
	   offset_sf->ecoef[index], offset_sf->ecoef[index + 1],
	   offset_sf->ecoef[index + 2]);
    index = (offset_sf->in2 - 1)*(offset_sf->in1)*(offset_sf->idim);
    printf("s1365 : (u_min, v_max)-corner=%lf %lf %lf \n",
	   offset_sf->ecoef[index], offset_sf->ecoef[index + 1],
	   offset_sf->ecoef[index + 2]);
    index = ((offset_sf->in2)*(offset_sf->in1) - 1)*(offset_sf->idim);
    printf("s1365 : (u_max, v_max)-corner=%lf %lf %lf \n",
	   offset_sf->ecoef[index], offset_sf->ecoef[index + 1],
	   offset_sf->ecoef[index + 2]);



    SplineSurface* offset_go_sf = SISLSurf2Go(offset_sf);

    ofstream outfile2(argv[4]); // Write the offset surface to output on Go-format.
    offset_go_sf->writeStandardHeader(outfile2);
    offset_go_sf->write(outfile2);

}
