//===========================================================================
//                                                                           
// File: testOffsetSurf.C                                                        
//                                                                           
// Created: Mon Dec 10 13:36:14 2001                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@math.sintef.no>
//                                                                           
// Revision: $Id: testOffsetSurf.C,v 1.6 2006-04-19 11:08:56 jbt Exp $
//                                                                           
// Description: Test program for offsetting of surface. Calls sisl, sh1929.
//                                                                           
//===========================================================================

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
using std::shared_ptr;
using std::dynamic_pointer_cast;

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
