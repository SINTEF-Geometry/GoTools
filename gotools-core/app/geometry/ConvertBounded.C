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
#include <sstream>
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/GoTools.h"


using namespace std;
using namespace Go;


string col(int nmb)
{
  if (nmb == 0)
    return "255 0 255 255";
  else if (nmb == 1)
    return "0 255 0 255";
  else if ((nmb % 2) == 0)
    return "255 0 0 255";
  else
    return "0 0 255 255";
}


int main(int argc, char** argv)
{
  GoTools::init();

  if (argc != 5) {
    cout << "Usage:  " << argv[0] << " nmbSurfcs infile geofilecore parfilecore" << endl;
    return 1;
  }

  int nmb_surfs = atoi(argv[1]);
  vector<BoundedSurface*> surfs;

  ifstream ins(argv[2]);

  ObjectHeader header;
  for (int i = 0; i < nmb_surfs; ++i)
    {
      BoundedSurface* bs = new BoundedSurface();
      ins >> header >> *bs;
      surfs.push_back(bs);
    }
  ins.close();

  for (int i = 0; i < nmb_surfs; ++i)
    {
      shared_ptr<ParamSurface> surf = surfs[i]->underlyingSurface();

      shared_ptr<Cylinder> surf_cyl = dynamic_pointer_cast<Cylinder>(surf);
      shared_ptr<Sphere> surf_sph = dynamic_pointer_cast<Sphere>(surf);

      /*
      if (surf_cyl)
	{
	  CurveBoundedDomain cbd = surfs[i]->parameterDomain();
	  RectDomain rd = cbd.containingDomain();
	  if (rd.umin() < -1.0 || rd.umax() > 8.0)
	    surf_cyl->setParamBoundsU(rd.umin(), rd.umax());
	  if (rd.vmin() < -1.0 || rd.vmax() > 8.0)
	    surf_cyl->setParamBoundsV(rd.vmin(), rd.vmax());
	}
      */

      /*
      if (i == 5)
	surf_cyl->setParamBoundsV(-800.0, -500.0);
      if (i == 6)
	surf_cyl->setParamBoundsU(200.0, 600.0);
      if (i == 8)
	surf_cyl->setParamBoundsU(-700.0, 200.0);
      if (i == 11)
	surf_cyl->setParamBoundsV(-2000.0, -400.0);
      */

      stringstream sstm_geo;
      sstm_geo << argv[3] << ((i < 10) ? "0" : "") << i << ".g2";
      ofstream outs_geo(sstm_geo.str());
      outs_geo << surf->instanceType() << " 1 0 4 255 255 0 255" << endl;
      surf->write(outs_geo);

      stringstream sstm_par;
      sstm_par << argv[4] << ((i < 10) ? "0" : "") << i << ".g2";
      ofstream outs_par(sstm_par.str());
      outs_par << "200 1 0 4 255 255 0 255" << endl << "3 0" << endl;
      outs_par << "2 2" << endl << "0 0 1 1" << endl;
      outs_par << "2 2" << endl << "0 0 1 1" << endl;
      RectDomain dom = surf->containingDomain();
      outs_par << dom.umin() << " " << dom.vmin() << " 0" << endl;
      outs_par << dom.umax() << " " << dom.vmin() << " 0" << endl;
      outs_par << dom.umin() << " " << dom.vmax() << " 0" << endl;
      outs_par << dom.umax() << " " << dom.vmax() << " 0" << endl;
      outs_par << endl;

      if (surf_sph.get())
	{
	  outs_geo << "130 1 0 4 0 0 0 255" << endl << 3 << endl << surf_sph->getRadius() << endl;
	  outs_geo << surf_sph->getLocation() << endl;
	  Point x, y, z;
	  surf_sph->getCoordinateAxes(x, y, z);
	  outs_geo << y << endl;
	  outs_geo << z << endl;
	  outs_geo << "0 " << M_PI << endl << 0 << endl << endl;
	}

      if (surf_cyl.get())
	{
	  outs_geo << "120 1 0 4 0 0 0 255" << endl << 3 << endl;
	  Point x, y, z;
	  surf_cyl->getCoordinateAxes(x, y, z);
	  outs_geo << (surf_cyl->getLocation() + x * surf_cyl->getRadius()) << endl;
	  outs_geo << z << endl;
	  outs_geo << 1 << endl;
	  RectDomain syl_dom = surf_cyl->containingDomain();
	  if (i == 6 || i == 8)
	    outs_geo << dom.umin() << " " << dom.umax() << endl;
	  else
	    outs_geo << dom.vmin() << " " << dom.vmax() << endl;
	  outs_geo << 0 << endl;
	}

      vector<CurveLoop> loops = surfs[i]->allBoundaryLoops();
      for (int j = 0; j < loops.size(); ++j)
	{
	  int curve_cnt = 0;
	  vector<shared_ptr<ParamCurve> > parCurves = loops[j].getCurves();
	  for (int k = 0; k < parCurves.size(); ++k, ++curve_cnt)
	    {
	      shared_ptr<CurveOnSurface> cos = dynamic_pointer_cast<CurveOnSurface>(parCurves[k]);
	      shared_ptr<ParamCurve> geo_c = cos->spaceCurve();
	      outs_geo << geo_c->instanceType() << " 1 0 4 " << col(curve_cnt) << endl;
	      geo_c->write(outs_geo);

	      shared_ptr<ParamCurve> par_c = cos->parameterCurve();
	      if (par_c.get())
		{
		  shared_ptr<SplineCurve> spl_c = dynamic_pointer_cast<SplineCurve>(par_c);
		  shared_ptr<Line> line_c = dynamic_pointer_cast<Line>(par_c);
		  if (spl_c.get())
		    {
		      vector<double> newKnots;
		      bool second_coord = false;
		      for (vector<double>::const_iterator it = spl_c->coefs_begin(); it != spl_c->coefs_end(); ++it)
			{
			  newKnots.push_back(*it);
			  if (second_coord)
			    newKnots.push_back(0);
			  second_coord = !second_coord;
			}
		      SplineCurve spl_c_3d = SplineCurve(spl_c->basis(), newKnots.begin(), 3);
		      outs_par << "100 1 0 4 " << col(curve_cnt) << endl;
		      spl_c_3d.write(outs_par);
		    }
		  else if (line_c.get())
		    {
		      outs_par << "120 1 0 4 " << col(curve_cnt) << endl;
		      outs_par << 3 << endl;
		      outs_par << line_c->getPoint() << " 0" << endl;
		      outs_par << line_c->getDirection() << " 0" << endl;
		      outs_par << 1 << endl;
		      outs_par << line_c->startparam() << " " << line_c->endparam() << endl;
		      if (!line_c->isReversed())
			outs_par << "0" << endl;
		      else
			outs_par << "1" << endl;
		      outs_par << endl;
		    }
		}
	    }
	}
      outs_geo.close();
      outs_par.close();
    }
}

