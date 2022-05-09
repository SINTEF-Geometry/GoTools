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

#include "GoTools/lrsplines2D/TrimUtils.h"
#include "GoTools/lrsplines2D/TrimSurface.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/utils/Point.h"
#include "GoTools/utils/Array.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/RectDomain.h"
#include "GoTools/geometry/LoopUtils.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/lrsplines2D/TrimCrvUtils.h"
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

using namespace Go;

using std::cout;
using std::endl;
using std::vector;
using std::string;

//#define DEBUG

int main(int argc, char* argv[])
{
  if (argc != 5)
    {
      cout << "Usage: " << argv[0] << " sf.g2 point_cloud.g2 tightness bd_sf.g2" << endl;
      return -1;
    }

    
  std::ifstream filein_sf(argv[1]);
  std::ifstream filein_points(argv[2]);
  int tightness = atoi(argv[3]);
  std::ofstream fileout_bd_sf(argv[4]);

  GoTools go_tools;
  go_tools.init();
  Registrator<LRSplineSurface> r293;

  ObjectHeader header;
  header.read(filein_sf);
  shared_ptr<LRSplineSurface> sf;
  shared_ptr<GeomObject> geom_obj(Factory::createObject(header.classType()));
  geom_obj->read(filein_sf);
  sf = dynamic_pointer_cast<LRSplineSurface, GeomObject>(geom_obj);
  if (!sf.get())
    {
      std::cerr << "Input file contains no LR B-spline surface" << std::endl;
      exit(-1);
    }

  header.read(filein_points);
  PointCloud3D points;
  points.read(filein_points);

  // Check for correspondance between the point set and the parameter
  // domain of the surface

  BoundingBox box = points.boundingBox();
  Point low = box.low();
  Point high = box.high();

  RectDomain dom = sf->containingDomain();
  double umin = dom.umin();
  double umax = dom.umax();
  double vmin = dom.vmin();
  double vmax = dom.vmax(); 
  
  if (umin > low[0] || umax < high[0] ||  vmin > low[1] || vmax < high[1])
    {
      std::cout << " Point cloud not corresponding to surface domain. Exiting " << std::endl;
      exit(-1);
    }

  // Translate to origo
  Vector3D vec(-0.5*(umin+umax), -0.5*(vmin+vmax), 0.0);
  points.translate(vec);

#ifdef DEBUG
  std::ofstream of1("translated_cloud.g2");
  points.writeStandardHeader(of1);
  points.write(of1);

  std::ofstream ofp("projected_cloud.g2");
  points.writeStandardHeader(ofp);
  int num = points.numPoints();
  ofp << num << std::endl;
  for (int ka=0; ka<num; ++ka)
    {
      Vector3D curr = points.point(ka);
      ofp << curr[0] << " " << curr[1] << " " << 0.0<< std::endl;
    }
#endif
  
  // Update parameter domain
  Point tmp_vec(vec.begin(), vec.end());
  TrimCrvUtils::translateSurfaceDomain(sf.get(), tmp_vec);

  if (sf->dimension() == 3)
    {
      sf->translate(tmp_vec);
    }

  // Set parameters for computations of trimming sequence
  int max_rec;
  int nmb_div;
  if (tightness <= 2)
    {
      max_rec = 1;
      nmb_div = (tightness == 2) ? 20 : 15;
    }
  else if (tightness <= 5)
    {
      max_rec = 2;
      nmb_div = (tightness == 3) ? 8 : ((tightness == 4) ? 12 : 15);
    }
  else
    {
      max_rec = 3;
      nmb_div = (tightness == 6) ? 10 : ((tightness == 7) ? 12 : 15);
    }

  // Compute trimming seqence
  vector<double> points2(points.rawData(), 
			 points.rawData()+points.dimension()*points.numPoints());
  RectDomain dom2 = sf->containingDomain();
  double domain[4];
  domain[0] = dom2.umin();
  domain[1] = dom2.umax();
  domain[2] = dom2.vmin();
  domain[3] = dom2.vmax();
  vector<vector<vector<double> > > seqs;
  TrimUtils trimutil(&points2[0], points.numPoints(), 1);
  trimutil.computeAllTrimSeqs(max_rec, nmb_div, seqs);
  
  double udel, vdel;
  trimutil.getDomainLengths(udel, vdel);
  std::cout << "Minimum sub domains, diag = "<< trimutil.getDomainDiag();
  std::cout << ", udel = " << udel << ", vdel = " << vdel << std::endl;

#ifdef DEBUG
  std::ofstream ofsf("trimsfs.g2");
#endif
  double eps = std::max(udel, vdel);
  bool isotrim[4];
  isotrim[0] = isotrim[1] = isotrim[2] = isotrim[3] = false;
  for (size_t kh=0; kh<seqs.size(); ++kh)
    {
      // Compute domain
      double umin, umax, vmin, vmax;
      umin = umax = seqs[kh][0][0];
      vmin = vmax = seqs[kh][0][1];
      vector<vector<vector<double> > > loop_seqs(seqs[kh].size());
      for (size_t kr=0; kr<seqs[kh].size(); ++kr)
	{
	  loop_seqs[kr].push_back(seqs[kh][kr]);
	  for (size_t kj=0; kj<seqs[kh][kr].size(); kj+=2)
	    {
	      umin = std::min(umin, seqs[kh][kr][kj]);
	      umax = std::max(umax, seqs[kh][kr][kj]);
	      vmin = std::min(vmin, seqs[kh][kr][kj+1]);
	      vmax = std::max(vmax, seqs[kh][kr][kj+1]);
	    }
	}
      double udel = 0.1*(umax - umin);
      double vdel = 0.1*(vmax - vmin);
      umin = std::max(sf->startparam_u(), umin-udel);
      umax = std::min(sf->endparam_u(), umax+udel);
      vmin = std::max(sf->startparam_v(), vmin-vdel);
      vmax = std::min(sf->endparam_v(), vmax+vdel);
      shared_ptr<BoundedSurface> trim_surf;
      shared_ptr<ParamSurface> sf2 = shared_ptr<ParamSurface>(sf->subSurface(umin, vmin,
									     umax, vmax,
									     1.0e-2));
      bool found = TrimSurface::defineBdSurface(sf2, domain, isotrim, eps,
						loop_seqs, trim_surf);
#ifdef DEBUG
       trim_surf->writeStandardHeader(ofsf);
      trim_surf->write(ofsf);
#endif
     if (found)
	{
	  // Translate back
	  trim_surf->setParameterDomain(umin-vec[0], umax-vec[0], vmin-vec[1], vmax-vec[1]);
	}

      trim_surf->writeStandardHeader(fileout_bd_sf);
      trim_surf->write(fileout_bd_sf);
      int stop_break = 1;
    }
  int stop_break2 = 1;
}
