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

#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/qualitymodule/FaceSetQuality.h"
#include "GoTools/geometry/PointOnCurve.h"
#include <fstream>

using std::vector;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 2) {
    std::cout << "Input parameters : Input file on g2 format" << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

    double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  shared_ptr<CompositeModel> model = shared_ptr<CompositeModel>(factory.createFromG2(file1));
  shared_ptr<SurfaceModel> sfmodel = dynamic_pointer_cast<SurfaceModel,CompositeModel>(model);

  FaceSetQuality quality(gap, kink, approxtol);
  quality.attach(sfmodel);

  vector<shared_ptr<PointOnCurve> > sing_points;
  vector<pair<shared_ptr<PointOnCurve>, shared_ptr<PointOnCurve> > > sing_curves;
  quality.vanishingCurveTangent(sing_points, sing_curves);

  std::cout << "Number of singular points: " << sing_points.size() << std::endl;
  std::cout << "Number of singular curves: " << sing_curves.size() << std::endl;

  std::ofstream out_file("crv_singularity.g2");
  size_t ki;
  if (sing_points.size() > 0)
  {
      out_file << "400 1 0 4 255 0 0 255" << std::endl;
      out_file << sing_points.size() << std::endl;
      for (ki=0; ki<sing_points.size(); ki++)
      {
	  Point pnt = sing_points[ki]->getPos();
	  out_file << pnt[0] << "  ";
	  out_file << pnt[1] << "  ";
	  out_file << pnt[2] << "  " << std::endl;
	  out_file << std::endl;
      }
  }
      
  for (ki=0; ki<sing_curves.size(); ki++)
  {
      shared_ptr<ParamCurve> crv = 
	  shared_ptr<ParamCurve>(sing_curves[ki].first->getCurve()->subCurve(sing_curves[ki].first->getPar(),
								 sing_curves[ki].second->getPar()));
      crv->writeStandardHeader(out_file);
      crv->write(out_file);
  }
      
  vector<shared_ptr<PointOnCurve> > sing_points2;
  vector<pair<shared_ptr<PointOnCurve>, shared_ptr<PointOnCurve> > > sing_curves2;
  quality.vanishingCurveTangent(sing_points2, sing_curves2);

  std::ofstream out_file2("crv_singularity2.g2");
  if (sing_points.size() > 0)
  {
      out_file2 << "400 1 0 4 255 0 0 255" << std::endl;
      out_file2 << sing_points2.size() << std::endl;
      for (ki=0; ki<sing_points2.size(); ki++)
      {
	  Point pnt = sing_points2[ki]->getPos();
	  out_file2 << pnt[0] << "  ";
	  out_file2 << pnt[1] << "  ";
	  out_file2 << pnt[2] << "  " << std::endl;
	  out_file2 << std::endl;
      }
  }
      
  for (ki=0; ki<sing_curves2.size(); ki++)
  {
      shared_ptr<ParamCurve> crv = 
	  shared_ptr<ParamCurve>(sing_curves2[ki].first->getCurve()->subCurve(sing_curves2[ki].first->getPar(),
								 sing_curves2[ki].second->getPar()));
      crv->writeStandardHeader(out_file2);
      crv->write(out_file2);
  }
      


}
