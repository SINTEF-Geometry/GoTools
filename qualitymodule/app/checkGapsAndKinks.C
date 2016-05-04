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
#include <fstream>

using std::vector;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 2) {
    std::cout << "Input parameters : Input file on g2 format," << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  double gap = 0.01; //0.00001;
  double neighbour = 0.1;
  double kink = 0.001;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  shared_ptr<CompositeModel> model = shared_ptr<CompositeModel>(factory.createFromG2(file1));
  shared_ptr<SurfaceModel> sfmodel = dynamic_pointer_cast<SurfaceModel,CompositeModel>(model);

  FaceSetQuality quality(gap, kink, approxtol);
  quality.attach(sfmodel);

  vector<pair<ftEdge*, ftEdge*> > gaps;
  vector<pair<ftEdge*, ftEdge*> > kinks;
  quality.facePositionDiscontinuity(gaps);  
  quality.faceTangentDiscontinuity(kinks);

  std::cout << "Number of gaps: " << gaps.size() << std::endl;
  std::cout << "Number of kinks: " << kinks.size() << std::endl;

  std::ofstream out_file("gaps.g2");
  size_t ki;
  for (ki=0; ki<gaps.size(); ++ki)
  {
      shared_ptr<ParamCurve> crv1 = gaps[ki].first->geomCurve();
      shared_ptr<CurveOnSurface> sfcv = 
	dynamic_pointer_cast<CurveOnSurface,ParamCurve>(crv1);
      if (sfcv.get())
	{
	  sfcv->spaceCurve()->writeStandardHeader(out_file);
	  sfcv->spaceCurve()->write(out_file);
	}
      else
	{
	  crv1->writeStandardHeader(out_file);
	  crv1->write(out_file);
	}
  }

  std::ofstream out_file2("kinks.g2");
  for (ki=0; ki<kinks.size(); ++ki)
  {
      shared_ptr<ParamCurve> crv1 = kinks[ki].first->geomCurve();
      shared_ptr<CurveOnSurface> sfcv = 
	dynamic_pointer_cast<CurveOnSurface,ParamCurve>(crv1);
      if (sfcv.get())
	{
	  sfcv->spaceCurve()->writeStandardHeader(out_file2);
	  sfcv->spaceCurve()->write(out_file2);
	}
      else
	{
	  crv1->writeStandardHeader(out_file2);
	  crv1->write(out_file2);
	}
  }

  vector<pair<ftEdge*, ftEdge*> > gaps2;
  vector<pair<ftEdge*, ftEdge*> > kinks2;
  quality.facePositionDiscontinuity(gaps2);  
  quality.faceTangentDiscontinuity(kinks2);

  std::ofstream out_filen("gaps2.g2");
  for (ki=0; ki<gaps2.size(); ++ki)
  {
      shared_ptr<ParamCurve> crv1 = gaps2[ki].first->geomCurve();
      shared_ptr<CurveOnSurface> sfcv = 
	dynamic_pointer_cast<CurveOnSurface,ParamCurve>(crv1);
      if (sfcv.get())
	{
	  sfcv->spaceCurve()->writeStandardHeader(out_filen);
	  sfcv->spaceCurve()->write(out_filen);
	}
      else
	{
	  crv1->writeStandardHeader(out_filen);
	  crv1->write(out_filen);
	}
  }

  std::ofstream out_filen2("kinks2.g2");
  for (ki=0; ki<kinks2.size(); ++ki)
  {
      shared_ptr<ParamCurve> crv1 = kinks2[ki].first->geomCurve();
      shared_ptr<CurveOnSurface> sfcv = 
	dynamic_pointer_cast<CurveOnSurface,ParamCurve>(crv1);
      if (sfcv.get())
	{
	  sfcv->spaceCurve()->writeStandardHeader(out_file);
	  sfcv->spaceCurve()->write(out_filen2);
	}
      else
	{
	  crv1->writeStandardHeader(out_filen2);
	  crv1->write(out_file);
	}
  }
}

