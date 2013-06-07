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
#include "GoTools/qualitymodule/FaceSetRepair.h"
#include <fstream>

using std::cout;
using std::cin;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 2 && argc != 3) {
    std::cout << "Input parameters : Input file on g2 format, (Insert knots)" << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

    double gap = 0.0001;
    double neighbour = 0.001;
  double kink = 0.01;
  double approxtol = 0.01;
  int insert = 0;
  if (argc == 3)
    insert = atoi(argv[2]);

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  shared_ptr<CompositeModel> model = shared_ptr<CompositeModel>(factory.createFromG2(file1));
  shared_ptr<SurfaceModel> sfmodel = dynamic_pointer_cast<SurfaceModel,CompositeModel>(model);

  if (insert)
    {
      cout << "Number of surfaces: " << sfmodel->nmbEntities() << std::endl;
      cout << "Surface to refine: ";
      int idx;
      cin >> idx;

      cout << "Parameter direction: ";
      int dir;
      cin >> dir;

      cout << "Number of knots: ";
      int nmb;
      cin >> nmb;

      cout << "Knots: ";
      vector<double> knots(nmb);
      for (int ki=0; ki<nmb; ++ki)
	cin >> knots[ki];

      shared_ptr<ParamSurface> srf = sfmodel->getSurface(idx);
      shared_ptr<SplineSurface> s1 = 
	dynamic_pointer_cast<SplineSurface,ParamSurface>(srf);
      if (s1.get())
	{
	  if (dir == 0)
	    s1->insertKnot_u(knots);
	  else
	    s1->insertKnot_v(knots);
	}
    }
      
  shared_ptr<FaceSetQuality> quality = 
    shared_ptr<FaceSetQuality>(new FaceSetQuality(gap, kink, approxtol));
  quality->attach(sfmodel);

  shared_ptr<FaceSetRepair> repair = 
    shared_ptr<FaceSetRepair>(new FaceSetRepair(quality));
 
  vector<pair<ftEdge*, ftEdge*> > gaps;
  quality->facePositionDiscontinuity(gaps);  

  std::cout << "Number of gaps: " << gaps.size() << std::endl;

  std::ofstream out_file("gaps_1.g2");
  size_t ki;
  for (ki=0; ki<gaps.size(); ++ki)
  {
    shared_ptr<ParamCurve> cv1 = gaps[ki].first->geomCurve();
    shared_ptr<ParamCurve> cv2 = gaps[ki].second->geomCurve();
    shared_ptr<SplineCurve> scv1 = shared_ptr<SplineCurve>(cv1->geometryCurve());
    shared_ptr<SplineCurve> scv2 = shared_ptr<SplineCurve>(cv2->geometryCurve());
    scv1->writeStandardHeader(out_file);
    scv1->write(out_file);
    scv2->writeStandardHeader(out_file);
    scv2->write(out_file);
  }

  repair->mendEdgeDistance();  

  std::ofstream out_model0("sfmodel0.g2");
  int nmb = sfmodel->nmbEntities();
  for (int ki=0; ki<nmb; ki++)
    {
      shared_ptr<ParamSurface> surf = sfmodel->getSurface(ki);

      surf->writeStandardHeader(out_model0);
      surf->write(out_model0);
    }

  repair->mendGaps();

  std::ofstream out_model("sfmodel.g2");
  nmb = sfmodel->nmbEntities();
  for (int ki=0; ki<nmb; ki++)
    {
      shared_ptr<ParamSurface> surf = sfmodel->getSurface(ki);

      surf->writeStandardHeader(out_model);
      surf->write(out_model);
    }

			      
  vector<pair<ftEdge*, ftEdge*> > gaps2;
  quality->facePositionDiscontinuity(gaps2);  

  std::cout << "Number of gaps after mending: " << gaps2.size() << std::endl;

  std::ofstream out_filen("gaps_2.g2");
  for (ki=0; ki<gaps2.size(); ++ki)
  {
    shared_ptr<ParamCurve> cv1 = gaps2[ki].first->geomCurve();
    shared_ptr<ParamCurve> cv2 = gaps2[ki].second->geomCurve();
    shared_ptr<SplineCurve> scv1 = shared_ptr<SplineCurve>(cv1->geometryCurve());
    shared_ptr<SplineCurve> scv2 = shared_ptr<SplineCurve>(cv2->geometryCurve());
    scv1->writeStandardHeader(out_filen);
    scv1->write(out_filen);
    scv2->writeStandardHeader(out_filen);
    scv2->write(out_filen);
  }

}

