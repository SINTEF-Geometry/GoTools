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
#include "GoTools/igeslib/IGESconverter.h"
#include <fstream>

using std::istream;
using std::vector;
using namespace Go;

ftMessage readg2(istream& is, vector<shared_ptr<ftSurface> >& faces)
//===========================================================================
{

  ftMessage status;
  if (is.bad())
    {
      status.setError(FT_BAD_INPUT_FILE);
      return status;
    }

  IGESconverter conv;
  try
    {
      conv.readgo(is);
    }
  catch (...)
    {
      status.setError(FT_ERROR_IN_READ_IGES);
      return status;
    }

  //      std::ofstream outfile("debug.out");
  //      conv.writedisp(outfile);

  vector<shared_ptr<GeomObject> > gogeom = conv.getGoGeom();

  int nmbgeom = (int)gogeom.size();
  faces.reserve(nmbgeom); // May be too much, but not really important
  int face_count = 0;


  // Remaining geometry.

  for (int i=0; i<nmbgeom; i++)
    {

      if (gogeom[i]->instanceType() == Class_SplineCurve)
	{
	  if (conv.getGroup().size() == 0)
	    {
	      // One mesh surface expected.

	    }
	  else
	    // Not expected. Ignore the current curve.
	    status.addWarning(FT_UNEXPECTED_INPUT_OBJECT_IGNORED);  
	}
      else if (gogeom[i]->instanceType() == Class_SplineSurface ||
	       gogeom[i]->instanceType() == Class_BoundedSurface)
	{
	  shared_ptr<GeomObject> lg = gogeom[i];
	  shared_ptr<ParamSurface> gosf =
	    dynamic_pointer_cast<ParamSurface, GeomObject>(lg);
	  shared_ptr<ftSurface> ftsf(new ftSurface(gosf, face_count++));
	  faces.push_back(ftsf);

	}
    }

  return status;

}



int main( int argc, char* argv[] )
{
  if (argc != 2) {
    std::cout << "Input parameters : Input file on g2 format," << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

    double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double approxtol = 0.01;

  vector<shared_ptr<ftSurface> > faces;
  ftMessage status = readg2(file1, faces);
  std::cout << "Read g2 file. Status message : " << status.getMessage() << std::endl;

  SurfaceModel model1(approxtol, gap, neighbour, kink, 10.0*kink, faces);
  int nmb_bd1 = model1.nmbBoundaries();
  std::cout << "No of boundaries (1): " << nmb_bd1 << std::endl;

  SurfaceModel model2(approxtol, gap, neighbour, kink, 10.0*kink, faces, true);
  int nmb_bd2 = model2.nmbBoundaries();
  std::cout << "No of boundaries (2): " << nmb_bd2 << std::endl;

  exit(0);
}



		      

