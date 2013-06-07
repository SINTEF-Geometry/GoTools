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

#include "GoTools/parametrization/PrExplicitConnectivity.h"
//#include "GoTools/parametrization/PrDijkstra.h"
#include "GoTools/parametrization/PrTexture.h"

//----------------------------------------------------------------------------
void printTexture(std::ostream& os, PrExplicitConnectivity& triang )
//-----------------------------------------------------------------------------
{
  int npts = triang.getNumNodes();
  vector<int> neighbours, face;


  os << "#Inventor V2.1 ascii \n";
  os << "\n";

  // TEXTURE COORDINATES VERTICES IN 2D 
  
  std::cout << "\nPlease enter the name of the file containing the texture: ";
  char texture_filename[300];
  std::cin >> texture_filename;
  std::cout << "Writing file for the texture mapping ..." << std::endl;
  
  os << "Separator {\n";
  os << "Texture2 {\n";
  os << "    filename	\"./" << texture_filename << "\" " << std::endl;
  os << "    model	MODULATE\n";
  os << "}\n";
  os << "TextureCoordinate2 {\n";
  os << "    point	[ \n";

  int i;  
  for (i=0; i<npts; i++)
  { 
     os << "        " << triang.getU(i) << " " << triang.getV(i);
     // transforme un point de [-1,1]x[-1,1] en un point de [0,1]x[0,1]
     // lorsque la parametrisation est sur le cercle unite
     // os << "        " << (triang.getU(i)+1)/2;
     // os << " " << (triang.getV(i)+1)/2;

     if (i<npts-1) 
      {
	  os << "," ;
      } 
      os << std::endl;
  }
  os << "    ]\n";
  os << "}\n";
  

  os << "    TextureCoordinateBinding {\n";
  os << "      value	PER_VERTEX_INDEXED\n";
  os << "    }";

  // CORRESPONDING VERTICES IN 3D

  os << "    ShapeHints {\n";
  os << "	vertexOrdering	COUNTERCLOCKWISE\n";
  os << "	shapeType	SOLID\n";
  os << "	faceType	CONVEX\n";
  os << "    }\n";
 
  int j,k,deg;
  Vector3D node;

  os << "    Separator {\n";
  os << "	Coordinate3 {\n";
  os << "	    point	[\n";
  
  for(i=0; i< npts; i++)
  {
    node = triang.get3dNode(i);
    node.write(os);
    if (i<npts-1)
      os << ",\n";
  }
  os << "]\n";
  os << "}\n";

  // FACES

  os << " IndexedFaceSet {\n";
  os << "    coordIndex [\n\t";
  for(i=0; i< npts; i++)
  {
    triang.getNeighbours(i,neighbours); 
    deg = (int)neighbours.size();
    for(j=0; j<(deg-1); j++)
    {
      triang.findFace(i,j,neighbours,face);
      if(triang.isMinimum(i,face))
      {
	os << "       ";
	for(k=0; k<int(face.size()); k++) 
	  os << face[k] << ", ";
        os << "-1,\n";
      }
    }
    if(!triang.isBoundary(i))
    {
      triang.findFace(i,j,neighbours,face);
      if(triang.isMinimum(i,face))
      {
	os << "       ";
        for(k=0; k<int(face.size()); k++) 
	  os << face[k] << ", ";
	os << "-1,\n";
      }
    }
  }

  os << "       ]\n";
  os << "    }\n";
  os << "  }\n";
  os << "}\n";

}

