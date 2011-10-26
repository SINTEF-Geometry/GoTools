/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1997 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/


/********************************************************************
 FILENAME    : PrTexture.C
 AUTHOR      : Valerie PHAM-TRONG, SINTEF
 DATE        : March 2002
 DESCRIPTION : Implementation of methods of texture mapping
 in the class PrOrganizePoints.
 CHANGE LOG  :
*********************************************************************/

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

