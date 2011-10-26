//===========================================================================
//                                                                           
// File: PrExplicitConnectivity.C                                            
//                                                                           
// Created: Wed Mar 29 16:52:57 2006                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: PrExplicitConnectivity.C,v 1.3 2007-03-02 16:35:52 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/parametrization/PrExplicitConnectivity.h"
#include<vector>

using namespace std;

//----------------------------------------------------------------------------
int PrExplicitConnectivity::findNumFaces() const
//-----------------------------------------------------------------------------
{
  int npts = getNumNodes();
  vector<int> neighbours, face;

  int i,j,deg,it = 0;
  for(i=0; i< npts; i++)
  {
    getNeighbours(i,neighbours);
    deg = (int)neighbours.size();
    for(j=0; j<deg-1; j++)
    {
      findFace(i,j,neighbours,face);
      if(isMinimum(i,face)) it++;
    }
    if(!isBoundary(i))
    {
      findFace(i,j,neighbours,face);
      if(isMinimum(i,face)) it++;
    }
  }
  return it;
}

//----------------------------------------------------------------------------
void PrExplicitConnectivity::findFace(int i, int j,
				      vector<int>& neighbours, vector<int>& face) const
//-----------------------------------------------------------------------------
{
/* Find the face of the graph which
// lies to the left of the directed edge whose end points
// are the node i and its j-th neighbour.
// The nodes of the face, starting with i, will be filled out in the
// list vector<int> face. The vector<int> neighbours is used for temporarily
// storing neighbours (for efficiency when calling this routine many times).
//
//   0-----------0
//   |            \ 
//   |             \ 
//   |             0 j-th neighbour of i
//   |            /
//   |           /
//   |          /
//   O---------0 i
// (j+1)-th
// neighbour of i
//
*/

  face.clear();
  face.push_back(i);
  getNeighbours(i,neighbours);
  int node = i;
  int ind  = j;

  while(neighbours[ind] != i)
  {
    findNextEdgeInFace(node,ind,neighbours);
    face.push_back(node);
  }

}

//----------------------------------------------------------------------------
int PrExplicitConnectivity::findGenus() const
//-----------------------------------------------------------------------------
{
  return getNumNodes() - findNumEdges() + findNumFaces();
}

//----------------------------------------------------------------------------
bool PrExplicitConnectivity::isTriangulation() const
//-----------------------------------------------------------------------------
{
  int npts = getNumNodes();
  vector<int> neighbours, face;

  int i,j,deg = 0;
  for(i=0; i< npts; i++)
  {
    getNeighbours(i,neighbours);
    deg = (int)neighbours.size();
    for(j=0; j<(deg-1); j++)
    {
      findFace(i,j,neighbours,face);
      if(face.size() != 3) return false;
    }
    if(!isBoundary(i))
    {
      findFace(i,j,neighbours,face);
      if(face.size() != 3) return false;
    }
  }
  return true;
}

//----------------------------------------------------------------------------
void PrExplicitConnectivity::findNextEdgeInFace(int& i, int& j,
						vector<int>& neighbours) const
//-----------------------------------------------------------------------------
{
  int nj = neighbours[j];
  getNeighbours(nj,neighbours);
  if(neighbours[0] == i)
  {
    j = (int)neighbours.size() - 1;
    i = nj;
    return;
  }
  for(size_t k=1; k<neighbours.size(); k++)
    if(neighbours[k] == i)
    {
	j = (int)k-1;
      i = nj;
      break;
    }
}


//----------------------------------------------------------------------------
void PrExplicitConnectivity::printXYZFaces(std::ostream& os) const
//-----------------------------------------------------------------------------
{
  int npts = getNumNodes();
  vector<int> neighbours, face;

  int i,j,k,deg;
  for(i=0; i< npts; i++)
  {
    getNeighbours(i,neighbours);
    deg = (int)neighbours.size();
    for(j=0; j<(deg-1); j++)
    {
      findFace(i,j,neighbours,face);
      if(isMinimum(i,face))
      {
        for(k=0; k<int(face.size()); k++)
          get3dNode(face[k]).write(os);
//	    os << get3dNode(face[k]);
        os << std::endl;
      }
    }
    if(!isBoundary(i))
    {
      findFace(i,j,neighbours,face);
      if(isMinimum(i,face))
      {
        for(k=0; k<int(face.size()); k++)
          get3dNode(face[k]).write(os);
//	    os << get3dNode(face[k]);
        os << "\n";
      }
    }
  }
}

//----------------------------------------------------------------------------
void PrExplicitConnectivity::printXYZFacesML(std::ostream& os) const
//-----------------------------------------------------------------------------
{
  int npts = getNumNodes();
  vector<int> neighbours, face;
  // VERTICES
  os << "v = [\n";

  int i,j,k,deg;
  Vector3D node;
  for(i=0; i< npts; i++)
  {
    node = get3dNode(i);
//    os << node;
    node.write(os);
  }
  os << "];\n";

  // FACES
  os << "f = [\n";
  for(i=0; i< npts; i++)
  {
    getNeighbours(i,neighbours);
    deg = (int)neighbours.size();
    for(j=0; j<(deg-1); j++)
    {
      findFace(i,j,neighbours,face);
      if(isMinimum(i,face))
      {
        for(k=0; k<int(face.size()); k++) os << face[k] << " ";
        os << ";\n";
      }
    }
    if(!isBoundary(i))
    {
      findFace(i,j,neighbours,face);
      if(isMinimum(i,face))
      {
        for(k=0; k<int(face.size()); k++) os << face[k] << " ";
        os << ";\n";
      }
    }
  }
  os << "];\n";

  // PRINTING
  os << "view(3);\n";
  os << "light('Position',[0 0 1]);\n";
  os << "light('Position',[-1 -1 0]);\n";
  os << "patch('Vertices', v, 'Faces', f, 'FaceColor', 'y');\n";
}

//----------------------------------------------------------------------------
void PrExplicitConnectivity::printXYZFacesVRML(std::ostream& os) const
//-----------------------------------------------------------------------------
{
  int npts = getNumNodes();
  vector<int> neighbours, face;

  os << "#VRML V1.0 ascii \n";
  os << "\n";
  os << "Separator {\n";
  os << "    ShapeHints {\n";
  os << "	vertexOrdering	COUNTERCLOCKWISE\n";
  os << "	shapeType	SOLID\n";
  os << "	faceType	CONVEX\n";
  os << "    }\n";
  os << "    Separator {\n";
  os << "	Coordinate3 {\n";
  os << "	    point	[\n";
  
  
  int i,j,k,deg;
  Vector3D node;
  for(i=0; i< npts; i++)
  {
    node = get3dNode(i);
    node.write(os);
    if (i<npts-1)
      os << ",\n";
  }
  os << "]\n";
  os << "}\n";
  os << " IndexedFaceSet {\n";
  os << "    coordIndex [\n\t";
  
  // FACES
  for(i=0; i< npts; i++)
  {
    getNeighbours(i,neighbours); 
    deg = (int)neighbours.size();
    for(j=0; j<(deg-1); j++)
    {
      findFace(i,j,neighbours,face);
      if(isMinimum(i,face))
      {
	os << "       ";
	for(k=0; k<int(face.size()); k++) 
	  os << face[k] << ", ";
        os << "-1,\n";
      }
    }
    if(!isBoundary(i))
    {
      findFace(i,j,neighbours,face);
      if(isMinimum(i,face))
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

//----------------------------------------------------------------------------
void PrExplicitConnectivity::printUVFaces(std::ostream& os) const
//-----------------------------------------------------------------------------
{
  int npts = getNumNodes();
  vector<int> neighbours, face;

  int i,j,k,deg;
  for(i=0; i< npts; i++)
  {
    getNeighbours(i,neighbours);
    deg = (int)neighbours.size();
    for(j=0; j<(deg-1); j++)
    {
      findFace(i,j,neighbours,face);
      if(isMinimum(i,face))
      {
        for(k=0; k<int(face.size()); k++)
	    os << getU(face[k]) << ' ' << getV(face[k]) << std::endl;
        os << "\n";
      }
    }
    if(!isBoundary(i))
    {
      findFace(i,j,neighbours,face);
      if(isMinimum(i,face))
      {
        for(k=0; k<int(face.size()); k++)
	    os << getU(face[k]) << ' ' << getV(face[k]) << std::endl;
        os << "\n";
      }
    }
  }
}

//----------------------------------------------------------------------------
void PrExplicitConnectivity::printTexture(std::ostream& os) const
//-----------------------------------------------------------------------------
{

  int npts = getNumNodes();
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
     os << "        " << getU(i) << " " << getV(i);
     // transforme un point de [-1,1]x[-1,1] en un point de [0,1]x[0,1]
     // lorsque la parametrisation est sur le cercle unite
     // os << "        " << (getU(i)+1)/2;
     // os << " " << (getV(i)+1)/2;

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
    node = get3dNode(i);
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
    getNeighbours(i,neighbours); 
    deg = (int)neighbours.size();
    for(j=0; j<(deg-1); j++)
    {
      findFace(i,j,neighbours,face);
      if(isMinimum(i,face))
      {
	os << "       ";
	for(k=0; k<int(face.size()); k++) 
	  os << face[k] << ", ";
        os << "-1,\n";
      }
    }
    if(!isBoundary(i))
    {
      findFace(i,j,neighbours,face);
      if(isMinimum(i,face))
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

//----------------------------------------------------------------------------
void PrExplicitConnectivity::printInfo(std::ostream& os) const
//-----------------------------------------------------------------------------
{
  os << "Number of nodes = " << getNumNodes() << std::endl;
  os << "Number of edges = " << findNumEdges() << std::endl;
  os << "Number of faces = " << findNumFaces() << std::endl;
  os << "Genus = " << findGenus() << std::endl;
  os << "Number of boundary nodes = " << findNumBdyNodes() << std::endl;
  if(isTriangulation()) os << "It is a triangulation" << std::endl;
  else os << "It is not a triangulation" << std::endl;
  os << "Number of connected components = ";
  os << findNumComponents() << std::endl;
  os << "Number of connected boundary components = ";
  os << findNumBdyComponents() << std::endl;
}





