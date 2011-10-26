/********************************************************************
 FILENAME    : PrGeomview.C
 AUTHOR      : Valerie PHAM-TRONG, SINTEF
 DATE        : Mai 2002
 DESCRIPTION : Computation of a geodesic path between two vertices
               in the class PrOrganizePoints.
 CHANGE LOG  :
*********************************************************************/

#include <iostream>

#include "GoTools/parametrization/PrGeodesics.h"
#include "GoTools/parametrization/PrPathTriangleSeq.h"
#include "GoTools/parametrization/PrInterface.h"

#include "GoTools/parametrization/PrOrganizedPoints.h"
#include "GoTools/parametrization/PrTriangulation_OP.h"
#include "GoTools/parametrization/PrTriangle.h"
#include "GoTools/parametrization/PrDijkstra.h"
#include "GoTools/utils/Array.h"

using namespace Go;
using namespace std;

//----------------------------------------------------------------------------
void print(vector<int>& v)
//----------------------------------------------------------------------------
{
    std::vector<int>::iterator iter;
    for (iter = v.begin(); iter != v.end(); iter++)
    {
	cout << *iter << " - ";
    }
    cout << endl;
}

//----------------------------------------------------------------------------
void print(vector<double>& v)
//----------------------------------------------------------------------------
{
    std::vector<double>::iterator iter;
    for (iter = v.begin(); iter != v.end(); iter++)
    {
	cout << *iter << " - ";
    }
    cout << endl;
}

//----------------------------------------------------------------------------
void print(vector<Vector3D>& v)
//----------------------------------------------------------------------------
{
    std::vector<Vector3D>::iterator iter;
    for (iter = v.begin(); iter !=  v.end(); iter++)
    {
	std::cout << "  " << (*iter) ;
    }
    cout << endl;
}

//----------------------------------------------------------------------------
void print(vector<Vector2D>& v)
//----------------------------------------------------------------------------
{
    std::vector<Vector2D>::iterator iter;
    for (iter = v.begin(); iter != v.end(); iter++)
    {
	std::cout << "  " << (*iter) ;
    }
    cout << endl;
}


//----------------------------------------------------------------------------
void print(vector<PrTriangle>& t)
//----------------------------------------------------------------------------
{
    std::vector<PrTriangle>::iterator iter;
    for (iter = t.begin(); iter != t.end(); iter++)
    {
	cout << "\n" << iter->n1();
	cout << " " << iter->n2();
	cout << " " << iter->n3();
    }
    cout << endl;
}

//----------------------------------------------------------------------------
void printGeomviewTriangulation(PrTriangulation_OP& t)
//-----------------------------------------------------------------------------
/*
  Writing file for the triangulation visualisation with geomview
 */
{  
    std::ofstream os("geom_triangulation.off");

    std::cout <<
      "Writing file for the triangulation visualisation with geomview..."
	 << std::endl;
    int nb_nodes = t.getNumNodes();
    int nb_tri = t.findNumFaces();
    int i;

    os.precision(5) ;
    os << "#Geomview Triangulation \n";
    os << "OFF\n";
    os << nb_nodes << " " << nb_tri+1;
    os << "\n";

    os << "\n#Vertices\n\n";
    for (i=0; i<nb_nodes; i++)	
    {
	os << t.get3dNode(i); 
    }	
    os << "\n";
    
    os << "\n#Faces\n\n";	
    os << "3 0 0 0 \n";
    for (i=0; i<nb_tri; i++)
    {	
	os << "3 ";
	os << t.getPrTriangle(i).n1() << " " ;
	os << t.getPrTriangle(i).n2() << " " ;
	os << t.getPrTriangle(i).n3() << "\n";
    }
}

//----------------------------------------------------------------------------
void printGeomviewTriangulation(std::ofstream& os, PrTriangulation_OP& t)
//-----------------------------------------------------------------------------
/*
  Writing file for the triangulation visualisation with geomview
  problem in giving the name function as parameter.
 */
{  
  std::cout << "Writing file for the triangulation visualisation with geomview ...";
  std::cout << std::endl;
  int nb_nodes = t.getNumNodes();
  int nb_tri = t.findNumFaces();
  int i;

  os.precision(5) ;
  os << "#Geomview Triangulation \n";
  os << "OFF\n";
  os << nb_nodes << " " << nb_tri+1;
  os << "\n";

  os << "\n#Vertices\n\n";
  for (i=0; i<nb_nodes; i++)	
  {
      os << t.get3dNode(i); 
  }	
  os << "\n";

  os << "\n#Faces\n\n";	
  os << "3 0 0 0 \n";
  for (i=0; i<nb_tri; i++)
  {	
      os << "3 ";
      os << t.getPrTriangle(i).n1() << " " ;
      os << t.getPrTriangle(i).n2() << " " ;
      os << t.getPrTriangle(i).n3() << "\n";
  }
}
//----------------------------------------------------------------------------
void printGeomviewTriangleSequence(PrTriangulation_OP& t, vector<int> tr_seq)
//----------------------------------------------------------------------------
/*
  Writing file for the triangle sequence visualisation with geomview
 */
{
    std::ofstream os("geom_sequence.off");

    int nb_tri = (int) tr_seq.size();
    int i;
 
    std::cout << "Writing file for the triangle sequence visualisation with geomview ...";
    std::cout << std::endl;

    os.precision(5) ;
    os << "#Geomview Triangle Sequence \n";
    os << "\n";
    os << "OFF" << "\n";
    os << 3*nb_tri << " " << nb_tri;
    os << "\n\n";

    os << "#Vertices\n\n" ;
    for (i=0; i<nb_tri; i++)	
    {
	(t.get3dNode(t.getPrTriangle(tr_seq[i]).n1())).write(os); 
	(t.get3dNode(t.getPrTriangle(tr_seq[i]).n2())).write(os); 
	(t.get3dNode(t.getPrTriangle(tr_seq[i]).n3())).write(os); 
	os << "\n"; 
    }	
    os << "\n";
    
    os <<  "#Faces\n\n" ;	
    os << "3 0 0 0 \n";
    for (i=0; i<nb_tri; i++)
    {	
	os << "3 ";
	os << 3*i << " " ;
	os << 3*i+1 << " " ;
	os << 3*i+2 << "\n";
    }
}

//----------------------------------------------------------------------------
void printGeomviewTriangleSequence(std::ofstream& os, PrTriangulation_OP& t, vector<int> tr_seq)
//-----------------------------------------------------------------------------
/*
  Writing file for the triangle sequence visualisation with geomview
  problem in giving the name function as parameter.
 */
{
    int nb_tri = (int) tr_seq.size();
    int i;
    
    std::cout << "Writing file for the triangle sequence visualisation with geomview ...";
    std::cout << std::endl;

    os.precision(5) ;
    os << "#Geomview Triangle Sequence \n";
    os << "\n";
    os << "OFF" << "\n";
    os << 3*nb_tri << " " << nb_tri;
    os << "\n\n";
    
    os << "#Vertices\n\n" ;
    for (i=0; i<nb_tri; i++)	
    {
	(t.get3dNode(t.getPrTriangle(tr_seq[i]).n1())).write(os); 
	(t.get3dNode(t.getPrTriangle(tr_seq[i]).n2())).write(os); 
	(t.get3dNode(t.getPrTriangle(tr_seq[i]).n3())).write(os); 
	os << "\n"; 
    }	
    os << "\n";

    os <<  "#Faces\n\n" ;	
    os << "3 0 0 0 \n";
    for (i=0; i<nb_tri; i++)
    {	
	os << "3 ";
	os << 3*i << " " ;
	os << 3*i+1 << " " ;
	os << 3*i+2 << "\n";
    }

    return;
}

//----------------------------------------------------------------------------
void printGeomviewPath(vector<Vector3D>& path)
//-----------------------------------------------------------------------------
/*
  Writing file for the shortest path visualisation with geomview
 */
{  
    std::ofstream os("geom_path.vect");

    int nb = (int) path.size();
    std::cout << "Writing file for the shortest path visualisation with geomview ...";
    std::cout << std::endl;

    os.precision(5);
    os << "#Geomview Shortest Path " << endl;
    os << "\n";
    os << "VECT" << "\n";
    os << "1 " << nb << " 1" << endl;
    os <<  nb << endl;
    os <<  "1" << endl;

    
    int i;
    for (i=0; i<nb; i++)	
    {
	path[i].write(os);
    }
    os << endl; 
    os << "0 0 1 1 " << endl;
}


//----------------------------------------------------------------------------
void printGeomviewPath(std::ofstream& os, vector<Vector3D>& path)
//-----------------------------------------------------------------------------
/*
  Writing file for the shortest path visualisation with geomview
  problem in giving the name function as parameter.
 */
{  
    int nb = (int) path.size();
    std::cout << "Writing file for the shortest path visualisation with geomview ...";
    std::cout << std::endl;

    os.precision(5);
    os << "#Geomview Shortest Path " << endl;
    os << "\n";
    os << "VECT" << "\n";
    os << "1 " << nb << " 1" << endl;
    os <<  nb << endl;
    os <<  "1" << endl;

    
    int i;
    for (i=0; i<nb; i++)	
    {
	path[i].write(os);
    }
    os << endl; 
    os << "0 0 1 1 " << endl;
}

//----------------------------------------------------------------------------
void print_edges_lengths(vector<int> tr_seq, PrTriangulation_OP& t)
//----------------------------------------------------------------------------
{
    std::vector<int>::iterator iter;
    Vector3D a, b, c;
    cout << "\nedges lengths in the 3d sequence : " << endl;
    for (iter = tr_seq.begin(); iter != tr_seq.end(); iter++)
    {
	a = t.get3dNode(t.getPrTriangle(*iter).n1());
	b = t.get3dNode(t.getPrTriangle(*iter).n2());
	c = t.get3dNode(t.getPrTriangle(*iter).n3());
	cout << "triangle " << *iter << " : ";
	cout <<  dist3D(a, b); 
	cout << " - ";
	cout <<  dist3D(b, c); 
	cout << " - ";
	cout <<  dist3D(c, a); 
	cout << endl;
    }
    cout << endl;
    
}

//----------------------------------------------------------------------------
void print_lengths(const vector<double>& lengths)
//----------------------------------------------------------------------------
{
    std::vector<double>::const_iterator iter;
    int cpt = 1;
    cout << "\nLength evolution during the " << lengths.size() 
	 << " iterations :\n"; 
    for (iter = lengths.begin(); iter != lengths.end(); iter ++)
    {
	cout << *iter << " > ";
	if (cpt ++ == 6)
	{
	    cout << endl; 
	    cpt = 1;
	}
    }
    cout << endl; 
}
