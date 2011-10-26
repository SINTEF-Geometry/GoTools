#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include <fstream>

using std::cout;
using std::cin;
using std::vector;
using std::shared_ptr;
using std::dynamic_pointer_cast;
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
  int insert = 0;
  if (argc == 3)
    insert = atoi(argv[2]);

    double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  CompositeModel *model = factory.createFromG2(file1);

  SurfaceModel *sfmodel = dynamic_cast<SurfaceModel*>(model);

  if (sfmodel && insert)
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

  if (sfmodel)
    {
      bool splines = sfmodel->allSplines();
      std::cout << "All surfaces are splines: " << splines << std::endl;

      if (splines)
	{
	  bool corner = sfmodel->isCornerToCorner();
	  std::cout << "Corner-to-corner configuaration: " << corner << std::endl;
	  if (!corner)
	    {
	      sfmodel->makeCornerToCorner();
	      std::ofstream out_file("corner2corner.g2");
	      int nmb = sfmodel->nmbEntities();
	      for (int ki=0; ki<nmb; ki++)
		{
		  shared_ptr<ParamSurface> surf = sfmodel->getSurface(ki);

		  surf->writeStandardHeader(out_file);
		  surf->write(out_file);
		}
	    }

	  sfmodel->makeCommonSplineSpaces();
	  std::ofstream out_file("splinespace.g2");
	  int nmb = sfmodel->nmbEntities();
	  int ki;
	  for (ki=0; ki<nmb; ki++)
	    {
	      shared_ptr<ParamSurface> surf = sfmodel->getSurface(ki);

	      surf->writeStandardHeader(out_file);
	      surf->write(out_file);
	    }
	    
	  std::ofstream out_file2("adjacency.dat");
	  for (ki=0; ki<nmb; ++ki)
	    {
	      shared_ptr<ftSurface> face = sfmodel->getFace(ki);

	      vector<int> free_bd;
	      bool found_free = face->getFreeBoundaryInfo(gap, free_bd);
	      out_file2 << "Surface " << ki << ", result: " << found_free;
	      out_file2 << ", number: " << free_bd.size() << std::endl;
	      size_t kj;
	      for (kj=0; kj<free_bd.size(); ++kj)
		out_file2 << free_bd[kj] << "  ";
	      out_file2 << std::endl << std::endl;

	      for (kj=0; kj<free_bd.size(); ++kj)
		{
		  vector<int> bd_enumeration;
		  found_free = face->getBoundaryCoefEnumeration(free_bd[kj],
								bd_enumeration);
		  out_file2 << "Boundary: " << kj << std::endl;
		  for (size_t kr=0; kr<bd_enumeration.size(); ++kr)
		    out_file2 << bd_enumeration[kr] << " ";
		  out_file2 << std::endl;
		}
	      out_file2 << std::endl;

  	      vector<ftSurface*> neighbours;
	      face->getAdjacentFaces(neighbours);
	      for (size_t kj=0; kj<neighbours.size(); ++kj)
		{
		  vector<pair<int,int> > enumeration;
		  bool found;
		  found = face->getCorrCoefEnumeration(neighbours[kj],
							    gap, enumeration);
		  int idx = sfmodel->getIndex(neighbours[kj]);
		  out_file2 << "Faces: " << ki << " " << idx << std::endl;
		  for (size_t kr=0; kr<enumeration.size(); ++kr)
		    {
		      out_file2 << enumeration[kr].first << " ";
		      out_file2 << enumeration[kr].second << std::endl;
		    }
		}

	      for (int kj=ki+1; kj<nmb; ++kj)
		{
		  shared_ptr<ftSurface> face2 = sfmodel->getFace(kj);
		  int bd1, bd2;
		  bool same;
		  bool adjacent = face->getAdjacencyInfo(face2.get(), gap,
							 bd1, bd2, same);
		  if (adjacent)
		    {
		      vector<pair<int,int> > enumeration;
		      bool found;
		      found = face->getCorrCoefEnumeration(face2.get(),
								gap, enumeration);
		      out_file2 << "Adjacent faces: " << ki << " " << kj << std::endl;
		      out_file2 << "Edges: " << bd1 << " " << bd2 << std::endl;
		      for (size_t kr=0; kr<enumeration.size(); ++kr)
			{
			  out_file2 << enumeration[kr].first << " ";
			  out_file2 << enumeration[kr].second << std::endl;
			}
		    }
		}
		      
	    }
	}
    }
}

	  
