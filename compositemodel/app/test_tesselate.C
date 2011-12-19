#ifdef __BORLANDC__
#include <vcl.h>
#endif

#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/geometry/LineCloud.h"
#include <fstream>
#include <stdlib.h> // For atof()

using namespace std;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 6) {
    std::cout << "Input parameters : Input file, IGES or g2 (1/0), n, m, density"  << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  // Removed from corresponding use of fairingToolbox
  double approx = 0.001;
  int useIGES = atoi(argv[2]);
  int n = atoi(argv[3]);
  int m = atoi(argv[4]);
  double density = atof(argv[5]);

  CompositeModelFactory factory(approx, gap, neighbour, kink, 10.0*kink);

  CompositeModel *model;
  if (useIGES)
      model = factory.createFromIges(file1);
  else
      model = factory.createFromG2(file1);

  SurfaceModel *sfmodel = dynamic_cast<SurfaceModel*>(model);
  CompositeCurve *cvmodel = dynamic_cast<CompositeCurve*>(model);


  std::vector<shared_ptr<GeneralMesh> > meshes;
  int res[2];
  res[0] = n;
  res[1] = m;
  model->tesselate(res, meshes);

  std::vector<shared_ptr<GeneralMesh> > meshes2;
  model->tesselate(density, meshes2);

  if (sfmodel)
    {
  std::vector<shared_ptr<GeneralMesh> > meshes3;
  sfmodel->tesselate(n*m, meshes3);

  size_t ki;
  std::ofstream out1("triangles1.g2");
  for (ki=0; ki<meshes.size(); ++ki)
    {
      double *nodes = meshes[ki]->vertexArray();
      int num_triang = meshes[ki]->numTriangles();
      int num_triang_nodes = 3*num_triang;
      unsigned int *triang_idx = meshes[ki]->triangleIndexArray();
      vector<double> seg;
      for (int kj=0; kj<num_triang_nodes; kj+=3)
	{
	  Point pt1(&nodes[3*triang_idx[kj]],&nodes[3*(triang_idx[kj]+1)]);
	  Point pt2(&nodes[3*triang_idx[kj+1]],&nodes[3*(triang_idx[kj+1]+1)]);
	  Point pt3(&nodes[3*triang_idx[kj+2]],&nodes[3*(triang_idx[kj+2]+1)]);
	  seg.insert(seg.end(), pt1.begin(), pt1.end());
	  seg.insert(seg.end(), pt2.begin(), pt2.end());
	  seg.insert(seg.end(), pt2.begin(), pt2.end());
	  seg.insert(seg.end(), pt3.begin(), pt3.end());
	  seg.insert(seg.end(), pt3.begin(), pt3.end());
	  seg.insert(seg.end(), pt1.begin(), pt1.end());
	}
      LineCloud line_seg(&seg[0], (int)seg.size()/6);
      line_seg.writeStandardHeader(out1);
      line_seg.write(out1);
    }
  
  std::ofstream out2("triangles2.g2");
  for (ki=0; ki<meshes2.size(); ++ki)
    {
      double *nodes = meshes2[ki]->vertexArray();
      int num_triang = meshes2[ki]->numTriangles();
      int num_triang_nodes = 3*num_triang;
      unsigned int *triang_idx = meshes2[ki]->triangleIndexArray();
      vector<double> seg;
      for (int kj=0; kj<num_triang_nodes; kj+=3)
	{
	  Point pt1(&nodes[3*triang_idx[kj]],&nodes[3*(triang_idx[kj]+1)]);
	  Point pt2(&nodes[3*triang_idx[kj+1]],&nodes[3*(triang_idx[kj+1]+1)]);
	  Point pt3(&nodes[3*triang_idx[kj+2]],&nodes[3*(triang_idx[kj+2]+1)]);
	  seg.insert(seg.end(), pt1.begin(), pt1.end());
	  seg.insert(seg.end(), pt2.begin(), pt2.end());
	  seg.insert(seg.end(), pt2.begin(), pt2.end());
	  seg.insert(seg.end(), pt3.begin(), pt3.end());
	  seg.insert(seg.end(), pt3.begin(), pt3.end());
	  seg.insert(seg.end(), pt1.begin(), pt1.end());
	}
      LineCloud line_seg(&seg[0], (int)seg.size()/6);
      line_seg.writeStandardHeader(out2);
      line_seg.write(out2);
     }
  
  std::ofstream out3("triangles3.g2");
  for (ki=0; ki<meshes3.size(); ++ki)
    {
      double *nodes = meshes3[ki]->vertexArray();
      int num_triang = meshes3[ki]->numTriangles();
      int num_triang_nodes = 3*num_triang;
      unsigned int *triang_idx = meshes3[ki]->triangleIndexArray();
      vector<double> seg;
      for (int kj=0; kj<num_triang_nodes; kj+=3)
	{
	  Point pt1(&nodes[3*triang_idx[kj]],&nodes[3*(triang_idx[kj]+1)]);
	  Point pt2(&nodes[3*triang_idx[kj+1]],&nodes[3*(triang_idx[kj+1]+1)]);
	  Point pt3(&nodes[3*triang_idx[kj+2]],&nodes[3*(triang_idx[kj+2]+1)]);
	  seg.insert(seg.end(), pt1.begin(), pt1.end());
	  seg.insert(seg.end(), pt2.begin(), pt2.end());
	  seg.insert(seg.end(), pt2.begin(), pt2.end());
	  seg.insert(seg.end(), pt3.begin(), pt3.end());
	  seg.insert(seg.end(), pt3.begin(), pt3.end());
	  seg.insert(seg.end(), pt1.begin(), pt1.end());
	}
      LineCloud line_seg(&seg[0], (int)seg.size()/6);
      line_seg.writeStandardHeader(out3);
      line_seg.write(out3);
     }
    }

  if (cvmodel)
    {
      std::ofstream out1("lines1.g2");
      for (size_t ki=0; ki<meshes.size(); ++ki)
	{
	  LineCloud lines;
	  vector<double> tmp_lines;
	  int nmb_vert = meshes[ki]->numVertices();
	  double *vertices = meshes[ki]->vertexArray();
	  for (int kj=0; kj<nmb_vert-1; ++kj)
	    tmp_lines.insert(tmp_lines.end(), vertices+kj*3, vertices+(kj+2)*3);

	  lines.setCloud(&tmp_lines[0], nmb_vert-1);
	  lines.writeStandardHeader(out1);
	  lines.write(out1);
	}

      std::ofstream out2("lines2.g2");
      for (size_t ki=0; ki<meshes2.size(); ++ki)
	{
	  LineCloud lines;
	  vector<double> tmp_lines;
	  int nmb_vert = meshes2[ki]->numVertices();
	  double *vertices = meshes2[ki]->vertexArray();
	  for (int kj=0; kj<nmb_vert-1; ++kj)
	    tmp_lines.insert(tmp_lines.end(), vertices+kj*3, vertices+(kj+2)*3);

	  lines.setCloud(&tmp_lines[0], nmb_vert-1);
	  lines.writeStandardHeader(out2);
	  lines.write(out2);
	}

     }

  std::vector<shared_ptr<LineCloud> > ctr_pol;
  model->tesselatedCtrPolygon(ctr_pol);

  std::ofstream out4("ctrl_pol.g2");
  for (size_t ki=0; ki<ctr_pol.size(); ++ki)
    {
      ctr_pol[ki]->writeStandardHeader(out4);
      ctr_pol[ki]->write(out4);
    }


  int break_point;
  break_point = 1;
}

