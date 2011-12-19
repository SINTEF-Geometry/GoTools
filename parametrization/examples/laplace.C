/*----------------------------------------------*/
/* Test program for testing approximations to solutions
   to Laplace equation u_xx + u_yy = 0 over the unit square
   with Dirichlet boundary condition equal to he function xy.
   We triangluate the unit square with a Type I triangulation.
   The FE method in this case reduces to the finite difference
   method: the (4, -1, -1, -1, -1) ifive-point stencil.
   The mean value approximation on the other hand has 6 non-zero weights,
   so is a seven-point stencil.
   The "uv"' of the mesh parameterization are used to store
   the height values (as u). v is set to zero.
   Michael Floater, April 2002.    */
/*----------------------------------------------*/

#include "GoTools/parametrization/PrTriangulation_OP.h"
#include "GoTools/parametrization/PrParametrizeInt.h"
#include "GoTools/parametrization/PrParametrizeBdy.h"
#include "GoTools/parametrization/PrPrmShpPres.h"
//#include "GoTools/parametrization/PrPrmLShape.h"
#include "GoTools/parametrization/PrPrmUniform.h"
#include "GoTools/parametrization/PrPrmLeastSquare.h"
#include "GoTools/parametrization/PrPrmEDDHLS.h"
#include "GoTools/parametrization/PrPrmMeanValue.h"
#include <memory>
#include <cstring>
#include <fstream>


using std::cout;
using std::endl;
using std::strcmp;


int main(int argc, const char** argv)
{
  // char inf[80];
  //std::string inf2;
  int intparam_type = 1;
  int bdyparam_type = 1;
  int bdy_method = 1;

  // For choosing grid. Number of squares in [0,1] * [0,1].
  int m_grid = 50; 
  int n_grid = 50; 
  int undefCmd = 0;

  for (int i=1; i<argc; i=i+2)
   {
     if( strcmp(argv[i],"-intpar") == 0 )
            intparam_type = atoi(argv[i+1]);
     else if( strcmp(argv[i],"-bdypar") == 0 )
            bdyparam_type = atoi(argv[i+1]);
     else if( strcmp(argv[i],"-bdymeth") == 0 )
            bdy_method = atoi(argv[i+1]);
     else if( strcmp(argv[i],"-m_grid") == 0 )
            m_grid = atoi(argv[i+1]);
     else if( strcmp(argv[i],"-n_grid") == 0 )
            n_grid = atoi(argv[i+1]);
     else
      undefCmd = 1;
   }

  if(undefCmd)
  {
    cout << "Error in command line" << endl;
    return -1;
  }

  // Make type II triangulation of unit square,
  // dividing into m_grid * n_grid rectangles.
  // Put boundary of uv points in too (the rest are dummies
  // to be solved for later).

  int np = (m_grid+1) * (n_grid+1);
  double *xyz_points = new double[3*np];
  double *uv_points = new double[2*np];
  int nt = 2 * m_grid * n_grid;
  int *triangles = new int[3*nt];

  for(int j=0; j<=n_grid; j++)
    for(int i=0; i<=m_grid; i++)
    {
      int ii = j*(m_grid+1) + i;
      xyz_points[3*ii] = (double)i / (double)m_grid;
      xyz_points[3*ii+1] = (double)j / (double)n_grid;
      xyz_points[3*ii+2] = 0.0;
      uv_points[2*ii] = xyz_points[3*ii] * xyz_points[3*ii+1];
      uv_points[2*ii+1] = 0.0;
    }

  int k=0;
  for(int j=0; j<n_grid; j++)
    for(int i=0; i<m_grid; i++)
    {
      int ii = j*(m_grid+1) + i;
      triangles[k] = ii;
      triangles[k+1] = ii+1;
      triangles[k+2] = ii+m_grid+2;
      triangles[k+3] = ii;
      triangles[k+4] = ii+m_grid+2;
      triangles[k+5] = ii+m_grid+1;
      k += 6;
    }

  // Create the triangulation class.
  shared_ptr<PrTriangulation_OP>
    pr_triang(new PrTriangulation_OP(xyz_points,uv_points,np,triangles,nt));
  pr_triang->printInfo(cout);

  std::ofstream pout("triang.out");
  pr_triang->printRawData(pout);

  std::ofstream xyz_nodes_file("xyz_nodes");
  pr_triang->printXYZNodes(xyz_nodes_file);
  std::ofstream xyz_edges_file("xyz_edges");
  pr_triang->printXYZEdges(xyz_edges_file);
  std::ofstream xyz_triangles_file("xyz_triangles");
  pr_triang->printXYZTriangles(xyz_triangles_file);

  std::ofstream xyz_face_file("faces.m");
  pr_triang->printXYZFacesML(xyz_face_file);

  int no_comps = pr_triang->findNumComponents();
  int genus = pr_triang->findGenus();

  if(no_comps == 1 && genus == 1)
  {
    cout << "Parametrizing interior..." << endl;

    PrParametrizeInt *pi;
    switch(intparam_type)
    {
      case 1: pi = new PrPrmShpPres; break;
      case 2: pi = new PrPrmUniform; break;
      case 3: pi = new PrPrmLeastSquare; break;
      case 4: pi = new PrPrmEDDHLS; break;
      // case 5: pi = new PrPrmLShape; break;
      case 6: pi = new PrPrmMeanValue; break;
    }
  
    pi->attach(pr_triang);
    pi->parametrize();
    delete pi;

    std::ofstream u_nodes_file("u_nodes");
    u_nodes_file << n_grid + 1 << " " << m_grid + 1 << std::endl;
    std::ofstream error_file("error");
    double max_error = 0.0;
    for(int j=0; j<=n_grid; j++)
    {
      for(int i=0; i<=m_grid; i++)
      {
        int ii = j*(m_grid+1) + i;
        u_nodes_file << pr_triang->getU(ii) << " ";
	if (i == m_grid) u_nodes_file << '\n';
        double error = pr_triang->getU(ii) - 
                          ( (double)i / (double)m_grid
                            * (double)j / (double)n_grid );
        if(fabs(error) > max_error) max_error = fabs(error);
          // true solution is u(x,y) = xy.
        error_file << error << "\n";
      } 
      u_nodes_file << "\n";
      error_file << "\n";
    }

    cout << "Max error = " << max_error << "\n";
  }
  else return 0;
}
