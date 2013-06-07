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
