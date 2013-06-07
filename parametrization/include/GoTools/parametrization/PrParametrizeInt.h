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

#ifndef PRPARAMETRIZEINT_H
#define PRPARAMETRIZEINT_H

#include "GoTools/parametrization/PrOrganizedPoints.h"
#include <memory>

/*<PrParametrizeInt-syntax: */

enum PrParamStartVector {
  PrBARYCENTRE            = 1,
  PrFROMUV                = 2
};

/** This class implements an algorithm for creating a
 * parametrization in \f$R^2\f$ of the interior of
 * a given embedding of a planar graph in \f$R^3\f$.
 * The method is described in the paper: M. S. Floater,
 * "Parametrization and smooth approximation of
 * surface triangulations", CAGD 14 (1997), 231-250.
 * The method sets u and v values to each interior
 * node of the graph.
 * The class PrParametrizeInt is an abstract base class.
 * One must call one of its derived classes in order
 * to choose a particular method of parametrization.
 * We recommend the shape preserving parametrization
 * implemented in PrPrmShpPres.
 */
class PrParametrizeInt
{
protected:

  double                 tolerance_;
  PrParamStartVector   startvectortype_;

  shared_ptr<PrOrganizedPoints> g_;

  vector<int>            neighbours_;
  vector<double>           weights_;

  vector< vector<double> > allWeights_;
  vector< vector<int> >    allNeighbours_;

// PRIVATE METHODS

// Used for parametrizing the interior.
  int          getNumIntNghrs(int i);
  int          getNumInt2Nghrs(int i, vector<int>&);
  void         findBarycentre(double& ucentre, double& vcentre);
  bool         isFixed (int, vector<int>&);

  virtual bool makeWeights(int i) = 0;

  const vector< vector<double> >& getAllWeights() const
  {
      return allWeights_;
  }
  const vector< vector<int> >&    getAllNeighbours() const
  {
      return allNeighbours_;
  }

public:
  /// Default constructor
  PrParametrizeInt();
  /// Empty default destructor
  virtual ~PrParametrizeInt();

  /// Set the graph.
  void attach(shared_ptr<PrOrganizedPoints> graph);

  /// Choose how to initialize the solution vector when running the 
  /// internal solver.  Choices are PrBARYCENTRE or PrFROMUV.
  void setStartVectorKind(PrParamStartVector svtype = PrBARYCENTRE)
    {startvectortype_ = svtype;}

  /// Set tolerance for Bi-CGSTAB.
  void setBiCGTolerance(double tolerance = 1.0e-6) {tolerance_ = tolerance;}

  /// Parametrize the given planar graph.
  bool parametrize();

  /** Parametrize the nodes of the 3D graph g_ except those
   * with indices in "fixedPnts". The parameterization is done
   * in 3D and the result is returned as vector "uvw"
   */
  bool parametrize3d(vector<int>&, vector<double>&);

  /** Parametrize the nodes of the 3D graph g_ except those
   * with indices in "fixedPnts". The parameterization is done
   * in 3D and the result is returned as vector "uvw"
   */
  bool new_parametrize3d(vector<int>&, vector<double>&);

  /** This is a simple routine which finds the indices fixedPnts\f$[0...3]\f$
   * of the four vertices of the graph
   * whose (x,y,z) points are the furthest in the directions
   * (-1/-1/1), (1/1/1), (-1/1/-1), and (1/-1/-1) in that order.
   * The indices can be used as fixed points for parametrising in 3D.
   */
  void findFixedPntsFromXYZ(vector<int>& fixedPnts);

  /// performs "nmb" Gauss-Seidel smoothing steps on the sphere
  void smooth(int nmb, vector<int>& fixedPnts); 

  /// computes and stores the weights for each node
  void computeWeights();

};

/*>PrParametrizeInt-syntax: */

/*Class:PrParametrizeInt

Name:              PrParametrizeInt
Syntax:	           @PrParametrizeInt-syntax
Keywords:
Description:       This class implements an algorithm for creating a
                   parametrization in $R^2$ of the interior of
                   a given embedding of a planar graph in $R^3$.
                   The method is described in the paper: M. S. Floater,
                   "Parametrization and smooth approximation of
                   surface triangulations", CAGD 14 (1997), 231-250.
                   The method sets u and v values to each interior
                   node of the graph.
                   The class PrParametrizeInt is an abstract base class.
                   One must call one of its derived classes in order
                   to choose a particular method of parametrization.
                   We recommend the shape preserving parametrization
                   implemented in PrPrmShpPres.
Member functions:
                   "attach(PrOrganizedPoints& graph)" --\\
                   Set the graph.

                   "setBiCGTolerance()" --\\
                   Set tolerance for Bi-CGSTAB.

                   "parametrize()" --\\
                   Parametrize the given planar graph.

		   "smooth(int nmb)" --\\
		   performs "nmb" Gauss-Seidel smoothing steps on the sphere

Constructors:
Files:
Example:

| #include "GoTools/parametrization/PrTriangulation_OP.h"
| #include "GoTools/parametrization/PrPrmShpPres.h"
| 
| main(int argc, char* argv[])
| {
|   // Initialize LA tools because we're using PrParametrizeInt.
|   initDIFFPACK(argc,argv,true);
| 
|   // Read in a triangulation called "gjoevik_triang" which
|   // should be in the directory "src/app/param_pr_triang".
|   // See PrTriangulation_OP.h for examples of the file format.
|   Is infile ("gjoevik_triang", INFILE);
|   PrTriangulation_OP pr_triang;
|   pr_triang.scanRawData(infile);
| 
|   pr_triang.printInfo(s_o);
|   Os xyz_nodes_file("xyz_nodes",NEWFILE);
|   pr_triang.printXYZNodes(xyz_nodes_file);
|   Os xyz_edges_file("xyz_edges",NEWFILE);
|   pr_triang.printXYZEdges(xyz_edges_file);
|   Os xyz_triangles_file("xyz_triangles",NEWFILE);
|   pr_triang.printXYZTriangles(xyz_triangles_file);
| 
|   int no_comps = pr_triang.findNumComponents();
|   int genus = pr_triang.findGenus();
| 
|   if(no_comps == 1 && genus == 1)
|   {
|     PrPrmShpPres p; // PrPrmShpPres is derived from PrParametrizeInt
|     p.setBdyParam(PrCHORDLENGTHBDY); // this is the default anyway
|     p.attach(pr_triang);
| 
|     int cor[4];
|     p.findCornersFromXYZ(cor);
|     s_o << "The corner indices and their xyz values are\n";
|     CgPoint3d pt;
|     for(int k=0; k<4; k++)
|     {
|       pt = pr_triang.get3dNode(cor[k]);
|       s_o << oform("%d  %lf %lf %lf\n",cor[k],pt.x(),pt.y(),pt.z());
|     }
|     p.parametrize(cor[0],cor[1],cor[2],cor[3]);
|   }
| 
|   Os uv_nodes_file("uv_nodes",NEWFILE);
|   pr_triang.printUVNodes(uv_nodes_file);
|   Os uv_edges_file("uv_edges",NEWFILE);
|   pr_triang.printUVEdges(uv_edges_file);
|   Os uv_triangles_file("uv_triangles",NEWFILE);
|   pr_triang.printUVTriangles(uv_triangles_file);
| 
|   Os uvxyz_nodes_file("uvxyz_nodes",NEWFILE);
|   pr_triang.printUVXYZNodes(uvxyz_nodes_file);
|   Os triangulation_file("triangulation",NEWFILE);
|   pr_triang.print(triangulation_file);
| 
|   return 0;
| }

See also:
Developed by:      SINTEF Applied Mathematics, Oslo, Norway
Author:	           Michael Floater, SINTEF
Date:              Mar. 97
*/

#endif // PRPARAMETRIZEINT_H
