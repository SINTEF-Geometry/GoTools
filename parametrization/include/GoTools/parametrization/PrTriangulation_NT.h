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

#ifndef PRTRIANGULATION_NT_H
#define PRTRIANGULATION_NT_H

#include "GoTools/parametrization/PrNestedTriangulation.h"
#include "GoTools/parametrization/PrLevelTriangulation_OP.h"
#include "GoTools/parametrization/PrParamTriangulation.h"

/*<PrTriangulation_NT-syntax: */

/** PrTriangulation_NT -  This class represents a nested sequence of 
 * triangulations of points in \f$R^3\f$. It implements the virtual functions
 * in PrNestedTriangulations. 
 */
class PrTriangulation_NT : public PrNestedTriangulation
{
private:
  vector<PrNestedNode> node_;
  vector<PrLevelTriangulation_OP*> triang_; // an array of triangulations
        // levels are 0,1,2,...,k
        // so there are k+1 triangulations in total.

  //int getNghrTriangle(int n1, int n2, vector<int>& tlist);
  //void buildTopology();

public:
  /// Default constructor
  PrTriangulation_NT() {}
  /** Constructor. Construct the PrTriangulation_NT from a PrTriangulation_OP,
   * using just one level in the hierarchy. This can later be refined by refine().
   */
  PrTriangulation_NT(PrTriangulation_OP& t);
  /// Empty destructor
  ~PrTriangulation_NT() {}

  /** @name Derived from base class */
  //@{
    virtual int      getFinestLevel() {return (int)triang_.size()-1;}
  virtual int      getNumNodes(int jlev)
                     {return triang_[jlev]->getNumNodes(); }
  virtual double   getX(int i) {return node_[i].x(); }
  virtual double   getY(int i) {return node_[i].y(); }
  virtual double   getZ(int i) {return node_[i].z(); }
  virtual void     setX(int i, const double& x) {node_[i].x() = x; }
  virtual void     setY(int i, const double& y) {node_[i].y() = y; }
  virtual void     setZ(int i, const double& z) {node_[i].z() = z; }
  virtual void     getNeighbours(int i, int jlev, vector<int>& neighbours)
                    {triang_[jlev]->getNeighbours(i,neighbours); }
  virtual bool   isBoundary(int i)
                    {return triang_[node_[i].level()]->isBoundary(i); }
  //@}

  /** @name Other functions */
  //@{

  /** Construct the PrTriangulation_NT from a PrTriangulation_OP,
   * using just one level in the hierarchy.
   * This can later be refined by refine().
   */
  void refine();

  /// Lift the vertices of the nested triangulation to some
  /// surface, using the parameterization defined by pt
  void lift(PrParamTriangulation* pt);

  /// Get pointer to triangulation at level 'i'.
  PrLevelTriangulation_OP* getLevel (int i) { return triang_[i]; }

  /// Get pointer to the finest triangulation
  PrLevelTriangulation_OP* getTopLevel () { return triang_[triang_.size()-1]; }
  //@}
};

/*>PrTriangulation_NT-syntax: */

/*Class:PrTriangulation_NT

Name:              PrTriangulation_NT
Syntax:	           @PrTriangulation_NT-syntax
Keywords:
Description:       This class represents a nested sequence of triangulations
                   of points in R^3.
                   It implements the virtual functions in
                   PrNestedTriangulations.
Member functions:

Constructors:
Files:
Example:


See also:
Developed by:      SINTEF Applied Mathematics, Oslo, Norway
Author:	           Michael Floater, SINTEF
Date:              Nov. 2000
*/

#endif // PRNESTEDPRTRIANGULATION_NT_H
