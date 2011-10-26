/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 2000 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

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
