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

#ifndef PRPARAMETRIZEBDY_H
#define PRPARAMETRIZEBDY_H

#include "GoTools/parametrization/PrOrganizedPoints.h"
#include <memory>

/*<PrParametrizeBdy-syntax: */


enum PrBdyParamKind {
  PrCHORDLENGTHBDY        = 1,
  PrCENTRIPETAL           = 2,
  PrUNIFBDY               = 3
};

/** This class implements an algorithm for
 * parametrizing the boundary of a planar graph in \f$R^3\f$.
 * Typical choices are uniform and chord length,
 * either round the whole boundary or along each "side".
 */
class PrParametrizeBdy
{
protected:

  PrBdyParamKind       bdyparamtype_;
  shared_ptr<PrOrganizedPoints> g_;
  vector<int>            neighbours_;

// PRIVATE METHODS

  /// Find the bounday node with the least index.
  int          findBdyNode();

  /// Given a boundary node i, return the next boundary node
  /// in an anticlockwise direction around the boundary.
  int          getNextBdyNode(int i);

  ///   Return "length" of chord between two points in \f$R^3\f$.
  double         chord(const Vector3D& a, const Vector3D& b);


public:
  /// Default constructor
  PrParametrizeBdy();
  /// Empty default constructor
  virtual ~PrParametrizeBdy();

  /// Set the graph.
  virtual void attach(shared_ptr<PrOrganizedPoints> graph);

  /** Set type of parametrization to be used along the boundary.
   * Unless you have a strong reason to do so you should
   * set the parametrization to "PrCHORDLENGTHBDY".
   * (Other choices are PrCENTRIPETAL and PrUNIFBDY).
   */
  void setParamKind(PrBdyParamKind bdyparamtype = PrCHORDLENGTHBDY)
                   {bdyparamtype_ = bdyparamtype;}

  /** Calculate the "length" of the section of the xyz boundary
   * between nodes boundary i and j according to the kind of parametrization
   * specified by setParamKind().  Choices are PrCHORDLENGTHBDY, PrCENTRIPETAL
   * and PrUNIFBDY. 
   */
  double         boundaryLength(int i1, int i2);

  /** Parametrize the boundary of the planar graph.
   * The parameter domain will be the unit circle.
   */
  bool         parametrize();

  /** Parametrize the boundary nodes between two boundary nodes i1,i2
   * of a given planar graph by chord length.
   * The two nodes i1 and i2 are assumed to be have u and v values.
   * The boundary parameter points between i1 and i2, starting from i1
   * in an anticlockwise direction around the boundary,
   * will be placed along the straight line between the two parameter points.
   */
    // M.F. Mar. 97.
  bool         parametrizeSide(int i1, int i2);

  /** Parametrize the boundary of the given planar graph.
   * The parameter domain will be the rectangle
   * [umin,umax]*[vmin,vmax]. Node c1 will be
   * mapped to (umin,vmin), c2 to (umax,vmin),
   * c3 to (umax,vmax) and c4 to (umin,vmax).
   * The nodes c1,c2,c3,c4 should be boundary nodes and
   * should be in anticlockwise sequence around the boundary.
   * The kind of parametrization along the four sides
   * of the square will be specified by setParamKind().
   */
  bool         parametrize(int c1, int c2, int c3, int c4,
                           double umin = 0.0, double umax = 1.0,  
                           double vmin = 0.0, double vmax = 1.0);

  /** Find four corner nodes c[0],c[1],c[2],c[3] by taking
   * c[0] to be an arbitrary boundary node and choosing
   * c[1],c[2],c[3] so that the lengths of four boundary
   * curves delimited by c[0],c[1],c[2],c[3] are as equal
   * in "chord" length (defined by setParamKind()) as possible.
   * The array c should be allocated outside with length 4.
   */
  void         findCornersFromXYZ(int* c);

  /** Find the indices c[0],c[1],c[2],c[3] of the four vertices
   * of the graph whose (u,v) points are the furthest
   * SW, SE, NE, and NW in that order.
   * The array c should be allocated outside with length 4.
   */
  void         findCornersFromUV(int* c);
};


/*>PrParametrizeBdy-syntax: */

/*Class:PrParametrizeBdy

Name:              PrParametrizeBdy
Syntax:	           @PrParametrizeBdy-syntax
Keywords:
Description:       This class implements an algorithm for
                   parametrizing the boundary of a planar graph in $R^3$.
                   Typical choices are uniform and chord length,
                   either round the whole boundary or along each "side".
Member functions:
                   "attach(GoHandle<PrOrganizedPoints> graph)" --\\
                   Set the graph.

                   "setParamKind()" --\\
                   Set type of parametrization to be used along the boundary.
                   Unless you have a strong reason to do so you should
                   set the parametrization to "PrCHORDLENGTHBDY".

                   "parametrize()" --\\
                   Parametrize the boundary of the planar graph.
                   The parameter domain will be the unit circle.

                   "parametrize(int& c1, int& c2, int& c3, int& c4)" --\\
                   Parametrize the boundary of the given planar graph.
                   The parameter domain will be the rectangle
                   [umin,umax]*[vmin,vmax]. Node c1 will be
                   mapped to (umin,vmin), c2 to (umax,vmin),
                   c3 to (umax,vmax) and c4 to (umin,vmax).
                   The nodes c1,c2,c3,c4 should be boundary nodes and
                   should be in anticlockwise sequence around the boundary.
                   The kind of parametrization along the four sides
                   of the square will be specified by setParamKind().

                   "findCornersFromXYZ(int* c)" --\\
                   Find four corner nodes c[0],c[1],c[2],c[3] by taking
                   c[0] to be an arbitrary boundary node and choosing
                   c[1],c[2],c[3] so that the lengths of four boundary
                   curves delimited by c[0],c[1],c[2],c[3] are as equal
                   in "chord" length (defined by setParamKind()) as possible.
                   The array c should be allocated outside with length 4.

                   "findCornersFromUV(int* c)" --\\
                   Find the indices c[0],c[1],c[2],c[3] of the four vertices
                   of the graph whose (u,v) points are the furthest
                   SW, SE, NE, and NW in that order.
                   The array c should be allocated outside with length 4.

Constructors:
Files:
Example:


See also:
Developed by:      SINTEF Applied Mathematics, Oslo, Norway
Author:	           Michael Floater, SINTEF
Date:              Nov. 98 (modified from previous PrParametrizer).
*/

#endif // PRPARAMETRIZEBDY_H
