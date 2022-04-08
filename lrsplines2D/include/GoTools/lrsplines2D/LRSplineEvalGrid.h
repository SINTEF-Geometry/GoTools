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

#ifndef _LRSPLINEEVALGRID_H
#define _LRSPLINEEVALGRID_H


#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/Element2D.h"
#include "GoTools/lrsplines2D/Mesh2D.h"

#include <vector>


namespace Go
{


// struct simpleElement {
// glm::vec2 &low() {return m_low;}
// const glm::vec2 &low() const {return m_low;}
// glm::vec2 &high() {return m_high;}
// const glm::vec2 &high() const {return m_high;}
// glm::vec2 m_low;
// glm::vec2 m_high;
// simpleElement(glm::vec2 newLow,	glm::vec2 newHigh) :
// m_low(newLow), m_high(newHigh) {;}
// };

typedef Element2D simpleElement;


  /// Grid evaluation of the elements of an LR B-spline surface.
// =============================================================================
class LRSplineEvalGrid
// =============================================================================
{
public:

  // Constructor
  LRSplineEvalGrid();

  typedef double param_float_type;

//    LRSplineEvalGrid(int numElements);

  /// Constructor given an LR B-spline surface
  LRSplineEvalGrid(LRSplineSurface& lr_spline);

  /// Iterator to the first element of the LR B-spline surface
  std::vector<Element2D>::iterator elements_begin()// const
    {
      return elements_.begin();
    }

  /// Iterator past the last element of the LR B-spline surface
  std::vector<Element2D>::iterator elements_end()// const
    {
      return elements_.end();
    }

  void evaluate(Element2D &elem, double u, double v, double *res) const
    {
		double scaledU = u *(orig_dom_.umax()-orig_dom_.umin());
		scaledU += orig_dom_.umin();
		double scaledV = v *(orig_dom_.vmax()-orig_dom_.vmin());
		scaledV += orig_dom_.vmin();

		assert(dim_ ==1 || dim_ == 3);

      Point result(dim_);
      result.setValue(0.0);

      const std::vector<LRBSpline2D*> covering_B_functions = elem.getSupport();

      for (auto b = covering_B_functions.begin();
	   b != covering_B_functions.end(); ++b)
	{
	  const bool u_on_end = (scaledU == (*b)->umax());
	  const bool v_on_end = (scaledV == (*b)->vmax());

	  result += (*b)->eval(scaledU, 
			       scaledV, 
			       0, // No derivs.
			       0, // No derivs.
			       u_on_end, 
			       v_on_end);
	}

      if (dim_ == 3)
	{
#if 1
	  res[0] = result[0];
	  res[1] = result[1];
      res[2] = result[2];
//      std::cout << "res: " << result << std::endl;
#else // @@sbr201301 Setting first two params to parameter domain.
	  res[0] = u;
	  res[1] = v;
      res[2] = 0.001*result[2];
#endif
	}
      else
	{
	  res[0] = scaledU;
	  res[1] = scaledV;
	  res[2] = result[0];
	  //std::cout << "res: " << result << std::endl;
	}
    }

  int numElements() const
    {
      return elements_.size();
    }

  int dim() const
    {
      return 3; // If dim is 1 we use the parameter domain as the first two dimensions.
    }

  /// Order (polynomial degree + 1) of the surface in the first parameter direction
  int orderU() const
    {
      return order_u_;
    }

  /// Order (polynomial degree + 1) of the surface in the second parameter direction
  int orderV() const
    {
      return order_v_;
    }

  /// Lower left corner of element
  void low(const Element2D &e, double &u, double &v)
    {
      u = e.umin();
      v = e.vmin();
	  u -= orig_dom_.umin();
	  u /= orig_dom_.umax()-orig_dom_.umin();
	  v -= orig_dom_.vmin();
	  v /= orig_dom_.vmax()-orig_dom_.vmin();
    }

  /// Upper right corner of element
  void high(const Element2D &e, double &u, double &v)
    {
      u = e.umax();
      v = e.vmax();
	  u -= orig_dom_.umin();
	  u /= orig_dom_.umax()-orig_dom_.umin();
	  v -= orig_dom_.vmin();
	  v /= orig_dom_.vmax()-orig_dom_.vmin();
    }

    // Copy and paste from code in r2gl.
    void testCoefComputation();

  /// Compute Bezier coefficients of element based on a grid of sample points obtained
  /// by evaluating the surface in the element. The number of sample points depends on
  /// the polynomial degree
    inline void computeBezCoefs(int dim, const double *points, int orderU, int orderV, double *coefs)
	{
	    if ((orderU != 3 && orderU != 4) || (orderV != 3 && orderV != 4)) {
		throw std::runtime_error("LRViz only supports quadratic and cubic surfs");
	    }
	    double scratchVec[4*4*4]; // maxOrderU*MaxOrderV*maxDim
	    double M3[9] = // interpolation matrix for quadratic splines
		{
		    1, 0, 0,
		    -.5, 2, -.5,
		    0, 0, 1
		};
	    double M4[16] = // interpolation matrix for cubic splines
		{
		    1.0, 0, 0, 0,
		    -5/6.0, 18/6.0, -9/6.0, 2/6.0,
		    2/6.0, -9/6.0, 18/6.0, -5/6.0,
		    0, 0, 0, 1
		};
	    double *p;
	    if (orderU == 3) {
		for(int row=0; row<orderV; row++) {
		    for(int i=0; i<dim; i++) {
			scratchVec[(3*row+0)*dim+i] = points[(3*row+0)*dim+i];
			scratchVec[(3*row+1)*dim+i] = M3[3]*points[(3*row+0)*dim+i]+M3[4]*points[(3*row+1)*dim+i]+M3[5]*points[(3*row+2)*dim+i];
			scratchVec[(3*row+2)*dim+i] = points[(3*row+2)*dim+i];
		    }
		}
	    } else {
		for(int row=0; row<orderV; row++) {
		    for(int i=0; i<dim; i++) {
			scratchVec[(4*row+0)*dim+i] = points[(4*row+0)*dim+i];
			scratchVec[(4*row+1)*dim+i] = M4[4]*points[(4*row+0)*dim+i]+M4[5]*points[(4*row+1)*dim+i]+M4[6]*points[(4*row+2)*dim+i]+M4[7]*points[(4*row+3)*dim+i];
			scratchVec[(4*row+2)*dim+i] = M4[8]*points[(4*row+0)*dim+i]+M4[9]*points[(4*row+1)*dim+i]+M4[10]*points[(4*row+2)*dim+i]+M4[11]*points[(4*row+3)*dim+i];
			scratchVec[(4*row+3)*dim+i] = points[(4*row+3)*dim+i];
		    }
		}
	    }
	    if (orderV == 3) {
		for(int column=0; column<orderU; column++) {
		    for(int i=0; i<dim; i++) {
			coefs[(0*orderU+column)*dim+i] = scratchVec[(0*orderU+column)*dim+i];
			coefs[(1*orderU+column)*dim+i] = M3[3]*scratchVec[(0*orderU+column)*dim+i]+M3[4]*scratchVec[(1*orderU+column)*dim+i]+M3[5]*scratchVec[(2*orderU+column)*dim+i];
			coefs[(2*orderU+column)*dim+i] = M3[8]*scratchVec[(2*orderU+column)*dim+i];
		    }
		}
	    } else {
		for(int column=0; column<orderU; column++) {
		    for(int i=0; i<dim; i++) {
			coefs[(0*orderU+column)*dim+i] = scratchVec[(0*orderU+column)*dim+i];
			coefs[(1*orderU+column)*dim+i] = M4[ 4]*scratchVec[(0*orderU+column)*dim+i]+M4[5]*scratchVec[(1*orderU+column)*dim+i]+M4[6 ]*scratchVec[(2*orderU+column)*dim+i]+M4[7 ]*scratchVec[(3*orderU+column)*dim+i];
			coefs[(2*orderU+column)*dim+i] = M4[ 8]*scratchVec[(0*orderU+column)*dim+i]+M4[9]*scratchVec[(1*orderU+column)*dim+i]+M4[10]*scratchVec[(2*orderU+column)*dim+i]+M4[11]*scratchVec[(3*orderU+column)*dim+i];
			coefs[(3*orderU+column)*dim+i] = M4[15]*scratchVec[(3*orderU+column)*dim+i];
		    }
		}
	    }
	}

    template <class V>
    void evaluateGrid(V &element, double *points)
	{
	    int order_U = orderU();
	    int order_V = orderV();
	    double ll_x, ll_y;
	    double ur_x, ur_y;
// float ll_x, ll_y;
// float ur_x, ur_y;
	    low(element, ll_x, ll_y);
	    high(element, ur_x, ur_y);
	    double *p=points;
	    auto v=ll_y;
	    auto du = (ur_x - ll_x) / (order_U - 1);
	    auto dv = (ur_y - ll_y) / (order_V - 1);
	    for (int i=0; i<order_V; i++, v += dv) {
		if (i==order_V-1) {
		    v=ur_y;
		}
		auto u=ll_x;
		for (int j=0; j<order_U; j++, u += du) {
		    if (j==order_U-1) {
			u=ur_x;
		    }
		    evaluate(element, u, v, p);
		    p+=dim_;
		}
	    }
	}

    RectDomain origDom() {
      return orig_dom_;
    }


private:
	RectDomain orig_dom_;
  std::vector<Element2D> elements_;
  int order_u_;
  int order_v_;
  int dim_;
  Mesh2D mesh_;


};

} // end namespace Go


#endif // _LRSPLINEEVALGRID_H

