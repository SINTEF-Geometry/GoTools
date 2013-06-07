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
#define _LRSPLINEEVAVGRID_H


#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/Element2D.h"

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


// =============================================================================
class LRSplineEvalGrid
// =============================================================================
{
public:

  LRSplineEvalGrid();

  typedef double param_float_type;

//    LRSplineEvalGrid(int numElements);

  LRSplineEvalGrid(LRSplineSurface& lr_spline);

  std::vector<Element2D>::iterator elements_begin()// const
    {
      return elements_.begin();
    }

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
#else // @@sbr201301 Setting first two params to parameter domain.
	  res[0] = u;
	  res[1] = v;
      res[2] = 0.001*result[2];
#endif
	}
      else
	{
	  res[0] = u;
	  res[1] = v;
	  res[2] = result[0];
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

  int orderU() const
    {
      return order_u_;
    }

  int orderV() const
    {
      return order_v_;
    }

  void low(const Element2D &e, double &u, double &v)
    {
      u = e.umin();
      v = e.vmin();
	  u -= orig_dom_.umin();
	  u /= orig_dom_.umax()-orig_dom_.umin();
	  v -= orig_dom_.vmin();
	  v /= orig_dom_.vmax()-orig_dom_.vmin();
    }

  void high(const Element2D &e, double &u, double &v)
    {
      u = e.umax();
      v = e.vmax();
	  u -= orig_dom_.umin();
	  u /= orig_dom_.umax()-orig_dom_.umin();
	  v -= orig_dom_.vmin();
	  v /= orig_dom_.vmax()-orig_dom_.vmin();
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

