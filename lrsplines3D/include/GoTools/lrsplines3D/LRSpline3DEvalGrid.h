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

#ifndef _LRSPLINE3DEVALGRID_H
#define _LRSPLINE3DEVAVGRID_H


#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/lrsplines3D/Element3D.h"
#include "GoTools/lrsplines3D/Mesh3D.h"

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

typedef Element3D simpleElement3D;


  /// Grid evaluation of the elements of an LR B-spline volume.
// =============================================================================
class LRSpline3DEvalGrid
// =============================================================================
{
public:

  // Constructor
  LRSpline3DEvalGrid();

  typedef double param_float_type;
  
  // LRSpline3DEvalGrid(int numElements);

  /// Constructor given an LR B-spline volume
  LRSpline3DEvalGrid(LRSplineVolume& lr_spline);

  /// Iterator to the first element of the LR B-spline volume
  std::vector<Element3D>::iterator elements_begin()// const
    {
      return elements_.begin();
    }

  /// Iterator past the last element of the LR B-spline volume
  std::vector<Element3D>::iterator elements_end()// const
    {
      return elements_.end();
    }

  /// Perform grid evaluation
    template <class V>
      void evaluateGrid(V &element, double *points)
      {
	int order_U = orderU();
	int order_V = orderV();
	int order_W = orderW();
	double ll_x=element.umin(), ll_y=element.vmin(), ll_z=element.wmin();
	double ur_x=element.umax(), ur_y=element.vmax(), ur_z=element.wmax();
	/* low(element, ll_x, ll_y, ll_z); */
	/* high(element, ur_x, ur_y, ur_z); */
	
	auto du = (ur_x - ll_x) / (order_U - 1);
	auto dv = (ur_y - ll_y) / (order_V - 1);
	auto dw = (ur_z - ll_z) / (order_W - 1);

	std::fill(points, points+order_U*order_V*order_W*dim_, 0.0);
	const std::vector<LRBSpline3D*> bfunctions = element.getSupport();
	size_t bsize = bfunctions.size();

	// Evaluate univariate B-splines
	std::vector<double> val1(bsize*order_U);
	std::vector<double> val2(bsize*order_V);
	std::vector<double> val3(bsize*order_W);
	for (size_t ki=0; ki<bsize; ++ki)
	  {
	    size_t kj;
	    const BSplineUniLR* uni1 =  bfunctions[ki]->getUnivariate(XDIR);
	    const BSplineUniLR* uni2 =  bfunctions[ki]->getUnivariate(YDIR);
	    const BSplineUniLR* uni3 =  bfunctions[ki]->getUnivariate(ZDIR);
	    for (kj=0; kj<ki; ++kj)
	      if (uni1 == bfunctions[kj]->getUnivariate(XDIR))
		break;
	    if (kj < ki)
	      std::copy(val1.begin()+kj*order_U, val1.begin()+(kj+1)*order_U,
			val1.begin()+ki*order_U);
	    else
	      {
		double par = ll_x;
		for (int ka=0; ka<order_U; ++ka, par+=du)
		  {
		    if (ka == order_U-1)
		      par = ur_x;
		    const bool on_end = (par == orig_dom_[1]); //bfunctions[ki]->umax());
		    val1[ki*order_U+ka] = uni1->evalBasisFunction(par, 0, on_end);
		  }
	      }
		
	    for (kj=0; kj<ki; ++kj)
	      if (uni2 == bfunctions[kj]->getUnivariate(ZDIR))
		break;
	    if (kj < ki)
	      std::copy(val2.begin()+kj*order_V, val2.begin()+(kj+1)*order_V,
			val2.begin()+ki*order_V);
	    else
	      {
		double par = ll_y;
		for (int ka=0; ka<order_V; ++ka, par+=dv)
		  {
		    if (ka == order_V-1)
		      par = ur_y;
		    const bool on_end = (par == orig_dom_[3]); //bfunctions[ki]->vmax());
		    val2[ki*order_V+ka] = uni2->evalBasisFunction(par, 0, on_end);
		  }
	      }
		
	    for (kj=0; kj<ki; ++kj)
	      if (uni3 == bfunctions[kj]->getUnivariate(ZDIR))
		break;
	    if (kj < ki)
	      std::copy(val3.begin()+kj*order_W, val3.begin()+(kj+1)*order_W,
			val3.begin()+ki*order_W);
	    else
	      {
		double par = ll_z;
		for (int ka=0; ka<order_W; ++ka, par+=dw)
		  {
		    if (ka == order_W-1)
		      par = ur_z;
		    const bool on_end = (par == orig_dom_[5]); //bfunctions[ki]->wmax());
		    val3[ki*order_W+ka] = uni3->evalBasisFunction(par, 0, on_end);
		  }
	      }
	  }

	Point pt;
	double *p=points;
	int ix2=0;
	for (int i=0; i<order_W; ++i)
	  for (int j=0; j<order_V; ++j)
	    for (int k=0; k<order_U; ++k, p+=dim_)
	      for (size_t ix=0; ix<bsize; ++ix)
		{
		  pt = bfunctions[ix]->coefTimesGamma();
		  for (int ka=0; ka<dim_; ++ka)
		    p[ka] +=
		      pt[ka]*val1[ix*order_U+k]*
		      val2[ix*order_V+j]*
		      val3[ix*order_W+i];
		  //(*p) += pt[ka]*val1[k*bsize+ix]*val2[j*bsize+ix]*val3[i*bsize+ix];
		}
		  
	    /* auto w=ll_z; */
	    /* for (int i=0; i<order_W; i++, w += dw) { */
            /*   if (i==order_W-1) { */
            /*     w=ur_z; */
            /*   } */
            /*   auto v=ll_y; */
            /*   for (int j=0; j<order_V; j++, v += dv) { */
	    /* 	if (j==order_V-1) { */
            /*       v=ur_y; */
	    /* 	} */
	    /* 	auto u=ll_x; */
	    /* 	for (int k=0; k<order_U; k++, u += du) { */
            /*       if (k==order_U-1) { */
            /*         u=ur_x; */
            /*       } */
            /*       evaluate(element, u, v, w, p); */
            /*       p+=dim_; */
	    /* 	} */
            /*   } */
            /* } */
      }

  // It appears that this function expects u,v,w \in [0,1]
  void evaluate(Element3D &elem, double u, double v, double w, double *res) const
    {
      double scaledU = u *(orig_dom_[1]-orig_dom_[0]);
      scaledU += orig_dom_[0];
      double scaledV = v *(orig_dom_[3]-orig_dom_[2]);
      scaledV += orig_dom_[2];
      double scaledW = w *(orig_dom_[5]-orig_dom_[4]);
      scaledW += orig_dom_[4];
      
      //assert(dim_ == 1 || dim_ == 3);

      Point result(dim_);
      result.setValue(0.0);

      const std::vector<LRBSpline3D*> covering_B_functions = elem.getSupport();

      for (auto b = covering_B_functions.begin();
	   b != covering_B_functions.end(); ++b)
	{
	  const bool u_on_end = (scaledU == (*b)->umax());
	  const bool v_on_end = (scaledV == (*b)->vmax());
	  const bool w_on_end = (scaledW == (*b)->wmax());

	  result += (*b)->eval(scaledU, 
			       scaledV, 
                               scaledW, 
			       0, // No derivs.
			       0, // No derivs.
                               0, // No derivs.
			       u_on_end, 
                               v_on_end, 
			       w_on_end);
	}

      for (int ix=0; ix!=dim_; ++ix) {
        res[ix] = result[ix];
      }
    }

  int numElements() const
    {
      return elements_.size();
    }

  int dim() const
    {
      return dim_; // If dim is 1 we use the parameter domain as the first two dimensions.
    }

  /// Order (polynomial degree + 1) of the volume in the first parameter direction
  int orderU() const
    {
      return order_u_;
    }

  /// Order (polynomial degree + 1) of the volume in the second parameter direction
  int orderV() const
    {
      return order_v_;
    }

  /// Order (polynomial degree + 1) of the volume in the third parameter direction
  int orderW() const
    {
      return order_w_;
    }

  /// Lower left front corner of element
  void low(const Element3D &e, double &u, double &v, double &w)
    {
      u = e.umin();
      v = e.vmin();
      w = e.wmin();
      u -= orig_dom_[0];
      u /= orig_dom_[1]-orig_dom_[0];
      v -= orig_dom_[2];
      v /= orig_dom_[3]-orig_dom_[2];
      w -= orig_dom_[4];
      w /= orig_dom_[5]-orig_dom_[4];
    }

  /// Upper right back corner of element
  void high(const Element3D &e, double &u, double &v, double &w)
    {
      u = e.umax();
      v = e.vmax();
      w = e.wmax();
      u -= orig_dom_[0];
      u /= orig_dom_[1]-orig_dom_[0];
      v -= orig_dom_[2];
      v /= orig_dom_[3]-orig_dom_[2];
      w -= orig_dom_[4];
      w /= orig_dom_[5]-orig_dom_[4];
    }

    // Copy and paste from code in r2gl.
    //void testCoefComputation();
  Array<double,6> origDom() {
    return orig_dom_;
  }


private:
    // umin, umax, vmin, vmax, wmin, wmax
    Array<double,6> orig_dom_;
    std::vector<Element3D> elements_;
    int order_u_;
    int order_v_;
    int order_w_;
    int dim_;
    Mesh3D mesh_;

    /*
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
        }*/

};

} // end namespace Go


#endif // _LRSPLINE3DEVALGRID_H

