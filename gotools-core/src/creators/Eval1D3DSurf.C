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




#include "GoTools/creators/Eval1D3DSurf.h"


#include <vector>
#include <assert.h>

using std::vector;
using std::pair;

namespace Go
{

    //===========================================================================
  Eval1D3DSurf::Eval1D3DSurf(shared_ptr<ParamSurface> sf)
    //===========================================================================
        : sf_(sf)
    {
    }


    //===========================================================================
    Eval1D3DSurf::~Eval1D3DSurf()
    //===========================================================================
    {
    }


    //===========================================================================
    Point Eval1D3DSurf::eval(double u, double v) const
    //===========================================================================
    {
        Point pt = sf_->point(u, v);

	Point res;
	if (pt.dimension() == 1)
	  res = Point(u, v, pt[0]);
	else
	  res = pt;
	return res;
    }


    //===========================================================================
    void Eval1D3DSurf::eval(double u, double v, int n, Point der[]) const
    //===========================================================================
    {
        if (n == 0)
        {
	  der[0] = eval(u, v);
	  return;
        }

	vector<Point> pts((n+1)*(n+2)/2);
	sf_->point(pts, u, v, n);
	if (sf_->dimension() == 1)
	  {
	    der[0] = Point(u, v, pts[0][0]);
	    if (n >= 1)
	      {
		der[1] = Point(1, 0, pts[1][0]);
		der[2] = Point(0, 1, pts[2][0]);
		for (size_t ki=3; ki<pts.size(); ++ki)
		  der[ki] = Point(0, 0, pts[ki][0]);
	      }
	  }
	else
	  {
	    for (size_t ki=0; ki<pts.size(); ++ki)
	      der[ki] = pts[ki];
	  }
    }


    //===========================================================================
    double Eval1D3DSurf::start_u() const
    //===========================================================================
    {
        RectDomain rect_dom = sf_->containingDomain();
        double start_u = rect_dom.umin();
        
        return start_u;
    }



    //===========================================================================
    double Eval1D3DSurf::start_v() const
    //===========================================================================
    {
        RectDomain rect_dom = sf_->containingDomain();
        double start_v = rect_dom.vmin();
        
        return start_v;
    }


    //===========================================================================
    double Eval1D3DSurf::end_u() const
    //===========================================================================
    {
        RectDomain rect_dom = sf_->containingDomain();
        double end_u = rect_dom.umax();
        
        return end_u;
    }


    //===========================================================================
    double Eval1D3DSurf::end_v() const
    //===========================================================================
    {
        RectDomain rect_dom = sf_->containingDomain();
        double end_v = rect_dom.vmax();
        
        return end_v;
    }


    //===========================================================================
    int Eval1D3DSurf::dim() const
    //===========================================================================
    {
        int dim = sf_->dimension();
        
	if (dim == 1)
	  return 3;
	else
	  return dim;
    }


    //===========================================================================
    bool Eval1D3DSurf::approximationOK(double par_u, double par_v, Point approxpos,
					     double tol1, double tol2) const
    //===========================================================================
    {

        Point eval_pos = eval(par_u, par_v);
        double dist = eval_pos.dist(approxpos);

        bool appr_ok = (dist < tol1);        
        if (dist > tol1)
        {
            ;//std::cout << "dist: " << dist << std::endl;
        }

        return appr_ok;
    }


} // namespace Go
