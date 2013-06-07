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

#include "GoTools/intersections/Line2DInt.h"


namespace Go {


//===========================================================================
Line2DInt::Line2DInt(Point point, Point dir)
    : AlgObj2DInt(1), point_(point), dir_(dir)
//===========================================================================
{
    // All values in terms_ initialized to 0.0.
    double a, b, c;
    computeConstants(point_, dir_, a, b, c);
    terms_.push_back(Alg2DElem(a, 1, 0));
    terms_.push_back(Alg2DElem(b, 0, 1));
    terms_.push_back(Alg2DElem(c, 0, 0));
}


//===========================================================================
Line2DInt::Line2DInt(double a, double b, double c)
    : AlgObj2DInt(1)
//===========================================================================
{
    // terms_ should be empty.
    terms_.push_back(Alg2DElem(a, 1, 0));
    terms_.push_back(Alg2DElem(b, 0, 1));
    terms_.push_back(Alg2DElem(c, 0, 0));
    computePoints();
}


//===========================================================================
double Line2DInt::a()
//===========================================================================
{
    int ki;
    double a = 0.0;
    int sum_order;
    for (ki = 0; ki < int(terms_.size()); ++ki) {
	sum_order = terms_[ki].degrees_[0] + terms_[ki].degrees_[1];
	ASSERT(sum_order <= 1);
	if (terms_[ki].degrees_[0] == 1)
	    a += terms_[ki].factor_;
    }

    return a;
}


//===========================================================================
double Line2DInt::b()
//===========================================================================
{
    int ki;
    double b = 0.0;
    int sum_order;
    for (ki = 0; ki < int(terms_.size()); ++ki) {
	sum_order = terms_[ki].degrees_[0] + terms_[ki].degrees_[1];
	ASSERT(sum_order <= 1);
	if (terms_[ki].degrees_[1] == 1)
	    b += terms_[ki].factor_;
    }

    return b;
}


//===========================================================================
double Line2DInt::c()
//===========================================================================
{
    int ki;
    double c = 0.0;
    for (ki = 0; ki < int(terms_.size()); ++ki) {
	int sum_order = terms_[ki].degrees_[0] + terms_[ki].degrees_[1];
	ASSERT(sum_order <= 1);
	if (sum_order == 0)
	    c += terms_[ki].factor_;
    }

    return c;
}


//===========================================================================
void Line2DInt::computePoints()
//===========================================================================
{
    // Setting x = 0.0 we get our ref point.
    double a_fac = a();
    double b_fac = b();
    double c_fac = c();
    point_ = (b_fac == 0.0) ? Point(0.0, 0.0) : Point(0.0, -c_fac/b_fac);

    dir_ = Point(-b_fac, a_fac);
    dir_.normalize();
}


//===========================================================================
void Line2DInt::computeConstants(Point point, Point dir,
			       double& a, double& b, double& c)
//===========================================================================
{
    // n*(x - x_0, y - y_0, z - z_0) = 0.
    a = dir_[1];
    b = -dir_[0];
    c = -(a*point_[0] + b*point_[1]);
}


} // namespace Go
