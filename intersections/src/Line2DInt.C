//===========================================================================
//                                                                           
// File: Line2DInt.C                                                           
//                                                                           
// Created: Mon Jan 24 12:28:23 2005                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: Line2DInt.C,v 1.3 2006-03-06 12:11:41 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


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
