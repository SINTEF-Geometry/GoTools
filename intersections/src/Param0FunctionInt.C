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

#include "GoTools/intersections/Param0FunctionInt.h"
#include "GoTools/utils/CompositeBox.h"


using std::vector;


namespace Go {


//===========================================================================
Param0FunctionInt::Param0FunctionInt(double C)
    : C_(C)
//===========================================================================
{
    dim_ = 0;
    parentcurve_ = 0;
}


//===========================================================================
Param0FunctionInt::Param0FunctionInt(double C,
				     ParamFunctionInt *parent)
    : C_(C)
//===========================================================================
{
    dim_ = 0;
    parentcurve_ = parent;
}


//===========================================================================
Param0FunctionInt::~Param0FunctionInt()
//===========================================================================
{
}


//===========================================================================
Param0FunctionInt* Param0FunctionInt::getParam0FunctionInt()
//===========================================================================
{ 
    return this;
}


//===========================================================================
int Param0FunctionInt::numParams() const
//===========================================================================
{
    return 0;
}


//===========================================================================
void Param0FunctionInt::getLengthAndWiggle(double *length, double *wiggle)
//===========================================================================
{
    *length = 0.0;
    *wiggle = 0.0; // Does not really apply to a point.
}


//===========================================================================
bool Param0FunctionInt::hasInnerKnots(int pardir) const
//===========================================================================
{
    return false;
}


//===========================================================================
vector<double> Param0FunctionInt::
getInnerKnotVals(int pardir, bool sort) const
//===========================================================================
{
    vector<double> knots;
    return knots;
}


//===========================================================================
bool Param0FunctionInt::hasCriticalVals(int pardir) const
//===========================================================================
{
    return false;
}


//===========================================================================
vector<double> Param0FunctionInt::getCriticalVals(int pardir) const
//===========================================================================
{
    vector<double> vals;
    return vals;
}


//===========================================================================
bool Param0FunctionInt::hasCriticalValsOrKnots(int pardir) const
//===========================================================================
{
    return false;
}


//===========================================================================
vector<double> Param0FunctionInt::getCriticalValsAndKnots(int pardir) const
//===========================================================================
{
    vector<double> vals;
    return vals;
}


//===========================================================================
bool Param0FunctionInt::canDivide(int pardir)
//===========================================================================
{
    return false;
}


//===========================================================================
int Param0FunctionInt::getMeshSize(int dir)
//===========================================================================
{
    return 0;
}


//===========================================================================
double Param0FunctionInt::paramFromMesh(int dir, int idx)
//===========================================================================
{
    THROW("Not applicable!");
}


//===========================================================================
vector<double>::iterator Param0FunctionInt::getMesh()
//===========================================================================
{
    vector<double> mesh;
    return mesh.begin();
}


//===========================================================================
double Param0FunctionInt::startParam(int pardir) const
//===========================================================================
{
    return 0.0;
}


//===========================================================================
double Param0FunctionInt::endParam(int pardir) const
//===========================================================================
{
    return 0.0;
}


//===========================================================================
bool Param0FunctionInt::boundaryPoint(const double* par, double eps) const
//===========================================================================
{
    return false; // Object not bounded.
}


//===========================================================================
void Param0FunctionInt::
subdivide(int pardir, double par, 
	  vector<shared_ptr<ParamFunctionInt> >& subdiv_objs,
	  vector<shared_ptr<ParamFunctionInt> >& bd_objs)
//===========================================================================
{
    return;
}


// //==========================================================================
// shared_ptr<Param0FunctionInt> 
//    Param0FunctionInt::makeIntFunction(double C)
// //==========================================================================
// {
//   shared_ptr<Param0FunctionInt> const_int =
//       shared_ptr<Param0FunctionInt>(new Param0FunctionInt(C, this));
//   return const_int;
// }


//===========================================================================
CompositeBox Param0FunctionInt::compositeBox() const
//===========================================================================
{
    Point c_pt(1);
    c_pt[0] = C_;
    CompositeBox box(c_pt, c_pt);
    return box;
}


//===========================================================================
bool Param0FunctionInt::monotone(Point& dir, double tol) const
//===========================================================================
{
    return false; // I guess the definition does not apply to a point.
}


//===========================================================================
void Param0FunctionInt::
getBoundaryObjects(vector<shared_ptr<BoundaryFunctionInt> >& bd_objs)
//===========================================================================
{
    return;
}


//===========================================================================
int Param0FunctionInt::dimension()
//===========================================================================
{
    return dim_;
}


//===========================================================================
void Param0FunctionInt::point(Point& res, const double *par) const
//===========================================================================
{
    res.resize(1);
    res[0] = C_;
}


// //=========================================================================
// void Param0FunctionInt::assureInRange(double& t)
// //=========================================================================
// {
//     t = 0.0;
// }


//===========================================================================
int Param0FunctionInt::knotIntervalFuzzy(double& t, double tol) const
//===========================================================================
{
    if (fabs(t) < tol)
	t = 0.0; // Well, we really do not care about the parameter.
    return -1;
}


//===========================================================================
double Param0FunctionInt::nextSegmentVal(double par, bool forward) const
//===========================================================================
{
    return 0.0;
}


} // namespace Go

