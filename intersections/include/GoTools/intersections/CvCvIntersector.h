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

#ifndef _CVCVINTERSECTOR_H
#define _CVCVINTERSECTOR_H


#include "GoTools/intersections/Intersector2Obj.h"


namespace Go {


/// Class that performs intersection between two parametric curves.

class CvCvIntersector : public Intersector2Obj {
public:

    /// Constructor.
    /// The last two variables pertains only if the parent has one
    /// more parameter than the Intersector to be constructed.
    /// \param curve1 of type ParamCurveInt.
    /// \param curve2 of type ParamCurveInt.
    /// \param epsge the associated tolerance.
    /// \param prev the "parent" Intersector (0 if there is no
    /// parent).
    /// \param eliminated_parameter the index of the parameter that
    /// was removed from the parent \a prev.
    /// \param eliminated_value the value of the parameter that was
    /// removed from the parent \a prev.
    CvCvIntersector(shared_ptr<ParamGeomInt> curve1, 
		    shared_ptr<ParamGeomInt> curve2,
		    double epsge,
		    Intersector* prev = 0,
		    int eliminated_parameter = -1,
		    double eliminated_value = 0);

    /// Constructor.
    /// The last two variables pertains only if the parent has one
    /// more parameter than the Intersector to be constructed.
    /// \param curve1 of type ParamCurveInt.
    /// \param curve2 of type ParamCurveInt.
    /// \param epsge the associated tolerance.
    /// \param prev the "parent" Intersector (0 if there is no
    /// parent).
    /// \param eliminated_parameter the index of the parameter that
    /// was removed from the parent \a prev.
    /// \param eliminated_value the value of the parameter that was
    /// removed from the parent \a prev.
    CvCvIntersector(shared_ptr<ParamGeomInt> curve1, 
		    shared_ptr<ParamGeomInt> curve2,
		    shared_ptr<GeoTol> epsge, 
		    Intersector* prev = 0,
		    int eliminated_parameter = -1,
		    double eliminated_value = 0);

    /// Destructor.
    virtual ~CvCvIntersector();

//     // Validation of given intersection results
//     virtual void validate(int level, ValidationStat status);

    /// Return the number of parameter directions for the
    /// intersection.
    /// \return The number of parameter directions for the
    /// intersection.
    virtual int numParams() const
    { return 2; }
	
protected:
    // Data members

    virtual shared_ptr<Intersector> 
    lowerOrderIntersector(shared_ptr<ParamGeomInt> obj1,
			  shared_ptr<ParamGeomInt> obj2, 
			  Intersector* prev = 0,
			  int eliminated_parameter = -1,
			  double eliminated_value = 0);

    virtual int performRotatedBoxTest(double eps1, double eps2);

    virtual bool foundIntersectionNearBoundary();

    virtual int simpleCase2(Point& axis1, Point& axis2);

    virtual int checkCoincidence();

    virtual void microCase();
  
    virtual int updateIntersections();

    virtual int repairIntersections()
    { return 0; }

    virtual int linearCase();

    virtual int doSubdivide();

    virtual void writeOut();

private:
    void doIterate(double& par1, double& par2,
		   double& dist, double *guess=0);

    double distInCandidatePar(double par, int dir, const double* seed);

    int sortParameterDirections(int perm[]);

    SubdivisionClassification getSubdivisionParameter(int dir,
						      double& par);

};


} // namespace Go


#endif  // _CVCVINTERSECTOR_H
