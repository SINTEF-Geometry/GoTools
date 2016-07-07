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

#ifndef _ELEMENTARYCURVE_H
#define _ELEMENTARYCURVE_H


#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/geometry/SplineCurve.h"


namespace Go {


/// \brief ElementaryCurve is a base class for elementary curves like
/// lines and circles. Such curves have natural parametrizations and
/// ElementaryCurve is therefore a subclass of ParamCurve.
class ElementaryCurve : public ParamCurve
{
public:
    ElementaryCurve();

    /// Virtual destructor, enables safe inheritance.
    virtual ~ElementaryCurve();

    // --- Functions inherited from GeomObject ---
    virtual ElementaryCurve* clone() const = 0;

    // --- Functions inherited from ParamCurve ---
    virtual void reverseParameterDirection(bool switchparam = false);

    // --- Functions specific to ElementaryCurve ---
    virtual SplineCurve* createSplineCurve() const = 0;

    /// Set bounds for the parametrization of the curve
    /// \param startpar start parameter
    /// \param endpar end parameter
    virtual void setParamBounds(double startpar, double endpar) = 0;

    // Translate the curve along a given vector
    virtual void translateCurve(const Point& dir) = 0;

    virtual bool isReversed() const;

    // Helper function for reverseParameterDirection() when switchparam
    // is true. Used when curve is a parameter curve and x and y coordinates
    // should be swapped.
    virtual void swapParameters2D() = 0;

    virtual ElementaryCurve* subCurve(double from_par, double to_par,
				      double fuzzy = DEFAULT_PARAMETER_EPSILON) const = 0;
protected:
    // Returns reversed parameter in [tmin, tmax] if isReversed_ is true
    void getReversedParameter(double& t) const;

    bool isReversed_;
};


} // namespace Go


#endif // _ELEMENTARYCURVE_H

