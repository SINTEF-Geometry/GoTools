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

#ifndef _GEOTOL_H
#define _GEOTOL_H


#include <fstream>


namespace Go {

// @@@ VSK 0904, As long as this class only contains a constant tolerance,
// then it is OK just to pass it down into Intersectors with fewer parameter
// directions. When/if the tolerance gets parameter dependent, then this
// is not good enough. Then also the tolerance class needs a reduction
// in the number of parameter directions. Thus, the constructors in the
// Intersector tree is not final.

// @@@jbt Some functions are left un-Doxygen-documented, since this
// class is partly under development


/// Class handling various tolerances

class GeoTol {
public:
    /// Constructor
    GeoTol(double epsge, double rel_par_res = 1.0e-12, 
	   double numerical_tol = 1.0e-16);

    /// Constructor
    GeoTol(GeoTol *epsge);

    /// Destructor
    ~GeoTol();

    /// Get constant tolerance
    const double& getEpsge() const 
    { return epsge_; }

    const double& getMinEpsge() const 
    { return epsge_; }  // For the time being until a varying
			// tolerance is implemented
 
    /// Get relative parameter resolution
    const double& getRelParRes() const 
    { return rel_par_res_;}

    const double& getMaxEpsge() const 
    { return epsge_; }   // For the time being

    /// Get "bracket" tolerance used to specify confidence intervals
    /// in the parameter domain
    const double& getEpsBracket() const 
    { return eps_bracket_;}

    const double& getRefAng() const 
    { return ref_ang_; }

    /// Get angular tolerance
    const double& getAngleTol() const 
    { return ang_tol_; }
    
    /// Get numerical tolerance
    const double& getNumericalTol() const
    { return numerical_tol_;  }

    /// Write to output stream
    void write(std::ostream& os) const;
    /// Read from input stream
    void read(std::istream& is);
    
private:
    double epsge_;  // Constant tolerance
    double eps_bracket_;
    double ref_ang_;
    double ang_tol_;

    // It should also be room for a varying tolerance for use together
    // with functions and a directional dependent tolerance for use
    // with subdivision and boundary curves and points.

    double rel_par_res_;
    double numerical_tol_;
};


} // namespace Go


#endif // _GEOTOL_H
