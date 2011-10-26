//===========================================================================
//                                                                           
// File: GeoTol.h
//                                                                           
// Created: 
//                                                                           
// Author: 
//                                                                           
// Revision: $Id: GeoTol.h,v 1.10 2006-02-23 11:42:07 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

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
