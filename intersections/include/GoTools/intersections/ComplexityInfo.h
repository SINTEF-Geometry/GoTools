//==========================================================================
//                                                                          
// File: ComplexityInfo.h                                                    
//                                                                          
// Created: Wed Sep 27 14:44:37 2006                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: ComplexityInfo.h,v 1.1 2006-09-27 12:52:13 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================

#ifndef _COMPLEXITYINFO_H
#define _COMPLEXITYINFO_H


namespace Go {


/// This class contains statistical information used to check whether
/// an intersection problem is simplified during recursive subdivision

class ComplexityInfo {
public:
    /// Default constructor
    ComplexityInfo()
    {
	cone_overlap_ = min_box_overlap_ = -1.0;
	nmb_intpts_ = nmb_sing_intpts_ = 0;
	complex_ = false;
    }

    /// Destructor
    ~ComplexityInfo() { }

    /// Set the overlapping angle between two cones
    /// \param overlap the overlapping angle
    void setConeOverlap(double overlap)
    { cone_overlap_ = overlap; }

    /// Get the overlapping angle between two cones
    /// \return the overlapping angle
    double getConeOverlap()
    { return cone_overlap_; }

    /// Set the overlap between two boxes
    /// \param overlap the box overlap
    void setBoxOverlap(double overlap)
    { min_box_overlap_ = overlap; }

    /// Get the overlap between two boxes
    /// \return the box overlap
    double getBoxOverlap()
    { return min_box_overlap_; }

    /// Set the number of intersection points between two objects
    /// \param nmbpts the number of intersection points
    void setNmbIntpts(int nmbpts)
    { nmb_intpts_ = nmbpts; }

    /// Get the number of intersection points between two objects
    /// \return the number of intersection points
    int getNmbIntpts()
    { return nmb_intpts_; }

    /// Set the number of singular points
    /// \param nmbpts the number of singular points
    void setNmbSingpts(int nmbpts)
    { nmb_sing_intpts_ = nmbpts; }

    /// Get the number of singular points
    /// \return the number of singular points
    int getNmbSingpts()
    { return nmb_sing_intpts_; }

    /// Set flag for complex case
    void setComplex()
    { complex_ = true; }

    /// Get flag for complex case
    /// \return \c true if complex, \c false otherwise
    bool isComplex()
    { return complex_; }

private:
    double cone_overlap_; // A negative number implies not set or no
			  // overlap
    double min_box_overlap_;  // Negative number = not set or not
			      // overlap
    int nmb_intpts_;
    int nmb_sing_intpts_;
    bool complex_;

};


} // namespace Go


#endif // _COMPLEXITYINFO_H

