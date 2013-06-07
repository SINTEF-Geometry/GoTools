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

