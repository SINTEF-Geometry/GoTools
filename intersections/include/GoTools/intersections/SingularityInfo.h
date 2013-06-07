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

#ifndef _SINGULARITYINFO_H
#define _SINGULARITYINFO_H


#include "GoTools/utils/Point.h"
#include <vector>
#include <memory>


namespace Go {


/// This class contains information about singularities in an
/// intersection problem

class SingularityInfo {
public:
    /// Default constructor
    SingularityInfo();

    /// Inheritance constructor
    /// \param previous shared pointer to a SingularityInfo object
    /// \param use_previous flag indicating if \a previous comes from
    /// a previous intersection problem
    SingularityInfo(shared_ptr<SingularityInfo> previous,
		    bool use_previous = false);

    /// Inheritance constructor
    /// \param previous shared pointer to a SingularityInfo object
    /// \param missing_dir if the dimension of the intersection
    /// problem has been reduced from previous, the index tells us
    /// which parameter direction that was removed. Indexing starts at
    /// 0.
    SingularityInfo(shared_ptr<SingularityInfo> previous,
		    int missing_dir);

    /// Check if singular point has been found
    /// \return \c true if a singular point has been found, \c false
    /// otherwise
    bool hasPoint()
    { return (param_sing_.size() > 0); }

    /// Get number of singular points
    /// \param nmb_par the number of parameters
    /// \return the number of sinular points
    int getNmbPoints(int nmb_par)
    { return (int)param_sing_.size()/nmb_par; }

    /// Get a parameter for th last found singularity or closest point
    /// \param dir the index of the parameter direction
    /// \return the parameter value in the \a dir direction
    double getParam(int dir);

    /// Increment the iteration count
    void addIterationCount()
    { 
	nmb_iter_++; 
	iteration_done_ = true;
    }

    /// Increment the number of times the intersection objects are
    /// characterized as simple.
    /// \param simple1 if \c true, increment for the first
    /// intersection object
    /// \param simple2 if \c true, increment for the second
    /// intersection object
    void addSimpleCount(bool simple1, bool simple2)
    {
	if (simple1)
	    nmb_simple1_++;
	if (simple2)
	    nmb_simple2_++;
    }

    /// Set a singular point
    /// \param par parameter array of the singular point
    /// \param nmb_par number of parameters in \a par
    void setSingularPoint(double* par, int nmb_par);

    /// Add a singular point
    /// \param par parameter array of the singular point
    /// \param nmb_par number of parameters in \a par
    void addSingularPoint(double* par, int nmb_par);

    /// Get a singular point
    /// \param idx index to the singular point in question
    /// \param nmb_par number of parameters
    /// \return the singular Point
    Point getPoint(int idx, int nmb_par);

    /// Check if the iteration is done
    /// \return \c true if the iteration is done, \c false otherwise
    bool iterationDone()
    { return iteration_done_; }

    /// Check if the iteration succeeded
    /// \return \c true if the iteration succeeded, \c false otherwise
    bool iterationSucceed()
    { return iteration_succeed_; }

    /// Get the number of times the first intersection object has been
    /// characterized as simple.
    /// \return the number of times
    int nmbSimple1()
    { return nmb_simple1_; }

    /// Get the number of times the second intersection object has
    /// been characterized as simple.
    /// \return the number of times
    int nmbSimple2()
    { return nmb_simple2_; }

    /// Set a high priority singular point
    /// \param par parameter array of the singular point
    /// \param nmb_par number of parameters in \a par
    void setHighPriSing(double* par, int nmb_par);

    /// Set the high priority singular point type
    /// \param type the type of the singular point
    void setHighPriSingType(SingularityClassification type)
    { high_pri_type_ = type; }
		  
    /// Get info on high priority singularities
    /// \return the singularity classification
    SingularityClassification hasHighPriSing()
    { return (high_pri_sing_.size() == 0) ? NO_SING : high_pri_type_; }

    /// Get the high priority singularity
    /// \param idx index specifying a parameter of the singular point
    /// \return the parameter of the singular point in the \a idx
    /// direction
    double getHighPriSing(int idx)
    { return high_pri_sing_[idx]; }

    /// Get the high priority singularity
    /// \return the vector of parameters of the singular point
    std::vector<double> getHighPriSing()
    { return high_pri_sing_; }

    // Remove singularities outside the current domain
    void cleanUp(std::vector<double> start,
		 std::vector<double> end,
		 double tol);

private:
    // Data members

    // Number of iterations for singularity or closest point performed
    // downward in the recursion
    int nmb_iter_;   

    // How often this iteration has succeeded
    int nmb_success_;

    // How often the first intersection object is characterized as
    // simple (only relevant in surface-surface intersection)
    int nmb_simple1_;

    // How often the second intersection object (if any) is
    // characterized as simple (only relevant in surface-surface
    // intersection)
    int nmb_simple2_;

    // Last found singularity or closest point (parameter values)
    std::vector<double> param_sing_;

    // If any iteration for singularity or closest point is performed
    // at this recursion level
    bool iteration_done_;

    // If the last iteration found a singularity
    bool iteration_succeed_;

    // Existence of high priority singularity
    SingularityClassification high_pri_type_;

    // Parameter of high priority singularity
    std::vector<double> high_pri_sing_;
};


} // namespace Go


#endif // _SINGULARITYINFO_H

