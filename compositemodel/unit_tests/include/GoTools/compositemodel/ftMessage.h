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

#ifndef _FTMESSAGE_H
#define _FTMESSAGE_H

//===========================================================================
//===========================================================================

#include <vector>

namespace Go
{

enum ftmessages
{
  FT_NOT_IMPLEMENTED = -99,
  FT_NO_DATA,
  FT_BAD_INPUT_FILE,
  FT_ERROR_IN_READ_IGES,
  FT_DISCONTINUITY,
  FT_UNEXPECTED_DATA_TYPE,
  FT_NON_4_SIDED_SURF,
  FT_WRONG_NO_OF_BOUNDARIES,
  FT_INSUFFICIENT_MESH_CURVES,
  FT_ERROR_IN_TOPOLOGY_ANALYSIS,
  FT_TOPOLOGY_PROBLEM,
  FT_ERROR_IN_SURFACE_CREATION,
  FT_ERROR_IN_PARAMETERIZE,
  FT_ERROR_IN_SURFACE_TRIMMING,
  FT_NOT_SUPPORTED,
  FT_NOT_SPLINE_SURF,
  FT_ERROR_IN_SURFACE_GRIDDING,
  FT_NEIGHBOUR_GRID_SIZE_DIFFER,
  FT_OK = 0,
  FT_APPROX_ERROR_TOO_LARGE,
  FT_COULD_NOT_REFINE_SURFACE,
  FT_FEATURE_ERROR_TOO_LARGE,
  FT_DEGENERATE_PATCH_CREATED,
  FT_SURFACE_ALREADY_CREATED,
  FT_DEGENERATE_SURFACE,
  FT_UNEXPECTED_INPUT_OBJECT_IGNORED,
  FT_FACES_OVERLAP,
  FT_COULD_NOT_SMOOTH_SURFACE,
  FT_COULD_NOT_SMOOTH_SURFACESET,
  FT_NO_SMOOTHING,
  FT_FAILED_CREATING_GRAPH,
  FT_GRID_NOT_CREATED,
  FT_SURFACE_NOT_MODIFIED
};

/** ftMessage -  Error and warning messages. Used only in some contexts
 * 
 */
  
class ftMessage
{
  public:
  
  /// Constructor
       ftMessage() 
	 : message_(FT_OK)
	 {}

  /// Constructor
       ftMessage(ftmessages message) 
	 : message_(message)
	 {}

  /// Destructor
       ~ftMessage(){}
       
       /// Set message
       void setError(ftmessages error)
	 { message_ = error; }

       /// Add a warning
       void addWarning(ftmessages warning)
	 { warnings_.push_back(warning); }

  // Inquiry
       /// Check for error and warning messages
       bool isOK()
	 { return message_ == FT_OK; }

       /// Return current message
       ftmessages getMessage()
	 { return message_; }

       /// Check number of warnings
       int noOfWarnings()
    { return (int)warnings_.size(); }

       /// Return requested warning
       ftmessages getWarning(int i)
	 { 
	   if (i < (int) warnings_.size())
	     return warnings_[i]; 
	   else
	     return FT_OK;
	 }


 private:
  ftmessages message_;
  std::vector<ftmessages> warnings_;
};

} // namespace Go



#endif //  _FTMESSAGE_H

       
