//===========================================================================
//                                                                           
// File: ftMessage.h                                                   
//                                                                           
// Created: Mon Mar 19 2001                                         
//                                                                           
// Author: Vibeke Skytt 
//                                                                           
// Revision: 
//                                                                           
// Description: 
//
//                                                                           
//===========================================================================

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

       
