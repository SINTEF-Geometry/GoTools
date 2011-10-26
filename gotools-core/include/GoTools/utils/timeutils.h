//===========================================================================
//                                                                           
// File: timeutils.h                                                         
//                                                                           
// Created: Thu Dec  6 14:56:13 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: timeutils.h,v 1.4 2005-06-09 07:29:50 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _TIMEUTILS_H
#define _TIMEUTILS_H

namespace Go {

/// Number of seconds since some (probably system-dependent) epoch.
double getCurrentTime();

/// Sleep for sleep_time seconds
void systemSleep(double sleep_time);

}; // end namespace Go

#endif // _TIMEUTILS_H

