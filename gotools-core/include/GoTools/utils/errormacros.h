//===========================================================================
//                                                                           
// File: errormacros.h                                                       
//                                                                           
// Created: Tue May 14 12:22:16 2002                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: errormacros.h,v 1.7 2006-10-27 16:14:59 vsk Exp $
//                                                                           
// Description: Error macros. In order to use some of them, you must also
//              include <iostream> or <exception>. The compile defines NDEBUG
//              and NVERBOSE control the behaviour of these macros.
//                                                                           
//===========================================================================

#ifndef _ERRORMACROS_H
#define _ERRORMACROS_H

#include <iostream>

/// Usage: REPORT;
/// Usage: MESSAGE("Message string.");
/// Usage: THROW("Error message string.");
#ifdef NVERBOSE // Not verbose mode
#  ifndef REPORT
#    define REPORT
#  endif
#  ifndef MESSAGE
// Cannot be empty because of the comma operator in THROW(x)
#    define MESSAGE(x) 0
#  endif
#  ifndef MESSAGE_IF
#    define MESSAGE_IF(cond, m)
#  endif
#  ifndef THROW
#    define THROW(x) throw std::exception()
#  endif
#else // Verbose mode
#  ifndef REPORT
#    define REPORT std::cerr << "\nIn file " << __FILE__ << ", line " << __LINE__ << std::endl
#  endif
#  ifndef MESSAGE
#    define MESSAGE(x) std::cerr << "\nIn file " << __FILE__ << ", line " << __LINE__ << ": " << x << std::endl
#  endif
#  ifndef MESSAGE_IF
#    define MESSAGE_IF(cond, m) do {if(cond) MESSAGE(m);} while(0)
#  endif
#  ifndef THROW
#    define THROW(x) MESSAGE(x), throw std::exception()
#  endif
#endif

#ifndef GO_NO_CHECKS
#define GO_NO_CHECKS
#endif


#define ALWAYS_ERROR_IF(condition, message) do {if(condition){ THROW(message);}} while(0)

/// Usage: ASSERT(condition)
/// Usage: ASSERT2(condition, "Error message string.")
/// Usage: DEBUG_ERROR_IF(condition, "Error message string.");
#ifdef NDEBUG // Not in debug mode
#  ifndef ASSERT
#    define ASSERT(x)
#  endif
#  ifndef ASSERT2
#    define ASSERT2(cond, x)
#  endif
#  ifndef DEBUG_ERROR_IF
#    define DEBUG_ERROR_IF(cond, x)
#  endif
#else // Debug mode
#  ifndef ASSERT
#    define ASSERT(cond) if (!(cond)) THROW("Assertion \'" #cond "\' failed.")
#  endif
#  ifndef ASSERT2
#    define ASSERT2(cond, x) do { if (!(cond)) THROW(x);} while(0)
#  endif
#  ifndef DEBUG_ERROR_IF
//#    define DEBUG_ERROR_IF(cond, x) if (cond) THROW(x) 
#    define DEBUG_ERROR_IF(cond, x) do { if (cond) THROW(x); } while(0)
#  endif
#endif


#endif // _ERRORMACROS_H




