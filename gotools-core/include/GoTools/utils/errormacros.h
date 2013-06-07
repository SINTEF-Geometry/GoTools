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




