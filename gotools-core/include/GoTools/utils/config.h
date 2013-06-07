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

#ifndef _CONFIG_H
#define _CONFIG_H

#ifdef __BORLANDC__
# define GO_API __declspec(package)
#elif defined(MICROSOFT) || defined(_MSC_VER)
# if defined(__DLL__) || defined(_DLL)
#  define GO_API __declspec(dllexport)
# else
#  define GO_API __declspec(dllimport)
# endif // __DLL__
#else
# define GO_API
#endif // __BORLANDC__

// The following pragma is not optimal, but it's a workaround to
// getting rid of warning C4251 in Visual Studio
#ifdef _MSC_VER
#pragma warning( disable: 4251 )
#endif // _MSC_VER

#ifdef USE_BOOST
#include <boost/shared_ptr.hpp>
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using boost::const_pointer_cast;
using boost::static_pointer_cast;
#include <boost/static_assert.hpp>
#define static_assert(x, msg) BOOST_STATIC_ASSERT(x)
#include <boost/type_traits.hpp>
using boost::is_floating_point;
#else
#include <memory>
using std::shared_ptr;
using std::dynamic_pointer_cast;
using std::const_pointer_cast;
using std::static_pointer_cast;
#include <type_traits>
using std::is_floating_point;
#endif

#endif // _CONFIG_H
