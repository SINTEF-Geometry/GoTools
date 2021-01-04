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

#ifndef STREAM_UTILS_H_
#define STREAM_UTILS_H_

#include <iostream>
#include <vector>

namespace { // anonymous, local namespace
  const char separator = ' ';
};

// =============================================================================
// Generic templates for sending objects to and retreaving objects from stream
// =============================================================================

// Generic write
template <typename T> void object_to_stream(std::ostream&  os, const T& obj) { os << obj << separator;}
template <typename T> void object_to_stream(std::wostream& os, const T& obj) { os << obj << separator;}

// Generic read
template <typename T> void object_from_stream(std::istream&  is, T& obj) { is >> obj; }
template <typename T> void object_from_stream(std::wistream& is, T& obj) { is >> obj; }

// =============================================================================
// SPECIALIZED TEMPLATES FOR CONTAINERS/OTHER PARTICULAR OBJECTS
// =============================================================================

// =============================================================================
// Write specialization for STL vectors
template<typename T>
void object_to_stream(std::ostream& os, const std::vector<T>& v)
// =============================================================================
{ 
  os << v.size() << separator;
  for (auto i = v.begin(); i != v.end(); ++i) 
object_to_stream(os, *i);
  os << '\n';
}

// =============================================================================
// Read specialization for STL vectors
template<typename T>
void object_from_stream(std::istream& is, std::vector<T>& v)
// =============================================================================
{
  size_t size;
  is >> size;
  v.resize(size);
  for (auto i = v.begin(); i != v.end(); ++i) { object_from_stream(is, *i);}
}

// =============================================================================
// Write specialization for STL vectors
template<typename T>
void object_to_stream(std::wostream& os, const std::vector<T>& v)
// =============================================================================
{ 
  os << v.size() << separator;
  for (auto i = v.begin(); i != v.end(); ++i ) object_to_stream(os, *i);
  os << '\n';
}

// =============================================================================
// Read specialization for STL vectors
template<typename T>
void object_from_stream(std::wistream& is, std::vector<T>& v)
// =============================================================================
{
  size_t size;
  is >> size;
  v.resize(size);
  for (auto i = v.begin(); i != v.end(); ++i) { object_from_stream(is, *i);}
}



#endif
