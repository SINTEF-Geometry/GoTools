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
  for (auto i = v.begin(); i != v.end(); ++i) object_to_stream(os, *i);
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
