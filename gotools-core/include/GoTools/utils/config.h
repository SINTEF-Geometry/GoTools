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
#else
#include <memory>
using std::shared_ptr;
using std::dynamic_pointer_cast;
using std::const_pointer_cast;
using std::static_pointer_cast;
#endif

#endif // _CONFIG_H
