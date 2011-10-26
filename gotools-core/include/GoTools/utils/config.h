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

#endif // _CONFIG_H
