#ifndef _VALUES_H_
#define _VALUES_H_


#include <limits>


//#ifndef MICROSOFT
//#include <values.h>
//#endif

#ifndef MAXINT
// #define MAXINT 2147483647;
const int MAXINT = std::numeric_limits<int>::max();
#endif

#ifndef MAXDOUBLE
// #define MAXDOUBLE   1.79769313486231570e+308
const double MAXDOUBLE = std::numeric_limits<double>::max();
#endif

#ifndef M_PI
const double M_PI = 3.14159265358979323846;
#endif

#ifndef DEFAULT_SPACE_EPSILON
#define DEFAULT_SPACE_EPSILON 1e-10
#endif 

#ifndef DEFAULT_PARAMETER_EPSILON
//#define DEFAULT_PARAMETER_EPSILON 1e-12
#define DEFAULT_PARAMETER_EPSILON 1e-10
#endif


#endif

/** @file Values.h
 * Defines the constants  MAXINT, MAXDOUBLE and M_PI if they are not
 * defined by system.
 */
