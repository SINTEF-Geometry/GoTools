//===========================================================================
//                                                                           
// File: timeutils.C                                                         
//                                                                           
// Created: Thu Dec  6 14:56:38 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: timeutils.C,v 1.5 2005-06-09 07:29:53 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/utils/timeutils.h"
#include "GoTools/utils/sleep.h"



// Time includes (crossplatform implementation)
#ifdef WIN32
#include <sys/timeb.h>
#include <stdlib.h>
#endif

#ifdef _WIN32_WCE
#include <afx.h>
#include <stdlib.h>
#endif

#ifdef __GNUC__
#include <sys/time.h>
#endif

#ifndef _WIN32_WCE
#include <time.h>
#endif

#ifndef WIN32
#include <unistd.h>
#endif

namespace Go {

// Time code (crossplatform implementation)
#ifdef WIN32
#ifndef _WIN32_WCE
#ifdef __BORLANDC__
#  define _timeb timeb
#  define _ftime ftime
#  define _sleep sleep
#endif

    struct _timeb  timeb_ptr_;       // Struct for holding time info, MS.
// Converts (Microsoft specific?) tstruct to seconds and milliseconds.
//-----------------------------------------------------------------------------
    inline double timeb2seconds(struct _timeb* timeb_ptr)
//-----------------------------------------------------------------------------
    {
	return (timeb_ptr->time + 0.001*timeb_ptr->millitm);
    }
#endif
#endif

#ifdef __GNUC__
// Removed...
#endif

#ifdef SGI
    inline double timespec2seconds( struct timespec* time_spec );
    // Converts from timespec to seconds
//-----------------------------------------------------------------------------
    double timespec2seconds(struct timespec* time_spec)
//-----------------------------------------------------------------------------
    {
	double value = (int)(time_spec->tv_sec)
	    + ((double)(time_spec->tv_nsec)/1e9);
	return value;
    }
#endif

#ifdef HPUX
    inline double timespec2seconds( struct timespec* time_spec );
    // Converts from timespec to seconds
//-----------------------------------------------------------------------------
    double timespec2seconds(struct timespec* time_spec)
//-----------------------------------------------------------------------------
    {
	double value = (int)(time_spec->tv_sec)
	    + ((double)(time_spec->tv_nsec)/1e9);
	return value;
    }
#endif


//-----------------------------------------------------------------------------
    double getCurrentTime()
//-----------------------------------------------------------------------------
    {
#ifdef WIN32
#ifndef _WIN32_WCE
	_ftime(&timeb_ptr_);
	return timeb2seconds(&timeb_ptr_);
#else
	// Might overrun every month...
	SYSTEMTIME s;
	GetSystemTime(&s);
	return 0.001*s.wMilliseconds + s.wSecond + 60.0*s.wMinute
	    + 3600.0*s.wHour + 3600.0*24.0*s.wDay;
#endif
#else
#ifdef __GNUC__
	timeval t;
	gettimeofday(&t,0);
	return (double)t.tv_sec + 1e-6*(double)t.tv_usec;
#else
	timespec time_spec_;
	clock_gettime(CLOCK_REALTIME, &time_spec_);
	return timespec2seconds(&time_spec_);
#endif
#endif
    }


//-----------------------------------------------------------------------------
    void systemSleep(double sleep_time)
//-----------------------------------------------------------------------------
    {
	// Sleep if sleep time is positive
	if( sleep_time > 0 )
	{
#ifdef WIN32
#ifndef _WIN32_WCE
	    Sleep((unsigned long)(sleep_time*1000));
#else
	    // Do not sleep...
#endif
#else
#ifdef __GNUC__
	    usleep((__useconds_t)(sleep_time*1e6));
#else
	    usleep(sleep_time*1e6);
#endif
#endif
	}
    }

} // end namespace Go
