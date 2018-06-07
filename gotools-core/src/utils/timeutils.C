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
#if(__GNUC__ >= 4 && __GNUC_MINOR__ >= 3) 
	    usleep((useconds_t)(sleep_time*1e6));
#else
	    usleep(sleep_time*1e6);
#endif
#endif
	}
    }

} // end namespace Go
