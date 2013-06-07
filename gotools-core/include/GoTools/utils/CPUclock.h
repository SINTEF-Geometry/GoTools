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

#ifndef CPUCLOCK_H
#define CPUCLOCK_H

#include <string>
#include "GoTools/utils/timeutils.h"

namespace Go {
/*<CPUclock:*/
///  A class for measuring CPU time in programs
class CPUclock
{
  double    last;  // previous clock value
  double    now;   // present clock value
  double    diff;  // diff = now - last

  void   swap();// last = now
  void   subt();// subtract: now - last, store in diff

public:
  /// The constructor initializes the time of last call.
  CPUclock ();

  /// The increase in time (seconds) since last call to  getInterval()
  /// or to  initTime ().
  double getInterval  ();  // returns user time
    //  double getInterval2 ();  // returns user time + system time

  /// Start measuring the length of an interval (typically called
  /// prior to the computational job to be measured).
  void   initTime ()     { getInterval(); }

    //  string report (const string& message);

  /// Returns time in seconds since January 1, 1970
  double getTime ();
};

};// namespace Go
/*>CPUclock:*/

/*Class:CPUclock

NAME:  CPUclock - measures the CPU time in C++ programs

SYNTAX: @CPUclock


KEYWORDS:

  CPU time, clock


DESCRIPTION:

  One can get the user time in absolute seconds, or the length of intervals
  can be measured (in seconds).

  "getTime" - returns time in seconds since January 1, 1970

  "initTime" - start measuring the length of an interval (typically called
               prior to the computational job to be measured). In fact,
               "initTime" is just a call to "getInterval" - the function
               was for making programs easier to read.

  "getInterval" - returns the length of the user time interval (in seconds)
                between the present call to getInterval and the last call
                to initTime or to getInterval

  "getInterval2" - as "getInterval", but the user time plus the system time
                   is returned.

  "report" - returns a string ("String") containing a message (f.ex. where in
             the program the "report" function was called) and the user and
             system time since last call to "initTime" or "getInterval"
             (or "getInterval2").

SEE ALSO:

  "time.h", "sys/times.h"

DEVELOPED BY:

                SINTEF Applied Mathematics, Oslo, Norway, and
                University of Oslo, Dept. of Mathematics, Norway

AUTHOR:

  	        Hans Petter Langtangen, SINTEF/UiO

End:
*/




#endif



/* LOG HISTORY of this file:

$Log: not supported by cvs2svn $
Revision 1.6  2004/04/02 06:06:38  bsp
Bug fix in getTime()

Revision 1.5  2004/01/14 11:19:19  bsp
Added doxygen code

Revision 1.4  2002/08/15 15:05:26  sbr
WinNT compilation fix.

Revision 1.3  2002/04/26 09:30:59  ers
Moved CPUclock from parametrization to utils

Revision 1.1  2000/10/20 10:28:38  afr
First commit.

Revision 1.1  2000/07/06 07:17:27  afr
Added a few new files and changed GoHandle a bit.

 * Revision 1.14  1996/11/15  10:58:31  job
 * Version 2.4.0
 *
*/
