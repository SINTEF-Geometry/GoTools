                                                                 // -*- C++ -*-
/*****************************************************************************/
/*                                                                           */
/* (c) Copyright 1994, 95, 96                                                */
/*     SINTEF, Oslo, Norway.                                                 */
/*     All rights reserved.                                                  */
/*                                                                           */
/*     See the file Copyright.h for further details.                         */
/*                                                                           */
/*****************************************************************************/

//#include<Copyright.h>
// $Id: CPUclock.C,v 1.5 2005-06-09 07:29:51 oan Exp $

#include "GoTools/utils/CPUclock.h"

namespace Go {

CPUclock:: CPUclock ()
{
    last = getCurrentTime();
}

void CPUclock:: swap ()
{
    last = now;
}

void CPUclock:: subt ()
{
    diff = now - last;
}


double CPUclock:: getInterval ()
{
    now = getCurrentTime();
    subt();
    swap();
    return diff;
}

double CPUclock:: getTime ()
{
    now = getCurrentTime();
    subt();
    swap();
    return now;
}

}; // namespace Go

/* LOG HISTORY of this file:

$Log: 
*/
