//===========================================================================
//                                                                           
// File: sleep.h                                                             
//                                                                           
// Created: Wed Mar 19 17:01:18 2003                                         
//                                                                           
// Author: Wu Yongwei
//                                                                           
// Revision: $Id: sleep.h,v 1.2 2004-01-14 11:19:19 bsp Exp $
//                                                                           
// Description: Defines cross-platform sleep, usleep, etc.
//                                                                           
//===========================================================================

/// Defines cross-platform sleep, usleep, etc.

#ifndef _SLEEP_H
#define _SLEEP_H
#ifdef _WIN32
# if defined(_NEED_SLEEP_ONLY) && (defined(_MSC_VER) || defined(__MINGW32__))
#  include <stdlib.h>
#  define sleep(t) _sleep((t) * 1000)
# else
#  include <windows.h>
#  define sleep(t)  Sleep((t) * 1000)
# endif
# ifndef _NEED_SLEEP_ONLY
#  define msleep(t) Sleep(t)
#  define usleep(t) Sleep((t) / 1000)
# endif
#else
# include <unistd.h>
# ifndef _NEED_SLEEP_ONLY
#  define msleep(t) usleep((t) * 1000)
# endif
#endif
#endif /* _SLEEP_H */
