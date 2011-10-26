//===========================================================================
//                                                                           
// File: Point.C                                                           
//                                                                           
// Created: Tue Jun  6 15:41:32 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: Point.C,v 1.5 2005-07-14 11:30:27 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/utils/Point.h"
#include <iostream> // @@ remove

using namespace Go;

//===========================================================================
void Point::read(std::istream& is)
//===========================================================================
{
    ALWAYS_ERROR_IF (n_ == 0,
		     "Trying to read into an 0-dimensional (empty) point.");
    for (int i = 0; i < n_; ++i)
	is >> pstart_[i];
}


//===========================================================================
void Point::write(std::ostream& os) const
//===========================================================================
{
    std::streamsize prev = os.precision(16);
    for (int i = 0; i < n_-1; ++i)
	os << pstart_[i] << ' ';
    os << pstart_[n_-1];
    os.precision(prev);   // Reset precision to it's previous value
}

