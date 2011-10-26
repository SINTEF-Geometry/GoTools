//===========================================================================
//                                                                           
// File: AlgObjectInt.h
//                                                                           
// Created: Thu Jan 27 10:38:07 2005                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: AlgObjectInt.h,v 1.3 2006-02-22 14:52:04 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#ifndef _ALGOBJECTINT_H
#define _ALGOBJECTINT_H


namespace Go {


/// This class is a purely abstract base class, providing an interface
/// to the algebraic intersection objects.

class AlgObjectInt {

public:
    /// Constructor.
    AlgObjectInt(){};

    /// Destructor.
    virtual ~AlgObjectInt(){};

};


} // end namespace Go


#endif // _ALGOBJECTINT_H

