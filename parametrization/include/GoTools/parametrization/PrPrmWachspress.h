//===========================================================================
//                                                                           
// File: PrPrmWachspress.h                                                   
//                                                                           
// Created: Mon Jun 30 15:18:21 2003                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: PrPrmWachspress.h,v 1.5 2007-03-02 16:18:19 jbt Exp $
//                                                                           
//===========================================================================

#ifndef _PRPRMWACHSPRESS_H
#define _PRPRMWACHSPRESS_H

#include "GoTools/parametrization/PrParametrizeInt.h"


/** PrPrmWachspress -  Implement the Wachspress parametrization
 * by implementing the virtual function makeWeights  
 */
class PrPrmWachspress : public PrParametrizeInt
{
protected:
    virtual bool makeWeights(int i);
    double  tanThetaOverTwo(Vector3D& a, Vector3D& b, Vector3D& c);

public:
    /// Default constructor
    PrPrmWachspress();
    /// Empty destructor
    virtual ~PrPrmWachspress();

};


#endif // _PRPRMWACHSPRESS_H

