//===========================================================================
//                                                                           
// File: PrPrmSymMeanValue.h                                                 
//                                                                           
// Created: Tue Jun 24 15:28:02 2003                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: PrPrmSymMeanValue.h,v 1.4 2007-03-02 16:18:19 jbt Exp $
//                                                                           
//===========================================================================

#ifndef _PRPRMSYMMEANVALUE_H
#define _PRPRMSYMMEANVALUE_H


#include "GoTools/parametrization/PrParametrizeInt.h"


/** PrPrmSymMeanValue -  Implements the virtual function makeWeights.
 * This is an attempt to improve the mean value method of
 * PrPrmMeanValue by making it symmetric.
 */
class PrPrmSymMeanValue : public PrParametrizeInt
{
protected:
    virtual bool makeWeights(int i);
    double  tanThetaOverTwo(Vector3D& a, Vector3D& b, Vector3D& c);

public:
    /// Default constructor
    PrPrmSymMeanValue();
    /// Empty destructor
    ~PrPrmSymMeanValue();

};



#endif // _PRPRMSYMMEANVALUE_H

