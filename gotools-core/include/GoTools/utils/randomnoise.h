//===========================================================================
//                                                                           
// File: randomnoise.h                                                       
//                                                                           
// Created: Thu Feb 12 09:55:04 2004                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: randomnoise.h,v 1.2 2005-06-06 09:32:01 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _RANDOMNOISE_H
#define _RANDOMNOISE_H

namespace Go {

//===========================================================================
//                    FUNCTIONS FOR RANDOM DATA
//===========================================================================
//! Gives a certain number of random samples drawn from the normal distribution.
//! \param res a pointer to the array where the resulting samples should be written.
//! \param mean_err the sigma parameter to the normal distribution.
//! \param num_samples the desired number of samples (should also be the size of the
//! array pointed to by \em res.
void normalNoise(double* res, double mean_err, int num_samples);

//! Gives a certain number of random samples drawn from the uniform distribution.
//! \param res a pointer to the array where the resulting samples should be written.
//! \param lval lower bound of the range from which the samples can take their values.
//! \param uval upper bound of the range from which the samples can take their values.
//! \param num_samples the desired number of samples (should also be the size of the
//! array pointed to by \em res.
void uniformNoise(double* res, double lval, double uval, int num_samples);

}; // namespace Go


#endif // _RANDOMNOISE_H

