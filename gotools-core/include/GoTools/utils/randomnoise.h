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

