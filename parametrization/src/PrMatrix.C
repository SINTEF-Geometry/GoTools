/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1998 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

/********************************************************************
 FILENAME    : PrMatrix.C
 AUTHOR      : Michael Floater, SINTEF
 DATE        : Dec 98
 DESCRIPTION : Implementation of methods in the class PrMatrix.
 CHANGE LOG  :
*********************************************************************/

#include "GoTools/parametrization/PrMatrix.h"
#include <iostream>

//-----------------------------------------------------------------------------
PrMatrix::~PrMatrix()
//-----------------------------------------------------------------------------
{
}

//-----------------------------------------------------------------------------
void PrMatrix::print(std::ostream& os)
//-----------------------------------------------------------------------------
{
  for(int i=0; i<rows(); i++)
  {
    for(int j=0; j<rows(); j++)
    {
      os << (*this)(i,j) << " ";
    }
    os << std::endl;
  }
}

