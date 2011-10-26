/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 2000 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/


/********************************************************************
 FILENAME    : PrNestedNode.C
 AUTHOR      : Michael Floater, SINTEF
 DATE        : Nov. 2000
 DESCRIPTION : Implementation of methods in the class PrNestedNode.
 CHANGE LOG  :
*********************************************************************/

#include "GoTools/parametrization/PrNestedNode.h"
#include <iostream>

// MEMBER FUNCTIONS

//-----------------------------------------------------------------------------
void PrNestedNode::print(std::ostream& os)
//-----------------------------------------------------------------------------
{
    os << x() << ' ' << y() << ' ' << z() << ' ' << u_ << ' ' << v_
       << ' ' << level_ << std::endl;
    for(size_t i=0; i<tr_.size(); i++) os << tr_[i] << " ";
    os << std::endl;
}

//-----------------------------------------------------------------------------
void PrNestedNode::printXYZ(std::ostream& os)
//-----------------------------------------------------------------------------
{
    os << x() << ' ' << y() << ' ' << z() << std::endl;
}

//-----------------------------------------------------------------------------
void PrNestedNode::printUV(std::ostream& os)
//-----------------------------------------------------------------------------
{
    os << u_ << ' ' << v_ << std::endl;
}
