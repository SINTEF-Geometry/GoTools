//===========================================================================
//
// File : BdConditionType.h
//
// Created: Thu Feb 25 13:55:37 2010
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================


#ifndef __BDCONDITIONTYPE_H
#define __BDCONDITIONTYPE_H



namespace Go
{

  enum BdConditionType
  {
    UNKNOWN = 0,
    ZERO_DIRICHLET,
    CONSTANT_DIRICHLET,
    DIRICHLET,
    ZERO_NEUMANN,
    NEUMANN,
    SYMMETRY
  };

}  // End namespace Go


#endif    // #ifndef __BDCONDITIONTYPE_H
