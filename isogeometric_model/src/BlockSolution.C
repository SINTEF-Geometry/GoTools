//===========================================================================
//
// File : BlockSolution.C
//
// Created: Mon Mar  1 08:54:51 2010
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================




// #include "GoTools/isogeometric_model/SfSolution.h"
#include "GoTools/isogeometric_model/BlockSolution.h"


namespace Go
{

  //===========================================================================
  BlockSolution::~BlockSolution()
  //===========================================================================
  {
  }


  //===========================================================================
  SfSolution* BlockSolution::asSfSolution()
  //===========================================================================
  {
    return NULL;
  }


  //===========================================================================
  VolSolution* BlockSolution::asVolSolution()
  //===========================================================================
  {
    return NULL;
  }


  //===========================================================================
  int BlockSolution::nmbCoefs(int pardir) const
  //===========================================================================
  {
    return basis(pardir).numCoefs();
  }


  //===========================================================================
  int BlockSolution::degree(int pardir) const
  //===========================================================================
  {
    return basis(pardir).order() - 1;
  }


  //===========================================================================
  vector<double> BlockSolution::knots(int pardir) const
  //===========================================================================
  {
    BsplineBasis bas = basis(pardir);
    vector<double> result(bas.order() + bas.numCoefs());
    vector<double>::const_iterator
      coef_begin = bas.begin(),
      coef_end = bas.end();
    copy(coef_begin, coef_end, result.begin());
    return result;
  }


  //===========================================================================
  vector<double> BlockSolution::distinctKnots(int pardir) const
  //===========================================================================
  {
    vector<double> result;
    basis(pardir).knotsSimple(result);
    return result;
  }


} // end namespace Go
