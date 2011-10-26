//===========================================================================
//
// File : IsogeometricBlock.C
//
// Created: Fri Feb 26 15:11:49 2010
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================




#include "GoTools/isogeometric_model/IsogeometricBlock.h"
#include <algorithm>


namespace Go
{

  //===========================================================================
  IsogeometricBlock::~IsogeometricBlock()
  //===========================================================================
  {
  }


  //===========================================================================
  IsogeometricSfBlock* IsogeometricBlock::asIsogeometricSfBlock()
  //===========================================================================
  {
    return NULL;
  }


  //===========================================================================
  IsogeometricVolBlock* IsogeometricBlock::asIsogeometricVolBlock()
  //===========================================================================
  {
    return NULL;
  }


  //===========================================================================
  tpTolerances IsogeometricBlock::getTolerances() const
  //===========================================================================
  {
    return model_->getTolerances();
  }


  //===========================================================================
  int IsogeometricBlock::nmbCoefs(int pardir) const
  //===========================================================================
  {
    return basis(pardir).numCoefs();
  }


  //===========================================================================
  int IsogeometricBlock::degree(int pardir) const
  //===========================================================================
  {
    return basis(pardir).order() - 1;
  }


  //===========================================================================
  vector<double> IsogeometricBlock::knots(int pardir) const
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
  vector<double> IsogeometricBlock::distinctKnots(int pardir) const
  //===========================================================================
  {
    vector<double> result;
    basis(pardir).knotsSimple(result);
    return result;
  }


} // end namespace Go
