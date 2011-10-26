//===========================================================================
//                                                                           
// File: GeoTol.C 
//                                                                           
// Created: 
//                                                                           
// Author: 
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/intersections/GeoTol.h"
#include <iostream>

using namespace std;
namespace Go
{
//===========================================================================
GeoTol::GeoTol(double epsge, double rel_par_res, double numerical_tol)
    : epsge_(epsge), 
      eps_bracket_(0.1), 
      ang_tol_(0.01), 
      rel_par_res_(rel_par_res),
      numerical_tol_(numerical_tol)
//===========================================================================
  {
      if (epsge <= 10e-5)
	  ref_ang_ = 0.01;
      else if (epsge <= 10e-3)
	  ref_ang_ = 0.1;
      else
	  ref_ang_ = 0.5;
  }

//===========================================================================
GeoTol::GeoTol(GeoTol *epsge)
{
    epsge_ = epsge->epsge_;
    eps_bracket_ = epsge->eps_bracket_;
    ref_ang_ = epsge->ref_ang_;
    ang_tol_ = epsge->ang_tol_;
    rel_par_res_ = epsge->rel_par_res_;
    numerical_tol_ = epsge->numerical_tol_;
}

//===========================================================================
  GeoTol::~GeoTol()
//===========================================================================
{
}

//===========================================================================
void GeoTol::write(ostream& os) const
//===========================================================================
{
    os << epsge_ << ' ';
    os << eps_bracket_ << ' ';
    os << ref_ang_ << ' ';
    os << ang_tol_ << ' ';

    os << rel_par_res_ << ' ';
    os << numerical_tol_ << ' ';
}

//===========================================================================
void GeoTol::read(istream& is)
//===========================================================================
{
    is >> epsge_;
    is >> eps_bracket_;
    is >> ref_ang_;
    is >> ang_tol_;
    is >> rel_par_res_;
    is >> numerical_tol_;
}

};
