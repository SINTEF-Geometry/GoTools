//===========================================================================
//                                                                           
// File: CurveTesselator.C                                                 
//                                                                           
// Created: Thu Nov 29 12:52:31 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: CurveTesselator.C,v 1.1 2008-03-27 16:14:05 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/tesselator/CurveTesselator.h"

namespace Go
{


//===========================================================================
CurveTesselator::~CurveTesselator()
//===========================================================================
{
}

//===========================================================================
void CurveTesselator::tesselate()
//===========================================================================
{
    Point pt(3);
    int n = mesh_->numVertices();
    for (int i = 0; i < n; ++i) {
	double rt = double(i)/double(n-1);
	double ta = curve_.startparam();
	double tb = curve_.endparam();

	// Ensure a finite extension in case of an unbounded curve
// 	double limit = 1.0e3;
// 	ta = std::max(ta, -1.0*limit);
// 	tb = std::min(tb, 1.0*limit);

	double t = ta*(1.0 - rt) + rt*tb;
	curve_.point(pt, t);
	mesh_->vertexArray()[i*3] = pt[0];
	mesh_->vertexArray()[i*3 + 1] = pt[1];
	mesh_->vertexArray()[i*3 + 2] = pt[2];
    }
}

} // namespace Go

