//===========================================================================
//                                                                           
// File: tpFace.C                                                           
//                                                                           
// Created: Thu Jul 11 13:14:35 2002                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@math.sintef.no>
//                                                                           
// Revision: $Id: tpFace.C,v 1.7 2002-07-12 08:02:53 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/topology/tpFace.h"
#include "GoTools/topology/tpEdge.h"

using namespace Go;

// //---------------------------------------------------------------------------
// void tpFace::turnFace(vector<tpFace*>& turned)
// //---------------------------------------------------------------------------
// {
//   int ki;
//   for (ki=0; ki<turned.size(); ki++)
//     if (turned[ki] == this)
//       return;   // This face is already turned.

//   //  std::cout << "Turning face: " << this << std::endl;
//   turnOrientation();  // Turn the current face
//   turned.push_back(this);

//   // Reverse the direction of all edges, and turn neighbouring surfaces.
//   vector<shared_ptr<tpEdge> > start_edges = startEdges();
//   for (ki=0; ki<start_edges.size(); ki++)
//     {
//       tpEdge *e0 = start_edges[ki].get();
//       bool finished = false;
//       while (!finished)
// 	{
// 	  e0->turnOrientation();
// 	  if (e0->twin())
// 	    e0->twin()->face()->turnFace(turned);
// 	  e0 = e0->next();
// 	  if (e0 == start_edges[ki].get())
// 	    finished = true;
// 	}
//     }
// }
