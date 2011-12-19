//===========================================================================
//                                                                           
// File: ftFaceBase.C                                                        
//                                                                           
// Created: Mon Jul  8 15:36:10 2002                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@math.sintef.no>
//                                                                           
// Revision: $Id: ftFaceBase.C,v 1.3 2008-12-01 14:02:55 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/compositemodel/ftFaceBase.h"


namespace Go
{

//===========================================================================
ftFaceBase::ftFaceBase()
//===========================================================================
{
    id_ = -1;
    //is_turned_ = false;
}

//===========================================================================
ftFaceBase::ftFaceBase(int id/*, bool is_turned*/)
//===========================================================================
{
    id_ = id;
    //is_turned_ = is_turned;
}

//===========================================================================
ftFaceBase::~ftFaceBase()
//===========================================================================
{
}


//===========================================================================
ftSurface* ftFaceBase::asFtSurface()
{
  return 0;
}
//===========================================================================


//===========================================================================
void ftFaceBase::setId(int id)
//===========================================================================
{
    id_ = id;
}

//===========================================================================
int ftFaceBase::getId()
//===========================================================================
{
    return id_;
}

// //---------------------------------------------------------------------------
// void ftFaceBase::turnFace(vector<ftFaceBase*>& turned)
// //---------------------------------------------------------------------------
// {
//   int ki;
//   for (ki=0; ki<(int)turned.size(); ki++)
//     if (turned[ki] == this)
//       return;   // This face is already turned.

//   //  std::cout << "Turning face: " << this << std::endl;
//   turnOrientation();  // Turn the current face
//   turned.push_back(this);

//   // Reverse the direction of all edges, and turn neighbouring surfaces.
//   vector<shared_ptr<ftEdgeBase> > start_edges = startEdges();
//   for (ki=0; ki<(int)start_edges.size(); ki++)
//     {
//       ftEdgeBase *e0 = start_edges[ki].get();
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

//---------------------------------------------------------------------------
void ftFaceBase::updateBoundaryLoops(shared_ptr<ftEdgeBase> new_edge)
//---------------------------------------------------------------------------
{
    // Nothing to do
}


} // namespace Go
