//===========================================================================
//                                                                           
// File: ttlTriang.C
//                                                                           
// Created: March 1 2001
//                                                                           
// Author: Øyvind Hjelle <oyvind.hjelle@math.sintef.no>,
//         ............. <............ @math.sintef.no>
//                                                                           
// Revision: $Id: ttlTriang.C,v 1.8 2008-07-10 17:57:17 kfp Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================
// Copyright (c) 2000 SINTEF Applied Mathematics
//===========================================================================
#include "GoTools/compositemodel/ttlTraits.h"
#include "GoTools/compositemodel/ttlTriang.h"
#include "ttl/ttl.h"
using namespace hetriang;



#include <iostream> 
#ifndef TTL_USE_OLD_STD
#include <algorithm>
#include <fstream>
#else
#include <algorithm.h>
#include <fstream.h>
#endif
using namespace Go;
using std::vector;
using std::list;
using std::min;
using std::max;
using std::ofstream;
using std::endl;


Triangulation* TTLtraits::triang_ = NULL;
//int Node::id_count = 0;

#ifdef DEBUG_HE
  static void errorAndExit(char* message) {
  cout << "\n!!! ERROR: "<< message << " !!!\n" << endl; exit(-1);}
#endif


//---------------------------------------------------------
static Edge* getLeadingEdgeInTriangle(Edge* e) {
  
  Edge* edge = e;
  
  // Code: 3EF (assumes triangle)
  if (!edge->isLeadingEdge()) {
    edge = edge->getNextEdgeInFace();
    if (!edge->isLeadingEdge())
      edge = edge->getNextEdgeInFace();
  }
  
  if (!edge->isLeadingEdge()) {
    return NULL;
  }
  
  return edge;
}

//------------------------------------------------------------------------------
static void getLimits(vector<ttlPoint*>::iterator first,
                      vector<ttlPoint*>::iterator last,
                      double& xmin, double& ymin,
                      double& xmax, double& ymax) {
  
  xmin = ymin = 1.0e+30;
  xmax = ymax = -1.0e+30;
  
  vector<ttlPoint*>::iterator it;
  for (it = first; it != last; ++it) {
    xmin = min(xmin, (*it)->x());
    ymin = min(ymin, (*it)->y());
    xmax = max(xmax, (*it)->x());
    ymax = max(ymax, (*it)->y());
  }
}

Edge*
Triangulation::initTwoEnclosingTriangles(vector<ttlPoint*>::iterator first,
                                         vector<ttlPoint*>::iterator last) {
  
  double xmin, ymin, xmax, ymax;
  getLimits(first, last, xmin, ymin, xmax, ymax);
  
  // Add 10% of range:  
  double fac = 10.0;
  double dx  = (xmax-xmin)/fac;
  double dy  = (ymax-ymin)/fac;
  
  double zval = 0;
  shared_ptr<Node> n1(new Node(xmin-dx,ymin-dy,zval));
  shared_ptr<Node> n2(new Node(xmax+dx,ymin-dy,zval));
  shared_ptr<Node> n3(new Node(xmax+dx,ymax+dy,zval));
  shared_ptr<Node> n4(new Node(xmin-dx,ymax+dy,zval));
  
  // diagonal
  Edge* e1d = new Edge; // lower
  Edge* e2d = new Edge; // upper, the twin edge
  
  // lower triangle
  Edge* e11 = new Edge;
  Edge* e12 = new Edge;
  
  // upper triangle
  Edge* e21 = new Edge; // upper upper
  Edge* e22 = new Edge;
  
  // lower triangle
  e1d->setSourceNode(n3);
  e1d->setNextEdgeInFace(e11);
  e1d->setTwinEdge(e2d);
  e1d->setAsLeadingEdge();
  addLeadingEdge(e1d);
  
  e11->setSourceNode(n1);
  e11->setNextEdgeInFace(e12);
  
  e12->setSourceNode(n2);
  e12->setNextEdgeInFace(e1d);
  
  // upper triangle
  e2d->setSourceNode(n1);
  e2d->setNextEdgeInFace(e21);
  e2d->setTwinEdge(e1d);
  e2d->setAsLeadingEdge();
  addLeadingEdge(e2d);
  
  e21->setSourceNode(n3);
  e21->setNextEdgeInFace(e22);
  
  e22->setSourceNode(n4);
  e22->setNextEdgeInFace(e2d);
  
  return e11;
  //dc1.init(e11);
  //dc2.init(e12);
  //dc3.init(e21);
  //dc4.init(e22);
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void Triangulation::createDelaunay(vector<ttlPoint*>::iterator first,
                                   vector<ttlPoint*>::iterator last) {
  
  TTLtraits::triang_ = this;
  cleanAll();
    
  Edge* bedge = initTwoEnclosingTriangles(first, last);
  Dart dc(bedge);
  
  Dart d_iter = dc;
  
  vector<ttlPoint*>::iterator it;
  for (it = first; it != last; ++it) {
      ttl::insertNode<TTLtraits>(d_iter, **it);
  }
  
  // In the general case the initial dart may have been changed.
  // It is the users responsibility to get a valid boundary dart here.
  // (A dart at the boundary can also be found by trying to locate a
  // triangle "outside" the triangulation.)

  // Assumes rectangular domain
  ttl::removeRectangularBoundary<TTLtraits>(dc);
}


//------------------------------------------------------------------------------
void Triangulation::removeTriangle(Edge& edge) {
  
  Edge* e1 = getLeadingEdgeInTriangle(&edge);

#ifdef DEBUG_HE
  if (e1 == NULL)
    errorAndExit("Triangulation::removeTriangle: could not find leading edge");
#endif
  
  removeLeadingEdgeFromList(e1);
  // cout << "No leading edges = " << leadingEdges_.size() << endl;  
  // Remove the triangle
  Edge* e2 = e1->getNextEdgeInFace();
  Edge* e3 = e2->getNextEdgeInFace();
  
  if (e1->getTwinEdge())
    e1->getTwinEdge()->setTwinEdge(NULL);
  if (e2->getTwinEdge())
    e2->getTwinEdge()->setTwinEdge(NULL);
  if (e3->getTwinEdge())
    e3->getTwinEdge()->setTwinEdge(NULL);
  
  delete e1;
  delete e2;
  delete e3;
}


/* OK, but not needed
//------------------------------------------------------------------------------
void Triangulation::removeBoundaryNode(Dart& dart) {
  
  // Assumes that swapping has been done in the Delaunay case
  // No dart given as output
  // Assumes that the dart represents a boundary edge
  // Note, convex boundary is not necessarily preserved.
  
  // Removing boundary node but NOT swap Delaunay
  
  Dart d_iter = dart;
  Dart dnext = dart;
  bool bend = false;
  while (bend == false) {
    dnext.alpha1().alpha2();
    if (ttl::isBoundaryEdge(dnext))
      bend = true; // Stop when boundary
    
    removeTriangle(*d_iter.getEdge());
    d_iter = dnext;
  }
}
*/

//------------------------------------------------------------------------------
void Triangulation::reverse_splitTriangle(Edge& edge) {
  
  // Reverse operation of splitTriangle
  
  Edge* e1 = edge.getNextEdgeInFace();
  Edge* le = getLeadingEdgeInTriangle(e1);
#ifdef DEBUG_HE
  if (le == NULL)
    HeUtil::errorAndExit("Triangulation::removeTriangle: could not find leading edge");
#endif
  removeLeadingEdgeFromList(le);
  
  Edge* e2 = e1->getNextEdgeInFace()->getTwinEdge()->getNextEdgeInFace();
  le = getLeadingEdgeInTriangle(e2);
#ifdef DEBUG_HE
  if (le == NULL)
    HeUtil::errorAndExit("Triangulation::removeTriangle: could not find leading edge");
#endif
  removeLeadingEdgeFromList(le);
  
  
  Edge* e3= edge.getTwinEdge()->getNextEdgeInFace()->getNextEdgeInFace();
  le = getLeadingEdgeInTriangle(e3);
#ifdef DEBUG_HE
  if (le == NULL)
    HeUtil::errorAndExit("Triangulation::removeTriangle: could not find leading edge");
#endif
  removeLeadingEdgeFromList(le);
  
  // The three triangles at the node have now been removed
  // from the triangulation, but the arcs have not been deleted.
  // Next delete the 6 half edges radiating from the node
  // The node is maintained by handle and need not be deleted explicitly
  
  Edge* estar = &edge;
  Edge* enext = estar->getTwinEdge()->getNextEdgeInFace();
  delete estar->getTwinEdge();
  delete estar;  
  
  estar = enext;
  enext = estar->getTwinEdge()->getNextEdgeInFace();
  delete estar->getTwinEdge();
  delete estar;
  
  delete enext->getTwinEdge();
  delete enext;
  
  // Create the new triangle
  e1->setNextEdgeInFace(e2);
  e2->setNextEdgeInFace(e3);
  e3->setNextEdgeInFace(e1);
  addLeadingEdge(e1);
}


//------------------------------------------------------------------------------
// "Template" for iterating the boundary    
/*
static void iterateBoundary2(const Dart& dart) {
cout << "Iterate boundary 2" << endl;
// input is a dart at the boundary

  Dart dart_iter = dart;
  do {
  if (ttl::isBoundaryEdge(dart_iter))
  dart_iter.alpha0().alpha1();
  else
  dart_iter.alpha2().alpha1();
  
    } while(dart_iter != dart);
    }
*/


// -----------------------------------
Dart Triangulation::createDart() { 
  
  // Return an arbitrary CCW dart
  return Dart(*leadingEdges_.begin());
}

//------------------------------------------------------------------------------
bool Triangulation::removeLeadingEdgeFromList(Edge* leadingEdge) {
  
  // Remove the edge from the list of leading edges,
  // but don't delete it.
  // Also set flag for leading edge to false.
  // Must search from start of list. Since edges are added to the
  // start of the list during triangulation, this operation will
  // normally be fast (when used in the triangulation algorithm)
  list<Edge*>::iterator it;
  for (it = leadingEdges_.begin(); it != leadingEdges_.end(); ++it) {
    
    Edge* edge = *it;
    if (edge == leadingEdge) {
      
      edge->setAsLeadingEdge(false);
      leadingEdges_.erase(it);
      
      break;
    }
  }
  
  if (it == leadingEdges_.end())
    return false;
  
  return true;
}

//--------------------------------------------------------------------------
/*
 * Picks out the edge with the node equal the spesified point
 * Function useless????
 * Just construct the Edges instead???
 */
/*Edge* Triangulation::getConstrainedEdges(ttlPoint point) const{
  cout << "getConstrainedEdges 1" << endl;//debugline
	const list<Edge*>& leadingEdges = getLeadingEdges();//picks out the list and makes it constant
	list<Edge*>::const_iterator it;
	cout << leadingEdges.size() << endl;//debugline
	cout << "getConstrainedEdges 2" << endl;//debugline
  for (it = leadingEdges.begin(); it != leadingEdges.end(); ++it) {
    cout << "getConstrainedEdges" << endl;//debugline
    Edge* edge = *it;
    if (edge->getSourceNode()->x() == point.x() && 
			edge->getSourceNode()->y() == point.y() && 
			edge->getSourceNode()->z() == point.z()) {
			cout << "fant en node" << endl;//debugline
      return edge;
		}//endif
  }//endfor
  return NULL;//
};*/



//------------------------------------------------------------------------------
void Triangulation::cleanAll() {
  
  list<Edge*>::const_iterator it;  
  for (it = leadingEdges_.begin(); it != leadingEdges_.end(); ++it) {
    Edge* e1 = *it;
    Edge* e2 = e1->getNextEdgeInFace();
    Edge* e3 = e2->getNextEdgeInFace();
    
    delete e1;
    delete e2;
    delete e3;
  }
  
  leadingEdges_.clear();
}


// for all nodes template, but multiple
/*
void 
Triangulation::flagNodes(bool flag) const {

  list<Edge*>::const_iterator it;
  for (it = leadingEdges_.begin(); it != leadingEdges_.end(); ++it) {
  Edge* edge = *it;
  
    for (int i = 0; i < 3; ++i) {
    edge->getSourceNode()->flag_ = flag;
    edge = edge->getNextEdgeInFace();
    }
    }
    }
    
      list<Edge*>*
      Triangulation::getNodes() const {
      
        flagNodes(false);
        list<Edge*>* nodeList = new list<Edge*>;
        
          list<Edge*>::const_iterator it;
          for (it = leadingEdges_.begin(); it != leadingEdges_.end(); ++it) {
          Edge* edge = *it;
          
            for (int i = 0; i < 3; ++i) {
            shared_ptr<Node> node = edge->getSourceNode();
            
              if (node->flag_ == false) {
              nodeList->push_back(edge);
              node->flag_ = true;
              }
              edge = edge->getNextEdgeInFace();
              }
              }
              return nodeList;
              }
*/

//------------------------------------------------------------------------------
list<Edge*>* Triangulation::getEdges(bool skip_boundary_edges) const {
  
  // collect all arcs (one half edge for each arc)
  // (boundary edges are also collected).
  
  list<Edge*>::const_iterator it;
  list<Edge*>* elist = new list<Edge*>;
  for (it = leadingEdges_.begin(); it != leadingEdges_.end(); ++it) {
    Edge* edge = *it;
    for (int i = 0; i < 3; ++i) {
      Edge* twinedge = edge->getTwinEdge();
      // only one of the half-edges
      
      if ( (twinedge == NULL && !skip_boundary_edges) ||
           (twinedge != NULL && ((long int)edge > (long int)twinedge)) )
        elist->push_front(edge);
      
      edge = edge->getNextEdgeInFace();
    }
  }
  return elist;
}

//------------------------------------------------------------------------------
Edge* Triangulation::splitTriangle(Edge& edge, const ttlPoint& point3d) {
  
  // Add a node by just splitting a triangle into three triangles
  // Assumes the half edge is located in the triangle
  // Returns a half edge with source node as the new node
  
  double x = point3d.x();
  double y = point3d.y();
  double z = point3d.z();
  PointIter pnt_iter = point3d.pnt_iter();
  
  // e#_n are new edges
  // e# are existing edges
  // e#_n and e##_n are new twin edges
  // e##_n are edges incident to the new node
  
  // Add the node to the structure
  shared_ptr<Node> node(new Node(pnt_iter,x,y,z));
  shared_ptr<Node> n1 = edge.getSourceNode();
  Edge* e1 = &edge;
  
  Edge* e2 = edge.getNextEdgeInFace();
  shared_ptr<Node> n2 = e2->getSourceNode();
  
  Edge* e3 = e2->getNextEdgeInFace();
  shared_ptr<Node> n3 = e3->getSourceNode();
  
  Edge* e1_n  = new Edge;
  Edge* e11_n = new Edge;
  Edge* e2_n  = new Edge;
  Edge* e22_n = new Edge;
  Edge* e3_n  = new Edge;
  Edge* e33_n = new Edge;
  
  e1_n->setSourceNode(n1);
  e11_n->setSourceNode(node);
  e2_n->setSourceNode(n2);
  e22_n->setSourceNode(node);
  e3_n->setSourceNode(n3);
  e33_n->setSourceNode(node);
  
  e1_n->setTwinEdge(e11_n);
  e11_n->setTwinEdge(e1_n);
  e2_n->setTwinEdge(e22_n);
  e22_n->setTwinEdge(e2_n);
  e3_n->setTwinEdge(e33_n);
  e33_n->setTwinEdge(e3_n);
  
  e1_n->setNextEdgeInFace(e33_n);
  e2_n->setNextEdgeInFace(e11_n);
  e3_n->setNextEdgeInFace(e22_n);
  
  e11_n->setNextEdgeInFace(e1);
  e22_n->setNextEdgeInFace(e2);
  e33_n->setNextEdgeInFace(e3);
  
  
  // and update old's next edge
  e1->setNextEdgeInFace(e2_n);
  e2->setNextEdgeInFace(e3_n);
  e3->setNextEdgeInFace(e1_n);
  
  // add the three new leading edges, 
  // Must remove the old leading edge from the list.
  // Use the field telling if an edge is a leading edge
  // NOTE: Must search in the list!!!
  
  
  Edge* leadingEdge;
  if (e1->isLeadingEdge())
    leadingEdge = e1;
  else if (e2->isLeadingEdge())
    leadingEdge = e2;
  else if(e3->isLeadingEdge())
    leadingEdge = e3;
  else
    return NULL;
  
  removeLeadingEdgeFromList(leadingEdge);
  
  
  addLeadingEdge(e1_n);
  addLeadingEdge(e2_n);
  addLeadingEdge(e3_n);
  
  // Return a half edge incident to the new node (with the new node as source node)
  
  return e11_n;
}


//------------------------------------------------------------------------------
void Triangulation::swapEdge(Edge& diagonal) {
  
  // Note that diagonal is both input and output and it is always
  // kept in counterclockwise direction (this is not required by all 
  // finctions in ttl:: now)
  
  // Swap by rotating counterclockwise
  // Use the same objects - no deletion or new objects
  Edge* eL   = &diagonal;
  Edge* eR   = eL->getTwinEdge();
  Edge* eL_1 = eL->getNextEdgeInFace();
  Edge* eL_2 = eL_1->getNextEdgeInFace();
  Edge* eR_1 = eR->getNextEdgeInFace();
  Edge* eR_2 = eR_1->getNextEdgeInFace();
  
  // avoid node to be dereferenced to zero and deleted
  shared_ptr<Node> nR = eR_2->getSourceNode();
  shared_ptr<Node> nL = eL_2->getSourceNode();
  
  eL->setSourceNode(nR);
  eR->setSourceNode(nL);
  
  // and now 6 1-sewings
  eL->setNextEdgeInFace(eL_2);
  eL_2->setNextEdgeInFace(eR_1);
  eR_1->setNextEdgeInFace(eL);
  
  eR->setNextEdgeInFace(eR_2);
  eR_2->setNextEdgeInFace(eL_1);
  eL_1->setNextEdgeInFace(eR);
  
  Edge* leL = 0;
  if (eL->isLeadingEdge())
    leL = eL;
  else if (eL_1->isLeadingEdge())
    leL = eL_1;
  else if (eL_2->isLeadingEdge())
    leL = eL_2;
  
  Edge* leR = 0;
  if (eR->isLeadingEdge())
    leR = eR;
  else if (eR_1->isLeadingEdge())
    leR = eR_1;
  else if (eR_2->isLeadingEdge())
    leR = eR_2;
  
  removeLeadingEdgeFromList(leL);
  removeLeadingEdgeFromList(leR);
  addLeadingEdge(eL);
  addLeadingEdge(eR);
}

// -----------------------------------------------------
//static void printEdge(const Dart& dart, ostream& ofile) {
//  
//  Dart d0 = dart;
//  d0.alpha0();
//  
//  ofile << dart.x() << " " << dart.y() << endl;
//  ofile << d0.x() << " " << d0.y() << endl;
//}

//------------------------------------------
bool Triangulation::checkDelaunay() const {
  
  // ???? outputs !!!!
  // ofstream os("qweND.dat");
  const list<Edge*>& leadingEdges = getLeadingEdges();
  
  list<Edge*>::const_iterator it;
  bool ok = true;
  int noNotDelaunay = 0;
  
  for (it = leadingEdges.begin(); it != leadingEdges.end(); ++it) {
    Edge* edge = *it;
    
    for (int i = 0; i < 3; ++i) {
      Edge* twinedge = edge->getTwinEdge();
      
      // only one of the half-edges
      if (twinedge == NULL || (long int)edge > (long int)twinedge) {
        Dart dart(edge);
        if (ttl::swapTestDelaunay<TTLtraits, Dart>(dart, false)) {
          noNotDelaunay++;
          
          //printEdge(dart,os); os << "\n";
          ok = false;
          //cout << "............. not Delaunay .... " << endl;
        }
      }
      edge = edge->getNextEdgeInFace();
    }
  }
  
  
  //if (ok)
    //cout << "---> Triangulation is OK acc. to Delaunay swap test\n" << endl;
  //else
    //cout << "!!! Triangulation is NOT Delaunay: " << noNotDelaunay << " edges\n" << endl;
  
  return ok;
}

// -----------------------------------------------------
void Triangulation::optimizeDelaunay() {

  // This function is also present in ttl where it is implemented
  // generically.
  // The implementation below is tailored for the half-edge data structure,
  // and is thus more efficient

  // Collect all interior edges (one half edge for each arc)
  bool skip_boundary_edges = true;
  list<Edge*>* elist = getEdges(skip_boundary_edges);
  
  // Assumes that elist has only one half-edge for each arc.
  bool cycling_check = true;
  bool optimal = false;
  list<Edge*>::const_iterator it;
  while(!optimal) {
    optimal = true;
    for (it = elist->begin(); it != elist->end(); ++it) {
      Edge* edge = *it;
      
      Dart dart(edge);
      if (ttl::swapTestDelaunay<TTLtraits>(dart, cycling_check)) {
        optimal = false;
        swapEdge(*edge);
      }
    }
  }
  delete elist;
}

// -----------------------------------------------------
Edge* Triangulation::getInteriorNode() const {
  
  const list<Edge*>& leadingEdges = getLeadingEdges();
  list<Edge*>::const_iterator it;
  for (it = leadingEdges.begin(); it != leadingEdges.end(); ++it) {
    Edge* edge = *it;
    
    // multiple checks, but only until found
    for (int i = 0; i < 3; ++i) {
      if (edge->getTwinEdge() != NULL) {
        
        if (!ttl::isBoundaryNode(Dart(edge)))
          return edge;
      }
      edge = edge->getNextEdgeInFace();
    }
  }
  return NULL; // no boundary nodes
}


// -----------------------------------------------------
static Edge* getBoundaryEdgeInTriangle(Edge* edge) {
  
  if (ttl::isBoundaryEdge(Dart(edge)))
    return edge;  
  edge = edge->getNextEdgeInFace();
  if (ttl::isBoundaryEdge(Dart(edge)))
    return edge;
  edge = edge->getNextEdgeInFace();
  if (ttl::isBoundaryEdge(Dart(edge)))
    return edge;
  
  return NULL;
}

Edge* Triangulation::getBoundaryEdge() const {

  // Get an arbitrary (CCW) boundary edge
  // If the triangulation is closed, NULL is returned

  const list<Edge*>& leadingEdges = getLeadingEdges();
  list<Edge*>::const_iterator it;
  Edge* edge;
  
  for (it = leadingEdges.begin(); it != leadingEdges.end(); ++it) {
    edge = getBoundaryEdgeInTriangle(*it);
    
    if (edge)
      return edge;
  }
  return NULL;
};




// ----------------------------------------------------------
void Triangulation::printEdges(ofstream& os) const {
      
  // Print source node and target node for each edge face by face,
  // but only one of the half-edges.
  
  const list<Edge*>& leadingEdges = getLeadingEdges();
	//cout << leadingEdges.size() << endl;//debugline
  list<Edge*>::const_iterator it;
  for (it = leadingEdges.begin(); it != leadingEdges.end(); ++it) {
    Edge* edge = *it;
    //cout << "printEdges" << endl;//debugline
    for (int i = 0; i < 3; ++i) {
      Edge* twinedge = edge->getTwinEdge();
      
      // Print only one edge (the highest value of the pointer)
      if (twinedge == NULL || (long int)edge > (long int)twinedge) {
        // Print source node and target node
        shared_ptr<Node> node = edge->getSourceNode();
        os << node->x() << " " << node->y() << " " << node->z() << endl;
        node = edge->getTargetNode();
        os << node->x() << " " << node->y() << " " << node->z() << endl;
        os << '\n'; // blank line
      }
      edge = edge->getNextEdgeInFace();
    }
  }
}
