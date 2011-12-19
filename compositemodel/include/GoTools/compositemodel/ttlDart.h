//===========================================================================
//                                                                           
// File: ttlDart.h
//                                                                           
// Created: March 1 2001
//                                                                           
// Author: Øyvind Hjelle <oyvind.hjelle@math.sintef.no>,
//         ............. <............ @math.sintef.no>
//                                                                           
// Revision: $Id: ttlDart.h,v 1.3 2004-01-21 08:29:55 bsp Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================
// Copyright (c) 2000 SINTEF Applied Mathematics
//===========================================================================
#ifndef _DART_
#define _DART_

#include "GoTools/compositemodel/ttlTriang.h"

namespace hetriang{
  //--------------------------------------------
  // Dart class for the half-edge data structure
  //--------------------------------------------
  /** \class Dart
  * \brief \b Dart class for the half-edge data structure.
  *
  * See \ref api for a detailed description of how the member functions
  * should be implemented.
  *
  *  This is how the dart class is implemented for the half-edge data structure:
  *  \include HeDart.h
  */
  class Dart {
    Edge* edge_;
    bool dir_; // true if dart is counterclockwise in face    
  public:
    /// Default constructur
    Dart() {edge_ = NULL; dir_ = true;}
    /// Constructur
    Dart(Edge* edge, bool dir = true) {edge_ = edge; dir_ = dir;}
    /// Copy constructor
    Dart(const Dart& dart) {edge_ = dart.edge_; dir_ = dart.dir_;}
    /// Destructor
    ~Dart() {}
		
    /// Assignment operator
    Dart& operator = (const Dart& dart) {
      if (this == &dart)
        return *this;
      edge_ = dart.edge_;
      dir_  = dart.dir_;
      return *this;
    }
    /// Comparing dart objects
    bool operator==(const Dart& dart) const {
      if (dart.edge_ == edge_ && dart.dir_ == dir_)
        return true;
      return false;
    }
    /// Comparing dart objects
    bool operator!=(const Dart& dart) const {
      return !(dart==*this);
    }
    /// Maps the dart to a different node
    Dart& alpha0(){dir_ = !dir_; return *this;}
    /// Maps the dart to a different edge
    Dart& alpha1() {
      if (dir_) {
        edge_ = edge_->getNextEdgeInFace()->getNextEdgeInFace();
        dir_ = false;
      }
      else {
        edge_ = edge_->getNextEdgeInFace();
        dir_ = true;
      }
      return *this;
    }
    /// Maps the dart to a different triangle. \b Note: the dart is not changed if it is at the boundary!
    Dart& alpha2() {
      if (edge_->getTwinEdge()) {
        edge_ = edge_->getTwinEdge();
        dir_ = !dir_;
      }
      // else, the dart is at the boundary and should not be changed
      return *this;
    }
    
    /** @name Utilities not required by TTL */
    //@{
    void init(Edge* edge, bool dir = true) {edge_ = edge; dir_ = dir;}
		
    double x() const {return getNode()->x();} // x-coordinate of source node
    double y() const {return getNode()->y();} // y-coordinate of source node
    bool isCounterClockWise() const {return dir_;}
		
    shared_ptr<Node> getNode() const {return dir_ ? edge_->getSourceNode() : edge_->getTargetNode();}
    shared_ptr<Node> getOppositeNode() const {return dir_ ? edge_->getTargetNode() : edge_->getSourceNode();}
    Edge* getEdge() const {return edge_;}

    /**
     * controllConstDart is a function who are suppost to check if the dart representing dstart or dend has been changed
     * during a swap. It's implemented so that this is compared to a spesified copy \e dart.
     * In Half-Edge datastrukture this is an ompossible situation, and the function is therefor empty
     */
    void controllConstDart ( Dart& dart ) {}

    Dart makeCopy () {
      Dart* dart = new Dart(edge_,dir_);
      return *dart; 
    }

    //debugfunction
    void printDart(){
      //not implemented
    }

    //@}
  };
}; // end half_edge_triang namespace

#endif
