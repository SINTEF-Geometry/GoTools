/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1997 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#ifndef PRHEAP_H
#define PRHEAP_H


/*<PrHeap-syntax: */

#include <iostream>

/** PrHeap - This is a class for minimum heaps. The PrHeap's functionality
 * is to store prioritized elements. It has a function for adding
 * elements, a function for returning (and removing from the heap) the
 * element with the "best" priority (least key) and a function for modifying
 * elements. These functions are implemented in a way making sure that
 * they are very fast in run time. All functions run in O(log n) time.
 */
class PrHeap
{
private:
  int     maxsize_; 
  int     size_;            
  double*   keys_;            // priority values
  int*   elements_;        // eg. indexes 
  int*   pos_of_elements_; // Back-pointers
 

public:
  /// Default constructor
  PrHeap  ( );

  /// Constructor. Creates a heap wich is able to store maxsize elements.
  /// maxsize will typically be set to the total number of elements 
  /// (here total number of vertices).  
  PrHeap  ( int maxsize);

  /// Constructor. Not implemented.
  PrHeap  ( int* elms, int n);

  /// Copy constructor.Not implemented, generates a warning.      
  PrHeap  ( const PrHeap& heap);

  /// Destructor. Deletes arrays.
  ~PrHeap ( ); 
  
  void redim( int size);

  /// return max. size for the heap.
  int  getMaxSize() const { return maxsize_;}

  /// return size for the heap.   
  int  getSize()    const { return size_;}

  /// Inserts an element with a key in the heap.
  void  push   ( double key,     int element);

  /// Modifies an elements key and restores the heap property. 
  void  modify ( double new_key, int element);

  /// Extracts the best element from the heap.
  int  pop( );
  
  /// Extracts the best element from the heap.
  void  pop( double& key, int& element );

  /// Restarts the heap.
  void  emptyHeap();

  /// print to stream
  void print( std::ostream& os) const;

private:
  void upHeap   ( int k);
  void downHeap ( int k);
};
/*>PrHeap-syntax: */
 

/*Class:PrHeap

Name:              "PrHeap" -

Syntax:	           @PrHeap-syntax

Keywords:          PrHeap

Description: 
  This is a class for minimum heaps. The PrHeap's functionality
  is to store prioritized elements. It has a function for adding
  elements, a function for returning (and removing from the heap)
  the element with the "best" priority (least key) and a function for modifying
  elements. These functions are implemented in a way making sure that
  they are very fast in run time. All functions run in O(log n)
  time.

See also:

Member functions:



  PUBLIC MEMBER FUNCTIONS:\\

  "Constructor" --\\
  Creates a heap wich is able to store maxsize elements.  maxsize will 
  typically be set to the total number of elements (here total number 
  of vertices). 

  "Destructor" --\\
  Deletes arrays.

  "Copy constructor" --\\
  Not implemented, generates a warning.

  "push" --\\
  Inserts an element with a key in the heap.

  "pop" --\\
  Extracts the best element from the heap.

  "modify" --\\
  Modifies an elements key and restores the heap property. 

  "emptyHeap" --\\
  Restarts the heap.

  "getMaxSize" --\\
  "getCurrentSize" --\\
  Sizes for the heap.

PRIVATE MEMBER FUNCTIONS: \\

  "upHeap"  --\\ Heap maintainance function.

  "downHeap" --\\ Heap maintainance function.

  PrHeap(int)   heap(size);

Example:


Developed by:      SINTEF Applied Mathematics, Oslo, Norway
Author:	           Martin Reimers, SINTEF
                   Modified to avoid BasicTools by Michael Floater, Oct 2000
Date:              1997
*/

#endif // PRHEAP_H
