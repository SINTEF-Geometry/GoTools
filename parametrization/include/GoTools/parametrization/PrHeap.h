/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

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
