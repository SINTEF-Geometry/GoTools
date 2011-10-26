/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1997 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

/********************************************************************
 FILENAME    : PrHeap.C
 AUTHOR      : Martin Reimers, SINTEF
 DATE        : 1997
 DESCRIPTION : Implementation of methods in the class PrHeap.
 CHANGE LOG  : Made independent of BasicTools by Michael FLoater, Oct. 2000
*********************************************************************/

#include "GoTools/parametrization/PrHeap.h"
#include <iostream>


// -------------------------------------------------------------------------
PrHeap::PrHeap  ()
// -------------------------------------------------------------------------
{  
  maxsize_  = 0;
  size_     = 0;
}

// -------------------------------------------------------------------------
PrHeap::PrHeap (int maxsize)
// -------------------------------------------------------------------------
{
  size_            = 0;
  maxsize_         = maxsize;
  keys_            = new double[maxsize+1];
  elements_        = new int[maxsize+1];

  pos_of_elements_ = new int[maxsize+1];
  for (int i = 1 ; i<= maxsize; i++)
      pos_of_elements_[i] = -1;
}
 

// -------------------------------------------------------------------------
PrHeap::PrHeap (const PrHeap& )
// -------------------------------------------------------------------------
{
  std::cout << "PrHeap::PrHeap(const PrHeap& heap)" <<
    "copy constructor called" << std::endl;
}


// -------------------------------------------------------------------------
PrHeap::~PrHeap ()
// -------------------------------------------------------------------------
{
  delete[] keys_;
  delete[] elements_;  
  delete[] pos_of_elements_;
}


// -------------------------------------------------------------------------
void PrHeap::redim( int maxsize)
// -------------------------------------------------------------------------
{
  // todo: hva med redim av pos_of_elements_ ??
  // not using [0] => size+1
  double*   keys_temp     = new double[maxsize+1];
  int*   elements_temp = new int[maxsize+1];

  int*   pos_of_elements_temp = new int[maxsize+1];

  for ( int j=1 ; j <= size_; j++) {
    keys_temp[j]     = keys_[j];
    elements_temp[j] = elements_[j];

    pos_of_elements_temp[j] = pos_of_elements_[j];
  }

  delete[] keys_;
  delete[] elements_;  
  delete[] pos_of_elements_;

  keys_            = keys_temp;
  elements_        = elements_temp;

  maxsize_ = maxsize;
  pos_of_elements_ = pos_of_elements_temp;
  for( int i = 1; i<=maxsize_; i++)
    pos_of_elements_[i] = -1;

}


// -------------------------------------------------------------------------
void PrHeap::upHeap( int k)
// -------------------------------------------------------------------------
{  
  if ( ( k > size_) || (k < 1))
   std::cout << "PrHeap::upHeap, index" << k <<
           "out of range, size = " << size_ << std::endl;

  double  target_key     = keys_[k];
  int  target_element = elements_[k];

  while (k > 1 && keys_[k/2] > target_key) // parent > child
    {
      keys_[k]     = keys_[k/2];     // k/2 down to k
      elements_[k] = elements_[k/2]; // k/2 down to k 

      pos_of_elements_[elements_[k/2]] = k; // update position

      k = k/2;
    }
  keys_[k]     = target_key;
  elements_[k] = target_element;

  pos_of_elements_[target_element] = k;
}

// -------------------------------------------------------------------------
void PrHeap::downHeap( int k)
// -------------------------------------------------------------------------
{  
  if ( ( (k > size_) && (size_ != 0) ) || (k < 1))
   std::cout << "PrHeap::downHeap, index" << k <<
           "out of range, size = " << size_ << std::endl;

  if ( k+k > size_) return; // no children, finished

  int j;
  double  target_key     = keys_[k];
  int  target_element = elements_[k];

  while ( k <= size_/2 ) { // at least one child
    j = k+k;   // j index to left child

    if ( (j < size_) && // l & r child
	 ( keys_[j] > keys_[j+1])) 
      {	j++;} // right child smallest
              // else left child smallest
    if ( target_key <= keys_[j]) break; // finished
    keys_[k]     = keys_[j];            // 2k up to k
    elements_[k] = elements_[j];        // 2k up to k

    pos_of_elements_[elements_[j]] = k; // update position

    k = j;
  }
  keys_[k]     = target_key;
  elements_[k] = target_element;

  pos_of_elements_[target_element] = k; // update position

}

// ---------------------------------------------------------------------------
void PrHeap::modify( double  key, int element)
// ---------------------------------------------------------------------------
{
  if ( key < keys_[pos_of_elements_[element]]) 
  {
    keys_[pos_of_elements_[element]] = key;
    upHeap( pos_of_elements_[element]);
  }
  else if ( key > keys_[pos_of_elements_[element]]) 
  {
    keys_[pos_of_elements_[element]] = key;
    downHeap( pos_of_elements_[element]);
  }
}


// ---------------------------------------------------------------------------
void PrHeap::push( double key, int element) 
// ---------------------------------------------------------------------------
{
  if( size_ >= maxsize_) // todo when stable?: redim( 2*maxsize_);
    std::cout << "PrHeap::insert, Heap full" << std::endl;

  if ( pos_of_elements_[element] != -1) // != -1:  in heap
    modify(key, element);
  else {                                // == -1 : not in heap
    // push element
    size_++;
    keys_[ size_]             = key;
    elements_[size_]          = element;
    pos_of_elements_[element] = size_; // update position
    upHeap (size_);
  }

}


// ---------------------------------------------------------------------------
int PrHeap::pop()
// ---------------------------------------------------------------------------
{
  if(size_ < 1) return -1;     // empty

  int     top = elements_[1];
  keys_[1]     = keys_[size_];
  elements_[1] = elements_[size_];

  pos_of_elements_[top] = -1;  // mark as not in heap
  pos_of_elements_[elements_[size_]] = 1; // update position

  size_--;
  downHeap(1);
  return top;
}

// ---------------------------------------------------------------------------
void PrHeap::pop( double& key, int& element )
// ---------------------------------------------------------------------------
{
  if(size_ < 1) {
    std::cout << "PrHeapIndirect::pop : Nothing to pop\n" << std::endl;
    element = -1;
    key     = -1.; 
    return;
  }

  key     = keys_[1];
  element = pop();
}


// ---------------------------------------------------------------------------
void PrHeap::emptyHeap()
// ---------------------------------------------------------------------------
{
  size_ = 0;
  for (int i = 1; i<= maxsize_; i++)
    pos_of_elements_[i] = -1;
}


// ---------------------------------------------------------------------------
void PrHeap::print( std::ostream& os) const
// ---------------------------------------------------------------------------
{
  int i;
  os << "\nPrint from class PrHeap:\n"
     << "-------------\n"
     << "Max size " << maxsize_  << std::endl
     << "Size     " << size_     << std::endl;

  os << "Keys: ";
  for ( i = 1; i <= size_ ; i++)
    os << keys_[i] << "\t";

  os << std::endl;

  os << "Elements: ";
  for ( i = 1; i <= size_ ; i++)
    os << elements_[i] << "\t";

  os << "Position of elements in heap: ";
  for ( i = 1; i <= size_ ; i++)
    os << pos_of_elements_[i] << "\t";

  os << "\n-------------\n";
}

