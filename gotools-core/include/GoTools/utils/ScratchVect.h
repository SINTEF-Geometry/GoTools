//===========================================================================
//                                                                           
// File: ScratchVect.h                                                           
//                                                                           
// Created: 
//                                                                           
// Author: Jon Mikkelsen <jami@sintef.no>
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _SCRATCHVECT_H
#define _SCRATCHVECT_H


#include <vector>
//#include "GoTools/utils/Utils.h"


namespace Go {

    /** A template Vector class that behaves much like the std::vector template,
     *  but stores its elements on the stack rather than on the heap as long as the
     *  total size stays less than 'N' (template argument).
     */

template <typename T, size_t N>
class ScratchVect {

private:
    T scratch_[N];
    std::vector<T> stl_vect_;
    T* begin_;
    T* end_;

public:
    typedef T* iterator;
    typedef const T*  const_iterator;

    /// Create a new, empty ScratchVect.
    ScratchVect() :
      stl_vect_(0),
      begin_(scratch_),
      end_(scratch_) {}

    /// Create a new ScratchVect with 'size' (uninitialized) elements
    ScratchVect(size_t size) :
      stl_vect_((size>N)?size:0) {
      if (size>N) {
	begin_ = &stl_vect_[0];
	end_ = (&stl_vect_[0])+size; 
      } else {
	begin_=scratch_;
	end_=(&scratch_[size]);
      }
    }

    /// Create a new ScratchVect with 'size' elements initialized to the
    /// value 'elem'.
    ScratchVect(size_t size, T elem) :
      stl_vect_((size>N)?size:0) {
      if (size>N) {
	begin_ = &stl_vect_[0];
	end_ = (&stl_vect_[0])+size; 
      } else {
	begin_=scratch_;
	end_=(&scratch_[size]);
      }
      std::fill(begin_, end_, elem);
    }



    /// Set the contents of this ScratchVect to be equal
    /// to that of the STL vector 'vec'.
    ScratchVect(const std::vector<T>& vec) {
	assign(vec.begin(), vec.end());
    }

    /// Set the contents of this ScratchVect to be equal to the contents
    /// of the range specified by 'beg' and 'end'.
    template<typename RAIter>
    ScratchVect(RAIter beg, RAIter end) {
	assign(beg, end);
    }

    /// Set the contents of this ScratchVector to be equal to that of 
    /// the ScratchVector 'orig'.  (A templatized copy constructor, if 
    /// you like).
    template <int M>
    ScratchVect(const ScratchVect<T,M> &orig) {
	assign(orig.begin(), orig.end());
    }

    /// Assignment operator
    ScratchVect& operator = (const ScratchVect &orig) {
	assign(orig.begin(), orig.end());
	return *this;
    }

    /// resize the ScratchVect
    void resize(size_t new_size) {
	const int cur_size = (const int)size();
	if ((int)new_size == cur_size) {
	    return;
	}
	if ((int)new_size < cur_size) {
	    if (cur_size < (int)N) {
		end_ -= (cur_size - new_size);		
	    } else if (new_size > (int)N) {
		// both cur_size and new_size > N
		end_ -= (cur_size - new_size);		
	    } else {
		// cur_size > N but new_size < N
		begin_ = scratch_;
		end_ = begin_ + new_size;
		stl_vect_.clear(); // release memory
	    }
	}
	// new_size > size
	if (new_size < N) {
	    end_ = begin_ + new_size;
	} else {
	    stl_vect_.resize(new_size);
	    begin_ = &stl_vect_[0];
	    end_ = begin_ + new_size;
	}
    }

    /// Set the contents of this ScratchVect to be equal to the
    /// contents of the range specified by 'beg' and 'end'.
    template<typename RAIter>
    void assign(RAIter beg, RAIter end) {
	size_t sz = end - beg;
	stl_vect_.clear();
	if (sz>N) {
	    stl_vect_.insert(stl_vect_.end(), beg, end);
	    begin_ = &stl_vect_[0];
	    end_ = (&stl_vect_[0])+sz; 
	} else {
	    begin_=scratch_;
	    end_=(&scratch_[sz]);
	    std::copy(beg, end, begin_);
	}
    }
     
    /// Get iterator to beginning of range
    inline iterator begin() { return begin_;}

    /// Get constant iterator to beginning of range
    inline const_iterator begin() const { return begin_;}

    /// Get iterator to one-past-end of range
    inline iterator end() { return end_;}

    /// Get constant iterator to one-past-end of range
    inline const_iterator end() const { return end_;}

    /// Get the number of elements in the ScratchVect
    inline size_t size() const {return (end_-begin_);}

    /// Element access operators (const and not-const)
    inline const T& operator [] (int i) const { return begin_[i]; }
    inline T& operator [] (int i)       { return begin_[i]; }
    
    /// Add an element to the back of the range.
    void push_back(const T &e)
    {
      size_t n=size();
      if (n<N)
	*(end_++)=e;
      else if (n>N)
      {
	stl_vect_.push_back(e);
	begin_=&stl_vect_[0];
	end_=&stl_vect_[0]+stl_vect_.size();
      }
      else
      {
	stl_vect_.reserve(2*n);
	stl_vect_.insert(stl_vect_.end(), begin_, end_);
	stl_vect_.insert(stl_vect_.end(), e);
	begin_=&stl_vect_[0];
	end_=&stl_vect_[0]+stl_vect_.size();
      }
    }
    
    /// Clear ScratchVect
    void clear()
    {
      if (begin_==scratch_)
	end_=scratch_;
      else
      {
	stl_vect_.clear();
	begin_=end_=scratch_;
      }
    }

};
}
#endif

