//===========================================================================
//                                                                           
// File: GpuMatrix.hpp                                                
//                                                                           
// Created: 2009-01-10                                     
//                                                                           
// Authors: Andr\'e Rigland Brodtkorb
//                                                                           
// Revision: $Id: GpuMatrix.hpp,v 1.1 2009-02-10 11:54:25 babrodtk Exp $
//                                                                           
// Description: A matrix represented in GPU memory
//                                                                           
//===========================================================================

#ifndef GPUMATRIX_H_
#define GPUMATRIX_H_

#include <iostream>


namespace Go {
    
    template <typename T>
    class GpuMatrix {
        public:
        /**
         * Constructs a matrix on the GPU
         */
        GpuMatrix(int rows, int cols);
        
        /**
         * Constructs a matrix on the GPU using data as initial values
         */
        GpuMatrix(int rows, int cols, T* data);
        
        /**
         * Constructs a matrix on the GPU setting all values to initializer.
         */
        GpuMatrix(int rows, int cols, T initializer);
        
        /**
         * Deletes data allocated on the gpu.
         */
        ~GpuMatrix();
      
        /**
         * Uploads data from cpu memory to gpu memory. Using "pinned memory"
         * will give best performance.
         */
        void upload(T* data);
        
        /**
         * Downloads data from gpu memory to cpu memory. Using "pinned memory"
         * will give best performance.
         */
        void download(T* data);
        
        int getRows();
        int getCols();
        
        /**
         * Prints out the contents of the matrix m.
         */
        template <typename U> 
        friend std::ostream& operator<< (std::ostream& os, const GpuMatrix<U>& m);
        
        private:
            T* h_ptr;
            T* d_ptr;
            int rows;
            int cols;
            int pitch;
                   
    };
    
}

#include "GoTools/implicitization/GpuMatrix_impl.hpp"

#endif /*GPUMATRIX_H_*/
