//==========================================================================
//                                                                          
// File: SubdivisionClassification.h                                         
//                                                                          
// Created: Wed Sep 27 14:37:20 2006                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: SubdivisionClassification.h,v 1.1 2006-09-27 12:52:13 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================

#ifndef _SUBDIVISIONCLASSIFICATION_H
#define _SUBDIVISIONCLASSIFICATION_H


namespace Go {


/// This enum classifies the type of subdivision

enum SubdivisionClassification {
    CANNOT_SUBDIVIDE = 0,
    DIVIDE_DEG,
    DIVIDE_CRITICAL,
    DIVIDE_HIGH_SING,
    DIVIDE_SING,
    DIVIDE_KNOT,
    DIVIDE_INT,
    DIVIDE_PAR
};


} // namespace Go


#endif // _SUBDIVISIONCLASSIFICATION_H

