//==========================================================================
//                                                                          
// File: SingularityClassification.h                                         
//                                                                          
// Created: Wed Sep 27 14:39:53 2006                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: SingularityClassification.h,v 1.1 2006-09-27 12:52:13 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================

#ifndef _SINGULARITYCLASSIFICATION_H
#define _SINGULARITYCLASSIFICATION_H


namespace Go {


/// This enum classifies the type of singularity

enum SingularityClassification {
    NO_SING = 0,
    INIT_SING,
    DIVIDED_SING,
    KEEP_SING
};


} // namespace Go

	    
#endif // _SINGULARITYCLASSIFICATION_H

