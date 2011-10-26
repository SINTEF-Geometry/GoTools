//==========================================================================
//                                                                          
// File: LinkType.h                                                          
//                                                                          
// Created: Tue Sep 26 13:16:40 2006                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: LinkType.h,v 1.4 2006-10-10 15:19:51 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================

#ifndef _LINKTYPE_H
#define _LINKTYPE_H


namespace Go {


/// A classification of IntersectionLink objects. The classification
/// is strictly based on how it was made within the algorithm.

enum LinkType {
    /// 0
    LINK_UNDEFINED = 0,
    /// 1
    SIMPLE_CONE,
    /// 2
    SIMPLE_MONOTONE,
    /// 3
    SIMPLE_IMPLICIT,
    /// 4
    SIMPLE_TWO_POINTS,
    /// 5
    COINCIDENCE_CVPT,
    /// 6
    COINCIDENCE_CVCV,
    /// 7
    COINCIDENCE_SFCV,
    /// 8
    COINCIDENCE_SFPT,
    /// 9
    MICRO_CVCV,
    /// 10
    MICRO_SFPT,
    /// 11
    MICRO_SFCV,
    /// 12
    MICRO_SFSF,
    /// 13
    MICRO_PAR1FUNC,
    /// 14
    MICRO_PAR2FUNC,
    /// 15
    LINEAR_CVCV,
    /// 16
    LINEAR_SFCV,
    /// 17
    LINEAR_SFSF,
    /// 18
    MERGED_UNDEFINED,
    /// 19
    MERGED_SIMPLE_CONE,
    /// 20
    MERGED_SIMPLE_MONOTONE,
    /// 21
    MERGED_SIMPLE_CONE_MONOTONE,
    /// 22
    MERGED_COINCIDENCE_SFCV_SFCV,
    /// 23
    DEG_TRIANGLE,
    /// 24
    POST_ITERATE,
    /// 25
    BRANCH_CONNECTION,
    /// 26
    COMPLEX_SFSF,
    /// 27
    SPLIT_LINK,
    /// 28
    REPAIRED_MISSING_LINK,
    /// 29
    INSIDE_OUTSIDE_SINGULARITY_BOX

};


} // namespace Go


#endif // _LINKTYPE_H

