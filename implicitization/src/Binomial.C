//==========================================================================
//                                                                          
// File: Binomial.C                                                          
//                                                                          
// Created: Mon Nov  4 17:21:18 2002                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: Binomial.C,v 1.4 2003-06-18 13:21:44 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include "GoTools/implicitization/Binomial.h"


using namespace std;


namespace Go {


//===========================================================================
void Binomial::init(int n)
//===========================================================================
{
    pascals_triangle.resize(n+1);
    for (int nr = 0; nr < n+1; ++nr) {
	pascals_triangle[nr].resize(nr+1);
	pascals_triangle[nr][0] = 1.0;
	for (int j = 1; j < nr; ++j) {
	    pascals_triangle[nr][j]
		= pascals_triangle[nr-1][j-1]
		+ pascals_triangle[nr-1][j];
	}
	pascals_triangle[nr][nr] = 1.0;
    }

    return;
}


//===========================================================================
void Binomial::expand(int n)
//===========================================================================
{
    int old_size = (int)pascals_triangle.size();
    pascals_triangle.resize(n+1);

    for (int nr = old_size; nr < n+1; ++nr) {
	pascals_triangle[nr].resize(nr+1);
	pascals_triangle[nr][0] = 1.0;
	for (int j = 1; j < nr; ++j) {
	    pascals_triangle[nr][j]
		= pascals_triangle[nr-1][j-1]
		+ pascals_triangle[nr-1][j];
	}
	pascals_triangle[nr][nr] = 1.0;
    }

    return;
}


//===========================================================================


} // namespace Go
