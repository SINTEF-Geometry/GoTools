/** Given a coarse vertex v of degree deg whose index is i
 * and two neighbouring fine vertices u1 and u2 which are the (j+1)-th and
 * (k+1)-th neighbours resp., in an anticlockwise direction,
 * (i.e.\ j and k begin theie numbering at 0)
 * (which may possibly be equal), return the function
 * theta(u1,u2,v), described in the paper: <br>
 * M. S. Floater and E. G. Quak,
 * Piecewise Linear Prewavelets on Arbitrary Triangulations,
 * Numer. Math. \b 82 (1999), 221--252.
 */

#ifndef _PR_WAVELET_UTIL
#define _PR_WAVELET_UTIL

double theta(int j, int k, int i, int deg, bool isBoundary);

#endif // _PR_WAVELET_UTIL

