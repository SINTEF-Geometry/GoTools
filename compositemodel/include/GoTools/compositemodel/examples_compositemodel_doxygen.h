#ifndef _EXAMPLES_COMPOSITEMODEL_DOXYMAIN_H
#define _EXAMPLES_COMPOSITEMODEL_DOXYMAIN_H

/**

\example createSplitDisc createSplitDisc.C 
\verbatim
\endverbatim
This programs creates a face set representing a disc as two trimmed
surfaces: The trimmed disc and a rectangular surface lying inside this
disc. This is a starting point for the creation of a disc as a multi patch
model with spline surfaces and no degeneracies.
The construction uses planar, rectangular surfaces and a truncated sylinder,
but the operations performed using these surfaces, do not depend on that
level of regularity.

\example createBlockStructuredDisc createBlockStructuredDisc.C
\verbatim
\endverbatim
This example file creates a block structured set of spline surfaces 
from a face set
consisting of possibly trimmed surfaces with arbitrary topology
(no corner-to-corner conditions). 

\example createVolumeBoundaries createVolumeBoundaries.C
\verbatim
\endverbatim
This programs creates a set of B-spline surfaces intended as the boundary
surfaces for a spline volume. A number of different methods are used in the
surface construction. Thus, this expample program illustrates some of the
possibilities for surface construction.

\example face2splineset face2splineSet.C
\verbatim
\endverbatim
This program demonstrates how to create a set of spline surfaces,
meeting in a corner-to-corner configuration and with corresponding
coefficients at common boundaries, from one possibly trimmed face
///
The program reads a bounded surface from a file, splits this surface
into several bounded surfaces where each surface has (at most) 4 boundary
curves. Finally, each bounded surface is approximated by a spline surface
within a given tolerance and C0 continuities at common boundaries is
ensured.
*/

#endif // _EXAMPLES_COMPOSITEMODEL_DOXYMAIN_H
