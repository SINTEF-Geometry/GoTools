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

#ifndef _REFINE_STRATEGIES_DOXYMAIN_H
#define _REFINE_STRATEGIES_DOXYMAIN_H
/**
 \page lrsplines2d_refine  LRsplines2D, strategies for refinement

Given an element to refine, a knot line segment to insert in the LR mesh
must be defined. According to the rules, this segment must cover the
support of at least one B-spline. The segment itself can cover the support or
it can extend an existing knot line segment partially covering the support
such that the combined segments will lead to splitting the B-spline. There
exists several strategies for defining knot line segments. One can
select a B-spline with the given element in its support and define one or more
segments covering the support of this B-spline or one can base the definition
of the knot line segment(s) on the element itself. The main strategies are:

- <b>Full span</b> Given an element selected for refinement, all B-splines overlapping this element
    are split in a parameter contained in the element.
- <b>Minimum span</b> The shortest possible knotline overlapping a chosen element that split
    at least one B-spline, is selected. Several candidate B-splines may exist and even if additional
    selection criteria are added, the candidate may not be unique. The selected knotline is not necessarily
    symmetric with respect to the element.
- <b>Structured mesh</b> Choose a B-spline and refine all knot intervals in this B-spline.
    No elements with an unbalanced aspect ratio will be created if the split is performed in the
    middle of the knot intervals.
- <b> Restricted mesh </b> Refine some knot intervals in a selected B-spline.
The strategy is mostly relevant in an approximation context. The knot 
intervals are selected by length and the approximation accuracy in the the
elements associated to a knot interval.

\subsection IGA Refinement strategies for IGA
Structured mesh is the preferred approximation strategy in IGA. It results
in a well behaved mesh, but the data size of the surface increases rapidly.
Structured mesh is known not to result in linearly dependent B-splines for
polynomial degrees one, two and three.

\subsection approx Refinement strategies for scattered data approximation
Refinement and reapproximation is performed iteratively. The procedure is 
as follows:
<ul>
<li> Define initial surface
<li> Check approximation accuracy
<li> Until a tolerance is reached or a maximum number of iterations is performed:
<ul>
<li> Identify elements or B-splines where the approximation is too poor
<li> Define new knot line segments
<li> Refine
<li> Approximate
<li> Check approximation accuracy
</ul>
</ul>

In approximation context, it is important to obtained the best possbiel 
accuracy with the least number of coefficients. Furthermore, the selected
refinement strategy should provide a stable good result. Experience shows
that the Full span strategy is a good choice in this context. Further
reduction in the number of coefficients can be obtained by refining one
parameter direction only in one iteration and alternating the directions.
This approach results in a sligthly less balanced mesh.

A further reduction on the speed in which new coeffients are introduced
can be obtained by not refining all elements or B-splines where the
approximation accuracy is not met. Parameters to guide the refinement process
can be set in LRSurfApprox::setRefinementStrategy, 
see \link Go::LRSurfApprox \endlink .

The Minimum span strategy leads to the least introduction of new coefficients
for each identified element to refine and can give good accuracy with few
coefficients in the beginning of the iterative approximation process. 
However, it has a tendency to get stuck adding coefficients without 
improving the accuracy later in the process. Moreover, this strategy is the
one most prone to lead to linear dependencies.

**/
#endif // _REFINE_STRATEGIES_DOXYMAIN_H
