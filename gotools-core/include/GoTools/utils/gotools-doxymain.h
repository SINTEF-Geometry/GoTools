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

#ifndef _DOXYMAIN_H
#define _DOXYMAIN_H


/// \namespace Go
/// The Go namespace is the common namespace for all GoTools modules.



//===========================================================================
//                        DOCUMENTATION ON MAIN PAGE
//===========================================================================

/**
\mainpage GoTools Library

\section intromain Introduction
The newest version of GoTools is available from <a
href="http://www.sintef.no/Geometry-Toolkits">http://www.sintef.no/Geometry-Toolkits</a>.

GoTools is licensed under the <a
href="http://www.gnu.org/copyleft/gpl.html">GNU General Public
License </a>

GoTools is the group name of many interdependent C++ software
modules developed by the geometry group at SINTEF ICT, Dept. of
Applied Mathematics.  GoTools software has been developed for a
range of different applications in many different projects.
However, a few key modules are used by almost all the others;
these have been grouped together in \ref geometry_doc.

At the moment, the following GoTools modules are offered with GPL license:
\li \ref geometry_doc. Parametric curves and surfaces including construction 
methods and operations on these entities
\li \ref parametrization.  Parametrization of scattered data points
\li \ref implicitization. Approximating spline curves and surfaces by 
implicitely defined algebraic entities
\li \ref intersections. Intersection functionality involving spline curves and 
surfaces and computations of self intersections.
\li \ref igeslib. Read from and write to iges format. 
\li \ref trivariate. Spline volumes and elementary volumes
including some creation methods
\li \ref topology. Adjacency analysis for surface sets
\li \ref compositemodel. Representation of a surface set including topological
entities and operations on a set of surfaces as one unit.
\li \ref trivariatemodel. A volume model including topology and some operations 
on the model
\li \ref viewlib. A utility viewer to visualize curves and surfaces
\li \ref qualitymodule. A set of tools to check the quality of CAD models
\li \ref isogeometric_model. A set of tools related to isogeometric analysis
\li \ref lrsplines2d. LR spline surfaces
\li \ref lrsplines3d. LR spline volumes

GoTools depends on:
\li <a href="http://www.sintef.no/SISL">SISL</a> (available with
GPL license), SINTEF's spline library, for various spline related functionality,
\li <a href="http://www.sintef.no/Projectweb/Geometry-Toolkits/TTL">TTL</a> (available with GPL license), a generic triangulation library,
\li <a href="http://www.robertnz.net">newmat</a> (available with a permissive license), for various matrix operations.

For convenience, these libraries are included in the GPL version of GoTools.

\section building Building GoTools
This GoTools package uses CMake to generate a Makefile (on Linux)
or MS Visual Studio project file (on Windows).

For information on using CMake, see <a href="http://www.cmake.org">www.cmake.org</a>.

As a Quick Start Guide, on Linux, make a build directory somewhere:

\verbatim
$ cd some_dir
$ mkdir build
$ cd build
$ ccmake <path_to_source_code>
\endverbatim

Follow the instructions of 'ccmake' - the CMake "GUI". Then:

\verbatim
$ make
$ sudo make install
\endverbatim

On Windows, add a new build folder somewhere. Start the CMake
executable and fill in the paths to the source and build folders. When
you run CMake, a Visual Studio project solution file will be generated
in the build folder.

\subsection compilers Compilers
The code uses certain features of the new C++ standard C++11, most
notably the smart pointer \c std::shared_ptr. It has been tested
on GCC 4.6.1 on Linux and Visual Studio 2010 on Windows.

A set of options to control the build can be accessed in CMake
(names starting with \c GoTools). For example, you can turn on/off
building the various modules by checking/unchecking \c
GoTools_COMPILE_MODULE_<modulename>.

Also provided is the option \c GoTools_USE_BOOST. If this option
is turned on, the building process uses \c boost::shared_ptr
instead of \c std::shared_ptr. If a C++11 compliant compiler is
not available, you may try this option and see if it
works. Requires Boost: <a
href="http://www.boost.org">www.boost.org</a>.
*/


#endif // _DOXYMAIN_H

