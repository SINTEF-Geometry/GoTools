//===========================================================================
//                                                                           
// File: doxymain.h                                                          
//                                                                           
// Created: Wed Nov 16 11:07:09 2005                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: doxymain.h,v 1.2 2005/12/20 12:13:40 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _DOXYMAIN_H
#define _DOXYMAIN_H

//===========================================================================
//                        DOCUMENTATION ON MAIN PAGE
//===========================================================================

/// \mainpage GoTools Library
/// \section d0 Introduction
///
/// The newest version of GoTools is available from <a
/// href="http://www.sintef.no/Geometry-Toolkits">http://www.sintef.no/Geometry-Toolkits</a>.
///
/// GoTools is licensed under the <a
/// href="http://www.gnu.org/copyleft/gpl.html">GNU General Public
/// License </a>
///
/// GoTools is the group name of many interdependent C++ software
/// modules developed by the geometry group at SINTEF ICT, Dept. of
/// Applied Mathematics.  GoTools software has been developed for a
/// range of different applications in many different projects.
/// However, a few key modules are used by almost all the others;
/// these have been grouped together in \ref geometry_doc.
///
/// At the moment, the following GoTools modules are offered with GPL license:
/// \li \ref geometry_doc. Parametric curves and surfaces including construction 
/// methods and operations on these entities
/// \li \ref parametrization.  Parametrization of scattered data points
/// \li \ref implicitization. Approximating spline curves and surfaces by 
/// implicitely defined algebraic entities
/// \li \ref intersections. Intersection functionality involving spline curves and 
/// surfaces and computations of self intersections.
/// \li \ref igeslib. Read from and write to iges format. 
/// \li \ref trivariate. Spline volumes and elementary volumes
/// including some creation methods
/// \li \ref topology. Adjacency analysis for surface sets
/// \li \ref compositemodel. Representation of a surface set including topological
/// entities and operations on a set of surfaces as one unit.
/// \li \ref trivariatemodel. A volume model including topology and some operations 
/// on the model
/// \li \ref viewlib. A utility viewer to visualize curves and surfaces
///
/// GoTools depends on:
/// \li <a href="http://www.sintef.no/SISL">SISL</a> (available with
/// GPL license), SINTEF's spline library, for various spline related functionality,
/// \li <a href="http://www.sintef.no/Projectweb/Geometry-Toolkits/TTL">TTL</a> (available with GPL license), a generic triangulation library,
/// \li <a href="http://www.robertnz.net">newmat</a> (available with a permissive license), for various matrix operations.
///
/// For convenience, these libraries are included in the GPL version of GoTools.
///
/// \section building Building GoTools
///
/// This GoTools package uses CMake to generate a Makefile (on Linux)
/// or MS Visual Studio project file (on Windows).
/// 
/// For information on using CMake, see <a href="http://www.cmake.org">www.cmake.org</a>.
/// 
/// As a Quick Start Guide, on Linux, make a build directory somewhere:
/// 
/// \verbatim
$ cd some_dir
$ mkdir build
$ cd build
$ ccmake <path_to_source_code>
\endverbatim
///
/// Follow the instructions of 'ccmake' - the CMake "GUI". Then:
/// 
/// \verbatim
$ make
$ sudo make install
\endverbatim
///
/// On Windows, add a new build folder somewhere. Start the CMake
/// executable and fill in the paths to the source and build folders. When
/// you run CMake, a Visual Studio project solution file will be generated
/// in the build folder.
/// 
/// \subsection compilers Compilers
///
/// The code uses certain features of the new C++ standard C++11, most
/// notably the smart pointer \c std::shared_ptr. It has been tested
/// on GCC 4.6.1 on Linux and Visual Studio 2010 on Windows.
///
/// A set of options to control the build can be accessed in CMake
/// (names starting with \c GoTools). For example, you can turn on/off
/// building the various modules by checking/unchecking \c
/// GoTools_COMPILE_MODULE_<modulename>.
///
/// Also provided is the option \c GoTools_USE_BOOST. If this option
/// is turned on, the building process uses \c boost::shared_ptr
/// instead of \c std::shared_ptr. If a C++11 compliant compiler is
/// not available, you may try this option and see if it
/// works. Requires Boost: <a
/// href="http://www.boost.org">www.boost.org</a>.



#endif // _DOXYMAIN_H

