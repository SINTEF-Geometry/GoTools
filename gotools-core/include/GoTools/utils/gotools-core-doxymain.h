//==========================================================================
//                                                                          
// File: gotools-core-doxymain.h                                              
//                                                                          
// Created: Fri Nov 11 13:33:54 2011                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
//==========================================================================

#ifndef _GOTOOLSCORE-DOXYMAIN_H
#define _GOTOOLSCORE-DOXYMAIN_H


//===========================================================================
//                        DOCUMENTATION ON GoTools Core PAGE
//===========================================================================

/// \page gotools-core GoTools Core Library
/// \section d0 Introduction
///
/// The GoTools Core Library is mainly a library for parametrical
/// curves and surfaces, usually spline based.  For clarity it is
/// divided into four parts. The two main parts are <em> 'geometry'
/// </em> and <em> 'utility' </em>.  The <em>'geometry' </em> part
/// contains classes and funcitonality for representing, storing and
/// manipulating parameterized geometrical objects, whereas the
/// 'utility' part contains general, low-level
/// functionality. <em>'creators'</em> and <em>'tesselator'</em>
/// contains functionality for generating curves and surfaces by
/// approximation, blending, etc., and for making tesselations,
/// respectively.
///
/// \section d1 The 'geometry' module
/// \subsection s1d1 Main, streamable classes
///
/// See \beginlink \link geometry_doc geometry\endlink
///
/// The main classes in the 'geometry' module are the \link
/// Go::SplineCurve SplineCurve \endlink and the \link
/// Go::SplineSurface SplineSurface \endlink, which are used for
/// reading, storing, manipulating NURBS curves and surfaces
/// respectively.  In addition to this, we have the \link
/// Go::BoundedSurface BoundedSurface \endlink, which represents a
/// <em>trimmed </em> surface (surface whose parametric domain has an
/// arbitrary shape and topology).  Other geometrical objects include
/// \link Go::PointCloud PointCloud \endlink representing a set of
/// (potentially very many) points, and \link Go::LineCloud LineCloud
/// \endlink, representing a set of lines.  All the abovementioned
/// objects are \link Go::Streamable Streamable \endlink, which means
/// that they can be written to and read from a stream (typically a
/// file) in a uniform way.  Similarly, they can be created in a
/// uniform way by the \link Go::Factory Factory \endlink, which is
/// useful when you want to generate objects whose exact kind are
/// unknown at compile time.
///
/// \subsection s2d1 Description of some of the core functionality
/// (non-exhaustive)
///
/// All the spline objects have member functions for basic tasks such
/// as evaluating points, tangents and higher-order derivatives,
/// changing parameterization (insertion of knots, rescaling or
/// reversing of parameter domains, etc.), bounding box specification,
/// order raising, conversion to Bezier form, reading and writing,
/// closeness to a given point in space, etc.  In addition,
/// namespace-level functions are providing for tasks such as:
///
/// <ul>
///    <li> 
///       closest point calculations 
///       <ul>
///          <li> \link Go::closestPtCurves2D between two curves in the plane \endlink </li>
///          <li> \link Go::closestPtCurves between two curves in space \endlink </li>
///          <li> \link Go::closestPtCurveSurf between a curve and a surface \endlink </li>
///          <li> \link Go::closestPtSurfSurfPlane between two surfaces and a plane \endlink </li>
///       </ul>
///    </li>
///    <li> computation of intersections </li>
///         <ul>
///          <li> \link Go::intersect2Dcurves between two curves in the plane \endlink </li>
///          <li> \link Go::intersectcurves between two curves in space \endlink </li>
///          <li> \link Go::BoundedUtils::intersectWithPlane intersect a surface and a plane \endlink</li>
///          <li> \link Go::closestPtCurveSurf between a curve and a surface \endlink </li>
///          <li> \link Go::closestPtSurfSurfPlane between two surfaces and a plane \endlink </li>
///         </ul>
///    <li> affine transformations </li>
///         <ul> 
///            <li> rotate \link Go::rotateSplineCurve curves \endlink, \link Go::rotateSplineSurf surfaces 
///                 \endlink, \link Go::BoundedUtils::rotateBoundedSurf bounded surfaces \endlink, \link
///                 Go::rotatePoint points \endlink and \link Go::rotateLineCloud clouds \endlink. </li>
///            <li> translate  \link Go::translateSplineCurve curves \endlink, \link Go::translateSplineSurf
///                 surfaces \endlink, \link Go::BoundedUtils::translateBoundedSurf bounded surfaces \endlink
///                 and \link Go::translateLineCloud clouds \endlink. </li>
///         </ul>
///    </li>
///    <li> splitting and merging </li>
///         <ul>
///            <li> \link Go::splitCurveIntoSegments split a curve into Bezier segments \endlink </li>
///            <li> \link Go::splitSurfaceIntoPatches split a surface into Bezier patches \endlink </li>
///            <li> \link Go::joinPatches merge individual patches into a spline surface \endlink </li>
///         </ul>
///    </li>
/// </ul>
///
/// \section d2 The 'creators' module
///
/// The 'creators' module contains various methods for modifying and creating
/// spline curves and surfaces. In addition to standard functionality such as
/// projecting and lifting a spline curve, the class incorporates the
/// following methods:
/// <ul>
///    <li> 
///       Approximation
///       <ul>
///          <li> Approximate parametrized points with a spline curve within
///               a prescribed tolerance, whilst fulfilling boundary 
///               conditions ( \link Go::ApproxCurve ApproxCurve \endlink ).
///               </li>
///          <li> Approximate parametrized points with a spline surface within
///               a prescribed tolerance, whilst fulfilling boundary 
///               conditions ( \link Go::ApproxSurf ApproxSurf
///               \endlink ). </li>
///       </ul>
///    </li>
///    <li>
///       Smoothing
///       <ul>
///          <li> Perform smoothing on a spline curve while fulfilling
///               various constraints
///               (\link Go::SmoothCurve SmoothCurve \endlink). </li>
///          <li> Perform smoothing on a spline surface while fulfilling
///               various constraints
///               (\link Go::SmoothSurf SmoothSurf \endlink). </li>
///          <li> Perform smoothing on a set of spline surfaces while
///               fulfilling various constraints
///               (\link Go::SmoothSurfSet SmoothSurfSet \endlink). </li>
///       </ul>
///    </li>
///    <li>
///       Skinning
///       <ul>
///          <li> Given input of a set of curves, and possibly the
///               corresponding cross tangent curves, create a spline surface
///               interpolating the curves
///               (\link Go::CoonsPatchGen::loftSurface loftSurface \endlink).
///               </li>
///       </ul>
///    </li>
///    <li>
///       Coons patch
///       <ul>
///          <li> Given input of four boundary curves, and possibly their
///               corresponding cross tangent curves, create a spline surface
///               interpolating the curves
///               (\link Go::CoonsPatchGen::createCoonsPatch createCoonsPatch
///               \endlink). </li>
///       </ul>
///    </li>
///    <li>
///       Gordons surface
///       <ul>
///          <li> The generalization of 'Coons patch'
///               to also include inner curves as well as their corresponding
///               cross tangent curve
///               (\link Go::CoonsPatchGen::createGordonSurface
///               createGordonSurface \endlink). </li>
///       </ul>
///    </li>
///    <li>
///       Hahns surface
///       <ul>
///          <li> The generalization of 'Coons patch' to cover input of
///               3 to 6 boundary curves, and possibly their corresponding
///               cross tangent curve. The number of surfaces created 
///               equals the number of boundary curves. For the case
///               of 4 boundary curves Coons patch is usually a better 
///               hole-filling strategy ( \link Go::HahnsSurfaceGen::constructHahnsSurface HahnsSurface \endlink ). </li>
///       </ul>
///    </li>
/// </ul>
///
/// \section d3 The 'utility' module
/// The 'utility' module was developed in parallel with the 'geometry' module,
/// and contains various useful classes and functions related to math. / 
/// computational geometry, while not being directly related to splines.
///
/// \subsection s1d3 Data structures
/// You can find templated arrays and optimised vectors.
///
/// \subsection s2d3 Geometrical objects
/// The module contain some objects that are much used in computational geometry,
/// like
/// <ul>
/// <li> \link Go::Point points \endlink </li>
/// <li> \link Go::BoundingBox bounding boxes \endlink, axis-aligned or 
///      \link Go::RotatedBox rotated \endlink, eventually 
///      \link Go::CompositeBox composite \endlink</li>
/// <li> \link Go::Rational rational \endlink numbers </li>
/// <li> \link Go::BaryCoordSystem barycentric \endlink coordinate systems </li>
/// <li> \link Go::DirectionCone direction cones \endlink </li>
/// </ul>
/// \subsection s3d3 Misc.
/// <ul>
/// <li> Compile-time computation of \link Go::Factorial factorials \endlink </li>
/// <li> A general, nonlinear \link Go::FunctionMinimizer function minimizer \endlink  
///      in an arbitrary number of variables.  This can be useful in a great number
///      of different settings! </li>
/// </ul>
///
/// \section d4 The 'tesselator' module
/// This module is not properly documented yet...


//===========================================================================
//                        OTHER DOCUMENTATION
//===========================================================================

/// \namespace Go
/// The Go namespace is the common namespace for all GoTools modules.


//===========================================================================

//FUNCTIONALITY THAT SHOULD BE MENTIONED IN UTILS

// * Factorial

// * Arrays (template, scratch vectors, etc.)
// * scratchvect

// * point
// * barycentric coordinate sysems
// * rational
// * bounding boxes (composite box, rotated box)
// * direction cones

// * Function minimizer





//FUNCTIONALITY THAT SHOULD BE MENTIONED IN GEOMETRY

//----------classes-------------
// * ApproxCurve
// * ApproxSUrf

// * BsplineBasis

// * Factory
// * Interpolator (HermiteInterpolator, SplineInterpolator, SplineApproximator)
// * Point- and Line clouds

// * ParamSurface
// * ParamCurve
// * SplineCurve
// * CurveOnSurface
// * SplineSurface
// * BoundedSurface

// * Domain
// * CurveLoop

//----------functions-----------

// * 'rotate' functions, translate functions


#endif // _GOTOOLSCORE-DOXYMAIN_H

