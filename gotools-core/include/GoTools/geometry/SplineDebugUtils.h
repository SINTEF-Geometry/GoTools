//===========================================================================
//                                                                           
// File: SplineDebugUtils.h                                                  
//                                                                           
// Created: Tue Sep 27 11:02:34 2005                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: SplineDebugUtils.h,v 1.2 2009-05-13 07:32:30 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _SPLINEDEBUGUTILS_H
#define _SPLINEDEBUGUTILS_H


#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include <memory>
#include "GoTools/utils/config.h"
#include "GoTools/geometry/BoundedSurface.h"


namespace Go
{

    /// For debugging. Writes a parameter curve in the xy-plane. Remove when
    /// GoViewer handles 2D curves.

    /// For debugging.  Writes a parameter curve (2D) in the xy-plane, for a given z-value.
    /// It will be written (with header) to the specified stream as a 3D curve.
    /// \param pcurve the parameter curve we want to write as a 3D curve
    /// \param os the stream to which we want to write the 3D curve
    /// \param z the constant z-value for the generated curve
    void GO_API writeSpaceParamCurve(const SplineCurve& pcurve,
                                     std::ostream& os, double z = 0.0);

    void GO_API writeTrimmedInfo(BoundedSurface& bd_sf,
				 std::ostream& os, double z = 0.0);

    /// writes the geometric object (with header) to the specified file name.
    /// \param geom_obj the object to write to file.
    /// \param to_file the file name to which the object will be written.
    void GO_API objToFile(GeomObject* geom_obj, char *to_file);

    /// writes the geometric objects (with header) to the specified file name.
    /// \param geom_objs the objects to write to file.
    /// \param to_file the file name to which the objects will be written.
    void GO_API objsToFile(std::vector<shared_ptr<GeomObject> >& geom_objs,
                           char *to_file);

    /// Write a SplineCurve to a stream using the SISL file format (not the Go format).
    /// \param spline_cv the curve to write to a stream
    /// \param os the stream to which the curve will be written (in SISL format)x
    void GO_API writeSISLFormat(const SplineCurve& spline_cv, std::ostream& os);

} // End of namespace Go


#endif // _SPLINEDEBUGUTILS_H

