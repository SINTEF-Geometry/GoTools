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

#ifndef _SPLINEDEBUGUTILS_H
#define _SPLINEDEBUGUTILS_H


#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/Line.h"
#include "GoTools/geometry/Circle.h"
#include <memory>
#include "GoTools/utils/config.h"
#include "GoTools/geometry/BoundedSurface.h"


namespace Go
{

    /// For debugging. Writes a parameter curve in the xy-plane. Remove when
    /// GoViewer handles 2D curves.

namespace SplineDebugUtils
{
    /// For debugging.  Writes a parameter curve (2D) in the xy-plane, for a given z-value.
    /// It will be written (with header) to the specified stream as a 3D curve.
    /// \param pcurve the parameter curve we want to write as a 3D curve
    /// \param os the stream to which we want to write the 3D curve
    /// \param z the constant z-value for the generated curve
    void GO_API writeSpaceParamCurve(const SplineCurve& pcurve,
                                     std::ostream& os, double z = 0.0);

    void GO_API writeSpaceParamCurve(shared_ptr<ParamCurve> pcurve,
                                     std::ostream& os, double z = 0.0);

    void GO_API writeSpaceParamCurve(const Line& pline,
                                     std::ostream& os, double z = 0.0);

    void GO_API writeSpaceParamCurve(const Circle& pcirc,
                                     std::ostream& os, double z = 0.0);

    void GO_API writeSpace1DCurve(const SplineCurve& pcurve,
                                     std::ostream& os, double z = 0.0);

    /// Write the parameter curve (if existing) and the space curve (if
    /// existing) to the output stream. Both curves are written as 3D curves
    /// extending 2D curves with the given z-value.
    void GO_API writeTrimmedInfo(BoundedSurface& bd_sf,
				 std::ostream& os, double z = 0.0);

    void GO_API writeOuterBoundaryLoop(ParamSurface& sf,
				       std::ostream& os);

     void GO_API writeBoundary(BoundedSurface& sf, std::ostream& os);

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

    // Assuming the curves form a loop.
    void GO_API writeCvsOnSf(const std::vector<shared_ptr<Go::ParamCurve> >& loop_cvs,
			     double epsgeo,
			     std::ofstream& fileout);

    // Assuming the curves form a loop.
    void GO_API writeCvsOnSf(const std::vector<shared_ptr<Go::CurveOnSurface> >& loop_cvs,
			     double epsgeo,
			     std::ofstream& fileout);

    void GO_API writeSeamInfo(Go::BoundedSurface& bd_sf,
                              std::ofstream& fileout);


} // End of namespace SplineDebugUtils

} // End of namespace Go


#endif // _SPLINEDEBUGUTILS_H

