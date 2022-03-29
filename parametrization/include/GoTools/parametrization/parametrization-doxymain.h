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

#ifndef _PARAMETRIZATION_DOXYMAIN_H
#define _PARAMETRIZATION_DOXYMAIN_H

//===========================================================================
//                        DOCUMENTATION ON MAIN PAGE
//===========================================================================
/// \page parametrization GoTools Parametrization
///
/// \section d0 Purpose
/// The purpose of the GoTools Parametrization Library is to compute a 
/// reasonable parametrization for a discrete geometrical object.  This 
/// geometrical object can be a planar graph embedded in \f$R^2\f$ or
/// \f$R^3\f$ (triangulation, rectangular grid, etc.) or a point cloud without 
/// any explicit connectivity, but with an identified border.
///
/// There are many contexts for which such a parametrization is useful, for
/// example when mapping textures onto polygonal meshes or when triangulating
/// a surface from a set of scattered 3D points.
/// The algorithms contained in this software are based on the work of 
/// Dr. Michael Floater.
///
/// \section d1 How it can be used
/// This library has three degrees of freedom for specifying the 
/// parametrization:
/// <ul>
/// <li> choice of <em> data structure </em></li>
/// <li> choice of <em> boundary parametrization </em></li>
/// <li> choice of <em> parametrization of the interior </em> </li>
/// </ul>
/// \subsection d11 Choice of data structure
/// As mentioned in the last section, the data to be parametrized can be of
/// various kinds, from general polygonal meshes (3D points and 
/// connectivity) to point clouds (no connectivity information, but boundary
/// nodes must be indentified and ordered).  The various possible data types
/// are given as sub-classes of \link PrOrganizedPoints \endlink.
/// You should choose the one that corresponds to the data that you want
/// to parametrize.  Many of these classes have read/write functionality for
/// a generic format; you might want to reimplement these to correspond with
/// the format of the data you are using.  
/// The list of currently available data structures are:
/// <ul>
/// <li> \link PrUnorganized_OP \endlink </li>
/// <li> \link PrFastUnorganized_OP  \endlink </li>
/// <li> \link PrTriangulation_OP \endlink </li>
/// <li> \link PrRectangularGrid_OP \endlink </li>
/// <li> \link PrPlanarGraph_OP \endlink </li>
/// </ul>
/// You can of course write a totally <em>new</em> specialization of the 
/// PrOrganizedPoints class, as long as you implement all the required member
/// functions.
/// \subsection d12 Choice of boundary parametrization
/// The boundary of the surface must be parametrized before the surface
/// interior.  This is done by the class \link PrParametrizeBdy \endlink.
/// There are several ways of parametrizing the boundary.
/// You can specify the method you want to use through the 
/// PrParametrizeBdy::setParamKind() member function.  There are currently
/// three options: 
/// <ul>
/// <li>\c PrCHORDLENGTHBDY  - parametrize by chordlength (generally preferred)</li>
/// <li>\c PrCENTRIPETAL - parametrize using square root of boundary point distances</li>
/// <li>\c PrUNIFBDY - parametrize by imposing a fixed distance between each boundary point</li>
/// </ul>
/// Another thing that must be specified when determining the parametrization of
/// the boundary is the <em> shape </em> of the parametrical domain.  It can be 
/// a circle (default), or a rectangle (corners must be specified).
/// \subsection d13 Choice of interior parametrization
/// The parametrization of the interior can also be done in several ways, which are
/// each implemented as a sub-class of PrParametrizeInt.  Currently available methods are:
/// (click on the links for individual explanations)
/// <ul>
/// <li>PrPrmMeanValue (generally preferred)</li>
/// <li>PrPrmEDDHLS</li>
/// <li>PrPrmLeastSquare</li>
/// <li>PrPrmShpPres</li>
/// <li>PrPrmUniform</li>
/// </ul>
/// \subsection d14 Putting it all together
/// The general procedure for parametrizing a data set can be summarized as follows:
///   - Prepare the data structure\n
///      - establish an object of the appropriate datastructure 
///        (sub-class of PrOrganizedPoints).  This object should be owned by a 
///        shared pointer (as found in the <a href="http://www.boost.org">boost</a> library.
///      - read your data into the structure, either by using one of the provided member 
///        functions like PrTriangulation_OP::scanRawData(), PrTriangulation_OP::scan(),
///        or by adding your own I/O-routines.
///   - Parametrize the boundary\n
///      - Establish a PrParametrizeBdy object.
///      - Associate it with your data by using the PrParametrizeBdy::attach() function.
///      - Choose parametrization type by calling the PrParametrizeBdy::setParamKind() 
///        function.
///      - Execute the parametrization by calling PrParametrizeBdy().  When calling this 
///        function without arguments, the parametrical domain will be a circle.  
///        A rectangular parameter domain can also be used.  This is specified by additional
///        arguments to the function.
///   - Parametrize the interior\n
///      - When the boundary has been successfully parametrized, the interior can be
///        treated.  First, you need to instantiate an object of the desired sub-class of 
///        PrParametrizeInt.
///      - Associate your data to this object by using the PrParametrizeInt::attach()
///        function.
///      - Carry out the parametrization by calling the PrParametrizeInt::parametrize()
///        function.
///   - Accessing the result
///      - The parametrization can be accessed in several ways, through the PrOrganizedPoints
///        interface of your data structure:
///        - Several print routines can be used to dump all nodes with/without coordinates
///          and associated parameters (PrOrganizedPoints::printUVNodes, 
///          PrOrganizedPoints::printUVXYZNodes, etc.).
///        - The parameters of an individual node can be accessed by the 
///          PrOrganizedPoints::getU() and PrOrganizedPoints::getV() functions.
///      - Most of the sub-classes of PrOrganizedPoints also have print() and scan() 
///        functions that can be used for reading and writing their data to a stream.
///      - You can write your own customised data access and streaming functions.
///       
/// \section d5 Example program
/// The procedure mentioned above can be observed in a working program located in the
/// \c examples folder.  The source file name is demo.C.  The program is run from the 
/// command line, where you can
/// also specify the data.  It can handle both point clouds and triangulations.  A sample
/// triangulation is provided in the \c data/ folder.
///
/// For instance, try
/// \verbatim $ ./demo -i data/gjoevik_triang \endverbatim
/// (There are also more options you can specify; to see the list of options, run the 
/// program without any arguments).
///
/// If the input data is a triangulation (the case above), the resulting, parametrized
/// triangulation will be saved in the file \c triangulation.  In any case, the 
/// parametrization for each node/point will be dumped to the file \c uv_nodes.
///
/// Good luck!


#endif // _PARAMETRIZATION_DOXYMAIN_H

