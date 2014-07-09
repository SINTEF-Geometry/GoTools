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

#ifndef _COMPOSITEMODELFACTORY_H
#define _COMPOSITEMODELFACTORY_H

#include "GoTools/utils/config.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeCurve.h"
#include "GoTools/igeslib/IGESconverter.h"
#include "GoTools/utils/Point.h"
#include <vector>

struct SISLSurf;

namespace Go
{

  class ftSurface;
  class ParamCurve;

//===========================================================================
/** Factory class for creating children of CompositeModel
 */
//===========================================================================

class GO_API CompositeModelFactory
{
 public:
  /// Constructor
  CompositeModelFactory(double approxtol,
			double gap,   // Gap between adjacent surfaces
			double neighbour,  // Threshold for whether surfaces are adjacent
			double kink,  // Kink between adjacent surfaces 
			double bend); // Intended G1 discontinuity between adjacent surfaces

  /// Destructor
  ~CompositeModelFactory();

  /// Create an empty model
  SurfaceModel* createEmpty();

  /// Read IGES file
  /// One model is produced
  CompositeModel* createFromIges(std::istream& is, bool use_filetol=false,
				 bool prefer_surfacemodel=true);

  /// Get all connected models in the IGES file
  std::vector<shared_ptr<CompositeModel> > 
    getModelsFromIges(std::istream& is, bool use_filetol=false);

  /// Read G2 file
  /// One model is produced
  CompositeModel* createFromG2(std::istream& is, bool prefer_surfacemodel=true);

  /// Get all connected models in the g2 file
  std::vector<shared_ptr<CompositeModel> > 
    getModelsFromG2(std::istream& is, bool use_filetol=false);

  /// Read a vector of sisl surfaces
  SurfaceModel* createFromSisl(std::vector<SISLSurf*>& surfaces);

  /// Read a vector of sisl curves
  CompositeCurve* createFromSisl(std::vector<SISLCurve*>& curves);

  /// Make surface model from a box
  /// Input is one box corner, the vector from this corner towards another corner in the
  /// same box side, yet another vector in the plane defining this box side and the lengts
  /// of the box sides. The first length corresponds to the vector defining the box side,
  /// the second to the other side in the defined plane and the third defines the depth of
  /// the box. The lengths may be negative.
  SurfaceModel* createFromBox(Point corner, Point side_vec, Point plane_vec,
				double side1_length, double side2_length, 
				double side3_length);

  /// Make surface model from sphere
  /// Input is sphere center and radius
  SurfaceModel* createFromSphere(Point centre, double radius);

  /// Make surface model from octants of a sphere
  /// Input is sphere center, axis toward norht pole, vector from centre to
  /// point on equator and the number of octants. The length on the equator vector
  /// gives the sphere radius
  /// latitude = 1 : Octants in northern hemisphere
  /// latitude = 2 : Octants in both hemispheres
  /// longitude = 1 : Octants in first quadrant
  /// longitude = 2 : Octants in first and second quadrant
  /// longitude = 3 : Octants in first, second and third quadrant
  /// longitude = 4 : Octants in all quadrants
  SurfaceModel* createFromSphere(Point centre, Point axis, 
				   Point equator, int latitude, int longitude);

  /// Make surface model from truncated cylinder
  /// bottom_pos : centre of cylinder at bottom
  /// cylinder_axis : Cylinder axis. The length of the vector defines the cylinder height
  /// major_axis : Major axis in ellipse at cylinder bottom. The vector length defines the major radius
  /// minor_axis : Minor axis in ellipse at cylinder bottom. The vector length defines the minor radius
  SurfaceModel* createFromCylinder(Point bottom_pos, Point cylinder_axis,
				     Point major_axis, Point minor_axis);

  /// Interpolate a set of positional curves. Automatic parameterization. The
  /// parameter values corresponding to the curves are given as output
  SurfaceModel* 
    interpolateCurves(const std::vector<shared_ptr<SplineCurve> >& curves,
		      std::vector<double>& parvals, 
		      int open, int degree = 3);

  /// Interpolate a set of positional curves. Parameterization is given.
  /// If a closed surface is specified, one extra parameter value must be
  /// specified
  SurfaceModel* 
    interpolateCurves2(const std::vector<shared_ptr<SplineCurve> >& curves,
		       std::vector<double>& param,
		       int open, int degree = 3);

  /// Interpolate a set of positional curves where a tangent curve may be added. 
  //  Automatic parameterization. The parameter values corresponding to the 
  /// positional curves are given as output
  /// NB! Rational input curves are not handled
  SurfaceModel* 
    interpolateCurves(const std::vector<shared_ptr<SplineCurve> >& curves,
		      std::vector<int>& crv_type,
		      std::vector<double>& parvals,
		       int open, int degree = 3);

  /// Interpolate a set of positional curves where a tangent curve may be added. 
  /// Parameterization is given and corresponds to the positional curves.
  /// NB! Rational input curves are not handled
  SurfaceModel* 
    interpolateCurves2(const std::vector<shared_ptr<SplineCurve> >& curves,
		       std::vector<int>& crv_type,
		       std::vector<double>& param,
		       int open, int degree = 3);


  /// Create a circular arc
  /// centre : Circle centre
  /// start_pt : Start point of arc
  /// angle : Angle of arc
  /// axis : Axis of the plane in which the circle lies
  CompositeCurve* createCircularArc(Point centre, Point start_pt, double angle,
				    Point axis);

  /// Create an elliptic arc
  /// centre : Ellipse centre
  /// direction : Direction of major ellips axis
  /// r1 : Major radius
  /// r2 : Minor radius
  /// startpar : Start parameter of segment
  /// angle : Angle of arc
  /// axis : Axis of the plane in which the ellipse lies
  CompositeCurve* createEllipticArc(Point centre, Point direction, 
				    double r1, double r2,
				    double startpar, double angle,
				    Point axis);

  /// Create line segment bewteen the two points startpt and endpt
  CompositeCurve* createLineSegment(Point startpt, Point endpt);

 private:
  double approxtol_;
  double gap_;        // Gap between adjacent surfaces
  double neighbour_;  // Threshold for whether surfaces are adjacent
  double kink_;       // Kink between adjacent surfaces 
  double bend_;       // Intended G1 discontinuity between adjacent surfaces

  // Read geometry from file converter
  CompositeModel* getGeometry(IGESconverter& conv, bool use_filetol,
			      bool prefer_surfacemodel);

  void getAllEntities(IGESconverter& conv, 
		      std::vector<shared_ptr<ftSurface> >& faces,
		      std::vector<shared_ptr<ParamCurve> >& curves);

  // Make spline surface from control points
  SplineSurface* fromKnotsAndCoefs(int order1, std::vector<double> knots1, int order2,
				   std::vector<double> knots2, vector<Point> coefs);

  void replaceElementaryCurves(shared_ptr<CurveOnSurface> sf_cv);
};

} // namespace Go

#endif // _COMPOSITEMODELFACTORY_H
