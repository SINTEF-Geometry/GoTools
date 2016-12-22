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

#ifndef _CURVEONVOLUME_H
#define _CURVEONVOLUME_H

#include "GoTools/geometry/ParamSurface.h"
#include <memory>
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/trivariate/ParamVolume.h"
#include "GoTools/utils/config.h"

namespace Go
{

    /** \brief A curve living on a parametric volume. It either has got 
	information about the curve in geometry space and in the parameter
	domain of the volume or the ability to compute the other representation 
	given one.
     *
     */

class GO_API CurveOnVolume : public ParamCurve
{
public:
    /// Define an empty CurveOnVolume that can be assigned or \ref read()
    /// into
    CurveOnVolume();

    /// Construct a CurveOnVolume by specifying the volume and either a space curve
    /// or a curve in the parameter plane.  Does not clone any of the input, just 
    /// sets (smart) pointers.
    /// \param vol pointer to the underlying volume
    /// \param curve pointer to the curve specifying the CurveOnVolume.  This
    ///              curve may either be in the parametric domain of the volume
    ///              (3D curve), or a space curve that the user has assured 
    ///              to be coincident with the volume.  'preferparameter' specifies
    ///              which kind of curve this is.
    /// \param preferparameter if this is set to 'true', then 'curve' is assumed
    ///                        to be a curve in the parametric domain of the surface.
    ///                        Otherwise, it is assumed to be a space (3D) curve.
  CurveOnVolume(shared_ptr<ParamVolume> vol,
		shared_ptr<ParamCurve> curve,
		bool preferparameter);
    

  CurveOnVolume(shared_ptr<ParamVolume> vol,
		shared_ptr<ParamCurve> pcurve,
		shared_ptr<ParamCurve> spacecurve,
		bool preferparameter);
    

    /// Copy constructor.  The copy constructor will not clone() the underlying
    /// volume, but it will clone() both the parametric and the spatial curve.
    /// \param volume_curve the CurveOnVolume to copy into 'this' CurveOnVolume.
    CurveOnVolume(const CurveOnVolume& volume_curve);

    /// Assignment operator.  Like the copy constructor, the assignment operator
    /// clone()s the curves, and not the volume.
    /// \param other the CurveOnVolume to copy into 'this' CurveOnVolume.
    CurveOnVolume& operator= (const CurveOnVolume& other);
    
    /// Destructor.
    /// Trivial because memory is managed by shared_ptr.
    virtual ~CurveOnVolume();


    // inherited from Streamable
    virtual void read (std::istream& is);
    virtual void write (std::ostream& os) const;

    // inherited from GeomObject
    /// Axis align box surrounding this object
    /// Computed with respect to the space curve if this one exists,
    /// otherwise the underlying surface
    virtual BoundingBox boundingBox() const;

    /// Cone surrounding the set of tangent directions corresponding to
    /// the surface. Only computed if the space curve exists
    virtual DirectionCone directionCone() const;

    /// Dimension of geometry space
    virtual int dimension() const;

    /// Type of current object. 
    virtual ClassType instanceType() const;
    static ClassType classType()
    { return Class_CurveOnVolume; }
    
    // The function clone() calls the copy constructor,
    // so clone() also makes deep copies of the curves,
    // but not the underlying surface.
    virtual CurveOnVolume* clone() const
    { return new CurveOnVolume(*this); }

    // inherited from ParamCurve
    virtual void point(Point& pt, double tpar) const;

    /// Inherited from \ref ParamCurve. Only works for 'derivs' = 0 or 1.
    /// \see ParamCurve::point()
    virtual void point(std::vector<Point>& pts, 
		       double tpar,
		       int derivs, bool from_right = true) const;
    
    // inherited from ParamCurve
    virtual double startparam() const;

    // inherited from ParamCurve
    virtual double endparam() const;

    // inherited from ParamCurve
    virtual void reverseParameterDirection(bool switchparam = false);

    // inherited from ParamCurve
    virtual void setParameterInterval(double t1, double t2);

    // inherited from ParamCurve
    virtual SplineCurve* geometryCurve();

    // inherited from ParamCurve
    virtual bool isDegenerate(double degenerate_epsilon);

    // inherited from ParamCurve
    virtual CurveOnVolume* subCurve(double from_par, double to_par,
				       double fuzzy =
				       DEFAULT_PARAMETER_EPSILON) const;

    /// Split curve in a specified parameter value
    virtual
      std::vector<shared_ptr<ParamCurve> > 
      split(double param,
	    double fuzzy = DEFAULT_PARAMETER_EPSILON) const; 

   // inherited from ParamCurve
    virtual void closestPoint(const Point& pt,
			      double         tmin,
			      double         tmax,
			      double&        clo_t,
			      Point&       clo_pt,
			      double&        clo_dist,
			      double const   *seed = 0) const;

    // Inherited from ParamCurve
    // Compute the total length of this curve
    virtual double length(double tol);

    // inherited from ParamCurve.  NB: does not check whether the resulting ParamCurve
    // stays inside parameter domain (or that the space curve stays on surface).
    virtual void appendCurve(ParamCurve* cv, bool reparam=true);

    // inherited from ParamCurve.  NB: does not check whether the resulting ParamCurve
    // stays inside parameter domain (or that the space curve stays on surface).
    virtual void appendCurve(ParamCurve* cv,
			     int continuity, double& dist, bool reparam=true);

    /// Set the underlying surface to the one pointed to by the argument
    /// \param surface the pointer to the surface we will set as underlying for this
    ///                CurveOnVolume.
    void setUnderlyingVolume(shared_ptr<ParamVolume> volume)
    {volume_ = volume;}

    /// Replace the curves describing the curve on volume. The curve preference
    /// is not changed
    void setCurves(shared_ptr<ParamCurve> spacecurve,
		   shared_ptr<ParamCurve> parametercurve)
    {
      spacecurve_ = spacecurve;
      pcurve_ = parametercurve;
    }

    /// Replace the space curve corresponding to this curve on volume curve.
    /// Use with care!
    void setSpaceCurve(shared_ptr<ParamCurve> spacecurve)
    {
      spacecurve_ = spacecurve;
    }

    /// Replace the parameter curve corresponding to this curve on volume curve.
    /// Used for instance in relation to reparameterizations of the related volume.
    /// Use with care!
    void setParameterCurve(shared_ptr<ParamCurve> parametercurve)
    {
      pcurve_ = parametercurve;
    }

    /// Remove parameter curve information
    void unsetParameterCurve()
    {
      pcurve_.reset();
      prefer_parameter_ = false;
    }

    /// Remove space curve information
    void unsetSpaceCurve()
    {
      spacecurve_.reset();
      prefer_parameter_ = true;
    }

    /// Inherited from \ref ParamCurve.  If the parametric curve is set to be the
    /// 'prefered' one, this function will return the next segment value for the 
    /// parametric curve; otherwise it will return the next segment value for the 
    /// spatial 3D curve.
    /// See also \ref ParamCurve::nextSegmentVal()
    virtual double nextSegmentVal(double par, bool forward, double tol) const;


    /// Get a shared pointer to the underlying volume
    /// \return a shared pointer to the underlying volume
    shared_ptr<ParamVolume> underlyingVolume()
    { return volume_; }

    /// Get a shared pointer to the curve in the parameter domain.
    /// \return a shared pointer to the curve in the parameter domain
    shared_ptr<ParamCurve> parameterCurve()
      { return pcurve_; }

    /// Get a shared pointer to the space curve
    /// \return a shared pointer to the space curve.
    shared_ptr<ParamCurve> spaceCurve()
    { return spacecurve_; }

    /// Get a constant, shared pointer to the underlying volume
    /// \return a const-pointer to the underlying volume
    shared_ptr<const ParamVolume> underlyingVolume() const
    { return volume_; }

    /// Get a constant, shared pointer to the curve in the parameter domain.
    /// \return a const-pointer to the curve in the parameter domain.
    shared_ptr<const ParamCurve> parameterCurve() const
    { return pcurve_; }

    /// Get a constant, shared pointer to the space curve
    /// \return a const-shared pointer to the space curve.
    shared_ptr<const ParamCurve> spaceCurve() const
    { return spacecurve_; }
    
    /// Query whether the parameter curve or the space curve is prefered for computation
    /// in this object.
    bool parPref() const
    { return prefer_parameter_; }

    void setParPref(bool prefer) 
    { prefer_parameter_ = prefer; }

     /// Get the rectangle enclosing the underlying volume's parametric domain.
    /// \return the RectDomain for the underlying volume.
    RectDomain containingDomain() const;

    /// Fetch the volume parameter corresponding to a curve parameter
    Point volumeParameter(double crv_par,
			  const RectDomain* domain_of_interest = NULL) const;


protected:
    /// The underlying volume
    shared_ptr<ParamVolume> volume_;
    /// The 2D curve in the parameter domain of the volume.
    /// May point to null.
    shared_ptr<ParamCurve> pcurve_;
    /// An instance of the curve in the volume. May point to null.
    shared_ptr<ParamCurve> spacecurve_;
    /// Which representation to prefer if both exist
    bool prefer_parameter_;


};


} // namespace Go

#endif // _CURVEONVOLUME_H

