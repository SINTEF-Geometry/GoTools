//===========================================================================
//                                                                           
// File: BoundedCurve.h                                                      
//                                                                           
// Created: Fri Aug 28 17:08:23 2009                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _BOUNDEDCURVE_H
#define _BOUNDEDCURVE_H


#include "GoTools/geometry/ParamCurve.h"


// Both parameter values and end points may be given to define the
// boundaries. Assuming that both points prefer parameter, or both
// points prefer points.

namespace Go
{

class GO_API BoundedCurve : public ParamCurve
{
public:

    /// Default constructor. Constructs an uninitialized Line which
    /// can only be assigned to or read into.
    BoundedCurve()
    {};

    /// Constructor. Input is start point and end point. Assumed to
    /// lie on curve (or at least close to it).
    BoundedCurve(shared_ptr<ParamCurve> curve, bool prefer_bd_par,
		 double start_par, double end_par,
		 Point start_pt, Point end_pt);

    BoundedCurve(shared_ptr<ParamCurve> curve,
		 Point start_pt, Point end_pt);

    BoundedCurve(shared_ptr<ParamCurve> curve,
		 double start_par, double end_par);

    /// virtual destructor - ensures safe inheritance
    virtual ~BoundedCurve();

    /// Read object from stream
    /// \param is stream from which object is read
    virtual void read (std::istream& is);
    /// Write object to stream
    /// \param os stream to which object is written
    virtual void write (std::ostream& os) const;


    // --- Functions inherited from GeomObject ---

    virtual BoundingBox boundingBox() const;
    
    virtual int dimension() const;

    virtual ClassType instanceType() const;

    static ClassType classType();

    virtual BoundedCurve* clone() const;

    // --- Functions inherited from ParamCurve ---

    virtual void point(Point& pt, double tpar) const;

    virtual void point(std::vector<Point>& pts, 
		       double tpar,
		       int derivs,
		       bool from_right = true) const;

    virtual double startparam() const;
    virtual double endparam() const;

    virtual void reverseParameterDirection(bool switchparam = false);
    
    virtual void setParameterInterval(double t1, double t2);

    virtual SplineCurve* geometryCurve();

    virtual bool isDegenerate(double degenerate_epsilon);

    virtual 
      BoundedCurve* subCurve(double from_par, double to_par,
			     double fuzzy = DEFAULT_PARAMETER_EPSILON) const;

    virtual DirectionCone directionCone() const;
 
    virtual void appendCurve(ParamCurve* cv, bool reparam=true);

    virtual void appendCurve(ParamCurve* cv,
			     int continuity, double& dist, bool reparam=true);

    virtual void closestPoint(const Point& pt,
			      double tmin,
			      double tmax,
			      double& clo_t,
			      Point& clo_pt,
			      double& clo_dist,
			      double const *seed = 0) const;

    virtual double length(double tol);

    /// Set bounds for the parametrization of the Line.
    /// \param startpar start parameter
    /// \param endpar end parameter
    void setParamBounds(double startpar, double endpar);

 private:
    shared_ptr<ParamCurve> curve_;

    bool prefer_parameter_; // As opposed to points.

    double startparam_;
    double endparam_;
    Point start_pt_;
    Point end_pt_;

//     // Also give an orientation?
//     bool opp_dir_;

};


} // namespace Go


#endif // _BOUNDEDCURVE_H

