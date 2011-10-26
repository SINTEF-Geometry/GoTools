//===========================================================================
//                                                                           
// File: BoundedCurve.C                                                      
//                                                                           
// Created: Fri Aug 28 17:08:37 2009                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/geometry/BoundedCurve.h"
#include "GoTools/geometry/Factory.h"
#include <vector>


using std::vector;
using std::shared_ptr;
using std::endl;


namespace Go {


//===========================================================================
BoundedCurve::BoundedCurve(shared_ptr<ParamCurve> curve,
			   bool prefer_parameter,
			   double startpar, double endpar,
			   Point start_pt, Point end_pt)
    : curve_(curve), prefer_parameter_(prefer_parameter),
      startparam_(startpar), endparam_(endpar),
      start_pt_(start_pt), end_pt_(end_pt)
//===========================================================================
{
}


// Constructor. Input is point and direction
//===========================================================================
BoundedCurve::BoundedCurve(shared_ptr<ParamCurve> curve,
			   Point start_pt,
			   Point end_pt)
    : curve_(curve), prefer_parameter_(false),
      start_pt_(start_pt), end_pt_(end_pt)
//===========================================================================
{
    startparam_ = curve_->startparam();
    endparam_ = curve_->endparam();
}


//===========================================================================
BoundedCurve::BoundedCurve(shared_ptr<ParamCurve> curve,
			   double startpar, double endpar)
  : curve_(curve), prefer_parameter_(true),
    startparam_(startpar), endparam_(endpar)
//===========================================================================
{
}

/// virtual destructor - ensures safe inheritance
//===========================================================================
BoundedCurve::~BoundedCurve()
//===========================================================================
{
}


/// Read object from stream
//===========================================================================
void BoundedCurve::read(std::istream& is)
//===========================================================================
{
    bool is_good = is.good();
    if (!is_good) {
	THROW("Invalid geometry file!");
    }

    int curve_type;
    is >> curve_type;

    is >> prefer_parameter_;
    is >> startparam_;
    is >> endparam_;

    int dim;
    is >> dim;
    start_pt_.resize(dim);
    end_pt_.resize(dim);
    int i;
    for (i=0; i<dim; ++i)
	is >> start_pt_[i];    
    for (i=0; i<dim; ++i)
	is >> end_pt_[i];

    ClassType type = ClassType(curve_type); // Needs this conversion
    shared_ptr<GeomObject> goobject(Factory::createObject(type));
    shared_ptr<ParamCurve> tmp_crv =
	std::dynamic_pointer_cast<ParamCurve, GeomObject>(goobject);
    ALWAYS_ERROR_IF(tmp_crv.get() == 0,
		    "Can not read this instance type");
    tmp_crv->read(is);
    curve_ = tmp_crv;
    
}


/// Write object to stream
//===========================================================================
void BoundedCurve::write(std::ostream& os) const
//===========================================================================
{
    int i;
    //os << setprecision(15);
    os << curve_->instanceType();
    os << endl;

    if (prefer_parameter_)
	os << 1 << " ";
    else 
	os << 0 << " ";
    os << startparam_ << " " << endparam_ << endl;

    int dim = start_pt_.dimension();
    os << dim << "  ";
    for (i=0; i<dim; ++i)
	os << start_pt_[i] << "  ";
    for (i=0; i<dim; ++i)
	os << end_pt_[i] << "  ";
    os << endl;

    curve_->write(os);

}


//===========================================================================
BoundingBox BoundedCurve::boundingBox() const
//===========================================================================
{
    shared_ptr<ParamCurve> sub_cv(curve_->subCurve(startparam_, endparam_));
    return sub_cv->boundingBox();
}
    

//===========================================================================
int BoundedCurve::dimension() const
//===========================================================================
{
    return curve_->dimension();
}


//===========================================================================
ClassType BoundedCurve::instanceType() const
//===========================================================================
{
    return classType();
}


//===========================================================================
ClassType BoundedCurve::classType()
//===========================================================================
{
    return Class_BoundedCurve;
}


//===========================================================================
BoundedCurve* BoundedCurve::clone() const
//===========================================================================
{
    BoundedCurve* bd_cv = new BoundedCurve(curve_, prefer_parameter_,
					   startparam_, endparam_,
					   start_pt_, end_pt_);

    return bd_cv;
}


//===========================================================================
void BoundedCurve::point(Point& pt, double tpar) const
//===========================================================================
{
    curve_->point(pt, tpar);
}
 

//===========================================================================
void BoundedCurve::point(std::vector<Point>& pts, 
			 double tpar,
			 int derivs,
			 bool from_right) const
//===========================================================================
{
    curve_->point(pts, tpar, derivs, from_right);
}

 
//===========================================================================
double BoundedCurve::startparam() const
//===========================================================================
{
    return startparam_;
}


//===========================================================================
double BoundedCurve::endparam() const
//===========================================================================
{
    return endparam_;
}

 
//===========================================================================
void BoundedCurve::reverseParameterDirection(bool switchparam)
//===========================================================================
{
    curve_->reverseParameterDirection(switchparam);
    std::swap(start_pt_, end_pt_);
}
    

//===========================================================================
void BoundedCurve::setParameterInterval(double t1, double t2)
//===========================================================================
{
    setParamBounds(t1, t2);
}


//===========================================================================
SplineCurve* BoundedCurve::geometryCurve()
//===========================================================================
{
    shared_ptr<ParamCurve> sub_cv(curve_->subCurve(startparam_, endparam_));

    return sub_cv->geometryCurve();
}


//===========================================================================
bool BoundedCurve::isDegenerate(double degenerate_epsilon)
//===========================================================================
{
    shared_ptr<ParamCurve> sub_cv(curve_->subCurve(startparam_, endparam_));
    return sub_cv->isDegenerate(degenerate_epsilon);
}

     
//===========================================================================
BoundedCurve* BoundedCurve::subCurve(double from_par, double to_par,
				     double fuzzy) const
//===========================================================================
{
    BoundedCurve* bd_cv = clone();
    bd_cv->setParamBounds(from_par, to_par);

    return bd_cv;
}


//===========================================================================
DirectionCone BoundedCurve::directionCone() const
//===========================================================================
{
    shared_ptr<ParamCurve> sub_cv(curve_->subCurve(startparam_, endparam_));
    return sub_cv->directionCone();
}
 

//===========================================================================
void BoundedCurve::appendCurve(ParamCurve* cv, bool reparam)
//===========================================================================
{
    MESSAGE("Not implemented!");
}


//===========================================================================
void BoundedCurve::appendCurve(ParamCurve* cv,
			       int continuity, double& dist, bool reparam)
//===========================================================================
{
    MESSAGE("Not implemented!");
}


//===========================================================================
void BoundedCurve::closestPoint(const Point& pt,
				double tmin,
				double tmax,
				double& clo_t,
				Point& clo_pt,
				double& clo_dist,
				double const *seed) const
//===========================================================================
{
    // Check and fix the parameter bounds
    if (tmin < startparam_) {
	tmin = startparam_;
	MESSAGE("tmin too small. Using startparam_.");
    }
    if (tmax > endparam_) {
	tmax = endparam_;
	MESSAGE("tmax too large. Using endparam_.");
    }

    shared_ptr<ParamCurve> sub_cv(curve_->subCurve(tmin, tmax));
    sub_cv->ParamCurve::closestPoint(pt, clo_t, clo_pt, clo_dist);
}


//===========================================================================
double BoundedCurve::length(double tol)
//===========================================================================
{
    return ParamCurve::length(tol, startparam_, endparam_);
}


//===========================================================================
void BoundedCurve::setParamBounds(double startpar, double endpar)
//===========================================================================
{
    if (startpar >= endpar)
	THROW("First parameter must be strictly less than second.");

    startparam_ = startpar;
    endparam_ = endpar;
}


} // namespace Go
