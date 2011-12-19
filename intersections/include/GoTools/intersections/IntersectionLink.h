//===========================================================================
//                                                                           
// File: IntersectionLink.h                                                  
//                                                                           
// Created: Fri Nov 26 12:18:05 2004                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: IntersectionLink.h,v 1.25 2007-11-01 14:31:36 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _INTERSECTIONLINK_H
#define _INTERSECTIONLINK_H


#include "GoTools/intersections/IntersectionPoint.h"
#include "GoTools/intersections/LinkType.h"
#include <memory>
#include <fstream> // remove after debugging
#include "GoTools/geometry/LineCloud.h" // remove after debugging


namespace Go {


/// This class represents a link between two intersection points.  It
/// is owned by both of the points in question, and can therefore refer
/// to them only by means of ordinary pointers (no shared_ptr stuff
/// here...)

class IntersectionLink {
public:
    /// Constructor
    IntersectionLink(IntersectionPoint* p1, 
		     IntersectionPoint* p2)
	: p1_(p1), p2_(p2), delimits_partial_coincidence_(false)
    { 
	ASSERT(p1_->getObj1() == p2->getObj1());
	ASSERT(p1_->getObj2() == p2->getObj2());
	for (int i = 0; i < 4; iso_[i++] = false);
    }

    /// This function will set the 'p1' and 'p2' to point to the
    /// IntersectionPoints participating in this link.
    void getIntersectionPoints(IntersectionPoint*& p1,
			       IntersectionPoint*& p2)
    {
	p1 = p1_;
	p2 = p2_;
    }

    /// This function will set the 'p1' and 'p2' to point to the
    /// IntersectionPoints participating in this link.
    void getIntersectionPoints(const IntersectionPoint*& p1,
			       const IntersectionPoint*& p2) const
    {
	p1 = p1_;
	p2 = p2_;
    }

    /// When one of the two points in the link is input, this function
    /// returns the other point.
    IntersectionPoint* getOtherPoint(const IntersectionPoint* p1)
    {
	return (p1 == p1_) ? p2_ : ((p1 == p2_) ? p1_ : 0);
    }

    /// Returns true if the link is attached to the input point.
    bool isAttachedTo(const IntersectionPoint* ip)
    {
	return (ip == p1_) || (ip == p2_);
    }

    /// Set all meta-information (all private info except
    /// IntersectionPoint-pointers) equal to the one in 'rhs'
    void copyMetaInformation(const IntersectionLink& rhs)
    {
	delimits_partial_coincidence_ = rhs.delimits_partial_coincidence_;
	for (int i = 0; i < 4; i++) {
	    iso_[i] = rhs.iso_[i];
	}
    }

    /// Find out whether this link is participating in defining a
    /// partial coincidence area
    bool delimitsPAC() const 
    {
	return delimits_partial_coincidence_;
    }

    void setPAC(bool value) 
    {
	delimits_partial_coincidence_ = value;
    }
    
    void setIsoparametricIn(int dir, bool iso) 
    {
	ASSERT(dir >= 0 && dir < (p1_->numParams1() + p1_->numParams2()));
	iso_[dir] = iso;
	
	// debug purposes
	if (iso_[0] && iso_[1]) {
	    //MESSAGE("Too many iso parametric directions");
// 	    ASSERT(false); // a link cannot be isoparametric in both
// 			   // directions
	} 
	if (iso_[2] && iso_[3]) {
	    //MESSAGE("Too many iso parametric directions");
//	    ASSERT(false);
	}
    }

    bool isIsoparametricIn(int dir)  const 
    {
	ASSERT(dir >= 0 && dir < numParams());
	return iso_[dir];
    }

    bool isIsoparametric() const 
    {
	int num_param = numParams();
	for (int i = 0; i < num_param; ++i) {
	    if (iso_[i]) {
		return true;
	    }
	}
	return false;
    }
    
    int numParams() const 
    {
	return p1_->numParams1() + p1_->numParams2();
    }

    /// Test for crossing between 'this' and 'link'.
    int crosses(const shared_ptr<IntersectionLink>& link) const
	// 
	// Return value : 2 - crossing links
	//                1 - close links
	//                0 - distant links
    {
	// Define two algebraic functions, a*x + b*y - c = 0, for
	// 'this' and 'link', respectively. Check if the two points on
	// one link lies on each side of the other, and vice versa.

	int close = 0;
	const IntersectionPoint *p1, *p2, *q1, *q2;
	getIntersectionPoints(p1, p2);
	link->getIntersectionPoints(q1, q2);

	double aeps = p1->getTolerance()->getEpsge();
	double tol = 100.0*aeps;

	// Only implemented for two surfaces
	if (p1->numParams1() != 2 || p1->numParams2() != 2) {
	    MESSAGE("IntersectionLink::crosses() is only implemented "
		    "for two surfaces.");
	    return 0;
	}

	// If the two links share a point, they don't cross and is not classified
	// as close
	if (p1 == q1 || p1 == q2 || p2 == q1 || p2 == q2)
	    return 0;

	// Don't consider links shorther than the tolerance
	if (p1->getPoint().dist(p2->getPoint()) < aeps ||
	    q1->getPoint().dist(q2->getPoint()) < aeps)
	    return 0;

	// Loop over first and second surface (i=0 or i=2)
	for (int i = 0; i < 4; i += 2) 
	{
	    // The points belonging to link two lies on either side
	    Point a = (i<2) ? p1->getPar1Point() : p1->getPar2Point();
	    Point b = (i<2) ? p2->getPar1Point() : p2->getPar2Point();
	    Point c = (i<2) ? q1->getPar1Point() : q1->getPar2Point();
	    Point d = (i<2) ? q2->getPar1Point() : q2->getPar2Point();
	    Point n1(-(b[1]-a[1]), b[0]-a[0]);
	    Point n2(-(d[1]-c[1]), d[0]-c[0]);
	    n1.normalize();
	    n2.normalize();
	    double s1 = (a - c)*n1;
	    double s2 = (a - d)*n1;
	    double s3 = (c - a)*n2;
	    double s4 = (c - b)*n2;
	    if (s1*s2 < 0.0 && s3*s4 < 0.0)
	    {
		return 2;
	    }
	    else if (fabs(s1) < tol || fabs(s2) < tol || fabs(s3) < tol ||
		     fabs(s4) < tol)
		close++;
	    // of the infinite line through the points belonging to link one
/* 	    // First link */
/* 	    double a0 = p1->getPar(i+1) - p2->getPar(i+1); */
/* 	    double b0 = p2->getPar(i) - p1->getPar(i); */
/* 	    double c0 = p1->getPar(i+1) * p2->getPar(i) */
/* 		- p1->getPar(i) * p2->getPar(i+1); */
/* 	    double q1val = a0*q1->getPar(i) + b0*q1->getPar(i+1) - c0; */
/* 	    double q2val = a0*q2->getPar(i) + b0*q2->getPar(i+1) - c0; */
/* 	    bool q_on_each_side = (q1val > 0.0 && q2val < 0.0) */
/* 		|| (q1val < 0.0 && q2val > 0.0); */
/* 	    if (!q_on_each_side) */
/* 		continue; */
/* 	    // Second link */
/* 	    double a1 = q1->getPar(i+1) - q2->getPar(i+1); */
/* 	    double b1 = q2->getPar(i) - q1->getPar(i); */
/* 	    double c1 = q1->getPar(i+1) * q2->getPar(i) */
/* 		- q1->getPar(i) * q2->getPar(i+1); */
/* 	    double p1val = a1*p1->getPar(i) + b1*p1->getPar(i+1) - c1; */
/* 	    double p2val = a1*p2->getPar(i) + b1*p2->getPar(i+1) - c1; */
/* 	    bool p_on_each_side = (p1val > 0.0 && p2val < 0.0) */
/* 		|| (p1val < 0.0 && p2val > 0.0); */
/* 	    if (p_on_each_side) { */
/* 		return true; */
/* 	    } */
	}

	return (close == 2) ? 1 : 0;
    }

    /// Returns the link type of the link
    /// \retval a LinkType enum
    const LinkType& linkType() const
    { return link_type_; }

    /// Returns the link type of the link
    /// \retval a LinkType enum
    LinkType& linkType()
    { return link_type_; }

    /// Writes info about the link to standard output
    void writeInfo() const
    {
	// Parameters
	p1_->writeParams(std::cout);
	std::cout << "  --  ";
	p2_->writeParams(std::cout);
	std::cout << std::endl << "\t";

	// Cosines. cos1 refers to the angle between the tangent in p1
	// and the vector from p1 to p2. Similarly with cos2.
	Point linkdir = p2_->getPoint() - p1_->getPoint();
	double len = linkdir.length();
	bool has_nonzero_link = (len != 0.0);
	if (has_nonzero_link) {
	    linkdir.normalize();
	    bool has_tangent;
	    has_tangent
		= (p1_->getSingularityType() == ORDINARY_POINT
		   || p1_->getSingularityType() == TANGENTIAL_POINT);
	    if (has_tangent) {
		Point tangent = p1_->getTangent();
		tangent.normalize();
		double cos = tangent * linkdir;
		std::cout << "cos1 = ";
		std::cout.width(std::cout.precision()+4);
		std::cout << cos << "  ";
	    }
	    else {
		std::cout << "cos1 = no tangent  ";
	    } 

	    has_tangent
		= (p2_->getSingularityType() == ORDINARY_POINT
		   || p2_->getSingularityType() == TANGENTIAL_POINT);
	    if (has_tangent) {
		Point tangent = p2_->getTangent();
		tangent.normalize();
		double cos = -tangent * linkdir;
		std::cout << "cos2 = ";
		std::cout.width(std::cout.precision()+4);
		std::cout << cos << "  ";
	    }
	    else {
		std::cout << "cos2 = no tangent  ";
	    }
	}
	else {
	    std::cout << "link length = 0.0                    ";
	}

	// Link type
	std::cout << "type = ";
	std::cout.width(2);
	std::cout << link_type_;
	
	return;
    }

    void dumpToStream(std::ostream& os) // debug reasons
    {
	double data[12];
	Point tmp = p1_->getPoint1();
	std::copy(tmp.begin(), tmp.end(), data);
	tmp = p2_->getPoint1();
	std::copy(tmp.begin(), tmp.end(), data+3);
	tmp = p1_->getPoint2();
	std::copy(tmp.begin(), tmp.end(), data+6);
	tmp = p2_->getPoint2();
	std::copy(tmp.begin(), tmp.end(), data+9);
	LineCloud lc(data,2);
	lc.writeStandardHeader(os);
	lc.write(os);
    }

private:
    // the IntersectionLink should not be able to change or delete the
    // pointers to the IntersectionPoints it refers to. Hence the
    // 'const' keyword.  (It will however be able to call non-const
    // member functions of p1 and p2, so we cannot declare them to
    // point to 'const'-objects).
    IntersectionPoint* const p1_;
    IntersectionPoint* const p2_;
    
    bool delimits_partial_coincidence_;
    bool iso_[4];  // confirmed isoparametric intersection curve in a
		   // parameter

    LinkType link_type_;

};


}; // end namespace Go


#endif // _INTERSECTIONLINK_H

