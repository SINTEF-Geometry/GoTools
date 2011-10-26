//===========================================================================
//                                                                           
// File: ObjectHeader.C                                                    
//                                                                           
// Created: Wed Nov  8 14:24:37 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: ObjectHeader.C,v 1.5 2005-09-19 09:24:42 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/Utils.h"

namespace Go
{

//===========================================================================
ObjectHeader::~ObjectHeader()
//===========================================================================
{
}   

//===========================================================================
void ObjectHeader::read (std::istream& is)
//===========================================================================
{
    // We should verify that the object is valid.
    bool is_good = is.good();
    if (!is_good) {
	THROW("Invalid object header!");
    }
    int dummy;
    is >> dummy;
    class_type_ = static_cast<ClassType>(dummy);
    is >> major_version_;
    is >> minor_version_;
    int auxsize;
    is >> auxsize;
    is_good = is.good();
    if (!is_good) {
	THROW("Invalid object header!");
    }
    auxillary_data_.resize(auxsize);
    for (int i = 0; i < auxsize; ++i) {
	is >> auxillary_data_[i];
    }
}   

//===========================================================================
void ObjectHeader::write (std::ostream& os) const
//===========================================================================
{
    os << class_type_ << ' ';
    os << major_version_ << ' ';
    os << minor_version_ << ' ';
    os << auxillary_data_.size();
    for (size_t i = 0; i < auxillary_data_.size(); ++i) {
	os << ' ' << auxillary_data_[i];
    }
    os << '\n' << std::endl;
    
}

} // namespace Go
