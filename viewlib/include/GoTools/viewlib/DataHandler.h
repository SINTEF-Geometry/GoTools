//===========================================================================
//                                                                           
// File: dataHandler.h                                                       
//                                                                           
// Created: Wed Nov 28 17:17:21 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: DataHandler.h,v 1.1 2007-04-17 12:25:30 sbr Exp $
//                                                                           
// Description: This class connects the four inheritance trees descending
//              from gvTesselator, gvPaintable, gvPropertySheet and
//              GeomObject.
//                                                                           
//===========================================================================

#ifndef _DATAHANDLER_H
#define _DATAHANDLER_H


#include "GoTools/viewlib/gvColor.h"
#include "GoTools/tesselator/Tesselator.h"

#include "GoTools/utils/config.h"


//class Tesselator;
class gvPaintable;
class gvPropertySheet;
namespace Go
{
    class GeomObject;
}

/** This class connects the four inheritance trees descending
 * from gvTesselator, gvPaintable, gvPropertySheet and GeomObject.
 */

class DataHandler
{
public:
    /// Default constructor.
    /// In inherited classes, it should Initialize factory
    /// (by registering classes).
    DataHandler();
    virtual ~DataHandler();

    virtual void create(shared_ptr<Go::GeomObject> obj,
			const gvColor& col, int id) = 0;

    shared_ptr<Go::Tesselator> tesselator()
    {
	return tesselator_;
    }
    shared_ptr<gvPaintable> paintable()
    {
	return paintable_;
    }
    shared_ptr<gvPropertySheet> propertySheet()
    {
	return property_sheet_;
    }

    void clear()
    {
        tesselator_.reset();
	paintable_.reset();
	property_sheet_.reset();
    }

protected:
    shared_ptr<Go::Tesselator> tesselator_;
    shared_ptr<gvPaintable> paintable_;
    shared_ptr<gvPropertySheet> property_sheet_;
};


#endif // _DATAHANDLER_H

