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

#ifndef _GVAPPLICATION_H
#define _GVAPPLICATION_H


#include <QWidget>
#include <QString>
//Added by qt3to4:
// #include <Q3PopupMenu>
#include "GoTools/viewlib/gvData.h"
#include "GoTools/geometry/LineCloud.h"
#include <QMainWindow>

#include "GoTools/utils/config.h"

class gvView;
// class Q3PopupMenu;
class Q3ButtonGroup;
class QMenuBar;

/** gvApplication:
 *  etc
 */

// moc not too happy about long identifiers.
// #ifdef USE_BOOST
typedef std::vector<shared_ptr<Go::GeomObject> > ObjContainer;
typedef std::vector<shared_ptr<gvColor> >  ColContainer;
// #else
// typedef std::vector<shared_ptr<Go::GeomObject> > ObjContainer;
// typedef std::vector<shared_ptr<gvColor> >  ColContainer;
// #endif

class gvApplication : public QWidget
{

Q_OBJECT

public:
    /// Constructor.
    gvApplication(std::auto_ptr<DataHandler> dh,
		  QWidget* parent = 0,
		  const char* name = 0,
		  Qt::WFlags f = 0);
    /// The destructor
    virtual ~gvApplication();

    virtual QSize sizeHint() const;


public slots:
    void open();
    void reload_last_opened_file();
    void save_selection_as();
    void close_document();
    void quit();

    void about();

    virtual void view_reset();
    void view_reset_visible();
    void view_wireframe();
    void view_axis();
    void view_cull();
    void view_specular();
    void view_orthographic();
    void toggle_blending_mode();
    void view_focus_point();
    void view_focus_point_cb(int x, int y);
    //    void change_resolution_dialog();
    void display_object_properties();
    void assign_texture();
    void set_curve_resolutions();
    virtual void set_surface_resolutions();
    void enable_objects();
    void disable_objects();
    void toggle_enable();

    void select_all();
    void select_none();
    void select_inverse();
    void select_all_surfaces();
    void select_all_curves();
    void select_all_visible();
    void toggle_selection_mode();
    void toggle_select_object(unsigned int name);
    void toggle_multiselect_object(unsigned int* names, int numnames);

    /// Put selected objects in a group (user required to name group).
    void group_selected();
    virtual void dismiss_selections();

    //    void change_resolution(int u, int v);  

    // For a cv or surface we plot the 3D control net. If surface is
    // trimmed we plot the net of the underlying sf.
    void show_control_nets();

    void set_random_color();

protected:
    void buildGUI();
    Q3ButtonGroup* createObjectToggleBox();

   // Selected objects are extracted from data_ and returned in vector.
   void getSelectedObjects(std::vector< shared_ptr< Go::GeomObject > >&
			   sel_objs,
			   std::vector< shared_ptr< Go::GeomObject > >&
			   not_sel_objs);

    gvData data_;
    gvView* view_;
    QString curr_file_type_;
    QString last_file_name_;
    QString app_name_; // Used when setting caption.

    QMenuBar* menu_;
//     QMenuBar* menu_;
    QMenu* view_menu_;
    QMenu* select_menu_;
    QMenu* group_menu_;
    QMenu* object_menu_;
//     Q3PopupMenu* view_menu_;
//     Q3PopupMenu* select_menu_;
//     Q3PopupMenu* group_menu_;
    shared_ptr <QWidget> actionForm;

protected slots:
    void changeCurveResolutions(int new_res); // Change resolution of
					      // all selected curves.
    virtual void
    changeSurfaceResolutions(int new_u_res,
			     int new_v_res); // Change resolution of
					     // all selected sfs.
  virtual void add_group(std::vector<int>& members, QString name);

    void add_objects(ObjContainer& new_objs,
		     ColContainer& new_colors);

private:

  // Expecting spline cvs or spline sfs (possibly trimmed sfs, in
  // which case we'll use underlying spline sf).
  shared_ptr<Go::LineCloud>
  getLineCloud(shared_ptr<Go::GeomObject>& obj);

 private slots:


};


#endif // _GVAPPLICATION_H

