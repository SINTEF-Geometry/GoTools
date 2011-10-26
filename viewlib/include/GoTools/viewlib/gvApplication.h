//===========================================================================
//                                                                           
// File: gvApplication.h                                                     
//                                                                           
// Created: Wed Jun 27 16:17:13 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: gvApplication.h,v 1.5 2008-03-07 13:03:51 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _GVAPPLICATION_H
#define _GVAPPLICATION_H


#include <QWidget>
#include <QString>
//Added by qt3to4:
// #include <Q3PopupMenu>
#include "GoTools/viewlib/gvData.h"
#include "GoTools/geometry/LineCloud.h"
#include <QMainWindow>

class gvView;
// class Q3PopupMenu;
class Q3ButtonGroup;
class QMenuBar;

/** gvApplication:
 *  etc
 */

// moc not too happy about long identifiers.
typedef std::vector<std::shared_ptr<Go::GeomObject> > ObjContainer;
typedef std::vector<std::shared_ptr<gvColor> >  ColContainer;

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

    void view_reset();
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

protected:
    void buildGUI();
    Q3ButtonGroup* createObjectToggleBox();

   // Selected objects are extracted from data_ and returned in vector.
   void getSelectedObjects(std::vector< std::shared_ptr< Go::GeomObject > >&
			   sel_objs,
			   std::vector< std::shared_ptr< Go::GeomObject > >&
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
//     Q3PopupMenu* view_menu_;
//     Q3PopupMenu* select_menu_;
//     Q3PopupMenu* group_menu_;
    std::shared_ptr <QWidget> actionForm;

protected slots:
    void changeCurveResolutions(int new_res); // Change resolution of
					      // all selected curves.
    virtual void
    changeSurfaceResolutions(int new_u_res,
			     int new_v_res); // Change resolution of
					     // all selected sfs.
  virtual void add_group(std::vector<int>& members, QString name);

private:

  // Expecting spline cvs or spline sfs (possibly trimmed sfs, in
  // which case we'll use underlying spline sf).
  std::shared_ptr<Go::LineCloud>
  getLineCloud(std::shared_ptr<Go::GeomObject>& obj);

 private slots:
 void add_objects(ObjContainer& new_objs,
		  ColContainer& new_colors);


};


#endif // _GVAPPLICATION_H

