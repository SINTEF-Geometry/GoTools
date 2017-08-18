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

#include "GoTools/viewlib/gvApplication.h"
#include "GoTools/viewlib/gvView.h"
#include "GoTools/viewlib/gvObjectList.h"
#include "GoTools/viewlib/gvPropertySheet.h"
#include "GoTools/viewlib/gvGroupPropertySheet.h"
#include "GoTools/viewlib/CurveResolutionSheet.h"
#include "GoTools/viewlib/SurfaceResolutionSheet.h"
#include "GoTools/tesselator/CurveTesselator.h"
#include "GoTools/viewlib/gvTexture.h"
#include "GoTools/tesselator/RectangularSurfaceTesselator.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/tesselator/ParametricSurfaceTesselator.h"
#include "GoTools/utils/errormacros.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/LineCloud.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/SplineCurve.h"
//#include "GoTools/viewlib/gvResolutionDialog.h"
//#include "GoTools/viewlib/DataHandler.h"

#include <QMenu>
#include <QMenuBar>
#include <QLayout>
#include <QApplication>
#include <QMessageBox>
#include <fstream>
#include <utility>
#include <QSlider>
#include <QLCDNumber>
#include <QString>
#include <vector>
// #include <QtGui/>
#include <QFileDialog>
#include <QtCore/QEvent>
#include <QSettings>
#include <assert.h>
#include <time.h>
#include <QScreen>

using std::vector;
using std::auto_ptr;
using std::ifstream;
using std::string;
using std::exception;
// using std::shared_ptr;
// using dynamic_pointer_cast;
using namespace Go;
class DataHandler;
//class ParametricSurfaceTesselator;



#define MAX_COLORS 12
int colors[MAX_COLORS][3] = {
  {255, 0, 0},
  {0, 255, 0},
  {0, 0, 255},
  {255, 255, 0},
  {255, 0, 255},
  {0, 255, 255},
  {128, 255, 0},
  {255, 128, 0},
  {128, 0, 255},
  {255, 0, 128},
  {0, 128, 255},
  {0, 255, 128},
};


string getFileName(const string& s) {

   char sep = '/';

#ifdef _WIN32
   sep = '\\';
#endif

   size_t i = s.rfind(sep, s.length());
   if (i != string::npos) {
      return(s.substr(i+1, s.length() - i));
   }

   return("");
}

//===========================================================================
gvApplication::gvApplication(auto_ptr<DataHandler> dh,
			     QWidget* parent,
			     const char* name,
			     Qt::WindowFlags f)
//===========================================================================
//   : QWidget(parent, name, f), data_(dh), curr_file_type_(".g2")
  : QWidget(parent, f), data_(dh), curr_file_type_("GO files (*.g2)"),
    last_file_name_("data")
{
//     app_name_ = QWidget::name();
    app_name_ = QWidget::accessibleName();
    // We're not interested in the full path.
    QChar sep('/');
//     int last_sep_occ = app_name_.findRev(sep, -1);
    int last_sep_occ = app_name_.lastIndexOf(sep, -1);
    if (last_sep_occ != -1)
      app_name_ = app_name_.remove(0, last_sep_occ + 1);

    buildGUI();
}


//===========================================================================
gvApplication::~gvApplication()
//===========================================================================
{
}


//===========================================================================
QSize gvApplication::sizeHint() const
//===========================================================================
{
    return QSize(500, 600);
}


//===========================================================================
void gvApplication::open()
//===========================================================================
{
//   curr_file_type_ = QString("IGES files (*.igs)");
//   string curr_file_type = ".g2";
  std::string text = string("GO files (*.g2);;") +
    string("SISL Surfs (*.srf);;") +
    string("SISL Curves (*.crv);;") +
    string("IGES files (*.igs);;") +
    string("All files (*)");
    try {
      // Get filename by dialog
	QSettings settings("SINTEF", "goview");
	QString last_open_dir = settings.value("last_open_dir").toString();
      QString s(QFileDialog::getOpenFileName
		(this,
		 tr("Open file"),
		 last_open_dir,
// 		 "Current file type (*" && ".g2" && ");;" &&
// 		 "SISL Surfs (*.srf);;" &&
// 		 "SISL Curves (*.crv);;" &&
// 		 "IGES files (*.igs);;" &&
// 		 "All files (*)")));//,
		 tr("GO files (*.g2);;"
		    "SISL Surfs (*.srf);;"
		    "SISL Curves (*.crv);;"
		    "IGES files (*.igs);;"
		    "All files (*)"),//));//,
		 &curr_file_type_));//,
// 		 QFileDialog::DontUseNativeDialog));
//       QString s(QFileDialog::getOpenFileName(0, "Current file type (*" +
// 				curr_file_type_ + ")\n"
// 				"Go files (*.g2)\n"
// 				"IGES files (*.igs)\n"
// 				"SISL surfs (*.srf)\n"
// 				"SISL curves (*.crv)\n"
// 				"all files (*)", 
// 				this));

	if (s.isNull())
	    return;

	ifstream infile(s.toLocal8Bit());
	ALWAYS_ERROR_IF(!infile, "No such file or file corrupt!");

	last_file_name_ = s;

	// Save the last open dir in settings
	last_open_dir = QFileInfo(last_file_name_).absolutePath();
	settings.setValue("last_open_dir", last_open_dir);

	int suffix_length = 4;//5;
	if (s.right(4) == ".igs") {
	  data_.readIges(infile);
	  curr_file_type_ = QString("IGES files (*.igs)");
	} else if (s.right(4) == ".srf") {
	  data_.readSislSrfs(infile);
	  curr_file_type_ = "SISL Surfs (*.srf)";
	} else if (s.right(4) == ".crv") {
	  data_.readSislCrvs(infile);
	  curr_file_type_ = "SISL Curves (*.crv)";
	} else if (s.right(3) == ".g2") {
	  data_.readGo(infile);
	  curr_file_type_ = "GO files (*.g2)";
	  --suffix_length;
	} else {
	  MESSAGE("Unexpected file suffix, returning.");
	  return;
	}

//  	cout << "Finished open()" << endl;
//  	curr_file_type_ = s.right(suffix_length);
// 	curr_file_type_ = "*"+s.right(suffix_length);

	// We're only interested in the file name (not the full path)
	QChar sep('/');
// 	int last_sep_occ = s.findRev(sep, -1);
	int last_sep_occ = s.lastIndexOf(sep, -1);
	QString local_file_name;
	if (last_sep_occ == -1)
	  local_file_name = s;
	else
	  local_file_name = s.remove(0, last_sep_occ + 1);
	QString caption_string = app_name_;
	caption_string.append(" (last opened file: " + local_file_name + ")");
 	setWindowTitle(caption_string);
//  	QEvent::WindowTitleChange(string);
    }
    catch (std::exception& e) {
        const char* what = e.what();
        MESSAGE("Continuing after failure: " << what);
    }
    catch (...) {
	MESSAGE("Continuing after failure.");
    }
}


//===========================================================================
void gvApplication::reload_last_opened_file()
//===========================================================================
{
    try {
	close_document();

	int suffix_length = 4;//5;
	ifstream infile(last_file_name_.toLocal8Bit());
	ALWAYS_ERROR_IF(!infile, "No such file or file corrupt!");
	if (last_file_name_.right(4) == ".igs") {
	    data_.readIges(infile);
	    curr_file_type_ = QString("IGES files (*.igs)");
	} else if (last_file_name_.right(4) == ".srf") {
	    data_.readSislSrfs(infile);
	    curr_file_type_ = "SISL Surfs (*.srf)";
	} else if (last_file_name_.right(4) == ".crv") {
	    data_.readSislCrvs(infile);
	    curr_file_type_ = "SISL Curves (*.crv)";
	} else if (last_file_name_.right(3) == ".g2") {
	    data_.readGo(infile);
	    curr_file_type_ = "GO files (*.g2)";
	    --suffix_length;
	} else {
	    MESSAGE("Unexpected file suffix, returning.");
	    return;
	}
    }
    catch (exception& e) {
        const char* what = e.what();
	MESSAGE("Continuing after failure: " << what);
    }
    catch (...) {
	MESSAGE("Continuing after failure.");
    }
}


//===========================================================================
void gvApplication::save_selection_as()
//===========================================================================
{
  try {
    // Get filename by dialog
    QString s(QFileDialog::getSaveFileName(this, "Go files (*.g2)\n"
					   "IGES files (*.igs)\n"
// 					   "SISL files (*.srf)\n"
// 					   "all files (*)", 
					   ));
    if (s.isNull())
	return;
  
    std::ofstream outfile(qPrintable(s));
    ALWAYS_ERROR_IF(!outfile, "File creation failed!");

    if (s.right(4) == ".igs")
      data_.writeSelectedIges(outfile);
    else if (s.right(3) == ".g2")    
      data_.writeSelectedGo(outfile);
    else
      {
	MESSAGE("Unexpected file suffix, returning.");
	return;
      }
  } catch (...) {
    MESSAGE("Continuing after failure.");
  }
}

//===========================================================================
void gvApplication::print_screenshot()
//===========================================================================
{
    string basename = getFileName(last_file_name_.toStdString());
    size_t lastindex = basename.find_last_of("."); 
    string rawname = basename.substr(0, lastindex); 
    string savename = "tmp/" + rawname + ".png";
    QImage qimage = view_->QGLWidget::grabFrameBuffer();//(0);
    QFile file(QString::fromStdString(savename));
    file.open(QIODevice::WriteOnly);
    qimage.save(&file, "PNG");

#if 0
    // Testing rotation and creation of video.
    const double fps = 15.0;
    const double sec = 1.0;
    int num_screenshots = 1 + ceiling(fps*sec);
    const double rot_rad = 2.0*M_PI/(double)num_screenshots;
    for (int ki = 0; ki < num_screenshots; ++ki)
    {
        
    }
#endif
}

//===========================================================================
void gvApplication::close_document()
//===========================================================================
{
    data_.clear();
    int nmb_objs = data_.numObjects(); // 0, I guess.
    data_.recreateDataStructure(nmb_objs);
//     QDialog::setCaption(app_name_);
}

//===========================================================================
void gvApplication::quit()
//===========================================================================
{
    QApplication::exit(0);
}


//===========================================================================
void gvApplication::about()
//===========================================================================
{
    // We print the date and time for the compilation of the app.
    QString date = __DATE__;
    QString time = __TIME__;
    QString build_date_info = "Build date: " + date + ", " + time;

    // We also store the build mode.
#ifdef QT_DEBUG
    QString build_mode = "debug";
#else
    QString build_mode = "release";
#endif
    QString build_mode_info = "Build mode: " + build_mode;

    QMessageBox msgBox;
    // The date and time variables are set during compile time, not
    // link time, hence not what we want (which is link date and time
    // for application).
//    msgBox.setText(build_date_info);
    msgBox.setInformativeText(build_mode_info);
    msgBox.setStandardButtons(QMessageBox::Ok);
    msgBox.setDefaultButton(QMessageBox::Ok);
    int ret = msgBox.exec();
}

//===========================================================================
void gvApplication::view_reset()
//===========================================================================
{
    view_->focusOnBox();
    view_->updateGL();
}

//===========================================================================
void gvApplication::view_reset_visible()
//===========================================================================
{
    try {
        view_->focusOnVisible();
    } catch (...) {
        MESSAGE("Failed calling focusOnVisible().");
    }
    view_->updateGL();
}


//===========================================================================
void gvApplication::view_wireframe()
//===========================================================================
{
    if (view_->wireframe()) {
// 	view_menu_->setItemChecked(2, false);
	view_->setWireframe(false);
    } else {
// 	view_menu_->setItemChecked(2, true);
	view_->setWireframe(true);
    }
}


//===========================================================================
void gvApplication::view_axis()
//===========================================================================
{
    if (view_->axis()) {
// 	view_menu_->setItemChecked(3, false);
	view_->setAxis(false);
    } else {
// 	view_menu_->setItemChecked(3, true);
	view_->setAxis(true);
    }
}

//===========================================================================
void gvApplication::view_cull()
//===========================================================================
{
    if (view_->backCull()) {
// 	view_menu_->setItemChecked(4, false);
	view_->setBackCull(false);
    } else {
// 	view_menu_->setItemChecked(4, true);
	view_->setBackCull(true);
    }
}

//===========================================================================
void gvApplication::view_specular()
//===========================================================================
{
    if (view_->specular()) {
// 	view_menu_->setItemChecked(5, false);
	view_->setSpecular(false);
    } else {
// 	view_menu_->setItemChecked(5, true);
	view_->setSpecular(true);
    }
}

//===========================================================================
void gvApplication::view_orthographic()
//===========================================================================
{
    if (view_->perspective()) {
// 	view_menu_->setItemChecked(6, true);
	view_->setPerspective(false);
    } else {
// 	view_menu_->setItemChecked(6, false);
	view_->setPerspective(true);
    }
}


//===========================================================================
void gvApplication::toggle_blending_mode()
//===========================================================================
{
    bool blending_mode = view_->blendingmode();
    view_->setBlendingmode(!blending_mode);
}

//===========================================================================

void gvApplication::view_focus_point()
//===========================================================================
{
  if (!view_->feedbackmode())
  {
     view_->setGetClickmode(true);
     connect(view_, SIGNAL(feedback(int, int)),
	     this, SLOT(view_focus_point_cb(int, int)));
  }
}

//===========================================================================
void gvApplication::view_focus_point_cb(int x, int y)
//===========================================================================
{
  disconnect(view_, SIGNAL(feedback(int, int)),
      this, SLOT(view_focus_point_cb(int, int)));
  view_->setCenter(x,y);
}
/*
//===========================================================================
void gvApplication::change_resolution_dialog()
//===========================================================================
{
    int ures = 10;
    int vres = 10;
    if (data_.numSurfaces() == 0) {
	QMessageBox::information (0, "Cannot use feature",
				  "You are not viewing any tesselated surfaces.", 0);
	return;
    }
    data_.getResSurface(0, ures, vres);
    gvResolutionDialog* resdia
	= new gvResolutionDialog(ures, vres, 2, 200, this);
    connect(resdia, SIGNAL(valuesChanged(int, int)),
	    this, SLOT(change_resolution(int, int)));
    resdia->show();
}
*/

//===========================================================================
void gvApplication::display_object_properties()
//===========================================================================
{
    // @@ Preliminary version. Works only with a single object selected.

    // Find the selected object's index
    int index = -1;
    for (int i = 0; i < data_.numObjects(); ++i) {
	if (data_.getSelectedStateObject(i)) {
	    if (index >= 0) {
		QMessageBox::warning( this, "Object properties",
				      "More than one object has been selected.",
				      QMessageBox::Ok, Qt::NoButton);
		return;
	    }
	    index = i;
	}
    }
    if (index == -1) {
	QMessageBox::warning( this, "Object properties",
			      "Please select an object.",
			      QMessageBox::Ok, Qt::NoButton);
	return;
    }


    // We have an object index. Let's get the property sheet
    gvPropertySheet* sh = data_.propertySheet(index).get();
    if (sh == 0) {
	QMessageBox::warning( this, "Object properties",
			      "No property sheet defined.",
			      QMessageBox::Ok, Qt::NoButton);
	return;
    }
    sh->createSheet(0, view_);
}

//===========================================================================
void gvApplication::assign_texture()
//===========================================================================
{
    // @@ Preliminary version. Works only with a single object selected.

    // Find the selected object's index
    int index = -1;
    for (int i = 0; i < data_.numObjects(); ++i) {
	if (data_.getSelectedStateObject(i)) {
	    if (index >= 0) {
		QMessageBox::warning( this, "Assign texture",
				      "More than one object has been selected.",
				      QMessageBox::Ok, Qt::NoButton);
		return;
	    }
	    index = i;
	}
    }
    if (index == -1) {
	QMessageBox::warning( this, "Assign texture",
			      "Please select an object.",
			      QMessageBox::Ok, Qt::NoButton);
	return;
    }

    // Open a dialog box to get a texture file
     try {
	// Get filename by dialog
	QString s(QFileDialog::getOpenFileName(this,"SGI rgb files (*.rgb)\n"
					       "all files (*)"
					       ));
	if (s.isNull())
	    return;
	shared_ptr<gvTexture> texture(new gvTexture(string(s.toLatin1())));
	texture->setEnvMode(envModulate);
	data_.setTexture(index, texture);
    }
    catch (exception& e) {
        const char* what = e.what();
	MESSAGE("Continuing after failure: " << what);
    }
   
}

//===========================================================================
void gvApplication::set_curve_resolutions()
//===========================================================================
{
    // We run through selected cvs and extract lowest resolution.
    int low = -1;
    int nmb_selected = 0;
    for (int ki = 0; ki < data_.numObjects(); ++ki) {
	if (data_.getSelectedStateObject(ki)) {
	    ++nmb_selected;
	    // Hmm ... Either a rectangular or parametric. Casting.
	    CurveTesselator* tess =
		dynamic_cast<CurveTesselator*>(data_.tesselator(ki).get());
	    int res = -1;
	    if (tess != 0) {
		tess->getRes(res);
	    }
	    if (res < 0) {
		MESSAGE("Changing curve resolutions: Unknown curve type.");
		// QMessageBox::warning( this, "Changing curve resolutions: ",
		// 		      "Unknown curve type.",
		// 		      QMessageBox::Ok, Qt::NoButton);
		res = 500;
	    }

	    if (low < 0 || res < low)
		low = res;
	}
    }

    if (nmb_selected == 0) {
	QMessageBox::warning( this, "Changing curve resolutions: ",
			      "No object has been selected.",
			      QMessageBox::Ok, Qt::NoButton);
	return;
    }

    CurveResolutionSheet* sh = new CurveResolutionSheet(low);
    sh->createSheet(this, view_);

    connect(sh, SIGNAL(return_value(int)),
	    this, SLOT(changeCurveResolutions(int)));
}

//===========================================================================
void gvApplication::set_surface_resolutions()
//===========================================================================
{
    // We run through selected sfs and extract lowest resolution.
    int low_u = -1;
    int low_v = -1;
    int nmb_selected = 0;
    for (int ki = 0; ki < data_.numObjects(); ++ki) {
	if (data_.getSelectedStateObject(ki)) {
	    ++nmb_selected;
	    // Hmm ... Either a rectangular or parametric. Casting.
	    RectangularSurfaceTesselator* tess =
		dynamic_cast<RectangularSurfaceTesselator*>(data_.tesselator(ki).get());
	    int ures = -1, vres = -1;
	    if (tess != 0) {
		tess->getRes(ures, vres);
	    } else {
		ParametricSurfaceTesselator* tess =
		    dynamic_cast<ParametricSurfaceTesselator*>(data_.tesselator(ki).get());
		if (tess != 0) {
		    tess->getRes(ures, vres);
		}
	    }

	    if (ures < 0 && vres < 0) {
		MESSAGE("Changing surface resolutions: Unknown surface type (or something else). Type: "
                        << data_.object(ki)->instanceType());
		// QMessageBox::warning( this, "Changing surface resolutions: ",
		// 		      "Unknown surface type.",
		// 		      QMessageBox::Ok, Qt::NoButton);
		vres = ures = 20;
	    }

	    if (low_u < 0 || ures < low_u)
		low_u = ures;
	    if (low_v < 0 || vres < low_v)
		low_v = vres;
	}
    }

    if (nmb_selected == 0) {
	QMessageBox::warning( this, "Changing surface resolutions: ",
			      "No object has been selected.",
			      QMessageBox::Ok, Qt::NoButton);
	return;
    }

    SurfaceResolutionSheet* sh = new SurfaceResolutionSheet(low_u, low_v);
    sh->createSheet(this, view_);

    connect(sh, SIGNAL(return_value(int, int)),
	    this, SLOT(changeSurfaceResolutions(int, int)));
}

//===========================================================================
void gvApplication::enable_objects()
//===========================================================================
{
  for (int i = 0; i < data_.numObjects(); ++i) {
    if (data_.getSelectedStateObject(i))
    {
      data_.setVisibleStateObject(i, true);
    }
  }
}

//===========================================================================
void gvApplication::disable_objects()
//===========================================================================
{
  for (int i = 0; i < data_.numObjects(); ++i) {
    if (data_.getSelectedStateObject(i))
    {
      data_.setVisibleStateObject(i, false);
    }
  }
}

//===========================================================================
void gvApplication::toggle_enable()
//===========================================================================
{
  for (int i = 0; i < data_.numObjects(); ++i) {
    if (data_.getSelectedStateObject(i))
    {
      bool state=data_.getVisibleStateObject(i);
      data_.setVisibleStateObject(i, !state);
    }
  }
}

//===========================================================================
void gvApplication::select_all()
//===========================================================================
{
    data_.disableUpdates();
    for (int i = 0; i < data_.numObjects(); ++i) {
	data_.setSelectedStateObject(i, true);
    }
    data_.enableUpdates();
    data_.updateObservers();
}

//===========================================================================
void gvApplication::select_all_surfaces()
//===========================================================================
{
    data_.disableUpdates();
    for (int i = 0; i < data_.numObjects(); ++i) {
	shared_ptr<Go::ParamSurface> surf =
	    dynamic_pointer_cast<Go::ParamSurface, Go::GeomObject>(data_.object(i));
	if (surf.get() != 0)
	    data_.setSelectedStateObject(i, true);
    }
    data_.enableUpdates();
    data_.updateObservers();
}

//===========================================================================
void gvApplication::select_all_curves()
//===========================================================================
{
    data_.disableUpdates();
    for (int i = 0; i < data_.numObjects(); ++i) {
	shared_ptr<Go::ParamCurve> crv =
	    dynamic_pointer_cast<Go::ParamCurve, Go::GeomObject>(data_.object(i));
	if (crv.get() != 0)
	    data_.setSelectedStateObject(i, true);
    }
    data_.enableUpdates();
    data_.updateObservers();
}

//===========================================================================
void gvApplication::select_all_visible()
//===========================================================================
{
  data_.disableUpdates();
  for (int i = 0; i < data_.numObjects(); ++i) {
    data_.setSelectedStateObject(i, data_.getVisibleStateObject(i));
  }
  data_.enableUpdates();
  data_.updateObservers();
}

//===========================================================================
void gvApplication::select_none()
//===========================================================================
{
    data_.disableUpdates();
    for (int i = 0; i < data_.numObjects(); ++i) {
	data_.setSelectedStateObject(i, false);
    }
    data_.enableUpdates();
    data_.updateObservers();
}

//===========================================================================
void gvApplication::select_inverse()
//===========================================================================
{
    data_.disableUpdates();
    for (int i = 0; i < data_.numObjects(); ++i) {
	bool current_state = data_.getSelectedStateObject(i);
	data_.setSelectedStateObject(i, !current_state);
    }
    data_.enableUpdates();
    data_.updateObservers();
}


//===========================================================================
void gvApplication::toggle_selection_mode()
//===========================================================================
{
    if (view_->selectionmode()) {
	view_->setSelectionmode(false);
// 	if (!view_->selectionmode())
// 	    select_menu_->setItemChecked(7, false);
    } else {
// 	select_menu_->setItemChecked(7, true);
	view_->setSelectionmode(true);
    }
}

//===========================================================================
void gvApplication::toggle_select_object(unsigned int i)
//===========================================================================
{
    data_.disableUpdates();
    bool current_state = data_.getSelectedStateObject(i);
    data_.setSelectedStateObject(i, !current_state);
    data_.enableUpdates();
    data_.updateObservers();
}

//===========================================================================
void gvApplication::toggle_multiselect_object(unsigned int* names, int numnames)
//===========================================================================
{
    data_.disableUpdates();
    for (int i = 0; i < numnames; ++i) {
	bool current_state = data_.getSelectedStateObject(names[i]);
	data_.setSelectedStateObject(names[i], !current_state);
    }
    data_.enableUpdates();
    data_.updateObservers();
}

//===========================================================================
void gvApplication::group_selected()
//===========================================================================
{
    vector<int> selected;
    for (int i = 0; i < data_.numObjects(); ++i) {
	if (data_.getSelectedStateObject(i)) {
	    selected.push_back(i);
	}
    }

    if (selected.size() > 0) {
	int nmb_groups = data_.nmbGroups() + 1;
	char buffer[8]; // Allows 999999 # of groups, should be enough.
	sprintf(buffer, "gp%d", nmb_groups);
	QString def_name(buffer);
	gvGroupPropertySheet* group_sheet = new gvGroupPropertySheet(selected, def_name);
	group_sheet->createSheet(0, view_);

	// If user clicks OK, group is created and pushed back.
	QObject::connect(group_sheet,
			 SIGNAL(value_changed(std::vector<int>&, QString)),
			 this, SLOT(add_group(std::vector<int>&, QString)));
    }
}

//===========================================================================
void gvApplication::dismiss_selections()
//===========================================================================
{
    data_.clearGroup();
}

//===========================================================================
void gvApplication::show_control_nets()
//===========================================================================
{
    // We extract all objects currently selected.
    vector<shared_ptr<Go::GeomObject> > sel_objs, not_sel_objs;
    getSelectedObjects(sel_objs, not_sel_objs);
    // We then extract those which are of the required types.
    vector<shared_ptr<GeomObject> > sel_geoms;
    for (size_t ki = 0; ki < sel_objs.size(); ++ki) {
	if (sel_objs[ki]->instanceType() == Class_SplineSurface) {
	    sel_geoms.push_back(sel_objs[ki]);
	} else if (sel_objs[ki]->instanceType() == Class_BoundedSurface) {
	    shared_ptr<BoundedSurface> bd_sf =
		dynamic_pointer_cast<BoundedSurface, GeomObject>(sel_objs[ki]);
	    if (bd_sf->underlyingSurface()->instanceType() ==
		Class_SplineSurface) {
		sel_geoms.push_back(bd_sf->underlyingSurface());
	    } else {
		not_sel_objs.push_back(sel_objs[ki]);
	    }
	} else if (sel_objs[ki]->instanceType() == Class_SplineCurve) {
	    sel_geoms.push_back(sel_objs[ki]);
	} else {
	    not_sel_objs.push_back(sel_objs[ki]);
	}
    }

    vector<shared_ptr<GeomObject> > new_objs;
    // We then proceed to update view with the control nets.
    for (size_t ki = 0; ki < sel_geoms.size(); ++ki) {
	// We create LineCloud
	shared_ptr<LineCloud> lc = getLineCloud(sel_geoms[ki]);
	new_objs.push_back(lc);
	if (sel_geoms[ki]->instanceType() == Class_SplineCurve) {
	    vector<double> pts;
	    int dim = lc->dimension();
            if (dim != 3) {
                MESSAGE("Expecting LineCloud (from SplineCurve) to be 3D.");
                continue;
            }
	    pts.insert(pts.end(), lc->rawData(), lc->rawData() + dim);
	    for (int kj = 0; kj < lc->numLines(); ++kj) {
		pts.insert(pts.end(), lc->rawData() + kj*dim*2 + dim,
			   lc->rawData() + kj*dim*2 + 2*dim);
	    }
            // If spline curve is 2D, the linecloud was raised to 3D.
	    shared_ptr<PointCloud3D> pt_cl
		(new PointCloud3D(pts.begin(), (int)pts.size()/dim));
	    new_objs.push_back(pt_cl);
	}
    }

    // Finally we send the new objs to the tesselator.
    vector<shared_ptr<gvColor> > new_cols(new_objs.size());
    add_objects(new_objs, new_cols);
}

//===========================================================================
void gvApplication::set_random_color()
//===========================================================================
{
    srand ( time(NULL) );
    int col_id = rand() % MAX_COLORS;
    for (int ki = 0; ki < data_.numObjects(); ++ki)
    {
	if (data_.getSelectedStateObject(ki))
	{
	    col_id = col_id%MAX_COLORS;
	    shared_ptr<gvColor> col = data_.color(ki);
	    if (col.get() == NULL)
	    {
		col = shared_ptr<gvColor>(new gvColor());

	    }
	    for (int kj = 0; kj < 3; ++kj)
		col->rgba[kj] = (float)colors[col_id][kj]/(float)255;
	    ++col_id;
	}
    }

    int num_objs = data_.numObjects();
    try {
        data_.recreateDataStructure(num_objs); // @@sbr A bit overkill perhaps ...
    } catch (...) {
        MESSAGE("Failed recreating data structure!");
    }
}

/*
//===========================================================================
void gvApplication::change_resolution(int u, int v)
//===========================================================================
{
    data_.disableUpdates();
    for (int i = 0; i < data_.numSurfaces(); ++i) {
	data_.setResSurface(i, u, v);
    }
    data_.enableUpdates();
    data_.updateObservers();
}
*/
//===========================================================================
void gvApplication::buildGUI()
//===========================================================================
{

  
    //=======================================================================
    //--------------- on top: a menu ----------------------------------------
    //=======================================================================

     menu_ = new QMenuBar(this);
//      menu_ = QLayout::menuBar();
    {
	//---------------------------------------------------------------------
	//------------ first menu item: File ----------------------------------
	//---------------------------------------------------------------------

	QMenu* file_menu = new QMenu( "&File" );
 	menu_->addMenu(file_menu);
	file_menu->addAction( "Open", this, SLOT(open()), Qt::CTRL+Qt::Key_O );
	file_menu->addAction( "Reload (last opened file)", this,
			       SLOT(reload_last_opened_file()),
			       Qt::CTRL+Qt::Key_R );
	file_menu->addAction( "Save selection as...", this,
			       SLOT(save_selection_as()),
			       Qt::CTRL+Qt::SHIFT+Qt::Key_S );
	file_menu->addAction( "Print screenshot (to tmp/<filename>.png)", this,
			       SLOT(print_screenshot()),
			       Qt::CTRL+Qt::SHIFT+Qt::Key_P );
	file_menu->addAction( "Close", this,
			       SLOT(close_document()),
			       Qt::CTRL+Qt::Key_W );
	file_menu->addSeparator();
	file_menu->addAction( "Exit", this, SLOT(quit()), Qt::Key_Q );

	//---------------------------------------------------------------------
	//------------ second menu item: View ------------------------------
	//---------------------------------------------------------------------
	
	view_menu_ = new QMenu( "&View" );
// 	view_menu_->setCheckable(true);
 	menu_->addMenu(view_menu_);
	view_menu_->addAction("Reset", this, 
			      SLOT(view_reset()), Qt::Key_R);//, 1);
	view_menu_->addAction("Wireframe mode", this, 
			      SLOT(view_wireframe()), Qt::Key_W);//, 2);
	view_menu_->addAction("Axis on/off", this, 
			      SLOT(view_axis()), Qt::Key_X);//, 3);
	view_menu_->addAction("Cullface mode", this, 
			      SLOT(view_cull()), Qt::Key_B);//, 4);
	view_menu_->addAction("Specular highlight mode", this, 
			      SLOT(view_specular()), Qt::Key_H);//, 5);
	view_menu_->addAction("Orthographic projection mode", this, 
			      SLOT(view_orthographic()), Qt::Key_O);//, 6);
	view_menu_->addAction("Focus on enabled", this, 
			      SLOT(view_reset_visible()), 
			      Qt::SHIFT+Qt::Key_R);//, 7);
	view_menu_->addAction("Set focus point", this, 
			      SLOT(view_focus_point()), 
			      Qt::SHIFT+Qt::Key_C);//, 9);
	view_menu_->addAction("Toggle blending mode", this, 
			      SLOT(toggle_blending_mode()));//, 0, 10);
//  	view_menu_->insertItem("Change resolution...", this, 
//  			      SLOT(change_resolution_dialog()), Qt::CTRL+Qt::Key_R, XX);
	

	//---------------------------------------------------------------------
	//------------ third menu item: Select --------------------------------
	//---------------------------------------------------------------------

  	select_menu_ = new QMenu( "&Select" );
 	menu_->addMenu(select_menu_);
	select_menu_->addAction("Select all", this, 
				SLOT(select_all()), Qt::CTRL+Qt::Key_A);//, 1);
	select_menu_->addAction("Select none", this, 
				SLOT(select_none()), Qt::CTRL+Qt::Key_N);//, 2);
	select_menu_->addAction("Inverse selection", this, 
				SLOT(select_inverse()), Qt::CTRL+Qt::Key_I);//, 3);
	select_menu_->addAction("Select all surfaces", this, 
				SLOT(select_all_surfaces()));//, 0, 4);
	select_menu_->addAction("Select all curves", this, 
				SLOT(select_all_curves()));//, 0, 5);
	select_menu_->addAction("Select all enabled", this, 
				SLOT(select_all_visible()));//, 0, 6);
	select_menu_->addAction("Selection mode", this, 
				SLOT(toggle_selection_mode()), Qt::Key_S);//, 7);

	//---------------------------------------------------------------------
	//------------ fourth menu item: Group -------------------------------
	//---------------------------------------------------------------------

	group_menu_ = new QMenu( "&Group" );
 	menu_->addMenu(group_menu_);
	group_menu_->addAction( "Group selected...", this,
				 SLOT(group_selected()),
				Qt::CTRL+Qt::Key_G);//, 1 );
	group_menu_->addAction( "Dismiss selections.", this, 
				SLOT(dismiss_selections()));//,
// 				 0, 2 );


	//---------------------------------------------------------------------
	//------------ fifth menu item: Object -------------------------------
	//---------------------------------------------------------------------

  	object_menu_ = new QMenu( "&Object" );
 	menu_->addMenu(object_menu_);
	object_menu_->addAction("Properties...", this, 
				SLOT(display_object_properties()),
			       Qt::CTRL+Qt::Key_P);//, 1);
	object_menu_->addAction("Assign texture...", this, 
				SLOT(assign_texture()),
			       Qt::Key_T);//, 2);
	object_menu_->addAction("Curve Resolutions...", this, 
			       SLOT(set_curve_resolutions()));//,
// 				0, 3);
	object_menu_->addAction("Surface Resolutions...", this, 
			       SLOT(set_surface_resolutions()));//,
// 				0, 4);
	object_menu_->addAction("Show control nets...", this, 
			       SLOT(show_control_nets()));//,
// 				 0, 5 );
	object_menu_->addAction("Set random color...", this, 
			       SLOT(set_random_color()));//,
// 				 0, 5 );
	object_menu_->addAction("Toggle enable", this, 
				SLOT(toggle_enable()),
			       Qt::CTRL+Qt::Key_T);//, 6);
	object_menu_->addAction("Enable objects", this, 
				SLOT(enable_objects()),
			       Qt::CTRL+Qt::Key_E);//, 7);
	object_menu_->addAction("Disable objects", this, 
				SLOT(disable_objects()),
			       Qt::CTRL+Qt::Key_D);//, 8);

	// object_menu_->addAction("Translate to origin", this, 
	// 			SLOT(translate_to_origin()));



	//---------------------------------------------------------------------
	//------------ last menu item: Help -----------------------------------
	//---------------------------------------------------------------------

	QMenu* help = new QMenu( "&Help" );
	menu_->addSeparator();
	menu_->addMenu(help);
	help->addAction( "About...", this, SLOT(about()) );
    }

    //=======================================================================
    //--------------- under the menubar: a gvView widget --------------------
    //=======================================================================

    view_ = new gvView(data_, this);
    view_->setMinimumSize(400, 400);
    // In order to handle mouse click events correctly
    setFocusProxy(view_);
    connect(view_, SIGNAL(objectPicked(unsigned int)),
	    this, SLOT(toggle_select_object(unsigned int)));
    connect(view_, SIGNAL(objectsPicked(unsigned int*, int)),
	    this, SLOT(toggle_multiselect_object(unsigned int*, int)));

    view_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

    //=======================================================================
    //--------------- to the right: object buttons  --------------------
    //=======================================================================
    gvObjectList* ol = new gvObjectList(data_, this);

    // Set up a vertical arrangement
    QVBoxLayout* vlayout = new QVBoxLayout(this);
    vlayout->setMenuBar(menu_);
    QHBoxLayout* hlayout = new QHBoxLayout();
//     hlayout->addItem(vlayout);
    hlayout->addWidget(view_, 1);
    hlayout->addWidget(ol, 0);
    vlayout->addItem(hlayout);
    //    vlayout->activate();
}

//===========================================================================
void gvApplication::changeCurveResolutions(int new_res)
//===========================================================================
{
    for (int i = 0; i < data_.numObjects(); ++i)
	if (data_.getSelectedStateObject(i)) {
	    CurveTesselator* tess =
	      dynamic_cast<CurveTesselator*>(data_.tesselator(i).get());
	    if (tess != 0)
		tess->changeRes(new_res);
	}
}

//===========================================================================
void gvApplication::changeSurfaceResolutions(int new_u_res, int new_v_res)
//===========================================================================
{
    for (int i = 0; i < data_.numObjects(); ++i)
	if (data_.getSelectedStateObject(i)) {
	    RectangularSurfaceTesselator* tess =
		dynamic_cast<RectangularSurfaceTesselator*>(data_.tesselator(i).get());
	    int u_res = -1, v_res = -1;
	    if (tess != 0) {
		tess->getRes(u_res, v_res);
		if ((u_res != new_u_res) || (v_res != new_v_res)) {
		    tess->changeRes(new_u_res, new_v_res);
		}
	    } else {
		ParametricSurfaceTesselator* tess =
		    dynamic_cast<ParametricSurfaceTesselator*>(data_.tesselator(i).get());
		if (tess != 0) {
		    tess->getRes(u_res, v_res);
		    if ((u_res != new_u_res) || (v_res != new_v_res)) {
			tess->changeRes(new_u_res, new_v_res);
		    }
		}
	    }
	}
}

//===========================================================================
void gvApplication::add_group(vector<int>& members, QString name)
//===========================================================================
{
    data_.addGroup(gvGroup(members, name));
}

//===========================================================================
void
gvApplication::getSelectedObjects(vector<shared_ptr<GeomObject> >& sel_objs,
				  vector<shared_ptr<GeomObject> >& not_sel_objs)
//===========================================================================
{
    sel_objs.clear();
    not_sel_objs.clear();
    int i;
    for (i = 0; i < data_.numObjects(); ++i) {
	if (data_.getSelectedStateObject(i))
	    sel_objs.push_back(data_.object(i));
	else
	    not_sel_objs.push_back(data_.object(i));
    }
}

//===========================================================================
shared_ptr<LineCloud> gvApplication::getLineCloud(shared_ptr<GeomObject>& obj)
//===========================================================================
{
    int dim = obj->dimension();
    ASSERT(dim == 2 || dim == 3);
    int ki;
    vector<double> lines;
    if (obj->instanceType() == Class_SplineSurface ||
	obj->instanceType() == Class_BoundedSurface) {
	shared_ptr<SplineSurface> spline_sf =
	    (obj->instanceType() == Class_SplineSurface) ?
	    dynamic_pointer_cast<SplineSurface, GeomObject>(obj) :
	    dynamic_pointer_cast<SplineSurface, ParamSurface>
	    ((dynamic_pointer_cast<BoundedSurface, GeomObject>(obj))->
	     underlyingSurface());
	if (spline_sf.get() != 0) {
	    int kj;
	    int in1 = spline_sf->numCoefs_u();
	    int in2 = spline_sf->numCoefs_v();
	    vector<double>::iterator iter = spline_sf->coefs_begin();
	    for (ki = 0; ki < in1; ++ki) {
		for (kj = 0; kj < in2; ++kj) {
		    if (ki < in1 - 1) {
			lines.insert(lines.end(), iter + (kj*in1 + ki)*dim,
				     iter + (kj*in1 + ki + 1)*dim);
			if (dim == 2)
			  lines.push_back(0.0);
			lines.insert(lines.end(), iter + (kj*in1 + ki + 1)*dim,
				     iter + (kj*in1 + ki + 2)*dim);
			if (dim == 2)
			  lines.push_back(0.0);
		    }
		    if (kj < in2 - 1) {
			lines.insert(lines.end(), iter + (kj*in1 + ki)*dim,
				     iter + (kj*in1 + ki + 1)*dim);
			if (dim == 2)
			  lines.push_back(0.0);
			lines.insert(lines.end(),
				     iter + ((kj + 1)*in1 + ki)*dim,
				     iter + ((kj + 1)*in1 + ki + 1)*dim);
			if (dim == 2)
			  lines.push_back(0.0);
		    }
		}
	    }
	}
    } else if (obj->instanceType() == Class_SplineCurve) {
	shared_ptr<SplineCurve> spline_cv =
	    dynamic_pointer_cast<SplineCurve, GeomObject>(obj);
	int in = spline_cv->numCoefs();
	vector<double>::iterator iter = spline_cv->coefs_begin();
	lines.insert(lines.end(), iter, iter + dim);
        if (dim == 2)
            lines.push_back(0.0);
	for (ki = 0; ki < in - 1; ++ki) {
	    lines.insert(lines.end(), iter + ki*dim, iter + (ki + 1)*dim);
	    if (dim == 2)
	      lines.push_back(0.0);
	    lines.insert(lines.end(), iter + ki*dim, iter + (ki + 1)*dim);
	    if (dim == 2)
	      lines.push_back(0.0);
	}
	lines.insert(lines.end(), iter + ki*dim, iter + (ki + 1)*dim);
	if (dim == 2)
	  lines.push_back(0.0);
    }

    int nmb_lines = (int)(lines.size())/(2*dim);
    shared_ptr<LineCloud> line_cloud(new LineCloud(lines.begin(), nmb_lines));

    return line_cloud;
}


//===========================================================================
void gvApplication::add_objects(vector<shared_ptr<GeomObject> >& new_objs,
				vector<shared_ptr<gvColor> >& new_colors)
//===========================================================================
{
    dismiss_selections();
    data_.readGo(new_objs, new_colors); // Objects are appended to existing.
}
