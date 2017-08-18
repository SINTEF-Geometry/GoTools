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

#define QT_CLEAN_NAMESPACE

#include "GoTools/viewlib/gvView.h"
#include "GoTools/viewlib/gvPainter.h"
#include "GoTools/viewlib/gvData.h"
#include "GoTools/viewlib/gvUtilities.h"
#include "GoTools/geometry/ParamSurface.h"
#include <QPainter>
#include <QImage>
#include <QtGui/QMouseEvent>
#include <QScreen>

//#include "GoTools/viewlib/PBuffer.h"

using namespace std;

//===========================================================================
gvView::gvView(gvData& data,
	       QWidget* parent, const char* name,
	       const QGLWidget * shareWidget, Qt::WindowFlags f)
//===========================================================================
//     : QGLWidget(QGLFormat(AlphaChannel), parent, name, shareWidget, f),
    : QGLWidget(QGLFormat(QGL::AlphaChannel), parent, shareWidget, f),
      gl_initialized_(false),
      no_data_(true),
      selection_mode_(false),
      selecting_(false),
      feedback_mode_(false),
      get_click_mode_(false),
      mouse_is_active_(false),
//       keyboard_is_active_(false),
      blending_mode_(false),
      data_(data),
      camera_(Go::Vector3D(0.0, 0.0, 0.0), 0, 0, width(), height()),
      lights_camera_(Go::Vector3D(0.0, 0.0, 0.0), 0, 0, width(), height()),
      base_axis_size_(5.0),
      mouse_button_(Qt::NoButton),
//       keyboard_modifier_(Qt::NoModifier),
//       keyboard_event_(NULL),
      wireframe_(false),
      axis_(false),
      backcull_(false),
      specular_(false),
      coarseTex_(0),
      fineTex_(0),
      painter_(0),
      focus_on_origin_(false)
{
    data.registerObserver(this);
    no_data_ = (data.numObjects() == 0);
//     adjustSize();
//     resize(800, 800);
    setAutoFillBackground(true);
}

//===========================================================================
gvView::gvView(const QGLFormat &format, gvData& data,
	       QWidget* parent, const char* name,
	       const QGLWidget * shareWidget, Qt::WindowFlags f)
//===========================================================================
//     : QGLWidget(format, parent, name, shareWidget, f),
    : QGLWidget(format, parent, shareWidget, f),
      gl_initialized_(false),
      no_data_(true),
      selection_mode_(false),
      selecting_(false),
      feedback_mode_(false),
      get_click_mode_(false),
      mouse_is_active_(false),
//       keyboard_is_active_(false),
      blending_mode_(false),
      data_(data),
      camera_(Go::Vector3D(0.0, 0.0, 0.0), 0, 0, width(), height()),
      lights_camera_(Go::Vector3D(0.0, 0.0, 0.0), 0, 0, width(), height()),
      base_axis_size_(5.0),
      mouse_button_(Qt::NoButton),
//       keyboard_event_(NULL),
      wireframe_(false),
      axis_(false),
      backcull_(false),
      specular_(false),
      coarseTex_(0),
      fineTex_(0),
      painter_(0),
      focus_on_origin_(false)
{
    data.registerObserver(this);
    no_data_ = (data.numObjects() == 0);
    adjustSize();
//     resize(800, 800);
    setAutoFillBackground(true);
}

//===========================================================================
gvView::~gvView()
//===========================================================================
{
}

//===========================================================================
void gvView::initializeGL()
//===========================================================================
{
    gl_initialized_ = true;
    // Setting the background color.
//     glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // black
//     glClearColor(0.95f, 0.95f, 0.95f, 1.0f); // soft white
    // glClearColor(0.7f, 0.7f, 0.7f, 1.0f); // light gray
    glClearColor(0.5f, 0.5f, 0.5f, 1.0f); // dark gray
 
    glShadeModel( GL_SMOOTH );

    glDepthFunc(GL_LESS);
    glEnable(GL_DEPTH_TEST);

//      GLfloat white[] = { 1.0, 1.0, 1.0, 1.0 };
//      GLfloat gray[] = { 0.3, 0.3, 0.3, 1.0 };
//      GLfloat red[] = { 1.0, 0.0, 0.0, 1.0 };
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHTING);

    glLineWidth(1.0); // @@@ var 2 081209
    glPointSize(2.0);

    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1); // Default value = 0.
    // @bjornc: We get a better performance with two lights than with
    // GL_LIGHT_MODEL_TWO_SIDE. WHY?  
    //  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

    focusOnBox();

    // We need to trap the keyboard commands.
    setFocusPolicy( Qt::StrongFocus );
    setFocus();
}

//===========================================================================
void gvView::resizeGL(int w, int h)
//===========================================================================
{
    camera_.setViewPort(0, 0, w, h);
    updateGL();

//     int side = qMin(w, h);
//     glViewport((w - side) / 2, (h - side) / 2, side, side);
}


//===========================================================================
void gvView::paintGL()
//===========================================================================
{
//     if (painter_) delete painter_;
//     painter_ = new QPainter(this);

    GLfloat light_position0[] = { 0, 0, 1.0, 0 };
//     GLfloat Qt::white[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat white[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat dark_gray[] = { 0.1f, 0.1f, 0.1f, 1.0 };
    GLfloat black[] = { 0.0, 0.0, 0.0, 1.0 };
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    // Color for the axis
    if (backcull_)
      glEnable(GL_CULL_FACE);
    else
      glDisable(GL_CULL_FACE);
    float purple[] = { 0.7f, 0.0, 0.9f, 1.0};
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, purple);
    camera_.use();
    glPushMatrix();
    lights_camera_.useModelView();
    glLightfv(GL_LIGHT0, GL_POSITION, light_position0);
    glPopMatrix();
    glLightfv(GL_LIGHT0, GL_AMBIENT, dark_gray);
    if (specular_) {
	glLightfv(GL_LIGHT0, GL_SPECULAR, white);
    } else {
	glLightfv(GL_LIGHT0, GL_SPECULAR, black);
    }

    data_.painter().drawScene(this);
//     static int p1=0;
//     if (p1)
//     {
//        int mousex, mousey;
//        shared_ptr<const Go::ParamSurface> obj;
//        double tex_u; double tex_v;

//        getObjAndParam(10, 10, obj, tex_u, tex_v);
//     }
    drawOverlay();

}



//===========================================================================
void gvView::saveSnapshot(int w, int h, const QString &filename)
//===========================================================================
{
#if 1
    {
        QImage qimage = QGLWidget::grabFrameBuffer();//(0);
        QFile file(filename);
        file.open(QIODevice::WriteOnly);
        qimage.save(&file, "PNG");
    }
#else
  // temporarly removed.
   bool withAlpha=true;
   gvCamera oldcam=camera_;

   camera_.setViewPort(0, 0, w, h);

   PBuffer pb("rgb alpha depth");
   pb.Initialize(w, h, false, true);   
   pb.Activate();

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glShadeModel( GL_SMOOTH );

    glDepthFunc(GL_LESS);
    glEnable(GL_DEPTH_TEST);

//      GLfloat white[] = { 1.0, 1.0, 1.0, 1.0 };
//      GLfloat gray[] = { 0.3, 0.3, 0.3, 1.0 };
//      GLfloat red[] = { 1.0, 0.0, 0.0, 1.0 };
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHTING);

    glLineWidth(2.0);
    glPointSize(2.0);

    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1); // Default value = 0.
   

    GLfloat light_position0[] = { 0, 0, 1.0, 0 };
    GLfloat white[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat dark_gray[] = { 0.1, 0.1, 0.1, 1.0 };
    GLfloat black[] = { 0.0, 0.0, 0.0, 1.0 };
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    // Color for the axis
    if (backcull_)
      glEnable(GL_CULL_FACE);
    else
      glDisable(GL_CULL_FACE);
    float purple[] = { 0.7, 0.0, 0.9, 1.0};
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, purple);
    camera_.use();
    glPushMatrix();
    lights_camera_.useModelView();
    glLightfv(GL_LIGHT0, GL_POSITION, light_position0);
    glPopMatrix();
    glLightfv(GL_LIGHT0, GL_AMBIENT, dark_gray);
    if (specular_) {
	glLightfv(GL_LIGHT0, GL_SPECULAR, white);
    } else {
	glLightfv(GL_LIGHT0, GL_SPECULAR, black);
    }
    data_.painter().drawScene();


   QImage res;
   res = QImage( w, h, 32 );
   glReadPixels( 0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, res.bits() );
   if ( QImage::systemByteOrder() == QImage::BigEndian ) {
      // OpenGL gives RGBA; Qt wants ARGB
      uint *p = (uint*)res.bits();
      uint *end = p + w*h;
      if ( withAlpha && format().alpha() ) {
	 while ( p < end ) {
	    uint a = *p << 24;
	    *p = (*p >> 8) | a;
	    p++;
	 }
      }
      else {
	 while ( p < end )
	    *p++ >>= 8;
      }
   }
   else {
      // OpenGL gives ABGR (i.e. RGBA backwards); Qt wants ARGB
      res = res.swapRGB();
   }
   res.setAlphaBuffer( withAlpha && format().alpha() );
   pb.Deactivate();
   camera_=oldcam;

   res.mirror().save(filename, "PNG", 100);
#endif
}

//===========================================================================
void gvView::pick(int mousex, int mousey)
//===========================================================================
{
//   std::cout << "Picking object ... " << std::endl;

    const int bufsize = 512;
    unsigned int select_buffer[bufsize];
    int hits;
    glSelectBuffer(bufsize, select_buffer);
    glRenderMode(GL_SELECT);
    glInitNames();
    glPushName(0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    camera_.pick(mousex - 1, mousey - 1, 2, 2);
    data_.painter().drawScene();
    hits = glRenderMode(GL_RENDER);
    //std::cout << "Hits: " << hits << std::endl;
    unsigned int name, numnames, z1, z2;
    unsigned int min_ind = 0; // Arbitrarily initialized in order to
			      // avoid complaints from the
			      // compiler. @@@jbt
    unsigned int min_z1 = select_buffer[1] + 1; // If hits == 0, this
						// makes no sense, but
						// at least it's
						// initialized with an
						// otherwise
						// meaningful value.
    unsigned int* bufp = select_buffer;
    // We're picking the nearest object.
    for (int i = 0; i < hits; ++i) {
	numnames = *bufp++;
	z1 = *bufp++;
	if (z1 < min_z1) {
	    min_z1 = z1;
	    min_ind = (unsigned int)(bufp - select_buffer + 1);
	}
	z2 = *bufp++;
	while(numnames--){
	    name = *bufp++;
	}
    }
//     std::cout << "Hit name: " << select_buffer[3] << std::endl;
    if (hits > 0)
	emit objectPicked(select_buffer[min_ind]);
    updateGL();
}


//===========================================================================
void gvView::pickRegion(int mousex, int mousey, int w, int h)
//===========================================================================
{
  makeCurrent();
    const int bufsize = 4096;
    unsigned int select_buffer[bufsize];
    int hits;
    glSelectBuffer(bufsize, select_buffer);
    glRenderMode(GL_SELECT);
    glInitNames();
    glPushName(0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    camera_.pick(mousex+w/2, mousey+h/2, w, h);
    data_.painter().drawScene();
    hits = glRenderMode(GL_RENDER);
    //std::cout << "Hits: " << hits << std::endl;
    
    vector<unsigned int> hitnames;
    unsigned int name, numnames, z1, z2;
    unsigned int* bufp = select_buffer;

    // [0x6]
    for(int i = 0; i < hits; ++i) {
	numnames = *bufp++;
	z1 = *bufp++;
	z2 = *bufp++;
	while(numnames--){
	    name = *bufp++;
	    hitnames.push_back(name);
	}
    }

    unsigned int* names = (hitnames.empty() ? NULL : &hitnames[0]);
    emit objectsPicked(names, (int)hitnames.size());
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    updateGL();
}


//===========================================================================
void gvView::getWindowCoords(const Vector3D& pt, int& mousex, int& mousey) const
//===========================================================================
{
   Vector3D window_coord = camera_.objectCoordsToWindowCoords(pt);
   mousex = (int)window_coord[0];
   mousey = (int)window_coord[1];
}


//===========================================================================
gvTexture* gvView::makeCoarseCheckImage()
//===========================================================================
{
    const int coarse_width = 256;
    const int coarse_height = 256;
//     GLuint texName;
    GLubyte check_image[coarse_width][coarse_height][4];
  
    for (int i=0; i<coarse_width; i++) {
	for (int j=0; j<coarse_height; j++) {
	    check_image[i][j][0]= (GLubyte) j;
	    check_image[i][j][1]= (GLubyte) 0;
	    check_image[i][j][2]= (GLubyte) i;
	    check_image[i][j][3]= (GLubyte) 0;
	}
    }

    gvTexture *new_tex = new gvTexture(coarse_height, coarse_width, check_image[0][0]);
    new_tex->setEnvMode(envReplace);
    new_tex->setWrapMode(wrapClamp);
    new_tex->setMinFilter(minNearest);
    new_tex->setMagFilter(magNearest);
    return new_tex;
}

//===========================================================================
gvTexture *gvView::makeFineCheckImage()
//===========================================================================
{
    const int fine_width = 256;
    const int fine_height = 256;
//     GLuint texName;
    GLubyte check_image[fine_width][fine_height][4];
  
    for (int i=0; i<fine_width; i++) {
	for (int j=0; j<fine_height; j++) {
	    check_image[i][j][0]= (GLubyte) 0;
	    check_image[i][j][1]= (GLubyte) j;
	    check_image[i][j][2]= (GLubyte) 0;
	    check_image[i][j][3]= (GLubyte) i;
	}
    }

    gvTexture *new_tex=new gvTexture(fine_height, fine_width, check_image[0][0]);
    new_tex->setEnvMode(envReplace);
    new_tex->setWrapMode(wrapRepeat);
    new_tex->setMinFilter(minNearest);
    new_tex->setMagFilter(magNearest);
    return new_tex;
}


//===========================================================================
void gvView::getObjAndParam(int mousex, int mousey, 
		    shared_ptr<const Go::ParamSurface> &obj,
		    double &tex_u, double &tex_v)
//===========================================================================
{
  obj.reset();
  shared_ptr<Go::GeomObject> object;
  shared_ptr<const Go::ParamSurface> paramobj;
  double tm[16];

  glClear(GL_DEPTH_BUFFER_BIT);
  if (wireframe_) 
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  
  GLboolean has_lighting=glIsEnabled(GL_LIGHTING);
  glDisable(GL_LIGHTING);
  
  std::vector<int > sel_objs;
  for (int i = 0; i < data_.numObjects(); ++i) {
    object=data_.object(i);
    paramobj = 
      dynamic_pointer_cast<const Go::ParamSurface, Go::GeomObject>(object);
    if (data_.getVisibleStateObject(i) &&
	data_.getSelectedStateObject(i) &&
	(paramobj.get()!=0))
      {
	sel_objs.push_back(i);
      }
  }

  int cx, cy, cw, ch;
  camera_.getViewPort(cx, cy, cw, ch);
  if (sel_objs.size()>0 &&
      mousex>=0 && mousex<cw && mousey>=0 && mousey<ch)
    {
      glMatrixMode(GL_PROJECTION);
      glPushMatrix();
      glMatrixMode(GL_MODELVIEW);
      glPushMatrix();

      //glEnable(GL_SCISSOR_TEST);
      //glScissor(mousex,ch-mousey-1, 1, 1);

      //drawObjs(sel_objs);
      //glDisable(GL_SCISSOR_TEST);
      
      if (coarseTex_==0)
	coarseTex_=makeCoarseCheckImage();
      if (fineTex_==0)
	fineTex_=makeFineCheckImage();
      
      glClear(GL_DEPTH_BUFFER_BIT);
      if (wireframe_) 
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      // @afr: What was the code below doing here?
      // Seems to me that it is dangerously overlapping with
      // the code above. It only worked since it was in a
      // different scope.
//       GLboolean has_lighting=glIsEnabled(GL_LIGHTING);
//       glDisable(GL_LIGHTING);

      glEnable(GL_TEXTURE_2D);
      glMatrixMode(GL_TEXTURE);
      glPushMatrix();
      glLoadIdentity();

      static int s1=1;
      static int s2=1;
      static int s3=1;

      for (size_t ki=0; ki<sel_objs.size(); ki++)
      {
	double min_s, max_s, min_t, max_t;
	texCoordInt((int)sel_objs.size(), (int)ki, min_s, max_s, min_t, max_t);
	glMatrixMode(GL_TEXTURE);
	glPushMatrix();
	if (s1)
	   {
	       glTranslatef((GLfloat)min_s, (GLfloat)min_t, 0.0f); 
	       glScalef((GLfloat)(max_s-min_s), (GLfloat)(max_t-min_t), 1.0f); 
	glGetDoublev(GL_TEXTURE_MATRIX, tm);
	coarseTex_->setTextureMatrix(tm);
	   }
	glMatrixMode(GL_MODELVIEW);
	data_.painter().getPaintable(sel_objs[ki]).paint(coarseTex_);
	glMatrixMode(GL_TEXTURE);
	glPopMatrix();
      }

      glMatrixMode(GL_TEXTURE);
      glPopMatrix();

      glEnable(GL_BLEND);
      glBlendFunc(GL_ONE, GL_ONE);

      glPushMatrix();
      glLoadIdentity();
      if (s2)
	 {
      glScalef(256.0, 256.0, 1.0f); 
	 }
      glClear(GL_DEPTH_BUFFER_BIT);

      for (size_t ki=0; ki<sel_objs.size(); ki++)
      {
	double min_s, max_s, min_t, max_t;
	texCoordInt((int)sel_objs.size(), (int)ki, min_s, max_s, min_t, max_t);
	glMatrixMode(GL_TEXTURE);
	glPushMatrix();
	if (s3){
	glTranslatef((GLfloat)min_s, (GLfloat)min_t, 0.0f); 
	glScalef((GLfloat)(max_s-min_s), (GLfloat)(max_t-min_t), 1.0f); 
	glGetDoublev(GL_TEXTURE_MATRIX, tm);
	fineTex_->setTextureMatrix(tm);
	}
	glMatrixMode(GL_MODELVIEW);
	
	data_.painter().getPaintable(sel_objs[ki]).paint(fineTex_);
	glMatrixMode(GL_TEXTURE);
	glPopMatrix();
      }
      glDisable(GL_BLEND);
    }
  glDisable(GL_TEXTURE_2D);
	
  glMatrixMode(GL_TEXTURE);
  glPopMatrix();

  glMatrixMode(GL_MODELVIEW);

  if (has_lighting)
    glEnable(GL_LIGHTING);
  
  if (wireframe_) 
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  
  GLfloat pixels[4*9];
  glReadPixels(mousex, ch-mousey-1, 1, 1, GL_RGBA, GL_FLOAT, pixels);

  GLfloat depth[9];
  glReadPixels(mousex, ch-mousey-1, 1, 1, GL_DEPTH_COMPONENT, 
	       GL_FLOAT, depth);


  if (depth[0]==1.0)
  {
    return;
  }

  double min_s, max_s, min_t, max_t;
  size_t ki = 0;
  for (; ki<sel_objs.size(); ki++)
  {
      texCoordInt((int)sel_objs.size(), (int)ki, min_s, max_s, min_t, max_t);

    if (pixels[0]>=min_s && pixels[0]<(max_s==1.0?1.1:max_s)
	&& pixels[2]>=min_t && pixels[2]<(max_t==1.0?1.1:max_t))
      break;
  }
      
  if (ki >=sel_objs.size())
  {
    std::cout << pixels[0]*255.0 << " " << pixels[2]*255.0 
	      << "                  oops" << std::endl;
    return;
  }

  object=data_.object(sel_objs[ki]);
  obj = 
    dynamic_pointer_cast<const Go::ParamSurface, Go::GeomObject>(object);
  
  double du_fine;
  int interval_u = (int)((min_s+pixels[0]*(max_s-min_s))*255.0);
  interval_u = (int)(pixels[0]*255.0);
  du_fine=(pixels[1])/257.0;
  tex_u=(((double)interval_u)/256.0+du_fine+1.0/256.0/256.0 - min_s)/
    (max_s-min_s);

  double dv_fine;
  int interval_v;
  interval_v = (int)(pixels[2]*255.0);
  dv_fine=(pixels[3])/257.0;

  tex_v=(((double)interval_v)/256.0+dv_fine+1.0/256.0/256.0 - min_t)/
    (max_t-min_t);
}

//===========================================================================
void gvView::texCoordInt(int num_objs, int ind, double &min_s, double &max_s,
			 double &min_t, double &max_t)
//===========================================================================
{
  min_s=0.0;
  max_s=1.0;
  min_t=0.0;
  max_t=1.0;

  int i;
  int mask;
  int ki;
  int num1=num_objs-1;

  int numB;
  for (numB=1, mask=1; (num1 >> (numB)) != 0; numB++)
    {
      mask <<=1;
    }

  for (ki=0; mask!=0; ki++, mask >>=1)
    {
      i=(ind) & mask;
      if (ki%2==0) // modify t
	{
	  if (i)
	    {
	      min_t = (min_t+max_t)/2.0;
	    }
	  else 
	    {
	      max_t = (min_t+max_t)/2.0;
	    }
	}
      else
	{
	  if (i)
	    {
	      min_s = (min_s+max_s)/2.0;
	    }
	  else 
	    {
	      max_s = (min_s+max_s)/2.0;
	    }
	}
    }
}

//===========================================================================
bool gvView::get3Dpoint(int mousex, int mousey, Vector3D &objpt) 
//===========================================================================
{
   bool wireframe=wireframe_;
   bool axis=axis_;
   bool success=false;
   if (wireframe || axis)
   {
      setWireframe(false);
      setAxis(false);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      paintGL();
   }
  int cx, cy, cw, ch;
  camera_.getViewPort(cx, cy, cw, ch);
  if (mousex>=0 && mousex<cw && mousey>=0 && mousey<ch)
  {
     GLfloat depth[1];
     glReadPixels(mousex, ch-mousey-1, 1, 1, GL_DEPTH_COMPONENT, 
		  GL_FLOAT, depth);

     if (depth[0]!=1.0)
     {
	Vector3D eyept(mousex, ch-mousey-1, depth[0]);
	objpt=camera_.eyeCoordsToObjectCoords(eyept);
	success=true;
     }
  }
  if (wireframe || axis)
  {
     setWireframe(wireframe);
     setAxis(axis);
  }
  return success;
}



//===========================================================================
void gvView::mousePressEvent(QMouseEvent* e)
//===========================================================================
{
//   std::cout << "A mouse key was pressed!" << std::endl;
    if (selection_mode_) {
//       std::cout << "Selection mode enabled." << std::endl;
	selecting_ = true;
    }
    if (mouse_button_ == Qt::NoButton)
	mouse_button_ = e->button();
//     if (!selecting_ && (e->state() & Qt::ControlModifier)) {
//     std::cout << "selecting: " << selecting_ << std::endl;
//     std::cout << "keyboard_event_: " << keyboard_event_ << std::endl;
//     std::cout << "keyboard_is_active_: " << keyboard_is_active_ << std::endl;
//     if (keyboard_event_ != NULL) {
//       std::cout << "keyboard_mod: " << keyboard_event_->modifiers() << std::endl;
//     }
//     if (!selecting_ && (e->buttons() & Qt::ControlModifier)) {
    if (!selecting_ && (e->modifiers() & Qt::ControlModifier)) {
//       std::cout << "We should be picking an object now ..." << std::endl;
      pick(e->x(), e->y());
    }
    last_mouse_pos_ = e->pos();
    starting_mouse_pos_ = e->pos();
    mouse_is_active_ = true;

    if (get_click_mode_) {
       get_click_mode_ = false;
       emit feedback(e->x(), e->y());
    }
}


//===========================================================================
void gvView::mouseReleaseEvent(QMouseEvent* e)
//===========================================================================
{
    if (mouse_button_ == e->button()) {
	mouse_button_ = Qt::NoButton;
	if (selection_mode_) {
	    selecting_ = false;
	    // Compute picking rectangle
	    int xmin = min(last_mouse_pos_.x(), starting_mouse_pos_.x());
	    int xmax = max(last_mouse_pos_.x(), starting_mouse_pos_.x());
	    int ymin = min(last_mouse_pos_.y(), starting_mouse_pos_.y());
	    int ymax = max(last_mouse_pos_.y(), starting_mouse_pos_.y());
	    if (xmax == xmin) xmax = xmin + 2;
	    if (ymax == ymin) ymax = ymin + 2;
	    pickRegion(xmin, ymin, xmax - xmin, ymax - ymin);
	}
    }
    mouse_is_active_ = false;
}


//===========================================================================
void gvView::mouseMoveEvent(QMouseEvent* e)
//===========================================================================
{
    if (!mouse_is_active_) return;
    // Handling mouse movement while selecting first,
    // so that we don't have to touch the rest of this function.
    if (selecting_) {
	last_mouse_pos_ = e->pos();
	updateGL();
	return;
    }

    if (feedback_mode_) {
	last_mouse_pos_ = e->pos();
	emit feedback(last_mouse_pos_.x() - starting_mouse_pos_.x(),
		      last_mouse_pos_.y() - starting_mouse_pos_.y());
	updateGL();
	return;	
    }

    if (e->buttons() & Qt::LeftButton) {
	if (e->buttons() & Qt::ShiftModifier) {
	    // Transversal rotation.
	    int dy = e->pos().y() - last_mouse_pos_.y();
	    camera_.rotateTransversal(-0.2*dy);
	} else if (e->modifiers() & Qt::AltModifier) {
	    // Rotate camera.
	    // Rotate back to state at start of mousemove,
	    // if not there already.
	    // This is to avoid hysteresis in rotations.
	    if (last_mouse_pos_ != starting_mouse_pos_) {
		// Current camera state is NOT identical to the
		// state at the start of the drag
		lights_camera_.rotate(unity_, unitx_, -0.2*draglength_);
	    }
	    // Rotate the scene
	    if (e->pos() != starting_mouse_pos_) {
		int dx = e->pos().x() - starting_mouse_pos_.x();
		int dy = e->pos().y() - starting_mouse_pos_.y();
		draglength_ = sqrt(double(dx*dx + dy*dy));
		unitx_ = double(dx)/draglength_;
		unity_ = double(dy)/draglength_;
		lights_camera_.rotate(unity_, unitx_, 0.2*draglength_);
	    }
	} else {
	    // Normal rotation.
	    // Rotate back to state at start of mousemove,
	    // if not there already.
	    // This is to avoid hysteresis in rotations.
	    if (last_mouse_pos_ != starting_mouse_pos_) {
		// Current camera state is NOT identical to the
		// state at the start of the drag
		camera_.rotate(unity_, unitx_, -0.2*draglength_);
	    }
	    // Rotate the scene
	    if (e->pos() != starting_mouse_pos_) {
		int dx = e->pos().x() - starting_mouse_pos_.x();
		int dy = e->pos().y() - starting_mouse_pos_.y();
		draglength_ = sqrt(double(dx*dx + dy*dy));
		unitx_ = double(dx)/draglength_;
		unity_ = double(dy)/draglength_;
		camera_.rotate(unity_, unitx_, 0.2*draglength_);
	    }
	}
    }
    if (e->buttons() & Qt::MidButton) {
	// Zoom the scene
	int amount = e->pos().y() - last_mouse_pos_.y();
	double dist = 0.0;
	camera_.getDistance(dist);
	camera_.setDistance(dist*(1.0 - 0.01 * amount));
    }
    if (e->buttons() & Qt::RightButton) {
	if (e->buttons() & Qt::ShiftModifier) {
	    // Depth panning
	    double dist = 0.0;
	    camera_.getDistance(dist);
	    int dy = e->pos().y() - last_mouse_pos_.y();
	    // @@ The next line is bad, but works ok most of the time.
	    camera_.moveFocalPointRelative(Vector3D(0, 0, -dy*0.1/(pow(dist,1.5))));
	} else {
	    // Pan the scene (move the focal point)
	    QPoint diff = e->pos() - last_mouse_pos_;
	    camera_.moveFocalPointRelative(Vector3D(diff.x(),
						      -diff.y(), // y is flipped
						      0.0));     // near-plane
	}
    }

    last_mouse_pos_ = e->pos();
    updateGL();
}


// //===========================================================================
// void gvView::keyPressEvent(QKeyEvent* e)
// //===========================================================================
// {
//     std::cout << "A keyboard key was pressed!" << std::endl;
// //     if (keyboard_modifier_ == Qt::NoModifier) {
//     if (keyboard_event_ == NULL) {
// //       keyboard_modifier_ = e->modifiers();
//       keyboard_event_ = e;
//     }
//     keyboard_is_active_ = true;
// }


// //===========================================================================
// void gvView::keyReleaseEvent(QKeyEvent* e)
// //===========================================================================
// {
//     std::cout << "A keyboard key was released!" << std::endl;
// //     if (keyboard_modifier_ & e->modifiers()) {
//     if (keyboard_event_ != NULL) {
// //       keyboard_modifier_ = Qt::NoModifier;
//       keyboard_event_ = NULL;
//     }
//     keyboard_is_active_ = false;
// }



//===========================================================================
void gvView::focusOnBox()
//===========================================================================
{
    double dist = 10.0;
    if ((!focus_on_origin_) && (data_.numObjects() != 0)) {
	Go::Point mid = 0.5*box_.low() + 0.5*box_.high();
	Go::Point size = box_.high() - box_.low();
	dist = max(size[0], max(size[1], size[2])) * 4.0;
	camera_.setFocalPoint(Vector3D(mid[0], mid[1], mid[2]));
// #ifndef NDEBUG
// 	MESSAGE("Focusing on the origin!");
// 	camera_.setFocalPoint(Vector3D(0.0, 0.0, 0.0));
// #endif
    } else {
	if (data_.numObjects() != 0)
	{
	    Go::Point size = box_.high() - box_.low();
	    dist = max(size[0], max(size[1], size[2])) * 4.0;
	}
	camera_.setFocalPoint(Vector3D(0.0, 0.0, 0.0));
    }
    //    std::cout << dist << std::endl;
    camera_.resetRotation();
    lights_camera_.resetRotation();
    camera_.setDistance(dist);
    base_axis_size_ = dist/5.0;
    if (axis_) {
	camera_.setAxisSize(base_axis_size_);
    }
}

//===========================================================================
void gvView::focusOnVisible()
//===========================================================================
{
   double dist = 10.0;
   std::vector<int> objs;
   for (int i=0; i<data_.numObjects(); i++)
   {
      if (data_.getVisibleStateObject(i))
      {
	 objs.push_back(i);
      }
   }
   if (objs.size()>0)
   {
      Go::BoundingBox box=data_.boundingBox(objs);
      Go::Point mid = 0.5*box.low() + 0.5*box.high();
      Go::Point size = box.high() - box.low();
      dist = max(size[0], max(size[1], size[2])) * 4.0;
      camera_.setFocalPoint(Vector3D(mid[0], mid[1], mid[2]));
   } else {
	camera_.setFocalPoint(Vector3D(0.0, 0.0, 0.0));
   }
   //    std::cout << dist << std::endl;
   camera_.resetRotation();
   lights_camera_.resetRotation();
   camera_.setDistance(dist);
   base_axis_size_ = dist/5.0;
   if (axis_) {
      camera_.setAxisSize(base_axis_size_);
   }
}

//===========================================================================
void gvView::drawOverlay()
//===========================================================================
{
     if (selecting_) {
        // @@ Paint a rectangle
        // Compute picking rectangle
        int xmin = min(last_mouse_pos_.x(), starting_mouse_pos_.x());
        int xmax = max(last_mouse_pos_.x(), starting_mouse_pos_.x());
        int ymin = min(last_mouse_pos_.y(), starting_mouse_pos_.y());
        int ymax = max(last_mouse_pos_.y(), starting_mouse_pos_.y());
        if (xmax == xmin) xmax = xmin + 2;
        if (ymax == ymin) ymax = ymin + 2;

        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadIdentity();
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glLoadIdentity();
        glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);

        // Set farge...
        float umin=(float)(-1.0 + 2.0*xmin/(float)width());
        float umax=(float)(-1.0 + 2.0*xmax/(float)width());
        float vmin=(float)( 1.0 - 2.0*ymin/(float)height());
        float vmax=(float)( 1.0 - 2.0*ymax/(float)height());
        glBegin(GL_LINE_STRIP);
	glColor4f(1.0f, 1.0f, 0.0f, 1.0f);
        glVertex3f(umin, vmin, 0.1f);
        glVertex3f(umin, vmax, 0.1f);
        glVertex3f(umax, vmax, 0.1f);
        glVertex3f(umax, vmin, 0.1f);
        glVertex3f(umin, vmin, 0.1f);
        glEnd();
        
        glMatrixMode(GL_PROJECTION);
        glPopMatrix();
        glMatrixMode(GL_MODELVIEW);
        glPopMatrix();
        glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);


/*
        QPainter painter(this);
        painter.setPen(Qt::yellow);
        painter.drawRect(xmin, ymin, xmax - xmin, ymax - ymin);
*/
     }
}

//===========================================================================
void gvView::setBlendingmode(bool mode)
//===========================================================================
{
    if (mode) {
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    } else {
	glDisable(GL_BLEND);
    }
    blending_mode_ = mode;
    paintGL();
    //    observedChanged();
}

//===========================================================================
void gvView::setOriginFocalPoint(bool origin)
//===========================================================================
{
//    MESSAGE("Under construction!");
    focus_on_origin_ = origin;
    if (focus_on_origin_)
    {
	Vector3D origin(0.0, 0.0, 0.0);
	camera_.setFocalPoint(origin);
    }
    else
    {
	focusOnBox();
    }
    updateGL();
}

//===========================================================================
void gvView::setSelectionmode(bool mode)
//===========================================================================
{
    if (selecting_) return;
    selection_mode_ = mode;
}

//===========================================================================
void gvView::setFeedbackmode(bool mode)
//===========================================================================
{
    feedback_mode_ = mode;
}

//===========================================================================
void gvView::setGetClickmode(bool mode)
//===========================================================================
{
    get_click_mode_ = mode;
}

//===========================================================================
void gvView::setCenter(int x, int y)
//===========================================================================
{
   Vector3D objpt;
   if (get3Dpoint(x, y, objpt))
   {
      camera_.setFocalPoint(objpt);
      updateGL();
   }
}


//===========================================================================
void gvView::setWireframe(bool mode)
//===========================================================================
{
    if (mode == wireframe_) return;
    if (mode) {
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    } else {
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    }
    wireframe_ = mode;
    updateGL();
}

//===========================================================================
void gvView::setAxis(bool mode)
//===========================================================================
{
    if (mode == axis_) return;
    if (mode) {
	camera_.setAxisSize(base_axis_size_);
    } else {
	camera_.setAxisSize(0.0);
    }
    axis_ = mode;
    updateGL();
}

//===========================================================================
void gvView::setBackCull(bool mode)
//===========================================================================
{
    if (mode == backcull_) return;
    backcull_ = mode;
    updateGL();
}

//===========================================================================
void gvView::setSpecular(bool mode)
//===========================================================================
{
    if (mode == specular_) return;
    specular_ = mode;
    updateGL();
}

//===========================================================================
void gvView::setPerspective(bool mode)
//===========================================================================
{
    if (mode == perspective()) return;
    camera_.setPerspectiveMode(mode);
    updateGL();
}

//===========================================================================
void gvView::observedChanged()
//===========================================================================
{
    box_ = data_.boundingBox();
    if (data_.numObjects()==0)
      no_data_=true; 
    if (no_data_ && data_.numObjects()>0) {
      focusOnBox();
      no_data_ = false;
    }
    if (gl_initialized_)
	updateGL();

}

//===========================================================================
QSize gvView::sizeHint() const
//===========================================================================
{
    return QSize(300, 300);
}
