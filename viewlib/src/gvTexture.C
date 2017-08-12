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

#include "GoTools/viewlib/gvTexture.h"
#include <stdio.h>
#include <stdlib.h> 
#include <string>
#include <iostream>
#include "GoTools/utils/errormacros.h"
#include <exception>

#include <QImage>
#include <QString>
#include "GoTools/viewlib/raster.h"
#include <QtOpenGL/QtOpenGL>
#ifdef _MSC_VER
#ifndef NOMINMAX
#define NOMINMAX
#endif
#    include <windows.h>
#endif
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

const int MIN_FILTER_LIST[] = {
    GL_NEAREST, GL_LINEAR,
    GL_NEAREST_MIPMAP_NEAREST, GL_NEAREST_MIPMAP_LINEAR,
    GL_LINEAR_MIPMAP_NEAREST, GL_LINEAR_MIPMAP_LINEAR
};
const int MAG_FILTER_LIST[] = {
    GL_NEAREST, GL_LINEAR
};
const int ENV_MODE_LIST[] = {
    GL_DECAL, GL_REPLACE, GL_MODULATE, GL_BLEND
};
const int WRAP_MODE_LIST[] = {
  GL_CLAMP, GL_REPEAT
};


using namespace std;

unsigned *
read_rgb_texture(const char *name, int *width, int *height, int *components);

gvTexture::gvTexture(bool mipmapped) :
  height_(0), width_(0), sx_(1.0), sy_(1.0),
  pixels_(0), texName_(0),
  min_filter_(mipmapped?GL_NEAREST_MIPMAP_NEAREST:GL_NEAREST), 
  mag_filter_(GL_NEAREST),
  env_mode_(GL_REPLACE), wrap_mode_(GL_CLAMP), 
  center_edge_textels_(false), mipmapped_(mipmapped)
{
   setIdentMatrix();
}


gvTexture::gvTexture(string filename, bool mipmapped) :
  height_(0), width_(0), sx_(1.0), sy_(1.0),
  pixels_(0), texName_(0),
  min_filter_(mipmapped?GL_NEAREST_MIPMAP_NEAREST:GL_NEAREST), 
  mag_filter_(GL_NEAREST),
  env_mode_(GL_REPLACE), center_edge_textels_(false), mipmapped_(mipmapped)
{
   setIdentMatrix();
  readFile(filename);
}


gvTexture::gvTexture(const gvTexture &other) :
  pixels_(0), texName_(0),
  min_filter_(other.min_filter_), 
  mag_filter_(other.mag_filter_),
  env_mode_(other.env_mode_), center_edge_textels_(false),
  mipmapped_(other.mipmapped_)
{
  init(other);
}

gvTexture::gvTexture(int height, int width, unsigned char *pixels,
		     bool mipmapped) :
  height_(height), width_(width), sx_(1.0), sy_(1.0),
  pixels_(0), texName_(0),
  min_filter_(mipmapped?GL_NEAREST_MIPMAP_NEAREST:GL_NEAREST), 
  mag_filter_(GL_NEAREST),
  env_mode_(GL_REPLACE), center_edge_textels_(false), 
  mipmapped_(mipmapped)
{
   setIdentMatrix();
  pixels_=(unsigned *)malloc(height_*width_*sizeof(unsigned));
  memcpy(pixels_, pixels, height_*width_*sizeof(unsigned));

  mirrorIfNessecary();
  genTexture();
}


void gvTexture::getQImage(QImage &res, bool withAlpha) const
{
   res = QImage( width_, height_, QImage::Format_ARGB32 );

//    if ( QImage::systemByteOrder() == QImage::BigEndian )
//      {
       // OpenGL gives RGBA; Qt wants ARGB
       uint *p = (uint*)res.bits();
       uint *p_orig = pixels_;

       uint *end = p + width_*height_;
       if ( withAlpha)
	 {
	   while ( p < end ) {
	     uint a = *p_orig << 24;
	     *p = (*p >> 8) | a;
	     p++;
	     p_orig++;
	   }
	 }
       else
	 {
	   memcpy(res.bits(), pixels_, sizeof(uint)*height_*width_);
	   while ( p < end )
	     *p++ >>= 8;
	 }
//      }
//    else
//      {
//        memcpy(res.bits(), pixels_, sizeof(uint)*height_*width_);
//        res = res.swapRGB();
//      }
}

void gvTexture::setIdentMatrix()
{
   for(int ki=0; ki<4; ki++)
      for(int kj=0; kj<4; kj++)
	 if (ki==kj)
	    tm_[ki*4+kj]=1.0;
	 else
	    tm_[ki*4+kj]=0.0;
}

void gvTexture::setTextureMatrix(const double tm[])
{
   for(int ki=0; ki<4; ki++)
      for(int kj=0; kj<4; kj++)
	 tm_[ki*4+kj]=tm[ki*4+kj];
}

gvTexture& gvTexture::operator= (const gvTexture &other)
{
  if (&other != this)
    {
      destroy();
      init(other);
    }
  return *this;
}

gvTexture::~gvTexture()
{
  destroy();
}

void gvTexture::init(const gvTexture &other)
{
  ASSERT(pixels_==0);

  height_=other.height_;
  width_=other.width_;
  sx_=other.sx_;
  sy_=other.sy_;
  pixels_=0; 
  texName_=0;
  min_filter_=other.min_filter_;
  mag_filter_=other.mag_filter_;
  env_mode_=other.env_mode_;
  center_edge_textels_ = other.center_edge_textels_;

  pixels_=(unsigned *)malloc(height_*width_*sizeof(unsigned));
  memcpy(pixels_, other.pixels_, height_*width_*sizeof(unsigned));

  genTexture();
}

void gvTexture::destroy()
{
  if (pixels_)
    free(pixels_);
  pixels_=0;
  if (glIsTexture(texName_))
    {
      glDeleteTextures(1, &texName_);
      texName_=0;
    }
}


void gvTexture::readFile(string filename)
{
    bool done=false;
    if (filename.size()>4)
    {
	string extension(filename.end()-4, filename.end());
	if (extension[0]=='.' &&
	    toupper(extension[1])=='R' &&
	    (toupper(extension[2])=='G') &&
	    (toupper(extension[3])=='B'))
	{	  // IRIX RGB format
	    done = true;
	    int components;
	    pixels_=read_rgb_texture(filename.c_str(), &width_, &height_, 
				     &components);
	    if (pixels_==0)
	    {
		THROW ("Unable to open file:"<<filename);
	    }
  	    mirrorIfNessecary();
	    genTexture();
	}
	// Other .extension - checks here.
	// else if...
    }

    QImage qtex1, qtex2;

    if (qtex1.load(filename.c_str()) ) {	// Load first image from file
       done=true;
       width_=qtex1.width();
       height_=qtex1.height();
       int maxSize;
       glGetIntegerv(GL_MAX_TEXTURE_SIZE, &maxSize);
       //       maxSize=min(512, maxSize);
       maxSize/=2;
       while (width_>maxSize)
       {
	  width_/=2;
	  qtex1=qtex1.scaled(width_, height_);
       }
       while (height_>maxSize)
       {
	  height_/=2;
	  qtex1=qtex1.scaled(width_, height_);
       }
//        qtex1.setAlphaBuffer(true); // Removed in qt4
       qtex2 = QGLWidget::convertToGLFormat( qtex1 );  // flipped 32bit RGBA
       pixels_=new unsigned int[height_*width_];
       memcpy(pixels_, qtex2.bits(), height_*width_*sizeof(unsigned int));
       mirrorIfNessecary();
       genTexture();
    }


    if(!done)	// Magic-number - checks here.
    {
	// Peek at the file to get magic number...

	// I know, I know... I should use those C++ streams,
	// but since stdio.h already was included, I just kept
	// to good old FILE* :)
	FILE* f = fopen(filename.c_str(), "r");
	if (f==0)
	{
	    THROW ("Unable to open file:"<<filename);
	}
	char buff[10];
	fgets(buff, sizeof(buff), f);
	fclose(f);

	if((buff[0] == 'P') && (buff[1] == '6'))
	{   // PPM binary/ascii, components=255
//	    pixels_ = (unsigned*)read_ppm_file(filename.c_str(), &width_, &height_);
	    // Adding alpha=0.
	    unsigned char* rgbbuff = read_ppm_file(filename.c_str(), &width_, &height_);
	    if (rgbbuff==0)
	    {
		THROW ("Unable to open file:"<<filename);
	    }
	    unsigned char* rgbabuff = new unsigned char[width_*height_*4];
	    for(int w = 0; w < width_; ++w)
		for(int h = 0; h < height_; ++h)
		{
		    int a = w + h*width_;
		    int b = w + (height_-h-1)*width_;
		    rgbabuff[4*b+0] = rgbbuff[3*a+0];
		    rgbabuff[4*b+1] = rgbbuff[3*a+1];
		    rgbabuff[4*b+2] = rgbbuff[3*a+2];
		    rgbabuff[4*b+3] = 0;
		}
	    pixels_ = (unsigned*)rgbabuff;
	    
	    mirrorIfNessecary();
	    genTexture();
	}
	// Other Magic-number - checks here.
	// else if...
    }
}

bool gvTexture::genTexture() const
{
  if (texName_==0)
    {
      if (height_<=0 || width_<=0)
	return false;
      glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
      glGenTextures(1, &texName_);
      if (texName_==0)
	return false;
      glBindTexture(GL_TEXTURE_2D, texName_);
      if (mipmapped_)
      {
	 gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGBA, width_, height_, 
			   GL_RGBA, GL_UNSIGNED_BYTE, pixels_);
      } else
      {
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA,
		   width_, height_, 0, GL_RGBA,
		   GL_UNSIGNED_BYTE, pixels_);
      }
      glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, (float)min_filter_); 
      glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, (float)mag_filter_); 
      glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, (float)env_mode_);
      if (!glIsTexture(texName_))
	  {
	      return false;
	  }
    }
  return true;
}


void gvTexture::setMinFilter(MinFilterSet filter, bool flush)
{
    min_filter_ = MIN_FILTER_LIST[filter];
    if(flush)
     {
	bind();
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
			(float)MIN_FILTER_LIST[filter]); 
     }
}

void gvTexture::setMagFilter(MagFilterSet filter, bool flush)
{
    mag_filter_ = MAG_FILTER_LIST[filter];
    if(flush)
     {
	bind();
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
			(float)MAG_FILTER_LIST[filter]); 
     }
}

void gvTexture::setEnvMode(EnvModeSet mode, bool flush)
{
    env_mode_ = ENV_MODE_LIST[mode];
    if(flush)
     {
	bind();
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE,
		  (float)ENV_MODE_LIST[mode]);
     }
}

void gvTexture::setWrapMode(WrapModeSet mode, bool flush)
{
  wrap_mode_ = WRAP_MODE_LIST[mode];
  if(flush)
   {
      bind();
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, WRAP_MODE_LIST[mode]);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, WRAP_MODE_LIST[mode]);
   }
}


MinFilterSet gvTexture::getMinFilter()
{
    switch(min_filter_) {
	case GL_NEAREST: return minNearest;
	case GL_LINEAR:  return minNearest;
	case GL_NEAREST_MIPMAP_NEAREST: return minNearestMipmapNearest;
	case GL_NEAREST_MIPMAP_LINEAR:  return minNearestMipmapLinear;
	case GL_LINEAR_MIPMAP_NEAREST:  return minLinearMipmapNearest;
	case GL_LINEAR_MIPMAP_LINEAR:   return minLinearMipmapLinear;
	default:  return minNearest;
    }
}

MagFilterSet gvTexture::getMagFilter()
{
    switch(mag_filter_) {
	case GL_NEAREST: return magNearest;
	case GL_LINEAR:  return magNearest;
	default:  return magNearest;
    }
}

EnvModeSet gvTexture::getEnvMode()
{
    switch(env_mode_) {
	case GL_DECAL:     return envDecal;
	case GL_REPLACE:   return envReplace;
	case GL_MODULATE:  return envModulate;
	case GL_BLEND:     return envBlend;
	default:  return envReplace;
    }
}



void gvTexture::setAlphaValue(int val)
{
    if (val < 0 || val > 255) {
	THROW("Argument out of range in gvTexture::setAlphaValue(...).");
    }
    int size = width() * height();
    for (int pos = 0; pos < size; ++pos) {
	reinterpret_cast<unsigned char*> (pixels_ + pos)[3] = (unsigned char)val;
    }
}

bool gvTexture::bind() const
{
  if (!glIsTexture(texName_))
      {
	  bool success = genTexture();
	  if (!success) {
	      return false;
	  }
      }
  glBindTexture(GL_TEXTURE_2D, texName_);

  glMatrixMode(GL_TEXTURE);
  glLoadIdentity();
  glScalef((GLfloat)sx_, (GLfloat)sy_, 1.0f); 
  glMultMatrixd(tm_);
  if(center_edge_textels_)
  {
      double onetextel_w = 1.0/width_;
      double onetextel_h = 1.0/height_;
      glTranslatef((GLfloat)onetextel_w/2, (GLfloat)onetextel_h/2, 0.0);
      glScalef((GLfloat)(1.0-onetextel_w), (GLfloat)(1.0-onetextel_h), 1.0);
  }

  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, (float)min_filter_); 
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, (float)mag_filter_); 
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, (float)env_mode_);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
		  (float)min_filter_); 
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
		  (float)mag_filter_); 
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE,
	    (float)env_mode_);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, wrap_mode_);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, wrap_mode_);

  glMatrixMode(GL_MODELVIEW);
  return true;
}

void gvTexture::mirrorIfNessecary()
{
  int h;
  int w;
  int i,j;
  unsigned *old_img;
  w=width_;
  while (w%2==0 && w>0)
    {
      w/=2;
    }
  if (w==1)
    {
      w=width_;
    } else
    {
      w=1;
      while (w<width_)
	{
	  w*=2;
	}
    }
  h= height_;
  while (h%2==0 && h>0)
    {
      h/=2;
    }
  if (h==1)
    {
      h=height_;
    } else
    {
      h=1;
      while (h<height_)
	h*=2;
    }
  if (h==height_ && w==width_) 
      return;

  old_img=pixels_;
  pixels_ = (unsigned *)malloc(w*h*sizeof(unsigned));
  
  for(i=0; i<height_; i++)
    {
      for(j=0; j<width_; j++)
	{
	  pixels_[i*w+j]=old_img[i*width_+j];
	}
    }


  for (i=0; i<height_; i++)
    {
      for(j=0; j<w-width_; j++)
	{
	  pixels_[i*w+width_+j]=old_img[i*width_+width_-j-1];
	}
    }

  for (i=0; i<h-height_; i++)
    {
      for(j=0; j<width_; j++)
	{
	  pixels_[(height_+i)*w+j]=old_img[(height_-i-1)*width_+j];
	}
    }

  for (i=0; i<h-height_; i++)
    {
      for(j=0; j<w-width_; j++)
	{
	  pixels_[(height_+i)*w+width_+j]=
	    old_img[(height_-i-1)*width_+width_-j-1];
	}
    }
  
  sx_=((float)width_) /((float)w);
  sy_=((float)height_) /((float)h);

  free(old_img);

  height_=h;
  width_=w;
}

void
bwtorgba(unsigned char *b,unsigned char *l,int n) {
    while(n--) {
	l[0] = *b;
	l[1] = *b;
	l[2] = *b;
	l[3] = 0xff;
	l += 4; b++;
    }
}

void
latorgba(unsigned char *b, unsigned char *a,unsigned char *l,int n) {
    while(n--) {
	l[0] = *b;
	l[1] = *b;
	l[2] = *b;
	l[3] = *a;
	l += 4; b++; a++;
    }
}

void
rgbtorgba(unsigned char *r,unsigned char *g,unsigned char *b,unsigned char *l,int n) {
    while(n--) {
	l[0] = r[0];
	l[1] = g[0];
	l[2] = b[0];
	l[3] = 0xff;
	l += 4; r++; g++; b++;
    }
}

void
rgbatorgba(unsigned char *r,unsigned char *g,unsigned char *b,unsigned char *a,unsigned char *l,int n) {
    while(n--) {
	l[0] = r[0];
	l[1] = g[0];
	l[2] = b[0];
	l[3] = a[0];
        l += 4; r++; g++; b++; a++;
    }
}


typedef struct _ImageRec {
    unsigned short imagic;
    unsigned short type;
    unsigned short dim;
    unsigned short xsize, ysize, zsize;
    unsigned int min, max;
    unsigned int wasteBytes;
    char name[80];
    unsigned long colorMap;
    FILE *file;
    unsigned char *tmp, *tmpR, *tmpG, *tmpB;
    unsigned long rleEnd;
    unsigned int *rowStart;
    int *rowSize;
} ImageRec;

static void
ConvertShort(unsigned short *array, long length) {
    unsigned b1, b2;
    unsigned char *ptr;

    ptr = (unsigned char *)array;
    while (length--) {
	b1 = *ptr++;
	b2 = *ptr++;
	*array++ = (short unsigned int)( (b1 << 8) | (b2) );
    }
}

static void
ConvertLong(unsigned *array, long length) {
    unsigned b1, b2, b3, b4;
    unsigned char *ptr;

    ptr = (unsigned char *)array;
    while (length--) {
	b1 = *ptr++;
	b2 = *ptr++;
	b3 = *ptr++;
	b4 = *ptr++;
	*array++ = (b1 << 24) | (b2 << 16) | (b3 << 8) | (b4);
    }
}

static ImageRec *ImageOpen(const char *fileName)
{
    union {
	int testWord;
	char testByte[4];
    } endianTest;
    ImageRec *image;
    int swapFlag;
    int x;

    endianTest.testWord = 1;
    if (endianTest.testByte[0] == 1) {
	swapFlag = 1;
    } else {
	swapFlag = 0;
    }

    image = (ImageRec *)malloc(sizeof(ImageRec));
    if (image == NULL) {
	fprintf(stderr, "Out of memory!\n");
	return 0;
    }
    if ((image->file = fopen(fileName, "rb")) == NULL) {
      return 0;
    }

    fread(image, 1, 12, image->file);

    if (swapFlag) {
	ConvertShort(&image->imagic, 6);
    }

    image->tmp = (unsigned char *)malloc(image->xsize*256);
    image->tmpR = (unsigned char *)malloc(image->xsize*256);
    image->tmpG = (unsigned char *)malloc(image->xsize*256);
    image->tmpB = (unsigned char *)malloc(image->xsize*256);
    if (image->tmp == NULL || image->tmpR == NULL || image->tmpG == NULL ||
	image->tmpB == NULL) {
	fprintf(stderr, "Out of memory!\n");
	return 0;
    }

    if ((image->type & 0xFF00) == 0x0100) {
	x = image->ysize * image->zsize * (int)sizeof(unsigned);
	image->rowStart = (unsigned *)malloc(x);
	image->rowSize = (int *)malloc(x);
	if (image->rowStart == NULL || image->rowSize == NULL) {
	    fprintf(stderr, "Out of memory!\n");
	    return 0;
	}
	image->rleEnd = 512 + (2 * x);
	fseek(image->file, 512, SEEK_SET);
	fread(image->rowStart, 1, x, image->file);
	fread(image->rowSize, 1, x, image->file);
	if (swapFlag) {
	    ConvertLong(image->rowStart, x/(int)sizeof(unsigned));
	    ConvertLong((unsigned *)image->rowSize, x/(int)sizeof(int));
	}
    } else {
	image->rowStart = NULL;
	image->rowSize = NULL;
    }
    return image;
}

static void
ImageClose(ImageRec *image) {
    fclose(image->file);
    free(image->tmp);
    free(image->tmpR);
    free(image->tmpG);
    free(image->tmpB);
    free(image->rowSize);
    free(image->rowStart);
    free(image);
}

static void
ImageGetRow(ImageRec *image, unsigned char *buf, int y, int z) {
    unsigned char *iPtr, *oPtr, pixel;
    int count;

    if ((image->type & 0xFF00) == 0x0100) {
	fseek(image->file, (long)image->rowStart[y+z*image->ysize], SEEK_SET);
	fread(image->tmp, 1, (unsigned int)image->rowSize[y+z*image->ysize],
	      image->file);

	iPtr = image->tmp;
	oPtr = buf;
	for (;;) {
	    pixel = *iPtr++;
	    count = (int)(pixel & 0x7F);
	    if (!count) {
		return;
	    }
	    if (pixel & 0x80) {
		while (count--) {
		    *oPtr++ = *iPtr++;
		}
	    } else {
		pixel = *iPtr++;
		while (count--) {
		    *oPtr++ = pixel;
		}
	    }
	}
    } else {
	fseek(image->file, 512+(y*image->xsize)+(z*image->xsize*image->ysize),
	      SEEK_SET);
	fread(buf, 1, image->xsize, image->file);
    }
}

unsigned *
read_rgb_texture(const char *name, int *width, int *height, int *components) {
    unsigned *base, *lptr;
    unsigned char *rbuf, *gbuf, *bbuf, *abuf;
    ImageRec *image;
    int y;

    image = ImageOpen(name);
    
    if(!image)
	return NULL;
    (*width)=image->xsize;
    (*height)=image->ysize;
    (*components)=image->zsize;
    base = (unsigned *)malloc(image->xsize*image->ysize*sizeof(unsigned));
    rbuf = (unsigned char *)malloc(image->xsize*sizeof(unsigned char));
    gbuf = (unsigned char *)malloc(image->xsize*sizeof(unsigned char));
    bbuf = (unsigned char *)malloc(image->xsize*sizeof(unsigned char));
    abuf = (unsigned char *)malloc(image->xsize*sizeof(unsigned char));
    if(!base || !rbuf || !gbuf || !bbuf)
      return NULL;
    lptr = base;
    for(y=0; y<image->ysize; y++) {
	if(image->zsize>=4) {
	    ImageGetRow(image,rbuf,y,0);
	    ImageGetRow(image,gbuf,y,1);
	    ImageGetRow(image,bbuf,y,2);
	    ImageGetRow(image,abuf,y,3);
	    rgbatorgba(rbuf,gbuf,bbuf,abuf,(unsigned char *)lptr,image->xsize);
	    lptr += image->xsize;
	} else if(image->zsize==3) {
	    ImageGetRow(image,rbuf,y,0);
	    ImageGetRow(image,gbuf,y,1);
	    ImageGetRow(image,bbuf,y,2);
	    rgbtorgba(rbuf,gbuf,bbuf,(unsigned char *)lptr,image->xsize);
	    lptr += image->xsize;
	} else if(image->zsize==2) {
	    ImageGetRow(image,rbuf,y,0);
	    ImageGetRow(image,abuf,y,1);
	    latorgba(rbuf,abuf,(unsigned char *)lptr,image->xsize);
	    lptr += image->xsize;
	} else {
	    ImageGetRow(image,rbuf,y,0);
	    bwtorgba(rbuf,(unsigned char *)lptr,image->xsize);
	    lptr += image->xsize;
	}
    }
    ImageClose(image);
    free(rbuf);
    free(gbuf);
    free(bbuf);
    free(abuf);

    return (unsigned *) base;
}
