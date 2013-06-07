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

#ifndef _GV_TEXTURE_H
#define _GV_TEXTURE_H

#include <string>

enum MinFilterSet {
    minNearest, minLinear,
    minNearestMipmapNearest, minNearestMipmapLinear,
    minLinearMipmapNearest, minLinearMipmapLinear
};
enum MagFilterSet {
    magNearest, magLinear
};
enum EnvModeSet {
    envDecal, envReplace, envModulate, envBlend
};
enum WrapModeSet {
  wrapClamp, wrapRepeat
};

class QImage;

/** Documentation ...
    etc
 */

class gvTexture
{
 public:
  gvTexture(bool mipmapped=false);
  gvTexture(const gvTexture &other);
  gvTexture& operator= (const gvTexture &other);
  /// the size of the pixels array should be
  /// height*width*4 in RGBA format
  gvTexture(int height, int width, unsigned char *pixels, 
	    bool mipmapped=false);

  gvTexture(std::string filename, bool mipmapped=false);
  virtual void readFile(std::string filename); 
  virtual ~gvTexture();
  bool genTexture() const;
  bool bind() const;

  void getQImage(QImage &res, bool withAlpha=false) const;
  
  int height() const {return height_;} 
  int width() const {return width_;}

  // Sets the alpha (transparency) value of all pixels in the 
  // texture equal to val, where val is in the range of 0..255.
  // Any other value of val will cause a runtime error.
  void setAlphaValue(int val);
   void setTextureMatrix(const double tm[]);
   void setIdentMatrix();

  // These set-methods take a GL-constant as parameter.
  // See OpenGL specs for valid values (glTexParameter and glTexEnv).
  // They accept all valid values for OpenGL 1.2.
  // flush=true makes the appropriate GL calls, remember to bind first!
  // flush=false only caches the calls, and readFile must be called afterwards.
  void setMinFilter(MinFilterSet filter, bool flush = true);
  void setMagFilter(MagFilterSet filter, bool flush = true);
  void setEnvMode(EnvModeSet mode, bool flush = true);
  void setWrapMode(WrapModeSet mode, bool flush = true);
  MinFilterSet getMinFilter();
  MagFilterSet getMagFilter();
  EnvModeSet getEnvMode();
  void setCenterEdgeTexels(bool enable) { center_edge_textels_ = enable; }


 private:
  void mirrorIfNessecary();
  int height_, width_;
  double sx_, sy_;
  unsigned int *pixels_;
  void init(const gvTexture &other);
  void destroy();
  mutable unsigned int texName_;
  int min_filter_, mag_filter_;
  int env_mode_;
  int wrap_mode_;
  bool center_edge_textels_;
   /// tm_ is the texture matrix without handling rectangular textures
   double tm_[16];
   bool mipmapped_;
};

#endif  // _GV_TEXTURE_H


