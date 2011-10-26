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


