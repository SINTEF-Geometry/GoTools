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

#ifndef _GVCOLOR_H
#define _GVCOLOR_H


//#include <boost/type_traits.hpp>
#include "GoTools/utils/config.h"
#include <cmath>

/// Represents a color by four floats, like in OpenGL.
struct gvColor
{
    float rgba[4];
    /// Default color is white
    gvColor()
    {
	rgba[0] = rgba[1] = rgba[2] = rgba[3] = 1.0f;
    }
    template <typename NumericType>
    gvColor(NumericType r, NumericType g, NumericType b)
    {
	setVal(Int2Type< is_floating_point<NumericType>::value >(),
	       r, g, b);
    }
    template <typename NumericType>
    gvColor(NumericType r, NumericType g, NumericType b,
	    NumericType alpha)
    {
	setVal(Int2Type< is_floating_point<NumericType>::value >(),
	       r, g, b, alpha);
    }

    template <typename NumericType>
    static gvColor hsva(NumericType h, NumericType s, NumericType v,
			NumericType alpha)
    {
       gvColor ret;
       ret.setHsva(Int2Type< is_floating_point<NumericType>::value >(),
		   h, s, v, alpha);
       return ret;
    }

//      gvColor(float r, float g, float b, float alpha = 1.0f)
//      {
//  	rgba[0] = r;
//  	rgba[1] = g;
//  	rgba[2] = b;
//  	rgba[3] = alpha;
//      }
//      gvColor(double r, double g, double b, double alpha = 1.0)
//      {
//  	rgba[0] = static_cast<float>(r);
//  	rgba[1] = static_cast<float>(g);
//  	rgba[2] = static_cast<float>(b);
//  	rgba[3] = static_cast<float>(alpha);
//      }
//      gvColor(unsigned char r, unsigned char g, unsigned char b,
//  	    unsigned char alpha = 255)
//      {
//  	rgba[0] = r/255.0f;
//  	rgba[1] = g/255.0f;
//  	rgba[2] = b/255.0f;
//  	rgba[3] = alpha/255.0f;
//      }
    gvColor(const gvColor& other)
    {
	rgba[0] = other.rgba[0];
	rgba[1] = other.rgba[1];
	rgba[2] = other.rgba[2];
	rgba[3] = other.rgba[3];
    }
    gvColor& operator=(const gvColor& other)
    {
	rgba[0] = other.rgba[0];
	rgba[1] = other.rgba[1];
	rgba[2] = other.rgba[2];
	rgba[3] = other.rgba[3];
	return *this;
    }

private:
    template <int v>
    struct Int2Type
    {
	enum { value = v };
    };

    template <typename NumericType>
    void setVal(Int2Type<0>,
	   NumericType r,
	   NumericType g,
	   NumericType b,
	   NumericType alpha = 255)
    {
	rgba[0] = r/255.0f;
	rgba[1] = g/255.0f;
	rgba[2] = b/255.0f;
	rgba[3] = alpha/255.0f;
    }

    template <typename NumericType>
    void setVal(Int2Type<1>,
	   NumericType r,
	   NumericType g,
	   NumericType b,
	   NumericType alpha = 1.0f)
    {
	rgba[0] = (float)r;
	rgba[1] = (float)g;
	rgba[2] = (float)b;
	rgba[3] = (float)alpha;
    }

   template <typename NumericType>
   void setHsva(Int2Type<0>, 
		NumericType h,
		NumericType s,
		NumericType v,
		NumericType alpha)
   {
      float fh=h/360.0f;
      float fs=s/100.0f;
      float fv=v/100.0f;
      float fa=alpha/100.0f;
      setHsva(Int2Type< is_floating_point<float>::value >(),
	      fh,fs,fv,fa);
   }

   template <typename NumericType>
   void setHsva(Int2Type<1>, 
		NumericType h,
		NumericType s,
		NumericType v,
		NumericType alpha)
   {
      h*=6.0f;
      int i=int(floor(h));
      float f=h-i;
      if ((i&1)==0) //i even
	 f=1.0f-f;

      float m=v*(1.0f-s);
      float n=v*(1.0f-s*f);

      switch (i)
      {
      case 6:
      case 0:
	 setVal(Int2Type<1>(), v, n, m, alpha);
	 break;
      case 1:
	 setVal(Int2Type<1>(), n, v, m, alpha);
	 break;
      case 2:
	 setVal(Int2Type<1>(),m, v, n, alpha);
	 break;
      case 3:
	 setVal(Int2Type<1>(),m, n, v, alpha);
	 break;
      case 4:
	 setVal(Int2Type<1>(),n, m, v, alpha);
	 break;
      case 5:
	 setVal(Int2Type<1>(),v, m, n, alpha);
	 break; 
     }
   }


};



#endif // _GVCOLOR_H

