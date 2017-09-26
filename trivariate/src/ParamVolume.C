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

#include "GoTools/trivariate/ParamVolume.h"
#include "GoTools/utils/Point.h"

using std::vector;
namespace Go
{


//===========================================================================
ParamVolume::~ParamVolume()

//===========================================================================
{
}

//===========================================================================
void ParamVolume::estimateVolSize(double& u_size, double& v_size, double& w_size,
				  int u_nmb, int v_nmb, int w_nmb)

//===========================================================================
{
  Array<double,6> dom = parameterSpan();
  double del_u = (dom[1] - dom[0])/(double)(u_nmb-1);
  double del_v = (dom[3] - dom[2])/(double)(v_nmb-1);
  double del_w = (dom[5] - dom[4])/(double)(w_nmb-1);

  int ki, kj, kr;
  double upar, vpar, wpar;
  vector<Point> pts(u_nmb*v_nmb*w_nmb);
  for (kr=0, wpar=dom[4]; kr<w_nmb; ++kr, wpar+=del_w)
    for (kj=0, vpar=dom[2]; kj<v_nmb; ++kj, vpar+=del_v)
      for (ki=0, upar=dom[0]; ki<u_nmb; ++ki, upar+=del_u)
	{
	  Point pos;
	  point(pos, upar, vpar, wpar);
	  pts[(kr*v_nmb+kj)*u_nmb+ki] = pos;
	}

  double acc_u=0.0, acc_v=0.0, acc_w=0.0;
  for (kr=0; kr<w_nmb; ++kr)
    for (kj=0; kj<v_nmb; ++kj)
      for (ki=1; ki<u_nmb; ++ki)
	acc_u += pts[(kr*w_nmb+kj)*u_nmb+ki-1].dist(pts[(kr*w_nmb+kj)*u_nmb+ki]);
  acc_u /= (double)(v_nmb*w_nmb);
  
  for (kr=0; kr<w_nmb; ++kr)
    for (ki=0; ki<u_nmb; ++ki)
      for (kj=1; kj<v_nmb; ++kj)
	acc_v += pts[(kr*w_nmb+kj-1)*u_nmb+ki].dist(pts[(kr*w_nmb+kj)*u_nmb+ki]);
  acc_v /= (double)(u_nmb*w_nmb);
  
  for (kj=0; kj<v_nmb; ++kj)
    for (ki=0; ki<u_nmb; ++ki)
      for (kr=1; kr<w_nmb; ++kr)
	acc_w += pts[((kr-1)*w_nmb+kj)*u_nmb+ki].dist(pts[(kr*w_nmb+kj)*u_nmb+ki]);
  acc_w /= (double)(u_nmb*v_nmb);
  
  u_size = acc_u;
  v_size = acc_v;
  w_size = acc_w;
}

} // namespace Go
