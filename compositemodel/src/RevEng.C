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

#include "GoTools/compositemodel/RevEng.h"
#include "GoTools/compositemodel/RevEngPoint.h"
#include "GoTools/compositemodel/RevEngEdge.h"
//#include "GoTools/compositemodel/RevEngRegion.h"
#include "GoTools/compositemodel/RevEngUtils.h"
#include "GoTools/compositemodel/HedgeSurface.h"
#include "GoTools/compositemodel/ImplicitApprox.h"
#include "GoTools/compositemodel/SurfaceModelUtils.h"
#include "GoTools/utils/DirectionCone.h"
#include "GoTools/utils/MatrixXD.h"
#include "GoTools/creators/CurveCreators.h"
#include "GoTools/creators/CreatorsUtils.h"
#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/Cone.h"
#include "GoTools/geometry/Plane.h"
#include "GoTools/geometry/Torus.h"
#include "GoTools/geometry/Sphere.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/geometry/ClosestPoint.h"
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/geometry/SplineDebugUtils.h"
#include "sislP.h"
#include <vector>
#include <set>
#include <fstream>
#include <iostream> // @@ debug

using namespace Go;
using std::vector;
using std::pair;
using std::istream;
using std::ostream;

typedef MatrixXD<double, 3> Matrix3D;

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

//#define DEBUG_DIV
//#define DEBUG_EDGE0
//#define DEBUG_BLEND
//#define DEBUG_MONGE
//#define DEBUG_ENHANCE
//#define DEBUG_SEG
//#define DEBUG
//#define DEBUGONE
//#define DEBUG_CHECK
//#define DEBUG_PLANAR
//#define DEBUG_AXIS
//#define DEBUG_GROW
//#define DEBUG_VALIDATE
//#define DEBUG_EDGE
//#define DEBUG_TRIANG
//#define DEBUG_TRIM
//#define DEBUG_MODEL
//#define DEBUG_SMALL

//===========================================================================
RevEng::RevEng(shared_ptr<ftPointSet> tri_sf)
  : tri_sf_(tri_sf)
//===========================================================================
{
  mean_edge_len_ = 0.0;
  int num = tri_sf_->size();
  if (num > 0.0)
    {
      double fac1 = 1.0/(double)num;
      for (int ki=0; ki<num; ++ki)
	{
	  double tmp_len = 0.0;
	  ftSamplePoint* curr = (*tri_sf_)[ki];
	  Vector3D xyz1 = curr->getPoint();
	  vector<ftSamplePoint*> adj = curr->getNeighbours();
	  if (adj.size() > 0)
	    {
	      double fac2 = 1.0/(double)adj.size();
	      for (size_t kj=0; kj<adj.size(); ++kj)
		{
		  Vector3D xyz2 = adj[kj]->getPoint();
		  tmp_len += fac2*xyz1.dist(xyz2);
		}
	    }
	  mean_edge_len_ += fac1*tmp_len;
	}
    }
	      
	  
  
  // Set default parameters
  model_character_ = ROUGH;
  initParameters();
  max_next_ = std::min(80, tri_sf_->size()/200);
  max_next_ = std::max(2*min_next_, max_next_);
}


//===========================================================================
RevEng::RevEng()
//===========================================================================
{
  // Empty infrastructure for reading stage
  model_character_ = ROUGH;
  initParameters();
}


//===========================================================================
RevEng::~RevEng()
//===========================================================================
{
}


int compare_x_par(const RevEngPoint* p1, const RevEngPoint* p2)
{
  return (p1->getPoint()[0] < p2->getPoint()[0]);
  // if (p1->getPoint()[0] < p2->getPoint()[0])
  //   return -1;
  // else if (p1->getPoint()[0] > p2->getPoint()[0])
  //   return 1;
  // else
  //   return 0;
}

int compare_y_par(const RevEngPoint* p1, const RevEngPoint* p2)
{
  return (p1->getPoint()[1] < p2->getPoint()[1]);
  // if (p1->getPoint()[1] < p2->getPoint()[1])
  //   return -1;
  // else if (p1->getPoint()[1] > p2->getPoint()[1])
  //   return 1;
  // else
  //   return 0;
}

int compare_z_par(const RevEngPoint* p1, const RevEngPoint* p2)
{
  return (p1->getPoint()[2] < p2->getPoint()[2]);
  // if (p1->getPoint()[2] < p2->getPoint()[2])
  //   return -1;
  // else if (p1->getPoint()[2] > p2->getPoint()[2])
  //   return 1;
  // else
  //   return 0;
}

//===========================================================================
void RevEng::setBoundingBox()
//===========================================================================
{
  int nmbpt = tri_sf_->size();
  vector<RevEngPoint*> all_pts(nmbpt);
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint* pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      Vector3D xyz = pt->getPoint();
      Point xyz2 = Point(xyz[0], xyz[1], xyz[2]);
      all_pts[ki] = pt;
      if (ki == 0)
	bbox_ = BoundingBox(xyz2, xyz2);
      else
	bbox_.addUnionWith(xyz2);
    }

}

//===========================================================================
void RevEng::enhancePoints()
//===========================================================================
{
  setBoundingBox();
#ifdef DEBUG_ENHANCE  
  std::cout << "Bounding box, min: " << bbox_.low() << ", max: " << bbox_.high() << std::endl;
#endif

  // Update parameters based on surface roughness
  updateParameters();
  int nmbpt = tri_sf_->size();

#ifdef DEBUG_TRIANG
  vector<vector<RevEngPoint*> > conn_groups;
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      if (pt->visited())
	continue;
      vector<RevEngPoint*> curr_group;
      pt->fetchConnected(0, nmbpt, curr_group);
      conn_groups.push_back(curr_group);
    }
  
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      pt->unsetVisited();
    }

  std::ofstream oftri("init_ptgroups.g2");
  for (size_t kj=0; kj<conn_groups.size(); ++kj)
    {
      oftri << "400 1 0 0" << std::endl;
      oftri << conn_groups[kj].size() << std::endl;
      for (size_t kr=0; kr<conn_groups[kj].size(); ++kr)
	oftri << conn_groups[kj][kr]->getPoint() << std::endl;
    }
#endif
  
#ifdef DEBUG_DIV
  int writepoints = 0;
  vector<double> tri_ang(nmbpt);
#endif
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);

      // Compute surface normal from triangulation
      pt->computeTriangNormal(100.0*mean_edge_len_);
      if (pt->getNmbNeighbour() == 0)
	pt->setOutlier();
      //double avlen = pt->getMeanEdgLen();
#ifdef DEBUG_DIV
      tri_ang[ki] = pt->getTriangAngle();
#endif
    }
#ifdef DEBUG_ENHANCE  
  std::sort(tri_ang.begin(), tri_ang.end());
  std::cout << "Triangle angles: " << tri_ang[0] << " " << tri_ang[nmbpt/4];
  std::cout << " " << tri_ang[nmbpt/2] << " " << tri_ang[3*nmbpt/4];
  std::cout << " " << tri_ang[nmbpt-1] << std::endl;
  std::cout << "norm_ang_lim_ : " << norm_ang_lim_ << std::endl;
#endif
  
  double wgt_nmb = 1.0/(double)nmbpt;
  double av_close = 0.0;
  int max_close = 0;
  int min_close = nmbpt;
  vector<double> lambda_3(nmbpt);
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      if (pt->nmbLocFunc() > 0)
	continue;  // Already enhanced

      // // Compute surface normal from triangulation
      // pt->computeTriangNormal(100.0*mean_edge_len_);
      if (pt->isOutlier())
	continue;

      // if (pt->getNmbNeighbour() == 0)
      // 	{
      // 	  pt->setOutlier();
      // 	  continue;
      // 	}

      //double avlen = pt->getMeanEdgLen();

      //Fetch nearby points
      vector<RevEngPoint*> nearpts;
      double local_len = pt->getMeanEdgLen(10.0*mean_edge_len_);
      double radius = rfac_*(local_len + mean_edge_len_);
      double radius2 = 0.5*radius;
      radius = std::min(radius, 20.0*mean_edge_len_);
      //radius *= 1.5; // TEST 
      //double radius = 0.5*rfac_*(local_len + mean_edge_len_);
      pt->fetchClosePoints2(radius, min_next_, max_next_, nearpts);

      Point mincvec(0.0, 0.0, 0.0), maxcvec(0.0, 0.0, 0.0);
      //      Point mincvec2(0.0, 0.0, 0.0), maxcvec2(0.0, 0.0, 0.0);

      av_close += wgt_nmb*(double)nearpts.size();
      max_close = std::max(max_close, (int)nearpts.size());
      min_close = std::min(min_close, (int)nearpts.size());
      
      if (nearpts.size() >= 3)
	{
#ifdef DEBUG_DIV
	  std::ofstream of("nearpts.g2");
	  if (writepoints)
	    {
	      of << "400 1 0 4 255 0 0 255" << std::endl;
	      of << "1" << std::endl;
	      of << pt->getPoint() << std::endl << std::endl;
	      of << "400 1 0 4 0 255 0 255" << std::endl;
	      of << nearpts.size() << std::endl;
	      for (size_t kr=0; kr<nearpts.size(); ++kr)
		of << nearpts[kr]->getPoint() << std::endl;
	    }
#endif
	  
	  // Compute eigenvectors and values of covariance matrix
	  nearpts.push_back(pt);
	  double lambda[3];
	  double eigenvec[3][3];
	  RevEngUtils::principalAnalysis(nearpts, lambda, eigenvec);
	  Point eigen1(eigenvec[0][0], eigenvec[0][1], eigenvec[0][2]);
	  Point eigen2(eigenvec[1][0], eigenvec[1][1], eigenvec[1][2]);
	  Point eigen3(eigenvec[2][0], eigenvec[2][1], eigenvec[2][2]);
	  lambda_3[ki] = lambda[2];
	  Point tnorm = pt->getTriangNormal();
	  if (tnorm.length() < 1.0e-10)
	    {
	      int stop_norm = 1;
	    }
	  else if (eigen3*tnorm <  0.0)
	    {
	      eigen2 *= -1;
	      eigen3 *= -1;
	    }

	  for (size_t kr=0; kr<nearpts.size(); ++kr)
	    {
	      if (pt->pntDist(nearpts[kr]) <= radius2)
		nearpts[kr]->addCovarianceEigen(eigen1, lambda[0], eigen2, lambda[1],
						eigen3, lambda[2]);
	    }
#ifdef DEBUG_DIV      
	  if (writepoints)
	    {
	      for (int ki=0; ki<3; ++ki)
		{
		  Vector3D vec(eigenvec[ki][0], eigenvec[ki][1], eigenvec[ki][2]);
		  of << "410 1 0 4 0 0 0 255" << std::endl;
		  of << "1" << std::endl;
		  Vector3D curr = pt->getPoint();
		  of << curr << " " << curr+0.1*vec << std::endl;
		}
	    }
#endif
	  // Compute normal and curvature using LocFunc patch
	  // Point normal;//, mincvec, maxcvec;
	  // double minc, maxc;
	  // double currdist, avdist;
	  // RevEngUtils::computeLocFunc(curr, nearpts, eigen1, eigen3, normal, mincvec, minc,
	  // 				maxcvec, maxc, currdist, avdist);
	  computeLocFunc(pt, nearpts, eigen1, eigen3, radius2);
	  // Orient vectors with respect to triangulation normal
	  // The normal vectors should be OK. Curvature vectors are not necessarily
	  // consistent with regard to orientation

	  int stop_break = 1;
	}
    }

#ifdef DEBUG
  vector<Vector3D> pts1;
  vector<Vector3D> lin1;
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]); 
      if (!pt->isOutlier())
	{
	  pts1.push_back(pt->getPoint());
	  vector<ftSamplePoint*> next = pt->getNeighbours();
	  for (size_t kj=0; kj<next.size(); ++kj)
	    {
	      RevEngPoint *pt2 = dynamic_cast<RevEngPoint*>(next[kj]);
	      if (!pt2->isOutlier())
		{
		  lin1.push_back(pt->getPoint());
		  lin1.push_back(pt2->getPoint());
		}
	    }
	}
    }
  std::ofstream ofout1("not_outliers.g2");
  ofout1 << "400 1 0 4 0 0 0 255" << std::endl;
  ofout1 << pts1.size() << std::endl;
  for (size_t kj=0; kj<pts1.size(); ++kj)
    ofout1 << pts1[kj] << std::endl;
  ofout1 << "410 1 0 4 55 100 100 255" << std::endl;
  ofout1 << lin1.size()/2 << std::endl;
  for (size_t kj=0; kj<lin1.size(); kj+=2)
    ofout1 << lin1[kj] << " " << lin1[kj+1] << std::endl;
#endif
  std::sort(lambda_3.begin(), lambda_3.end());
#ifdef DEBUG_ENHANCE  
  std::cout << "No close, min: " << min_close << ", max: " << max_close << ", average: " << av_close << std::endl;
  std::cout << "lambda3, min: " << lambda_3[0] << ", max: " << lambda_3[nmbpt-1] << ", medium: " << lambda_3[nmbpt/2] << std::endl;
#endif


  for (int ki=0; ki<nmbpt; ++ki)
    {
      // Check orientation of curvature
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]); 
      if (pt->isOutlier())
	continue;
      Point tnorm = pt->getTriangNormal();
      if (tnorm.length() < 1.0e-10)
	{
	  // Fetch triangulation normal in neighbour
	  vector<ftSamplePoint*> next = pt->getNeighbours();
	  double mindist = std::numeric_limits<double>::max();
	  int ix = -1;
	  for (size_t kr=0; kr<next.size(); ++kr)
	    {
	      double dist = pt->pntDist(next[kr]);
	      RevEngPoint *nextpt = dynamic_cast<RevEngPoint*>(next[kr]);
	      Point nextnorm = nextpt->getTriangNormal();
	      if (dist < mindist && nextnorm.length() > 1.0e-10)
		{
		  mindist = dist;
		  ix = (int)kr;
		}
	    }
	  if (ix >= 0)
	    {
	      RevEngPoint *nextpt = dynamic_cast<RevEngPoint*>(next[ix]);
	      tnorm = nextpt->getTriangNormal();
	      Point PCAnorm = pt->getPCANormal();
	      if (tnorm*PCAnorm < 0.0)
		pt->turnPCA();
	      Point LocFuncnorm = pt->getLocFuncNormal();
	      if (tnorm*LocFuncnorm < 0.0)
		pt->turnLocFuncNorm();
	    }
	  
	}
    }
  

#ifdef DEBUG_ENHANCE  
  std::cout << "Start curvature filter" << std::endl;
#endif
  curvatureFilter();
  
#ifdef DEBUG_ENHANCE  
  std::cout << "Finish curvature filter" << std::endl;
#endif 
  int stop_break = 1;

}

//===========================================================================
void RevEng::computeLocFunc(RevEngPoint* pt, std::vector<RevEngPoint*>& points,
			  Point& vec1, Point& vec2, double radius)
//===========================================================================
{
  // Transform points to coordinate system given by vec1 (x-axis) and vec2 (y-axis)
  Matrix3D mat1, mat2, rotmat;
  Vector3D vec1_2(vec1[0], vec1[1], vec1[2]);
  Vector3D vec2_2(vec2[0], vec2[1], vec2[2]);
  Vector3D xaxis(1, 0, 0);
  Vector3D zaxis(0, 0, 1);
  mat1.setToRotation(vec1_2, xaxis);
  Vector3D v1 = mat1*vec1_2;
  Vector3D vec2_3 = mat1*vec2_2;
  mat2.setToRotation(vec2_3, zaxis);
  Vector3D v2 = mat2*vec2_3;
  rotmat = mat2*mat1;
  //rotmat.identity();

  // Perform rotation and sort parameter values and z-value
  int nmbpts = (int)points.size();
  vector<double> par(2*nmbpts);
  vector<double> zval(nmbpts);
  Vector3D curr = pt->getPoint();
  for (int ki=0; ki<nmbpts; ++ki)
    {
      Vector3D dv = points[ki]->getPoint() - curr;
      Vector3D dvrot = rotmat*dv;
      //Vector3D dvrot = mat2*dvrot0;
      par[2*ki] = curr[0] + dvrot[0];
      par[2*ki+1] = curr[1] + dvrot[1];
      zval[ki] = curr[2] + dvrot[2];
    }

  // Approximate z-component by biquadratic Bezier function in x and y
  int order = 3;
  shared_ptr<SplineSurface> mongesf = RevEngUtils::surfApprox(zval, 1, par, order,
							      order, order, order);

  vector<double> coefs2(3*order*order);
  std::vector<double>::iterator cf = mongesf->coefs_begin();
  for (int ka=0; ka<order; ++ka)
    {
      double vpar = mongesf->basis_v().grevilleParameter(ka);
      for (int kb=0; kb<order; ++kb, ++cf)
	{
	  double upar = mongesf->basis_u().grevilleParameter(kb);
	  coefs2[(ka*order+kb)*3] = upar;
	  coefs2[(ka*order+kb)*3+1] = vpar;
	  coefs2[(ka*order+kb)*3+2] = *cf;
	}
    }
  shared_ptr<SplineSurface> tmp(new SplineSurface(order, order, order, order, 
						  mongesf->basis_u().begin(),
						  mongesf->basis_v().begin(), &coefs2[0], 3));
#ifdef DEBUG_DIV
  int writesurface = 0;
  if (writesurface)
    {
      std::ofstream of("approx_sf.g2");
      tmp->writeStandardHeader(of);
      tmp->write(of);
      of << "400 1 0 4 0 255 0 255" << std::endl;
      of << 1 << std::endl;
      of << curr << std::endl;
      of << "400 1 0 4 255 0 0 255" << std::endl;
      of << nmbpts << std::endl;
      for (int ka=0; ka<nmbpts; ++ka)
	{
	  Point tmppt(par[2*ka], par[2*ka+1], zval[ka]);
	  of << tmppt << std::endl;
	}
    }
#endif
  
  // Compute surface normal 
  double avdist = 0.0;
  for (int ki=0; ki<nmbpts; ++ki)
    {
      Point pos;
      mongesf->point(pos, par[2*ki], par[2*ki+1]);
      avdist += fabs(zval[ki] - pos[0]);
    }
  avdist /= (double)nmbpts;

#ifdef DEBUG_MONGE
  std::ofstream of2("LocFunc_curvature.g2");
  std::ofstream of3("LocFunc_curvature2.g2");
#endif
  vector<Point> monge1, monge2, monge3, monge4;
   for (size_t kr=0; kr<points.size(); ++kr)
    {
      if (pt->pntDist(points[kr]) > radius)
	continue;

      Point triang_norm = points[kr]->getTriangNormal();
      vector<Point> der(3);
      mongesf->point(der, par[2*kr], par[2*kr+1], 1);
      Vector3D norm(-der[1][0], -der[2][0], 1.0);
      norm.normalize();

      // Accuracy of approximation
      double currdist = fabs(zval[kr] - der[0][0]);
  
      // Compute principal curvatures in curr
      SISLSurf *sislsf = GoSurf2SISL(*mongesf, false);
      int left1 = 0, left2 = 0;
      int stat = 0;
      double minc, maxc;
      double d1[2], d2[2];
      s2542(sislsf, 0, 0, 0, &par[0], &left1, &left2, &minc, &maxc, d1, d2, &stat);
      Vector3D du(1.0, 0.0, der[1][0]);
      Vector3D dv(0.0, 1.0, der[2][0]);
      Vector3D cvec1 = d1[0]*du + d1[1]*dv;
      Vector3D cvec2 = d2[0]*du + d2[1]*dv;
      if (sislsf) freeSurf(sislsf);

      // Vector3D origin(par[0], par[1], zval[0]);
      // of << "410 1 0 4 0 0 0 255" << std::endl;
      // of << "1" << std::endl;
      // of << origin << " " << origin+norm << std::endl;

      // of << "410 1 0 4 0 55 155 255" << std::endl;
      // of << "1" << std::endl;
      // of << origin << " " << origin+cvec1 << std::endl;

  
      // of << "410 1 0 4 155 55 0 255" << std::endl;
      // of << "1" << std::endl;
      // of << origin << " " << origin+cvec2 << std::endl;

  
  
      // Transform results to original coordinate system
      Matrix3D mat3, mat4, rotmat2;
      mat4.setToRotation(zaxis, vec2_3);
      mat3.setToRotation(xaxis, vec1_2);
      rotmat2 = mat3*mat4;
      //rotmat2.identity();
      //Vector3D norm0 = mat4*norm;
      Vector3D norm2 = rotmat2*norm;
      Point normal = Point(norm2[0], norm2[1], norm2[2]);
      if (triang_norm.length() > 1.0e-10 && normal*triang_norm < 0.0)
	normal *= -1;
  
      Vector3D cvec3 = rotmat2*cvec1;
      Point mincvec = Point(cvec3[0], cvec3[1], cvec3[2]); 
      Vector3D cvec4 = rotmat2*cvec2;
      Point maxcvec = Point(cvec4[0], cvec4[1], cvec4[2]);
      points[kr]->addLocFuncInfo(normal, mincvec, minc, maxcvec, maxc, currdist, avdist);

      // Vector3D xyz = points[kr]->getPoint();
      // Point xyz2 = Point(xyz[0], xyz[1], xyz[2]);
      // Vector3D der2(der[0][0], der[0][1], der[0][2]);
      // Vector3D der3 = rotmat2*der2;
      // Point der4(der3[0], der3[1], der3[2]);  // Not a 3D point!!!
      // monge1.push_back(xyz2);
      // monge1.push_back(xyz2+mincvec);
      // monge2.push_back(xyz2);
      // monge2.push_back(xyz2+maxcvec);
      // monge3.push_back(der4);
      // monge3.push_back(der4+mincvec);
      // monge4.push_back(der4);
      // monge4.push_back(der4+maxcvec);
    }

   // int writeLocFunc = 0;

   // if (writeLocFunc)
   //   {
   //     of2 << "410 1 0 4 255 0 0 255" << std::endl;
   //     of2 << monge1.size()/2 << std::endl;
   //     for (size_t kr=0; kr<monge1.size(); kr+=2)
   // 	 of2 << monge1[kr] << " " << monge1[kr+1] << std::endl;
   //     of2 << "410 1 0 4 0 0 255 255" << std::endl;
   //     of2 << monge2.size()/2 << std::endl;
   //     for (size_t kr=0; kr<monge2.size(); kr+=2)
   // 	 of2 << monge2[kr] << " " << monge2[kr+1] << std::endl;
   
   //     of3 << "400 1 0 4 0 255 0 255" << std::endl;
   //     of3 << monge3.size()/2 << std::endl;
   //     for (size_t kr=0; kr<monge3.size(); kr+=2)
   // 	 of3 << monge3[kr] << std::endl;
   //     of3 << "410 1 0 4 255 0 0 255" << std::endl;
   //     of3 << monge3.size()/2 << std::endl;
   //     for (size_t kr=0; kr<monge3.size(); kr+=2)
   // 	 of3 << monge3[kr] << " " << monge3[kr+1] << std::endl;
   //     of3 << "410 1 0 4 0 0 255 255" << std::endl;
   //     of3 << monge4.size()/2 << std::endl;
   //     for (size_t kr=0; kr<monge4.size(); kr+=2)
   // 	 of3 << monge4[kr] << " " << monge4[kr+1] << std::endl;
   //   }
  int stop_break = 1;
}
 

//===========================================================================
void RevEng::curvatureFilter()
//===========================================================================
{
  int nmbpt = tri_sf_->size();
  double radius_fac = 0.2; //0.7;
  vector<vector<RevEngPoint*> > nearpts(nmbpt);
  bool smoothcurv = true; //false;
  if (smoothcurv)
    {
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      
      // Fetch nearby points
      double local_len = pt->getMeanEdgLen();
      double radius = 0.5*radius_fac*rfac_*(local_len + mean_edge_len_);
      radius = std::min(radius, 20.0*mean_edge_len_);
      pt->fetchClosePoints2(radius, min_next_/2, max_next_/2, nearpts[ki]);
      if (nearpts[ki].size() == 0)
	continue;
      vector<double> H0(nearpts[ki].size()+1);
      vector<double> K0(nearpts[ki].size()+1);
      H0[0] = pt->meanCurvature0();
      K0[0] = pt->GaussCurvature0();
      for (size_t kj=0; kj<nearpts[ki].size(); ++kj)
	{
	  H0[kj+1] = nearpts[ki][kj]->meanCurvature0();
	  K0[kj+1] = nearpts[ki][kj]->GaussCurvature0();
	}
      std::sort(H0.begin(), H0.end());
      std::sort(K0.begin(), K0.end());
      pt->setMeanCurvature(0.5*H0[H0.size()/2] + H0[(H0.size()+1)/2]);
      pt->setGaussCurvature(0.5*K0[K0.size()/2] + K0[(K0.size()+1)/2]);
    }
  
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      pt->updateCurvature();
    }

  int nmbsmooth = (model_character_ <= MEDIUM_ROUGH) ? 2 : 5;
  for (int ka=0; ka<nmbsmooth; ++ka)
    {
      for (int ki=0; ki<nmbpt; ++ki)
	{
	  RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      
	  // Fetch nearby points
	  double Hmean = pt->meanCurvature0();
	  double Kmean = pt->GaussCurvature0();
	  for (size_t kj=0; kj<nearpts[ki].size(); ++kj)
	    {
	      Hmean += nearpts[ki][kj]->meanCurvature0();
	      Kmean += nearpts[ki][kj]->GaussCurvature0();
	    }
	  Hmean /= (double)(nearpts[ki].size()+1);
	  Kmean /= (double)(nearpts[ki].size()+1);
	  pt->setMeanCurvature(Hmean);
	  pt->setGaussCurvature(Kmean);
	}
      for (int ki=0; ki<nmbpt; ++ki)
	{
	  RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
	  pt->updateCurvature();
	}
   }
    }
}


//===========================================================================
void RevEng::edgeClassification()
//===========================================================================
{
  int nmbpts = tri_sf_->size();
  vector<Vector3D> triangcorners;
  vector<Vector3D> curvaturecorners;
  for (int ki=0; ki<nmbpts; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      Vector3D xyz = pt->getPoint();
      
      if (pt->getTriangAngle() > norm_ang_lim_)
	triangcorners.push_back(xyz);
      
       // Curvature edge classification
      double avlen = pt->getMeanEdgLen();
      double maxpc = std::max(fabs(pt->maxPrincipalCurvature()),
			      fabs(pt->minPrincipalCurvature()));
      double crvrad = 1.0/maxpc; 
      int c1_edge = (crvrad < cfac_*avlen) ? C1_EDGE : C1_NOT_EDGE;
      if (c1_edge == C1_EDGE)
	curvaturecorners.push_back(xyz);


      // Store classification in point
      pt->setEdgeClassification(c1_edge);
    }
  
  // Specify/clean edge classification
  bool closeedge = false;
  int nmbedge = 2;
  for (int ka=0; ka<nmbpts; ++ka)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ka]);
      if (pt->isEdge(edge_class_type_))
	{
	  if (pt->isolatedEdge(edge_class_type_, nmbedge, closeedge))
	    pt->setEdgeUndef();

	  else 
	    // If the angular difference between triangle normals is less
	    // then the limit, classify the point as CLOSE_EDGE.
	    pt->adjustWithTriangNorm(norm_ang_lim_);
	}
    }
  
  vector<Vector3D> edgepts;
  for (int ka=0; ka<nmbpts; ++ka)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ka]);
      if (pt->isEdge(edge_class_type_))
	edgepts.push_back(pt->getPoint());
    }


#ifdef DEBUG_ENHANCE
  std::ofstream of2("triangcorners.g2");
  of2 << "400 1 0 4 0 0 0 255" << std::endl;
  of2 << triangcorners.size() << std::endl;
  for (size_t kj=0; kj<triangcorners.size(); ++kj)
    of2 << triangcorners[kj] << std::endl;

  std::ofstream of3("curvaturecorners.g2");
  of3 << "400 1 0 4 10 10 10 255" << std::endl;
  of3 << curvaturecorners.size() << std::endl;
  for (size_t kj=0; kj<curvaturecorners.size(); ++kj)
    of3 << curvaturecorners[kj] << std::endl;
#endif

#ifdef DEBUG_EDGE0
   if (edgepts.size() > 0)
    {
      std::ofstream ofedg("edgepts.g2");
      ofedg << "400 1 0 4 255 0 0 255" << std::endl;
      ofedg << edgepts.size() << std::endl;
      for (size_t kr=0; kr<edgepts.size();++kr)
	ofedg << edgepts[kr] << std::endl;
    }
#endif
   int stop_break = 1;
 }

//===========================================================================
void RevEng::classifyPoints()
//===========================================================================
{
  // First extract obvious edges
  edgeClassification();
  
  // Fetch relevant values for all points
  vector<vector<Vector3D> > class_pts(9);
  int nmbpts = tri_sf_->size();
  for (int ki=0; ki<nmbpts; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      Vector3D xyz = pt->getPoint();

      // Curvature surface classification
      int ctype = C1_UNDEF;
      double gausscurv = pt->GaussCurvature();
      double meancurv = pt->meanCurvature();
      if (meancurv < -zero_H_)
	{
	  if (gausscurv > zero_K_)
	    {
	      ctype = C1_PEAK;
	      class_pts[0].push_back(xyz);
	    }
	  else if (gausscurv < -zero_K_)
	    {
	      ctype = C1_SRIDGE;
	      class_pts[2].push_back(xyz);
	    }
	  else
	    {
	      ctype = C1_RIDGE;
	      class_pts[1].push_back(xyz);
	    }
	}
      else if (meancurv > zero_H_)
	{
	  if (gausscurv > zero_K_)
	    {
	      ctype = C1_PIT;
	      class_pts[6].push_back(xyz);
	    }
	  else if (gausscurv < -zero_K_)
	    {
	      ctype = C1_SVALLEY;
	      class_pts[8].push_back(xyz);
	    }
	  else
	    {
	      ctype = C1_VALLEY;
	      class_pts[7].push_back(xyz);
	    }
	}
      else
	{
	  if (gausscurv > zero_K_)
	    {
	      ctype = C1_NONE;
	      class_pts[3].push_back(xyz);
	    }
	  else if (gausscurv < -zero_K_)
	    {
	      ctype = C1_MINSURF;
	      class_pts[5].push_back(xyz);
	    }
	  else
	    {
	      ctype = C1_FLAT;
	      class_pts[4].push_back(xyz);
	    }
	}



      // Store classification in point
      pt->setClassification(ctype);
   }

#ifdef DEBUG_SEG
  std::ofstream of("curvature_segments.g2");
  for (int ka=0; ka<3; ++ka)
    for (int kb=0; kb<3; ++kb)
      {
	of << "400 1 0 4 ";
	for (int kc=0; kc<3; ++kc)
	  of << colors[3*ka+kb][kc] << " ";
	of << "255" << std::endl;
	of << class_pts[3*ka+kb].size() << std::endl;
	for (size_t kr=0; kr<class_pts[3*ka+kb].size(); ++kr)
	  of << class_pts[3*ka+kb][kr] << std::endl;
      }


 #endif  
  int stop_break = 1;
}

struct
{
  bool operator()(shared_ptr<RevEngRegion> a, shared_ptr<RevEngRegion> b)
  {
    return (a->numPoints() > b->numPoints());
    //   return 1;
    // else if (a->numPoints() == b->numPoints())
    //   return 2;
    // else
    //   return 3;
  }
} sort_region;

//===========================================================================
void RevEng::setApproxTolerance()
//===========================================================================
{
  double eps = getInitApproxTol();
#ifdef DEBUG
  std::cout << "Approx tol: " << eps << std::endl;
#endif
  // std::cout << "New tolerance: " << std::endl;
  // std::cin >> eps;
  setApproxTol(eps);
}

//===========================================================================
double RevEng::getInitApproxTol()
//===========================================================================
{
  int nmbpts = tri_sf_->size();
  vector<double> pointdist;
  double maxdist = 0.0;
  double avptdist = 0.0;
  double minptdist = std::numeric_limits<double>::max();
  double maxptdist = 0.0;
  for (int ki=0; ki<nmbpts; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      if (!pt->isEdge(edge_class_type_))
	{
	  double ptdist = pt->getPointDistance();
	  pointdist.push_back(ptdist);
	  minptdist = std::min(minptdist, ptdist);
	  maxptdist = std::max(maxptdist, ptdist);
	  avptdist += ptdist;
	}
    }
  if (pointdist.size() > 0)
    avptdist /= (double)pointdist.size();
  
  std::sort(pointdist.begin(), pointdist.end());
  
  std::sort(pointdist.begin(), pointdist.end());
  double dlim = 0.93; //0.75;
  int dix = (int)(dlim*(double)pointdist.size());
  double eps = pointdist[dix];
  if (model_character_ == MEDIUM_ROUGH)
    eps *= 1.5;
  else if (model_character_ == ROUGH)
    eps *= 2.0;

  // Just to test
  //eps = 0.1;

#ifdef DEBUG_DIV
  std::cout << "Maxptdist: " << maxdist << ", avdist: " << avptdist;
  std::cout << ", medptdist: " << pointdist[pointdist.size()/2];
  std::cout << ", eps: " << eps << std::endl;
  
  std::ofstream ofd("ptdist.g2");
#endif
  double ptd1 = minptdist;
  double ptd2 = maxptdist;
  int nmbd = 12;
  double pdel = (ptd2 - ptd1)/(double)nmbd;
  vector<vector<Vector3D> > ptrs(nmbd);
  for (int ki=0; ki<nmbpts; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      Vector3D xyz = pt->getPoint();
      double dd = pt->getPointDistance();
      int ix = (int)((dd-ptd1)/pdel);
      ix = std::min(ix, nmbd-1);
      ptrs[ix].push_back(xyz);
    }

#ifdef DEBUG_DIV
  for (int ki=0; ki<nmbd; ++ki)
    {
      if (ptrs[ki].size() > 0)
	{
	  ofd << "400 1 0 4 " << colors[ki][0] << " " << colors[ki][1] << " ";
	  ofd << colors[ki][2] << " 255"  << std::endl;
	  ofd << ptrs[ki].size();
	  for (size_t kr=0; kr<ptrs[ki].size(); ++kr)
	    ofd << ptrs[ki][kr] << std::endl;
	}
    }
#endif

  return eps;
 }

//===========================================================================
void RevEng::segmentIntoRegions()
//===========================================================================
{
  // Collect continous regions
  //int min_point_in = 50; //10; //20;  // Should be set depending on the total
  // number of points. Probably class parameter. Need to look at the use
  int nmbpts = tri_sf_->size();
  for (int ki=0; ki<nmbpts; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      if (pt->hasRegion())
	continue;
      if (pt->closeEdge(edge_class_type_))
	continue;

      if (pt->nmbSameClassification(classification_type_) == 0)
	single_points_.push_back(pt);
      else
	{
	  shared_ptr<RevEngRegion> region(new RevEngRegion(classification_type_,
							   edge_class_type_));
	  pt->setRegion(region.get());
	  region->collect(pt);
	  if (region->numPoints() == 1)
	    {
	      RevEngPoint *pt_single = region->getPoint(0);
	      pt_single->unsetRegion();
	      single_points_.push_back(pt_single);
	    }
	  else
	    regions_.push_back(region);
	}
    }

  std::sort(regions_.begin(), regions_.end(), sort_region);
#ifdef DEBUG
  if (regions_.size() > 0)
    {
      std::cout << "Regions 1" << std::endl;
      std::ofstream of("regions1.g2");
      std::ofstream ofm("mid_regions1.g2");
      std::ofstream ofs("small_regions1.g2");
      writeRegionStage(of, ofm, ofs);
    }
#endif

  // Sort regions according to number of points
  std::sort(regions_.begin(), regions_.end(), sort_region);

  min_point_region_ = setSmallRegionNumber();
#ifdef DEBUG
  std::cout << "Min point region: " << min_point_region_ << std::endl;
#endif
#ifdef DEBUG_PLANAR
  std::ofstream ofpc("cand_planar.g2");
#endif  
  double lim_cone = 0.1*M_PI;
  size_t numreg = regions_.size();
  for (size_t ki=0; ki<numreg; ++ki)
    {
      double avH, avK, MAH, MAK;
      regions_[ki]->getAvCurvatureInfo(avH, avK, MAH, MAK);
      if (regions_[ki]->numPoints() > min_point_region_ &&
	  (regions_[ki]->planartype() || MAH <= zero_H_) && 
	  (!regions_[ki]->feasiblePlane(zero_H_, zero_K_)))
	{
#ifdef DEBUG_PLANAR
	  regions_[ki]->writeRegionPoints(ofpc);
#endif
	  vector<vector<RevEngPoint*> > other_groups;
	  vector<RevEngPoint*> single;
	  regions_[ki]->splitPlanar(lim_cone, min_point_region_/2, other_groups, single);
	  if (other_groups.size() > 0)
	    {
	      vector<HedgeSurface*> dummy;
	      surfaceExtractOutput((int)ki, other_groups, dummy);
	    }
	  if (single.size())
	    single_points_.insert(single_points_.end(), single.begin(), single.end());

	  int stop_cand_plane = 1;
	}
    }
  
  // Set adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

#ifdef DEBUG_PLANAR
  std::ofstream ofp("planar_reg.g2");
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      if (regions_[ki]->feasiblePlane(zero_H_, zero_K_))
	regions_[ki]->writeRegionPoints(ofp);
    }
#endif
  // Integrate single points when appropriate
  vector<RevEngPoint*> remaining_single;
  for (int ka=0; ka<(int)single_points_.size(); ++ka)
    {
      bool merged = single_points_[ka]->mergeWithAdjacent(mean_edge_len_);
      if (!merged)
	{
	  single_points_[ka]->setOutlier();
	  remaining_single.push_back(single_points_[ka]);
	}
    }
  std::swap(single_points_, remaining_single);

  // Merge adjacent planar regions
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      if (regions_[ki]->feasiblePlane(zero_H_, zero_K_))
	{
	  vector<RevEngRegion*> grown_regions;
	  vector<HedgeSurface*> adj_surfs;
	  bool merged = regions_[ki]->mergePlanarReg(zero_H_, zero_K_, 
						     approx_tol_, mainaxis_,
						     grown_regions);
	  if (merged)
	    {
	      if (grown_regions.size() > 0 || adj_surfs.size() > 0)
		updateRegionsAndSurfaces(ki, grown_regions, adj_surfs);
	      regions_[ki]->updateRegionAdjacency();
	    }
	}
    }
  
  std::sort(regions_.begin(), regions_.end(), sort_region);
  
#ifdef DEBUG
  checkConsistence("Regions1_2");

  if (regions_.size() > 0)
    {
      std::cout << "Regions 1_2" << std::endl;
      std::ofstream of("regions1_2.g2");
      std::ofstream ofm("mid_regions1_2.g2");
      std::ofstream ofs("small_regions1_2.g2");
      writeRegionStage(of, ofm, ofs);
     }
#endif
  
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }
  
  // Simplify regions structure
#ifdef DEBUG
  std::cout << "Number of regions pre integrate: " << regions_.size() << std::endl;
#endif
  int num_reg = (int)regions_.size();
  for (int ka=num_reg-1; ka>=0; --ka)
    {
      bool to_be_removed = regions_[ka]->integrateInAdjacent(mean_edge_len_,
							     min_next_, max_next_,
							     approx_tol_, 0.5,
							     max_nmb_outlier_);
      if (to_be_removed)
	{
	  regions_[ka]->removeFromAdjacent();
	  regions_[ka]->clearRegionAdjacency();
	  if (ka < num_reg-1)
	    std::swap(regions_[ka], regions_[num_reg-1]);
	  num_reg--;
	}
    }
  if (num_reg < (int)regions_.size())
    regions_.erase(regions_.begin()+num_reg, regions_.end());
#ifdef DEBUG
  std::cout << "Number of regions post integrate: " << regions_.size() << std::endl;
#endif  
  std::sort(regions_.begin(), regions_.end(), sort_region);
  
#ifdef DEBUG
  checkConsistence("Regions2");

  if (regions_.size() > 0)
    {
      std::cout << "Regions 2" << std::endl;
      std::ofstream of("regions2.g2");
      std::ofstream ofm("mid_regions2.g2");
      std::ofstream ofs("small_regions2.g2");
      writeRegionStage(of, ofm, ofs);
     }
#endif
  
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }
  
  // Update minimum number of points in region for surface generation.
  // First sort regions according to number of points
  std::sort(regions_.begin(), regions_.end(), sort_region);

  min_point_region_ = setSmallRegionNumber();
#ifdef DEBUG
  std::cout << "Min point region (2): " << min_point_region_ << std::endl;
#endif
}


//===========================================================================
void RevEng::initialSurfaces()
//===========================================================================
{
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }
  
  // Create surfaces in simple regions (planes and cylinders) and extract
  // deviant points
  int min_point_in = 50; //10; //20;  // Should be set depending on the total
  // number of points. Probably class parameter. Need to look at the use
  double angtol = 5.0*anglim_;
  //int regsize = (int)regions_.size();
  std::sort(regions_.begin(), regions_.end(), sort_region);
#ifdef DEBUG
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      RevEngRegion *first = regions_[ki]->getPoint(0)->region();
      int num = regions_[ki]->numPoints();
      for (int ka=1; ka<num; ++ka)
	if (regions_[ki]->getPoint(ka)->region() != first)
	  std::cout << "Inconsistent region pointers, pre initPlaneCyl: " << ki << " " << ka << std::endl;
    }
#endif
  for (int kr=0; kr<(int)regions_.size(); ++kr)
    {
#ifdef DEBUG
      std::ofstream of0("init_reg.g2");
      regions_[kr]->writeRegionInfo(of0);
#endif
      
      if (regions_[kr]->numPoints() < min_point_region_)
	continue;

      vector<vector<RevEngPoint*> > out_groups;
      vector<RevEngPoint*> single;
      vector<shared_ptr<HedgeSurface> > sfs;
      vector<HedgeSurface*> prev_sfs;
      bool repeat = false;
      regions_[kr]->initPlaneCyl(min_point_in, min_point_region_,
				 approx_tol_, angtol, mainaxis_,
				 zero_H_, zero_K_, sfs, out_groups, single, repeat);
      if (single.size() > 0)
	single_points_.insert(single_points_.end(), single.begin(), single.end());
      if (out_groups.size() > 0 || prev_sfs.size() > 0)
	surfaceExtractOutput((int)kr, out_groups, prev_sfs);

#ifdef DEBUG_CHECK
      bool connect = regions_[kr]->isConnected();
      if (!connect)
	std::cout << "initPlaneCyl, disconnected region " << kr << std::endl;
#endif
      if (sfs.size() > 0)
	surfaces_.insert(surfaces_.end(), sfs.begin(), sfs.end());
      if (repeat)
	--kr;
    }
  
  std::sort(regions_.begin(), regions_.end(), sort_region);

#ifdef DEBUG
  checkConsistence("Regions3");

  if (regions_.size() > 0)
    {
      std::cout << "Regions 3" << std::endl;
      std::ofstream of("regions3.g2");
      std::ofstream ofm("mid_regions3.g2");
      std::ofstream ofs("small_regions3.g2");
      writeRegionStage(of, ofm, ofs);
      std::ofstream of0("regions3_helix.g2");
      for (int ki=0; ki<(int)regions_.size(); ++ki)
	{
	  if (regions_[ki]->getSurfaceFlag() == PROBABLE_HELIX)
	    {
	      regions_[ki]->writeRegionInfo(of0);
	      if (regions_[ki]->hasSurface())
		regions_[ki]->writeSurface(of0);
	    }
	}
    }
  
  if (single_points_.size() > 0)
    {
      std::ofstream of_single("single_pts3.g2");
      of_single << "400 1 0 4 0 0 0 255" << std::endl;
      of_single << single_points_.size() << std::endl;
      for (size_t kr=0; kr<single_points_.size(); ++kr)
	of_single << single_points_[kr]->getPoint() << std::endl;
     }

  if (surfaces_.size() > 0)
    {
      std::ofstream of("regsurf3.g2");
      writeRegionWithSurf(of);
    }
#endif
  
#ifdef DEBUG
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      std::set<RevEngPoint*> tmpset(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
      if ((int)tmpset.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch, regions 3. " << ki << " " << tmpset.size() << " " << regions_[ki]->numPoints() << std::endl;
      int num = regions_[ki]->numPoints();
      for (int ka=0; ka<num; ++ka)
	if (regions_[ki]->getPoint(ka)->region() != regions_[ki].get())
	  std::cout << "Inconsistent region pointers, post initPlaneCyl: " << ki << " " << ka << std::endl;
    }
#endif

  // Integrate single points when appropriate
  vector<RevEngPoint*> remaining_single;
  for (int ka=0; ka<(int)single_points_.size(); ++ka)
    {
      if (single_points_[ka]->isOutlier())
	continue;
      bool merged = single_points_[ka]->mergeWithAdjacent(mean_edge_len_);
      if (!merged)
	{
	  single_points_[ka]->setOutlier();
	  remaining_single.push_back(single_points_[ka]);
	}
    }
  std::swap(single_points_, remaining_single);

  // Sort regions according to number of points
  std::sort(regions_.begin(), regions_.end(), sort_region);

#ifdef DEBUG
  if (surfaces_.size() > 0)
    {
      std::ofstream of("regsurf4_0.g2");
      writeRegionWithSurf(of);
    }
#endif
	  
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }
}

//===========================================================================
void RevEng::growSurfaces()
//===========================================================================
{
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }
  
  //int min_point_in = 50; //10; //20;  // Should be set depending on the total
  // number of points. Probably class parameter. Need to look at the use
  double angtol = 5.0*anglim_;

  bool joinreg = true;
  if (joinreg)
    {
#ifdef DEBUG
      std::cout << "Pre join. Number of regions: " << regions_.size() << std::endl;
#endif
      // int num_reg = (int)regions_.size();
      for (int ki=0; ki<(int)regions_.size(); ++ki)
	{
#ifdef DEBUG
      std::ofstream of0("init_join.g2");
      regions_[ki]->writeRegionInfo(of0);
#endif
	  if (regions_[ki]->numPoints() < min_point_region_)
	    continue;   // Grow into larger
	  
	  if (regions_[ki]->hasSurface())
	    {
	      growSurface(ki);
#ifdef DEBUG_CHECK
	      int num = regions_[ki]->numPoints();
	      for (int ka=0; ka<num; ++ka)
		if (regions_[ki]->getPoint(ka)->region() != regions_[ki].get())
		  std::cout << "Inconsistent region pointers, post grow: " << ki << " " << ka << std::endl;
#endif
	    }
	  // else
	  //   {
	  //     vector<RevEngRegion*> adapted_regions;
	  //     regions_[ki]->joinRegions(mainaxis_, approx_tol_,
	  // 				angtol, adapted_regions);
	  //     for (size_t kj=0; kj<adapted_regions.size(); ++kj)
	  // 	{
	  // 	  size_t kr=0;
	  // 	  for (kr=0; kr<regions_.size(); ++kr)
	  // 	    if (adapted_regions[kj] == regions_[kr].get())
	  // 	      break;

	  // 	  if (kr < regions_.size())
	  // 	    {
	  // 	      // std::swap(regions_[kr], regions_[num_reg-1]);
	  // 	      // num_reg--;
	  // 	      regions_.erase(regions_.begin()+kr);
	  // 	    }
	  // 	}
	  //   }
#ifdef DEBUG
      std::ofstream of02("post_join.g2");
      regions_[ki]->writeRegionInfo(of02);
#endif
      int stop_grow = 1;
	}
      // if (num_reg < (int)regions_.size())
      //   regions_.erase(regions_.begin()+num_reg, regions_.end());
#ifdef DEBUG
      std::cout << "Post join. Number of regions: " << regions_.size() << std::endl;

      checkConsistence("Regions4");

      if (regions_.size() > 0)
	{
	  std::cout << "Regions 4" << std::endl;
	  std::ofstream of("regions4.g2");
	  std::ofstream ofm("mid_regions4.g2");
	  std::ofstream ofs("small_regions4.g2");
	  writeRegionStage(of, ofm, ofs);
	  if (surfaces_.size() > 0)
	    {
	      std::ofstream of("regsurf4.g2");
	      writeRegionWithSurf(of);
	    }
	  std::ofstream of0("regions4_helix.g2");
	  for (int ki=0; ki<(int)regions_.size(); ++ki)
	    {
	      if (regions_[ki]->getSurfaceFlag() == PROBABLE_HELIX)
		{
		  regions_[ki]->writeRegionInfo(of0);
		  if (regions_[ki]->hasSurface())
		    regions_[ki]->writeSurface(of0);
		}
	    }
	}
#endif
    }
#ifdef DEBUG_CHECK
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      int num = regions_[ki]->numPoints();
      for (int ka=0; ka<num; ++ka)
	if (regions_[ki]->getPoint(ka)->region() != regions_[ki].get())
	  std::cout << "Inconsistent region pointers, post grow: " << ki << " " << ka << std::endl;
    }

  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      std::set<RevEngPoint*> tmpset(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
      if ((int)tmpset.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch, regions 4. " << ki << " " << tmpset.size() << " " << regions_[ki]->numPoints() << std::endl;
    }
#endif
  
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }


  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      if (regions_[ki]->hasSurface())
	{
	  vector<RevEngRegion*> grown_regions;
	  vector<HedgeSurface*> adj_surfs;
	  vector<RevEngEdge*> adj_edgs;
	  regions_[ki]->mergeAdjacentSimilar(approx_tol_, angtol,
					     grown_regions, adj_surfs, adj_edgs);
	if (grown_regions.size() > 0 || adj_surfs.size() > 0)
	  updateRegionsAndSurfaces(ki, grown_regions, adj_surfs);
	for (size_t kr=0; kr<adj_edgs.size(); ++kr)
	  {
	    size_t kj;
	    for (kj=0; kj<edges_.size(); ++kj)
	      if (edges_[kj].get() == adj_edgs[kr])
		break;
	    if (kj < edges_.size())
	      edges_.erase(edges_.begin()+kj);
	  }
	}
    }
  
#ifdef DEBUG_CHECK
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      int num = regions_[ki]->numPoints();
      for (int ka=0; ka<num; ++ka)
	if (regions_[ki]->getPoint(ka)->region() != regions_[ki].get())
	  std::cout << "Inconsistent region pointers, post mergeAdjacentSimilar: " << ki << " " << ka << std::endl;
    }

  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      std::set<RevEngPoint*> tmpset(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
      if ((int)tmpset.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch, regions 4_2. " << ki << " " << tmpset.size() << " " << regions_[ki]->numPoints() << std::endl;
    }
#endif
  
#ifdef DEBUG
      std::cout << "Post merge similar. Number of regions: " << regions_.size() << std::endl;

      checkConsistence("Regions4_2");

      if (regions_.size() > 0)
	{
	  std::cout << "Regions 4_2" << std::endl;
	  std::ofstream of("regions4_2.g2");
	  std::ofstream ofm("mid_regions4_2.g2");
	  std::ofstream ofs("small_regions4_2.g2");
	  writeRegionStage(of, ofm, ofs);
	  std::ofstream of0("regions4_2_helix.g2");
	  for (int ki=0; ki<(int)regions_.size(); ++ki)
	    {
	      if (regions_[ki]->getSurfaceFlag() == PROBABLE_HELIX)
		{
		  regions_[ki]->writeRegionInfo(of0);
		  if (regions_[ki]->hasSurface())
		    regions_[ki]->writeSurface(of0);
		}
	    }
	}
      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf4_2.g2");
	  writeRegionWithSurf(of);
	}
#endif
      
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
#ifdef DEBUG
      std::ofstream ofp("init_growplane.g2");
      regions_[ki]->writeRegionInfo(ofp);
#endif
      if (regions_[ki]->numPoints() < min_point_region_)
	continue;   // Not a stable source
	  
      if (!regions_[ki]->hasSurface())
	continue;

      vector<RevEngRegion*> grown_regions;
      vector<HedgeSurface*> adj_surfs;
      vector<vector<RevEngPoint*> > added_groups;
      regions_[ki]->growPlaneOrCyl(mainaxis_, min_point_region_, approx_tol_,
				   angtol, grown_regions, adj_surfs,
				   added_groups);
      updateRegionsAndSurfaces(ki, grown_regions, adj_surfs);
      vector<HedgeSurface*> dummy_surfs;
      if (added_groups.size() > 0)
	surfaceExtractOutput((int)ki, added_groups, dummy_surfs);
      int stop_grow = 1;
}

#ifdef DEBUG
      std::cout << "Grow plane. Number of regions: " << regions_.size() << std::endl;

      checkConsistence("Regions4_3");

      if (regions_.size() > 0)
	{
	  std::cout << "Regions 4_3" << std::endl;
	  std::ofstream of("regions4_3.g2");
	  std::ofstream ofm("mid_regions4_3.g2");
	  std::ofstream ofs("small_regions4_3.g2");
	  writeRegionStage(of, ofm, ofs);
	  std::ofstream of0("regions4_3_helix.g2");
	  for (int ki=0; ki<(int)regions_.size(); ++ki)
	    {
	      if (regions_[ki]->getSurfaceFlag() == PROBABLE_HELIX)
		{
		  regions_[ki]->writeRegionInfo(of0);
		  if (regions_[ki]->hasSurface())
		    regions_[ki]->writeSurface(of0);
		}
	    }
	}
      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf4_3.g2");
	  writeRegionWithSurf(of);
	}
#endif
      
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

  // Simplify regions structure
#ifdef DEBUG
  std::cout << "Number of regions pre integrate: " << regions_.size() << std::endl;
#endif
  int num_reg = (int)regions_.size();
  for (int ka=num_reg-1; ka>=0; --ka)
    {
      HedgeSurface* hedge = (regions_[ka]->hasSurface()) ?
	regions_[ka]->getSurface(0) : 0;
      bool to_be_removed = regions_[ka]->integrateInAdjacent(mean_edge_len_,
							     min_next_, max_next_,
							     approx_tol_, 0.5,
							     max_nmb_outlier_);
      if (to_be_removed)
	{
	  regions_[ka]->removeFromAdjacent();
	  regions_[ka]->clearRegionAdjacency();
	  if (hedge)
	    {
	      size_t kr;
	      for (kr=0; kr<surfaces_.size(); ++kr)
		if (surfaces_[kr].get() == hedge)
		  break;
	      if (kr < surfaces_.size())
		surfaces_.erase(surfaces_.begin()+kr);
	    }

	  if (ka < num_reg-1)
	    std::swap(regions_[ka], regions_[num_reg-1]);
	  num_reg--;
	}
    }
  if (num_reg < (int)regions_.size())
    regions_.erase(regions_.begin()+num_reg, regions_.end());
#ifdef DEBUG
  std::cout << "Number of regions post integrate: " << regions_.size() << std::endl;
  
  checkConsistence("Regions5");

  if (regions_.size() > 0)
    {
      std::cout << "Regions 5" << std::endl;
      std::ofstream of("regions5.g2");
      std::ofstream ofm("mid_regions5.g2");
      std::ofstream ofs("small_regions5.g2");
      writeRegionStage(of, ofm, ofs);
      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf5.g2");
	  writeRegionWithSurf(of);
	}
     }
#endif
  
#ifdef DEBUG
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      std::set<RevEngPoint*> tmpset(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
      if ((int)tmpset.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch, regions 5. " << ki << " " << tmpset.size() << " " << regions_[ki]->numPoints() << std::endl;
    }
#endif
  std::sort(regions_.begin(), regions_.end(), sort_region);
  
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

  
  int stop_break2 = 1;
}


//===========================================================================
void RevEng::updateAxesAndSurfaces()
//===========================================================================
{
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }
  
  // // First extract low accuracy points from groups with surface
  double angtol = 5.0*anglim_;
  // for (size_t ki=0; ki<regions_.size(); ++ki)
  //   {
  //     if (!regions_[ki]->hasSurface())
  // 	continue;

  //     vector<vector<RevEngPoint*> > added_groups;
  //     vector<HedgeSurface*> dummy_surfs;
  //     regions_[ki]->removeLowAccuracyPoints(min_point_region_, 
  // 					    approx_tol_, angtol, added_groups);
  //     if (added_groups.size() > 0)
  // 	surfaceExtractOutput((int)ki, added_groups, dummy_surfs);
  //   }

  std::sort(regions_.begin(), regions_.end(), sort_region);
  
  vector<int> reg_size(surfaces_.size());
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    reg_size[ki] = surfaces_[ki]->numPoints();

  std::sort(reg_size.begin(), reg_size.end());

  Point axis[3];
  int min_num = reg_size[(int)reg_size.size()/4];
  min_num = std::min(min_num, reg_size[reg_size.size()-1]/10);
  min_num = std::max(min_num, reg_size[reg_size.size()-1]/100);
  double max_ang = 0.1*M_PI;

  Point plane_axis[3];
  int num_pts1[3];
  computeAxisFromPlane(mainaxis_, min_num, max_ang, plane_axis, num_pts1);

  Point cyl_axis[3];
  int num_pts2[3];
  computeAxisFromCylinder(plane_axis, min_num, max_ang, cyl_axis, num_pts2);

  // Update main axes. Prioritize information from planes
  for (int ka=0; ka<3; ++ka)
    {
      num_pts1[ka] *= 2;
      int all_pts = num_pts1[ka] + num_pts2[ka];
      if (all_pts == 0)
	continue;
      double fac1 = (double)num_pts1[ka]/(double)(all_pts);
      double fac2 = (double)num_pts2[ka]/(double)(all_pts);
      if (cyl_axis[ka]*plane_axis[ka] < 0.0)
	cyl_axis[ka] *= -1.0;
      mainaxis_[ka] = fac1*plane_axis[ka] + fac2*cyl_axis[ka];
      mainaxis_[ka].normalize();
    }

  // Ensure orthogonality
  for (int ka=0; ka<3; ++ka)
    for (int kb=ka+1; kb<3; ++kb)
      if (num_pts1[ka] + num_pts2[ka] < num_pts1[kb] + num_pts2[kb])
	{
	  std::swap(num_pts1[ka], num_pts1[kb]);
	  std::swap(num_pts2[ka], num_pts2[kb]);
	  std::swap(mainaxis_[ka], mainaxis_[kb]);
	}

  Point tmp_axis[3];
  tmp_axis[2] = mainaxis_[0].cross(mainaxis_[1]);
  tmp_axis[1] = (num_pts1[2]+num_pts2[2] > 0) ?
	    mainaxis_[2].cross(mainaxis_[0]) : mainaxis_[1];
  tmp_axis[0] = (num_pts1[1]+num_pts2[1] > 0 && num_pts1[2]+num_pts2[2] > 0) ?
		 mainaxis_[1].cross(mainaxis_[2]) : mainaxis_[0];
  for (int ka=0; ka<3; ++ka)
    {
      if (tmp_axis[ka]*mainaxis_[ka] < 0.0)
	tmp_axis[ka] *= -1.0;
      if (num_pts1[ka] + num_pts2[ka] > 0)
	tmp_axis[ka] = 0.5*(tmp_axis[ka] + mainaxis_[ka]);
      tmp_axis[ka].normalize();
    }

  if (num_pts1[0]+num_pts2[0] > 0 && num_pts1[1]+num_pts2[1] > 0 &&
      num_pts1[2]+num_pts2[2] > 0)
    {
      mainaxis_[0] = tmp_axis[1].cross(tmp_axis[2]);
      mainaxis_[1] = tmp_axis[2].cross(tmp_axis[0]);
      mainaxis_[2] = tmp_axis[0].cross(tmp_axis[1]);
      for (int ka=0; ka<3; ++ka)
	{
	  if (mainaxis_[ka]*tmp_axis[ka] < 0.0)
	    mainaxis_[ka] *= -1.0;
	  mainaxis_[ka] = 0.5*(tmp_axis[ka] + mainaxis_[ka]);
	  mainaxis_[ka].normalize();
	}
    }
  else
    {
      for (int ka=0; ka<3; ++ka)
	{
	  mainaxis_[ka] = tmp_axis[ka];
	}
    }
  // Point tmp_axis = mainaxis_[0].cross(mainaxis_[1]);
  // tmp_axis.normalize();
  // if (mainaxis_[2]*tmp_axis < 0.0)
  //   tmp_axis *= -1;
  // mainaxis_[2] = 0.5*(mainaxis_[2] + tmp_axis);
  // tmp_axis = mainaxis_[2].cross(mainaxis_[0]);
  // tmp_axis.normalize();
  // if (mainaxis_[1]*tmp_axis < 0.0)
  //   tmp_axis *= -1;
  // mainaxis_[1] = 0.5*(mainaxis_[1] + tmp_axis);
  // tmp_axis = mainaxis_[1].cross(mainaxis_[2]);
  // tmp_axis.normalize();
  // if (mainaxis_[0]*tmp_axis < 0.0)
  //   tmp_axis *= -1;
  // mainaxis_[0] = 0.5*(mainaxis_[0] + tmp_axis);

  mainaxis_[2] = mainaxis_[0].cross(mainaxis_[1]);
  mainaxis_[1] = mainaxis_[2].cross(mainaxis_[0]);
  for (int ka=0; ka<3; ++ka)
    mainaxis_[ka].normalize();

  // Update surfaces
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      if (regions_[ki]->hasSurface())// && (!regions_[ki]->hasRevEdges()))
	{
	  /*bool updated =*/ (void)axisUpdate(ki, max_ang, angtol);
	  // if (updated)
	  //   {
	  //     growSurface(ki);
	  //     int stop_break0 = 1;
	  //   }
	}
    }
  
#ifdef DEBUG
  std::ofstream of0("regions5_2_0_helix.g2");
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      if (regions_[ki]->getSurfaceFlag() == PROBABLE_HELIX)
	{
	  regions_[ki]->writeRegionInfo(of0);
	  if (regions_[ki]->hasSurface())
	    regions_[ki]->writeSurface(of0);
	}
    }
#endif
	  
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      if (regions_[ki]->hasSurface())
	{
	  vector<RevEngRegion*> grown_regions;
	  vector<HedgeSurface*> adj_surfs;
	  vector<RevEngEdge*> adj_edgs;
	  regions_[ki]->mergeAdjacentSimilar(approx_tol_, angtol,
					     grown_regions, adj_surfs, adj_edgs);
	if (grown_regions.size() > 0 || adj_surfs.size() > 0)
	  {
	    updateRegionsAndSurfaces(ki, grown_regions, adj_surfs);
	    for (size_t kr=0; kr<adj_edgs.size(); ++kr)
	      {
		size_t kj;
		for (kj=0; kj<edges_.size(); ++kj)
		  if (edges_[kj].get() == adj_edgs[kr])
		    break;
		if (kj < edges_.size())
		  edges_.erase(edges_.begin()+kj);
	      }
	    // if (!regions_[ki]->hasRevEdges())
	    //   {
		bool updated = axisUpdate(ki, max_ang, angtol);
		if (!updated)
		  int stop_break1 = 1;
	      // }
	  }
	}
    }
#ifdef DEBUG
  std::cout << "Post merge similar. Number of regions: " << regions_.size() << std::endl;

  checkConsistence("Regions5_2");

  if (regions_.size() > 0)
    {
      std::cout << "Regions 5_2" << std::endl;
      std::ofstream of("regions5_2.g2");
      std::ofstream ofm("mid_regions5_2.g2");
      std::ofstream ofs("small_regions5_2.g2");
      writeRegionStage(of, ofm, ofs);
      std::ofstream of0("regions5_2_helix.g2");
      for (int ki=0; ki<(int)regions_.size(); ++ki)
	{
	  if (regions_[ki]->getSurfaceFlag() == PROBABLE_HELIX)
	    {
	      regions_[ki]->writeRegionInfo(of0);
	      if (regions_[ki]->hasSurface())
		regions_[ki]->writeSurface(of0);
	    }
	}
      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf5_2.g2");
	  writeRegionWithSurf(of);
	}
    }
#endif
#ifdef DEBUG
  std::ofstream ofu1("unresolved.g2");
  for (size_t kr=0; kr<regions_.size(); ++kr)
    {
      if (regions_[kr]->hasSurface())
	continue;
      if (regions_[kr]->hasAssociatedBlend())
	continue;
      regions_[kr]->writeRegionPoints(ofu1);
    }
#endif
   int stop_break = 1;
}


//===========================================================================
bool RevEng::axisUpdate(int ix, double max_ang, double angtol)
//===========================================================================
{
  if (!regions_[ix]->hasSurface())
    return false;

  HedgeSurface *hsurf = regions_[ix]->getSurface(0);
  shared_ptr<ElementarySurface> elem =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(hsurf->surface());
  if (!elem.get())
    return false;

#ifdef DEBUG_AXIS
  std::ofstream of("axis_adapt.g2");
  elem->writeStandardHeader(of);
  elem->write(of);
  regions_[ix]->writeRegionPoints(of);
#endif
  vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*>  > adj_elem, adj_elem_base;
  regions_[ix]->getAdjacentElemInfo(adj_elem, adj_elem_base);
  Point adj_axis, adj_pos;
  double adj_ang = M_PI;
  Point vec = elem->direction();
  double pihalf = 0.5*M_PI;
  double min_perpen = M_PI;
  for (size_t ki=0; ki<adj_elem.size(); ++ki)
    {
      if (adj_elem[ki].first->instanceType() != Class_Plane &&
	  adj_elem[ki].first->instanceType() != Class_Cylinder)
	continue;
      if (adj_elem[ki].second->hasBlendEdge())
	continue; // Derived information
      int sfflag = adj_elem[ki].second->getSurfaceFlag();
      double anglim = (sfflag == ACCURACY_OK) ? 2.0*max_ang : max_ang;
      Point dir = adj_elem[ki].first->direction();
      double ang = vec.angle(dir);
      if (fabs(pihalf - ang) < min_perpen)
	min_perpen = fabs(pihalf - ang);
      ang = std::min(ang, M_PI-ang);
      if (ang < anglim && ang < adj_ang)
	{
	  adj_ang = ang;
	  adj_axis = dir;
	  if (adj_elem[ki].first->instanceType() == Class_Cylinder)
	    adj_pos = adj_elem[ki].first->location();
	}
    }

  int ka_min = -1;
  double ang_main = M_PI;
  for (int ka=0; ka<3; ++ka)
    {
      double ang = vec.angle(mainaxis_[ka]);
      ang = std::min(ang, M_PI-ang);
      if (ang < max_ang && ang < ang_main)
	{
	  ang_main = ang;
	  ka_min = ka;
	}
    }
  
  bool updated = false;
  if (adj_axis.dimension() != 0 || ka_min >= 0) 
    updated = regions_[ix]->updateSurfaceWithAxis(min_point_region_, adj_axis, mainaxis_,
						  ka_min, approx_tol_, angtol, adj_pos);
  // if (updated == false && elem->instanceType() == Class_Plane && min_perpen < anglim_)
  //   {
  //     regions_[ix]->setPlaneParam(min_point_region_, mainaxis_, approx_tol_, angtol);
  //   }
  
  return updated;
}

//===========================================================================
void RevEng::recognizeEdges(bool only_curve)
//===========================================================================
{
 // Ensure some limitation of surface size
  Point low = bbox_.low();
  Point high = bbox_.high();
  double diag = low.dist(high);
  double blendfac = 2.0;

  // Ensure bounded surfaces
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    surfaces_[ki]->limitSurf(diag);

  double angtol = 5.0*anglim_;
  double pihalf = 0.5*M_PI;
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      if (!regions_[ki]->hasSurface())
	continue;
      if (regions_[ki]->hasAssociatedBlend())
	continue;
      if (regions_[ki]->hasBlendEdge())
	continue;

      vector<RevEngEdge*> rev_edgs1 = regions_[ki]->getAllRevEdges();
      
      int code;
      int classtype = regions_[ki]->getSurface(0)->instanceType(code);
      if (classtype != Class_Plane && classtype != Class_Cylinder &&
	  classtype != Class_Cone && classtype != Class_Sphere)
	continue;  // Preliminary

      shared_ptr<ParamSurface> surf1 = regions_[ki]->getSurface(0)->surface();
      shared_ptr<ElementarySurface> elem1 =
	dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf1);
      Point dir1 = elem1->direction();
      if (elem1->instanceType() == Class_Plane &&
	  dir1*regions_[ki]->getPoint(0)->getTriangNormal() < 0.0)
	dir1 *= -1;
      for (size_t kj=ki+1; kj<regions_.size(); ++kj)
	{
	  bool only_curve2 = only_curve;
	  if (!regions_[kj]->hasSurface())
	    continue;
	  if (regions_[ki]->hasAssociatedBlend())
	    continue;
	  if (regions_[ki]->hasBlendEdge())
	    continue;
	  
	  vector<RevEngEdge*> rev_edgs2 = regions_[kj]->getAllRevEdges();
	  if (rev_edgs1.size() > 0 && rev_edgs2.size() > 0)
	    {
	      size_t kr, kh;
	      for (kr=0; kr<rev_edgs1.size(); ++kr)
		{
		  for (kh=0; kh<rev_edgs2.size(); ++kh)
		    if (rev_edgs1[kr] == rev_edgs2[kh])
		      break;
		  if (kh < rev_edgs2.size())
		    break;
		}
	      if (kr < rev_edgs1.size())
		continue;
	    }
	  
 	  int code2;
	  int classtype2 = regions_[kj]->getSurface(0)->instanceType(code2);
	  if (classtype2 != Class_Plane && classtype2 != Class_Cylinder &&
	      classtype2 != Class_Cone && classtype2 != Class_Sphere)
	    continue;  // Preliminary
	  if (classtype == Class_Sphere || classtype2 == Class_Sphere)
	    only_curve2 = true;
	  if (classtype == Class_Cylinder && classtype2 == Class_Cylinder)
	    continue;
	  if (classtype == Class_Cone && classtype2 == Class_Cone)
	    continue;
	  if ((classtype == Class_Sphere && classtype2 != Class_Plane) ||
	      (classtype2 == Class_Sphere && classtype != Class_Plane))
	    continue;

#ifdef DEBUG_EDGE
	  std::ofstream of1("adj_regs.g2");
	  regions_[ki]->writeRegionPoints(of1);
	  regions_[ki]->writeSurface(of1);
	  regions_[kj]->writeRegionPoints(of1);
	  regions_[kj]->writeSurface(of1);
#endif
	  shared_ptr<ParamSurface> surf2 = regions_[kj]->getSurface(0)->surface();
	  shared_ptr<ElementarySurface> elem2 =
	    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf2);
	  if (regions_[ki]->isAdjacent(regions_[kj].get())
	      || (only_curve2 == false &&
		  regions_[ki]->isNextToAdjacent(regions_[kj].get())))
	    {
	      Point dir2 = elem2->direction();
	      if (elem2->instanceType() == Class_Plane &&
		  dir2*regions_[kj]->getPoint(0)->getTriangNormal() < 0.0)
		dir2 *= -1;
	      double ang = dir1.angle(dir2);
	      ang = std::min(ang, M_PI-ang);
	      bool compute_edge = false;
	      if (classtype == Class_Plane && classtype2 == Class_Plane)
		{
		  if (ang > blendfac*angtol)
		    compute_edge = true;
		}
	      else if (classtype == Class_Plane || classtype2 == Class_Plane)
		{
		  if (ang < blendfac*angtol ||
		      (only_curve2 && fabs(pihalf-ang) > blendfac*angtol))
		    compute_edge = true;
		  else if (fabs(pihalf-ang) < blendfac*angtol)
		    {
		      // Check for near tangential cases
		      Point norm = (classtype == Class_Plane) ? elem1->direction() :
			elem2->direction();
		      double dlen = fabs((elem1->location() - elem2->location())*norm);
		      double rad = (classtype == Class_Plane) ? elem2->radius(0.0, 0.0) :
			elem1->radius(0.0, 0.0);   // Not appropriate for a cone
		      if (dlen > rad + blendfac*approx_tol_ ||
			  dlen < rad - blendfac*approx_tol_)
			compute_edge = true;
		    }
		}
	      else
		{
		  // One cone and one cylinder. Check for same axis and
		  // a significant cone angle
		  Point loc1 = elem1->location();
		  Point loc2 = elem2->location();
		  Point tmp = loc1 + ((loc2-loc1)*dir1)*dir1;
		  if (ang >= blendfac*angtol || loc2.dist(tmp) > approx_tol_)
		    compute_edge = false;
		  else
		    {
		      double phi = 0.0;
		      shared_ptr<Cone> cone1 =
			dynamic_pointer_cast<Cone,ElementarySurface>(elem1);
		      shared_ptr<Cone> cone2 =
			dynamic_pointer_cast<Cone,ElementarySurface>(elem2);
		      if (cone1.get())
			phi = cone1->getConeAngle();
		      else if (cone2.get())
			phi = cone2->getConeAngle();
		      double conefac = 4.0;
		      if (fabs(phi) > conefac*angtol)
			compute_edge = true;
		      else
			compute_edge = false;
		    }
		}
	      

	      if (compute_edge)
		{
		  // Make sure that cone domains do not cover the apex
		  int ka;
		  shared_ptr<ElementarySurface> elem;
		  shared_ptr<RevEngRegion> reg;
		  for (ka=0, elem=elem1, reg=regions_[ki]; ka<2;
		       ++ka, elem=elem2, reg=regions_[kj])
		    {
		      if (elem->instanceType() == Class_Cone)
			{
			  shared_ptr<Cone> cone =
			    dynamic_pointer_cast<Cone,ElementarySurface>(elem);
			  double apar;
			  int adir;
			  cone->getDegenerateParam(apar, adir);
			  if (adir > 0)
			    {
			      RectDomain dom = cone->getParameterBounds();
			      double dom2[4];
			      dom2[0] = dom.umin();
			      dom2[1] = dom.umax();
			      dom2[2] = dom.vmin();
			      dom2[3] = dom.vmax();
			      double dom3[4];
			      reg->getDomain(dom3);
			      adir--;
			      // Assumes that the relevant part of the cone
			      // does not cover the apex
			      double midp = 0.5*(dom3[2*adir]+dom3[2*adir+1]);
			      double del = 0.01*(dom3[2*adir+1]-dom3[2*adir]);
			      if (midp < apar)
				dom2[2*adir+1] = std::min(dom2[2*adir+1], apar-del);
			      else
				dom2[2*adir] = std::max(dom2[2*adir], apar+del);
			      cone->setParameterBounds(dom2[0], dom2[2],
							dom2[1], dom2[3]);
			    }
			}
		    }

		  vector<shared_ptr<RevEngEdge> > edges =
		    defineEdgesBetween(ki, elem1, dir1,kj, elem2, dir2,
				       only_curve2);
		  if (edges.size() > 0)
		    edges_.insert(edges_.end(), edges.begin(),edges.end());
		}
	    }
	}
    }
}

//===========================================================================
void RevEng::firstEdges()
//===========================================================================
{
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

  recognizeEdges();

#ifdef DEBUG
  checkConsistence("Regions6");

   if (regions_.size() > 0)
    {
      std::cout << "Regions6" << std::endl;
      std::ofstream of("regions6.g2");
      std::ofstream ofm("mid_regions6.g2");
      std::ofstream ofs("small_regions6.g2");
      writeRegionStage(of, ofm, ofs);
     }

   if (edges_.size() > 0)
     {
       std::ofstream ofe("edges6.g2");
       writeEdgeStage(ofe);
     }
  
#endif

  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

   
#ifdef DEBUG
   std::cout << "Extend blend region collection" << std::endl;
#endif
   for (size_t ki=0; ki<edges_.size(); ++ki)
     {
       extendBlendAssociation(ki);
     }
   
 #ifdef DEBUG
  checkConsistence("Regions6_1");
#endif
  
  // Just in case composed regions have been segmented
  int min_num_point = min_point_region_/10;
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      if (regions_[ki]->hasSurface())
	continue;
      if (regions_[ki]->hasAssociatedBlend())
	continue;
      if (regions_[ki]->numPoints() < min_num_point)
	continue;
#ifdef DEBUG_EDGE
      std::ofstream of("planar_merge_cand.g2");
      regions_[ki]->writeRegionInfo(of);
#endif
      if (regions_[ki]->feasiblePlane(zero_H_, zero_K_))
	{
	  vector<RevEngRegion*> grown_regions;
	  vector<HedgeSurface*> adj_surfs;
	  bool merged = regions_[ki]->mergePlanarReg(zero_H_, zero_K_, 
						     approx_tol_, mainaxis_,
						     grown_regions);
	  if (merged)
	    {
	      if (grown_regions.size() > 0 || adj_surfs.size() > 0)
		updateRegionsAndSurfaces(ki, grown_regions, adj_surfs);
	      regions_[ki]->updateRegionAdjacency();
	    }
	}
    }

  
 #ifdef DEBUG
  checkConsistence("Regions6_2");

   if (regions_.size() > 0)
    {
      std::cout << "Regions6_2" << std::endl;
      std::ofstream of("regions6_2.g2");
      std::ofstream ofm("mid_regions6_2.g2");
      std::ofstream ofs("small_regions6_2.g2");
      writeRegionStage(of, ofm, ofs);
     }

   if (edges_.size() > 0)
     {
       std::ofstream ofe("edges6_2.g2");
       writeEdgeStage(ofe);
     }
  
#endif
   int stop_break = 1;
}

//===========================================================================
void RevEng::extendBlendAssociation(size_t ix)
//===========================================================================
{
  // Regions in blend area
  vector<RevEngRegion*> blend_regs;
  edges_[ix]->getAllBlendRegs(blend_regs);
  size_t num_blend_regs = blend_regs.size();

#ifdef DEBUG_BLEND
  std::ofstream of1("blend_regs.g2");
  for (size_t kj=0; kj<blend_regs.size(); ++kj)
    blend_regs[kj]->writeRegionPoints(of1);
#endif
  // Intersection curve
  vector<shared_ptr<CurveOnSurface> > cvs;
  edges_[ix]->getCurve(cvs);
#ifdef DEBUG_BLEND
  for (size_t kj=0; kj<cvs.size(); ++kj)
    {
      cvs[kj]->spaceCurve()->writeStandardHeader(of1);
      cvs[kj]->spaceCurve()->write(of1);
    }
  RevEngRegion *adj1, *adj2;
  edges_[ix]->getAdjacent(adj1, adj2);
  adj1->writeRegionPoints(of1);
  adj2->writeRegionPoints(of1);
#endif
  
  // Width
  double width = edges_[ix]->getDistance();

  for (size_t ki=0; ki<blend_regs.size(); ++ki)
    {
      vector<RevEngRegion*> new_blends;
      blend_regs[ki]->neighbourBlends(cvs, 2.0*width, approx_tol_, new_blends);
#ifdef DEBUG_BLEND
      std::ofstream of2("blend_regs2.g2");
      for (size_t kj=0; kj<new_blends.size(); ++kj)
	new_blends[kj]->writeRegionPoints(of2);
#endif
      for (size_t kj=0; kj<new_blends.size(); ++kj)
	{
	  size_t kr;
	  for (kr=0; kr<blend_regs.size(); ++kr)
	    if (blend_regs[kr] == new_blends[kj])
	      break;
	  if (kr == blend_regs.size())
	    blend_regs.push_back(new_blends[kj]);
	}
      int stop_break = 1;
    }

  for (size_t kj=num_blend_regs; kj<blend_regs.size(); ++kj)
    {
      edges_[ix]->addBlendRegion(blend_regs[kj]);
      blend_regs[kj]->setAssociatedBlend(edges_[ix].get());
    }
  int stop_break2 = 1;
}

//===========================================================================
bool RevEng::setBlendEdge(size_t ix)
//===========================================================================
{
 Point low = bbox_.low();
  Point high = bbox_.high();
  double diag = low.dist(high);
  double blendfac = 2.0;

  //double pihalf = 0.5*M_PI;
  double angtol = 5.0*anglim_;
  vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> > adj_elem, adj_elem_base;
  regions_[ix]->getAdjacentElemInfo(adj_elem, adj_elem_base);
  for (size_t ki=0; ki<adj_elem.size(); ++ki)
    {
      ClassType type1 = adj_elem[ki].first->instanceType();
      if (type1 != Class_Plane && type1 != Class_Cylinder)
	continue;
      Point dir1 = adj_elem[ki].first->direction();
      if (type1 == Class_Plane &&
	  dir1*adj_elem[ki].second->getMeanNormalTriang() < 0.0)
	dir1 *= -1;
      for (size_t kj=ki+1; kj<adj_elem.size(); ++kj)
	{
	  ClassType type2 = adj_elem[kj].first->instanceType();
	  if (type2 != Class_Plane && type2 != Class_Cylinder)
	    continue;
	  if (type1 == Class_Cylinder && type2 == Class_Cylinder)
	    continue;
	  Point dir2 = adj_elem[kj].first->direction();
	  if (type2 == Class_Plane &&
	      dir2*adj_elem[kj].second->getMeanNormalTriang() < 0.0)
	    dir2 *= -1;
	  double ang = dir1.angle(dir2);
	  ang = std::min(ang, M_PI-ang);
	  bool compute_edge = false;
#ifdef DEBUG_BLEND
	  std::ofstream of("adj_reg.g2");
	  adj_elem[ki].second->writeRegionPoints(of);
	  adj_elem[ki].second->writeSurface(of);
	  adj_elem[kj].second->writeRegionPoints(of);
	  adj_elem[kj].second->writeSurface(of);
#endif
	  if (type1 == Class_Plane && type2 == Class_Plane)
	    {
	      if (ang > blendfac*angtol)
		compute_edge = true;
	    }
	  else
	    {
	      if (ang < blendfac*angtol)
		compute_edge = true;
	    }

	  if (compute_edge)
	    {
	      // Check if the two regions already has a common edge
	      vector<RevEngEdge*> edges1 = adj_elem[ki].second->getAllRevEdges();
	      vector<RevEngEdge*> edges2 = adj_elem[kj].second->getAllRevEdges();
	      size_t kr, kh;
	      for (kr=0; kr<edges1.size(); ++kr)
		{
		  for (kh=0; kh<edges2.size(); ++kh)
		    if (edges1[kr] == edges2[kh])
		      break;
		  if (kh < edges2.size())
		    break;
		}

	      if (kr < edges1.size())
		{
		  // An edge exist already. Extend
		  continue;  // For the time being
		}
	      else
		{
		  // Define new edge
		  size_t ix1, ix2;
		  for (ix1=0; ix1<regions_.size(); ++ix1)
		    if (regions_[ix1].get() == adj_elem[ki].second)
		      break;
		  for (ix2=0; ix2<regions_.size(); ++ix2)
		    if (regions_[ix2].get() == adj_elem[kj].second)
		      break;
		  if (ix1 == regions_.size() || ix2 == regions_.size())
		    continue;

		  // Make sure that the adjacent surfaces is bounded
		  if (!adj_elem[ki].first->isBounded())
		    regions_[ix1]->getSurface(0)->limitSurf(diag);
		  if (!adj_elem[kj].first->isBounded())
		    regions_[ix2]->getSurface(0)->limitSurf(diag);
		  vector<shared_ptr<RevEngEdge> > edges =
		    defineEdgesBetween(ix1, adj_elem[ki].first,
				       dir1, ix2, adj_elem[kj].first, dir2);
		  if (edges.size() > 0)
		    {
		      edges_.insert(edges_.end(), edges.begin(), edges.end());
		      return true;
		    }
		}
	    }
	}
    }
  return false;
}

//===========================================================================
vector<shared_ptr<RevEngEdge> >
RevEng::defineEdgesBetween(size_t ix1, shared_ptr<ElementarySurface>& surf1,
			   Point& dir1, size_t ix2,
			   shared_ptr<ElementarySurface>& surf2, Point& dir2,
			   bool only_curve, double lenlim0, bool check_common)
//===========================================================================
{
  vector<shared_ptr<RevEngEdge> > edges;
  if (regions_[ix1]->hasAssociatedBlend() || regions_[ix2]->hasAssociatedBlend())
    return edges;
  if (((regions_[ix1]->getSurfaceFlag() == ACCURACY_POOR  ||
	regions_[ix1]->getSurfaceFlag() == NOT_SET) &&
       (!regions_[ix1]->hasBlendEdge())) ||
      ((regions_[ix2]->getSurfaceFlag() == ACCURACY_POOR  ||
	regions_[ix2]->getSurfaceFlag() == NOT_SET) &&
       (!regions_[ix2]->hasBlendEdge())))
    return edges;

  //double angtol = 5.0*anglim_;
  double diag = bbox_.low().dist(bbox_.high());
  double blendlim = std::min(0.1*diag, 30.0*mean_edge_len_);
  double maxwidth = std::max(blendlim, 0.1*diag);
  double lenlim = 10.0*mean_edge_len_; //blendlim;
  if (lenlim0 > 0.0)
    lenlim= std::min(lenlim, lenlim0);
  double int_tol = 1.0e-6;
  //bool plane1 = (surf1->instanceType() == Class_Plane);
  //bool plane2 = (surf2->instanceType() == Class_Plane);
  bool adjacent = regions_[ix1]->isAdjacent(regions_[ix2].get());

  shared_ptr<BoundedSurface> bd1, bd2;
  vector<shared_ptr<CurveOnSurface> > int_cvs1, int_cvs2;
  BoundedUtils::getSurfaceIntersections(surf1, surf2, int_tol_,
					int_cvs1, bd1,
					int_cvs2, bd2);
  if (int_cvs1.size() == 0)
    return edges;

  bool keep_length = (int_cvs1.size() == 1 && (regions_[ix1]->hasBlendEdge() ||
						regions_[ix2]->hasBlendEdge()));
		  
  // Limit intersection curves to relevant intervals
  if (!keep_length)
    {
      for (size_t kr=0; kr<int_cvs1.size(); )
	{
	  vector<shared_ptr<CurveOnSurface> > tmp_int1, tmp_int2;
	  tmp_int1.push_back(int_cvs1[kr]);
	  tmp_int2.push_back(int_cvs2[kr]);
	  vector<pair<double,double> > t1_t2, t3_t4;
	  bool OK1 =
	    regions_[ix1]->getCurveRestriction(tmp_int1, approx_tol_,
					       anglim_, t1_t2);
	  bool OK2 =
	    regions_[ix2]->getCurveRestriction(tmp_int2, approx_tol_,
					       anglim_, t3_t4);
      
	  double t1 = std::max(t1_t2[0].first, t3_t4[0].first);
	  double t2 = std::min(t1_t2[0].second, t3_t4[0].second);
	  if (t2 > t1 && (t1 > int_cvs1[kr]->startparam() ||
			  t2 < int_cvs1[kr]->endparam()))
	    {
	      double pmin = std::max(t1, int_cvs1[kr]->startparam());
	      double pmax = std::min(t2, int_cvs1[kr]->endparam());
	      shared_ptr<CurveOnSurface> sub1(int_cvs1[kr]->subCurve(pmin,pmax));
	      int_cvs1[kr] = sub1;
	      shared_ptr<CurveOnSurface> sub2(int_cvs2[kr]->subCurve(pmin,pmax));
	      int_cvs2[kr] = sub2;
	    }

	  if (t2 <= t1 || int_cvs1[kr]->estimatedCurveLength() < lenlim)
	    {
	      int_cvs1.erase(int_cvs1.begin()+kr);
	      int_cvs2.erase(int_cvs2.begin()+kr);
	    }
	  else
	    ++kr;
	}
      if (int_cvs1.size() == 0)
	return edges;
    }
  
  vector<RevEngRegion*> common_reg =
    regions_[ix1]->commonAdjacent(regions_[ix2].get());
  for (size_t kj=0; kj<common_reg.size(); )
    {
      if (common_reg[kj]->hasAssociatedBlend() || common_reg[kj]->hasBlendEdge())
  	common_reg.erase(common_reg.begin()+kj);
      else
  	++kj;
    }

  if (!adjacent && check_common)
    {
      // Check if any of the regions between is more significant
      int num1 = regions_[ix1]->numPoints();
      int num2 = regions_[ix2]->numPoints();
      int fac = 2;
      for (size_t kj=0; kj<common_reg.size(); ++kj)
	{
	  if (!common_reg[kj]->hasSurface())
	    continue;
	  int num3 = common_reg[kj]->numPoints();
	  if (num3 > fac*std::min(num1, num2))
	    return edges;  // Could be a small passage. 
	}
    }
		  
#ifdef DEBUG_EDGE
  std::ofstream of1("adj_regs_cv.g2");
  for (size_t kr=0; kr<int_cvs1.size(); ++kr)
    {
      shared_ptr<ParamCurve> cv = int_cvs1[kr]->spaceCurve();
      cv->writeStandardHeader(of1);
      cv->write(of1);
    }

  if (common_reg.size() > 0)
    {
      std::ofstream of2("between_regs.g2");
      for (size_t kr=0; kr<common_reg.size(); ++kr)
	{
	  common_reg[kr]->writeRegionPoints(of2);
	  if (common_reg[kr]->hasSurface())
	    common_reg[kr]->writeSurface(of2);
	}
    }
#endif

  vector<RevEngRegion*> regs1;
  regs1.push_back(regions_[ix2].get());
  regs1.insert(regs1.end(), common_reg.begin(), common_reg.end());
  vector<RevEngPoint*> bd_pts1 =
    regions_[ix1]->extractBdPoints(); //regs1);

  vector<RevEngRegion*> regs2;
  regs2.push_back(regions_[ix1].get());
  regs2.insert(regs2.end(), common_reg.begin(), common_reg.end());
  vector<RevEngPoint*> bd_pts2 =
    regions_[ix2]->extractBdPoints();//regs2);
  if (bd_pts1.size() == 0 || bd_pts2.size() == 0)
    return edges;
  
#ifdef DEBUG_EDGE
  std::ofstream of3("bd_pts.g2");
  of3 << "400 1 0 4 255 0 0 255" << std::endl;
  of3 << bd_pts1.size() << std::endl;
  for (size_t kr=0; kr<bd_pts1.size(); ++kr)
    of3 << bd_pts1[kr]->getPoint() << std::endl;
  of3 << "400 1 0 4 255 0 0 255" << std::endl;
  of3 << bd_pts2.size() << std::endl;
  for (size_t kr=0; kr<bd_pts2.size(); ++kr)
    of3 << bd_pts2[kr]->getPoint() << std::endl;
#endif

  // if (int_cvs1.size() > 1)
  //   {
  //     // Ensure consistent curve parameterization. Sort
  //     vector<Point> startpt(int_cvs1.size()), endpt(int_cvs1.size());
  //     for (size_t kr=0; kr<int_cvs1.size(); ++kr)
  // 	{
  // 	  startpt[kr] = int_cvs1[kr]->ParamCurve::point(int_cvs1[kr]->startparam());
  // 	  endpt[kr] = int_cvs1[kr]->ParamCurve::point(int_cvs1[kr]->endparam());
  // 	}
  //     int stop_breakp = 1;
  //   }
  int num_in_lim1=0, num_in_lim2=0;
  vector<double> width2;
  vector<shared_ptr<CurveOnSurface> > cvs1, cvs2;
  for (size_t kr=0; kr<int_cvs1.size(); ++kr)
    {
      vector<pair<double,double> > tmin_tmax;
      vector<double> width;
      vector<shared_ptr<CurveOnSurface> > tmp_int1, tmp_int2;
      tmp_int1.push_back(int_cvs1[kr]);
      tmp_int2.push_back(int_cvs2[kr]);
      vector<pair<double, double> > t5_t6, t7_t8;
      vector<double> wwd1, wwd2;
      int numlim1, numlim2;
      regions_[ix1]->estimateBlendDimensions(tmp_int1, bd_pts1,
					     approx_tol_, blendlim,
					     t5_t6, wwd1, numlim1);

      regions_[ix2]->estimateBlendDimensions(tmp_int2, bd_pts2,
					     approx_tol_, blendlim,
					     t7_t8, wwd2, numlim2);
      if (t5_t6.size() == 0 || t7_t8.size() == 0)
	continue;
      num_in_lim1 += numlim1;
      num_in_lim2 += numlim2;

     // Unify intersection curve limitations
      if (keep_length)
	{
	  tmin_tmax.push_back(std::make_pair(int_cvs1[kr]->startparam(),
					     int_cvs1[kr]->endparam()));
	  double wwd = 0.0;
	  double w1min = std::numeric_limits<double>::max();
	  double w2min = std::numeric_limits<double>::max();
	  for (size_t kh=0; kh<wwd1.size(); ++kh)
	    {
	      wwd += wwd1[kh];
	      w1min = std::min(w1min, wwd1[kr]);
	    }
	  for (size_t kh=0; kh<wwd2.size(); ++kh)
	    {
	      wwd += wwd2[kh];
	      w2min = std::min(w2min, wwd2[kh]);
	    }
	  wwd /= (double)(wwd1.size()+wwd2.size());
	  wwd = std::max(wwd, std::max(w1min, w2min));
	  width.push_back(wwd);
	}
      else
	{
	  size_t kk1, kk2;
	  for (kk1=0, kk2=0; kk1<wwd1.size() && kk2<wwd2.size(); )
	    {
	      for (; kk1<wwd1.size() && t7_t8[kk2].first >= t5_t6[kk1].second; ++kk1);
	      if (kk1 == wwd1.size())
		break;
	      for (; kk2<wwd2.size() && t5_t6[kk1].first >= t7_t8[kk2].second; ++kk2);
	      tmin_tmax.push_back(std::make_pair(std::max(t5_t6[kk1].first,t7_t8[kk2].first),
						 std::min(t5_t6[kk1].second,t7_t8[kk2].second)));
	      width.push_back(0.5*(wwd1[kk1]+wwd2[kk2]));
	      if (t5_t6[kk1].second < t7_t8[kk2].second)
		++kk1;
	      else if (t7_t8[kk2].second < t5_t6[kk1].second)
		++kk2;
	      else
		{
		  ++kk1;
		  ++kk2;
		}
	    }
	}

      size_t ki;
      for (ki=0; ki<width.size(); ++ki)
	{
	  double tmin = tmin_tmax[ki].first;
	  double tmax = tmin_tmax[ki].second;
	  double width3 = width[ki];
	  if (tmax <= tmin+int_tol || width3 > maxwidth)
	    continue;
      
	  double tp1 = std::max(int_cvs1[kr]->startparam(), tmin);
	  double tp2 = std::min(int_cvs1[kr]->endparam(), tmax);
	  if (tp2 <= tp1+int_tol)
	    continue;
	  if (tp1 > int_cvs1[kr]->startparam()+int_tol ||
	      tp2 < int_cvs1[kr]->endparam()-int_tol)
	    {
	      shared_ptr<CurveOnSurface> sub1(int_cvs1[kr]->subCurve(tp1, tp2));
	      shared_ptr<CurveOnSurface> sub2(int_cvs2[kr]->subCurve(tp1, tp2));
	      cvs1.push_back(sub1);
	      cvs2.push_back(sub2);
	    }
	  else if (fabs(tp1-int_cvs1[kr]->startparam()) < int_tol &&
		   fabs(tp2-int_cvs1[kr]->endparam()) < int_tol)
	    {
	      cvs1.push_back(int_cvs1[kr]);
	      cvs2.push_back(int_cvs2[kr]);
	    }
	  width2.push_back(width3);
	}
    }
   
  if (width2.size() == 0)
    return edges;
  if (num_in_lim1 == 0 || num_in_lim2 == 0)
    return edges;
  if (cvs1.size() == 0)
    return edges;

#ifdef DEBUG_EDGE
  std::ofstream of1e("one_edgcv.g2");
  for (size_t kr=0; kr<cvs1.size(); ++kr)
    {
      cvs1[kr]->spaceCurve()->writeStandardHeader(of1e);
      cvs1[kr]->spaceCurve()->write(of1e);
    }
#endif
  for (size_t kj=0; kj<cvs1.size(); ++kj)
    {
      shared_ptr<RevEngEdge> edg = defineOneEdge(ix1, surf1, dir1, ix2,
						 surf2, dir2, cvs1[kj],
						 cvs2[kj], width2[kj],
						 common_reg, only_curve,
						 check_common);
      if (edg.get())
	edges.push_back(edg);
    }

  return edges;
}
	

//===========================================================================
shared_ptr<RevEngEdge> 
RevEng::defineOneEdge(size_t ix1, shared_ptr<ElementarySurface>& surf1,
		      Point& dir1, size_t ix2,
		      shared_ptr<ElementarySurface>& surf2, Point& dir2,
		      shared_ptr<CurveOnSurface>& int_cv1,
		      shared_ptr<CurveOnSurface>& int_cv2,
		      double width, vector<RevEngRegion*>& common_reg,
		      bool only_curve, bool check_common)
//===========================================================================
{
  shared_ptr<RevEngEdge> dummy_edg;
  bool plane1 = (surf1->instanceType() == Class_Plane);
  bool plane2 = (surf2->instanceType() == Class_Plane);
  double angtol = 5.0*anglim_;
  double tol10 = 10.0*approx_tol_;
  int min_pt_blend = 20;

  if (plane1)
    {
      // Make sure that the plane normal poins out
      Point avnorm = regions_[ix1]->getMeanNormal();
      if (dir1*avnorm < 0.0)
	dir1 *= -1;
    }
  if (plane2)
    {
      // Make sure that the plane normal poins out
      Point avnorm = regions_[ix2]->getMeanNormal();
      if (dir2*avnorm < 0.0)
	dir2 *= -1;
    }
  int state = (plane1 || plane2) ? 1 : 2;


  vector<vector<RevEngPoint*> > near_pts(common_reg.size()+2);
  vector<RevEngPoint*> curr_near1, curr_near2;
  double tmin1 = int_cv1->startparam();
  double tmax1 = int_cv1->endparam();
  if (!only_curve)
    regions_[ix1]->getNearPoints(int_cv1, tmin1, tmax1, width,
				 angtol, curr_near1);
  double tmin2 = int_cv2->startparam();
  double tmax2 = int_cv2->endparam();
  if (!only_curve)
    regions_[ix2]->getNearPoints(int_cv2, tmin2, tmax2, width,
				 angtol, curr_near2);
  if (curr_near1.size() == 0 && curr_near2.size() == 0 && (!only_curve))
    return dummy_edg;

  RevEngPoint *distant1 = 0, *distant2 = 0;
  vector<RevEngPoint*> reg1_pts, reg2_pts;
  if (only_curve)
    {
      reg1_pts = regions_[ix1]->getPoints();
      reg2_pts = regions_[ix2]->getPoints();
    }
  distant1 = getDistantPoint(int_cv1, std::max(tmin1, tmin2),
			     std::min(tmax1, tmax2), approx_tol_,
			     width, (curr_near1.size() > 0) ? curr_near1 :
			     reg1_pts);
  distant2 = getDistantPoint(int_cv2, std::max(tmin1, tmin2),
			     std::min(tmax1, tmax2), approx_tol_,
			     width, (curr_near2.size() > 0) ? curr_near2 :
			     reg2_pts);
  if (!distant1)
    {
      Point mid;
      int_cv1->point(mid, 0.5*(int_cv1->startparam()+int_cv1->endparam()));
      double distmid;
      distant1 = regions_[ix1]->closestPoint(mid, distmid);
    }
  if (!distant2)
    {
      Point mid;
      int_cv2->point(mid, 0.5*(int_cv2->startparam()+int_cv2->endparam()));
      double distmid;
      distant2 = regions_[ix2]->closestPoint(mid, distmid);
    }
  if ((!distant1) || (!distant2))
    return dummy_edg;

  vector<Point> der(2);
  Vector3D xyz1 = distant1->getPoint();
  Vector3D xyz2 = distant2->getPoint();
  Point loc1(xyz1[0], xyz1[1], xyz1[2]);
  Point loc2(xyz2[0], xyz2[1], xyz2[2]);
  Point distpt = (plane1) ? loc1 : loc2;
  double dist;
  double tpar;
  Point close;
  int_cv1->closestPoint(distpt, int_cv1->startparam(),
			int_cv1->endparam(), tpar, close, dist);
  int_cv1->point(der, tpar, 1);
		  
  Point loc1_2 = loc1 - ((loc1-der[0])*der[1])*der[1];
  Point loc2_2 = loc2 - ((loc2-der[0])*der[1])*der[1];
  Point vec = loc1_2 - loc2_2;
  bool outer1, outer2;
  Point pos2;
  if (plane1)
    outer1 = (dir1*vec < 0.0);
  else if (plane2)
    {
      // Cylinder or cone combined with plane
      Point pos = surf1->location();
       double vval = (der[0]-pos)*dir1;
     double rad = surf1->radius(0.0,vval);
      pos2 = pos + vval*dir1;
      Point loc2_3 = loc2 - ((loc2-der[0])*dir1)*dir1;
      outer1 = (pos2.dist(loc2_3) < rad);
    }
  else if (surf1->instanceType() == Class_Cylinder)
    {
      // Cylinder combined with cone
      Point axs = surf1->direction();
      Point tmp = loc1 + ((der[0]-loc1)*axs)*axs;
      Point vec2 = loc1 - tmp;
      Vector2D uv = distant2->getPar();
      double rad1 = surf2->radius(0.0, uv[1]);
      double rad2 = surf1->radius(0.0, 0.0);
      outer1 = (rad1 >= rad2);
      if (vec2*axs < 0.0)
	outer1 = (!outer1);
    }
  else
    {
      // Cone combined with cylinder
      Vector2D uv = distant1->getPar();
      double rad1 = surf1->radius(0.0, uv[1]);
      double rad2 = surf2->radius(0.0, 0.0);
      outer1 = (rad1 >= rad2);
    }
  
  if (plane2)
    outer2 = (dir2*vec > 0.0);
  else if (plane1)
    {
      Point pos = surf2->location();
      double vval = (der[0]-pos)*dir2;
      double rad = surf2->radius(0.0,vval);
      pos2 = pos + vval*dir2;
      Point loc1_3 = loc1 - ((loc1-der[0])*dir2)*dir2;
      outer2 = (pos2.dist(loc1_3) < rad);
    }
  else if (surf2->instanceType() == Class_Cylinder)
    {
      // Cylinder combined with cone
      Vector2D uv = distant1->getPar();
      double rad1 = surf1->radius(0.0, uv[1]);
      double rad2 = surf2->radius(0.0, 0.0);
      Point axs = surf2->direction();
      Point tmp = loc2 + ((der[0]-loc2)*axs)*axs;
      Point vec2 = loc2 - tmp;
      outer2 = (rad1 >= rad2);
      if (vec2*axs < 0.0)
	outer2 = (!outer2);
    }
  else
    {
      // Cone combined with cylinder
      Vector2D uv = distant2->getPar();
      double rad1 = surf2->radius(0.0, uv[1]);
      double rad2 = surf1->radius(0.0, 0.0);
      outer2 = (rad1 >= rad2);
    }

  double mindist1=0.0, mindist2=0.0;
  near_pts[0] = (curr_near1.size() == 0 || state == 2) ? curr_near1 :
    regions_[ix2]->removeOutOfSurf(curr_near1, tol10,
				   angtol, outer2, mindist1);
  near_pts[1] = (curr_near2.size() == 0 || state == 2) ? curr_near2 :
    regions_[ix1]->removeOutOfSurf(curr_near2, tol10,
				   angtol, outer1, mindist2);

  if ((near_pts[0].size() == 0 || near_pts[1].size() == 0) && (!only_curve))
    return dummy_edg;
  int nmb_near = (int)(near_pts[0].size() + near_pts[1].size());
  for (size_t kr=0; kr<common_reg.size(); ++kr)
    {
      vector<RevEngPoint*> adj_near1, adj_near2;
      double dummy_min = 0.0;
      double tmin3 = int_cv1->startparam();
      double tmax3 = int_cv1->endparam();
      common_reg[kr]->getNearPoints(int_cv1, tmin3, tmax3, width,
				    angtol, adj_near1);
      if (state == 1)
	{
	  adj_near2 =
	    regions_[ix1]->removeOutOfSurf(adj_near1, tol10,
					   angtol, outer1, dummy_min);
	  near_pts[2+kr] = 
	    regions_[ix2]->removeOutOfSurf(adj_near2, tol10,
					   angtol, outer2, dummy_min);
	}
      else
	near_pts[2+kr] = adj_near1;
      
      nmb_near += (int)near_pts[2+kr].size();
    }
  if (nmb_near < min_pt_blend && (!only_curve))
    return dummy_edg;
  
#ifdef DEBUG_EDGE
  std::ofstream of4("blend_pts.g2");
  for (size_t kr=0; kr<near_pts.size(); ++kr)
    {
      of4 << "400 1 0 0" << std::endl;
      of4 << near_pts[kr].size() << std::endl;
      for (size_t kh=0; kh<near_pts[kr].size(); ++kh)
	of4 << near_pts[kr][kh]->getPoint() << std::endl;
    }
#endif

  // Ensure connection between adjacent regions
  bool adjacent = only_curve ? true : false;
  bool connection = false;
  for (size_t kr=0; kr<near_pts[0].size(); ++kr)
    {
      for (size_t kh=0; kh<near_pts[1].size(); ++kh)
	if (near_pts[0][kr]->isNeighbour(near_pts[1][kh]))
	  {
	    adjacent = true;
	    break;
	  }
      if (adjacent)
	break;
    }
  
  for (size_t kr=2; kr<near_pts.size(); ++kr)
    {
      bool con1 = false, con2 = false;
      for (size_t kh=0; kh<near_pts[kr].size(); ++kh)
	{
	  if (kh != 0 && near_pts[kr][kh]->nextToRegion(regions_[ix1].get()))
	    con1 = true;
	  if (kh != 1 && near_pts[kr][kh]->nextToRegion(regions_[ix2].get()))
	    con2 = true;
	  if (con1 && con2)
	    {
	      connection = true;
	      break;
	    }
	}
      if (connection)
	break;
    }
  if ((!connection) && (!adjacent) && (!regions_[ix1]->hasBlendEdge()) &&
      (!regions_[ix2]->hasBlendEdge()))
    return dummy_edg;
	      
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > blend_groups;
  for (size_t kr=0; kr<near_pts.size(); ++kr)
    blend_groups.push_back(std::make_pair(near_pts[kr].begin(),
					  near_pts[kr].end()));
  double radius = 0.0, ylen = 0.0;
  if (!only_curve)
    {
      if (plane1 && plane2)
	{
	  Point lin1 = der[1].cross(dir1);
	  Point lin2 = der[1].cross(dir2);
	  lin1.normalize();
	  if (lin1*dir2 < 0.0)
	    lin1 *= -1;
	  lin2.normalize();
	  if (lin2*dir1 < 0.0)
	    lin2 *= -1;
	  radius = computeCylinderRadius(near_pts, width, der[0],
					 der[1], lin1, lin2);
	  Point vec = lin1 + lin2;
	  vec.normalize();
	  double alpha = lin1.angle(lin2);
	  double fac = 1.0/sin(0.5*alpha);
	  Point centre = der[0] - fac*radius*vec;
	  //double ang = lin1.angle(lin2);
	  double xlen = der[0].dist(centre);
	  ylen = sqrt(xlen*xlen - radius*radius);
	}
      else
	{
	  int sgn = 1;
	  if ((plane1 && surf1->direction()*dir1 < 0.0) ||
	      (plane2 && surf2->direction()*dir2 < 0.0))
	    sgn = -1;
	  double d2 = 0.0;
	  radius = computeTorusRadius(near_pts, int_cv1, surf1, surf2, 
				      width, outer1, outer2, sgn, d2);


	  Point centre, normal, Cx;
	  double Rrad;
	  bool OK = getTorusParameters(surf1, surf2, der[0], radius, d2, outer1, 
				       outer2, sgn, Rrad, centre, normal, Cx,
				       check_common);
	  if (!OK)
	    return dummy_edg;
	  double xlen = der[0].dist(centre);
	  ylen = fabs(xlen - Rrad);
	}
      ylen *= 1.5;  // Allow some slack
      ylen= std::max(ylen, 1.1*std::max(mindist1, mindist2));
    

      // Reduce near points from adjacent groups according to the updated width
      if (ylen < width)
	{
	  for (size_t kr=2; kr<near_pts.size(); ++kr)
	    {
	      vector<RevEngPoint*> near_pts2;
	      common_reg[kr-2]->getNearPoints2(near_pts[kr], int_cv1, ylen, near_pts2);
	      std::swap(near_pts[kr], near_pts2);
	    }
	}
    }

#ifdef DEBUG_EDGE
  std::ofstream of5("blend_pts2.g2");
  for (size_t kr=0; kr<near_pts.size(); ++kr)
    {
      of5 << "400 1 0 0" << std::endl;
      of5 << near_pts[kr].size() << std::endl;
      for (size_t kh=0; kh<near_pts[kr].size(); ++kh)
	of5 << near_pts[kr][kh]->getPoint() << std::endl;
    }
#endif

  // An intersection/blend curve is found
  // Segment identified adjacent groups according to curve
  vector<RevEngRegion*> adj_regs;
  vector<size_t> remove_ix;
  for (size_t kr=0; kr<common_reg.size(); ++kr)
    {
      if (near_pts[kr+2].size() == 0)
	continue;

      if ((int)near_pts[kr+2].size() < common_reg[kr]->numPoints())
	{
	  int num_init = common_reg[kr]->numPoints();
	  vector<vector<RevEngPoint*> > out_groups;
	  vector<HedgeSurface*> out_sfs;
	  vector<vector<RevEngPoint*> > near_groups;
	  common_reg[kr]->extractSpesPoints(near_pts[kr+2], near_groups);
	  common_reg[kr]->updateInfo();
	  common_reg[kr]->splitRegion(out_groups);
	  if (common_reg[kr]->hasSurface() && common_reg[kr]->numPoints() < num_init/2 &&
	      !(common_reg[kr]->hasRevEdges() || common_reg[kr]->hasTrimEdges()))
	    {
	      int num_sf = common_reg[kr]->numSurface();
	      for (int ka=0; ka<num_sf; ++ka)
		out_sfs.push_back(common_reg[kr]->getSurface(ka));
	      common_reg[kr]->clearSurface();
	    }

	  // Make new region
	  shared_ptr<RevEngRegion> curr_adj(new RevEngRegion(common_reg[kr]->getClassificationType(),
							     common_reg[kr]->getEdgeClassificationType(),
							     near_pts[kr+2]));
	  regions_.push_back(curr_adj);
	  adj_regs.push_back(curr_adj.get());

	  size_t kh=0;
	  for (kh=0; kh<regions_.size(); ++kh)
	    if (common_reg[kr] == regions_[kh].get())
	      break;
	  surfaceExtractOutput((int)kh, out_groups, out_sfs);
	}
      else
	{
	  adj_regs.push_back(common_reg[kr]);
	  remove_ix.push_back(kr);
	}
    }

  for (int ka=(int)remove_ix.size()-1; ka>=0; --ka)
    common_reg.erase(common_reg.begin() + remove_ix[ka]);

  for (size_t kr=0; kr<adj_regs.size(); ++kr)
    {
      vector<vector<RevEngPoint*> > curr_out;
      vector<HedgeSurface*> dummy_surfs;
      adj_regs[kr]->splitRegion(curr_out);
      if (curr_out.size() > 0)
	{
	  size_t kh=0;
	  for (kh=0; kh<regions_.size(); ++kh)
	    if (adj_regs[kr] == regions_[kh].get())
	      break;
	  surfaceExtractOutput((int)kh, curr_out, dummy_surfs);
	}
    }

  for (size_t kr=0; kr<adj_regs.size(); ++kr)
    adj_regs[kr]->updateRegionAdjacency();
  regions_[ix1]->updateRegionAdjacency();
  regions_[ix2]->updateRegionAdjacency();
  
  int edge_type = (only_curve) ? NOT_BLEND : BLEND_NOT_SET;
  vector<shared_ptr<CurveOnSurface> > int_cvs1(1, int_cv1);
  vector<shared_ptr<CurveOnSurface> > int_cvs2(1, int_cv2);
  shared_ptr<RevEngEdge> edge(new RevEngEdge(edge_type, regions_[ix1].get(),
					     int_cvs1, outer1, regions_[ix2].get(),
					     int_cvs2, outer2, radius,
					     std::min(width, ylen)));
  regions_[ix1]->addRevEdge(edge.get());
  regions_[ix2]->addRevEdge(edge.get());
  if (adj_regs.size() > 0)
    edge->addBlendRegions(adj_regs);
  for (size_t kr=0; kr<adj_regs.size(); ++kr)
    adj_regs[kr]->setAssociatedBlend(edge.get());

  return edge;
}

					  
//===========================================================================
RevEngPoint* RevEng::getDistantPoint(shared_ptr<CurveOnSurface>& cv,
				     double tmin, double tmax, double dd,
				     double width, vector<RevEngPoint*>& points)
//===========================================================================
{
  RevEngPoint *distant = 0;
  double ptdist = 0.0;
  double eps = 1.0e-9;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      if (points[ki]->getSurfaceDist() > dd)
	continue;
      
      double tpar, dist;
      Point close;
      Vector3D xyz = points[ki]->getPoint();
      Point pt(xyz[0], xyz[1], xyz[2]);
      cv->closestPoint(pt, cv->startparam(), cv->endparam(),
		       tpar, close, dist);
      if (tpar > tmin+eps && tpar < tmax-eps && dist > ptdist &&
	  dist < width)
	{
	  ptdist = dist;
	  distant = points[ki];
	}
    }
  return distant;
}



//===========================================================================
void RevEng::computeAxisFromCylinder(Point initaxis[3], int min_num, double max_ang,
				     Point axis[3], int num_points[3])
//===========================================================================
{
  vector<vector<pair<std::vector<RevEngPoint*>::iterator,
		     std::vector<RevEngPoint*>::iterator> > > points(3);
  num_points[0] =  num_points[1] = num_points[2] = 0;
  for (int ka=0; ka<3; ++ka)
    axis[ka] = initaxis[ka];

  Point axis0[3];
  for (int ka=0; ka<3; ++ka)
    axis0[ka] = initaxis[ka];
  double max_ang2 = 0.5*M_PI - 1.0e-9;
  for (int kc=0; kc<2; ++kc)
    {
      Point dummy(0.0, 0.0, 0.0);
      vector<double> min_angle(3, max_ang2);
      vector<Point> min_axis(3, dummy);
   
      for (size_t ki=0; ki<surfaces_.size(); ++ki)
	{
	  int sfcode;
	  ClassType type = surfaces_[ki]->instanceType(sfcode);
	  if ((type != Class_Cylinder) && (type != Class_Cone))
	    continue;
      
	  shared_ptr<ElementarySurface> curr =
	    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surfaces_[ki]->surface());
	  if (!curr.get())
	    continue;

	  int num_reg = surfaces_[ki]->numRegions();
	  if (num_reg == 0)
	    continue; // Should not happen

	  RevEngRegion *reg0 = surfaces_[ki]->getRegion(0);
	  if (reg0->hasBlendEdge())
	    continue;  // Derived information
      
	  if (reg0->hasAssociatedBlend())
	    continue;  // Outdated information
      
	  int num_pts = surfaces_[ki]->numPoints();
	  if (num_pts < min_num)
	    continue;

	  Point vec = curr->direction();
	  int ix = -1;
	  double min_angle0 = std::numeric_limits<double>::max();
	  for (int ka=0; ka<3; ++ka)
	    {
	      double ang = vec.angle(axis0[ka]);
	      ang = std::min(ang, M_PI-ang);
	      if (ang < min_angle0)
		{
		  min_angle0 = ang;
		  ix = ka;
		}
	      if (ang < max_ang)
		{
		  for (int kb=0; kb<num_reg; ++kb)
		    {
		      RevEngRegion *reg = surfaces_[ki]->getRegion(kb);
		      points[ka].push_back(std::make_pair(reg->pointsBegin(),
							  reg->pointsEnd()));
		    }
		  num_points[ka] += num_pts;
		}
	    }
	  if (ix >= 0)
	    {
	      min_angle[ix] = min_angle0;
	      min_axis[ix] = vec;
	    }
	}
      if (num_points[0]+num_points[1]+num_points[2] == 0 && kc == 0)
	{
	  for (int kb=0; kb<3; ++kb)
	    if (min_axis[kb].length() > 0.0)
	      axis0[kb] = min_axis[kb];
	}
      else if (kc == 0)
	break;
    }

  Point Cx, Cy;
  for (int ka=0; ka<3; ++ka)
    {
      if (points[ka].size() > 0)
	{
	  RevEngUtils::computeAxis(points[ka], axis[ka], Cx, Cy);
	}
      else
	axis[ka] = initaxis[ka];
    }

}

//===========================================================================
void RevEng::computeAxisFromPlane(Point initaxis[3], int min_num, double max_ang,
				  Point axis[3], int num_points[3])
//===========================================================================
{
  vector<vector<shared_ptr<Plane> > > planes(3);
  vector<vector<int> > num(3);
  vector<vector<double> > avd(3);
  num_points[0] =  num_points[1] = num_points[2] = 0;
  Point axis0[3];
  for (int ka=0; ka<3; ++ka)
    axis0[ka] = initaxis[ka];
  double max_ang2 = 0.5*M_PI - 1.0e-9;
  for (int kc=0; kc<2; ++kc)
    {
      Point dummy(0.0, 0.0, 0.0);
      vector<double> min_angle(3, max_ang2);
      vector<Point> min_axis(3, dummy);
      for (size_t ki=0; ki<surfaces_.size(); ++ki)
	{
	  int sfcode;
	  ClassType type = surfaces_[ki]->instanceType(sfcode);
	  if (type != Class_Plane)
	    continue;

	  shared_ptr<Plane> curr =
	    dynamic_pointer_cast<Plane,ParamSurface>(surfaces_[ki]->surface());
	  if (!curr.get())
	    continue;

	  int num_pts = surfaces_[ki]->numPoints();
	  if (num_pts < min_num)
	    continue;
      
	  // Check for a more accurate base surface
	  double avdist = 0.0;
	  int num_reg = surfaces_[ki]->numRegions();
	  double fac = 1.0/(double)num_pts;
	  if (num_reg > 0)
	    {
	      // Should always be the case
	      RevEngRegion* reg = surfaces_[ki]->getRegion(0);
	      shared_ptr<ParamSurface> base = reg->getBase();
	      if (base.get() && base->instanceType() == Class_Plane)
		{
		  double av1 = 0.0, av2 = 0.0;
		  //int nmb = 0;
		  for (int ka=0; ka<num_reg; ++ka)
		    {
		      reg = surfaces_[ki]->getRegion(ka);
		      double maxds, avds, maxdb, avdb;
		      int num_in, num2_in, num_inb, num2_inb;
		      reg->getAccuracy(maxds, avds, num_in, num2_in);
		      reg->getBaseDist(maxdb, avdb, num_inb, num2_inb);
		      int nn = reg->numPoints();
		      av1 += (double)nn*fac*avds;
		      av2 += (double)nn*fac*avdb;
		    }
		  if (av2 < av1)
		    {
		      curr = dynamic_pointer_cast<Plane,ParamSurface>(base);
		      avdist = av2;
		    }
		  else
		    avdist = av1;
		}
	    }

	  Point normal = curr->direction();
	  int ix = -1;
	  double min_angle0 = std::numeric_limits<double>::max();
	  for (int kb=0; kb<3; ++kb)
	    {
	      double ang = normal.angle(axis0[kb]);
	      ang = std::min(ang, M_PI-ang);
	      if (ang < min_angle0)
		{
		  min_angle0 = ang;
		  ix = kb;
		}
	      if (ang <= max_ang)
		{
		  planes[kb].push_back(curr);
		  num[kb].push_back(num_pts);
		  avd[kb].push_back(avdist);
		  num_points[kb] += num_pts;
		}
	    }
	  if (ix >= 0)
	    {
	      min_angle[ix] = min_angle0;
	      min_axis[ix] = normal;
	    }
	}
      if (num_points[0]+num_points[1]+num_points[2] == 0 && kc == 0)
	{
	  for (int kb=0; kb<3; ++kb)
	    if (min_axis[kb].length() > 0.0)
	      axis0[kb] = min_axis[kb];
	}
      else if (kc == 0)
	break;
    }
  
  for (int kb=0; kb<3; ++kb)
    {
      if (planes[kb].size() == 0)
	{
	  axis[kb] = initaxis[kb];
	  continue;
	}

      axis[kb] = Point(0.0, 0.0, 0.0);
      double fac = 1.0/(double)num_points[kb];
      for (size_t ki=0; ki<planes[kb].size(); ++ki)
	{
	  Point normal = planes[kb][ki]->direction();
	  if (normal*axis0[kb] < 0.0)
	    normal *= -1.0;

	  double wgt = fac*(1.0 - avd[kb][ki])*num[kb][ki];
	  wgt = std::max(wgt, 0.0);
	  axis[kb] += wgt*normal;
	}
      axis[kb].normalize_checked();
    }
}

//===========================================================================
bool RevEng::recognizeOneSurface(int& ix, int min_point_in, double angtol,
				 int pass)
//===========================================================================
{
  bool firstpass = (pass <= 2);
  if (regions_[ix]->tryOtherSurf(prefer_elementary_, firstpass) /*&&
								  regions_[ix]->feasiblePlane(zero_H_, zero_K_)*/)
    {
      vector<shared_ptr<HedgeSurface> > plane_sfs;
      vector<HedgeSurface*> prev_surfs;
      vector<vector<RevEngPoint*> > out_groups;
      bool found = regions_[ix]->extractPlane(mainaxis_,
					      approx_tol_, min_point_in,
					      min_point_region_, angtol, 
					      prefer_elementary_,
					      plane_sfs, prev_surfs, out_groups);

      if (out_groups.size() > 0 || prev_surfs.size() > 0)
	surfaceExtractOutput(ix, out_groups, prev_surfs);

      if (plane_sfs.size() > 0)
	surfaces_.insert(surfaces_.end(), plane_sfs.begin(), plane_sfs.end());
      if (found && regions_[ix]->getMaxSfDist() < 0.5*approx_tol_)
	return true;  // Result accepted
    }

  if (regions_[ix]->tryOtherSurf(prefer_elementary_, firstpass))
    {
      vector<shared_ptr<HedgeSurface> > cyl_sfs;
      vector<HedgeSurface*> prev_surfs;
      vector<vector<RevEngPoint*> > out_groups;
      bool repeat = false;
      
      // Cylinder from context information
      bool found =
	regions_[ix]->contextCylinder(mainaxis_, approx_tol_,
				      min_point_in, min_point_region_, 
				      angtol, prefer_elementary_,
				      cyl_sfs, prev_surfs,
				      out_groups);

      if (found == false && regions_[ix]->feasibleCylinder(zero_H_, zero_K_))
	(void)regions_[ix]->extractCylinder(approx_tol_, min_point_in,
					    min_point_region_, angtol,
					    prefer_elementary_,
					    cyl_sfs, prev_surfs,
					    out_groups, repeat);

      
      if (out_groups.size() > 0 || prev_surfs.size() > 0)
	surfaceExtractOutput(ix, out_groups, prev_surfs);
	
      if (repeat)
	{
	  --ix;
	  return false;  // New try
	}
	

      if (cyl_sfs.size() > 0)
	surfaces_.insert(surfaces_.end(), cyl_sfs.begin(), cyl_sfs.end());
      if (found && regions_[ix]->getMaxSfDist() < 0.5*approx_tol_)
	return true;
    }
  
  if (regions_[ix]->tryOtherSurf(prefer_elementary_, firstpass) && (!firstpass) &&
      regions_[ix]->hasSweepInfo() && regions_[ix]->sweepType() == 1)
    {
      vector<shared_ptr<HedgeSurface> > linsweep_sfs;
      vector<HedgeSurface*> prev_surfs;
      vector<vector<RevEngPoint*> > out_groups;
      bool found = regions_[ix]->extractLinearSweep(approx_tol_, min_point_in,
						 min_point_region_, angtol, 
						 prefer_elementary_,
						 linsweep_sfs, 
						 prev_surfs);
      if (out_groups.size() > 0 || prev_surfs.size() > 0)
	surfaceExtractOutput(ix, out_groups, prev_surfs);
	
      if (linsweep_sfs.size() > 0)
	surfaces_.insert(surfaces_.end(), linsweep_sfs.begin(), linsweep_sfs.end());
	   
      if (found && regions_[ix]->getMaxSfDist() < 0.5*approx_tol_)
	return true;
    }
      
  if (regions_[ix]->tryOtherSurf(prefer_elementary_, firstpass) && pass > 1)
    {
      vector<shared_ptr<HedgeSurface> > tor_sfs;
      vector<HedgeSurface*> prev_surfs;
      vector<vector<RevEngPoint*> > out_groups;

      // Torus from context information
      bool found =  regions_[ix]->contextTorus(mainaxis_, approx_tol_,
					       min_point_in, min_point_region_, 
					       angtol, prefer_elementary_,
					       tor_sfs, prev_surfs,
					       out_groups);
      if (regions_[ix]->numPoints() == 0)
	{
	  regions_.erase(regions_.begin()+ix);
	  --ix;
	  return false;
	}

      if (!found)
	found = regions_[ix]->extractTorus(mainaxis_, approx_tol_, min_point_in,
					   min_point_region_,
					   angtol, prefer_elementary_,
					   tor_sfs, prev_surfs,
					   out_groups);

      if (out_groups.size() > 0 || prev_surfs.size() > 0)
	surfaceExtractOutput(ix, out_groups, prev_surfs);

      if (tor_sfs.size() > 0)
	surfaces_.insert(surfaces_.end(), tor_sfs.begin(), tor_sfs.end());
      if (found && regions_[ix]->getMaxSfDist() < 0.5*approx_tol_)
	return true;
    }

  if (regions_[ix]->tryOtherSurf(prefer_elementary_, firstpass) && pass > 1)
    {
      vector<shared_ptr<HedgeSurface> > sph_sfs;
      vector<HedgeSurface*> prev_surfs;
      vector<vector<RevEngPoint*> > out_groups;
      bool found = regions_[ix]->extractSphere(mainaxis_, approx_tol_,
					       min_point_in, min_point_region_,
					       angtol, prefer_elementary_,
					       sph_sfs, prev_surfs,
					       out_groups);

      
      if (out_groups.size() > 0 || prev_surfs.size() > 0)
	surfaceExtractOutput(ix, out_groups, prev_surfs);
	
      if (sph_sfs.size() > 0)
	surfaces_.insert(surfaces_.end(), sph_sfs.begin(), sph_sfs.end());
      if (found && regions_[ix]->getMaxSfDist() < 0.5*approx_tol_)
	return true;
    }

  if (regions_[ix]->tryOtherSurf(prefer_elementary_, firstpass) && pass > 1)
    {
      vector<shared_ptr<HedgeSurface> > cone_sfs;
      vector<HedgeSurface*> prev_surfs;
      vector<vector<RevEngPoint*> > out_groups;
      bool found =
	regions_[ix]->extractCone(approx_tol_, min_point_in, min_point_region_, 
				  angtol, prefer_elementary_,
				  cone_sfs, prev_surfs,
				  out_groups);
	   
      if (out_groups.size() > 0 || prev_surfs.size() > 0)
	surfaceExtractOutput(ix, out_groups, prev_surfs);
	
      if (cone_sfs.size() > 0)
	surfaces_.insert(surfaces_.end(), cone_sfs.begin(), cone_sfs.end());
      if (found && regions_[ix]->getMaxSfDist() < 0.5*approx_tol_)
	return true;
    }

  if (regions_[ix]->tryOtherSurf(prefer_elementary_, firstpass) && (!firstpass)) 
    {
      vector<shared_ptr<HedgeSurface> > spl_sfs;
      vector<HedgeSurface*> prev_surfs;
      vector<vector<RevEngPoint*> > out_groups;
      bool found = 
	regions_[ix]->extractFreeform(approx_tol_, min_point_in,
				      min_point_region_, angtol,
				      prefer_elementary_,
				      spl_sfs, prev_surfs,
				      out_groups);

      if (out_groups.size() > 0 || prev_surfs.size() > 0)
	surfaceExtractOutput(ix, out_groups, prev_surfs);

      if (spl_sfs.size() > 0)
	surfaces_.insert(surfaces_.end(), spl_sfs.begin(), spl_sfs.end());
    }

  if (!regions_[ix]->hasSurface())
    {
      vector<shared_ptr<HedgeSurface> > adj_sfs;
      vector<HedgeSurface*> prev_surfs;
      vector<vector<RevEngPoint*> > out_groups;
       bool found =
	regions_[ix]->adjacentToCylinder(mainaxis_, approx_tol_, min_point_in,
					 min_point_region_, angtol,
					 prefer_elementary_,
					 adj_sfs, prev_surfs,
					 out_groups);
      
      if (out_groups.size() > 0 || prev_surfs.size() > 0)
	surfaceExtractOutput(ix, out_groups, prev_surfs);

      if (adj_sfs.size() > 0)
	surfaces_.insert(surfaces_.end(), adj_sfs.begin(), adj_sfs.end());
    }
#ifdef DEBUG
  if (regions_[ix]->hasSurface())
    {
      std::ofstream of("one_surface.g2");
      regions_[ix]->writeSurface(of);
    }
#endif
  return (regions_[ix]->hasSurface());
}

//===========================================================================
void RevEng::recognizeSurfaces(int min_point_in, int pass)
//===========================================================================
{
  double angfac = 5.0;
  double angtol = angfac*anglim_;
  int pass2 = pass + 1;
  std::sort(regions_.begin(), regions_.end(), sort_region);
#ifdef DEBUG
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      std::set<RevEngPoint*> tmpset(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
      if ((int)tmpset.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch, pre recognizeSurfaces. " << ki << " " << tmpset.size() << " " << regions_[ki]->numPoints() << std::endl;
      int num = regions_[ki]->numPoints();
      for (int ka=0; ka<num; ++ka)
	if (regions_[ki]->getPoint(ka)->region() != regions_[ki].get())
	  std::cout << "Inconsistent region pointers, pre recognizeSurfaces: " << ki << " " << ka << std::endl;
    }
#endif

  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
#ifdef DEBUG_CHECK
      std::set<RevEngPoint*> tmpset(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
      if ((int)tmpset.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch, pre. " << ki << " " << tmpset.size() << " " << regions_[ki]->numPoints() << std::endl;
      
 #endif
    }
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
#ifdef DEBUG_CHECK
      std::set<RevEngPoint*> tmpset(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
      if ((int)tmpset.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch. " << ki << " " << tmpset.size() << " " << regions_[ki]->numPoints() << std::endl;
      for (size_t kj=0; kj<regions_.size(); ++kj)
	{
	  int nump = regions_[kj]->numPoints();
	  for (int ka=0; ka<nump; ++ka)
	    if (regions_[kj]->getPoint(ka)->region() != regions_[kj].get())
	      std::cout << "Inconsistent region pointers, recognizeSurfaces: " << ki << " " << kj << " " << ka << std::endl;
	}
#endif
      if (regions_[ki]->numPoints() < min_point_region_)
	continue;
      if (regions_[ki]->hasAssociatedBlend())
	continue;  // Treated separately
#ifdef DEBUGONE
      //int classtype = regions_[ki]->getClassification();
      //std::cout << "No " << ki <<  ", classtype: " << classtype << std::endl;
      
      std::ofstream of1("region.g2");
      regions_[ki]->writeRegionInfo(of1);
      if (regions_[ki]->hasSurface())
	regions_[ki]->writeSurface(of1);
      std::ofstream of2("unitsphere.g2");
      regions_[ki]->writeUnitSphereInfo(of2);
#endif
      if (regions_[ki]->hasBlendEdge())
	continue;
      if (regions_[ki]->hasRevEdges())
	continue;
      if (regions_[ki]->getAdaptionHistory() > INITIAL)
	continue; // Do not redo
      
      bool pot_blend = false;
      if (!regions_[ki]->hasSurface())
	pot_blend = regions_[ki]->potentialBlend(angtol);
      if (pot_blend)
	{
#ifdef DEBUG
	  std::cout << "Potential blend" << std::endl;
#endif
	  bool done = setBlendEdge(ki);
	  if (done)
	    continue;
	}
      bool found = recognizeOneSurface(ki, min_point_in, angtol, pass2);
#ifdef DEBUG_CHECK
      bool connect = regions_[ki]->isConnected();
      if (!connect)
	std::cout << "recognizeOneSurface, disconnected region " << ki << std::endl;
#endif
        
#ifdef DEBUG_CHECK
  for (int kj=0; kj<(int)regions_.size(); ++kj)
    {
      std::set<RevEngPoint*> tmpset2(regions_[kj]->pointsBegin(), regions_[kj]->pointsEnd());
      if ((int)tmpset2.size() != regions_[kj]->numPoints())
	std::cout << "Point number mismatch 2. " << kj << " " << ki << " " << tmpset2.size() << " " << regions_[kj]->numPoints() << std::endl;
    }
#endif
  if ((!regions_[ki]->hasSurface()) ||
      (regions_[ki]->hasSurface() &&
       (regions_[ki]->getSurfaceFlag() == ACCURACY_POOR ||
	regions_[ki]->getSurfaceFlag() == NOT_SET)))
	{
	  // Still no surface. Try to divide composite regions into smaller
	  // pieces
	  bool split = segmentComposite(ki, min_point_in, angtol);
#ifdef DEBUG_CHECK
	  if (ki >= 0)
	    {
	      bool connect = regions_[ki]->isConnected();
	      if (!connect)
		std::cout << "segmentComposite, disconnected region " << ki << std::endl;
	    }
#endif
	  int stop_split = 1;
	}
   }

#ifdef DEBUG
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      std::set<RevEngPoint*> tmpset(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
      if ((int)tmpset.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch, post recognizeSurfaces. " << ki << " " << tmpset.size() << " " << regions_[ki]->numPoints() << std::endl;
      int num = regions_[ki]->numPoints();
      for (int ka=0; ka<num; ++ka)
	if (regions_[ki]->getPoint(ka)->region() != regions_[ki].get())
	  std::cout << "Inconsistent region pointers, post recognizeSurfaces: " << ki << " " << ka << std::endl;
    }
#endif

  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

}


//===========================================================================
bool RevEng::segmentComposite(int& ix, int min_point_in, double angtol)
//===========================================================================
{
  // To be updated
  double frac_norm_lim = 0.2; //0.75;
  bool segmented = false;
  
  if (regions_[ix]->numPoints() < min_point_region_)
    return false;

  if (regions_[ix]->hasSurface() &&
      (!(regions_[ix]->getSurfaceFlag() == ACCURACY_POOR ||
	 regions_[ix]->getSurfaceFlag() == NOT_SET)))
    return false;

  if (regions_[ix]->hasSweepInfo())
    return false;

#ifdef DEBUG_DIV
  std::cout << "Segment composite, ix=" << ix << std::endl;
  std::ofstream ofs("region_to_segm.g2");
  regions_[ix]->writeRegionInfo(ofs);
  std::ofstream of2("unitsphere_segm.g2");
  regions_[ix]->writeUnitSphereInfo(of2);
  Point origo(0.0, 0.0, 0.0);
  of2 << "410 1 0 4 0 0 0 255" << std::endl;
  of2 << "3" << std::endl;
  for (int ka=0; ka<3; ++ka)
    of2 << origo << " " << mainaxis_[ka] << std::endl;

  std::ofstream ofa("adjacent_to_segm.g2");
  regions_[ix]->writeAdjacentPoints(ofa);
  //regions_[ix]->sortByAxis(mainaxis_, min_point_in, approx_tol_);
#endif

  size_t num_points = regions_[ix]->numPoints();
  
  bool plane_grow = false;
  if (plane_grow && regions_[ix]->getFracNorm() > frac_norm_lim)
    {
      // Grow sub regions according to surface type
#ifdef DEBUG_DIV
      //std::cout << "Grow planes" << std::endl;
#endif
      segmented = segmentByPlaneGrow(ix, min_point_in, angtol);
#ifdef DEBUG_DIV
      if (segmented)
	std::cout << "Segmented by grow planes" << std::endl;
#endif
    }

  bool repeat = false;
  vector<shared_ptr<HedgeSurface> > hedgesfs;
  vector<shared_ptr<RevEngRegion> > added_reg;
  vector<vector<RevEngPoint*> > separate_groups;
  vector<RevEngPoint*> single_points;
  vector<HedgeSurface*> prevsfs;

  if (!segmented && regions_[ix]->hasDivideInfo())
    {
      for (int ki=0; ki<regions_[ix]->numDivideInfo(); ++ki)
	{
	  segmented = regions_[ix]->divideWithSegInfo(ki,
						      min_point_region_,
						      separate_groups,
						      single_points);
	  if (segmented)
	    break;
	}
    }
  
  if (!segmented)
    {
      // Check if a segmentation into several cylinder like
      // regions is feasible
      double avH, avK, MAH, MAK;
      regions_[ix]->getAvCurvatureInfo(avH, avK, MAH, MAK);
      double fac = 5.0;
      if (MAH > fac*MAK)
	{
	  segmented = regions_[ix]->extractCylByAxis(mainaxis_, min_point_in,
						     min_point_region_,
						     approx_tol_, angtol,
						     prefer_elementary_,
						     hedgesfs, added_reg,
						     separate_groups,
						     single_points);
	}
    }

#ifdef DEBUG_CHECK
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      std::set<RevEngPoint*> tmpset1(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
      if ((int)tmpset1.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch, composite 1. " << ki << " " << tmpset1.size() << " " << regions_[ki]->numPoints() << std::endl;
    }
#endif
    
  if (!segmented)
    {
      vector<RevEngRegion*> adj_planar = regions_[ix]->fetchAdjacentPlanar();
      if (adj_planar.size() > 1)
	{
	  segmented = regions_[ix]->segmentByPlaneAxis(mainaxis_, min_point_in,
						       min_point_region_,
						       approx_tol_, angtol,
						       prefer_elementary_,
						       adj_planar, hedgesfs,
						       added_reg,
						       prevsfs, separate_groups);
#ifdef DEBUG_CHECK
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      std::set<RevEngPoint*> tmpset2(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
      if ((int)tmpset2.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch, composite 2. " << ki << " " << tmpset2.size() << " " << regions_[ki]->numPoints() << std::endl;
    }
    
  for (size_t kr=0; kr<adj_planar.size(); ++kr)
    {
      std::set<RevEngPoint*> tmpset(adj_planar[kr]->pointsBegin(),
				    adj_planar[kr]->pointsEnd());
      if ((int)tmpset.size() != adj_planar[kr]->numPoints())
	std::cout << "Point number mismatch (ByPlaneAxis). " << kr << " " << tmpset.size() << " " << adj_planar[kr]->numPoints() << std::endl;
    }
	  
#endif
	}
      
      if (!segmented)
	{
	  // Extend with cylindrical
	  vector<RevEngRegion*> adj_cyl =
	    regions_[ix]->fetchAdjacentCylindrical();
	  if (adj_cyl.size() > 0)
	    adj_planar.insert(adj_planar.end(), adj_cyl.begin(),
			      adj_cyl.end());
	  if (adj_planar.size() > 0)
	    segmented =
	      regions_[ix]->segmentByAdjSfContext(mainaxis_, min_point_in,
						  min_point_region_,
						  approx_tol_, angtol,
						  adj_planar, separate_groups);
#ifdef DEBUG_CHECK
	  for (int ki=0; ki<(int)regions_.size(); ++ki)
	    {
	      std::set<RevEngPoint*> tmpset4(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
	      if ((int)tmpset4.size() != regions_[ki]->numPoints())
		std::cout << "Point number mismatch, composite 4. " << ki << " " << tmpset4.size() << " " << regions_[ki]->numPoints() << std::endl;
	    }
	  for (size_t kr=0; kr<adj_planar.size(); ++kr)
	    {
	      std::set<RevEngPoint*> tmpset(adj_planar[kr]->pointsBegin(),
					    adj_planar[kr]->pointsEnd());
	      if ((int)tmpset.size() != adj_planar[kr]->numPoints())
		std::cout << "Point number mismatch (ByAdjSfContext). " << kr << " " << tmpset.size() << " " << adj_planar[kr]->numPoints() << std::endl;
	    }
#endif
	}
    }
  
  if (added_reg.size() > 0)
    repeat = true;
  else if (separate_groups.size() > 0)
    {
      int num = 0;
      for (size_t ki=0; ki<separate_groups.size(); ++ki)
	num += (int)separate_groups[ki].size();
      if (num > (int)num_points/10)
	repeat = true;
    }
  
#ifdef DEBUG_CHECK
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      std::set<RevEngPoint*> tmpset5(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
      if ((int)tmpset5.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch, composite 5. " << ki << " " << tmpset5.size() << " " << regions_[ki]->numPoints() << std::endl;
    }
#endif
#ifdef DEBUG
  if (regions_[ix]->numPoints() < (int)num_points)
    {
      std::ofstream of("seg_by_context.g2");
      int num = regions_[ix]->numPoints();
      of << "400 1 0 4 255 0 0 255" << std::endl;
      of <<  num << std::endl;
      for (int ka=0; ka<num; ++ka)
	of << regions_[ix]->getPoint(ka)->getPoint() << std::endl;

      for (size_t ki=0; ki<separate_groups.size(); ++ki)
	{
	  of << "400 1 0 4 0 255 0 255" << std::endl;
	  of <<  separate_groups[ki].size() << std::endl;
	  for (int ka=0; ka<(int)separate_groups[ki].size(); ++ka)
	    of << separate_groups[ki][ka]->getPoint() << std::endl;
	}
      
      for (size_t ki=0; ki<added_reg.size(); ++ki)
	{
	  int num = added_reg[ki]->numPoints();
	  of << "400 1 0 4 0 0 255 255" << std::endl;
	  of <<  num << std::endl;
	  for (int ka=0; ka<num; ++ka)
	    of << added_reg[ki]->getPoint(ka)->getPoint() << std::endl;
	  if (added_reg[ki]->hasSurface())
	    added_reg[ki]->writeSurface(of);
	}
    }
#endif

  // Update adjacency information for current region
  regions_[ix]->updateRegionAdjacency();
  
#ifdef DEBUG_CHECK
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      std::set<RevEngPoint*> tmpset6(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
      if ((int)tmpset6.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch, composite 6. " << ki << " " << tmpset6.size() << " " << regions_[ki]->numPoints() << std::endl;
    }
#endif
  if (added_reg.size() > 0 && regions_[ix]->hasSurface() == false)
    {
      // Swap sequence of regions
      int max_pt = 0;
      int max_ix = -1;
      for (size_t ki=0; ki<added_reg.size(); ++ki)
	if (added_reg[ki]->numPoints() > max_pt)
	  {
	    max_pt = added_reg[ki]->numPoints();
	    max_ix = (int)ki;
	  }
      if (max_ix >= 0)
	std::swap(regions_[ix], added_reg[max_ix]);
    }
#ifdef DEBUG_CHECK
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      std::set<RevEngPoint*> tmpset6_2(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
      if ((int)tmpset6_2.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch, composite 6_2. " << ki << " " << tmpset6_2.size() << " " << regions_[ki]->numPoints() << std::endl;
    }
#endif
  if (separate_groups.size() > 0)
    {
      vector<HedgeSurface*> prev_surfs;
      surfaceExtractOutput(ix, separate_groups, prev_surfs);
    }

#ifdef DEBUG_CHECK
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      std::set<RevEngPoint*> tmpset7(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
      if ((int)tmpset7.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch, composite 7. " << ki << " " << tmpset7.size() << " " << regions_[ki]->numPoints() << std::endl;
    }
  if (!regions_[ix]->isConnected())
    std::cout << "Disconnected region (split), ix= " << ix << " " << regions_[ix].get() << std::endl;
#endif
  if (added_reg.size() > 0)
    regions_.insert(regions_.end(), added_reg.begin(), added_reg.end());
  if (hedgesfs.size() > 0)
    surfaces_.insert(surfaces_.end(), hedgesfs.begin(), hedgesfs.end());
  if (repeat && (!regions_[ix]->hasSurface()))
    {
#ifdef DEBUG_DIV
      std::cout << "Repeat ix = " << ix << std::endl;
#endif
      --ix;
    }
  
  return segmented;
}

//===========================================================================
void RevEng::adjustPointRegions(int min_point_in)
//===========================================================================
{
  double angfac = 5.0;
  double angtol = angfac*anglim_;
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      if (regions_[ki]->hasSurface())
	continue;
      if (regions_[ki]->hasAssociatedBlend())
	continue;
      if (regions_[ki]->numPoints() < min_point_in)
	continue;
      
      vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*>  > adj_elem, adj_elem_base;
      regions_[ki]->getAdjacentElemInfo(adj_elem, adj_elem_base);
      if (adj_elem.size() == 0)
	continue;

      vector<RevEngRegion*> adj_groups(adj_elem.size());
      for (size_t kj=0; kj<adj_elem.size(); ++kj)
	adj_groups[kj] = adj_elem[kj].second;

      vector<vector<RevEngPoint*> > out_groups;
      (void)regions_[ki]->segmentByAdjSfContext(mainaxis_, min_point_in,
						min_point_region_,
						approx_tol_,
						angtol, adj_groups,
						out_groups);
      if (regions_[ki]->numPoints() == 0)
	{
	  vector<RevEngRegion*> adjacent;
	  regions_[ki]->getAdjacentRegions(adjacent);
	  for (size_t kj=0; kj<adjacent.size(); ++kj)
	    adjacent[kj]->removeAdjacentRegion(regions_[ki].get());

	  regions_.erase(regions_.begin()+ki);
	  --ki;
	}
      
      if (out_groups.size() > 0)
	{
	  vector<HedgeSurface*> dummy;
	  surfaceExtractOutput((int)ki, out_groups, dummy);
	}

    }
}

//===========================================================================
void RevEng::surfaceCreation(int pass)
//===========================================================================
{
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }
  
  int min_point_in = 50; //10; //20;
  adjustPointRegions(min_point_in);
  
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }
  
  // First pass. Recognize elementary surfaces
  recognizeSurfaces(min_point_in, pass);

#ifdef DEBUG
  checkConsistence("Regions9");

   if (regions_.size() > 0)
    {
      std::cout << "Regions9" << std::endl;
      std::ofstream of("regions9.g2");
      std::ofstream ofm("mid_regions9.g2");
      std::ofstream ofs("small_regions9.g2");
      writeRegionStage(of, ofm, ofs);
      std::ofstream of0("regions9_helix.g2");
      for (int ki=0; ki<(int)regions_.size(); ++ki)
	{
	  if (regions_[ki]->getSurfaceFlag() == PROBABLE_HELIX)
	    {
	      regions_[ki]->writeRegionInfo(of0);
	      if (regions_[ki]->hasSurface())
		regions_[ki]->writeSurface(of0);
	    }
	}
      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf9.g2");
	  writeRegionWithSurf(of);
	}
   
   if (edges_.size() > 0)
     {
       std::ofstream ofe("edges9.g2");
       writeEdgeStage(ofe);
     }

   std::cout << "Merge adjacent regions, regions: " << regions_.size() << ", surfaces: " << surfaces_.size() << std::endl;
#endif
   double angtol = 5.0*anglim_;
   for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      if (regions_[ki]->hasSurface())
	{
	  vector<RevEngRegion*> grown_regions;
	  vector<HedgeSurface*> adj_surfs;
	  vector<RevEngEdge*> adj_edgs;
	  regions_[ki]->mergeAdjacentSimilar(approx_tol_, angtol,
					     grown_regions, adj_surfs, adj_edgs);
	if (grown_regions.size() > 0 || adj_surfs.size() > 0)
	  updateRegionsAndSurfaces(ki, grown_regions, adj_surfs);
	for (size_t kr=0; kr<adj_edgs.size(); ++kr)
	  {
	    size_t kj;
	    for (kj=0; kj<edges_.size(); ++kj)
	      if (edges_[kj].get() == adj_edgs[kr])
		break;
	    if (kj < edges_.size())
	      edges_.erase(edges_.begin()+kj);
	  }
	}
    }
#ifdef DEBUG
  checkConsistence("Regions10");

   std::cout << "Finished merge adjacent regions, regions: " << regions_.size() << ", surfaces: " << surfaces_.size() << std::endl;
   if (regions_.size() > 0)
    {
      std::cout << "Regions10" << std::endl;
      std::ofstream of("regions10.g2");
      std::ofstream ofm("mid_regions10.g2");
      std::ofstream ofs("small_regions10.g2");
      writeRegionStage(of, ofm, ofs);
      std::ofstream of0("regions10_helix.g2");
      for (int ki=0; ki<(int)regions_.size(); ++ki)
	{
	  if (regions_[ki]->getSurfaceFlag() == PROBABLE_HELIX)
	    {
	      regions_[ki]->writeRegionInfo(of0);
	      if (regions_[ki]->hasSurface())
		regions_[ki]->writeSurface(of0);
	    }
	}
        if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf10.g2");
	  writeRegionWithSurf(of);
	}
     }
   
   if (edges_.size() > 0)
     {
       std::ofstream ofe("edges10.g2");
       writeEdgeStage(ofe);
     }
#endif
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

#ifdef DEBUG_CHECK
  for (int ka=0; ka<(int)regions_.size(); ++ka)
    {
      std::set<RevEngPoint*> tmpset(regions_[ka]->pointsBegin(), regions_[ka]->pointsEnd());
      if ((int)tmpset.size() != regions_[ka]->numPoints())
	std::cout << "Point number mismatch, pre grow. " << ka << " " << tmpset.size() << " " << regions_[ka]->numPoints() << std::endl;

      bool con = regions_[ka]->isConnected();
      if (!con)
	std::cout << "Disconnected region, ka= " << ka << " " << regions_[ka].get() << std::endl;
    }
#endif

  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      if (regions_[ki]->hasSurface() && regions_[ki]->getSurfaceFlag() < ACCURACY_POOR)
	{
#ifdef DEBUG_GROW
      std::cout << "ki=" << ki << ", nmb reg: " << regions_.size() << ", nmb surf: " << surfaces_.size() << std::endl;
#endif
      growSurface(ki, pass);
#ifdef DEBUG_CHECK
  if (!regions_[ki]->isConnected())
    std::cout << "Disconnected region (grow), ki= " << ki << " " << regions_[ki].get() << std::endl;
#endif
	}
      // for (int kh=0; kh<(int)surfaces_.size(); ++kh)
      // 	{
      // 	  int numreg = surfaces_[kh]->numRegions();
      // 	  for (int ka=0; ka<numreg; ++ka)
      // 	    {
      // 	      RevEngRegion *reg = surfaces_[kh]->getRegion(ka);
      // 	      size_t kr;
      // 	      for (kr=0; kr<regions_.size(); ++kr)
      // 		if (reg == regions_[kr].get())
      // 		  break;
      // 	      if (kr == regions_.size())
      // 		std::cout << "Region4, surface 1. Obsolete region pointer, ki=" << ki << ", kh=" << kh << ". Region: " << reg << ", surface: " << surfaces_[kh].get() << std::endl;
      // 	    }
      // 	}
    }
      
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

#ifdef DEBUG
  checkConsistence("Regions11");

   std::cout << "Finished grow with surf, regions: " << regions_.size() << ", surfaces: " << surfaces_.size() << std::endl;
   if (regions_.size() > 0)
    {
      std::cout << "Regions11" << std::endl;
      std::ofstream of("regions11.g2");
      std::ofstream ofm("mid_regions11.g2");
      std::ofstream ofs("small_regions11.g2");
      writeRegionStage(of, ofm, ofs);
      std::ofstream of0("regions11_helix.g2");
      for (int ki=0; ki<(int)regions_.size(); ++ki)
	{
	  if (regions_[ki]->getSurfaceFlag() == PROBABLE_HELIX)
	    {
	      regions_[ki]->writeRegionInfo(of0);
	      if (regions_[ki]->hasSurface())
		regions_[ki]->writeSurface(of0);
	    }
	}
        if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf11.g2");
	  writeRegionWithSurf(of);
	}
     }
   
   if (edges_.size() > 0)
     {
       std::ofstream ofe("edges11.g2");
       writeEdgeStage(ofe);
     }
#endif
}

//===========================================================================
void RevEng::manageBlends1()
//===========================================================================
{
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

  
#ifdef DEBUG
  std::ofstream ofu1("unresolved3.g2");
  for (size_t kr=0; kr<regions_.size(); ++kr)
    {
      if (regions_[kr]->hasSurface())
	continue;
      if (regions_[kr]->hasAssociatedBlend())
	continue;
      regions_[kr]->writeRegionPoints(ofu1);
    }
#endif

#ifdef DEBUG
  if (edges_.size() > 0)
    {
      std::ofstream ofe("edges10.g2");
      writeEdgeStage(ofe);
    }
#endif
  
#ifdef DEBUG
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      int num = regions_[ki]->numPoints();
      for (int ka=0; ka<num; ++ka)
	if (regions_[ki]->getPoint(ka)->region() != regions_[ki].get())
	  std::cout << "Inconsistent region pointers, pre createBlendSurface: " << ki << " " << ka << std::endl;
    }
#endif
  
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

#ifdef DEBUG
   std::cout << "Create blends" << std::endl;
#endif
   for (size_t ki=0; ki<edges_.size(); ++ki)
     (void)createBlendSurface((int)ki);

   
#ifdef DEBUG
    std::cout << "Finished create blends, regions: " << regions_.size() << ", surfaces: " << surfaces_.size() << std::endl;
   if (regions_.size() > 0)
    {
      std::cout << "Regions11_2" << std::endl;
      std::ofstream of("regions11_2.g2");
      std::ofstream ofm("mid_regions11_2.g2");
      std::ofstream ofs("small_regions11_2.g2");
      writeRegionStage(of, ofm, ofs);
      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf11_2.g2");
	  writeRegionWithSurf(of);
	}
     }
   
   if (edges_.size() > 0)
     {
       std::ofstream ofe("edges11_2.g2");
       writeEdgeStage(ofe);
     }
#endif

   // Update blend radii
   equalizeBlendRadii();
   
#ifdef DEBUG
   if (regions_.size() > 0)
    {
      std::cout << "Regions11_3_1" << std::endl;
      std::ofstream of("regions11_3_1.g2");
      std::ofstream ofm("mid_regions11_3_1.g2");
      std::ofstream ofs("small_regions11_3_1.g2");
      writeRegionStage(of, ofm, ofs);
      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf11_4.g2");
	  writeRegionWithSurf(of);
	}
    }
#endif
   int stop_break = 1;
}

//===========================================================================
void RevEng::manageBlends2()
//===========================================================================
{
#ifdef DEBUG
   std::cout << "Set blend boundaries" << std::endl;
#endif
   for (size_t ki=0; ki<surfaces_.size(); ++ki)
     {
       int nreg = surfaces_[ki]->numRegions();
       if (nreg != 1)
	 continue;  // Not an issue currently
       RevEngRegion *reg = surfaces_[ki]->getRegion(0);
       if (reg->hasBlendEdge())
	 setBlendBoundaries(reg);
     }
   
#ifdef DEBUG
   //std::cout << "Finished set blend boundaries, regions: " << regions_.size() << ", surfaces: " << surfaces_.size() << std::endl;
   if (regions_.size() > 0)
    {
      std::cout << "Regions11_3" << std::endl;
      std::ofstream of("regions11_3.g2");
      std::ofstream ofm("mid_regions11_3.g2");
      std::ofstream ofs("small_regions11_3.g2");
      writeRegionStage(of, ofm, ofs);
      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf11_3.g2");
	  writeRegionWithSurf(of);
	}

      vector<shared_ptr<ftEdge> > trim_edgs;
      for (size_t kr=0; kr<regions_.size(); ++kr)
	{
	  // if (regions_[kr]->numPoints() == 0)
	  //   std::cout << "Finished set blend boundaries, empty region, ki=" << kr << ", region: " << regions_[kr].get() << std::endl;
	  int num = regions_[kr]->numTrimEdges();
	  if (num > 0)
	    {
	      vector<shared_ptr<ftEdge> > curr_edgs = regions_[kr]->getTrimEdges();
	      trim_edgs.insert(trim_edgs.end(), curr_edgs.begin(), curr_edgs.end());
	    }
	}
      if (trim_edgs.size() > 0)
	{
	  std::ofstream oft("trim_edges11_3.g2");
	  for (size_t kr=0; kr<trim_edgs.size(); ++kr)
	    {
	      shared_ptr<ParamCurve> tmp = trim_edgs[kr]->geomCurve();
	      shared_ptr<CurveOnSurface> tmp2 =
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp);
#ifdef DEBUG_BLEND
	      bool same_orient = tmp2->sameOrientation();
	      bool same_trace = tmp2->sameTrace(approx_tol_);
	      bool same_cv = tmp2->sameCurve(approx_tol_);
	      if ((!same_orient) || (!same_trace) || (!same_cv))
		std::cout << "Surface curve mismatch " << kr << " " << same_orient << " " << same_trace << " " << same_cv << std::endl;
#endif
	      shared_ptr<ParamCurve> tmp3 = tmp2->spaceCurve();
	      tmp3->writeStandardHeader(oft);
	      tmp3->write(oft);
	    }
	}
     }
   
   if (edges_.size() > 0)
     {
       std::ofstream ofe("edges11_3.g2");
       writeEdgeStage(ofe);
     }
#endif

   for (int ka=0; ka<(int)regions_.size(); ++ka)
     {
       if (regions_[ka]->hasSurface() && regions_[ka]->hasBlendEdge())
   	 growBlendSurface(ka);
     }
   
#ifdef DEBUG
    std::cout << "Finished grow blend surface, regions: " << regions_.size() << ", surfaces: " << surfaces_.size() << std::endl;
   if (regions_.size() > 0)
    {
      std::cout << "Regions11_3_2" << std::endl;
      std::ofstream of("regions11_3_2.g2");
      std::ofstream ofm("mid_regions11_3_2.g2");
      std::ofstream ofs("small_regions11_3_2.g2");
      writeRegionStage(of, ofm, ofs);
      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf11_3_2.g2");
	  writeRegionWithSurf(of);
	}

      vector<shared_ptr<ftEdge> > trim_edgs;
      for (size_t kr=0; kr<regions_.size(); ++kr)
	{
	  if (regions_[kr]->numPoints() == 0)
	    std::cout << "Finished grow blend surface, empty region, ki=" << kr << ", region: " << regions_[kr].get() << std::endl;
	  int num = regions_[kr]->numTrimEdges();
	  if (num > 0)
	    {
	      vector<shared_ptr<ftEdge> > curr_edgs = regions_[kr]->getTrimEdges();
	      trim_edgs.insert(trim_edgs.end(), curr_edgs.begin(), curr_edgs.end());
	    }
	}
      if (trim_edgs.size() > 0)
	{
	  std::ofstream oft("trim_edges11_3_2.g2");
	  for (size_t kr=0; kr<trim_edgs.size(); ++kr)
	    {
	      shared_ptr<ParamCurve> tmp = trim_edgs[kr]->geomCurve();
	      shared_ptr<CurveOnSurface> tmp2 =
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp);
#ifdef DEBUG_BLEND
	      bool same_orient = tmp2->sameOrientation();
	      bool same_trace = tmp2->sameTrace(approx_tol_);
	      bool same_cv = tmp2->sameCurve(approx_tol_);
	      if ((!same_orient) || (!same_trace) || (!same_cv))
		std::cout << "Surface curve mismatch " << kr << " " << same_orient << " " << same_trace << " " << same_cv << std::endl;
#endif
	      shared_ptr<ParamCurve> tmp3 = tmp2->spaceCurve();
	      tmp3->writeStandardHeader(oft);
	      tmp3->write(oft);
	    }
	}
     }
   
   if (edges_.size() > 0)
     {
       std::ofstream ofe("edges11_3_2.g2");
       writeEdgeStage(ofe);
     }
#endif

    for (int ka=0; ka<(int)regions_.size(); ++ka)
     {
       if (regions_[ka]->hasSurface() && regions_[ka]->numTrimEdges() > 0 &&
	   (!regions_[ka]->hasBlendEdge()))
    	 {
    	   vector<vector<RevEngPoint*> > added_groups;
    	   vector<HedgeSurface*> dummy_surfs;
    	   double tol = 1.5*approx_tol_;
    	   double angtol = 5.0*anglim_;
    	   regions_[ka]->removeLowAccuracyPoints(min_point_region_, 
    						 tol, angtol, added_groups);
    	   if (added_groups.size() > 0)
    	     surfaceExtractOutput(ka, added_groups, dummy_surfs);
    	 }
     }
    
    for (int ka=0; ka<(int)regions_.size(); ++ka)
     {
       if (regions_[ka]->hasSurface() && regions_[ka]->numTrimEdges() > 0)
	 growMasterSurface(ka);
     }
   
#ifdef DEBUG
    std::cout << "Finished grow blend surface, regions: " << regions_.size() << ", surfaces: " << surfaces_.size() << std::endl;
   if (regions_.size() > 0)
    {
      std::cout << "Regions11_3_3" << std::endl;
      std::ofstream of("regions11_3_3.g2");
      std::ofstream ofm("mid_regions11_3_3.g2");
      std::ofstream ofs("small_regions11_3_3.g2");
      writeRegionStage(of, ofm, ofs);
      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf11_3_3.g2");
	  writeRegionWithSurf(of);
	}

      vector<shared_ptr<ftEdge> > trim_edgs;
      for (size_t kr=0; kr<regions_.size(); ++kr)
	{
	  if (regions_[kr]->numPoints() == 0)
	    std::cout << "Grow master surface, empty region, ki=" << kr << ", region: " << regions_[kr].get() << std::endl;
	  int num = regions_[kr]->numTrimEdges();
	  if (num > 0)
	    {
	      vector<shared_ptr<ftEdge> > curr_edgs = regions_[kr]->getTrimEdges();
	      trim_edgs.insert(trim_edgs.end(), curr_edgs.begin(), curr_edgs.end());
	    }
	}
      if (trim_edgs.size() > 0)
	{
	  std::ofstream oft("trim_edges11_3_3.g2");
	  for (size_t kr=0; kr<trim_edgs.size(); ++kr)
	    {
	      shared_ptr<ParamCurve> tmp = trim_edgs[kr]->geomCurve();
	      shared_ptr<CurveOnSurface> tmp2 =
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp);
#ifdef DEBUG_BLEND
	      bool same_orient = tmp2->sameOrientation();
	      bool same_trace = tmp2->sameTrace(approx_tol_);
	      bool same_cv = tmp2->sameCurve(approx_tol_);
	      if ((!same_orient) || (!same_trace) || (!same_cv))
		std::cout << "Surface curve mismatch " << kr << " " << same_orient << " " << same_trace << " " << same_cv << std::endl;
#endif
	      shared_ptr<ParamCurve> tmp3 = tmp2->spaceCurve();
	      tmp3->writeStandardHeader(oft);
	      tmp3->write(oft);
	    }
	}
     }
   
   if (edges_.size() > 0)
     {
       std::ofstream ofe("edges11_3_3.g2");
       writeEdgeStage(ofe);
     }
#endif

#ifdef DEBUG_BLEND
   std::ofstream ofbb("blend_branch.g2");
   vector<RevEngPoint*> bbpts;
   for (size_t kr=0; kr<regions_.size(); ++kr)
     {
       if (!regions_[kr]->hasBlendEdge())
	 continue;
       vector<RevEngPoint*> currbb = regions_[kr]->extractBranchPoints();
       if (currbb.size() > 0)
	 bbpts.insert(bbpts.end(), currbb.begin(), currbb.end());
     }
   if (bbpts.size() > 0)
     {
       ofbb << "400 1 0 4 0 0 0 255" << std::endl;
       ofbb << bbpts.size() << std::endl;
       for (size_t kr=0; kr<bbpts.size(); ++kr)
	 ofbb << bbpts[kr]->getPoint() << std::endl;
     }
#endif

  // TESTING. Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

#ifdef DEBUG
   std::cout << "Torus corners" << std::endl;
#endif

   for (size_t ki=0; ki<regions_.size(); ++ki)
     {
       if (regions_[ki]->toBeRemoved())
	 continue;
       
      if (!regions_[ki]->hasSurface())
	continue;

      if (!regions_[ki]->hasBlendEdge())
	continue;

      bool done = defineTorusCorner(ki);
#ifdef DEBUG_BLEND
      std::ofstream ofsfs("curr_sfs.g2");
      for (size_t kj=0; kj<regions_.size(); ++kj)
	{
	  if (regions_[kj]->hasSurface())
	    {
	      shared_ptr<ParamSurface> curr_sf = regions_[kj]->getSurface(0)->surface();
	      curr_sf->writeStandardHeader(ofsfs);
	      curr_sf->write(ofsfs);
	    }
	}
      int stop_tor = 1;
#endif
     }
   vector<RevEngRegion*> removereg;
   vector<HedgeSurface*> removehedge;
   for (size_t ki=0; ki<regions_.size(); ++ki)
     {
       if (regions_[ki]->toBeRemoved())
	 {
	   removereg.push_back(regions_[ki].get());
	   int num_sfs = regions_[ki]->numSurface();
	   for (int ka=0; ka<num_sfs; ++ka)
	     removehedge.push_back(regions_[ki]->getSurface(ka));
	 }
     }
   if (removereg.size() > 0)
     {
      int dummy_ix = 0;
      updateRegionsAndSurfaces(dummy_ix, removereg, removehedge);
     }
   
#ifdef DEBUG
    std::cout << "Finished torus corners, regions: " << regions_.size() << ", surfaces: " << surfaces_.size() << std::endl;
   if (regions_.size() > 0)
    {
      std::cout << "Regions11_4" << std::endl;
      std::ofstream of("regions11_4.g2");
      std::ofstream ofm("mid_regions11_4.g2");
      std::ofstream ofs("small_regions11_4.g2");
      writeRegionStage(of, ofm, ofs);
      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf11_4.g2");
	  writeRegionWithSurf(of);
	}
      
      vector<shared_ptr<ftEdge> > trim_edgs;
      for (size_t kr=0; kr<regions_.size(); ++kr)
	{
	  // if (regions_[kr]->numPoints() == 0)
	  //   std::cout << "Torus corner, ki=" << kr << ", region: " << regions_[kr].get() << std::endl;
	  int num = regions_[kr]->numTrimEdges();
	  if (num > 0)
	    {
	      vector<shared_ptr<ftEdge> > curr_edgs = regions_[kr]->getTrimEdges();
	      trim_edgs.insert(trim_edgs.end(), curr_edgs.begin(), curr_edgs.end());
	    }
	}
      if (trim_edgs.size() > 0)
	{
	  std::ofstream oft("trim_edges11_4.g2");
	  for (size_t kr=0; kr<trim_edgs.size(); ++kr)
	    {
	      shared_ptr<ParamCurve> tmp = trim_edgs[kr]->geomCurve();
	      shared_ptr<CurveOnSurface> tmp2 =
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp);
#ifdef DEBUG_BLEND
	      bool same_orient = tmp2->sameOrientation();
	      bool same_trace = tmp2->sameTrace(approx_tol_);
	      bool same_cv = tmp2->sameCurve(approx_tol_);
	      if ((!same_orient) || (!same_trace) || (!same_cv))
		std::cout << "Surface curve mismatch " << kr << " " << same_orient << " " << same_trace << " " << same_cv << std::endl;
#endif
	      shared_ptr<ParamCurve> tmp3 = tmp2->spaceCurve();
	      tmp3->writeStandardHeader(oft);
	      tmp3->write(oft);
	    }
	}
    }
#endif
   // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

  vector<RevEngRegion*> cand_corner_adj;
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      if (!regions_[ki]->hasBlendEdge())
	continue;
      if (regions_[ki]->getBlendEdge()->isClosed(approx_tol_))
	continue;
      int num_edg = regions_[ki]->numTrimEdges();
      if (num_edg < 4)
	cand_corner_adj.push_back(regions_[ki].get());
    }

#ifdef DEBUG_BLEND
  std::ofstream ofmb0("adj_candidate_blend_corner.g2");
  for (size_t ki=0; ki<cand_corner_adj.size(); ++ki)
    {
      cand_corner_adj[ki]->writeRegionPoints(ofmb0);
      cand_corner_adj[ki]->writeSurface(ofmb0);
    }
#endif
  if (cand_corner_adj.size() > 2)
    defineMissingCorner(cand_corner_adj);
   
   removereg.clear();
   removehedge.clear();
   for (size_t ki=0; ki<regions_.size(); ++ki)
     {
       if (regions_[ki]->toBeRemoved())
	 {
	   removereg.push_back(regions_[ki].get());
	   int num_sfs = regions_[ki]->numSurface();
	   for (int ka=0; ka<num_sfs; ++ka)
	     removehedge.push_back(regions_[ki]->getSurface(ka));
	 }
     }
   if (removereg.size() > 0)
     {
      int dummy_ix = 0;
      updateRegionsAndSurfaces(dummy_ix, removereg, removehedge);
     }

#ifdef DEBUG
    std::cout << "Finished missing corners, regions: " << regions_.size() << ", surfaces: " << surfaces_.size() << std::endl;
   if (regions_.size() > 0)
    {
      std::cout << "Regions11_5" << std::endl;
      std::ofstream of("regions11_5.g2");
      std::ofstream ofm("mid_regions11_5.g2");
      std::ofstream ofs("small_regions11_5.g2");
      writeRegionStage(of, ofm, ofs);
      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf11_5.g2");
	  writeRegionWithSurf(of);
	}
      
      vector<shared_ptr<ftEdge> > trim_edgs;
      for (size_t kr=0; kr<regions_.size(); ++kr)
	{
	  // if (regions_[kr]->numPoints() == 0)
	  //   std::cout << "Missing corner, ki=" << kr << ", region: " << regions_[kr].get() << std::endl;
	  int num = regions_[kr]->numTrimEdges();
	  if (num > 0)
	    {
	      vector<shared_ptr<ftEdge> > curr_edgs = regions_[kr]->getTrimEdges();
	      trim_edgs.insert(trim_edgs.end(), curr_edgs.begin(), curr_edgs.end());
	    }
	}
      if (trim_edgs.size() > 0)
	{
	  std::ofstream oft("trim_edges11_5.g2");
	  for (size_t kr=0; kr<trim_edgs.size(); ++kr)
	    {
	      shared_ptr<ParamCurve> tmp = trim_edgs[kr]->geomCurve();
	      shared_ptr<CurveOnSurface> tmp2 =
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp);
#ifdef DEBUG_BLEND
	      bool same_orient = tmp2->sameOrientation();
	      bool same_trace = tmp2->sameTrace(approx_tol_);
	      bool same_cv = tmp2->sameCurve(approx_tol_);
	      if ((!same_orient) || (!same_trace) || (!same_cv))
		std::cout << "Surface curve mismatch " << kr << " " << same_orient << " " << same_trace << " " << same_cv << std::endl;
#endif
	      shared_ptr<ParamCurve> tmp3 = tmp2->spaceCurve();
	      tmp3->writeStandardHeader(oft);
	      tmp3->write(oft);
	    }
	}
    }
#endif
      
  int stop_break_blend = 1;
 
}


//===========================================================================
void RevEng::equalizeBlendRadii()
//===========================================================================
{
  // Equialize adjacent blends between the same surfaces
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      if (!regions_[ki]->hasSurface())
	continue;
      
      if (!regions_[ki]->hasBlendEdge())
	continue;

      RevEngEdge *edge1 = regions_[ki]->getBlendEdge();
      RevEngRegion* adj1[2];
      edge1->getAdjacent(adj1[0], adj1[1]);

      
      for (size_t kj=ki+1; kj<regions_.size(); ++kj)
	{
	  if (!regions_[kj]->hasSurface())
	    continue;
      
	  if (!regions_[kj]->hasBlendEdge())
	    continue;
	  
	  RevEngEdge *edge2 = regions_[kj]->getBlendEdge();
	  RevEngRegion* adj2[2];
	  edge2->getAdjacent(adj2[0], adj2[1]);

	  if ((adj1[0] == adj2[0] || adj1[0] == adj2[1]) &&
	      (adj1[1] == adj2[0] || adj1[1] == adj2[1]))
	    {
#ifdef DEBUG_BLEND
	      std::ofstream of("adj_blend.g2");
	      regions_[ki]->writeRegionPoints(of);
	      regions_[kj]->writeRegionPoints(of);
#endif
	      double par1, par2;
	      bool is_adjacent = edge1->isAdjacent(edge2, approx_tol_, par1,
						   par2);
	      if (is_adjacent)
		equalizeAdjacent(ki, kj);
	    }
	}
  
    }
  
  // Collect information about blend radii
  vector<double> rad;
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      if (!regions_[ki]->hasSurface())
	continue;
      
      if (!regions_[ki]->hasBlendEdge())
	continue;

      shared_ptr<ParamSurface> surf = regions_[ki]->getSurface(0)->surface();
      shared_ptr<ElementarySurface> elem =
	dynamic_pointer_cast<ElementarySurface, ParamSurface>(surf);
      if (!elem.get())
	continue;

      double radius = elem->radius(0.0, 0.0);
      double radius2 = elem->radius2(0.0, 0.0);
      if (radius2 > 0)
	rad.push_back(radius2);
      else
	rad.push_back(radius);
    }

  if (rad.size() == 0)
    return;
  
  std::sort(rad.begin(), rad.end());
  vector<size_t> ixs;
  ixs.push_back(0);
  double tmean = rad[0];
  size_t prev = 0;
  int nn = 1;
  size_t ki=1;
  double fac = 0.5;
  double curr_mean, range;
  for (; ki<rad.size(); ++ki)
    {
      curr_mean = (tmean + rad[ki])/(double)(nn+1);
      range = rad[ki]-rad[prev];
      if (range < fac*curr_mean)
	{
	  tmean += rad[ki];
	  ++nn;
	}
      else
	{
	  ixs.push_back(ki);
	  prev = ki;
	  tmean = rad[ki];
	  nn = 1;
	}
      int stop_break0 = 1;
    }

  vector<pair<double,double> > rad_range;
  for (ki=1; ki<ixs.size(); ++ki)
    rad_range.push_back(std::make_pair(rad[ixs[ki-1]], rad[ixs[ki]-1]));
  rad_range.push_back(std::make_pair(rad[ixs[ixs.size()-1]], rad[rad.size()-1]));

  vector<double> mean_rad(rad_range.size());
  for (ki=0; ki<rad_range.size(); ++ki)
    {
      mean_rad[ki] = 0.5*(rad_range[ki].first + rad_range[ki].second);
      if (mean_rad[ki] <= 0.1)
	{
	  double low = 0.01*(double)((int)(100.0*mean_rad[ki]));
	  double high = 0.01*(double)((int)(100.0*mean_rad[ki]+1));
	  mean_rad[ki] = (mean_rad[ki]-low < high-mean_rad[ki]) ? low : high;
	}
      else if (mean_rad[ki] <= 1.0)
	{
	  double low = 0.1*(double)((int)(10.0*mean_rad[ki]));
	  double high = 0.1*(double)((int)(10.0*mean_rad[ki]+1));
	  mean_rad[ki] = (mean_rad[ki]-low < high-mean_rad[ki]) ? low : high;
	}
      else if (mean_rad[ki] <= 10.0)
	{
	  double low = (double)((int)(mean_rad[ki]));
	  double high = (double)((int)(mean_rad[ki]+1));
	  double mid = 0.5*(low+high);
	  mean_rad[ki] = (mean_rad[ki]-low < std::min(fabs(mean_rad[ki]-mid), high-mean_rad[ki])) ? 
	    low : ((high-mean_rad[ki]) < fabs(mean_rad[ki]-mid) ? high : mid);;
	}
      else
	{
	  double low = (double)((int)(mean_rad[ki]));
	  double high = (double)((int)(mean_rad[ki]+1));
	  mean_rad[ki] = (mean_rad[ki]-low < high-mean_rad[ki]) ? low : high;
	}
    }
  
  for (ki=0; ki<regions_.size(); ++ki)
    {
      if (!regions_[ki]->hasSurface())
	continue;
      
      if (!regions_[ki]->hasBlendEdge())
	continue;

      shared_ptr<ParamSurface> surf = regions_[ki]->getSurface(0)->surface();
      shared_ptr<ElementarySurface> elem =
	dynamic_pointer_cast<ElementarySurface, ParamSurface>(surf);
      if (!elem.get())
	continue;

      double radius = elem->radius(0.0, 0.0);
      double radius2 = elem->radius2(0.0, 0.0);
      if (radius2 > 0.0)
	radius = radius2;

      size_t kj;
      for (kj=0; kj<rad_range.size(); ++kj)
	if (radius >= rad_range[kj].first && radius <= rad_range[kj].second)
	  {
	    updateBlendRadius(ki, mean_rad[kj]);
	    break;
	  }
    }

    int stop_break = 1;

}

//===========================================================================
void RevEng::equalizeAdjacent(size_t ix1, size_t ix2)
//===========================================================================
{
  shared_ptr<ParamSurface> surf1 = regions_[ix1]->getSurface(0)->surface();
  shared_ptr<ParamSurface> surf2 = regions_[ix2]->getSurface(0)->surface();
  if ((!surf1) || (!surf2))
    return;

  if (surf1->instanceType() != surf2->instanceType())
    return;

  double angtol = 5.0*anglim_;
  shared_ptr<ElementarySurface> elem1 =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf1);
  shared_ptr<ElementarySurface> elem2 =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf2);
  Point axis1 = elem1->direction();
  Point axis2 = elem1->direction();
  double ang = axis1.angle(axis2);
  ang = std::min(ang, M_PI-ang);
  if (ang > angtol)
    return;
  Point pos1 = elem1->location();
  Point pos2 = elem2->location();
  double rad1 = elem1->radius(0.0, 0.0);
  double rad2 = elem2->radius(0.0, 0.0);
  double fac = 0.1;
  double fac2 = 0.25;
  if (fabs(rad1-rad2) > fac*std::max(rad1, rad2))
    return;  // Too large difference
  double rad = 0.5*(rad1 + rad2);
  
  if (axis1*axis2 < 0.0)
    axis2 *= -1;
  Point axis = 0.5*(axis1 + axis2);
  axis.normalize();

  Point Cx;
  Point Cx1 = elem1->direction2();
  Point Cx2 = elem2->direction2();
  double ang2 = Cx1.angle(Cx2);
  if (ang2 <= angtol || M_PI-ang2 <= angtol)
    {
      if (Cx1*Cx2 < 0.0)
	Cx2 *= -1;
      Cx = 0.5*(Cx1 + Cx2);
      Point Cy = axis.cross(Cx);
      Cx = Cy.cross(axis);
      Cx.normalize();
    }
  else
    {
      int ka = -1;
      double minang = std::numeric_limits<double>::max();
      for (int kb=0; kb<3; ++kb)
	{
	  double ang3 = mainaxis_[kb].angle(axis);
	  ang3 = std::min(ang3, M_PI-ang3);
	  if (ang3 < minang)
	    {
	      minang = ang3;
	      ka = kb;
	    }
	}
      if (ka < 0)
	return;
      int kb = (ka > 0) ? ka - 1 : 2;
      Cx = mainaxis_[kb].cross(axis);
    }

  RectDomain dom1 = elem1->getParameterBounds();
  RectDomain dom2 = elem2->getParameterBounds();
  shared_ptr<ElementarySurface> upd1, upd2;
  bool cyllike = true;
  if (surf1->instanceType() == Class_Cylinder)
    {
      Point pos2_2 = pos1 + ((pos2 - pos1)*axis)*axis;
      if (pos2.dist(pos2_2) > approx_tol_)
	return;

     Point pos1_2 = pos2 + ((pos1 - pos2)*axis)*axis;
      if (pos1.dist(pos1_2) > approx_tol_)
	return;

      Point pos3 = 0.5*(pos1_2 + pos2_2);  // Point on updated axis
      Point pos3_1 = pos3 + ((pos1 - pos3)*axis)*axis;
      Point pos3_2 = pos3 + ((pos2 - pos3)*axis)*axis;

      double rad = 0.5*(rad1 + rad2);
      upd1 = shared_ptr<Cylinder>(new Cylinder(rad, pos3_1, axis, Cx));
      upd2 = shared_ptr<Cylinder>(new Cylinder(rad, pos3_2, axis, Cx));

      double upar1 = 0.5*(dom1.umin() + dom2.umin());  // Could be a problem if
      // the direction of Cx is changed significantly
      double upar2 = 0.5*(dom1.umax() + dom2.umax());
      upd1->setParameterBounds(upar1, dom1.vmin(), upar2, dom1.vmax());
      upd2->setParameterBounds(upar1, dom2.vmin(), upar2, dom2.vmax());
      cyllike = true;
    }
  else if (surf1->instanceType() == Class_Torus)
    {
      Point pos2_2 = pos1 + ((pos2 - pos1)*axis)*axis;
      if (pos2_2.dist(pos2) > approx_tol_)
	return;

      double minrad1 = elem1->radius2(0.0, 0.0);
      double minrad2 = elem2->radius2(0.0, 0.0);
      if (fabs(minrad1-minrad2) > fac2*std::max(minrad1, minrad2))
	return;

      Point centre = 0.5*(pos1 + pos2);
      double minrad = 0.5*(minrad1 + minrad2);
      upd1 = shared_ptr<Torus>(new Torus(rad, minrad, centre, axis, Cx));
      upd2 = shared_ptr<Torus>(new Torus(rad, minrad, centre, axis, Cx));
      double vpar1 = 0.5*(dom1.vmin() + dom2.vmin());
      double vpar2 = 0.5*(dom1.vmax() + dom2.vmax());
      upd1->setParameterBounds(dom1.umin(), vpar1, dom1.umax(), vpar2);
      upd2->setParameterBounds(dom2.umin(), vpar1, dom2.umax(), vpar2);
      cyllike = false;
    }
  else
    return;  // Not supported
#ifdef DEBUG_BLEND
  std::ofstream of("updated_adjacent.g2");
  upd1->writeStandardHeader(of);
  upd1->write(of);
  upd2->writeStandardHeader(of);
  upd2->write(of);
#endif

  // Parameterize points
  double maxd1, avd1;
  int nmb_in1, nmb2_in1;
  vector<RevEngPoint*> in1, out1;
  vector<pair<double,double> > dist_ang1;
  vector<double> parvals1;
  RevEngUtils::distToSurf(regions_[ix1]->pointsBegin(), regions_[ix1]->pointsEnd(),
			  upd1, approx_tol_, maxd1, avd1, nmb_in1, nmb2_in1, in1, out1,
			  parvals1, dist_ang1, angtol);
  int sf_flag1 = regions_[ix1]->defineSfFlag(0, approx_tol_, nmb_in1, nmb2_in1, avd1,
					     cyllike);
  int num_pts1 = regions_[ix1]->numPoints();
  for (int ka=0; ka<num_pts1; ++ka)
    {
      RevEngPoint *curr = regions_[ix1]->getPoint(ka);
      curr->setPar(Vector2D(parvals1[2*ka],parvals1[2*ka+1]));
      curr->setSurfaceDist(dist_ang1[ka].first, dist_ang1[ka].second);
    }
  regions_[ix1]->updateInfo(approx_tol_, angtol);
  regions_[ix1]->setSurfaceFlag(sf_flag1);

  // Replace Surface
  regions_[ix1]->getSurface(0)->replaceSurf(upd1);
  
  double maxd2, avd2;
  int nmb_in2, nmb2_in2;
  vector<RevEngPoint*> in2, out2;
  vector<pair<double,double> > dist_ang2;
  vector<double> parvals2;
  RevEngUtils::distToSurf(regions_[ix2]->pointsBegin(), regions_[ix2]->pointsEnd(),
			  upd1, approx_tol_, maxd2, avd2, nmb_in2, nmb2_in2, in2, out2,
			  parvals2, dist_ang2, angtol);
  int sf_flag2 = regions_[ix2]->defineSfFlag(0, approx_tol_, nmb_in2, nmb2_in2, avd2,
					     cyllike);
  int num_pts2 = regions_[ix2]->numPoints();
  for (int ka=0; ka<num_pts2; ++ka)
    {
      RevEngPoint *curr = regions_[ix2]->getPoint(ka);
      curr->setPar(Vector2D(parvals2[2*ka],parvals2[2*ka+1]));
      curr->setSurfaceDist(dist_ang2[ka].first, dist_ang2[ka].second);
    }
  regions_[ix2]->updateInfo(approx_tol_, angtol);
  regions_[ix2]->setSurfaceFlag(sf_flag2);

  // Replace Surface
  regions_[ix2]->getSurface(0)->replaceSurf(upd2);
  
}

//===========================================================================
void RevEng::updateBlendRadius(size_t ix, double radius)
//===========================================================================
{
  double angtol = 5.0*anglim_;
  double eps = 1.0e-6;
  RevEngEdge *edge = regions_[ix]->getBlendEdge();
  RevEngRegion* adj[2];
  edge->getAdjacent(adj[0], adj[1]);

  if ((!adj[0]->hasSurface()) || (!adj[1]->hasSurface()))
    return;  // Something wrong

  bool out1 = false, out2 = false;
  edge->getOuterInfo(out1, out2);

  // Intersection curve
  vector<shared_ptr<CurveOnSurface> > cvs;
  edge->getCurve(cvs);
  vector<Point> der(2);
  cvs[0]->point(der, 0.5*(cvs[0]->startparam()+cvs[0]->endparam()), 1);

  shared_ptr<ParamSurface> surf = regions_[ix]->getSurface(0)->surface();
  shared_ptr<Cylinder> init_cyl =
    dynamic_pointer_cast<Cylinder, ParamSurface>(surf);
  shared_ptr<Torus> init_tor =
    dynamic_pointer_cast<Torus, ParamSurface>(surf);
  if ((!init_cyl.get()) && (!init_tor.get()))
    return;
  shared_ptr<ParamSurface> surf1 = adj[0]->getSurface(0)->surface();
  shared_ptr<ParamSurface> surf2 = adj[1]->getSurface(0)->surface();
  shared_ptr<ElementarySurface> elem1 =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf1);
  shared_ptr<ElementarySurface> elem2 =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf2);
  shared_ptr<ElementarySurface> upd_blend;
  double ylen;
  double init_rad;
  Point dir1 = elem1->direction();
  Point norm1 = adj[0]->getMeanNormalTriang();
  if (elem1->instanceType() == Class_Plane && dir1*norm1 < 0.0)
    dir1 *= -1;
  Point dir2 = elem2->direction();
  Point norm2 = adj[1]->getMeanNormalTriang();
  if (elem2->instanceType() == Class_Plane && dir2*norm2 < 0.0)
    dir2 *= -1;
  if (init_cyl.get())
    {
      if (elem1->instanceType() != Class_Plane || elem2->instanceType() != Class_Plane)
	return; 

      Point lin1, lin2;
      Point dir1_2 = dir1, dir2_2 = dir2;
      if (elem1->instanceType() == Class_Plane)
	lin1 = der[1].cross(dir1);
      else
	{
	  double clo_u, clo_v, clo_dist;
	  Point clo;
	  surf1->closestPoint(der[0], clo_u, clo_v, clo, clo_dist, eps);
	  surf1->normal(dir1_2, clo_u, clo_v);
	  lin1 = der[1].cross(dir1_2);
	}
      
      if (elem2->instanceType() == Class_Plane)
	lin2 = der[1].cross(dir2);
      else
	{
	  double clo_u, clo_v, clo_dist;
	  Point clo;
	  surf2->closestPoint(der[0], clo_u, clo_v, clo, clo_dist, eps);
	  surf2->normal(dir2_2, clo_u, clo_v);
	  lin2 = der[1].cross(dir2_2);
	}
     
      // Create new cylinder
      Point axis = init_cyl->direction();
      Point Cx = init_cyl->direction2();
      init_rad = init_cyl->getRadius();
      lin1.normalize();
      if (lin1*dir2_2 < 0.0)
	lin1 *= -1;
      lin2.normalize();
      if (lin2*dir1_2 < 0.0)
	lin2 *= -1;

      Point vec = lin1 + lin2;
      vec.normalize();
      double alpha = lin1.angle(lin2);
      double fac = 1.0/sin(0.5*alpha);
      int sgn = (out1 && out2) ? 1 : -1;
      Point loc = der[0] + sgn*fac*radius*vec;
      double xlen = der[0].dist(loc);
      ylen = sqrt(xlen*xlen - radius*radius);

      upd_blend = shared_ptr<Cylinder>(new Cylinder(radius, loc, axis, Cx));
    }
  else if (init_tor.get())
    {
      // int plane_ix = 0;
      // if (elem2->instanceType() == Class_Plane)
      // 	{
      // 	  std::swap(elem1, elem2);
      // 	  plane_ix = 1;
      // 	}
      // if (elem1->instanceType() != Class_Plane)
      // 	return; 
      // if (elem2->instanceType() != Class_Cylinder &&
      // 	  elem2->instanceType() != Class_Cone)
      // 	return;
      
      // init_rad = init_tor->radius2(0.0, 0.0);

      // bool rot_out = (plane_ix == 0) ? out2 : out1;
      // bool plane_out = (plane_ix == 0) ? out1 : out2;
      // int sgn1 = plane_out ? -1 : 1;
      // int sgn2 = rot_out ? -1 : 1;
      // Point normal0 = elem1->direction();
      // Point norm1 = adj[plane_ix]->getMeanNormalTriang();
      // if (normal0*norm1 < 0.0)
      // 	sgn1 *= -1;

      bool plane1 = (elem1->instanceType() == Class_Plane);
      bool plane2 = (elem2->instanceType() == Class_Plane);
      shared_ptr<Cone> cone1 = dynamic_pointer_cast<Cone,ElementarySurface>(elem1);
      shared_ptr<Cone> cone2 = dynamic_pointer_cast<Cone,ElementarySurface>(elem2);
      shared_ptr<ElementarySurface> rotational;
      if (plane1)
	rotational = elem2;
      else if (plane2)
	rotational = elem1;
      else if (cone1.get())
	rotational = elem2;  // elem2 is expected to be a cylinder
      else if (cone2.get())
	rotational = elem1;
      double alpha1 = 0.0, alpha2 = 0.0;
      if (cone1.get())
	alpha1 = cone1->getConeAngle();
      if (cone2.get())
	alpha2 = cone2->getConeAngle();
      double alpha = fabs(alpha1) + fabs(alpha2);
      double beta = (plane1 || plane2) ? 0.5*M_PI + alpha : M_PI-fabs(alpha);
      double phi2 = 0.5*(M_PI - beta);
      Point axis = rotational->direction();
      Point loc = rotational->location();
      double rd = (der[0]-loc)*axis;
      Point centr = loc + rd*axis;
      double hh = init_tor->location().dist(centr);
      init_rad = init_tor->radius2(0.0, 0.0);
      double d2 = hh/sin(phi2) - init_rad;
      int sgn = 1;
      if ((elem1->instanceType() == Class_Plane && elem1->direction()*dir1 < 0.0) ||
	  (elem2->instanceType() == Class_Plane && elem2->direction()*dir2 < 0.0))
	sgn = -1;
      double Rrad;
      Point centre, normal, Cx;
      bool OK = getTorusParameters(elem1, elem2, der[0], radius, d2, out1, out2, 
				   sgn, Rrad, centre, normal, Cx);
      if (!OK)
	return;  // Don't change
     
      upd_blend = shared_ptr<Torus>(new Torus(Rrad, radius, centre, normal, Cx));

      // bool setbounds = false;
      // if ((plane1 && out2) || (plane2 && out1) ||
      // 	  (plane1 && false && plane2 == false &&
      // 	   ((cone1.get() && out1 == false) || (cone2.get() && out2 == false))))
      // 	setbounds = true;
	
      // if (setbounds)
      // 	{
	  RectDomain dom = init_tor->getParameterBounds();
	  upd_blend->setParameterBounds(dom.umin(), dom.vmin()-M_PI,
					dom.umax(), dom.vmax()-M_PI);
	// }
      double xlen = der[0].dist(centre);
      ylen = fabs(xlen - Rrad);
     }
  
  // RectDomain dom = init_cyl->getParameterBounds();
  // cyl->setParameterBounds(dom.umin(), dom.vmin(), dom.umax(), dom.vmax());
#ifdef DEBUG_BLEND
  std::ofstream of("blend2.g2");
  regions_[ix]->writeRegionPoints(of);
  surf->writeStandardHeader(of);
  surf->write(of);

  adj[0]->writeRegionPoints(of);
  adj[1]->writeRegionPoints(of);
  upd_blend->writeStandardHeader(of);
  upd_blend->write(of);
  RectDomain dom2 = upd_blend->getParameterBounds();
  vector<shared_ptr<ParamCurve> > tmp_cvs = upd_blend->constParamCurves(dom2.vmin(), true);
  tmp_cvs[0]->writeStandardHeader(of);
  tmp_cvs[0]->write(of);
#endif

  if (radius < init_rad)
    {
      // Move outside points from the blend surface to the adjacent
      // surfaces. If the new radius is larger than the initial one,
      // moving of points will be performed by the next function
      vector<RevEngPoint*> points(regions_[ix]->pointsBegin(), regions_[ix]->pointsEnd());
      vector<vector<RevEngPoint*> > out_pts(2);
      adj[0]->sortBlendPoints(points, cvs, ylen, adj[1], out_pts[0], out_pts[1]);
      vector<RevEngPoint*> all_out;
      if (out_pts[0].size() > 0)
	all_out.insert(all_out.end(), out_pts[0].begin(), out_pts[0].end());
      if (out_pts[1].size() > 0)
	all_out.insert(all_out.end(), out_pts[1].begin(), out_pts[1].end());
      if (all_out.size() > 0)
	regions_[ix]->removePoints(all_out);
      
      for (int ka=0; ka<2; ++ka)
	{
	  vector<RevEngPoint*> to_blend;  // Not here
	  adj[ka]->updateWithPointsInOut(to_blend, out_pts[ka], approx_tol_, angtol);
	}
    }

  // Parameterize remaining points
  double maxd, avd;
  int nmb_in, nmb2_in;
  vector<RevEngPoint*> in, out;
  vector<pair<double,double> > dist_ang;
  vector<double> parvals;
  RevEngUtils::distToSurf(regions_[ix]->pointsBegin(), regions_[ix]->pointsEnd(),
			  upd_blend, approx_tol_, maxd, avd, nmb_in, nmb2_in, in, out,
			  parvals, dist_ang, angtol);
  int sf_flag = regions_[ix]->defineSfFlag(0, approx_tol_, nmb_in, nmb2_in, avd, true);
  int num_pts = regions_[ix]->numPoints();
  for (int ka=0; ka<num_pts; ++ka)
    {
      RevEngPoint *curr = regions_[ix]->getPoint(ka);
      curr->setPar(Vector2D(parvals[2*ka],parvals[2*ka+1]));
      curr->setSurfaceDist(dist_ang[ka].first, dist_ang[ka].second);
    }
  regions_[ix]->setSurfaceFlag(sf_flag);

  // Replace Surface
  regions_[ix]->getSurface(0)->replaceSurf(upd_blend);
}


//===========================================================================
bool RevEng::defineTorusCorner(size_t ix)
//===========================================================================
{
  double pihalf = 0.5*M_PI;
  double angtol = 5.0*anglim_;
  double blendfac = 2.0;
  
  shared_ptr<ParamSurface> surf1 = regions_[ix]->getSurface(0)->surface();
  if (surf1->instanceType() != Class_Cylinder)
    return false;

  if (!regions_[ix]->hasBlendEdge())
    return false;

  // Bound cylinder in the length direction
  double diag = bbox_.low().dist(bbox_.high());
  if (!surf1->isBounded())
    {
      shared_ptr<Cylinder> cyl = dynamic_pointer_cast<Cylinder,ParamSurface>(surf1);
      cyl->setParamBoundsV(-0.5*diag, 0.5*diag);
    }

  RevEngEdge* blend_edge = regions_[ix]->getBlendEdge();
  RevEngRegion *adj1=0, *adj2=0;
  blend_edge->getAdjacent(adj1, adj2);
  
  vector<RevEngEdge*> revedg1 = regions_[ix]->getAllRevEdges();
  shared_ptr<ElementarySurface> elem1 =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf1);
  Point dir1 = elem1->direction();

#ifdef DEBUG_BLEND
  std::ofstream of0("pot_tor_adj.g2");
  surf1->writeStandardHeader(of0);
  surf1->write(of0);
#endif  
  // Looking for a plane
  bool found_torus = false;
  for (size_t kj=0; kj<regions_.size(); ++kj)
    {
      if (kj == ix)
	continue;
      if (regions_[kj]->toBeRemoved())
	continue;
      if (!regions_[kj]->hasSurface())
	continue;
      if (regions_[kj].get() == adj1 || regions_[kj].get() == adj2)
	continue;

      shared_ptr<ParamSurface> surf2 = regions_[kj]->getSurface(0)->surface();
      if (surf2->instanceType() != Class_Plane)
	continue;

      if (!(regions_[ix]->isAdjacent(regions_[kj].get())
	    || regions_[ix]->isNextToAdjacent(regions_[kj].get())))
	continue;

      // Check if the plane and the cylinder already has a common RevEngEdge
      vector<RevEngEdge*> revedg2 = regions_[kj]->getAllRevEdges();
      size_t kr, kh;
      for (kr=0; kr<revedg1.size(); ++kr)
	{
	  for (kh=0; kh<revedg2.size(); ++kh)
	    if (revedg1[kr] == revedg2[kh])
	      break;
	  if (kh < revedg2.size())
	    break;
	}
      if (kr < revedg1.size())
	continue;
      
#ifdef DEBUG_BLEND
      std::ofstream of1("torus_blend1.g2");
      regions_[ix]->writeRegionPoints(of1);
      regions_[ix]->writeSurface(of1);
      regions_[kj]->writeRegionPoints(of1);
      regions_[kj]->writeSurface(of1);
#endif
      
      shared_ptr<ElementarySurface> elem2 =
	dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf2);
      if (!elem2->isBounded())
	elem2->setParameterBounds(-0.5*diag, -0.5*diag, 0.5*diag, 0.5*diag);
      Point dir2 = elem2->direction();
      double ang = dir1.angle(dir2);
      ang = fabs(pihalf - ang);
      if (ang > blendfac*angtol)
	{
	  double usz, vsz;
	  elem1->estimateSfSize(usz, vsz);
	  double lenlim = 0.9*usz;

	  vector<shared_ptr<RevEngEdge> > edges =
	    defineEdgesBetween(ix, elem1, dir1, kj, elem2, dir2, false,
			       lenlim, false);
	  if (edges.size() > 0)
	    {
	      size_t ix2 = edges_.size();
	      edges_.insert(edges_.end(), edges.begin(), edges.end());
	      bool found = createTorusBlend(ix2);
	      // if (!found)
	      // 	found = createBlendSurface((int)ix2);
	      if (found)
		found_torus = true;
	      else
		edges_.erase(edges_.begin()+ix2, edges_.end());
	    }
	}
    }

  return found_torus;
}

//===========================================================================
void RevEng::defineMissingCorner(vector<RevEngRegion*>& cand_adj)
//===========================================================================
{
  double pihalf = 0.5*M_PI;
  double angtol = 5.0*anglim_;
  double diag = bbox_.low().dist(bbox_.high());
  
  // Get candidate opposite regions
  vector<RevEngRegion*> opposite_reg;
  for (size_t ki=0; ki<cand_adj.size(); ++ki)
    {
      RevEngEdge* edg1 = cand_adj[ki]->getBlendEdge();
      if (!edg1)
	continue;
      RevEngRegion *adj1, *adj2;
      edg1->getAdjacent(adj1, adj2);
      if ((!adj1) || (!adj2))
	continue;
      for (size_t kj=ki+1; kj<cand_adj.size(); ++kj)
	{
	  RevEngEdge* edg2 = cand_adj[kj]->getBlendEdge();
	  if (!edg2)
	    continue;
	  RevEngRegion *adj3, *adj4;
	  edg2->getAdjacent(adj3, adj4);
	  if ((!adj3) || (!adj4))
	    continue;
	  if (adj1 == adj3 || adj1 == adj4)
	    opposite_reg.push_back(adj1);
	  if (adj2 == adj3 || adj2 == adj4)
	    opposite_reg.push_back(adj2);
	}
    }
#ifdef DEBUG_BLEND
  if (opposite_reg.size() > 0)
    {
      std::ofstream of1("opposite_reg.g2");
      for (size_t ki=0; ki<opposite_reg.size(); ++ki)
	opposite_reg[ki]->writeRegionPoints(of1);
    }
#endif
  
  double blendfac = 2.0;
  for (size_t ki=0; ki<cand_adj.size(); ++ki)
    {
      RevEngEdge* edg1 = cand_adj[ki]->getBlendEdge();
      RevEngRegion *adj1, *adj2;
      edg1->getAdjacent(adj1, adj2);
      if ((!adj1) || (!adj2))
	continue;
      shared_ptr<ParamSurface> surf1 = cand_adj[ki]->getSurface(0)->surface();
      if (surf1->instanceType() != Class_Cylinder)
	continue;
      if (!surf1->isBounded())
	{
	  shared_ptr<Cylinder> cyl =
	    dynamic_pointer_cast<Cylinder,ParamSurface>(surf1);
	  cyl->setParamBoundsV(-0.5*diag, 0.5*diag);
	}
      shared_ptr<ElementarySurface> elem1 =
	dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf1);
      Point dir1 = elem1->direction();
      
      for (size_t kj=0; kj<opposite_reg.size(); ++kj)
	{
	  if (!opposite_reg[kj]->hasSurface())
	    continue;
	  if (opposite_reg[kj] == adj1 || opposite_reg[kj] == adj2)
	    continue;
	  shared_ptr<ParamSurface> surf2 =
	    opposite_reg[kj]->getSurface(0)->surface();
	  if (surf2->instanceType() != Class_Plane)
	    continue;
	  shared_ptr<ElementarySurface> elem2 =
	    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf2);
	  if (!elem2->isBounded())
	    elem2->setParameterBounds(-0.5*diag, -0.5*diag, 0.5*diag, 0.5*diag);
	  Point dir2 = elem2->direction();
	  double ang = dir1.angle(dir2);
	  ang = fabs(pihalf - ang);
	  if (ang > blendfac*angtol)
	    {
	      double usz, vsz;
	      elem1->estimateSfSize(usz, vsz);
	      double lenlim = 0.9*usz;

	      size_t kr, kh;
	      for (kr=0; kr<regions_.size(); ++kr)
		if (regions_[kr].get() == cand_adj[ki])
		  break;
	      if (kr == regions_.size())
		continue;
	      for (kh=0; kh<regions_.size(); ++kh)
		if (regions_[kh].get() == opposite_reg[kj])
		  break;
	      if (kh == regions_.size())
		continue;
#ifdef DEBUG_BLEND
	      std::ofstream of2("corner_adj.g2");
	      regions_[kr]->writeRegionPoints(of2);
	      regions_[kr]->writeSurface(of2);
	      regions_[kh]->writeRegionPoints(of2);
	      regions_[kh]->writeSurface(of2);
#endif
	      vector<shared_ptr<RevEngEdge> > edges =
		defineEdgesBetween(kr, elem1, dir1, kh, elem2, dir2, false,
				   lenlim, false);
	      if (edges.size() > 0)
		{
		  size_t ix2 = edges_.size();
		  edges_.insert(edges_.end(), edges.begin(), edges.end());
		  bool found = createTorusBlend(ix2);
		  if (!found)
		    edges_.erase(edges_.begin()+ix2, edges_.end());
		}
	    }
	  
	}
    }
   int stop_break = 1;
}


//===========================================================================
bool getAdjacentToTorus(RevEngEdge* edge, vector<RevEngEdge*>& rev_edgs,
			double tol, RevEngEdge*& adj_edg1,
			RevEngEdge*& adj_edg2, double& rad1, double& rad2)
//===========================================================================
{
  adj_edg1 = adj_edg2 = 0;
  rad1 = rad2 = -1.0;
  
  Point pos1, pos2;
  edge->getCrvEndPoints(pos1, pos2);

  double mindist1 = std::numeric_limits<double>::max();
  double mindist2 = std::numeric_limits<double>::max();
  int min_ix1 = -1, min_ix2 = -1;
  for (size_t ki=0; ki<rev_edgs.size(); ++ki)
    {
      RevEngRegion *reg = rev_edgs[ki]->getBlendRegSurf();
      if (!reg)
	continue;

      double tpar1, tpar2, dist1, dist2;
      Point close1, close2;
      rev_edgs[ki]->closestPoint(pos1, tpar1, close1, dist1);
      rev_edgs[ki]->closestPoint(pos2, tpar2, close2, dist2);
      if (dist1 < mindist1)
	{
	  mindist1 = dist1;
	  min_ix1 = (int)ki;
	}
       if (dist2 < mindist2)
	{
	  mindist2 = dist2;
	  min_ix2 = (int)ki;
	}
      int stop_break = 1;
    }
  if (min_ix1 < 0 && min_ix2 < 0)
    return false;
  if (mindist1 > tol && mindist2 > tol)
    return false;   // Might need to tune tolerance

  if (min_ix1 >= 0 && mindist1 <= tol)
    {
      adj_edg1 = rev_edgs[min_ix1];
      RevEngRegion *reg = adj_edg1->getBlendRegSurf();
      if (reg->hasSurface())
	{
	  shared_ptr<ParamSurface> bsurf = reg->getSurface(0)->surface();
	  shared_ptr<ElementarySurface> belem =
	    dynamic_pointer_cast<ElementarySurface,ParamSurface>(bsurf);
	  if (bsurf->instanceType() == Class_Cylinder /*||
							bsurf->instanceType() == Class_Cone*/)
	    rad1 = belem->radius(0.0, 0.0);
	  else if (bsurf->instanceType() == Class_Torus)
	    rad1 = belem->radius2(0.0, 0.0);
	}
    }
  if (min_ix2 >= 0 && mindist2 <= tol)
    {
      adj_edg2 = rev_edgs[min_ix2];
      RevEngRegion *reg = adj_edg2->getBlendRegSurf();
      if (reg->hasSurface())
	{
	  shared_ptr<ParamSurface> bsurf = reg->getSurface(0)->surface();
	  shared_ptr<ElementarySurface> belem =
	    dynamic_pointer_cast<ElementarySurface,ParamSurface>(bsurf);
	  if (bsurf->instanceType() == Class_Cylinder /*||
							bsurf->instanceType() == Class_Cone*/)
	    rad2 = belem->radius(0.0, 0.0);
	  else if (bsurf->instanceType() == Class_Torus)
	    rad2 = belem->radius2(0.0, 0.0);
	}
    }

  if (rad1 < 0 && rad2 < 0)
    return false;

  return true;
}

//===========================================================================
bool RevEng::createTorusBlend(size_t ix)
//===========================================================================
{
  double eps = 1.0e-6;
  double angtol = 5.0*anglim_;
  double pihalf = 0.5*M_PI;
  double tol2 = 2*approx_tol_;  // Due to possible inaccuracies
  double tol5 = 5*approx_tol_;  // Due to possible inaccuracies

  // Adjacent regions
  RevEngRegion* adj[2];
  edges_[ix]->getAdjacent(adj[0], adj[1]);

  if ((!adj[0]->hasSurface()) || (!adj[1]->hasSurface()))
    return false;  // Something wrong
  
  vector<shared_ptr<CurveOnSurface> > intcv1, intcv2;
  edges_[ix]->getCurve(intcv1, true);
  edges_[ix]->getCurve(intcv2, false);
  
  bool out1 = false, out2 = false;
  edges_[ix]->getOuterInfo(out1, out2);
  
  shared_ptr<ParamSurface> surf[2];
  surf[0] = adj[0]->getSurface(0)->surface();
  surf[1] = adj[1]->getSurface(0)->surface();
  ClassType type[2];
  type[0] = surf[0]->instanceType();
  type[1] = surf[1]->instanceType();
  int jx1 = (type[0] == Class_Plane) ? 0 :
    ((type[1] == Class_Plane) ? 1 : -1);
  if (jx1 == -1)
    return false;
  int jx2 = 1 - jx1;
  if (type[jx2] != Class_Cylinder)
    return false;
  bool outp = (jx1 == 0) ? out1 : out2;
  bool outr = (jx1 == 0) ? out2 : out1;
  shared_ptr<Cylinder> cyl = dynamic_pointer_cast<Cylinder,ParamSurface>(surf[jx2]);
  if (!cyl.get())
    return false;  // Should not happen
  double radius1 = cyl->getRadius();

  // Get adjacent blend edges and corresponding radii
  vector<RevEngEdge*> rev_edgs = adj[jx1]->getAllRevEdges();
  if (rev_edgs.size() == 0)
    return false;

  double rad1 = -1.0, rad2 = -1.0;
  RevEngEdge *revedg1, *revedg2;
  bool OK = getAdjacentToTorus(edges_[ix].get(), rev_edgs, tol5, 
			       revedg1, revedg2, rad1, rad2);
// #ifdef DEBUG_BLEND
//   if (!OK)
//     std::cout << "getAdjacentToTorus not OK" << std::endl;
//   if (!revedg1)
//     std::cout << "getAdjacentToTorus, revedg1 missing" << std::endl;
//   if (!revedg2)
//     std::cout << "getAdjacentToTorus, revedg2 missing" << std::endl;
// #endif
  if (!OK)
    return false;
  if ((!revedg1) || (!revedg2))
    return false;
  bool possible_suitcase =
    (fabs(rad1-rad2) > std::min(fabs(radius1-rad1), fabs(radius1-rad2)));

  double radius2 = (rad1 > 0.0 && rad2 > 0.0) ? 0.5*(rad1 + rad2) :
    ((rad1 > 0.0) ? rad1 : rad2);  // Should do extra checking if rad1 very
  // different from rad2 and both are larger than zero

  Point loc1 = cyl->getLocation();
  Point axis1 = cyl->getAxis();
  Point Cx = cyl->direction2();

  // The torus location coincides with the intersection point between the
  // cylinder axis and the plane
  shared_ptr<Plane> plane = dynamic_pointer_cast<Plane,ParamSurface>(surf[jx1]);
  Point dirp = plane->direction();
  Point avnorm = adj[jx1]->getMeanNormal();
  if (dirp*avnorm < 0.0)
    dirp *= -1;
  if (axis1*dirp < 0.0)
    axis1 *= -1;
  double alpha = axis1.angle(dirp);
  Point locp = plane->location();
  Point loct0 = loc1 + ((locp - loc1)*axis1)*axis1;
  double d2 = locp.dist(loct0);
  double xlen = d2*atan(alpha);
  int sgn1 = ((locp - loct0)*dirp < 0.0) ? -1 : 1;
  loct0 += sgn1*xlen*dirp;
  int sgn2 = outp ? 1 : -1;
  Point loct = loct0 + sgn2*radius2*dirp;
  int sgn3 = outr ? -1 : 1;
  double radiust = radius1 + sgn3*radius2;
  if (radius2 >= radiust)
    return false;  // Not expecting degenerate torus

  shared_ptr<Torus> torus(new Torus(radiust, radius2, loct, dirp, Cx));

  // Regions in blend area
  vector<RevEngRegion*> blend_regs;
  edges_[ix]->getAllBlendRegs(blend_regs);
#ifdef DEBUG_BLEND
  std::ofstream of("torus_blend2.g2");
  for (size_t ki=0; ki<blend_regs.size(); ++ki)
    blend_regs[ki]->writeRegionPoints(of);
#endif
  
  // Parameterize on torus
  vector<RevEngPoint*> blend_pts;
  for (size_t ki=0; ki<blend_regs.size(); ++ki)
    blend_pts.insert(blend_pts.end(), blend_regs[ki]->pointsBegin(),
		     blend_regs[ki]->pointsEnd());
  double maxd, avd;
  int num_in, num2_in;
  vector<double> parvals;
  vector<pair<double,double> > dist_ang;
  vector<RevEngPoint*> inpt, outpt;
  RevEngUtils::distToSurf(blend_pts.begin(), blend_pts.end(), torus,
			  approx_tol_, maxd, avd, num_in, num2_in, inpt, outpt,
			  parvals, dist_ang, angtol);
  
  RectDomain dom = cyl->getParameterBounds();
  if (sgn3 < 0)
    torus->setParameterBounds(dom.umin(), 0.0, dom.umax(), pihalf);
  else
    torus->setParameterBounds(dom.umin(), M_PI, dom.umax(), 1.5*M_PI);
  // What if the cylinder and the plane are not perpendicular? It is not a
  // torus, but it is anyway not handled
  
#ifdef DEBUG_BLEND
  torus->writeStandardHeader(of);
  torus->write(of);
  vector<shared_ptr<ParamCurve> > c_cvs = torus->constParamCurves(0.0, true);
  c_cvs[0]->writeStandardHeader(of);
  c_cvs[0]->write(of);
#endif

  // Define longitudinal boundary edges of torus
  RectDomain tordom = torus->getParameterBounds();
  double torlim[4];
  torlim[0] = tordom.umin();
  torlim[1] = tordom.umax();
  torlim[2] = tordom.vmin();
  torlim[3] = tordom.vmax();
  vector<shared_ptr<CurveOnSurface> > torbound(4);
  for (int ka=0; ka<2; ++ka)
    {
      shared_ptr<Circle> circle = torus->getMajorCircle(torlim[2+ka]);
      Point mid(torlim[0], torlim[2+ka]);
      Point vec = Point(torlim[1],torlim[2+ka]) - Point(torlim[0],torlim[2+ka]);
      vec.normalize();
      shared_ptr<Line> line(new Line(mid, vec));
      line->setParamBounds(0.0, torlim[1]-torlim[0]);
      torbound[2+ka] =
	shared_ptr<CurveOnSurface>(new CurveOnSurface(torus, line, circle, false,
						      3, 2, torlim[2+ka], 2+ka, true));
    }

  // Restrict adjacent cylinders
  // The cylinder and the edges does not necessarily have the same
  // parameterization
  Point pt1, pt2;
  double midpar = 0.5*(tordom.vmin() + tordom.vmax());
  torus->point(pt1, tordom.umin(), midpar);
  torus->point(pt2, tordom.umax(), midpar);
  vector<RevEngRegion*> adj_reg(2, 0);
  bool upper2[2];
  double lim[2];
  vector<shared_ptr<ElementarySurface> > adj_surf(2);
  vector<double> parbound(8);
  int kx[2];
  for (int ka=0; ka<2; ++ka)
    {
      if ((ka == 0 && rad1 <= 0.0) || (ka == 1 && rad2 <= 0.0))
	{
	  kx[ka] = -1;
	  continue;
	}

      RevEngRegion *reg = (ka == 0) ? revedg1->getBlendRegSurf() :
	revedg2->getBlendRegSurf();
      shared_ptr<ParamSurface> bsurf = reg->getSurface(0)->surface();
      shared_ptr<ElementarySurface> belem =
	dynamic_pointer_cast<ElementarySurface,ParamSurface>(bsurf);
      adj_reg[ka] = reg;
      adj_surf[ka] = belem;
      
     double reg_dom[4];
      reg->getDomain(reg_dom);
      Point pt3[2];
      int kx1 = (bsurf->instanceType() == Class_Cylinder) ? 0 : 1;
      for (int kb=0; kb<2; ++kb)
	{
	  double par[2];
	  par[kx1] = 0.5*(reg_dom[2*kx1] + reg_dom[2*kx1+1]);
	  par[1-kx1] = reg_dom[2*(1-kx1)+kb];
	  pt3[kb] = bsurf->point(par[0], par[1]);
	}
      
      double upar1, vpar1, bdist1, upar2, vpar2, bdist2;
      Point bclose1, bclose2;
      bsurf->closestPoint(pt1, upar1, vpar1, bclose1, bdist1, eps);
      bsurf->closestPoint(pt2, upar2, vpar2, bclose2, bdist2, eps);
      double uvpar;
      if (belem->instanceType() == Class_Cylinder)
	uvpar = (bdist1 <= bdist2) ? vpar1 : vpar2;
      else
	uvpar = (bdist1 <= bdist2) ? upar1 : upar2;
      Point tor_pt = (bdist1 <= bdist2) ? pt1 : pt2;
      lim[ka] = uvpar;
      kx[ka] = (bdist1 <= bdist2) ? 0 : 1;

      // Assume the part of the cylinder to remove is the smaller one
      RectDomain bdom = bsurf->containingDomain();
      parbound[4*ka] = bdom.umin();
      parbound[4*ka+1] = bdom.vmin();
      parbound[4*ka+2] = bdom.umax();
      parbound[4*ka+3] = bdom.vmax();
      if (tor_pt.dist(pt3[1]) < tor_pt.dist(pt3[0]))
	{
	  if (belem->instanceType() == Class_Cylinder)
	    parbound[4*ka+3] = uvpar;
	  else
	    parbound[4*ka+2] = uvpar;
	  upper2[ka] = true;
	}
      else 
	{
	  if (belem->instanceType() == Class_Cylinder)
	    parbound[4*ka+1] = uvpar;
	  else
	    parbound[4*ka] = uvpar;
	  upper2[ka] = false;
	}
#ifdef DEBUG_BLEND
      bsurf->writeStandardHeader(of);
      bsurf->write(of);
#endif
    }

  Point vec1 = (adj_surf[0].get()) ? adj_surf[0]->location() - plane->location() : Point(0.0, 0.0, 0.0);
  Point vec2 = (adj_surf[1].get()) ? adj_surf[1]->location() - plane->location() : Point(0.0, 0.0, 0.0);
  double sc1 = vec1*dirp;
  double sc2 = vec2*dirp;
  if (sc1*sc2 < 0.0) //kx[0] == kx[1] && kx[0]>=0)
    {
      // Need to wait with "setParameterBounds"
      vector<RevEngRegion*> blends(3);
      blends[0] = adj[jx2];
      blends[1] = revedg1->getBlendRegSurf();
      blends[2] = revedg2->getBlendRegSurf();
      return suitcaseCorner(blends, edges_[ix].get());
    }
  else if (possible_suitcase)
    return false;  // Not an expected configuration

  if (kx[0] == kx[1])
    return false;   // Not an expected configuration

  for (int ka=0; ka<2; ++ka)
    {
      if (adj_surf[ka].get())
	{
	  adj_surf[ka]->setParameterBounds(parbound[4*ka], parbound[4*ka+1],
					   parbound[4*ka+2], parbound[4*ka+3]);
	  adj_reg[ka]->adaptEdges();
	}
   }

  // Bound blend cylinder
  Point pt3;
  torus->point(pt3, 0.5*(dom.umin()+dom.vmin()), 0);
  double upar3, vpar3, dist3;
  Point close3;
  bool upper = false;
  cyl->closestPoint(pt3, upar3, vpar3, close3, dist3, eps);
  if (fabs(vpar3-dom.vmin()) < fabs(dom.vmax()-vpar3))
    cyl->setParamBoundsV(vpar3, dom.vmax());
  else
    {
      cyl->setParamBoundsV(dom.vmin(), vpar3);
      upper = true;
    }
  adj[jx2]->adaptEdges();
  
 #ifdef DEBUG_BLEND
  cyl->writeStandardHeader(of);
  cyl->write(of);
#endif

  shared_ptr<Circle> bcircle = cyl->getCircle(vpar3);
  shared_ptr<ElementaryCurve> bpar =
    cyl->getElementaryParamCurve(bcircle.get(), approx_tol_);
  shared_ptr<CurveOnSurface> blendbound(new CurveOnSurface(cyl, bpar, bcircle,
							   false, 3, 2, vpar3,
							   (upper) ? 3 : 2, true));
 
  // Define traverse boundary edges of torus and adjacent cylinders
  vector<shared_ptr<CurveOnSurface> > cylbound(2);
  for (int ka=0; ka<2; ++ka)
    {
      if (!adj_reg[ka])
	continue;

      bool udir = (adj_surf[ka]->instanceType() == Class_Cylinder) ? true : false;
      shared_ptr<Circle> circle = torus->getMinorCircle(torlim[ka]);
      Point mid(torlim[ka], torlim[2]);
      Point vec = Point(torlim[ka],torlim[3]) - Point(torlim[ka],torlim[2]);
      vec.normalize();
      shared_ptr<Line> line(new Line(mid, vec));
      line->setParamBounds(0.0, torlim[3]-torlim[2]);
      torbound[ka] =
	shared_ptr<CurveOnSurface>(new CurveOnSurface(torus, line, circle,
						      false, 3, 1, torlim[ka],
						      ka, true));
      vector<shared_ptr<ParamCurve> > tmpcv =
	adj_surf[ka]->constParamCurves(lim[ka], udir);
      if (tmpcv.size() != 1)
	continue;  // Should not happen
      shared_ptr<ElementaryCurve> space =
	dynamic_pointer_cast<ElementaryCurve,ParamCurve>(tmpcv[0]);
      shared_ptr<ElementaryCurve> par =
	adj_surf[ka]->getElementaryParamCurve(space.get(), approx_tol_);
      int bd;
      if (udir)
	bd = (upper2[ka]) ? 3 : 2;
      else
	bd = (upper2[ka]) ? 1 : 0;
      shared_ptr<ParamCurve> par2;
      if (!par.get())
	{
	  Point pt1(2), pt2(2);
	  if (udir)
	    {
	      pt1[1] = pt2[1] = lim[ka];
	      pt1[0] = parbound[4*ka];
	      pt2[0] = parbound[4*ka+2];
	    }
	  else
	    {
	      pt1[0] = pt2[0] = lim[ka];
	      pt1[1] = parbound[4*ka+1];
	      pt2[1] = parbound[4*ka+3];
	    }
	  par2 = shared_ptr<ParamCurve>(new SplineCurve(pt1, space->startparam(),
							pt2, space->endparam()));
	}
      cylbound[ka] =
	shared_ptr<CurveOnSurface>(new CurveOnSurface(adj_surf[ka],
						      par.get() ? par : par2,
						      space,
						      false, 3, udir ? 1 : 2,
						      lim[ka],
						      bd, true));
    }

  // Define planar boundary curve
  shared_ptr<CurveOnSurface> planebound;
  shared_ptr<ParamCurve> pspace1 = torbound[3]->spaceCurve();
  shared_ptr<Circle> pspace = dynamic_pointer_cast<Circle,ParamCurve>(pspace1);
  shared_ptr<ElementaryCurve> ppar =
    plane->getElementaryParamCurve(pspace.get(), 10.0*tol2);  // Just to test
  planebound = shared_ptr<CurveOnSurface>(new CurveOnSurface(plane, ppar,
							     pspace, false, 1));

#ifdef DEBUG_BLEND
  std::ofstream ofc("torus_trim.g2");
  for (int ka=0; ka<4; ++ka)
    {
      if (!torbound[ka].get())
	continue;
      shared_ptr<ParamCurve> tmp = torbound[ka]->spaceCurve();
      tmp->writeStandardHeader(ofc);
      tmp->write(ofc);
    }
  for (int ka=0; ka<2; ++ka)
    {
      if (!cylbound[ka].get())
	continue;
      shared_ptr<ParamCurve> tmp = cylbound[ka]->spaceCurve();
      tmp->writeStandardHeader(ofc);
      tmp->write(ofc);
    }
  if (planebound.get())
    {
      shared_ptr<ParamCurve> tmp = planebound->spaceCurve();
      tmp->writeStandardHeader(ofc);
      tmp->write(ofc);
    }

  std::ofstream ofp("torus_par.g2");
  for (int ka=0; ka<4; ++ka)
    {
      if (!torbound[ka].get())
	continue;
      shared_ptr<ParamCurve> tmp = torbound[ka]->parameterCurve();
      if (!tmp.get())
	continue;
      SplineDebugUtils::writeSpaceParamCurve(tmp, ofp, 0.0);
    }
#endif
  
  // Move points as appropriate
  // From torus
  vector<vector<RevEngPoint*> > move1(4);
  vector<RevEngPoint*> torus_pts;
  for (size_t ki=0; ki<blend_pts.size(); ++ki)
    {
      if (parvals[2*ki+1] > M_PI)
	move1[0].push_back(blend_pts[ki]);  // To blend cylinder
      else if (parvals[2*ki+1] > pihalf)
	move1[1].push_back(blend_pts[ki]);  // To plane
      else if (parvals[2*ki] < dom.umin())
	move1[2].push_back(blend_pts[ki]);  // First adjacent cylinder
      else if (parvals[2*ki] > dom.umax())
	move1[3].push_back(blend_pts[ki]);  // Second adjacent cylinder
      else
	{
	  blend_pts[ki]->setPar(Vector2D(parvals[2*ki],parvals[2*ki+1]));
	  blend_pts[ki]->setSurfaceDist(dist_ang[ki].first, dist_ang[ki].second);
	  torus_pts.push_back(blend_pts[ki]);
	}
    }

  // From blend cylinder
  vector<RevEngPoint*> move_cyl;
  int num_pts_cyl = adj[jx2]->numPoints();
  for (int ka=0; ka<num_pts_cyl; ++ka)
    {
      RevEngPoint* curr = adj[jx2]->getPoint(ka);
      Vector2D par = curr->getPar();
      if ((upper && par[1] > vpar3) || (upper == false && par[1] < vpar3))
	move_cyl.push_back(curr);
    }

  // Remove identified points from cylinder
  if (move_cyl.size() > 0)
    adj[jx2]->removeAndUpdatePoints(move_cyl);

  // From adjacent cylinders
  vector<vector<RevEngPoint*> > move_adj(2);
  for (int kb=0; kb<2; ++kb)
    {
      if (!adj_reg[kb])
	continue;
      int num_pts_cyl = adj_reg[kb]->numPoints();
      for (int ka=0; ka<num_pts_cyl; ++ka)
	{
	  RevEngPoint* curr = adj_reg[kb]->getPoint(ka);
	  Vector2D par = curr->getPar();
	  if ((upper2[kb] && par[1] > lim[kb]) ||
	      (upper2[kb] == false && par[1] < lim[kb]))
	    move_adj[kb].push_back(curr);
	}
      if (move_adj[kb].size() > 0)
	adj_reg[kb]->removeAndUpdatePoints(move_adj[kb]);
    }

  // From adjacent plane
  vector<RevEngPoint*> move_plane;
  adj[jx1]->extractOutOfEdge(planebound, (jx1 == 0) ? intcv1 : intcv2,
			     radius2, approx_tol_, angtol, move_plane);

  // To plane
  bool OK1, OK2, OK3;
  if (move1[1].size() > 0)
    OK1 = adj[jx1]->addPointsToGroup(move1[1], approx_tol_, angtol);

  // To blend cylinder
  if (move1[0].size() > 0)
    OK2 = adj[jx2]->addPointsToGroup(move1[0], approx_tol_, angtol);

  // To adjacent cylinders
  vector<vector<RevEngPoint*> > added_groups;
  for (int ka=0; ka<2; ++ka)
    if (move1[2+ka].size() > 0)
      {
	if (adj_reg[ka])
	  OK3 = adj_reg[ka]->addPointsToGroup(move1[2+ka], approx_tol_, angtol);
	else
	  added_groups.push_back(move1[2+ka]);
      }


  // To torus
  if (move_plane.size() > 0)
    torus_pts.insert(torus_pts.end(), move_plane.begin(), move_plane.end());
  if (move_adj[0].size() > 0)
    torus_pts.insert(torus_pts.end(), move_adj[0].begin(), move_adj[0].end());
  if (move_adj[1].size() > 0)
    torus_pts.insert(torus_pts.end(), move_adj[1].begin(), move_adj[1].end());
  if (move_cyl.size() > 0)
    torus_pts.insert(torus_pts.end(), move_cyl.begin(), move_cyl.end());
  if (torus_pts.size() == 0)
    return false;

  // Define torus blend
  shared_ptr<RevEngRegion> blendreg(new RevEngRegion(classification_type_,
						     edge_class_type_,
						     torus_pts));
  blendreg->setRegionAdjacency();
  regions_.push_back(blendreg);
  shared_ptr<HedgeSurface> hedge;
  shared_ptr<ParamSurface> torus_tmp = torus;
  blendreg->setAssociatedSurface(torus_tmp, approx_tol_, angtol,
				 min_point_region_, hedge);
  if (hedge.get())
    surfaces_.push_back(hedge);
  
  if (added_groups.size() > 0)
    {
      vector<HedgeSurface*> dummy_hedge;
      surfaceExtractOutput((int)regions_.size()-1, added_groups, dummy_hedge);
    }

#ifdef DEBUG_BLEND
  std::ofstream oft("torus_corner.g2");
  blendreg->writeRegionPoints(oft);
  adj[jx1]->writeRegionPoints(oft);
  adj[jx2]->writeRegionPoints(oft);
  if (adj_reg[0])
    adj_reg[0]->writeRegionPoints(oft);
  if (adj_reg[1])
    adj_reg[1]->writeRegionPoints(oft);
#endif
  
  // Add edges to regions
  int status = 0;
  vector<shared_ptr<ftEdge> > tor_edg(4);
  for (int ka=0; ka<4; ++ka)
    {
      if (torbound[ka].get())
	{
	  tor_edg[ka] = shared_ptr<ftEdge>(new ftEdge(hedge.get(), torbound[ka],
						      torbound[ka]->startparam(),
						      torbound[ka]->endparam()));
	  blendreg->addTrimEdge(tor_edg[ka]);
	}
    }

  shared_ptr<ftEdge> plane_edg(new ftEdge(adj[jx1]->getSurface(0), planebound,
					  planebound->startparam(),
					  planebound->endparam()));
  adj[jx1]->addTrimEdge(plane_edg);
  plane_edg->setReversed(true);
  tor_edg[3]->connectTwin(plane_edg.get(), status);


   shared_ptr<ftEdge> cyl_edg(new ftEdge(adj[jx2]->getSurface(0), blendbound,
					 blendbound->startparam(),
					 blendbound->endparam()));
  adj[jx2]->addTrimEdge(cyl_edg);
  cyl_edg->setReversed(true);
  tor_edg[2]->connectTwin(cyl_edg.get(), status);

  for (int ka=0; ka<2; ++ka)
    {
      if (cylbound[ka].get())
	{
	  shared_ptr<ftEdge> adj_edg(new ftEdge(adj_reg[ka]->getSurface(0),
						cylbound[ka],
						cylbound[ka]->startparam(),
						cylbound[ka]->endparam()));
	  adj_reg[ka]->addTrimEdge(adj_edg);
	  adj_edg->setReversed(true);
	  if (kx[ka] >= 0 && tor_edg[kx[ka]].get())
	    tor_edg[kx[ka]]->connectTwin(adj_edg.get(), status);
	}
    }

  for (size_t ki=0; ki<blend_regs.size(); ++ki)
    {
      //blend_regs[ki]->removeFromAdjacent();
      blend_regs[ki]->setRemove();
    }
   edges_[ix]->clearBlendRegions();
 
  return true;
}

//===========================================================================
void RevEng::setBlendBoundaries(RevEngRegion *reg)
//===========================================================================
{
  double eps = 1.0e-6;
  if (!reg->hasSurface())
    return;

  shared_ptr<ParamSurface> surf = reg->getSurface(0)->surface();
  if (surf->instanceType() != Class_Cylinder &&
      surf->instanceType() != Class_Torus)
    return;   

  double angtol = 5.0*anglim_;
  double diag0 = bbox_.low().dist(bbox_.high());
  
  // Adjacent regions
  RevEngEdge *edge = reg->getBlendEdge();
  RevEngRegion* adj[2];
  edge->getAdjacent(adj[0], adj[1]);

  if ((!adj[0]->hasSurface()) || (!adj[1]->hasSurface()))
    return;  // Something wrong

  bool out1 = false, out2 = false;
  edge->getOuterInfo(out1, out2);

  vector<shared_ptr<CurveOnSurface> > intcv1, intcv2;
  edge->getCurve(intcv1, true);
  edge->getCurve(intcv2, false);
  //double eps1 = approx_tol_; //std::max(1.0e-4, 0.1*approx_tol_);
  // for (size_t ki=0; ki<intcv1.size(); ++ki)
  //   intcv1[ki]->fixMismatchCurves(eps1);
  // for (size_t ki=0; ki<intcv2.size(); ++ki)
  //   intcv2[ki]->fixMismatchCurves(eps1);
  
#ifdef DEBUG_BLEND
  std::ofstream of("blend.g2");
  reg->writeRegionPoints(of);
  surf->writeStandardHeader(of);
  surf->write(of);

  adj[0]->writeRegionPoints(of);
  adj[1]->writeRegionPoints(of);
#endif

  // Extract longitudial boundary curves
  shared_ptr<ElementarySurface> elem =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf);
  //double radius = elem->radius(0.0, 0.0);

  vector<shared_ptr<ElementarySurface> > adj_elem(2);
  for (int ka=0; ka<2; ++ka)
    {
      shared_ptr<ParamSurface> tmp = adj[ka]->getSurface(0)->surface();
      adj_elem[ka] = dynamic_pointer_cast<ElementarySurface,ParamSurface>(tmp);
    }

#ifdef DEBUG_BLEND
  std::ofstream ofsf("adj_sfs.g2");
  adj_elem[0]->writeStandardHeader(ofsf);
  adj_elem[0]->write(ofsf);
  shared_ptr<ParamCurve> tmp_cv1 = intcv1[0]->spaceCurve();
  tmp_cv1->writeStandardHeader(ofsf);
  tmp_cv1->write(ofsf);
  adj_elem[1]->writeStandardHeader(ofsf);
  adj_elem[1]->write(ofsf);
  shared_ptr<ParamCurve> tmp_cv2 = intcv2[0]->spaceCurve();
  tmp_cv2->writeStandardHeader(ofsf);
  tmp_cv2->write(ofsf);
#endif

  Point posi1, posi2, pos;
  intcv1[0]->point(posi1, intcv1[0]->startparam());
  double pari = (edge->isClosed(approx_tol_)) ?
    0.5*(intcv1[0]->startparam() + intcv1[intcv1.size()-1]->endparam()) :
    intcv1[intcv1.size()-1]->endparam();
  intcv1[intcv1.size()-1]->point(posi2, pari);
  pos = 0.5*(posi1 + posi2);
  RectDomain surfdom = surf->containingDomain();
  double regdom[4];
  reg->getDomain(regdom);
  double regfac = 2.0;
  double rad;
  double tpar1, tpar2;
  bool udir;
  int constdir;
  double delfac = 0.6;
  double seamfac = 0.1;
  bool plane1 =  (adj_elem[0]->instanceType() == Class_Plane);
  bool plane2 =  (adj_elem[1]->instanceType() == Class_Plane);
  if (surf->instanceType() == Class_Cylinder)
    {
      Point norm1 = adj_elem[0]->direction();
      Point norm2 = adj_elem[1]->direction();
      double ang = norm1.angle(norm2);
      ang = std::min(ang, M_PI-ang);
      
      tpar1 = M_PI - 0.5*ang; 
      tpar2 = tpar1 + ang;
      shared_ptr<Cylinder> cyl = dynamic_pointer_cast<Cylinder,ParamSurface>(surf);
      cyl->setParamBoundsU(tpar1, tpar2);
      double tdelreg = regdom[3] - regdom[2];
      if (tdelreg < regfac*approx_tol_)
	tdelreg = 0.5*diag0;
      if (surfdom.vmax() - surfdom.vmin() > regfac*tdelreg)
	cyl->setParamBoundsV(std::max(surfdom.vmin(), regdom[2]-delfac*tdelreg),
			     std::min(surfdom.vmax(), regdom[3]+delfac*tdelreg));
      rad = cyl->getRadius();
      udir = false;
      constdir = 0;
    }
  else
    {
      shared_ptr<Torus> tor = dynamic_pointer_cast<Torus,ParamSurface>(surf);
      RectDomain tor_dom = tor->getParameterBounds();
      if (plane1 || plane2)
	{
	  int plane_ix = (plane1) ? 0 : 1;
	  if (adj_elem[plane_ix]->instanceType() != Class_Plane)
	    return;  // Not as expected
	  rad = tor->getMinorRadius();
	  bool plane_out = (plane_ix == 0) ? out1 : out2;
	  bool rot_out = (plane_ix == 0) ? out2 : out1;
	  double phi = 0.5*M_PI;
	  Point norm = adj_elem[plane_ix]->direction();
	  Point norm2 = adj[plane_ix]->getMeanNormalTriang();
	  double beta = 0.5; //(norm*norm2 > 0.0) ? 0.5 : 0.0;
	  if (adj_elem[1-plane_ix]->instanceType() == Class_Cone)
	    {
	      shared_ptr<Cone> cone =
		dynamic_pointer_cast<Cone,ElementarySurface>(adj_elem[1-plane_ix]);
	      double alpha = cone->getConeAngle();
	      Point axis = adj_elem[1-plane_ix]->direction();
	      int sgn = (norm*axis < 0.0) ? -1 : 1;
	      phi += sgn*alpha;
	    }
	  if (rot_out)
	    {
	      tpar2 = tor_dom.vmin()+(1.0+beta)*M_PI;
	      tpar1 = tpar2 - phi;
	    }
	  else
	    {
	      tpar1 = tor_dom.vmin()+(1.0-beta)*M_PI;
	      tpar2 = tpar1 + phi;
	    }
	  if (plane_out)
	    {
	      int sgn = 1; //(norm*norm2 < 0.0 && rot_out == false) ? -1 : 1;
	      tpar1 += 0.5*sgn*M_PI;
	      tpar2 += 0.5*sgn*M_PI;
	    }
	}
      else
	{
	  int cyl_ix = (adj_elem[0]->instanceType() == Class_Cylinder) ? 0 : 1;
	  int cone_ix = 1 - cyl_ix;
	  if (adj_elem[cyl_ix]->instanceType() != Class_Cylinder ||
	      adj_elem[cone_ix]->instanceType() != Class_Cone)
	    return;

	  shared_ptr<Cone> cone =
	    dynamic_pointer_cast<Cone,ElementarySurface>(adj_elem[cone_ix]);
	  double alpha = cone->getConeAngle();
	  double beta = M_PI - fabs(alpha);
	  double phi = 0.5*(M_PI - beta);
	  rad = tor->getMinorRadius();
	  Point axis = adj_elem[cyl_ix]->direction();
	  Point loc = adj_elem[cyl_ix]->location();
	  double rd = (pos-loc)*axis;
	  Point centr = loc + rd*axis;
	  Point loc2 = tor->location();
	  Point loc2_2 = centr + ((loc2 - centr)*axis)*axis;
	  double hh = loc2_2.dist(centr);
	  double d2 = hh/sin(phi) - rad;
	  double xlen = sqrt((rad+d2)*(rad+d2) - rad*rad);
	  tpar1 = tor_dom.vmin() + M_PI - 2.0*xlen;
	  tpar2 = tpar1 + 2*xlen;
	}
      
      double tdelreg = regdom[1] - regdom[0];
      double upar1 = tor_dom.umin();
      double upar2 = tor_dom.umax();
      double umid = 0.5*(regdom[0]+regdom[1]);
      double eps2 = std::max(eps, 0.001*(upar2-upar1));

      double cp1, cp2;
      Point cpos1, cpos2;
      double seed[2];
      seed[0] = 0.5*(regdom[0]+regdom[1]);
      seed[1] = 0.5*(regdom[2]+regdom[3]);
      int seam1 = edge->closedSfAtEnd(approx_tol_, cp1, cpos1, true);
      int seam2 = edge->closedSfAtEnd(approx_tol_, cp2, cpos2, false);
      double dd = cpos1.dist(cpos2);
      double u1, v1, u2, v2, d1, d2;
      Point cl1, cl2;
      elem->closestPoint(cpos1, u1, v1, cl1, d1, eps, 0, seed);
      elem->closestPoint(cpos2, u2, v2, cl2, d2, eps, 0, seed);
      if (u2 < u1)
	std::swap(u1, u2);
      if ((u2-u1 < tdelreg-eps || umid > u2 || umid < u1) && (seam1 || seam2))
	{
	  // Check parameter at seam
	  if (seam1 && u1-upar1 < eps2)
	    u1 = upar2;
	  else if (seam2 && upar2-u2 < eps2)
	    u2 = upar1;
	  if (u2 < u1)
	    std::swap(u1, u2);
	}
	  
      if (seam1 != seam2 && dd > approx_tol_)
	{
	  bool close1, close2;
	  tor->isClosed(close1, close2);  // Expects close1 = true
	  if (close1)
	    {
	      if (tdelreg > u2-u1)
		{
		  int stop_seam = 1;
		}
	      
	      // Check for seam. Solution could be improved
	      if (fabs(u1-upar1) < eps2 && umid > u2) // && upar2 < u2)
		u1 = upar2;
	      else if (fabs(upar2-u1) < eps2 && umid < u2)
		u1 = upar1;
	      if (fabs(u2-upar1) < eps2 && umid > u1)
		u2 = upar2;
	      else if (fabs(upar2-u2) < eps2 && umid < u1) // && upar1 > u1)
		u2 = upar1;

	      if (u2 < u1)
		std::swap(u1, u2);
	    }
	      
	  upar1 = std::max(upar1, u1);
	  upar2 = std::min(upar2, u2);
	}
      
      if (upar2 - upar1 > regfac*tdelreg)
	{
	  upar1 = std::max(upar1, std::min(regdom[0]-delfac*tdelreg, u1));
	  upar2 = std::min(upar2, std::max(regdom[1]+delfac*tdelreg, u2));
	  if (upar1 < seamfac)
	    upar1 = 0.0;
	  if (2*M_PI-upar2 < seamfac)
	    upar2 = 2*M_PI;
	}
      tor->setParameterBounds(upar1, tpar1, upar2, tpar2);
      udir = true;
      constdir = 1;
    }


  vector<shared_ptr<ParamCurve> > bdcv1 = elem->constParamCurves(tpar1, udir); 
  vector<shared_ptr<ParamCurve> > bdcv2 = elem->constParamCurves(tpar2, udir);
  if (bdcv1.size() != 1 || bdcv2.size() != 1)
    return;  // Something strange

  shared_ptr<ElementaryCurve> space[2];
  space[0] = dynamic_pointer_cast<ElementaryCurve,ParamCurve>(bdcv1[0]);
  space[1] = dynamic_pointer_cast<ElementaryCurve,ParamCurve>(bdcv2[0]);
  if ((!space[0].get()) || (!space[1].get()))
    return;
  BoundingBox bb1 = adj[0]->getBbox();
  BoundingBox bb2 = adj[1]->getBbox();
  double diag = std::min(bb1.low().dist(bb1.high()), bb2.low().dist(bb2.high()));
  if (constdir == 0 && (space[0]->startparam() < -0.5*diag ||
			space[0]->endparam() > 0.5*diag))
    {
      for (int ka=0; ka<2; ++ka)
	space[ka]->setParamBounds(std::max(space[0]->startparam(), -0.5*diag),
				  std::min(space[0]->endparam(), 0.5*diag));
    }
  
#ifdef DEBUG_BLEND
  space[0]->writeStandardHeader(of);
  space[0]->write(of);
  space[1]->writeStandardHeader(of);
  space[1]->write(of);
#endif
  RectDomain dom = elem->getParameterBounds();
  Point parpt1(dom.umin(), dom.vmin()), parpt2(dom.umin(), dom.vmax());
  Point parpt3(dom.umax(), dom.vmin()), parpt4(dom.umax(), dom.vmax());
  Point lpos1 = (udir) ? Point(0.0, dom.vmin()) : Point(dom.umin(), 0.0);
  Point lpos2 = (udir) ? Point(0.0, dom.vmax()) : Point(dom.umax(), 0.0);
  Point pvec = (udir) ? Point(1.0, 0.0) : Point(0.0, 1.0);

  double t1 = space[0]->startparam();
  double t2 = space[0]->endparam();
  shared_ptr<ElementaryCurve> par1(new Line(lpos1, pvec));
  par1->setParamBounds(t1, t2);
  shared_ptr<ElementaryCurve> par2(new Line(lpos2, pvec));
  par2->setParamBounds(t1, t2);
#ifdef DEBUG_BLEND
  // Check
  double tm1 = 0.75*t1 + 0.25*t2;
  double tm2 = 0.25*t1 + 0.75*t2;
  Point pp1, pp2, pp3, pp4;
  Point sp1, sp2, sp3, sp4;
  Point ssp1, ssp2, ssp3, ssp4;
  space[0]->point(sp1, tm1);
  space[0]->point(sp2, tm2);
  space[1]->point(sp3, tm1);
  space[1]->point(sp4, tm2);
  par1->point(pp1, tm1);
  par1->point(pp2, tm2);
  par2->point(pp3, tm1);
  par2->point(pp4, tm2);
  surf->point(ssp1, pp1[0], pp1[1]);
  surf->point(ssp2, pp2[0], pp2[1]);
  surf->point(ssp3, pp3[0], pp3[1]);
  surf->point(ssp4, pp4[0], pp4[1]);
#endif

  shared_ptr<ElementaryCurve> adj_par[2];
  Point close1, close2;
  double upar1, upar2, vpar1, vpar2, dist1, dist2;
  Point pos1, pos2;
  int ix[2];
  ix[0] = 0;
  ix[1] = 1;
  
  // Check for curve matching
  space[0]->point(pos1, t1);
  adj_elem[0]->closestPoint(pos1, upar1, vpar1, close1, dist1, eps);
  adj_elem[1]->closestPoint(pos1, upar2, vpar2, close2, dist2, eps);
  if (dist2 < dist1)
    std::swap(ix[0], ix[1]);

  double tol5 = 5.0*approx_tol_;
  for (int ka=0; ka<2; ++ka)
    {
      adj_par[ix[ka]] = adj_elem[ix[ka]]->getElementaryParamCurve(space[ka].get(),
								  tol5);
      if (!adj_par[ix[ka]].get())
	{
	  std::cout << "No parameter curve found" << std::endl;
	}

      // space[ka]->point(pos1, t1);
      // space[ka]->point(pos2, t2);
      // adj_elem[ix[ka]]->closestPoint(pos1, upar1, vpar1, close1, dist1, eps);
      // adj_elem[ix[ka]]->closestPoint(pos2, upar2, vpar2, close2, dist2, eps);
      // Point pp1(upar1,vpar1), pp2(upar2,vpar2);
      // Point mpar = (t2*pp1 -t1*pp2)/(t2-t1);
      // Point mvec = pp2 - pp1;
      // mvec.normalize();

      // adj_par[ix[ka]] = shared_ptr<ElementaryCurve>(new Line(mpar, mvec));
      // adj_par[ix[ka]]->setParamBounds(t1, t2);
    }
  
  // Create edges. Orientation in adjacent regions is not set
  int status = 0;
  vector<shared_ptr<ftEdge> > bdedg(2);
  shared_ptr<CurveOnSurface> sfcv1(new CurveOnSurface(elem, par1, space[0], false, 3,
						      constdir+1, tpar1, 2*constdir, true));
  bdedg[0] = shared_ptr<ftEdge>(new ftEdge(reg->getSurface(0), sfcv1, t1, t2));
#ifdef DEBUG_BLEND
  bool same_orient = sfcv1->sameOrientation();
  bool same_trace = sfcv1->sameTrace(approx_tol_);
  bool same_cv = sfcv1->sameCurve(approx_tol_);
  if ((!same_orient) || (!same_trace) || (!same_cv))
    std::cout << "Surface curve 1 mismatch " << same_orient << " " << same_trace << " " << same_cv << std::endl;
#endif
  shared_ptr<CurveOnSurface> sfcv2(new CurveOnSurface(elem, par2, space[1], false, 3,
						      constdir+1, tpar2, 2*constdir+1, true));
  //sfcv2->reverseParameterDirection();
  bdedg[1] = shared_ptr<ftEdge>(new ftEdge(reg->getSurface(0), sfcv2, t1, t2));
#ifdef DEBUG_BLEND
  same_orient = sfcv2->sameOrientation();
  same_trace = sfcv2->sameTrace(approx_tol_);
  same_cv = sfcv2->sameCurve(approx_tol_);
  if ((!same_orient) || (!same_trace) || (!same_cv))
    std::cout << "Surface curve 2 mismatch " << same_orient << " " << same_trace << " " << same_cv << std::endl;
#endif
  reg->addTrimEdge(bdedg[0]);
  reg->addTrimEdge(bdedg[1]);

  shared_ptr<CurveOnSurface> adj_sfcv[2];
  shared_ptr<ftEdge> adj_edg[2];
  vector<RevEngPoint*> out[2];
  for (int ka=0; ka<2; ++ka)
    {
      adj_sfcv[ix[ka]] = shared_ptr<CurveOnSurface>(new CurveOnSurface(adj_elem[ix[ka]],
								       adj_par[ix[ka]],
								       space[ka],
								       false, 1));
      // if (!adj_sfcv[ix[ka]]->hasParameterCurve())
      // 	{
      // 	  vector<shared_ptr<CurveOnSurface> > tmp_cvs;
      // 	  tmp_cvs.push_back(adj_sfcv[ix[ka]]);
      // 	  vector<pair<double,double> > t1_t2;
      // 	  adj[ix[ka]]->getCurveRestriction(tmp_cvs, approx_tol_, anglim_, t1_t2);
      // 	  if (t1_t2.size() == 1)
      // 	    {
      // 	      if (t1_t2[0].first > adj_sfcv[ix[ka]]->startparam() ||
      // 		  t1_t2[0].second < adj_sfcv[ix[ka]]->endparam())
      // 		{
      // 		}
      // 	    }
      // 	}
      adj_edg[ix[ka]] = shared_ptr<ftEdge>(new ftEdge(adj[ix[ka]]->getSurface(0),
						      adj_sfcv[ix[ka]], t1, t2));

      adj[ix[ka]]->addTrimEdge(adj_edg[ix[ka]]);
      bdedg[ka]->setReversed(true);
      adj_edg[ix[ka]]->connectTwin(bdedg[ka].get(), status);
  
      // Identify points from the adjacent regions lying outside the corresponding
      // trimming curve
      adj[ix[ka]]->extractOutOfEdge(adj_sfcv[ix[ka]],
				    (ix[ka] == 0) ? intcv1 : intcv2,
				    rad, approx_tol_, angtol, out[ix[ka]]);
    }

#ifdef DEBUG_BLEND
  std::ofstream of2("out_points.g2");
  for (int ka=0; ka<2; ++ka)
    {
      if (out[ka].size() > 0)
	{
	  of2 << "400 1 0 4 0 255 0 255" << std::endl;
	  of2 << out[ka].size() << std::endl;
	  for (size_t kr=0; kr<out[ka].size(); ++kr)
	    of2 << out[ka][kr]->getPoint() << std::endl;
	}
    }
#endif


  // Add out-points to the blend regions.
  // NB! This can be too simple in a more complex configuration.
  // Let's wait for the problem to turn up
  if (out[0].size() > 0 || out[1].size() > 0)
    {
      vector<RevEngPoint*> points;
      for (int ka=0; ka<2; ++ka)
	if (out[ka].size() > 0)
	  points.insert(points.end(), out[ka].begin(), out[ka].end());
      bool integrated = reg->addPointsToGroup(points, approx_tol_, angtol);
      if (!integrated)
	{
	  vector<HedgeSurface*> dummy_sfs;
	  vector<RevEngPoint*> dummy_pts;
	  for (int ka=0; ka<2; ++ka)
	    {
	      if (out[ix[ka]].size() > 0)
		{
		  vector<vector<RevEngPoint*> > out_1;
		  adj[ix[ka]]->connectedGroups(out[ix[ka]], out_1, false, dummy_pts);
		  for (size_t kr=0; kr<regions_.size(); ++kr)
		    if (regions_[kr].get() == adj[ix[ka]])
		      {
			surfaceExtractOutput((int)kr, out_1, dummy_sfs);
			break;
		      }
		}
	    }
	}
    }

  // Identify points associated to the blend region that should be
  // moved to the adjacent regions
  vector<int> adj_ix(4);
  adj_ix[0] = 3;
  adj_ix[1] = 1;
  adj_ix[2] = 0;
  adj_ix[3] = 2;
  vector<vector<RevEngPoint*> > move2adj(4);
  vector<RevEngPoint*> remain;
  vector<RevEngPoint*> regpoints = reg->getPoints();
  extractOutPoints(regpoints, elem, adj_ix, approx_tol_, angtol,
		   move2adj, remain);

#ifdef DEBUG_BLEND
  std::ofstream of2_3("in_points.g2");
  for (int ka=0; ka<4; ++ka)
    {
      if (move2adj[ka].size() > 0)
	{
	  of2_3 << "400 1 0 4 100 155 0 255" << std::endl;
	  of2_3 << move2adj[ka].size() << std::endl;
	  for (size_t kr=0; kr<move2adj[ka].size(); ++kr)
	    of2_3 << move2adj[ka][kr]->getPoint() << std::endl;
	}
    }
#endif

  int kx = 2*(1-constdir);
  for (int ka=kx; ka<=kx+1; ++ka)
    {
      if (move2adj[ka].size() > 0)
	{
	  remain.insert(remain.end(), move2adj[ka].begin(), move2adj[ka].end());
	  move2adj[ka].clear();
	}
    }

  for (int ka=0; ka<=1; ++ka)
    {
      int kb = 2*constdir+ka;
      if (move2adj[kb].size() > 0)
	{
	  reg->removePoints(move2adj[kb]);
	  bool integrated = adj[ix[ka]]->addPointsToGroup(move2adj[kb],
							  approx_tol_, angtol);
	  if (!integrated)
	    MESSAGE("RevEng::setBlendBoundaris. Missing adjacent surface");
	}
    }
  
#ifdef DEBUG_BLEND
  std::ofstream of3("updated_blend.g2");
  reg->writeRegionPoints(of3);
  adj[0]->writeRegionPoints(of3);
  adj[1]->writeRegionPoints(of3);
#endif
 int stop_break = 1;
}

//===========================================================================

// Service functionality for suitcaseCorner

bool getBlendRegMatches(vector<RevEngRegion*>& adj_blends,
			vector<vector<shared_ptr<ftEdge> > >& trim_edgs,
			vector<vector<pair<double, double> > >& par_lim,
			vector<vector<size_t> >& match)
{
#ifdef DEBUG_BLEND
  std::ofstream of("int_pt.g2");
#endif
  for (size_t ki=0; ki<adj_blends.size(); ++ki)
    for (size_t kj=ki+1; kj<adj_blends.size(); ++kj)
      {
	size_t kr, kh;
	for (kr=0; kr<trim_edgs[ki].size(); ++kr)
	  {
	    ftEdgeBase *edg1 = trim_edgs[ki][kr]->twin();
	    ftFaceBase *face1 = (edg1->geomEdge()) ?
	      edg1->geomEdge()->face() : 0;
	    if (!face1)
	      return false;
	    for (kh=0; kh<trim_edgs[kj].size(); ++kh)
	      {
		ftEdgeBase *edg2 = trim_edgs[kj][kh]->twin();
		ftFaceBase *face2 = (edg2->geomEdge()) ?
		  edg2->geomEdge()->face() : 0;
		if (!face2)
		  return false;
		if (face1 == face2)
		  break;
	      }
	    if (kh < trim_edgs[kj].size())
	      break;
	  }
	if (kr == trim_edgs[ki].size() || kh == trim_edgs[kj].size())
	  continue;

	shared_ptr<ParamCurve> cv1 = trim_edgs[ki][kr]->geomCurve();
	shared_ptr<ParamCurve> cv2 = trim_edgs[kj][kh]->geomCurve();
	double par1, par2, dist;
	Point ptc1, ptc2;
	ClosestPoint::closestPtCurves(cv1.get(), cv2.get(), par1, par2,
				      dist, ptc1, ptc2);
#ifdef DEBUG_BLEND
	of << "400 1 0 4 155 100 0 255" << std::endl;
	of << "1" << std::endl;
	of << ptc1 << std::endl;
	of << "400 1 0 4 155 100 0 255" << std::endl;
	of << "1" << std::endl;
	of << ptc2 << std::endl;
#endif
	par_lim[ki][kr] = std::make_pair(par1, dist);
	par_lim[kj][kh] = std::make_pair(par2, dist);
	vector<size_t> match0{ki, kr, kj, kh};
	match.push_back(match0);
      }
  return true;
}


bool getTrimEdgeMidpoint(vector<vector<shared_ptr<ftEdge> > >& trim_edgs,
			 vector<vector<pair<double, double> > >& par_lim,
			 vector<vector<pair<double, double> > >& midp)
{
  for (size_t ki=0; ki<par_lim.size(); ++ki)
    {
      int num = 0;
      double par = 0.0;
      for (size_t kj=0; kj<par_lim[ki].size(); ++kj)
	if (par_lim[ki][kj].second >= 0.0)
	  {
	    par += par_lim[ki][kj].first;
	    ++num;
	  }
      if (num != 2)
	return false;
      par /= (double)num;
      midp[ki].resize(par_lim[ki].size());
      for (size_t kj=0; kj<par_lim[ki].size(); ++kj)
	{
	  if (par_lim[ki][kj].second < 0.0)
	    {
	      midp[ki][kj] = std::make_pair(0.0, -1.0);
	      continue;
	    }
	  shared_ptr<ParamCurve> cv = trim_edgs[ki][kj]->geomCurve();
	  Point pt1 = cv->point(par_lim[ki][kj].first);
	  Point pt2 = cv->point(par);
	  double dist = pt1.dist(pt2);
	  midp[ki][kj] = std::make_pair(par,dist);
	  //int stop_break = 1;
	}
    }
  return true;
}

bool computeCornerBoundaryCurves(vector<RevEngRegion*>& adj_blends,
				 vector<vector<shared_ptr<ftEdge> > >& trim_edgs,
				 vector<vector<pair<double, double> > >& par_lim,
				 vector<vector<pair<double, double> > >& midp,
				 vector<vector<size_t> >& match, double tol,
				 vector<shared_ptr<CurveOnSurface> >& blend_bd,
				 vector<RevEngRegion*>& adjreg)
{
  double tol1 = 0.5*tol;

  // Want four boundary curves if possible
  vector<double> midd(midp.size());
  for (size_t ki=0; ki<midp.size(); ++ki)
    {
      double middist = 0.0;
      for (size_t kj=0; kj<midp[ki].size(); ++kj)
	middist = std::max(middist, midp[ki][kj].second);
      midd[ki] = middist;
    }
  vector<double> midd2(midd.begin(), midd.end());
  
  std::sort(midd2.begin(), midd2.end());
  double tol2;
  if (midd2.size() > 4 && midd2.size() < 2)
    return false;  // Currently not handled
  else if (midd2.size() == 4)
    tol2 = 2.0*midd2[3];
  else if (midd2.size() == 2)
    tol2 = tol1;
  else
    tol2 = 0.5*(midd2[0] + midd2[1]);
  
  adjreg.insert(adjreg.end(), adj_blends.begin(), adj_blends.end());

  vector<pair<double,double> > cvpar(adj_blends.size());
  vector<bool> parset(adj_blends.size(), false);

  // Start defining iso-parametric curves due to close intersection points
  for (size_t ki=0; ki<adj_blends.size(); ++ki)
    {
      if (midd[ki] <= tol1)
	{
	  double par;
	  size_t kj;
	  for (kj=0; kj<midp[ki].size() && midp[ki][kj].second < 0.0; ++kj);
	  par = midp[ki][kj].first;
	  cvpar[ki] = std::make_pair(par, par);
	  parset[ki] = true;
	}
    }

  // Continue with straight curve between opposite intersection points
  for (size_t ki=0; ki<adj_blends.size(); ++ki)
    {
      if (midd[ki] <= tol2)
	{
	  double par1, par2;
	  size_t kj, kr;
	  for (kj=0; kj<par_lim[ki].size() && par_lim[ki][kj].second < 0.0; ++kj);
	  par1 = par_lim[ki][kj].first;
	  for (kr=kj+1; kr<par_lim[ki].size() && par_lim[ki][kr].second < 0.0; ++kr);
	  par2 = par_lim[ki][kr].first;
	  cvpar[ki] = std::make_pair(par1, par2);
	  parset[ki] = true;
	}
    }

  if (adj_blends.size() != 3)
    return false;  // Waiting for a test case

  // Set iso-curve parameters for adjacent blends
  size_t kj;
  vector<pair<size_t, size_t> > missing;
  vector<double> div_par;
  for (kj=0; kj<parset.size(); ++kj)
    {
      if (parset[kj])
	{
	  size_t kh, kr;
	  for (kh=0; kh<par_lim[kj].size() && par_lim[kj][kh].second < 0.0; ++kh);
	  for (kr=kh+1; kr<par_lim[kj].size() && par_lim[kj][kr].second < 0.0; ++kr);

	  int ix1=-1, ix2=-1;
	  double par1, par2;
	  for (size_t kv=0; kv<match.size(); ++kv)
	    {
	      if (match[kv][0] == kj)
		{
		  if (match[kv][1] == kh)
		    {
		      ix1 = (int)match[kv][2];
		      par1 = par_lim[ix1][match[kv][3]].first;
		      missing.push_back(std::make_pair(match[kv][2], 1-match[kv][3]));
		      div_par.push_back(par1);
		    }		      
		  else if (match[kv][1] == kr)
		    {
		      ix2 = (int)match[kv][2];
		      par2 = par_lim[ix2][match[kv][3]].first;
		      missing.push_back(std::make_pair(match[kv][2], 1-match[kv][3]));
		      div_par.push_back(par2);
		    }
		}
	      else if (match[kv][2] == kj)
		{
		  if (match[kv][3] == kh)
		    {
		      ix1 = (int)match[kv][0];
		      par1 = par_lim[ix1][match[kv][1]].first;
		      missing.push_back(std::make_pair(match[kv][0], 1-match[kv][1]));
		      div_par.push_back(par1);
		    }
		  else if (match[kv][3] == kr)
		    {
		      ix2 = (int)match[kv][0];
		      par2 = par_lim[ix2][match[kv][1]].first;
		      missing.push_back(std::make_pair(match[kv][0], 1-match[kv][1]));
		      div_par.push_back(par2);
		    }
		}
	    }
	  if (ix1 < 0 || ix2 < 0)
	    return false;
	  cvpar[ix1] = std::make_pair(par1, par1);
	  parset[ix1] = true;
	  cvpar[ix2] = std::make_pair(par2, par2);
	  parset[ix2] = true;
	  break;
	}
    }
  if (kj == parset.size())
    return false;
  
  double eps = 1.0e-9;
  for (size_t ki=0; ki<adj_blends.size(); ++ki)
    {
      if (!parset[ki])
	return false;  // Missing information
      if (fabs(cvpar[ki].first-cvpar[ki].second) < eps)
	{
	  // Iso-parametric curve
	  size_t kj;
	  for (kj=0; kj<midp[ki].size() && midp[ki][kj].second < 0.0; ++kj);
	  double par = 0.5*(cvpar[ki].first+cvpar[ki].second);
	  shared_ptr<ParamCurve> cv = trim_edgs[ki][kj]->geomCurve();
	  shared_ptr<CurveOnSurface> sfcv =
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv);
	  shared_ptr<ParamCurve> pcv = sfcv->parameterCurve();
	  shared_ptr<ParamSurface> surf = sfcv->underlyingSurface();
	  if ((!pcv.get()) || (!surf.get()))
	    return false;
	  double val;
	  int dir;
	  bool isconst = sfcv->isConstantCurve(tol1, dir, val);
	  if (!isconst)
	    return false;
	  RectDomain dom = surf->containingDomain();
	  double pmin = (dir == 1) ? dom.umin() : dom.vmin();
	  double pmax = (dir == 1) ? dom.umax() : dom.vmax();
	  shared_ptr<CurveOnSurface> sfcv2(new CurveOnSurface(surf, (3-dir),
							      par, pmin, pmax,
							      -1));
	  blend_bd.push_back(sfcv2);
	}
      else
	{
	  // Straight curve in the parameter domain
	  size_t kj, kr;
	  for (kj=0; kj<par_lim[ki].size() && par_lim[ki][kj].second < 0.0; ++kj);
	  for (kr=kj+1; kr<par_lim[ki].size() && par_lim[ki][kr].second < 0.0; +kr);
	  double par1 = cvpar[ki].first;
	  double par2 = cvpar[ki].second;
	  shared_ptr<ParamCurve> cv1 = trim_edgs[ki][kj]->geomCurve();
	  shared_ptr<ParamCurve> cv2 = trim_edgs[ki][kr]->geomCurve();
	  shared_ptr<CurveOnSurface> sfcv1 =
			dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv1);
	  shared_ptr<CurveOnSurface> sfcv2 =
			dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv2);
	  shared_ptr<ParamCurve> pcv1 = sfcv1->parameterCurve();
	  shared_ptr<ParamCurve> pcv2 = sfcv2->parameterCurve();
	  if ((!pcv1.get()) || (!pcv2.get()))
	    return false;
	  Point ppt1 = pcv1->point(par1);
	  Point ppt2 = pcv2->point(par2);
	  shared_ptr<SplineCurve> pcv3(new SplineCurve(ppt1, ppt2));
	  shared_ptr<CurveOnSurface> sfcv3(new CurveOnSurface(sfcv1->underlyingSurface(),
							      pcv3, true));
	  bool space = sfcv3->ensureSpaceCrvExistence(tol1);
	  if (!space)
	    return false;
	  blend_bd.push_back(sfcv3);
	}
    }

#ifdef DEBUG_BLEND
  std::ofstream of1("space_bd.g2");
  for (size_t kr=0; kr<blend_bd.size(); ++kr)
    {
      shared_ptr<ParamCurve> space = blend_bd[kr]->spaceCurve();
      shared_ptr<ParamCurve> parcv = blend_bd[kr]->parameterCurve();
      space->writeStandardHeader(of1);
      space->write(of1);
    }
#endif

  if (missing.size() > 0)
    {
      // Add remaining curves
      size_t ki, kj;
      for (ki=0; ki<missing.size(); ++ki)
	for (kj=ki+1; kj<missing.size(); ++kj)
	  {
	    size_t kr;
	    for (kr=0; kr<match.size(); ++kr)
	      if (missing[ki].first == match[kr][0] &&
		  missing[ki].second == match[kr][1] &&
		  missing[kj].first == match[kr][2] &&
		  missing[kj].second == match[kr][3])
		break;
	    if (kr == match.size())
	      continue; // Should not happen
	    size_t ki1 = missing[ki].first, kr1 = missing[ki].second;
	    size_t kj1 = missing[kj].first, kh1 = missing[kj].second;
	    ftEdgeBase *edg = trim_edgs[ki1][kr1]->twin();
	    ftSurface *face = edg->geomEdge()->face()->asFtSurface();
	    if (!face)
	      continue;
	    shared_ptr<ParamSurface> surf = face->surface();
	    shared_ptr<ParamCurve> cv1 = trim_edgs[ki1][kr1]->geomCurve();
	    shared_ptr<ParamCurve> cv2 = trim_edgs[kj1][kh1]->geomCurve();
	    vector<Point> der1(2), der2(2);
	    cv1->point(der1, div_par[ki], 1);
	    cv2->point(der2, div_par[kj], 1);
	    double dd = der1[0].dist(der2[0]);
	    der1[1].normalize();
	    der2[1].normalize();
	    if ((der2[0]-der1[0])*der1[1] < 0.0)
	      der1[1] *= -1;
	    if ((der2[0]-der1[0])*der2[1] < 0.0)
	      der2[1] *= -1;
#ifdef DEBUG_BLEND
	    std::ofstream of3("missing_cvs.g2");
	    of3 << "400 1 0 4 255 0 0 255" << std::endl;
	    of3 << "1" << std::endl;
	    of3 << der1[0] << std::endl;
	    of3 << "410 1 0 4 255 0 0 255" << std::endl;
	    of3 << "1" << std::endl;
	    of3 << der1[0] << " " << der1[0]+0.3*dd*der1[1] << std::endl;
	    of3 << "400 1 0 4 255 0 0 255" << std::endl;
	    of3 << "1" << std::endl;
	    of3 << der2[0] << std::endl;
	    of3 << "410 1 0 4 255 0 0 255" << std::endl;
	    of3 << "1" << std::endl;
	    of3 << der2[0] << " " << der2[0]+0.3*dd*der2[1] << std::endl;
#endif
	    vector<double> knots(8, 0.0);
	    for (int ka=0; ka<4; ++ka)
	      knots[4+ka] = dd;
	    dd /= 3.0;
	    vector<double> coefs;
	    coefs.insert(coefs.end(), der1[0].begin(), der1[0].end());
	    Point cf1 = der1[0] + dd*der1[1];
	    coefs.insert(coefs.end(), cf1.begin(), cf1.end());
	    Point cf2 = der2[0] - dd*der2[1];
	    coefs.insert(coefs.end(), cf2.begin(), cf2.end());
	    coefs.insert(coefs.end(), der2[0].begin(), der2[0].end());
	    shared_ptr<ParamCurve> mcv(new SplineCurve(4, 4, &knots[0],
							&coefs[0], 3));
#ifdef DEBUG_BLEND
	    mcv->writeStandardHeader(of3);
	    mcv->write(of3);
	    surf->writeStandardHeader(of3);
	    surf->write(of3);
#endif

	    shared_ptr<SplineCurve> space_proj, par_proj;
	    CurveCreators::projectCurve(mcv, surf, tol1, space_proj,
					par_proj);

	    shared_ptr<CurveOnSurface> sfcv_proj(new CurveOnSurface(surf,
								    par_proj,
								    space_proj,
								    true, 1));
#ifdef DEBUG_BLEND
	    space_proj->writeStandardHeader(of3);
	    space_proj->write(of3);
#endif
	    blend_bd.push_back(sfcv_proj);
	    HedgeSurface *regface =  dynamic_cast<HedgeSurface*>(face);
	    if (!regface)
	      return false;
	    adjreg.push_back(regface->getRegion(0));
	  }
    }
  return true;
}

bool getCoonsBoundaryInfo(vector<shared_ptr<CurveOnSurface> >& blend_bd,
			  double tol,
			  vector<shared_ptr<ParamCurve> >& bdcvs,
			  vector<shared_ptr<ParamCurve> >& crosscvs)
{
  bdcvs.resize(blend_bd.size());
  crosscvs.resize(blend_bd.size());
  
  for (size_t ki=0; ki<blend_bd.size(); ++ki)
    {
      shared_ptr<ParamCurve> tmp_cv(blend_bd[ki]->spaceCurve()->clone());
      if (tmp_cv->instanceType() == Class_Circle)
	{
	  shared_ptr<Circle> circ = dynamic_pointer_cast<Circle,ParamCurve>(tmp_cv);
	  bdcvs[ki] = shared_ptr<SplineCurve>(circ->createNonRationalSpline(tol));
	}
      else if (tmp_cv->instanceType() == Class_Line)
	{
	  shared_ptr<Line> line = dynamic_pointer_cast<Line,ParamCurve>(tmp_cv);
	  bdcvs[ki] = shared_ptr<SplineCurve>(line->createSplineCurve());
	}
      else if (tmp_cv->instanceType() == Class_SplineCurve)
	{
	  shared_ptr<SplineCurve> cv = dynamic_pointer_cast<SplineCurve,ParamCurve>(tmp_cv);
	  if (cv->rational())
	    {
	      shared_ptr<ParamCurve> tmp_par = blend_bd[ki]->parameterCurve();
	      shared_ptr<ParamSurface> tmp_sf = blend_bd[ki]->underlyingSurface();
	      bdcvs[ki] =
		shared_ptr<SplineCurve>(CurveCreators::liftParameterCurve(tmp_par,
									 tmp_sf,
									 tol));
	    }
	  else
	    bdcvs[ki] = tmp_cv;
	}
      else
	return false;
      shared_ptr<CurveOnSurface> tmp_sfcv(blend_bd[ki]->clone());
      tmp_sfcv->setSpaceCurve(bdcvs[ki]);
      crosscvs[ki] = CreatorsUtils::createCrossTangent(*tmp_sfcv);
    }

  // Ensure correct direction of cross tangent curves
  for (size_t ki=0; ki<crosscvs.size(); ++ki)
    {
      size_t kj = (ki == crosscvs.size()-1) ? 0 : ki+1;
      Point ctan = crosscvs[ki]->point(crosscvs[ki]->endparam());
      double tdel = bdcvs[kj]->endparam() - bdcvs[kj]->startparam();
      Point pos1 = bdcvs[kj]->point(bdcvs[kj]->startparam());
      Point pos2 = bdcvs[kj]->point(bdcvs[kj]->startparam() + 0.1*tdel);
      Point vec = pos2 - pos1;
      if (ctan*vec <  0.0)
	{
	  shared_ptr<SplineCurve> cross =
	    dynamic_pointer_cast<SplineCurve,ParamCurve>(crosscvs[ki]);
	  for (auto it=cross->coefs_begin(); it!=cross->coefs_end(); ++it)
	    (*it) *= -1.0;
	}
    }
  return true;
}

void pairOfRegionEdges(RevEngRegion* blendreg, HedgeSurface *hedge,
		       shared_ptr<ParamCurve>& bdcv1,
		       shared_ptr<CurveOnSurface>& bdcv2,
		       RevEngRegion* adjreg, Point& mid, double tol)
{
  int stat = 0;
  double eps = 1.0e-9;
  shared_ptr<ftEdge> blend_edge(new ftEdge(hedge, bdcv1,
					   bdcv1->startparam(),
					   bdcv1->endparam()));
  blendreg->addTrimEdge(blend_edge);

  int pdir;
  double pval;
  if (adjreg->hasBlendEdge() && bdcv2->isConstantCurve(tol, pdir, pval))
    {
      shared_ptr<ElementarySurface> adj_elem =
	dynamic_pointer_cast<ElementarySurface,ParamSurface>(adjreg->getSurface(0)->surface());
      if (adj_elem.get())
	{
	  RectDomain adjdom = adj_elem->getParameterBounds();
	  double regdom[4];
	  adjreg->getDomain(regdom);
	  // double tmin = regdom[2*(pdir-1)];
	  // double tmax = regdom[2*(pdir-1)+1];
	  double upar, vpar, dist;
	  Point close;
	  adj_elem->closestPoint(mid, upar, vpar, close, dist, eps);
	  if (pdir == 1)
	    {
	      if (upar > pval)
		adj_elem->setParameterBounds(adjdom.umin(), adjdom.vmin(),
					     pval, adjdom.vmax());
	      else
		adj_elem->setParameterBounds(pval, adjdom.vmin(),
					     adjdom.umax(), adjdom.vmax());
	    }
	  else
	    {
	      if (vpar > pval)
		adj_elem->setParameterBounds(adjdom.umin(), adjdom.vmin(), 
					     adjdom.umax(), pval);
	      else
		adj_elem->setParameterBounds(adjdom.umin(), pval, 
					     adjdom.umax(), adjdom.vmax());
	    }
	  // if (fabs(tmax-pval) > fabs(pval-tmin))
	  //   {
	  //     if (pdir == 1)
	  // 	adj_elem->setParameterBounds(pval, adjdom.vmin(),
	  // 				     adjdom.umax(), adjdom.vmax());
	  //     else
	  // 	adj_elem->setParameterBounds(adjdom.umin(), pval, 
	  // 				     adjdom.umax(), adjdom.vmax());
	  //   }
	  // else
	  //   {
	  //     if (pdir == 1)
	  // 	adj_elem->setParameterBounds(adjdom.umin(), adjdom.vmin(),
	  // 				     pval, adjdom.vmax());
	  //     else
	  // 	adj_elem->setParameterBounds(adjdom.umin(), adjdom.vmin(), 
	  // 				     adjdom.umax(), pval);
	  //   }
	}
      adjreg->adaptEdges();
    }

  HedgeSurface *other_hedge = adjreg->getSurface(0);
  shared_ptr<ftEdge> other_edge(new ftEdge(other_hedge, bdcv2,
					   bdcv2->startparam(),
					   bdcv2->endparam()));
  adjreg->addTrimEdge(other_edge);
  other_edge->setReversed(true);
  blend_edge->connectTwin(other_edge.get(), stat);
}

//===========================================================================
void RevEng::extractOutPoints(vector<RevEngPoint*>& points, shared_ptr<ParamSurface> surf,
			      vector<int>& cv_ix,
			      double tol, double angtol,
			      vector<vector<RevEngPoint*> >& move2adj,
			      vector<RevEngPoint*>& remain)
//===========================================================================
{
  double eps = 1.0e-9;
  double maxd, avd;
  int num_in, num2_in;
  vector<double> parvals;
  vector<pair<double,double> > dist_ang;
  vector<RevEngPoint*> inpt, outpt;
  RevEngUtils::distToSurf(points.begin(), points.end(), surf,
			  tol, maxd, avd, num_in, num2_in, inpt, outpt,
			  parvals, dist_ang, angtol);

  RectDomain dom = surf->containingDomain();
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      int bd=-1, bd2=-1;
      Vector2D par = Vector2D(parvals[2*ki], parvals[2*ki+1]);
      if (dom.isOnBoundary(par, eps, bd, bd2))
	{
	  if (bd2 >= 0)
	    {
	      // Check best fit
	    }
	  int move_ix = (bd == 0) ? cv_ix[3] :
	    ((bd == 1) ? cv_ix[1] : ((bd == 2) ? cv_ix[0] : cv_ix[2]));
	  move2adj[move_ix].push_back(points[ki]);
	}
      else
	remain.push_back(points[ki]);
    }
}

//===========================================================================
bool RevEng::suitcaseCorner(vector<RevEngRegion*>& adj_blends,
			    RevEngEdge* rev_edge)
//===========================================================================
{
  vector<vector<shared_ptr<ftEdge> > > trim_edgs(adj_blends.size());
  vector<vector<pair<double, double> > > par_lim(adj_blends.size());
  for (size_t ki=0; ki<adj_blends.size(); ++ki)
    {
      trim_edgs[ki] = adj_blends[ki]->getTrimEdges();
      par_lim[ki].resize(trim_edgs[ki].size());
      for (size_t kj=0; kj<par_lim[ki].size(); ++kj)
	par_lim[ki][kj] = std::make_pair(-1.0, -1.0);   // Dummy
    }

  vector<vector<size_t> > match;
  bool found1 = getBlendRegMatches(adj_blends, trim_edgs, par_lim, match);
  if (!found1)
    return false;

  vector<vector<pair<double, double> > > midp(par_lim.size());
  bool found2 = getTrimEdgeMidpoint(trim_edgs, par_lim, midp);
  if (!found2)
    return false;

  // Compute boundary curves
  vector<shared_ptr<CurveOnSurface> > blend_bd;
  vector<RevEngRegion*> adjreg;
  bool found3 = computeCornerBoundaryCurves(adj_blends, trim_edgs, par_lim,
					    midp, match, approx_tol_,
					    blend_bd, adjreg);
  if (!found3)
    return false;
					    

  // Ensure consistent curve sequence and direction
  vector<shared_ptr<CurveOnSurface> > blend_bd0(blend_bd.begin(), blend_bd.end());
  RevEngUtils::setLoopSeq(blend_bd);
  
#ifdef DEBUG_BLEND
  std::ofstream of4("space_bd2.g2");
  for (size_t kr=0; kr<blend_bd.size(); ++kr)
    {
      shared_ptr<ParamCurve> space = blend_bd[kr]->spaceCurve();
      space->writeStandardHeader(of4);
      space->write(of4);
      of4 << "400 1 0 4 255 0 0 255" << std::endl;
      of4 << "1" << std::endl;
      of4 << space->point(space->startparam()) << std::endl;
    }
#endif

  // Ensure non-rational spline boundary curves and extract cross parameter curves
  double tol1 = 0.5*approx_tol_;
  vector<shared_ptr<ParamCurve> > crosscvs;
  vector<shared_ptr<ParamCurve> > bdcvs;
  bool found4 = getCoonsBoundaryInfo(blend_bd, tol1, bdcvs, crosscvs);
  if (!found4)
    return false;

#ifdef DEBUG_BLEND
  std::ofstream of4_2("space_bd3.g2");
  for (size_t kr=0; kr<bdcvs.size(); ++kr)
    {
      bdcvs[kr]->writeStandardHeader(of4_2);
      bdcvs[kr]->write(of4_2);
    }
#endif
  
  // Create surface
  shared_ptr<SplineSurface> coons;
  if (bdcvs.size() == 4)
    {
      coons =
	shared_ptr<SplineSurface>(CoonsPatchGen::createCoonsPatch(bdcvs,
								  crosscvs,
								  tol1,
								  anglim_));
      RevEngUtils::smoothSurf(coons, 2);
    }
  else
    {
      MESSAGE("RevEng::suitcaseCorner. Only 4 boundary curves are supported.");
      return false;
    }

  if (!coons.get())
    return false;
  
#ifdef DEBUG_BLEND
  std::ofstream of5("coons_patch.g2");
  coons->writeStandardHeader(of5);
  coons->write(of5);
#endif

  // Set trimming curves
  double eps = 1.0e-9;
  CurveLoop cvloop = SurfaceTools::outerBoundarySfLoop(coons, eps);
  vector<shared_ptr<ParamCurve> > loopcvs = cvloop.getCurves();

  // Define match between boundary curves of corner surface and
  // adjacent surfaces
  vector<int> cv_ix(loopcvs.size());
  for (size_t ki=0; ki<loopcvs.size(); ++ki)
    {
      Point mid = loopcvs[ki]->point(0.5*(loopcvs[ki]->startparam()+
					  loopcvs[ki]->endparam()));
      int ix = -1;
      double mindist = std::numeric_limits<double>::max();
      for (size_t kj=0; kj<blend_bd0.size(); ++kj)
	{
	  double tpar, dist;
	  Point close;
	  blend_bd0[kj]->closestPoint(mid, blend_bd0[kj]->startparam(),
				     blend_bd0[kj]->endparam(), tpar,
				     close, dist);
	  if (dist < mindist)
	    {
	      ix = (int)kj;
	      mindist = dist;
	    }
	}
      cv_ix[ki] = ix;
    }

  // Move points from blend region as appropriate. First fetch blend points
    vector<RevEngRegion*> blend_regs;
  rev_edge->getAllBlendRegs(blend_regs);
  vector<RevEngPoint*> blend_pts;
  for (size_t ki=0; ki<blend_regs.size(); ++ki)
    blend_pts.insert(blend_pts.end(), blend_regs[ki]->pointsBegin(),
		     blend_regs[ki]->pointsEnd());

  // double angtol = 5.0*anglim_;
  // double maxd, avd;
  // int num_in, num2_in;
  // vector<double> parvals;
  // vector<pair<double,double> > dist_ang;
  // vector<RevEngPoint*> inpt, outpt;
  // RevEngUtils::distToSurf(blend_pts.begin(), blend_pts.end(), coons,
  // 			  approx_tol_, maxd, avd, num_in, num2_in, inpt, outpt,
  // 			  parvals, dist_ang, angtol);

  double angtol = 5.0*anglim_;
  vector<vector<RevEngPoint*> > move2adj(4);
  vector<RevEngPoint*> remain;
  extractOutPoints(blend_pts, coons, cv_ix, approx_tol_, angtol,
		   move2adj, remain);
  // RectDomain dom = coons->containingDomain();
  // for (size_t ki=0; ki<blend_pts.size(); ++ki)
  //   {
  //     int bd=-1, bd2=-1;
  //     Vector2D par = Vector2D(parvals[2*ki], parvals[2*ki+1]);
  //     if (dom.isOnBoundary(par, eps, bd, bd2))
  // 	{
  // 	  if (bd2 >= 0)
  // 	    {
  // 	      // Check best fit
  // 	    }
  // 	  int move_ix = (bd == 0) ? cv_ix[3] :
  // 	    ((bd == 1) ? cv_ix[1] : ((bd == 2) ? cv_ix[0] : cv_ix[2]));
  // 	  move2adj[move_ix].push_back(blend_pts[ki]);
  // 	}
  //     else
  // 	remain.push_back(blend_pts[ki]);
  //   }
  // if (remain.size() == 0 && remain.size() < blend_pts.size())
  //   return false;

  bool OK;
  for (size_t ki=0; ki<move2adj.size(); ++ki)
    if (move2adj[ki].size() > 0)
      OK = adjreg[ki]->addPointsToGroup(move2adj[ki], approx_tol_, angtol);

  // Delete current blend regions
  vector<HedgeSurface*> prev_sfs;
  for (size_t ki=0; ki<blend_regs.size(); ++ki)
    {
      //blend_regs[ki]->removeFromAdjacent();
      blend_regs[ki]->setRemove();
    }
  
  // Create blend region
  shared_ptr<RevEngRegion> blendreg(new RevEngRegion(classification_type_,
						     edge_class_type_,
						     remain));  
  blendreg->setRegionAdjacency();
  regions_.push_back(blendreg);
  shared_ptr<HedgeSurface> hedge;
  shared_ptr<ParamSurface> blend_tmp = coons;
  blendreg->setAssociatedSurface(blend_tmp, approx_tol_, angtol,
				 min_point_region_, hedge);
  rev_edge->setBlendRegSurf(blendreg.get());
  blendreg->setBlendEdge(rev_edge);
  if (hedge.get())
    surfaces_.push_back(hedge);
  
  // Add edges to regions
  RectDomain dom = coons->containingDomain();
  double umid = 0.5*(dom.umin()+dom.umax());
  double vmid = 0.5*(dom.vmin()+dom.vmax());
  Point mid = coons->ParamSurface::point(umid, vmid);
  for (size_t ki=0; ki<loopcvs.size(); ++ki)
    {
      shared_ptr<CurveOnSurface> other = blend_bd0[cv_ix[ki]];
      RevEngRegion *adj = adjreg[cv_ix[ki]];
      pairOfRegionEdges(blendreg.get(), hedge.get(), loopcvs[ki], other, adj, mid,
			tol1);
    }

  // Move points from adjacent regions to corner blend region as appropriate
  for (size_t ki=0; ki<adjreg.size(); ++ki)
    {
      blendreg->growInDomain(adjreg[ki], approx_tol_, angtol);
      if (adjreg[ki]->numPoints() == 0)
	adjreg[ki]->setRemove();
    }
  
  // Remove obsolete edge information
  RevEngRegion* adj[2];
  rev_edge->getAdjacent(adj[0], adj[1]);
  size_t kr;
  for (int ka=0; ka<2; ++ka)
    {
      for (kr=0; kr<adj_blends.size(); ++kr)
  	if (adj_blends[kr] == adj[ka])
  	  break;
      if (kr == adj_blends.size())
  	{
  	  rev_edge->eraseAdjacent(adj[ka]);
  	  break;
  	}
    }
  rev_edge->eraseCurves();
  
  return true;
}

//===========================================================================
bool RevEng::createBlendSurface(int ix)
//===========================================================================
{
  double angtol = 5.0*anglim_;
  double diag = bbox_.low().dist(bbox_.high());
  double eps = 1.0e-6;
  double blendlim = std::min(0.1*diag, 30.0*mean_edge_len_); //50.0*mean_edge_len_;
  double lenlim = 10.0*mean_edge_len_; //blendlim;
  
  // Adjacent regions
  RevEngRegion* adj[2];
  edges_[ix]->getAdjacent(adj[0], adj[1]);

  if ((!adj[0]) || (!adj[1]))
    return false;
  if ((!adj[0]->hasSurface()) || (!adj[1]->hasSurface()))
    return false;  // Something wrong

#ifdef DEBUG_BLEND
  std::ofstream of1("adj_groups.g2");
  adj[0]->writeRegionPoints(of1);
  adj[1]->writeRegionPoints(of1);
#endif

  if (adj[0]->hasAssociatedBlend() || adj[1]->hasAssociatedBlend())
    return false;

  // Update intersection curve if necessary
  bool updated_edge = false;
  if (edges_[ix]->getExtendCount() > 0)
    {
      vector<shared_ptr<RevEngRegion> > added_regions;
      vector<vector<RevEngPoint*> > extract_groups;
      vector<HedgeSurface*> out_sfs;
      updated_edge = edges_[ix]->extendCurve(int_tol_, approx_tol_, anglim_, diag, lenlim,
					    blendlim, added_regions, extract_groups, out_sfs);
      if (extract_groups.size() > 0 || out_sfs.size() > 0)
	surfaceExtractOutput(-1, extract_groups, out_sfs);
      for (size_t kj=0; kj<added_regions.size(); ++kj)
	{
	  added_regions[kj]->setRegionAdjacency();
	  regions_.push_back(added_regions[kj]);
	}

      if (updated_edge)
	{
	  for (int ka=ix+1; ka<(int)edges_.size(); ++ka)
	    {
	      // Check for overlap (simplified version)
	      bool embedded = edges_[ix]->contains(edges_[ka].get(), approx_tol_);
	      if (embedded)
		{
		  bool done = edges_[ix]->integrate(edges_[ka].get());
		  if (done)
		    edges_.erase(edges_.begin()+ka);
		}
	    }
	}
    }

  if (updated_edge == false && edges_[ix]->getSurfChangeCount() > 0)
    updated_edge = edges_[ix]->updateCurve(int_tol_, approx_tol_, diag);
  if (!updated_edge)
    edges_[ix]->fixMismatchCurves(approx_tol_);  // @@@ Should this be done always?

  if (edges_[ix]->getType() == NOT_BLEND)
    return false;  // Should not create blend surfaces
  
  // Regions in blend area
  vector<RevEngRegion*> blend_regs;
  edges_[ix]->getAllBlendRegs(blend_regs);
  
  // Intersection curve
  vector<shared_ptr<CurveOnSurface> > cvs;
  edges_[ix]->getCurve(cvs, true);

  vector<Point> der(2);
  cvs[0]->point(der, 0.5*(cvs[0]->startparam()+cvs[0]->endparam()), 1);

#ifdef DEBUG_BLEND
  for (size_t ki=0; ki<cvs.size(); ++ki)
    {
      shared_ptr<ParamCurve> tmp_cv = cvs[ki]->spaceCurve();
      tmp_cv->writeStandardHeader(of1);
      tmp_cv->write(of1);
    }
#endif

  // Width
  double width = edges_[ix]->getDistance();

  bool out1 = false, out2 = false;
  edges_[ix]->getOuterInfo(out1, out2);

  shared_ptr<ParamSurface> surf1 = adj[0]->getSurface(0)->surface();
  shared_ptr<ParamSurface> surf2 = adj[1]->getSurface(0)->surface();
  ClassType classtype1 = surf1->instanceType();
  ClassType classtype2 = surf2->instanceType();
  if (classtype1 != Class_Plane && classtype1 != Class_Cylinder &&
      classtype1 != Class_Cone)
    return false;
  if (classtype2 != Class_Plane && classtype2 != Class_Cylinder &&
      classtype2 != Class_Cone)
    return false;
  shared_ptr<ElementarySurface> elem1 =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf1);
  shared_ptr<ElementarySurface> elem2 =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf2);
  vector<shared_ptr<ElementarySurface> > elem_sfs(2);
  elem_sfs[0] = elem1;
  elem_sfs[1] = elem2;
  if (!(elem1.get() && elem2.get()))
    return false;
  Point dir1 = elem1->direction();
  Point norm1 = adj[0]->getMeanNormalTriang();
  if (elem1->instanceType() == Class_Plane && dir1*norm1 < 0.0)
    dir1 *= -1;
  Point dir2 = elem2->direction();
  Point norm2 = adj[1]->getMeanNormalTriang();
  if (elem2->instanceType() == Class_Plane && dir2*norm2 < 0.0)
    dir2 *= -1;

  if (classtype1 == Class_Plane && classtype2 == Class_Plane)
    {
      // Create cylinder
    }
  else if (classtype1 == Class_Plane || classtype2 == Class_Plane)
    {
      // Create torus
    }
  else if ((classtype1 == Class_Cylinder && classtype2 == Class_Cone) ||
	   (classtype2 == Class_Cylinder && classtype1 == Class_Cone))
    {
      // Create torus
    }
  else
    return false;  // Not supported

  // Collect points from adjacent surfaces
  vector<vector<RevEngPoint*> > blend_pts(3);
  double tmin = cvs[0]->startparam();
  double tmax = cvs[0]->endparam();
  for (int ka=0; ka<2; ++ka)
    adj[ka]->getNearPoints(cvs[0], tmin, tmax, 2*width, angtol, blend_pts[ka]);
 
  // Collect points from associated blend regions
  for (size_t ki=0; ki<blend_regs.size(); ++ki)
    blend_pts[2].insert(blend_pts[2].end(), blend_regs[ki]->pointsBegin(),
			blend_regs[ki]->pointsEnd());

#ifdef DEBUG_BLEND
  std::ofstream of2("blend_pts.g2");
  for (size_t ki=0; ki<3; ++ki)
    {
      if (blend_pts[ki].size() > 0)
	{
	  of2 << "400 1 0 4 0 255 0 255" << std::endl;
	  of2 << blend_pts[ki].size() << std::endl;
	  for (size_t kr=0; kr<blend_pts[ki].size(); ++kr)
	    of2 << blend_pts[ki][kr]->getPoint() << std::endl;
	}
    }
#endif

  shared_ptr<ElementarySurface> blend_surf;
  double ylen = 0.0;
  double axis_ang = dir1.angle(dir2);
  int ldir = -1;
  if ((classtype1 == Class_Plane && classtype2 == Class_Plane) ||
      ((classtype1 == Class_Plane || classtype2 == Class_Plane) &&
       fabs(0.5*M_PI - axis_ang) <= angtol))
    {
      // Create cylinder
      Point lin1, lin2;
      Point dir1_2 = dir1, dir2_2 = dir2;
      if (classtype1 == Class_Plane)
	lin1 = der[1].cross(dir1);
      else
	{
	  double clo_u, clo_v, clo_dist;
	  Point clo;
	  surf1->closestPoint(der[0], clo_u, clo_v, clo, clo_dist, eps);
	  surf1->normal(dir1_2, clo_u, clo_v);
	  lin1 = der[1].cross(dir1_2);
	}
      
      if (classtype2 == Class_Plane)
	lin2 = der[1].cross(dir2);
      else
	{
	  double clo_u, clo_v, clo_dist;
	  Point clo;
	  surf2->closestPoint(der[0], clo_u, clo_v, clo, clo_dist, eps);
	  surf2->normal(dir2_2, clo_u, clo_v);
	  lin2 = der[1].cross(dir2_2);
	}
      
      lin1.normalize();
      if (lin1*dir2_2 < 0.0)
	lin1 *= -1;
      lin2.normalize();
      if (lin2*dir1_2 < 0.0)
	lin2 *= -1;
      int sgn = (out1 && out2) ? 1 : -1;
      blend_surf = createCylinderBlend(blend_pts, width, der[0],
				       der[1], lin1, lin2, sgn);
      double radius = blend_surf->radius(0.0, 0.0);
      //double ang = dir1.angle(dir2);
      Point centre = blend_surf->location();
      double xlen = der[0].dist(centre);
      ylen = sqrt(xlen*xlen - radius*radius);
      
      ldir = 1;
    }
  else
    {
      int sgn = 1;
      if ((classtype1 == Class_Plane && dir1*elem1->direction() < 0.0) ||
	  (classtype2 == Class_Plane && dir2*elem2->direction() < 0.0))
	sgn = -1;

      blend_surf = torusBlend(blend_pts, cvs[0], elem1, elem2, width,
			      out1, out2, sgn);

      double Rrad = blend_surf->radius(0.0, 0.0);
      Point centre = blend_surf->location();
      double xlen = der[0].dist(centre);
      ylen = fabs(xlen - Rrad);

      ldir = 1;
    }
  
  if (!blend_surf.get())
    {
#ifdef DEBUG_BLEND
      std::cout << "No blend_surf" << std::endl;
#endif
      return false;
    }

#ifdef DEBUG_BLEND
  for (int ka=0; ka<2; ++ka)
    {
      int num = adj[ka]->numPoints();
      for (int kb=0; kb<num; ++kb)
	if (adj[ka]->getPoint(kb)->region() != adj[ka])
	  std::cout << "Inconsistent region pointers, pre sortBlendPoints: " << ka << " " << kb << std::endl;
    }
#endif
  
  // Associate blend points with the appropriate surface (region)
  // First remove blend points from the adjacent surfaces
  vector<vector<RevEngPoint*> > out_pts(2);
  for (int ka=0; ka<2; ++ka)
    {
      adj[ka]->sortBlendPoints(blend_pts[ka], cvs, ylen, true, out_pts[ka]);
      blend_pts[2].insert(blend_pts[2].end(), out_pts[ka].begin(),
			  out_pts[ka].end());
    }
#ifdef DEBUG_BLEND
  for (int ka=0; ka<2; ++ka)
    {
      int num = adj[ka]->numPoints();
      for (int kb=0; kb<num; ++kb)
	if (adj[ka]->getPoint(kb)->region() != adj[ka])
	  std::cout << "Inconsistent region pointers, post sortBlendPoints1: " << ka << " " << kb << std::endl;
    }
#endif
      
  // Move blend points to the adjacent surfaces if appropriate
  vector<vector<RevEngPoint*> > in_pts(2);
  adj[0]->sortBlendPoints(blend_pts[2], cvs, ylen, adj[1],
			  in_pts[0], in_pts[1]);

#ifdef DEBUG_BLEND
  for (int ka=0; ka<2; ++ka)
    {
      int num = adj[ka]->numPoints();
      for (int kb=0; kb<num; ++kb)
	if (adj[ka]->getPoint(kb)->region() != adj[ka])
	  std::cout << "Inconsistent region pointers, post sortBlendPoints2: " << ka << " " << kb << std::endl;
    }
#endif
      
  // Update adjacent surfaces with modified collection of points
  for (int ka=0; ka<2; ++ka)
    {
      adj[ka]->updateWithPointsInOut(out_pts[ka], in_pts[ka], approx_tol_, angtol);
    }

  // Delete current blend regions
  vector<HedgeSurface*> prev_sfs;
  for (size_t ki=0; ki<blend_regs.size(); ++ki)
    {
      blend_regs[ki]->removeFromAdjacent();
      blend_regs[ki]->clearRegionAdjacency();
      int num_sfs = blend_regs[ki]->numSurface();
      for (int ka=0; ka<num_sfs; ++ka)
	prev_sfs.push_back(blend_regs[ki]->getSurface(ka));
    }
  edges_[ix]->clearBlendRegions();
  if (blend_regs.size() > 0)
    {
      int dummy_ix = -1;
      updateRegionsAndSurfaces(dummy_ix, blend_regs, prev_sfs);
    }
  
  if (blend_pts[2].size() == 0)
    {
#ifdef DEBUG_BLEND
      std::cout << "No points for blend_surf" << std::endl;
#endif
      return false;
    }
  
  // Define blend region
  shared_ptr<RevEngRegion> blendreg(new RevEngRegion(classification_type_,
						     edge_class_type_,
						     blend_pts[2]));
#ifdef DEBUG_BLEND
  std::ofstream of3("updated_blend_pts.g2");
  adj[0]->writeRegionPoints(of3);
  adj[1]->writeRegionPoints(of3);
  blendreg->writeRegionPoints(of3);
#endif
  
  blendreg->setRegionAdjacency();
  regions_.push_back(blendreg);
  shared_ptr<HedgeSurface> hedge;
  shared_ptr<ParamSurface> blend_surf_tmp = blend_surf;
  blendreg->setAssociatedSurface(blend_surf_tmp, approx_tol_, angtol,
				 min_point_region_, hedge);
  if (hedge.get())
    surfaces_.push_back(hedge);

  // // Check blend points
  // vector<vector<RevEngPoint*> > out_groups;
  // blendreg->extractOutPoints(ldir, approx_tol_, angtol,
  // 			     1.1*angtol, out_groups);
  // if (out_groups.size() > 0)
  //   {
  //     vector<HedgeSurface*> dummy_sfs;
  //     surfaceExtractOutput(regions_.size()-1, out_groups, dummy_sfs);
  //   }
   
  // Update edge with blend region (surface)
  edges_[ix]->setBlendRegSurf(blendreg.get());
  blendreg->setBlendEdge(edges_[ix].get());
  for (int ka=0; ka<2; ++ka)
    adj[ka]->updateRegionAdjacency();
  blendreg->setRegionAdjacency();
  //edges_[ix]->setAltRadius(radius1);

  for (int ka=0; ka<2; ++ka)
    {
      vector<vector<RevEngPoint*> > sep_groups;
      adj[ka]->splitRegion(sep_groups);
      if (sep_groups.size() > 0)
	{
	  size_t kh;
	  for (kh=0; kh<regions_.size(); ++kh)
	    if (regions_[kh].get() == adj[ka])
	      break;
	  vector<HedgeSurface*> dummy_sfs;
	  surfaceExtractOutput((kh<regions_.size()) ? (int)kh : 0, sep_groups,
			       dummy_sfs);
	}
   }

  return true;
}

//===========================================================================
double
RevEng::computeTorusRadius(vector<vector<RevEngPoint*> >& blend_pts,
			   shared_ptr<CurveOnSurface>& cv,
			   const Point& locp, const Point& normal,
			   shared_ptr<ElementarySurface> rotational,
			   double width, bool plane_out, bool rot_out)
//===========================================================================
{
  double alpha = 0;
  shared_ptr<Cone> cone = dynamic_pointer_cast<Cone,ElementarySurface>(rotational);
  if (cone.get())
    alpha = cone->getConeAngle();
  Point axis = rotational->direction();
  
  // Compute minor radius of torus
  double lrad = 0.0;
  int lnmb = 0;
  double beta = 0.5*M_PI + alpha;
  //double phi = 0.5*M_PI - alpha;
  double fac = 1.0/sin(0.5*beta);
  double div = fac - 1.0;
  Point loc = rotational->location();
  double rd = (locp-loc)*axis;
  Point centr = loc + rd*axis;
  for (size_t kj=0; kj<blend_pts.size(); ++kj)
    {
      for (size_t ki=0; ki<blend_pts[kj].size(); ++ki)
	{
	  Vector3D xyz = blend_pts[kj][ki]->getPoint();
	  Point ptpos(xyz[0], xyz[1], xyz[2]);

	  // Localize with respect to the guiding curve(s)
	  double tpar, dist;
	  Point close;
	  cv->closestPoint(ptpos, cv->startparam(), cv->endparam(),
				tpar, close, dist);

	  // Define line in the point between the adjacent surfaces
	  vector<Point> der(2);
	  cv->point(der, tpar, 1);
	  der[1].normalize();
	  Point dir1 = der[1].cross(normal);
	  Point vec1 = der[0] - centr;
	  if ((dir1*vec1 < 0.0 && (!rot_out)) || (dir1*vec1 > 0.0 && rot_out))
	    dir1 *= -1;
	  Point dir2 = (plane_out) ? normal : -normal;
	  Point dir3;
	  if (alpha > 0)
	    {
	      Matrix3D mat;
	      mat.setToRotation(-alpha, der[1][0], der[1][1], der[1][2]);
	      Vector3D dir2_2(dir2[0], dir2[1], dir2[2]);
	      Vector3D dir3_2 = mat*dir2_2;
	      dir3 = Point(dir3_2[0], dir3_2[1], dir3_2[2]);
	    }
	  else
	    dir3 = dir2;
	  Point dir4 = 0.5*(dir1 + dir3);
#ifdef DEBUG_BLEND
	  std::ofstream of("midlin.g2");
	  of << "410 1 0 4 0 255 0 255" << std::endl;
	  of << "1" << std::endl;
	  of << der[0] << " " << der[0]+dir1 << std::endl;
	  of << "410 1 0 4 255 0 0 255" << std::endl;
	  of << "1" << std::endl;
	  of << der[0] << " " << der[0]+dir4 << std::endl;
	  of << "410 1 0 4 0 0 255 255" << std::endl;
	  of << "1" << std::endl;
	  of << der[0] << " " << der[0]+dir2 << std::endl;
	      
#endif
	  shared_ptr<Line> line(new Line(der[0], dir4));
	  line->setParameterInterval(-2*width, 2*width);
	  double lpar, ldist;
	  Point lclose;
	  line->closestPoint(ptpos, line->startparam(), line->endparam(),
			     lpar, lclose, ldist);
	  if (ldist <= approx_tol_)
	    {
	      double dd2 = close.dist(lclose);
	      double xlen = dd2/div;
	      lrad += xlen;
	      lnmb++;
	    }
	  int stop_break = 1;
	}
    }

  if (lnmb > 0)
    lrad /= (double)lnmb;

  return lrad;
}

//===========================================================================
void RevEng::getTorusParameters(shared_ptr<ElementarySurface> planar,
				shared_ptr<ElementarySurface> rotational,
				double radius, int sgn1, int sgn2, double& Rrad, 
				Point& centre, Point& normal, Point& Cx)
//===========================================================================
{
  Point locp = planar->location();
  normal = planar->direction();
  Point loc = rotational->location();
  Point axis = rotational->direction();
  double rd = (locp - loc)*axis;
  Point centre0 = loc + rd*axis;
  centre = centre0 - sgn1*radius*normal;
  Rrad = rotational->radius(0.0, rd);
  double alpha = 0.0;
  if (rotational->instanceType() == Class_Cone)
    alpha = ((Cone*)(rotational.get()))->getConeAngle();
  double phi = 0.5*M_PI - alpha;
  double sd = radius/sin(phi);
  Cx = rotational->direction2();
  Rrad += (sgn2*sd);
}

//===========================================================================
shared_ptr<Torus>
RevEng::torusBlend(vector<vector<RevEngPoint*> >& blend_pts,
		   vector<shared_ptr<CurveOnSurface> >& cvs,
		   const Point& locp, const Point& normal,
		   shared_ptr<ElementarySurface> rotational,
		   double width, bool plane_out, bool rot_out)
//===========================================================================
{
  shared_ptr<Torus> torus;

  double alpha = 0;
  shared_ptr<Cone> cone = dynamic_pointer_cast<Cone,ElementarySurface>(rotational);
  if (cone.get())
    alpha = cone->getConeAngle();
  Point axis = rotational->direction();
  
  // Compute minor radius of torus
  double lrad = 0.0;
  double d2 = 0.0;
  int lnmb = 0;
  double beta = 0.5*M_PI + alpha;
  double phi = 0.5*M_PI - alpha;
  double fac = 1.0/sin(0.5*beta);
  double div = fac - 1.0;
  Point loc = rotational->location();
  double rd = (locp-loc)*axis;
  Point centr = loc + rd*axis;
  for (size_t kj=0; kj<blend_pts.size(); ++kj)
    {
      for (size_t ki=0; ki<blend_pts[kj].size(); ++ki)
	{
	  Vector3D xyz = blend_pts[kj][ki]->getPoint();
	  Point ptpos(xyz[0], xyz[1], xyz[2]);

	  // Localize with respect to the guiding curve(s)
	  double tpar, dist=std::numeric_limits<double>::max();
	  Point close;
	  int ix = -1;
	  for (size_t kr=0; kr<cvs.size(); ++kr)
	    {
	      double tpar0, dist0;
	      Point close0;
	      cvs[kr]->closestPoint(ptpos, cvs[kr]->startparam(), cvs[kr]->endparam(),
				tpar0, close0, dist0);
	      if (dist0 < dist)
		{
		  tpar = tpar0;
		  dist = dist0;
		  close = close0;
		  ix = (int)kr;
		}
	      if (ix < 0)
		continue;

	      // Define line in the point between the adjacent surfaces
	      vector<Point> der(2);
	      cvs[ix]->point(der, tpar, 1);
	      der[1].normalize();
	      Point dir1 = der[1].cross(normal);
	      Point vec1 = der[0] - centr;
	      if ((dir1*vec1 < 0.0 && (!rot_out)) || (dir1*vec1 > 0.0 && rot_out))
		dir1 *= -1;
	      Point dir2 = (plane_out) ? normal : -normal;
	      Point dir3;
	      if (alpha > 0)
		{
		  Matrix3D mat;
		  mat.setToRotation(-alpha, der[1][0], der[1][1], der[1][2]);
		  Vector3D dir2_2(dir2[0], dir2[1], dir2[2]);
		  Vector3D dir3_2 = mat*dir2_2;
		  dir3 = Point(dir3_2[0], dir3_2[1], dir3_2[2]);
		}
	      else
		dir3 = dir2;
	      Point dir4 = 0.5*(dir1 + dir3);
#ifdef DEBUG_BLEND
	      std::ofstream of("midlin.g2");
	      of << "410 1 0 4 0 255 0 255" << std::endl;
	      of << "1" << std::endl;
	      of << der[0] << " " << der[0]+dir1 << std::endl;
	      of << "410 1 0 4 255 0 0 255" << std::endl;
	      of << "1" << std::endl;
	      of << der[0] << " " << der[0]+dir4 << std::endl;
	      of << "410 1 0 4 0 0 255 255" << std::endl;
	      of << "1" << std::endl;
	      of << der[0] << " " << der[0]+dir2 << std::endl;
	      
#endif
	      shared_ptr<Line> line(new Line(der[0], dir4));
	      line->setParameterInterval(-2*width, 2*width);
	      double lpar, ldist;
	      Point lclose;
	      line->closestPoint(ptpos, line->startparam(), line->endparam(),
				 lpar, lclose, ldist);
	      if (ldist <= approx_tol_)
		{
		  double dd2 = close.dist(lclose);
		  double xlen = dd2/div;
		  lrad += xlen;
		  d2 += dd2;
		  lnmb++;
		}
	      int stop_break = 1;
	    }
	}
    }

  if (lnmb > 0)
    {
      lrad /= (double)lnmb;
      d2 /= (double)lnmb;
    }

  int sgn1 = (plane_out) ? -1 : 1;
  int sgn2 = rot_out ? -1 : 1;
  Point pos = centr - sgn1*lrad*normal;
  double Rrad = rotational->radius(0.0, rd);
  double sd = lrad/sin(phi);
  Point Cx = rotational->direction2();
  torus = shared_ptr<Torus>(new Torus(Rrad+sgn2*sd, lrad, pos, normal, Cx));
  if (rot_out)
    {
      RectDomain dom = torus->getParameterBounds();
      try {
	torus->setParameterBounds(dom.umin(), dom.vmin()-M_PI,
				  dom.umax(), dom.vmax()-M_PI);
      }
      catch (...)
	{
	}
    }
  
#ifdef DEBUG_BLEND
  std::ofstream of2("tor_blend.g2");
  torus->writeStandardHeader(of2);
  torus->write(of2);
  RectDomain dom2 = torus->getParameterBounds();
  vector<shared_ptr<ParamCurve> > tmp_cvs = torus->constParamCurves(dom2.vmin(), true);
  tmp_cvs[0]->writeStandardHeader(of2);
  tmp_cvs[0]->write(of2);
  double tf = (sgn2 == 1) ? 0.5 : 1.5;
  vector<shared_ptr<ParamCurve> > tmp_cvs2 =
    torus->constParamCurves(dom2.vmin()+tf*M_PI, true);
  tmp_cvs2[0]->writeStandardHeader(of2);
  tmp_cvs2[0]->write(of2);
#endif
  return torus;
}

//===========================================================================
double
RevEng::computeTorusRadius(vector<vector<RevEngPoint*> >& blend_pts,
			   shared_ptr<CurveOnSurface>& cv,
			   shared_ptr<ElementarySurface> elem1,
			   shared_ptr<ElementarySurface> elem2,
			   double width, bool out1, bool out2, int sgn,
			   double& d2)
//===========================================================================
{
  d2 = 0.0;
  double alpha1 = 0.0, alpha2 = 0.0;
  shared_ptr<Cone> cone1 = dynamic_pointer_cast<Cone,ElementarySurface>(elem1);
  if (cone1.get())
    alpha1 = cone1->getConeAngle();
  shared_ptr<Cone> cone2 = dynamic_pointer_cast<Cone,ElementarySurface>(elem2);
  if (cone2.get())
    alpha2 = cone2->getConeAngle();
  if (cone1.get() && cone2.get())
    return 0.0;   // Two cones are not handled
  double alpha = fabs(alpha1) + fabs(alpha2);
  shared_ptr<Plane> plane;
  if (elem1->instanceType() == Class_Plane)
    plane = dynamic_pointer_cast<Plane,ElementarySurface>(elem1);
  else if (elem2->instanceType() == Class_Plane)
    plane = dynamic_pointer_cast<Plane,ElementarySurface>(elem2);
  int state = (plane.get()) ? 1 : 2;
  shared_ptr<ElementarySurface> rotational;
  if (elem1->instanceType() == Class_Plane)
    rotational = elem2;
  else if (elem2->instanceType() == Class_Plane)
    rotational = elem1;
  else if (cone1.get())
    rotational = elem2;  // elem2 is expected to be a cylinder
  else if (cone2.get())
    rotational = elem1;
  Point axis = rotational->direction();
  Point normal = (plane.get()) ? plane->direction() : axis;
  if (state == 1)
    normal *= sgn;
  if ((state == 1 && elem2->instanceType() == Class_Plane) ||
      (state == 2 && elem2->instanceType() == Class_Cylinder))
    std::swap(out1,out2);  // Call by value means this swap stays local
  
  // Compute minor radius of torus
  double lrad = 0.0;
  int lnmb = 0;
  double beta = (plane.get()) ? 0.5*M_PI + alpha : M_PI-fabs(alpha);
  //double phi = beta - 2.0*alpha;
  double fac = 1.0/sin(0.5*beta);
  double div = fac - 1.0;
  Point loc = rotational->location();
  Point locp = (plane.get()) ? plane->location() :
    cv->ParamCurve::point(0.5*(cv->startparam()+cv->endparam()));;
  double rd = (locp-loc)*axis;
  Point centr = loc + rd*axis;
  for (size_t kj=0; kj<blend_pts.size(); ++kj)
    {
      for (size_t ki=0; ki<blend_pts[kj].size(); ++ki)
	{
	  Vector3D xyz = blend_pts[kj][ki]->getPoint();
	  Point ptpos(xyz[0], xyz[1], xyz[2]);

	  // Localize with respect to the guiding curve(s)
	  double tpar, dist;
	  Point close;
	  cv->closestPoint(ptpos, cv->startparam(), cv->endparam(),
			   tpar, close, dist);

	  // Define line in the point between the adjacent surfaces
	  vector<Point> der(2);
	  cv->point(der, tpar, 1);
	  der[1].normalize();
	  Point dir1 = der[1].cross(normal);
	  Point vec1 = der[0] - centr;
	  if ((dir1*vec1 < 0.0 && (!out2)) || (dir1*vec1 > 0.0 && out2))
	    dir1 *= -1;
	  Point dir2 = (out1) ? normal : -normal;
	  Point dir3;
	  if (alpha > 0)
	    {
	      Matrix3D mat;
	      if (state == 1)
		mat.setToRotation(-alpha, der[1][0], der[1][1], der[1][2]);
	      else
		mat.setToRotation(0.5*beta, der[1][0], der[1][1], der[1][2]);
	      Vector3D dir2_2(dir2[0], dir2[1], dir2[2]);
	      Vector3D dir3_2 = mat*dir2_2;
	      dir3 = Point(dir3_2[0], dir3_2[1], dir3_2[2]);
	    }
	  else
	    dir3 = dir2;
	  Point dir4 = (state == 1 )? 0.5*(dir1 + dir3) : dir3;
	  if (state == 2 && (!out2))
	    dir4 *= -1;
#ifdef DEBUG_BLEND
	  std::ofstream of("midlin.g2");
	  of << "410 1 0 4 0 255 0 255" << std::endl;
	  of << "1" << std::endl;
	  of << der[0] << " " << der[0]+dir1 << std::endl;
	  of << "410 1 0 4 255 0 0 255" << std::endl;
	  of << "1" << std::endl;
	  of << der[0] << " " << der[0]+dir4 << std::endl;
	  of << "410 1 0 4 0 0 255 255" << std::endl;
	  of << "1" << std::endl;
	  of << der[0] << " " << der[0]+dir2 << std::endl;
	  of << "410 1 0 4 100 100 55 255" << std::endl;
	  of << "1" << std::endl;
	  of << der[0] << " " << der[0]+dir3 << std::endl;
	      
#endif
	  shared_ptr<Line> line(new Line(der[0], dir4));
	  line->setParameterInterval(-2*width, 2*width);
	  double lpar, ldist;
	  Point lclose;
	  line->closestPoint(ptpos, line->startparam(), line->endparam(),
			     lpar, lclose, ldist);
	  if (ldist <= approx_tol_)
	    {
	      double dd2 = close.dist(lclose);
	      double xlen = dd2/div;
	      lrad += xlen;
	      d2 += dd2;
	      lnmb++;
	    }
	  int stop_break = 1;
	}
    }


  if (lnmb > 0)
    {
      lrad /= (double)lnmb;
     d2 /= (double)lnmb;
    }

  return lrad;
}

//===========================================================================
bool
RevEng::getTorusParameters(shared_ptr<ElementarySurface> elem1,
			   shared_ptr<ElementarySurface> elem2, Point pos,
			   double radius, double d2, bool out1, bool out2, int sgn,
			   double& Rrad, Point& centre, Point& normal, Point& Cx,
			   bool check_common)
//===========================================================================
{
  double alpha1 = 0.0, alpha2 = 0.0;
  shared_ptr<Cone> cone1 = dynamic_pointer_cast<Cone,ElementarySurface>(elem1);
  if (cone1.get())
    alpha1 = cone1->getConeAngle();
  shared_ptr<Cone> cone2 = dynamic_pointer_cast<Cone,ElementarySurface>(elem2);
  if (cone2.get())
    alpha2 = cone2->getConeAngle();
  if (cone1.get() && cone2.get())
    return false;   // Two cones are not handled
  double alpha = fabs(alpha1) + fabs(alpha2);
  shared_ptr<Plane> plane;
  if (elem1->instanceType() == Class_Plane)
    plane = dynamic_pointer_cast<Plane,ElementarySurface>(elem1);
  else if (elem2->instanceType() == Class_Plane)
    plane = dynamic_pointer_cast<Plane,ElementarySurface>(elem2);
  int state = (plane.get()) ? 1 : 2;
  shared_ptr<ElementarySurface> rotational;
  if (elem1->instanceType() == Class_Plane)
    rotational = elem2;
  else if (elem2->instanceType() == Class_Plane)
    rotational = elem1;
  else if (cone1.get())
    rotational = elem2;  // elem2 is expected to be a cylinder
  else if (cone2.get())
    rotational = elem1;
  Point axis = rotational->direction();
  normal = (plane.get()) ? plane->direction() : axis;
  if (state == 1)
    normal *= sgn;
  Point loc = rotational->location();
  if ((state == 1 && elem2->instanceType() == Class_Plane) ||
      (state == 2 && elem2->instanceType() == Class_Cylinder))
    std::swap(out1,out2);  // Call by value means this swap stays local
  
  double rd = (pos-loc)*axis;
  Point centr = loc + rd*axis;
  double beta = (plane.get()) ? 0.5*M_PI + alpha : M_PI-fabs(alpha);
  double phi = beta - 2.0*alpha;
  double phi2 = 0.5*(M_PI - beta);
  int sgn1 = (out1) ? -1 : 1;
  int sgn2 = out2 ? -1 : 1;
  if (state == 2)
    {
      sgn2 *= -1;
    }
  double hh = (state == 1) ? radius : sgn2*(radius + d2)*sin(phi2);
  centre = centr - sgn1*hh*normal;
  Rrad = rotational->radius(0.0, rd);
  double sd = (state == 1) ? radius/sin(phi) : (radius + d2)*cos(phi2);
  Cx = rotational->direction2();
  Rrad += (sgn2*sd);
  if (radius > Rrad && check_common)
    return false;
  
#ifdef DEBUG_BLEND
  shared_ptr<Torus> torus(new Torus(Rrad, radius, centre, normal, Cx));
  // if (sgn2 < 0)
  //   {
  //     RectDomain dom = torus->getParameterBounds();
  //     torus->setParameterBounds(dom.umin(), dom.vmin()-M_PI,
  // 				dom.umax(), dom.vmax()-M_PI);
  //   }
  std::ofstream of2("tor_blend.g2");
  torus->writeStandardHeader(of2);
  torus->write(of2);
  RectDomain dom2 = torus->getParameterBounds();
  vector<shared_ptr<ParamCurve> > tmp_cvs = torus->constParamCurves(dom2.vmin(), true);
  tmp_cvs[0]->writeStandardHeader(of2);
  tmp_cvs[0]->write(of2);
  double tf = (sgn2 == 1) ? 0.5 : 1.5;
  vector<shared_ptr<ParamCurve> > tmp_cvs2 =
    torus->constParamCurves(dom2.vmin()+tf*M_PI, true);
  tmp_cvs2[0]->writeStandardHeader(of2);
  tmp_cvs2[0]->write(of2);
#endif
  return true;
}

//===========================================================================
shared_ptr<Torus>
RevEng::torusBlend(vector<vector<RevEngPoint*> >& blend_pts,
		   shared_ptr<CurveOnSurface>& cv,
		   shared_ptr<ElementarySurface> elem1,
		   shared_ptr<ElementarySurface> elem2,
		   double width, bool out1, bool out2, int sgn)
//===========================================================================
{
  shared_ptr<Torus> torus;

  double alpha1 = 0.0, alpha2 = 0.0;
  shared_ptr<Cone> cone1 = dynamic_pointer_cast<Cone,ElementarySurface>(elem1);
  if (cone1.get())
    alpha1 = cone1->getConeAngle();
  shared_ptr<Cone> cone2 = dynamic_pointer_cast<Cone,ElementarySurface>(elem2);
  if (cone2.get())
    alpha2 = cone2->getConeAngle();
  if (cone1.get() && cone2.get())
    return torus;   // Two cones are not handled
  double alpha = fabs(alpha1) + fabs(alpha2);
  shared_ptr<Plane> plane;
  if (elem1->instanceType() == Class_Plane)
    plane = dynamic_pointer_cast<Plane,ElementarySurface>(elem1);
  else if (elem2->instanceType() == Class_Plane)
    plane = dynamic_pointer_cast<Plane,ElementarySurface>(elem2);
  int state = (plane.get()) ? 1 : 2;
  shared_ptr<ElementarySurface> rotational;
  if (elem1->instanceType() == Class_Plane)
    rotational = elem2;
  else if (elem2->instanceType() == Class_Plane)
    rotational = elem1;
  else if (cone1.get())
    rotational = elem2;  // elem2 is expected to be a cylinder
  else if (cone2.get())
    rotational = elem1;
  Point axis = rotational->direction();
  Point normal = (plane.get()) ? plane->direction() : axis;
  if (state == 1)
    normal *= sgn;
  if ((state == 1 && elem2->instanceType() == Class_Plane) ||
      (state == 2 && elem2->instanceType() == Class_Cylinder))
    std::swap(out1,out2);  // Call by value means this swap stays local
  
  // Compute minor radius of torus
  double lrad = 0.0;
  double d2 = 0.0;
  int lnmb = 0;
  double beta = (plane.get()) ? 0.5*M_PI + alpha : M_PI-fabs(alpha);
  double phi = beta - 2.0*alpha;
  double fac = 1.0/sin(0.5*beta);
  double div = fac - 1.0;
  Point loc = rotational->location();
  Point locp = (plane.get()) ? plane->location() :
    cv->ParamCurve::point(0.5*(cv->startparam()+cv->endparam()));;
  double rd = (locp-loc)*axis;
  Point centr = loc + rd*axis;
  for (size_t kj=0; kj<blend_pts.size(); ++kj)
    {
      for (size_t ki=0; ki<blend_pts[kj].size(); ++ki)
	{
	  Vector3D xyz = blend_pts[kj][ki]->getPoint();
	  Point ptpos(xyz[0], xyz[1], xyz[2]);

	  // Localize with respect to the guiding curve(s)
	  double tpar, dist;
	  Point close;
	  cv->closestPoint(ptpos, cv->startparam(), cv->endparam(),
			   tpar, close, dist);

	  // Define line in the point between the adjacent surfaces
	  vector<Point> der(2);
	  cv->point(der, tpar, 1);
	  der[1].normalize();
	  Point dir1 = der[1].cross(normal);
	  Point vec1 = der[0] - centr;
	  if ((dir1*vec1 < 0.0 && (!out2)) || (dir1*vec1 > 0.0 && out2))
	    dir1 *= -1;
	  Point dir2 = (out1) ? normal : -normal;
	  Point dir3;
	  if (alpha > 0)
	    {
	      Matrix3D mat;
	      if (state == 1)
		mat.setToRotation(-alpha, der[1][0], der[1][1], der[1][2]);
	      else
		mat.setToRotation(0.5*beta, der[1][0], der[1][1], der[1][2]);
	      Vector3D dir2_2(dir2[0], dir2[1], dir2[2]);
	      Vector3D dir3_2 = mat*dir2_2;
	      dir3 = Point(dir3_2[0], dir3_2[1], dir3_2[2]);
	    }
	  else
	    dir3 = dir2;
	  Point dir4 = (state == 1 )? 0.5*(dir1 + dir3) : dir3;
	  if (state == 2 && (!out2))
	    dir4 *= -1;
#ifdef DEBUG_BLEND
	  std::ofstream of("midlin.g2");
	  of << "410 1 0 4 0 255 0 255" << std::endl;
	  of << "1" << std::endl;
	  of << der[0] << " " << der[0]+dir1 << std::endl;
	  of << "410 1 0 4 255 0 0 255" << std::endl;
	  of << "1" << std::endl;
	  of << der[0] << " " << der[0]+dir4 << std::endl;
	  of << "410 1 0 4 0 0 255 255" << std::endl;
	  of << "1" << std::endl;
	  of << der[0] << " " << der[0]+dir2 << std::endl;
	  of << "410 1 0 4 100 100 55 255" << std::endl;
	  of << "1" << std::endl;
	  of << der[0] << " " << der[0]+dir3 << std::endl;
	      
#endif
	  shared_ptr<Line> line(new Line(der[0], dir4));
	  line->setParameterInterval(-2*width, 2*width);
	  double lpar, ldist;
	  Point lclose;
	  line->closestPoint(ptpos, line->startparam(), line->endparam(),
			     lpar, lclose, ldist);
	  if (ldist <= approx_tol_)
	    {
	      double dd2 = close.dist(lclose);
	      double xlen = dd2/div;
	      lrad += xlen;
	      d2 += dd2;
	      lnmb++;
	    }
	  int stop_break = 1;
	}
    }


  if (lnmb > 0)
    {
      lrad /= (double)lnmb;
      d2 /= (double)lnmb;
    }

  int sgn1 = (out1) ? -1 : 1;
  int sgn2 = out2 ? -1 : 1;
  if (state == 2)
    {
      sgn2 *= -1;
    }
  double phi2 = 0.5*(M_PI - beta);
  double hh = (state == 1) ? lrad : sgn2*(lrad + d2)*sin(phi2);
  Point pos = centr - sgn1*hh*normal;
  double Rrad = rotational->radius(0.0, rd);
  double sd = (state == 1) ? lrad/sin(phi) : (lrad + d2)*cos(phi2);
  Point Cx = rotational->direction2();
  torus = shared_ptr<Torus>(new Torus(Rrad+sgn2*sd, lrad, pos, normal, Cx));
  if (sgn2 > 0)
    {
      RectDomain dom = torus->getParameterBounds();
      try {
	torus->setParameterBounds(dom.umin(), dom.vmin()-M_PI,
				  dom.umax(), dom.vmax()-M_PI);
      }
      catch (...)
	{
	}
    }
  
#ifdef DEBUG_BLEND
  std::ofstream of2("tor_blend.g2");
  torus->writeStandardHeader(of2);
  torus->write(of2);
  RectDomain dom2 = torus->getParameterBounds();
  vector<shared_ptr<ParamCurve> > tmp_cvs = torus->constParamCurves(dom2.vmin(), true);
  tmp_cvs[0]->writeStandardHeader(of2);
  tmp_cvs[0]->write(of2);
  double tf = (sgn2 == 1) ? 0.5 : 1.5;
  vector<shared_ptr<ParamCurve> > tmp_cvs2 =
    torus->constParamCurves(dom2.vmin()+tf*M_PI, true);
  tmp_cvs2[0]->writeStandardHeader(of2);
  tmp_cvs2[0]->write(of2);
#endif
  return torus;
}

//===========================================================================
double
RevEng::computeCylinderRadius(vector<vector<RevEngPoint*> > blend_pts,
			    double width, const Point& pos, const Point& axis,
			    const Point& dir1, const Point& dir2)
//===========================================================================
{
  double eps = 1.0e-6;
  double upar2, vpar2, dist2;
  Point close, surfnorm, close2;
  double alpha = dir1.angle(dir2);
  Point linedir = width*dir1 + width*dir2;
  Point planenorm = linedir.cross(axis);
  planenorm.normalize();
  shared_ptr<Plane> plane(new Plane(pos, planenorm, linedir));
  double lrad = 0.0;
  int lnmb = 0;
  double fac = 1.0/sin(0.5*alpha);
  double div = fac - 1.0;
  for (size_t kj=0; kj<blend_pts.size(); ++kj)
    {
      for (size_t ki=0; ki<blend_pts[kj].size(); ++ki)
	{
	  Vector3D xyz = blend_pts[kj][ki]->getPoint();
	  Point ptpos(xyz[0], xyz[1], xyz[2]);

	  plane->closestPoint(ptpos, upar2, vpar2, close2, dist2, eps);
	  if (dist2 <= approx_tol_)
	    {
	      Point pos2 = plane->ParamSurface::point(0.0, vpar2);
	      double dd2 = pos2.dist(close2);
	      double xlen = dd2/div;
	      lrad += xlen;
	      lnmb++;
	    }
	}
    }
  if (lnmb > 0)
    lrad /= (double)lnmb;
  return lrad;
}


//===========================================================================
shared_ptr<Cylinder>
RevEng::createCylinderBlend(vector<vector<RevEngPoint*> > blend_pts,
			    double rad1, const Point& pos, const Point& axis,
			    const Point& dir1, const Point& dir2, int sgn)
//===========================================================================
{
  double eps = 1.0e-6;
  double upar2, vpar2, dist2;
  Point close, surfnorm, close2;
  double alpha = dir1.angle(dir2);
  Point linedir = rad1*dir1 + rad1*dir2;
  Point planenorm = linedir.cross(axis);
  planenorm.normalize();
  shared_ptr<Plane> plane(new Plane(pos, planenorm, linedir));
  double lrad = 0.0;
  int lnmb = 0;
  double fac = 1/sin(0.5*alpha);
  double div = fac - 1.0;
  for (size_t kj=0; kj<blend_pts.size(); ++kj)
    {
      for (size_t ki=0; ki<blend_pts[kj].size(); ++ki)
	{
	  Vector3D xyz = blend_pts[kj][ki]->getPoint();
	  Point ptpos(xyz[0], xyz[1], xyz[2]);

	  plane->closestPoint(ptpos, upar2, vpar2, close2, dist2, eps);
	  if (dist2 <= approx_tol_)
	    {
	      Point pos2 = plane->ParamSurface::point(0.0, vpar2);
	      double dd2 = pos2.dist(close2);
	      double xlen = dd2/div;
	      lrad += xlen;
	      lnmb++;
	    }
	}
    }
  if (lnmb > 0)
    lrad /= (double)lnmb;
  
  Point Cx = sgn*(dir1 + dir2);
  Point Cy = axis.cross(Cx);
  // Point pos2;
  // double rad2;
  Point low = bbox_.low();
  Point high = bbox_.high();

  Point vec = dir1 + dir2;
  vec.normalize();
  Point pos3 = pos + sgn*fac*lrad*vec;
  shared_ptr<Cylinder> cyl(new Cylinder(lrad, pos3, axis, Cx));
  double diag = low.dist(high);
  cyl->setParamBoundsV(-0.5*diag,0.5*diag);
  
#ifdef DEBUG_BLEND
  std::ofstream of("blend_cyl2.g2");
  cyl->writeStandardHeader(of);
  cyl->write(of);
#endif

  
  return cyl;
}



//===========================================================================
void RevEng::surfaceExtractOutput(int idx,
				  vector<vector<RevEngPoint*> > out_groups,
				  vector<HedgeSurface*> prev_surfs)
//===========================================================================
{
  for (size_t kr=0; kr<out_groups.size(); ++kr)
    {
      for (size_t kh=0; kh<out_groups[kr].size(); ++kh)
	out_groups[kr][kh]->unsetRegion();  
    }
  
  int classtype = (idx >= 0) ? regions_[idx]->getClassification() :
    CLASSIFICATION_UNDEF;
  for (size_t kr=0; kr<out_groups.size(); ++kr)
    {
      shared_ptr<RevEngRegion> reg(new RevEngRegion(classtype,
						    edge_class_type_,
						    out_groups[kr]));
      if (idx >= 0)
	reg->setPreviousReg(regions_[idx].get());
      reg->setRegionAdjacency();
#ifdef DEBUG_CHECK
      bool connect = reg->isConnected();
      connect = reg->isConnected();
      if (!connect)
	{
	  std::cout << "surfaceExtractOutput, disconnected region " << idx << std::endl;
	}
#endif
      bool integrate = reg->integrateInAdjacent(mean_edge_len_,
						min_next_, max_next_,
						approx_tol_, 0.5,
						max_nmb_outlier_,
						(idx >= 0) ? regions_[idx].get() : 0);
      if (!integrate)
	regions_.push_back(reg);
    }	  
  for (size_t kr=0; kr<prev_surfs.size(); ++kr)
    {
      size_t kj;
      for (kj=0; kj<surfaces_.size(); ++kj)
	if (surfaces_[kj].get() == prev_surfs[kr])
	  break;
      if (kj < surfaces_.size())
	surfaces_.erase(surfaces_.begin()+kj);
    }
}

//===========================================================================
bool RevEng::segmentByPlaneGrow(int ix, int min_point_in, double angtol)
//===========================================================================
{
#ifdef DEBUG
  std::ofstream ofreg("segment_reg.g2");
  regions_[ix]->writeRegionInfo(ofreg);
#endif

  vector<shared_ptr<HedgeSurface> > plane_sfs;
  vector<HedgeSurface*> prev_surfs;
  vector<vector<RevEngPoint*> > out_groups;
  regions_[ix]->segmentByPlaneGrow(mainaxis_, approx_tol_, angtol, min_point_in, 
				   plane_sfs, prev_surfs, out_groups);
#ifdef DEBUG
  for (size_t ki=0; ki<plane_sfs.size(); ++ki)
    {
      plane_sfs[ki]->surface()->writeStandardHeader(ofreg);
      plane_sfs[ki]->surface()->write(ofreg);
    }
#endif
  
  if (out_groups.size() > 0 || prev_surfs.size() > 0)
    surfaceExtractOutput(ix, out_groups, prev_surfs);
  if (plane_sfs.size() > 0)
    surfaces_.insert(surfaces_.end(), plane_sfs.begin(), plane_sfs.end());

  bool segmented = (out_groups.size() > 0);
  return segmented;
}

//===========================================================================
bool RevEng::segmentByAxis(int ix, int min_point_in)
//===========================================================================
{
  double angtol = 5.0*anglim_;
  vector<shared_ptr<HedgeSurface> > hedgesfs;
  vector<shared_ptr<RevEngRegion> > added_reg;
  vector<vector<RevEngPoint*> > separate_groups;
  vector<RevEngPoint*> single_points;
  bool segmented = regions_[ix]->extractCylByAxis(mainaxis_, min_point_in,
						  min_point_region_,
					       approx_tol_, angtol,
					       prefer_elementary_,
					       hedgesfs, added_reg,
					       separate_groups,
					       single_points);
  if (segmented && single_points.size() > 0)
    single_points_.insert(single_points_.end(), single_points.begin(),
			  single_points.end());
#ifdef DEBUG
  if (segmented)
    {
      std::ofstream of("seg_by_axis.g2");
      int num = regions_[ix]->numPoints();
      of << "400 1 0 4 255 0 0 255" << std::endl;
      of <<  num << std::endl;
      for (int ka=0; ka<num; ++ka)
	of << regions_[ix]->getPoint(ka)->getPoint() << std::endl;

      for (size_t ki=0; ki<separate_groups.size(); ++ki)
	{
	  of << "400 1 0 4 0 255 0 255" << std::endl;
	  of <<  separate_groups[ki].size() << std::endl;
	  for (int ka=0; ka<(int)separate_groups[ki].size(); ++ka)
	    of << separate_groups[ki][ka]->getPoint() << std::endl;
	}
    }
#endif
  
  if (separate_groups.size() > 0)
    {
      vector<HedgeSurface*> prev_surfs;
      surfaceExtractOutput(ix, separate_groups, prev_surfs);
    }
  
  if (added_reg.size() > 0)
    regions_.insert(regions_.end(), added_reg.begin(), added_reg.end());
  if (hedgesfs.size() > 0)
    surfaces_.insert(surfaces_.end(), hedgesfs.begin(), hedgesfs.end());
  
  return segmented;
}

//===========================================================================
bool RevEng::segmentByContext(int ix, int min_point_in, double angtol, bool first)
//===========================================================================
{
#ifdef DEBUG_DIV
  vector<RevEngPoint*> branchpt = regions_[ix]->extractBranchPoints();
  if (branchpt.size() > 0)
    {
      std::ofstream ofb("branch_pts_seg.g2");
      ofb << "400 1 0 4 0 0 0 255" << std::endl;
      ofb << branchpt.size() << std::endl;
      for (size_t ki=0; ki<branchpt.size(); ++ki)
	ofb << branchpt[ki]->getPoint() << std::endl;
    }
  
#endif
  
  vector<vector<RevEngPoint*> > separate_groups;
  vector<RevEngRegion*> adj_planar = regions_[ix]->fetchAdjacentPlanar();
  vector<shared_ptr<HedgeSurface> > hedgesfs;
  vector<shared_ptr<RevEngRegion> > added_reg;
  vector<HedgeSurface*> prevsfs;
  bool segmented = false;
  if (adj_planar.size() > 1)
    {
      segmented = regions_[ix]->segmentByPlaneAxis(mainaxis_, min_point_in,
						   min_point_region_,
						   approx_tol_, angtol,
						   prefer_elementary_,
						   adj_planar, hedgesfs,
						   added_reg,
						   prevsfs, separate_groups);
    }

  if (!segmented)
    {
      // Extend with cylindrical
      vector<RevEngRegion*> adj_cyl = regions_[ix]->fetchAdjacentCylindrical();
      if (adj_cyl.size() > 0)
	adj_planar.insert(adj_planar.end(), adj_cyl.begin(), adj_cyl.end());
      if (adj_planar.size() > 0)
	segmented =
	  regions_[ix]->segmentByAdjSfContext(mainaxis_, min_point_in,
					      min_point_region_,
					       approx_tol_, angtol,
					       adj_planar, separate_groups);
    }
  
  if (!segmented)
    {
      // Search for context direction
      double angtol2 = 2.0*angtol;

      Point direction = regions_[ix]->directionFromAdjacent(angtol);
      vector<vector<RevEngPoint*> > separate_groups2;
      if (direction.dimension() == 3)
	segmented = regions_[ix]->segmentByDirectionContext(min_point_in, approx_tol_,
							    direction, angtol2,
							    separate_groups2);
      if (separate_groups2.size() > 0)
	separate_groups.insert(separate_groups.end(), separate_groups2.begin(),
			       separate_groups2.end());
      if (segmented && (!regions_[ix]->hasSurface()))
	{
	  double angtol = -1.0;
	  int pass = 1;
	  bool found = recognizeOneSurface(ix, min_point_in, angtol, pass);
	  int stop_break = 1;
	}
    }
  
#ifdef DEBUG
  if (segmented)
    {
      std::ofstream of("seg_by_context.g2");
      int num = regions_[ix]->numPoints();
      of << "400 1 0 4 255 0 0 255" << std::endl;
      of <<  num << std::endl;
      for (int ka=0; ka<num; ++ka)
	of << regions_[ix]->getPoint(ka)->getPoint() << std::endl;

      for (size_t ki=0; ki<separate_groups.size(); ++ki)
	{
	  of << "400 1 0 4 0 255 0 255" << std::endl;
	  of <<  separate_groups[ki].size() << std::endl;
	  for (int ka=0; ka<(int)separate_groups[ki].size(); ++ka)
	    of << separate_groups[ki][ka]->getPoint() << std::endl;
	}
    }
#endif
  
  if (separate_groups.size() > 0)
    {
      vector<HedgeSurface*> prev_surfs;
      surfaceExtractOutput(ix, separate_groups, prev_surfs);
    }

  if (added_reg.size() > 0)
    regions_.insert(regions_.end(), added_reg.begin(), added_reg.end());
  if (hedgesfs.size() > 0)
    surfaces_.insert(surfaces_.end(), hedgesfs.begin(), hedgesfs.end());
  
  return segmented;
}

//===========================================================================
void RevEng::growSurface(int& ix, int pass)
//===========================================================================
{
  vector<RevEngRegion*> grown_regions;
  int min_nmb = 5*min_point_region_;  // Should be set from distribution of how many
	  // points the regions have
  double angtol = 5.0*anglim_;
  vector<HedgeSurface*> adj_surfs;
  vector<RevEngEdge*> adj_edgs;
  regions_[ix]->growWithSurf(mainaxis_, min_point_region_,
			     approx_tol_, angtol, grown_regions,
			     adj_surfs, adj_edgs, (pass>1));
  updateRegionsAndSurfaces(ix, grown_regions, adj_surfs);
  for (size_t ki=0; ki<adj_edgs.size(); ++ki)
    {
      size_t kj;
      for (kj=0; kj<edges_.size(); ++kj)
	if (edges_[kj].get() == adj_edgs[ki])
	  break;
      if (kj < edges_.size())
	edges_.erase(edges_.begin()+kj);
    }
}

//===========================================================================
void RevEng::growBlendSurface(int& ix)
//===========================================================================
{
  // Collect associated blend surfaces
  RevEngEdge *edge = regions_[ix]->getBlendEdge();
  if (!edge)
    return;

#ifdef DEBUG_BLEND
  std::ofstream of("blend_grow.g2");
  regions_[ix]->writeRegionPoints(of);
#endif

  vector<RevEngRegion*> next_blend;
  RevEngRegion *adj[2];
  edge->getAdjacent(adj[0], adj[1]);
  for (int ka=0; ka<2; ++ka)
    {
      if (!adj[ka])
	continue;
      vector<RevEngEdge*> rev_edgs = adj[ka]->getAllRevEdges();
      for (size_t ki=0; ki<rev_edgs.size(); ++ki)
	{
	  RevEngRegion *blendreg = rev_edgs[ki]->getBlendRegSurf();
	  if (blendreg)
	    {
	      if (blendreg == regions_[ix].get())
		continue;
	      size_t kj=0;
	      for (kj=0; kj<next_blend.size(); ++kj)
		if (next_blend[kj] == blendreg)
		  break;
	      if (kj < next_blend.size())
		continue;
	      next_blend.push_back(blendreg);
	    }
	}
    }

#ifdef DEBUG_BLEND
  for (size_t kr=0; kr<next_blend.size(); ++kr)
    next_blend[kr]->writeRegionPoints(of);
#endif
  
  vector<RevEngRegion*> grown_regions;
  double angtol = 5.0*anglim_;
  vector<HedgeSurface*> adj_surfs;
  vector<vector<RevEngPoint*> > added_regs;
  regions_[ix]->growBlendSurf(next_blend, approx_tol_, angtol,
			      grown_regions, added_regs);

  updateRegionsAndSurfaces(ix, grown_regions, adj_surfs);
  if (added_regs.size() > 0)
    {
      vector<HedgeSurface*> prev_surfs;
      surfaceExtractOutput(ix, added_regs, prev_surfs);
    }
}

//===========================================================================
void RevEng::growMasterSurface(int& ix)
//===========================================================================
{

#ifdef DEBUG_BLEND
  std::ofstream of("master_grow.g2");
  regions_[ix]->writeRegionPoints(of);
#endif

  vector<RevEngRegion*> grown_regions;
  double angtol = 5.0*anglim_;
  vector<HedgeSurface*> adj_surfs;
  int small_lim = min_point_region_/20;
  regions_[ix]->joinToCurrent(approx_tol_, angtol, small_lim,
			      grown_regions);

  updateRegionsAndSurfaces(ix, grown_regions, adj_surfs);
}


//===========================================================================
void RevEng::updateRegionsAndSurfaces(int& ix, vector<RevEngRegion*>& grown_regions,
				      vector<HedgeSurface*>& adj_surfs)
//===========================================================================
{
  if (grown_regions.size() > 0)
    {
      for (size_t kr=0; kr<grown_regions.size(); ++kr)
	{
// #ifdef DEBUG_CHECK
// 	  bool connect = grown_regions[kr]->isConnected();
// 	  connect = grown_regions[kr]->isConnected();
// 	  if (!connect)
// 	    std::cout << "updateRegionsAndSurfaces, disconnected region " << ix << kr << std::endl;
// #endif
	  size_t kj;
	  for (kj=0; kj<regions_.size(); )
	    {
	      if ((int)kj == ix)
		{
		  ++kj;
		  continue;
		}

	      if (grown_regions[kr] == regions_[kj].get())
		{
		  regions_.erase(regions_.begin()+kj);
		  if ((int)kj < ix)
		    --ix;
		}
	      else
		++kj;
	    }
	}
      for (size_t kr=0; kr<adj_surfs.size(); ++kr)
	{
	  size_t kj;
	  for (kj=0; kj<surfaces_.size(); ++kj)
	    if (surfaces_[kj].get() == adj_surfs[kr])
	      break;
	  if (kj < surfaces_.size())
	    surfaces_.erase(surfaces_.begin()+kj);
	}
    }

  for (size_t kj=0; kj<regions_.size(); ++kj)
    regions_[kj]->setVisited(false);

// #ifdef DEBUG_DIV
//   std::ofstream ofpts("sfpoints2.g2");
//   std::ofstream ofpar("sfparpoints2.txt");
//   int numpt = regions_[ix]->numPoints();
//   ofpts << "400 1 0 0" << std::endl;
//   ofpts << numpt << std::endl;
//   for (int ka=0; ka<numpt; ++ka)
//     {
//       RevEngPoint *pt = regions_[ix]->getPoint(ka);
//       ofpar << pt->getPar() << " " << pt->getPoint() << std::endl;
//       ofpts << pt->getPoint() << std::endl;
//     }
// #endif
  int stop_break = 1;
  
}


//===========================================================================
void RevEng::smallRegionSurfaces()
//===========================================================================
{
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

  // Add possible missing edges
  recognizeEdges();
  
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

   
#ifdef DEBUG
   std::cout << "Extend blend region collection" << std::endl;
#endif
   for (size_t ki=0; ki<edges_.size(); ++ki)
     {
       extendBlendAssociation(ki);
     }


  int nmb_sfs = (int)surfaces_.size();
  defineSmallRegionSurfaces();

#ifdef DEBUG
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      RevEngRegion *first = regions_[ki]->getPoint(0)->region();
      int num = regions_[ki]->numPoints();
      for (int ka=1; ka<num; ++ka)
	if (regions_[ki]->getPoint(ka)->region() != first)
	  std::cout << "Inconsistent region pointers, post defineSmallRegionSurfaces: " << ki << " " << ka << std::endl;
    }
#endif

#ifdef DEBUG
   checkConsistence("Regions10_1");

   if (regions_.size() > 0)
    {
      std::cout << "Regions10_1" << std::endl;
      std::ofstream of("regions10_1.g2");
      std::ofstream ofm("mid_regions10_1.g2");
      std::ofstream ofs("small_regions10_1.g2");
      writeRegionStage(of, ofm, ofs);

      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf10_1.g2");
	  writeRegionWithSurf(of);
	}
    }
   if (edges_.size() > 0)
     {
       std::ofstream ofe("edges10_1.g2");
       writeEdgeStage(ofe);
     }

   if (surfaces_.size() > 0)
     {
       std::ofstream ofsf("surfaces10_1.g2");
       for (size_t kr=0; kr<surfaces_.size(); ++kr)
	 {
	   RevEngRegion *reg = surfaces_[kr]->getRegion(0);
	   reg->writeSurface(ofsf);
	   //reg->writeRegionPoints(ofsf);
	 }
     }
#endif

  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

  // Check if small regions surfaces really belongs to a blend
  for (size_t ki=nmb_sfs; ki<surfaces_.size();)
    {
      RevEngRegion *reg = surfaces_[ki]->getRegion(0);
      vector<HedgeSurface*> removed_sfs;
      reg->checkEdgeAssociation(approx_tol_, min_point_region_, removed_sfs);
      if (removed_sfs.size() > 0)
	{
	  for (size_t kr=0; kr<removed_sfs.size(); ++kr)
	    {
	      size_t kj;
	      for (kj=0; kj<surfaces_.size(); ++kj)
		if (surfaces_[kj].get() == removed_sfs[kr])
		  break;
	      if (kj < surfaces_.size())
		surfaces_.erase(surfaces_.begin()+kj);
	    }
	}
      else
	++ki;
    }
  
#ifdef DEBUG
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      int num = regions_[ki]->numPoints();
      for (int ka=0; ka<num; ++ka)
	if (regions_[ki]->getPoint(ka)->region() != regions_[ki].get())
	  std::cout << "Inconsistent region pointers, post define small region surfaces: " << ki << " " << ka << std::endl;
    }
#endif
   // if ((int)surfaces_.size() > nmb_sfs)
   //   recognizeEdges();
   
#ifdef DEBUG
   checkConsistence("Regions10_1_2");

   if (regions_.size() > 0)
    {
      std::cout << "Regions10_1_2" << std::endl;
      std::ofstream of("regions10_1_2.g2");
      std::ofstream ofm("mid_regions10_1_2.g2");
      std::ofstream ofs("small_regions10_1_2.g2");
      writeRegionStage(of, ofm, ofs);

      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf10_1.g2");
	  writeRegionWithSurf(of);
	}
    }
   if (edges_.size() > 0)
     {
       std::ofstream ofe("edges10_1_2.g2");
       writeEdgeStage(ofe);
     }
#endif
  
   double angtol = 5.0*anglim_;
   for (int ka=nmb_sfs; ka<(int)surfaces_.size(); ++ka)
     {
       growSmallRegionSurface(ka);
     }
   //#if 0
   // Dismiss too small surfaces. First check connectivity
   int min_sf_pts = min_point_region_/5;
   for (int ka=0; ka<(int)regions_.size(); ++ka)
     {
       if ((!regions_[ka]->hasSurface()) || regions_[ka]->hasBlendEdge())
	 continue;
       vector<vector<RevEngPoint*> > separate_groups;
       vector<HedgeSurface*> out_sfs;
       vector<RevEngEdge*> out_edgs;
       regions_[ka]->splitRegion(separate_groups);
       if (separate_groups.size() > 0)
	 regions_[ka]->updateInfo(approx_tol_, angtol);
       if (regions_[ka]->numPoints() < min_sf_pts)
	 {
	   out_sfs.push_back(regions_[ka]->getSurface(0));
	   vector<RevEngEdge*> rev_edgs = regions_[ka]->getAllRevEdges();
	   for (size_t kr=0; kr<rev_edgs.size(); ++kr)
	     {
	       RevEngRegion *adj1, *adj2;
	       rev_edgs[kr]->getAdjacent(adj1, adj2);
	       RevEngRegion *other = (adj1 == regions_[ka].get()) ? adj2 : adj1;
	       other->removeRevEngEdge(rev_edgs[kr]);
	     }
	   out_edgs.insert(out_edgs.end(), rev_edgs.begin(), rev_edgs.end());
	   regions_[ka]->clearSurface();
	 }

       if (separate_groups.size() > 0 || out_sfs.size() > 0)
	   surfaceExtractOutput(ka, separate_groups, out_sfs);
       
       for (size_t ki=0; ki<out_edgs.size(); ++ki)
	 {
	   size_t kj;
	   for (kj=0; kj<edges_.size(); ++kj)
	     if (edges_[kj].get() == out_edgs[ki])
	       break;
	   if (kj < edges_.size())
	     edges_.erase(edges_.begin()+kj);
	 }
 
     }
   //#endif 
   recognizeEdges();
 
   

#ifdef DEBUG
   checkConsistence("Regions10_2");

   if (regions_.size() > 0)
    {
      std::cout << "Regions10_2" << std::endl;
      std::ofstream of("regions10_2.g2");
      std::ofstream ofm("mid_regions10_2.g2");
      std::ofstream ofs("small_regions10_2.g2");
      writeRegionStage(of, ofm, ofs);

      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf10_2.g2");
	  writeRegionWithSurf(of);
	}
    }
   if (edges_.size() > 0)
     {
       std::ofstream ofe("edges10_2.g2");
       writeEdgeStage(ofe);
     }
   if (surfaces_.size() > 0)
     {
       std::ofstream ofsf("surfaces10_2.g2");
       for (size_t kr=0; kr<surfaces_.size(); ++kr)
	 {
	   RevEngRegion *reg = surfaces_[kr]->getRegion(0);
	   reg->writeSurface(ofsf);
	   //reg->writeRegionPoints(ofsf);
	 }
     }
#ifdef DEBUG
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      int num = regions_[ki]->numPoints();
      for (int ka=0; ka<num; ++ka)
	if (regions_[ki]->getPoint(ka)->region() != regions_[ki].get())
	  std::cout << "Inconsistent region pointers, post grow small regions surfaces: " << ki << " " << ka << std::endl;
    }
#endif
   std::cout << "Define small region surfaces, regions: " << regions_.size() << ", surfaces: " << surfaces_.size() << std::endl;
#endif
   int stop_break = 1;
}

//===========================================================================
// Service function for growSmallRegionSurface
int selectBestAccuracy(double tol, int num, double maxd[2], double avd[2],
		       int num_in[2], int num2_in[2])
{
  int ix = (avd[0] < avd[1]) ? 0 : 1;
  int ix2 = 1 - ix;
  double fac = 1.1;
  if ((avd[ix2] < fac*avd[ix] && num_in[ix2] > num_in[ix] && num2_in[ix2] > num_in[ix]))
    std::swap(ix,ix2);

  return ix;
}

//===========================================================================
void RevEng::growSmallRegionSurface(int& ix)
//===========================================================================
{
  double angtol = 5.0*anglim_;
  
  RevEngRegion *reg = surfaces_[ix]->getRegion(0);
  shared_ptr<ParamSurface> surf = surfaces_[ix]->surface();
  bool cyllike = (surf->instanceType() == Class_Cylinder ||
		  surf->instanceType() == Class_Cone);
  int num_group_points = reg->numPoints();

  // Collect feasible adjacent regions
  vector<RevEngRegion*> adj_surf, adj_nosurf;
  reg->getAdjCandMerge(adj_surf, adj_nosurf);

#ifdef DEBUG_SMALL
  std::ofstream of1("grow_small.g2");
  reg->writeRegionPoints(of1);
  reg->writeSurface(of1);

  for (size_t kr=0; kr<adj_surf.size(); ++kr)
    {
      adj_surf[kr]->writeRegionPoints(of1);
      adj_surf[kr]->writeSurface(of1);
    }
  
  for (size_t kr=0; kr<adj_nosurf.size(); ++kr)
    {
      adj_nosurf[kr]->writeRegionPoints(of1);
    }      
#endif

  // Decide on master regions with surface
  vector<RevEngRegion*> to_remove;
  for (size_t ki=0; ki<adj_surf.size(); ++ki)
    {
      // Check combined accuracy based on current region and adjacent region
      int numpt = reg->numPoints() + adj_surf[ki]->numPoints();
      double maxd[2], avd[2], maxd_loc[2], avd_loc[2];
      int num_in[2], num2_in[2], num_in_loc[2], num2_in_loc[2];
      vector<vector<double> > parvals(2);   // Only points from the other region
      vector<vector<pair<double,double> > > distang(2);
      int sf_flag_loc[2];
      sf_flag_loc[0] =
	reg->getGrowAccuracy(adj_surf[ki], approx_tol_, angtol, maxd[0], avd[0],
			     num_in[0], num2_in[0], maxd_loc[0], avd_loc[0],
			     num_in_loc[0], num2_in_loc[0], parvals[0], distang[0]);
      sf_flag_loc[1] = 
	adj_surf[ki]->getGrowAccuracy(reg, approx_tol_, angtol, maxd[1], avd[1], 
				      num_in[1], num2_in[1], maxd_loc[1], avd_loc[1],
				      num_in_loc[1], num2_in_loc[1], parvals[1], distang[1]);

      int m_ix = selectBestAccuracy(approx_tol_, numpt, maxd, avd, num_in, num2_in);
      if (m_ix < 0 || sf_flag_loc[m_ix] == ACCURACY_POOR || sf_flag_loc[m_ix] == NOT_SET)
      	continue;

      shared_ptr<ParamSurface> surf2 = adj_surf[ki]->getSurface(0)->surface();
      bool cyllike2 = (surf2->instanceType() == Class_Cylinder ||
		       surf2->instanceType() == Class_Cone);
      int sf_flag = reg->defineSfFlag(numpt, 0, approx_tol_, num_in[m_ix], num2_in[m_ix],
				      avd[m_ix], (m_ix==0) ? cyllike : cyllike2);
      if (sf_flag == ACCURACY_POOR || sf_flag == NOT_SET)
	continue;
      
      vector<RevEngRegion*> added_adjacent;
      if (m_ix == 0)
	{
	  reg->includeAdjacentRegion(adj_surf[ki], maxd_loc[m_ix], avd_loc[m_ix],
				     num_in_loc[m_ix], num2_in_loc[ix], parvals[m_ix],
				     distang[m_ix], added_adjacent);
	  to_remove.push_back(adj_surf[ki]);
	}
      else
	{
	  adj_surf[ki]->includeAdjacentRegion(reg, maxd_loc[m_ix], avd_loc[m_ix],
					      num_in_loc[m_ix], num2_in_loc[m_ix],
					      parvals[m_ix], distang[m_ix], added_adjacent);
	  to_remove.push_back(reg);
	  reg = adj_surf[ki];
	  cyllike = cyllike2;
	}

      reg->setSurfaceFlag(sf_flag);
      int stop_break0 = 1;
    }

  for (size_t ki=0; ki<adj_nosurf.size(); ++ki)
    {
      int numpt = reg->numPoints() + adj_nosurf[ki]->numPoints();
      double maxd, maxd_loc, avd, avd_loc;
      int num_in, num2_in, num_in_loc, num2_in_loc;
      vector<double> parvals;
      vector<pair<double,double> > distang;
      int sf_flag_loc =
	reg->getGrowAccuracy(adj_nosurf[ki], approx_tol_, angtol, maxd, avd,
			     num_in, num2_in, maxd_loc, avd_loc,
			     num_in_loc, num2_in_loc, parvals, distang);
      int sf_flag = reg->defineSfFlag(numpt, 0, approx_tol_, num_in, num2_in,
				      avd, cyllike);
      if (sf_flag == ACCURACY_POOR || sf_flag == NOT_SET ||
	  sf_flag_loc == ACCURACY_POOR || sf_flag_loc == NOT_SET)
	continue;
      
      vector<RevEngRegion*> added_adjacent;
      reg->includeAdjacentRegion(adj_nosurf[ki], maxd_loc, avd_loc,
				 num_in_loc, num2_in_loc, parvals,
				 distang, added_adjacent);
      to_remove.push_back(adj_nosurf[ki]);
      reg->setSurfaceFlag(sf_flag);
    }

  vector<HedgeSurface*> to_remove_sfs;
  vector<RevEngEdge*> to_remove_edg;
  for (size_t ki=0; ki<to_remove.size(); ++ki)
    {
      if (to_remove[ki]->hasSurface())
	{
	  int numsf = to_remove[ki]->numSurface();
	  for (int ka=0; ka<numsf; ++ka)
	    to_remove_sfs.push_back(to_remove[ki]->getSurface(ka));
	}
      vector<RevEngEdge*> rev_edgs = to_remove[ki]->getAllRevEdges();
      for (size_t kr=0; kr<rev_edgs.size(); ++kr)
	{
	  RevEngRegion *adj1, *adj2;
	  rev_edgs[kr]->getAdjacent(adj1, adj2);
	  RevEngRegion *other = (adj1 == to_remove[ki]) ? adj2 : adj1;
	  other->removeRevEngEdge(rev_edgs[kr]);
	}
      to_remove_edg.insert(to_remove_edg.end(), rev_edgs.begin(), rev_edgs.end());
      to_remove[ki]->removeFromAdjacent();
      to_remove[ki]->clearRegionAdjacency();
    }

  double fac = 1.5;
  if (reg->numPoints() > (int)(fac*num_group_points) ||
      to_remove_sfs.size() > 0)
    {
  vector<RevEngEdge*> rev_edges_curr = reg->getAllRevEdges();
  for (size_t kj=0; kj<rev_edges_curr.size(); ++kj)
    rev_edges_curr[kj]->increaseExtendCount();
    }

  for (size_t ki=0; ki<to_remove_sfs.size(); ++ki)
    {
      size_t kj;
      for (kj=0; kj<surfaces_.size(); ++kj)
	if (surfaces_[kj].get() == to_remove_sfs[ki])
	  break;
      if ((int)kj <= ix)
	--ix;
     }
  
  for (size_t ki=0; ki<to_remove_edg.size(); ++ki)
    {
      size_t kj;
      for (kj=0; kj<edges_.size(); ++kj)
	if (edges_[kj].get() == to_remove_edg[ki])
	  break;
      if (kj < edges_.size())
	edges_.erase(edges_.begin()+kj);
    }
 
  if (to_remove.size() > 0)
    {
      for (size_t ki=0; ki<to_remove.size(); ++ki)
	{
	  to_remove[ki]->removeFromAdjacent();
	  to_remove[ki]->clearRegionAdjacency();
	}
 
      int dummy_ix = -1;
      updateRegionsAndSurfaces(dummy_ix, to_remove, to_remove_sfs);
    }
  
#ifdef DEBUG_SMALL
  std::ofstream of2("updated_small.g2");
  reg->writeRegionPoints(of2);
  reg->writeSurface(of2);
#endif
  int stop_break = 1;
}

//===========================================================================
// Service functionality for defineSmallRegionSurfaces
//

void groupAdjacentRegions(vector<RevEngRegion*>& regs,
			  vector<vector<RevEngRegion*> >& groups,
			  vector<int>& num_group)
{
  for (size_t kh=0; kh<regs.size(); ++kh)
    {
      vector<RevEngRegion*> currgroup;
      if (regs[kh]->visited())
	continue;
      regs[kh]->setVisited(true);
      currgroup.push_back(regs[kh]);
      for (size_t kh2=0; kh2<currgroup.size(); ++kh2)
	{
	  RevEngRegion *curr=currgroup[kh2];
	  vector<RevEngRegion*> adjacent;
	  curr->getAdjacentRegions(adjacent);
	  for (size_t kh3=0; kh3<adjacent.size(); ++kh3)
	    {
	      if (adjacent[kh3]->visited())
		continue;
	      if (adjacent[kh3]->hasSurface())
		continue;
	      if (adjacent[kh3]->hasAssociatedBlend())
		continue;
	      size_t kh4;
	      for (kh4=0; kh4<regs.size(); ++kh4)
		if (adjacent[kh3] == regs[kh4])
		  break;
	      if (kh4 < regs.size())
		{
		  adjacent[kh3]->setVisited(true);
		  currgroup.push_back(adjacent[kh3]);
		}
	    }
	}
      groups.push_back(currgroup);
    }

  for (size_t kh=0; kh<regs.size(); ++kh)
    regs[kh]->setVisited(false);

  num_group.resize(groups.size(), 0);
  for (size_t kr=0; kr<groups.size(); ++kr)
    for (size_t kh=0; kh<groups[kr].size(); ++kh)
      num_group[kr] += groups[kr][kh]->numPoints();
}

void disconnectRegion(RevEngRegion* reg, vector<RevEngRegion*>& removed_regs,
		      vector<shared_ptr<RevEngRegion> >& added_regs,
		      vector<RevEngRegion*>& nosf_reg)
{
  reg->removeFromAdjacent();
  reg->clearRegionAdjacency();
 auto it = std::find(nosf_reg.begin(), nosf_reg.end(), reg);
  if (it != nosf_reg.end())
    nosf_reg.erase(it);
  size_t ki;
  for (ki=0; ki<added_regs.size(); ++ki)
    if (added_regs[ki].get() == reg)
      break;
  if (ki == added_regs.size())
    removed_regs.push_back(reg);
  else
    added_regs.erase(added_regs.begin()+ki);
 }


void sortAlongAxis(vector<RevEngPoint*>& points, const Point& pos,
		   const Point& axis, vector<double>& ppar, double delta,
		   vector<vector<RevEngPoint*> >& sorted, size_t nsort,
		   vector<pair<RevEngPoint*, double> >& remaining,
		   double& tmin, double& tmax)
 {
  tmin = std::numeric_limits<double>::max();
  tmax = std::numeric_limits<double>::lowest();
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      Vector3D xyz = points[ki]->getPoint();
      Point point(xyz[0], xyz[1], xyz[2]);
      double tpar = (point - pos)*axis;
      tmin = std::min(tmin, tpar);
      tmax = std::max(tmax, tpar);

      if (ppar.size() == 0)
	{
	  sorted[nsort].push_back(points[ki]);
	}
      else
	{
	  size_t kr;
	  for (kr=0; kr<ppar.size() && ppar[kr]<tpar; ++kr);
	  if (kr > 0 && kr == ppar.size())
	    --kr;
	  if (fabs(tpar-ppar[kr]) < delta || (kr>0 && fabs(tpar-ppar[kr-1]) < delta) ||
	      (kr < ppar.size()-1 && fabs(tpar-ppar[kr+1]) < delta))
	    remaining.push_back(std::make_pair(points[ki],tpar));
	  else if (tpar < ppar[0]-delta)
	    sorted[nsort].push_back(points[ki]);
	  else if (tpar > ppar[ppar.size()-1]+delta)
	    sorted[nsort+ppar.size()].push_back(points[ki]);
	  else
	    sorted[nsort+kr].push_back(points[ki]);
	}
    }
 }

BoundingBox bBoxGroup(vector<RevEngPoint*>& group)
{
  BoundingBox bbox(3);
  for (size_t ki=0; ki<group.size(); ++ki)
    {
      Vector3D xyz = group[ki]->getPoint();
      Point pnt(xyz[0], xyz[1], xyz[2]);
      bbox.addUnionWith(pnt);
    }
  return bbox;
}

//===========================================================================
void RevEng::defineSmallRegionSurfaces()
//===========================================================================
{
  if (model_axis_.size() == 0)
    return;
  
#ifdef DEBUG_SMALL
  std::ofstream of1("nosurf_reg.g2");
  std::ofstream of1_2("withsurf_reg.g2");
  std::ofstream of1_3("assblend_reg.g2");
#endif

  vector<RevEngRegion*> nosf_reg;
  for (size_t kr=0; kr<regions_.size(); ++kr)
    {
      if (regions_[kr]->hasSurface())
	{
#ifdef DEBUG_SMALL
	  regions_[kr]->writeRegionPoints(of1_2);
#endif
	  continue;
	}
      if (regions_[kr]->hasAssociatedBlend())
	{
#ifdef DEBUG_SMALL
	  regions_[kr]->writeRegionPoints(of1_3);
#endif
	  continue;
	}
#ifdef DEBUG_SMALL
      regions_[kr]->writeRegionPoints(of1);
#endif
      nosf_reg.push_back(regions_[kr].get());
    }

  // Sort according to identified planes
  double angtol = 5.0*anglim_;
  int num_pt_lim = min_point_region_/3;
  double axisang = 2.0*angtol;
  double planeang = angtol;

  double diag = bbox_.low().dist(bbox_.high());

  // Count surfaces and points associated to each model axis
  size_t num_model_axis = model_axis_.size();
  vector<int> n_msurf(num_model_axis, 0), n_mpts(num_model_axis, 0);
  for (size_t ki=0; ki<num_model_axis; ++ki)
    {
      n_msurf[ki] += ((int)model_axis_[ki].plane_loc_.size() +
		      (int)model_axis_[ki].rotational_loc_.size());
      for (size_t kj=0; kj<model_axis_[ki].plane_loc_.size(); ++kj)
	n_mpts[ki] += model_axis_[ki].plane_loc_[kj].second;
      for (size_t kj=0; kj<model_axis_[ki].rotational_loc_.size(); ++kj)
	n_mpts[ki] += model_axis_[ki].rotational_loc_[kj].second;
    }

  double allnsurf = 0.0;
  double allnpt = 0.0;
  double lim1 = 0.5, lim2 = 0.25, lim3 = 0.1;
  for (size_t ki=0; ki<num_model_axis; ++ki)
    {
      allnsurf += n_msurf[ki];
      allnpt += n_mpts[ki];
    }

  vector<size_t> axis_split;
  axis_split.push_back(0);
  size_t kia;
  for (kia=0; kia<n_msurf.size() && (double)n_msurf[kia] >= lim1*(double)allnsurf &&
	 (double)n_mpts[kia] >= lim1*(double)allnpt; ++kia);
  if (kia > axis_split[axis_split.size()-1])
    axis_split.push_back(kia);
 for (; kia<n_msurf.size() && (double)n_msurf[kia] >= lim2*(double)allnsurf &&
	 (double)n_mpts[kia] >= lim2*(double)allnpt; ++kia);
  if (kia > axis_split[axis_split.size()-1])
    axis_split.push_back(kia);
 for (; kia<n_msurf.size() && (double)n_msurf[kia] >= lim3*(double)allnsurf &&
	 (double)n_mpts[kia] >= lim3*(double)allnpt; ++kia);
  if (kia > axis_split[axis_split.size()-1])
    axis_split.push_back(kia);
  if (num_model_axis > axis_split[axis_split.size()-1])
    axis_split.push_back(num_model_axis);
  
  for (size_t ax=1; ax<axis_split.size(); ++ax)
    {
      size_t ax1 = axis_split[ax-1];
      size_t ax2 = axis_split[ax];
      
      // Extract points according to identified axes and differ between planar
      // and rotational candidates
      vector<Point> axis_dir(ax2-ax1);
      for (size_t ki=ax1; ki<ax2; ++ki)
	axis_dir[ki-ax1] = model_axis_[ki].axis_;

      vector<vector<RevEngPoint*> > axis_group1(2*axis_dir.size());
      vector<vector<RevEngPoint*> > axis_group2(axis_dir.size());
      vector<RevEngPoint*> remaining;
      for (size_t ki=0; ki<nosf_reg.size(); ++ki)
	{
	  vector<vector<RevEngPoint*> > curr_group1, curr_group2;
	  vector<RevEngPoint*> curr_remaining;
	  (void)nosf_reg[ki]->sortByAxis(axis_dir, approx_tol_, axisang, planeang,
					 curr_group1, curr_group2, curr_remaining);
	  for (size_t kr=0; kr<curr_group1.size(); ++kr)
	    if (curr_group1[kr].size() > 0)
	      axis_group1[kr].insert(axis_group1[kr].end(), curr_group1[kr].begin(),
				     curr_group1[kr].end());
	  for (size_t kr=0; kr<curr_group2.size(); ++kr)
	    if (curr_group2.size() > 0)
	      axis_group2[kr].insert(axis_group2[kr].end(), curr_group2[kr].begin(),
				     curr_group2[kr].end());
	  if (curr_remaining.size() > 0)
	    remaining.insert(remaining.end(), curr_remaining.begin(), curr_remaining.end());
	}

#ifdef DEBUG_SMALL
      std::ofstream of2("small_axis_sort.g2");
      for (size_t kr=0; kr<axis_group1.size(); kr+=2)
	{
	  if (axis_group1[kr].size() > 0)
	    {
	      of2 << "400 1 0 4 0 255 0 255" << std::endl;
	      of2 << axis_group1[kr].size() << std::endl;
	      for (size_t kw=0; kw<axis_group1[kr].size(); ++kw)
		of2 << axis_group1[kr][kw]->getPoint() << std::endl;
	    }
	  if (axis_group1[kr+1].size() > 0)
	    {
	      of2 << "400 1 0 4 100 155 0 255" << std::endl;
	      of2 << axis_group1[kr+1].size() << std::endl;
	      for (size_t kw=0; kw<axis_group1[kr+1].size(); ++kw)
		of2 << axis_group1[kr+1][kw]->getPoint() << std::endl;
	    }
	}

      for (size_t kr=0; kr<axis_group2.size(); ++kr)
	{
	  if (axis_group2[kr].size() > 0)
	    {
	      of2 << "400 1 0 4 255 0 0 255" << std::endl;
	      of2 << axis_group2[kr].size() << std::endl;
	      for (size_t kw=0; kw<axis_group2[kr].size(); ++kw)
		of2 << axis_group2[kr][kw]->getPoint() << std::endl;
	    }
	}

      if (remaining.size() > 0)
	{
	  of2 << "400 1 0 4 0 0 255 255" << std::endl;
	  of2 << remaining.size() << std::endl;
	  for (size_t kw=0; kw<remaining.size(); ++kw)
	    of2 << remaining[kw]->getPoint() << std::endl;
	}
    
#endif
  
      // Sort according to identified axis and planes
      vector<SmallSurface> small_surf;
      vector<shared_ptr<RevEngRegion> > small_sf_reg;
      vector<shared_ptr<HedgeSurface> > small_sf_hedge;
     for (size_t ki=ax1, ki2=0; ki<ax2; ++ki, ++ki2)
	{
	  if (n_mpts[ki] < min_point_region_)
	    continue;
	  size_t num_planes = model_axis_[ki].plane_loc_.size();
	  Point axis = model_axis_[ki].axis_;
	  vector<shared_ptr<Plane> > axis_planes;

	  if (num_planes > 0)
	    {
	      axis_planes.resize(num_planes);
	      for (size_t kr=0; kr<num_planes; ++kr)
		axis_planes[kr] = shared_ptr<Plane>(new Plane(model_axis_[ki].plane_loc_[kr].first,
							      axis));
#ifdef DEBUG_SMALL
	      std::ofstream of4("axis_planes.g2");
	      for (size_t kr=0; kr<num_planes; ++kr)
		{
		  shared_ptr<Plane> tmp_plane(axis_planes[kr]->clone());
		  tmp_plane->setParameterBounds(-0.5*diag, -0.5*diag, 0.5*diag, 0.5*diag);
		  tmp_plane->writeStandardHeader(of4);
		  tmp_plane->write(of4);
		}
#endif
	    }
	  int ix = -1;
	  double minang = std::numeric_limits<double>::max();
	  for (int kb=0; kb<3; ++kb)
	    {
	      double ang = mainaxis_[kb].angle(axis);
	      ang = std::min(ang, M_PI-ang);
	      if (ang < minang)
		{
		  minang = ang;
		  ix = kb;
		}
	    }
	  ix = (ix > 0) ? ix-1 : 2;
	  Point Cx = mainaxis_[ix].cross(axis);
	  Cx.normalize();
	
	  Point midp(0.0, 0.0, 0.0);
	  vector<double> ppar;
	  double delta = 0.0001;
	  if (num_planes > 0)
	    {
	      double fac = 1.0/(double)(num_planes);
	      for (size_t kj=0; kj<num_planes; ++kj)
		midp += fac*model_axis_[ki].plane_loc_[kj].first;

	      ppar.resize(num_planes);
	      for (size_t kj=0; kj<num_planes; ++kj)
		ppar[kj] = (model_axis_[ki].plane_loc_[kj].first - midp)*axis;

	      // Sort locations along axis
	      for (size_t kj=0; kj<ppar.size(); ++kj)
		for (size_t kr=kj+1; kr<ppar.size(); ++kr)
		  if (ppar[kr] < ppar[kj])
		    {
		      std::swap(ppar[kj], ppar[kr]);
		      std::swap(axis_planes[kj], axis_planes[kr]);
		    }
	      delta = std::max(delta, 0.01*(ppar[ppar.size()-1]-ppar[0]));
	    }
	  
	  size_t nsort = ppar.size()+1;
	  vector<vector<RevEngPoint*> > planar_cand(2*nsort);
	  vector<vector<RevEngPoint*> > rotational_cand(nsort);
	  vector<pair<RevEngPoint*,double> > at_planes;
	  double tmin, tmax;
	  double all_min = std::numeric_limits<double>::max();
	  double all_max = std::numeric_limits<double>::lowest();
	  sortAlongAxis(axis_group1[2*ki2], midp, axis_dir[ki2], ppar, delta,
			planar_cand, 0, at_planes, tmin, tmax);
	  all_min = std::min(all_min, tmin);
	  all_max = std::max(all_max, tmax);
	  sortAlongAxis(axis_group1[2*ki2+1], midp, axis_dir[ki2], ppar, delta,
			planar_cand, nsort, at_planes, tmin, tmax);
	  all_min = std::min(all_min, tmin);
	  all_max = std::max(all_max, tmax);
	  sortAlongAxis(axis_group2[ki2], midp, axis_dir[ki2], ppar, delta,
			rotational_cand, 0, at_planes, tmin, tmax);
	  all_min = std::min(all_min, tmin);
	  all_max = std::max(all_max, tmax);
	  ppar.insert(ppar.begin(), all_min-delta);
	  ppar.push_back(all_max+delta);

#ifdef DEBUG_SMALL
	  std::ofstream of14("at_planes.g2");
	  of14 << "400 1 0 4 100 0 155 255" << std::endl;
	  of14 << at_planes.size() << std::endl;
	  for (size_t kw=0; kw<at_planes.size(); ++kw)
	    of14 << at_planes[kw].first->getPoint() << std::endl;
#endif

	  int min_point_reg2 = std::max(min_point_region_/10, 200);
	  int min_point_reg3 = std::max(min_point_region_/20, 200);
	  for (size_t kr=0; kr<rotational_cand.size(); ++kr)
	    {
      
#ifdef DEBUG_SMALL
	      std::ofstream of3("small_plane_int.g2");
	      if (planar_cand[kr].size() > 0)
		{
		  of3 << "400 1 0 4 0 255 0 255" << std::endl;
		  of3 << planar_cand[kr].size() << std::endl;
		  for (size_t kw=0; kw<planar_cand[kr].size(); ++kw)
		    of3 << planar_cand[kr][kw]->getPoint() << std::endl;
		}

	      if (planar_cand[nsort+kr].size() > 0)
		{
		  of3 << "400 1 0 4 100 155 0 255" << std::endl;
		  of3 << planar_cand[nsort+kr].size() << std::endl;
		  for (size_t kw=0; kw<planar_cand[nsort+kr].size(); ++kw)
		    of3 << planar_cand[nsort+kr][kw]->getPoint() << std::endl;
		}

	      if (rotational_cand[kr].size() > 0)
		{
		  of3 << "400 1 0 4 255 0 0 255" << std::endl;
		  of3 << rotational_cand[kr].size() << std::endl;
		  for (size_t kw=0; kw<rotational_cand[kr].size(); ++kw)
		    of3 << rotational_cand[kr][kw]->getPoint() << std::endl;
		}
#endif
	      // Start with rotational
	      if (rotational_cand[kr].size() >= min_point_reg2) //num_pt_lim) //min_point_region_)
		{
		  // Group adjacent
		  vector<vector<RevEngPoint*> > rot_groups;
		  RevEngUtils::identifyConGroups(rotational_cand[kr], rot_groups);
#ifdef DEBUG_SMALL
		  std::ofstream of5("con_rotate.g2");
		  for (size_t kh=0; kh<rot_groups.size(); ++kh)
		    {
		      of5 << "400 1 0 0" << std::endl;
		      of5 << rot_groups[kh].size() << std::endl;
		      for (size_t kv=0; kv<rot_groups[kh].size(); ++kv)
			of5 << rot_groups[kh][kv]->getPoint() << std::endl;
		    }
#endif

		  for (size_t kh=0; kh<rot_groups.size(); ++kh)
		    {
		      BoundingBox bb = bBoxGroup(rot_groups[kh]);
		      Point high = bb.high();
		      Point low = bb.low();
		      for (size_t kv=0; kv<=model_axis_[ki].rotational_loc_.size(); ++kv)
			{
			  Point loc;
			  if (kv < model_axis_[ki].rotational_loc_.size())
			    loc = model_axis_[ki].rotational_loc_[kv].first;
			  // Point low2 = loc + ((low-loc)*axis)*axis;
			  // Point high2 = loc + ((high-loc)*axis)*axis;
			  // if ((low-low2)*(high-high2) > 0.0 &&
			  // 	  std::min(low.dist(low2), high.dist(high2)) > low.dist(high))
			  // 	continue;   // Large distance between points and axis

			  // Check if the group can be associated any of the existing
			  // surfaces
			  size_t kj;
			  for (kj=0; kj<small_surf.size(); ++kj)
			    {
			      if ((int)ki != small_surf[kj].axis_ix_ ||
				  (int)kr != small_surf[kj].lev_ix_ ||
				  (int)kv != small_surf[kj].pos_ix_)
				continue;  // Not the same axis or interval

			      if (small_surf[kj].type_ != 3)
				continue;  // Surface type  not compatible

			      // Test distance
			      int all_in = 0, all_in2 = 0;
			      for (size_t kw=0; kw<small_surf[kj].surfs_.size(); ++kw)
				{
				  double maxdist, avdist;
				  int num_in, num2_in;
				  vector<pair<double, double> > dist_ang;
				  shared_ptr<ParamSurface> sf = small_surf[kj].surfs_[kw];
				  RevEngUtils::distToSurf(rot_groups[kh], sf,
							  approx_tol_, angtol, maxdist, avdist,
							  num_in, num2_in, dist_ang);
				  all_in += num_in;
				  all_in2 += num2_in;
				}
			      if (all_in2 >= std::max((int)rot_groups[kh].size()/2, 1))
				break;  // Preliminary match found
			    }

			  if (kj < small_surf.size())
			    {
			      small_surf[kj].addPoints(rot_groups[kh], bb);
			      break;
			    }

			  // Try to adapt a rotational surface
			  if ((int)rot_groups[kh].size() < min_point_reg3)
			    continue;

			  vector<shared_ptr<ElementarySurface> > sfs;
			  bool found = identifySmallRotational(rot_groups[kh], midp, loc,
							       axis, Cx,
							       ppar[kr], ppar[kr+1], sfs);
			  if (found)
			    {
#ifdef DEBUG_SMALL
			      std::ofstream of6("rot_sfs.g2");
			      for (size_t kw=0; kw<sfs.size(); ++kw)
				{
				  sfs[kw]->writeStandardHeader(of6);
				  sfs[kw]->write(of6);
				}
#endif
		      
			      SmallSurface small((int)ki, (int)kv, (int)kr, 3, sfs);
			      small.addPoints(rot_groups[kh], bb);
			      small_surf.push_back(small);
			      break;
			    }
			  int stop_break = 1;
			}
		    }
		}
	  
	      // Planar
	      for (int ka=0; ka<2; ++ka)
		{
		  if (planar_cand[ka*nsort+kr].size() < min_point_reg2) //num_pt_lim) //min_point_region_)
		    continue;

		  // Group adjacent
		  vector<vector<RevEngPoint*> > pla_groups;
		  RevEngUtils::identifyConGroups(planar_cand[ka*nsort+kr], pla_groups);
#ifdef DEBUG_SMALL
		  std::ofstream of9("pla_groups.g2");
		  for (size_t kh=0; kh<pla_groups.size(); ++kh)
		    {
		      of9 << "400 1 0 0" << std::endl;
		      of9 << pla_groups[kh].size() << std::endl;
		      for (size_t kv=0; kv<pla_groups[kh].size(); ++kv)
			of9 << pla_groups[kh][kv]->getPoint() << std::endl;
		    }
#endif

		  for (size_t kh=0; kh<pla_groups.size(); ++kh)
		    {
		      BoundingBox bb = bBoxGroup(pla_groups[kh]);
		      Point high = bb.high();
		      Point low = bb.low();
		  
		      // Check if the group can be associated any of the existing
		      // surfaces
		      size_t kj;
		      for (kj=0; kj<small_surf.size(); ++kj)
			{
			  if ((int)ki != small_surf[kj].axis_ix_ ||
			      (int)kr != small_surf[kj].lev_ix_)
			    continue;  // Not the same axis or interval
		      
			  if (small_surf[kj].type_ != ka+1)
			    continue;  // Surface type  not compatible

			  // Test distance
			  int all_in = 0, all_in2 = 0;
			  for (size_t kw=0; kw<small_surf[kj].surfs_.size(); ++kw)
			    {
			      double maxdist, avdist;
			      int num_in, num2_in;
			      vector<pair<double, double> > dist_ang;
			      shared_ptr<ParamSurface> sf = small_surf[kj].surfs_[kw];
			      RevEngUtils::distToSurf(pla_groups[kh], sf,
						      approx_tol_, angtol, maxdist, avdist,
						      num_in, num2_in, dist_ang);
			      all_in += num_in;
			      all_in2 += num2_in;
			    }
			  if (all_in2 >= std::max((int)pla_groups[kh].size()/2, 1))
			    break;  // Prelimenary match found
			}
		  
		      if (kj < small_surf.size())
			{
			  small_surf[kj].addPoints(pla_groups[kh], bb);
			  continue;
			}
		  
		      // Try to adapt a planar surface
		      if ((int)pla_groups[kh].size() < min_point_reg3)
			continue;

		      vector<shared_ptr<ElementarySurface> > sfs;
		      bool found = identifySmallPlanar(pla_groups[kh], midp, axis, Cx,
						       ppar[kr], ppar[kr+1], delta, sfs);
		      if (found)
			{
#ifdef DEBUG_SMALL
			  std::ofstream of9("pla_sfs.g2");
			  for (size_t kw=0; kw<sfs.size(); ++kw)
			    {
			      sfs[kw]->writeStandardHeader(of9);
			      sfs[kw]->write(of9);
			    }
#endif
		      
			  SmallSurface small((int)ki, -1, (int)kr, ka+1, sfs);
			  small.addPoints(pla_groups[kh], bb);
			  small_surf.push_back(small);
			}
		      int stop_break = 1;
		    }
		}
#ifdef DEBUG_SMALL
	      std::ofstream of8("small_sfs_curr.g2");
	      for (size_t kv=0; kv<small_surf.size(); ++kv)
		{
		  if (small_surf[kv].lev_ix_ != (int)kr)
		    continue;
		  for (size_t kw=0; kw<small_surf[kv].surfs_.size(); ++kw)
		    {
		      small_surf[kv].surfs_[kw]->writeStandardHeader(of8);
		      small_surf[kv].surfs_[kw]->write(of8);
		    }
		  for (size_t kw=0; kw<small_surf[kv].assos_points_.size(); ++kw)
		    {
		      of8 << "400 1 0 0" << std::endl;
		      of8 << small_surf[kv].assos_points_[kw].size() << std::endl;
		      for (size_t kw2=0; kw2<small_surf[kv].assos_points_[kw].size(); ++kw2)
			of8 << small_surf[kv].assos_points_[kw][kw2]->getPoint() << std::endl;
		    }
		}
#endif
	      int stop_break3 = 1;
	    }

#ifdef DEBUG_SMALL
	  std::ofstream of7("small_sfs_axis.g2");
	  for (size_t kv=0; kv<small_surf.size(); ++kv)
	    {
	      for (size_t kw=0; kw<small_surf[kv].surfs_.size(); ++kw)
		{
		  small_surf[kv].surfs_[kw]->writeStandardHeader(of7);
		  small_surf[kv].surfs_[kw]->write(of7);
		}
	      for (size_t kw=0; kw<small_surf[kv].assos_points_.size(); ++kw)
		{
		  of7 << "400 1 0 0" << std::endl;
		  of7 << small_surf[kv].assos_points_[kw].size() << std::endl;
		  for (size_t kw2=0; kw2<small_surf[kv].assos_points_[kw].size(); ++kw2)
		    of7 << small_surf[kv].assos_points_[kw][kw2]->getPoint() << std::endl;
		}
	    }
#endif
	  // Associate near-plane points to existing planes or create new
	  // planes as appopriate
 	  for (size_t kj=0; kj<axis_planes.size(); ++kj)
	    {
	      double tpar = ppar[kj+1];
	      vector<RevEngPoint*> plane_pts;
	      for (size_t kr=0; kr<at_planes.size(); ++kr)
		if (fabs(at_planes[kr].second - tpar) <= delta)
		  plane_pts.push_back(at_planes[kr].first);

	      Point loc = axis_planes[kj]->location();
	      Point norm = axis_planes[kj]->direction();
	      vector<HedgeSurface*> hedge_surf;
	      for (size_t kr=0; kr<surfaces_.size(); ++kr)
		{
		  shared_ptr<ParamSurface> surf = surfaces_[kr]->surface();
		  if (surf->instanceType() != Class_Plane)
		    continue;
		  shared_ptr<Plane> plane = dynamic_pointer_cast<Plane,ParamSurface>(surf);
		  Point loc2 = plane->location();
		  Point norm2 = plane->direction();
		  double ang = norm.angle(norm2);
		  ang = std::min(ang, M_PI-ang);
		  double dist = (loc2 - loc)*norm;
		  if (ang <= angtol && fabs(dist) <= approx_tol_)
		    hedge_surf.push_back(surfaces_[kr].get());
		}
#ifdef DEBUG_SMALL
	      std::ofstream of15("at_curr_planes.g2");
	      of15 << "400 1 0 4 0 255 0 255" << std::endl;
	      of15 << plane_pts.size() << std::endl;
	      for (size_t kv=0; kv<plane_pts.size(); ++kv)
		of15 << plane_pts[kv]->getPoint() << std::endl;
	      for (size_t kv=0; kv<hedge_surf.size(); ++kv)
		{
		  RevEngRegion *reg = hedge_surf[kv]->getRegion(0);
		  reg->writeRegionPoints(of15);
		  reg->writeSurface(of15);
		}
#endif
	      size_t curr_nmb_small = small_sf_reg.size();
	      planarAtPlane(axis_planes[kj], plane_pts, hedge_surf, small_sf_reg,
			    small_sf_hedge);
#ifdef DEBUG_SMALL
	      std::ofstream of16("at_curr_planes2.g2");
	      for (size_t kv=0; kv<hedge_surf.size(); ++kv)
		{
		  if (!hedge_surf[kv])
		    continue;
		  RevEngRegion *reg = hedge_surf[kv]->getRegion(0);
		  reg->writeRegionPoints(of16);
		  if (reg->hasSurface())
		    reg->writeSurface(of16);
		}
	      for (size_t kv=curr_nmb_small; kv<small_sf_reg.size(); ++kv)
		{
		  small_sf_reg[kv]->writeRegionPoints(of16);
		  small_sf_reg[kv]->writeSurface(of16);
		}
#endif
	      
	      int stop_plane = 0;
	    }
      
	  int stop_break2 = 1;
	}

     vector<vector<RevEngPoint*> > remain_groups;
     RevEngUtils::identifyConGroups(remaining, remain_groups);
#ifdef DEBUG_SMALL
     std::ofstream ofr("con_remain.g2");
     for (size_t kh=0; kh<remain_groups.size(); ++kh)
       {
	 ofr << "400 1 0 0" << std::endl;
	 ofr << remain_groups[kh].size() << std::endl;
	 for (size_t kv=0; kv<remain_groups[kh].size(); ++kv)
	   ofr << remain_groups[kh][kv]->getPoint() << std::endl;
       }
#endif

     for (size_t kh=0; kh<remain_groups.size(); ++kh)
       {
	 if ((int)remain_groups[kh].size() <  num_pt_lim)
	   continue;
	 vector<RevEngPoint*> in_pts, out_pts;
	 BoundingBox bb;
	 shared_ptr<ElementarySurface> elem =
	   defineElemSurf(remain_groups[kh], in_pts, bb, out_pts);
	 if (elem.get())
	   {
	     vector<shared_ptr<ElementarySurface> > sfs;
	     sfs.push_back(elem);
#ifdef DEBUG_SMALL
	     std::ofstream ofe("elem_remain_sf.g2");
	     for (size_t kw=0; kw<sfs.size(); ++kw)
	       {
		 sfs[kw]->writeStandardHeader(ofe);
		 sfs[kw]->write(ofe);
	       }
#endif
		      
	     SmallSurface small(-1, -1, -1, 4, sfs);
	     small.addPoints(remain_groups[kh], bb);
	     small_surf.push_back(small);
	   }
	 int stop_break_elem = 1;
       }
     

     // Extract identified surfaces
      vector<RevEngPoint*> non_assigned_pts;
      for (size_t kj=0; kj<small_surf.size(); ++kj)
	{
	  extractSmallSurfs(small_surf[kj], small_sf_reg, small_sf_hedge,
			    nosf_reg, non_assigned_pts);
	}

     // Remove empty regions
      for (int ka=(int)nosf_reg.size()-1; ka>=0; --ka)
	if (nosf_reg[ka]->numPoints() == 0)
	  {
	    for (int kb=(int)regions_.size()-1; kb>=0; --kb)
	      if (regions_[kb].get() == nosf_reg[ka])
		{
		  regions_.erase(regions_.begin()+kb);
		  break;
		}
	    nosf_reg.erase(nosf_reg.begin()+ka);
	  }
#ifdef DEBUG_SMALL
      std::ofstream of10("all_small_sfs.g2");
      for (size_t kj=0; kj<small_sf_reg.size(); ++kj)
	{
	  small_sf_reg[kj]->writeRegionPoints(of10);
	  small_sf_reg[kj]->writeSurface(of10);
	}
      std::ofstream of11("nosurf_reg2.g2");
      for (size_t kj=0; kj<nosf_reg.size(); ++kj)
	nosf_reg[kj]->writeRegionPoints(of11);
#endif

      // Ensure connectivity of remaining small regions
      vector<shared_ptr<RevEngRegion> > nosf_add;
      size_t nmb_nosf = nosf_reg.size();
     for (size_t kj=0; kj<nosf_reg.size(); ++kj)
	{
	  vector<vector<RevEngPoint*> > tmp_add;
	  int classtype = nosf_reg[kj]->getClassificationType();
	  nosf_reg[kj]->splitRegion(tmp_add);
	  for (size_t kr=0; kr<tmp_add.size(); ++kr)
	    {
	      shared_ptr<RevEngRegion> tmp_reg(new RevEngRegion(classtype,
								edge_class_type_,
								tmp_add[kr]));
	      nosf_add.push_back(tmp_reg);
	    }
	}

      for (size_t kj=0; kj<nosf_add.size(); ++kj)
	nosf_reg.push_back(nosf_add[kj].get());
  
#ifdef DEBUG_SMALL
      std::ofstream of12("nosurf_reg3.g2");
      for (size_t kj=0; kj<nosf_reg.size(); ++kj)
	nosf_reg[kj]->writeRegionPoints(of12);
#endif

      vector<RevEngRegion*> include_reg;
      integrateInSmallSurfs(small_sf_reg, nosf_reg, include_reg);
      for (size_t kj=0; kj<include_reg.size(); ++kj)
	{
	  include_reg[kj]->removeFromAdjacent();
	  include_reg[kj]->clearRegionAdjacency();
	}
      for (int ka=(int)include_reg.size()-1; ka>=0; --ka)
	{
	  auto it = std::find(nosf_reg.begin(), nosf_reg.end(), include_reg[ka]);
	  if (it == nosf_reg.end())
	    continue; // Should not happen
	  int ix = it - nosf_reg.begin();
	  if (ix >= (int)nmb_nosf)
	    {
	      // nosf_add[ix-(int)nmb_nosf]->removeFromAdjacent();
	      // nosf_add[ix-(int)nmb_nosf]->clearRegionAdjacency();
	      nosf_add.erase(nosf_add.begin()+ix-(int)nmb_nosf);
	    }
	  else
	    {
	      for (int kb=(int)regions_.size()-1; kb>=0; --kb)
		if (regions_[kb].get() == include_reg[ka])
		  {
		    // regions_[kb]->removeFromAdjacent();
		    // regions_[kb]->clearRegionAdjacency();
		    regions_.erase(regions_.begin()+kb);
		    break;
		  }
	      nmb_nosf--;
	    }
	  nosf_reg.erase(nosf_reg.begin()+ix);
	}

#ifdef DEBUG_SMALL
      std::ofstream of13("all_small_sfs2.g2");
      for (size_t kj=0; kj<small_sf_reg.size(); ++kj)
	{
	  small_sf_reg[kj]->writeRegionPoints(of13);
	  small_sf_reg[kj]->writeSurface(of13);
	}
      std::ofstream of14("nosurf_reg4.g2");
      for (size_t kj=0; kj<nosf_reg.size(); ++kj)
	nosf_reg[kj]->writeRegionPoints(of14);
#endif

      if (small_sf_reg.size() > 0)
	regions_.insert(regions_.end(), small_sf_reg.begin(), small_sf_reg.end());
      if (nosf_add.size() > 0)
	regions_.insert(regions_.end(), nosf_add.begin(), nosf_add.end());
      if (small_sf_hedge.size() > 0)
	surfaces_.insert(surfaces_.end(), small_sf_hedge.begin(), small_sf_hedge.end());
    }
}

//===========================================================================
struct FittingResults
{
  double maxd, avd;
  int inside, inside_2;
  vector<double> param;
  vector<pair<double,double> > distang;
};

void collectAdjacentSurfRegs(vector<RevEngPoint*>& points,
			     vector<pair<RevEngRegion*, int> >& adj_regs)
{
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      vector<RevEngRegion*> adj = points[ki]->adjacentRegsWithSurf();
      for (size_t kj=0; kj<adj.size(); ++kj)
	{
	  size_t kr;
	  for (kr=0; kr<adj_regs.size(); ++kr)
	    if (adj_regs[kr].first == adj[kj])
	      break;
	  if (kr < adj_regs.size())
	    adj_regs[kr].second++;
	  else
	    adj_regs.push_back(std::make_pair(adj[kj], 1));
	}
    }
}

double mostCompatibleDir(const Point& curr, vector<Point>& vecs1,
			 vector<Point>& vecs2, Point& dir)
{
  double min_ang = M_PI;
  for (size_t ki=0; ki<vecs1.size(); ++ki)
    {
      double ang = curr.angle(vecs1[ki]);
      ang = std::min(ang, M_PI-ang);
      if (ang < min_ang)
	{
	  min_ang = ang;
	  dir = vecs1[ki];
	}
    }
  for (size_t ki=0; ki<vecs2.size(); ++ki)
    {
      double ang = curr.angle(vecs2[ki]);
      ang = std::min(ang, M_PI-ang);
      if (ang < min_ang)
	{
	  min_ang = ang;
	  dir = vecs2[ki];
	}
    }
  return min_ang;
}

void identifyDistPoints(vector<RevEngPoint*>& points, vector<pair<double,double> >& distang,
			double lim, vector<RevEngPoint*>& in_pts,
			vector<RevEngPoint*>& out_pts)
{
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      if (distang[ki].first <= lim)
	in_pts.push_back(points[ki]);
      else
	out_pts.push_back(points[ki]);
    }
}

shared_ptr<ElementarySurface>
RevEng::defineElemSurf(vector<RevEngPoint*>& points,
		       vector<RevEngPoint*>& in_points,
		       BoundingBox& bbox, vector<RevEngPoint*>& remain)
//===========================================================================
{
  shared_ptr<ElementarySurface> dummy_sf;
  
  // Collect direction information in adjacent surfaces
  vector<pair<RevEngRegion*,int> > adj_sf_regs;
  collectAdjacentSurfRegs(points, adj_sf_regs);

  vector<Point> vecs1;
  vector<int> num_reg_pts;
  for (size_t ki=0; ki<adj_sf_regs.size(); ++ki)
    {
      shared_ptr<ParamSurface> surf =
	adj_sf_regs[ki].first->getSurface(0)->surface();
      shared_ptr<ElementarySurface> elem =
	dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf);
      if (!elem.get())
	continue;
      vecs1.push_back(elem->direction());
      num_reg_pts.push_back(adj_sf_regs[ki].first->numPoints());
    }

  vector<Point> vecs2(model_axis_.size());
  vector<int> num_axis_pt(model_axis_.size(), 0);
  for (size_t ki=0; ki<model_axis_.size(); ++ki)
    {
      vecs2[ki] = model_axis_[ki].axis_;
      for (size_t kj=0; kj<model_axis_[ki].plane_loc_.size(); ++kj)
	num_axis_pt[ki] += model_axis_[ki].plane_loc_[kj].second;
      for (size_t kj=0; kj<model_axis_[ki].rotational_loc_.size(); ++kj)
	num_axis_pt[ki] += model_axis_[ki].rotational_loc_[kj].second;
    }

  // Prepare points for surface generation
  vector<RevEngPoint*> dummy_pts;
  vector<pair<vector<RevEngPoint*>::iterator,vector<RevEngPoint*>::iterator> > pts_it;
  pts_it.push_back(std::make_pair(points.begin(), points.end()));
  vector<Point> pts(points.size());
  bbox = BoundingBox(3);
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      Vector3D xyz = points[ki]->getPoint();
      pts[ki] = Point(xyz[0], xyz[1], xyz[2]);
      bbox.addUnionWith(pts[ki]);
    }
  Point low = bbox.low();
  Point high = bbox.high();
  
#ifdef DEBUG_SMALL
  std::ofstream of("remaining_group.g2");
  of << "400 1 0 4 55 100 100 255" << std::endl;
  of << points.size() << std::endl;
  for (size_t ki=0; ki<points.size(); ++ki)
    of << points[ki]->getPoint() << std::endl;
#endif
  // Collect information about feasible approximating surfaces
  double eps = 1.0e-9;
  double angtol = 5.0*anglim_;
  double tol2 = 2.0*approx_tol_;
  double avd_fac = 2.0;
  double in_fac = 0.2;
  double fac1 = 0.9;
  double fac2 = 0.9;
  double red_fac = 0.75;
  double ang_lim2 = 0.25*M_PI;

  vector<shared_ptr<ElementarySurface> > elem_sfs;
  vector<int> flag;
  vector<FittingResults> acc_fit;
  vector<vector<RevEngPoint*> > used_pts;
  vector<vector<RevEngPoint*> > not_used_pts;
  
  // Plane
  Point init_norm = points[0]->getTriangNormal();
  Point pos1, norm1, Cx1, Cy1;
  RevEngUtils::computePlane(pts, init_norm, mainaxis_, pos1, norm1, Cx1, Cy1);
  if (norm1.length() < 0.1)
    return dummy_sf;   // The length should be one, but could be zero if the
  // approximation fails
  shared_ptr<Plane> plane1(new Plane(pos1, norm1, Cx1));

#ifdef DEBUG_SMALL
  plane1->writeStandardHeader(of);
  plane1->write(of);
#endif
  FittingResults acc1;
  RevEngUtils::distToSurf(points.begin(), points.end(), plane1, approx_tol_,
			  acc1.maxd, acc1.avd, acc1.inside, acc1.inside_2,
			  acc1.param, acc1.distang, angtol);
  if (acc1.avd < avd_fac*approx_tol_ && (double)acc1.inside_2 > in_fac*(double)points.size())
    {
      int surf_flag1 = regions_[0]->defineSfFlag((int)points.size(), 0, approx_tol_,
						 acc1.inside, acc1.inside_2,
						 acc1.avd, false);
      // Possible surface. Check with alternative plane normal
      Point dir;
      double ang = mostCompatibleDir(plane1->direction(), vecs1, vecs2, dir);
      if (ang > eps && ang <= angtol)
	{
	  // Update plane with new normal
	  shared_ptr<Plane> plane2 =
	    RevEngUtils::planeWithAxis(points, dir, plane1->location(), mainaxis_);
	  FittingResults acc2;
	  RevEngUtils::distToSurf(points.begin(), points.end(), plane1, approx_tol_,
				  acc2.maxd, acc2.avd, acc2.inside, acc2.inside_2,
				  acc2.param, acc2.distang, angtol);
	  int surf_flag2 = regions_[0]->defineSfFlag((int)points.size(), 0, approx_tol_,
						     acc2.inside, acc2.inside_2,
						     acc2.avd, false);
 	  if (surf_flag2 <= surf_flag1 && acc2.avd <= fac1*acc1.avd &&
	      (double)acc2.inside_2 >= fac2*(double)acc1.inside_2)
	    {
#ifdef DEBUG_SMALL
	      plane2->writeStandardHeader(of);
	      plane2->write(of);
#endif
	      std::swap(acc1, acc2);
	      std::swap(plane1, plane2);
	    }
	}

      bool use_reduced = false;
      if (acc1.maxd > tol2)
	{
	  // Remove most distant points and make a refit
	  vector<RevEngPoint*> in_pts, out_pts;
	  identifyDistPoints(points, acc1.distang, tol2, in_pts, out_pts);
	  if ((double)in_pts.size() > red_fac*(double)points.size())
	    {
	      vector<Point> in_pts2(in_pts.size());
	      for (size_t ki=0; ki<in_pts.size(); ++ki)
		{
		  Vector3D xyz = in_pts[ki]->getPoint();
		  in_pts2[ki] = Point(xyz[0], xyz[1], xyz[2]);
		}
	      Point pos2, norm2, Cx2, Cy2;
	      RevEngUtils::computePlane(in_pts2, plane1->direction(), mainaxis_, pos2,
					norm2, Cx2, Cy2);
	      shared_ptr<Plane> plane2(new Plane(pos2, norm2, Cx2));
	  
	      FittingResults acc2;
	      RevEngUtils::distToSurf(in_pts.begin(), in_pts.end(), plane2, approx_tol_,
				      acc2.maxd, acc2.avd, acc2.inside, acc2.inside_2,
				      acc2.param, acc2.distang, angtol);
	      int surf_flag2 = regions_[0]->defineSfFlag((int)in_pts.size(), 0, approx_tol_,
							 acc2.inside, acc2.inside_2,
							 acc2.avd, false);
	      double frac1 = (double)acc1.inside_2/(double)points.size();
	      double frac2 = (double)acc2.inside_2/(double)in_pts.size();
	      if (surf_flag2 < ACCURACY_POOR &&
		  surf_flag2 <= surf_flag1 && acc2.avd < acc1.avd && frac2 > frac1)
		{
#ifdef DEBUG_SMALL
		  plane2->writeStandardHeader(of);
		  plane2->write(of);
#endif
		  use_reduced = true;
		  elem_sfs.push_back(plane2);
		  flag.push_back(surf_flag2);
		  acc_fit.push_back(acc2);
		  used_pts.push_back(in_pts);
		  not_used_pts.push_back(out_pts);
		}
	    }
	}
      if (surf_flag1 < ACCURACY_POOR && (!use_reduced))
	{
	  elem_sfs.push_back(plane1);
	  flag.push_back(surf_flag1);
	  acc_fit.push_back(acc1);
	  used_pts.push_back(points);
	  not_used_pts.push_back(dummy_pts);
	}
    }

  // Cylinder and cone
  Point axis2, Cx2, Cy2, pos2;
  double rad2;
  RevEngUtils::computeAxis(pts_it, axis2, Cx2, Cy2);
  RevEngUtils::computeCylPosRadius(pts_it, low, high, axis2, Cx2, Cy2, pos2, rad2);
  shared_ptr<ElementarySurface> csf1(new Cylinder(rad2, pos2, axis2, Cx2));
#ifdef DEBUG_SMALL
  csf1->writeStandardHeader(of);
  csf1->write(of);
#endif
  
  FittingResults acc2;
  RevEngUtils::distToSurf(points.begin(), points.end(), csf1, approx_tol_,
			  acc2.maxd, acc2.avd, acc2.inside, acc2.inside_2,
			  acc2.param, acc2.distang, angtol);

  shared_ptr<Cone> cone =
    RevEngUtils::coneWithAxis(points, csf1->direction(), low, high, mainaxis_);
  if (fabs(cone->getConeAngle()) > angtol)
    {
      FittingResults acc3;
      RevEngUtils::distToSurf(points.begin(), points.end(), cone, approx_tol_,
			  acc3.maxd, acc3.avd, acc3.inside, acc3.inside_2,
			  acc3.param, acc3.distang, angtol);
      if (acc3.avd < acc2.avd && acc3.inside_2 > acc2.inside_2)
	{
#ifdef DEBUG_SMALL
	  cone->writeStandardHeader(of);
	  cone->write(of);
#endif
	  std::swap(acc2, acc3);
	  csf1 = cone;
	}
    }
      
  if (acc2.avd < avd_fac*approx_tol_ && (double)acc2.inside_2 > in_fac*(double)points.size())
    {
      int surf_flag2 = regions_[0]->defineSfFlag((int)points.size(), 0, approx_tol_,
						 acc2.inside, acc2.inside_2,
						 acc2.avd, true);
      
      // Possible surface. Check with alternative plane normal
      Point dir;
      double ang = mostCompatibleDir(csf1->direction(), vecs1, vecs2, dir);
      if (ang > eps && ang < ang_lim2)
	{
	  // Update rotational surface with new axis
	  shared_ptr<ElementarySurface> csf2 =
	    RevEngUtils::cylinderWithAxis(points, dir, low, high, mainaxis_);
	  
	  FittingResults acc3;
	  RevEngUtils::distToSurf(points.begin(), points.end(), csf2, approx_tol_,
				  acc3.maxd, acc3.avd, acc3.inside, acc3.inside_2,
				  acc3.param, acc3.distang, angtol);
 				  
	  shared_ptr<Cone> cone2 =
	    RevEngUtils::coneWithAxis(points, csf2->direction(), low, high, mainaxis_);
	  if (fabs(cone2->getConeAngle()) > angtol)
	    {
	      FittingResults acc4;
	      RevEngUtils::distToSurf(points.begin(), points.end(), cone, approx_tol_,
				      acc4.maxd, acc4.avd, acc4.inside, acc4.inside_2,
				      acc4.param, acc4.distang, angtol);
	      if (acc4.avd < acc3.avd && acc4.inside_2 > acc3.inside_2)
		{
		  std::swap(acc3, acc4);
		  csf2 = cone2;
		}
	    }
	  
	  int surf_flag3 = regions_[0]->defineSfFlag((int)points.size(), 0, approx_tol_,
						     acc3.inside, acc3.inside_2,
						     acc3.avd, true);
 	  if (surf_flag3 <= surf_flag2 && acc3.avd <= fac1*acc2.avd &&
	      (double)acc3.inside_2 >= fac2*(double)acc2.inside_2)
	    {
#ifdef DEBUG_SMALL
	      csf2->writeStandardHeader(of);
	      csf2->write(of);
#endif
	      std::swap(acc2, acc3);
	      std::swap(csf1, csf2);
	    }
	}

      bool use_reduced = false;
      if (acc2.maxd > tol2)
	{
	  // Remove most distant points and make a refit
	  vector<RevEngPoint*> in_pts, out_pts;
	  identifyDistPoints(points, acc2.distang, tol2, in_pts, out_pts);
	  if ((double)in_pts.size() > red_fac*(double)points.size())
	    {
	      vector<pair<vector<RevEngPoint*>::iterator,vector<RevEngPoint*>::iterator> > pts_it2;
	      pts_it2.push_back(std::make_pair(in_pts.begin(), in_pts.end()));
	      Point axis3, Cx3, Cy3, pos3;
	      double rad3;
	      RevEngUtils::computeAxis(pts_it2, axis3, Cx3, Cy3);
	      RevEngUtils::computeCylPosRadius(pts_it2, low, high, axis3, Cx3, Cy3,
					       pos3, rad3);
	      shared_ptr<ElementarySurface> csf2(new Cylinder(rad3, pos3, axis3, Cx3));
	      
	  
	      FittingResults acc3;
	      RevEngUtils::distToSurf(in_pts.begin(), in_pts.end(), csf2, approx_tol_,
				      acc3.maxd, acc3.avd, acc3.inside, acc3.inside_2,
				      acc3.param, acc3.distang, angtol);
	      shared_ptr<Cone> cone2 =
		RevEngUtils::coneWithAxis(in_pts, csf2->direction(), low, high, mainaxis_);
	      if (fabs(cone2->getConeAngle()) > angtol)
		{
		  FittingResults acc4;
		  RevEngUtils::distToSurf(in_pts.begin(), in_pts.end(), cone2, approx_tol_,
				      acc4.maxd, acc4.avd, acc4.inside, acc4.inside_2,
				      acc4.param, acc4.distang, angtol);
		  if (acc4.avd < acc3.avd && acc4.inside_2 > acc3.inside_2)
		    {
		      std::swap(acc3, acc4);
		      csf2 = cone2;
		    }
		}
	      
	      int surf_flag3 = regions_[0]->defineSfFlag((int)in_pts.size(), 0, approx_tol_,
							 acc3.inside, acc3.inside_2,
							 acc3.avd, true);
	      double frac1 = (double)acc2.inside_2/(double)points.size();
	      double frac2 = (double)acc3.inside_2/(double)in_pts.size();
	      if (surf_flag3 < ACCURACY_POOR &&
		  surf_flag3 <= surf_flag2 && acc3.avd < acc2.avd && frac2 > frac1)
		{
#ifdef DEBUG_SMALL
		  csf2->writeStandardHeader(of);
		  csf2->write(of);
#endif
		  use_reduced = true;
		  elem_sfs.push_back(csf2);
		  flag.push_back(surf_flag3);
		  acc_fit.push_back(acc3);
		  used_pts.push_back(in_pts);
		  not_used_pts.push_back(out_pts);
		}
	    }
	}
      if (surf_flag2 < ACCURACY_POOR && (!use_reduced))
	{
	  elem_sfs.push_back(csf1);
	  flag.push_back(surf_flag2);
	  acc_fit.push_back(acc2);
	  used_pts.push_back(points);
	  not_used_pts.push_back(dummy_pts);
	}
    }

  // Sphere
  int ix = 0;
  for (size_t ki=1; ki<num_axis_pt.size(); ++ki)
    if (num_axis_pt[ki] > num_axis_pt[ix])
      ix = (int)ki;
  Point axis3 = vecs2[ix];

  shared_ptr<Sphere> sph1 = RevEngUtils::sphereWithAxis(points, axis3, mainaxis_);
  FittingResults acc3;
  RevEngUtils::distToSurf(points.begin(), points.end(), sph1, approx_tol_,
			  acc3.maxd, acc3.avd, acc3.inside, acc3.inside_2,
			  acc3.param, acc3.distang, angtol);
  
#ifdef DEBUG_SMALL
  sph1->writeStandardHeader(of);
  sph1->write(of);
#endif
  if (acc3.avd < avd_fac*approx_tol_ && (double)acc3.inside_2 > in_fac*(double)points.size())
    {
      int surf_flag3 = regions_[0]->defineSfFlag((int)points.size(), 0, approx_tol_,
						 acc3.inside, acc3.inside_2,
						 acc3.avd, false);
      bool use_reduced = false;
      if (acc3.maxd > tol2)
	{
	  // Remove most distant points and make a refit
	  vector<RevEngPoint*> in_pts, out_pts;
	  identifyDistPoints(points, acc3.distang, tol2, in_pts, out_pts);
	  if ((double)in_pts.size() > red_fac*(double)points.size())
	    {
	      shared_ptr<Sphere> sph2 = RevEngUtils::sphereWithAxis(in_pts, axis3, mainaxis_);
	      FittingResults acc4;
	      RevEngUtils::distToSurf(in_pts.begin(), in_pts.end(), sph2, approx_tol_,
				      acc4.maxd, acc4.avd, acc4.inside, acc4.inside_2,
				      acc4.param, acc4.distang, angtol);
	      int surf_flag4 = regions_[0]->defineSfFlag((int)in_pts.size(), 0, approx_tol_,
							 acc3.inside, acc3.inside_2,
							 acc3.avd, false);
	      double frac1 = (double)acc3.inside_2/(double)points.size();
	      double frac2 = (double)acc4.inside_2/(double)in_pts.size();
	      if (surf_flag4 < ACCURACY_POOR &&
		  surf_flag4 <= surf_flag3 && acc4.avd < acc3.avd && frac2 > frac1)
		{
#ifdef DEBUG_SMALL
		  sph2->writeStandardHeader(of);
		  sph2->write(of);
#endif
		  use_reduced = true;
		  elem_sfs.push_back(sph2);
		  flag.push_back(surf_flag4);
		  acc_fit.push_back(acc4);
		  used_pts.push_back(in_pts);
		  not_used_pts.push_back(out_pts);
		}
	    }
	}
      if (surf_flag3 < ACCURACY_POOR && (!use_reduced))
	{
	  elem_sfs.push_back(sph1);
	  flag.push_back(surf_flag3);
	  acc_fit.push_back(acc3);
	  used_pts.push_back(points);
	  not_used_pts.push_back(dummy_pts);
	}
    }

  // Select result
  ix = -1;
  int sf_flag = NOT_SET;
  double inside_fac = 0.1;
  for (size_t ki=0; ki<elem_sfs.size(); ++ki)
    {
      if ((double)acc_fit[ki].inside < inside_fac*(double)used_pts[ki].size())
	continue;
      if (flag[ki] < sf_flag)
	{
	  ix = (int)ki;
	  sf_flag = flag[ki];
	}
      else if ((flag[ki] == sf_flag ||
		(std::min(flag[ki],sf_flag) == 1 && std::max(flag[ki],sf_flag) == 2)) &&
	       acc_fit[ki].avd < acc_fit[ix].avd &&
	       acc_fit[ki].inside_2 > acc_fit[ix].inside_2)
	{
	  ix = (int)ki;
	  sf_flag = flag[ki];
	}
    }
      
  if (ix >= 0)
    {
      in_points = used_pts[ix];
      remain = not_used_pts[ix];
      return elem_sfs[ix];
    }
  else
    return dummy_sf;
}

//===========================================================================
void RevEng::planarAtPlane(shared_ptr<Plane> axis_plane,
			   vector<RevEngPoint*>& points,
			   vector<HedgeSurface*>& sfs,
			   vector<shared_ptr<RevEngRegion> >& plane_sf_reg,
			   vector<shared_ptr<HedgeSurface> >& plane_sf_hedge)
//===========================================================================
{
  if (sfs.size() == 0 || points.size() == 0)
    return;

  RevEngRegion *reg0 = sfs[0]->getRegion(0);
  Point axis = axis_plane->direction();
  Point loc = axis_plane->location();
  int min_num_pts = min_point_region_/5;
  double int_tol = 1.0e-4;

  // Check accuracy with respect to axis plane
  double angtol = 5.0*anglim_;
  double maxdist, avdist;
  int num_in, num2_in;
  vector<double> parvals;
  vector<pair<double,double> > distang;
  vector<RevEngPoint*> in, out;
  RevEngUtils::distToSurf(points.begin(), points.end(), axis_plane, approx_tol_,
			  maxdist, avdist, num_in, num2_in, in, out,
			  parvals, distang, angtol);
  int sf_flag = reg0->defineSfFlag((int)points.size(), 0, approx_tol_, num_in,
				  num2_in, avdist, false);

  vector<bool> same_plane(sfs.size(), true);
  double dtol = std::min(1.0e-4, 0.1*approx_tol_);
  double atol = 0.001;
  for (size_t ki=0; ki<sfs.size(); ++ki)
    {
      shared_ptr<ParamSurface> surf = sfs[ki]->surface();
      shared_ptr<Plane> plane = dynamic_pointer_cast<Plane,ParamSurface>(surf);
      if (!plane.get())
	continue;
      double ang = plane->direction().angle(axis);
      ang = std::min(ang, M_PI-ang);
      double dd = (plane->location() - loc)*axis;
      if (ang > atol || dd > dtol)
	same_plane[ki] = false;
    }

  // Check accuracy with respect to existing planar surfaces
  vector<int> flag_sfs(sfs.size(), sf_flag);
  vector<int> flag_axis(sfs.size(), NOT_SET);
  vector<bool> has_revedgs(sfs.size(), false);
  vector<int> num_points(sfs.size());
  vector<RevEngRegion*> sfs_reg(sfs.size());
  for (size_t ki=0; ki<sfs.size(); ++ki)
    {
      RevEngRegion *reg = sfs[ki]->getRegion(0);
      sfs_reg[ki] = reg;
      num_points[ki] = reg->numPoints();
      has_revedgs[ki] = reg->hasRevEdges();
      if (same_plane[ki])
	{
	  flag_axis[ki] = reg->getSurfaceFlag();
	  continue;
	}
      
      shared_ptr<ParamSurface> surf = sfs[ki]->surface();
      shared_ptr<Plane> plane = dynamic_pointer_cast<Plane,ParamSurface>(surf);
      if (!plane.get())
	continue;
      // To avoid bounded plane
      shared_ptr<Plane> plane2(new Plane(plane->location(), plane->direction()));
      double maxd, avd;
      int inside, inside2;
      vector<double> param;
      vector<pair<double,double> > d_a;
      vector<RevEngPoint*> in0, out0;
      RevEngUtils::distToSurf(points.begin(), points.end(), plane2, approx_tol_,
			      maxd, avd, inside, inside2, in0, out0,
			      param, d_a, angtol);
      flag_sfs[ki] = reg0->defineSfFlag((int)points.size(), 0, approx_tol_, inside,
					inside2, avd, false);
      if (sf_flag >= ACCURACY_POOR && flag_sfs[ki] >= ACCURACY_POOR)
	continue;  // Not compatible

      if (sf_flag < ACCURACY_POOR)
	{
	  double maxd2, avd2;
	  int inside_2, inside2_2;
	  vector<double> param2;
	  vector<pair<double,double> > d_a2;
	  vector<RevEngPoint*> in2, out2;
	  RevEngUtils::distToSurf(reg->pointsBegin(), reg->pointsEnd(), axis_plane,
				  approx_tol_, maxd2, avd2, inside_2, inside2_2, in2, out2,
				  param2, d_a2, angtol);
	  flag_axis[ki] = reg->defineSfFlag(0, approx_tol_, inside_2,
					    inside2_2, avd2, false);
	  int stop_break = 1;
	}
      int stop_break2 = 1;
    }

  // Check if the planar points can be integrated in any of the existing surfaces
  int ix = -1;
  for (size_t ki=0; ki<sfs.size(); ++ki)
    {
      if (flag_sfs[ki] < ACCURACY_POOR &&
	  (ix < 0 || flag_sfs[ki] < flag_sfs[ix] ||
	   (flag_sfs[ki] == flag_sfs[ix] && has_revedgs[ki] && (!has_revedgs[ix])) ||
	   (flag_sfs[ki] == flag_sfs[ix] && has_revedgs[ki] == has_revedgs[ix] &&
	    num_points[ki] > num_points[ix])))
	ix = ki;
    }

  if (ix < 0 && (sf_flag >= ACCURACY_POOR || (int)points.size() < min_num_pts))
    return;

  // Collect all points and parameterize on selected surface
  shared_ptr<ParamSurface> sel_sf = (ix < 0) ? axis_plane : sfs[ix]->surface();
  vector<RevEngPoint*> all_pts;
  for (size_t ki=0; ki<sfs.size(); ++ki)
    if (flag_axis[ki] < ACCURACY_POOR)
      all_pts.insert(all_pts.end(), sfs_reg[ki]->pointsBegin(),
		     sfs_reg[ki]->pointsEnd());
  all_pts.insert(all_pts.end(), points.begin(), points.end());
  double maxd_all, avd_all;
  int num_in_all, num2_in_all;
  vector<double> parvals_all;
  vector<pair<double,double> > distang_all;
  vector<RevEngPoint*> in_all, out_all;
  RevEngUtils::distToSurf(all_pts.begin(), all_pts.end(), sel_sf, approx_tol_,
			  maxd_all, avd_all, num_in_all, num2_in_all, in_all,
			  out_all, parvals_all, distang_all, angtol);
  int sf_flag_all = reg0->defineSfFlag((int)all_pts.size(), 0, approx_tol_, 
				       num_in_all, num2_in_all, avd_all, false);
  if (sf_flag_all >= ACCURACY_POOR)
    return;  // Do nothing

  // Assign parameter value and accuracy to points. This might have to be redone
  for (size_t ki=0; ki<all_pts.size(); ++ki)
    {
      all_pts[ki]->setPar(Vector2D(parvals_all[2*ki], parvals_all[2*ki+1]));
      all_pts[ki]->setSurfaceDist(distang_all[ki].first, distang_all[ki].second);
    }
  
  vector<vector<RevEngPoint*> > con_all;
  RevEngUtils::identifyConGroups(all_pts, con_all);
#ifdef DEBUG_SMALL
  std::ofstream of1("sep_pts_all.g2");
  for (size_t ki=0; ki<con_all.size(); ++ki)
    {
      of1 << "400 1 0 0" << std::endl;
      of1 << con_all[ki].size() << std::endl;
      for (size_t kj=0; kj<con_all[ki].size(); ++kj)
	of1 << con_all[ki][kj]->getPoint() << std::endl;
    }
#endif

  for (size_t ki=0; ki<con_all.size(); ++ki)
    {
      // Assemble accuracy information
      double maxdc = 0.0, avdc = 0.0;
      int num_inc = 0, num2_inc = 0, ang_in = 0;
      double fac = 1.0/(double)con_all[ki].size();
      for (size_t kj=0; kj<con_all[ki].size(); ++kj)
	{
	  double dd, ang;
	  con_all[ki][kj]->getSurfaceDist(dd, ang);
	  maxdc = std::max(maxdc, dd);
	  avdc += fac*dd;
	  if (dd <= approx_tol_ && ang <= angtol)
	    num_inc++;
	  if (dd <= approx_tol_)
	    num2_inc++;
	  if (ang <= angtol)
	    ang_in++;
	}
      int flagc = reg0->defineSfFlag((int)con_all[ki].size(), 0, approx_tol_,
				     num_inc, num2_inc, avdc, false);
      int min_nmb = std::max((int)con_all[ki].size()/4, 10);
      if (ang_in < min_nmb && num_inc < min_nmb)
	flagc = ACCURACY_POOR;
      if (flagc < ACCURACY_POOR && (int)con_all[ki].size() > min_num_pts/2)
	{
	  // Check connection to existing surface
	  vector<int> num_sf_points(sfs.size(), 0);
	  for (size_t kj=0; kj<con_all[ki].size(); ++kj)
	    {
	      RevEngRegion *curr_reg = con_all[ki][kj]->region();
	      for (size_t kr=0; kr<sfs.size(); ++kr)
		if (curr_reg == sfs_reg[kr])
		  {
		    num_sf_points[kr]++;
		    break;
		  }
	    }

	  // Decide on corresponding surface
	  int ix2 = -1;
	  int num_corr = 0;
	  for (size_t kr=0; kr<num_sf_points.size(); ++kr)
	    if (num_sf_points[kr] > num_corr)
	      {
		ix2 = (int)kr;
		num_corr = num_sf_points[kr];
	      }

	  // Remove the associated points from their previous region
	  // Collect points with respect to regions
	  vector<vector<RevEngPoint*> > reg_points;
	  for (size_t kj=0; kj<con_all[ki].size(); ++kj)
	    {
	      RevEngRegion *curreg = con_all[ki][kj]->region();
	      if ((ix2 >=0 && curreg != sfs_reg[ix2]) || ix2 < 0)
		{
		  size_t kr;
		  for (kr=0; kr<reg_points.size(); ++kr)
		      if (reg_points[kr][0]->region() == curreg)
			break;
		    if (kr < reg_points.size())
		      reg_points[kr].push_back(con_all[ki][kj]);
		    else
		      {
			vector<RevEngPoint*> curr_reg_pt;
			curr_reg_pt.push_back(con_all[ki][kj]);
			reg_points.push_back(curr_reg_pt);
		      }
		}
	    }
	  for (size_t kr=0; kr<reg_points.size(); ++kr)
	    {
	      RevEngRegion *currreg = reg_points[kr][0]->region();
	      currreg->removePoints(reg_points[kr]);
	      for (size_t kh=0; kh<reg_points[kr].size(); ++kh)
		reg_points[kr][kh]->unsetRegion();
	    }

	  if (ix2 < 0)
	    {
	      // Define new region and surface
	      shared_ptr<RevEngRegion> plane_reg(new RevEngRegion(reg0->getClassificationType(),
								  edge_class_type_,
								  con_all[ki]));
	      shared_ptr<Plane> plane2(axis_plane->clone());
	      shared_ptr<HedgeSurface> hedge(new HedgeSurface(sel_sf, plane_reg.get()));
	      plane_reg->setHedge(hedge.get());
	      plane_reg->updateInfo(approx_tol_, angtol);
	      plane_reg->setSurfaceFlag(flagc);
	      plane_sf_reg.push_back(plane_reg);
	      plane_sf_hedge.push_back(hedge);
	    }
	  else
	    {
	      // Enhance region with added points
	      vector<RevEngPoint*> added_pts;
	      for (size_t kj=0; kj<con_all[ki].size(); ++kj)
		{
		  if (con_all[ki][kj]->region() != sfs_reg[ix2])
		    added_pts.push_back(con_all[ki][kj]);
		}
	      if (added_pts.size() > 0)
		{
		  sfs_reg[ix2]->addPointsToGroup(added_pts, approx_tol_, angtol, false);
		  vector<RevEngEdge*> rev_edgs = sfs_reg[ix2]->getAllRevEdges();
		  for (size_t kj=0; kj<rev_edgs.size(); ++kj)
		    rev_edgs[kj]->increaseExtendCount();
		}
	      if (ix2 != ix)
		{
		  // Reset parameterization
		  sfs_reg[ix2]->parameterizePoints(approx_tol_, angtol);
		}
	    }
	  int stop_break5 = 1;
	}
      else
	{
	  // Reset current parameterization, if any
	  double eps = 1.0e-6;
	  double upar, vpar, dist;
	  Point close;
	  Point norm1, norm2, norm3;
	  double ang, ang2;
	  for (size_t kj=0; kj<con_all[ki].size(); ++kj)
	    {
	      RevEngPoint *currpt = con_all[ki][kj];
	      RevEngRegion *curreg = currpt->region();
	      if (curreg->hasSurface())
		{
		  shared_ptr<ParamSurface> surf = curreg->getSurface(0)->surface();
		  Vector3D xyz = currpt->getPoint();
		  Point pnt(xyz[0], xyz[1], xyz[2]);
		  surf->closestPoint(pnt, upar, vpar, close, dist, eps);
		  surf->normal(norm1, upar, vpar);
		  norm2 = currpt->getLocFuncNormal();
		  norm3 = currpt->getTriangNormal();
		  ang = norm1.angle(norm2);
		  ang2 = norm1.angle(norm3);
		  ang = std::min(std::min(M_PI-ang, ang), std::min(M_PI-ang2,ang2));
		  currpt->setPar(Vector2D(upar, vpar));
		  currpt->setSurfaceDist(dist, ang);
		}
	    }
	}
      int stop_break4 = 1;
    }

  vector<vector<RevEngPoint*> > separate_groups;
  vector<HedgeSurface*> removed_surfs;
  vector<RevEngEdge*> removed_edgs;
  for (size_t ki=0; ki<sfs_reg.size(); ++ki)
    {
     if (sfs_reg[ki]->numPoints() == 0)
	{
	  vector<RevEngEdge*> rev_edgs = sfs_reg[ki]->getAllRevEdges();
	  for (size_t kr=0; kr<rev_edgs.size(); ++kr)
	    {
	      size_t kj;
	      for (kj=0; kj<edges_.size(); ++kj)
		if (edges_[kj].get() == rev_edgs[kr])
		  break;
	      if (kj < edges_.size())
		edges_.erase(edges_.begin()+kj);
	    }
	  
	  // Remove from region pool
	  for (size_t kr=0; kr<regions_.size(); ++kr)
	    if (regions_[kr].get() == sfs_reg[ki])
	      {
		regions_[kr]->removeFromAdjacent();
		regions_[kr]->clearRegionAdjacency();
		regions_.erase(regions_.begin()+kr);
		break;
	      }
	  removed_surfs.push_back(sfs[ki]);
	  sfs[ki] = 0;
	}
     else
       {
	 // Split region if not connected
	 sfs_reg[ki]->splitRegion(separate_groups);
	 if (sfs_reg[ki]->numPoints() < min_num_pts)
	   {
	     vector<RevEngEdge*> rev_edgs = sfs_reg[ki]->getAllRevEdges();
	     if (rev_edgs.size() > 0)
	       removed_edgs.insert(removed_edgs.end(), rev_edgs.begin(), rev_edgs.end());
	     removed_surfs.push_back(sfs_reg[ki]->getSurface(0));
	     sfs_reg[ki]->clearSurface();
	   }
       }
    }
  surfaceExtractOutput(-1, separate_groups, removed_surfs);
  for (size_t ki=0; ki<removed_edgs.size(); ++ki)
    {
      size_t kj;
      for (kj=0; kj<edges_.size(); ++kj)
	if (edges_[kj].get() == removed_edgs[ki])
	  break;
      if (kj < edges_.size())
	edges_.erase(edges_.begin()+kj);
    }
      
  
  int stop_break3 = 1;
}

//===========================================================================
void RevEng::integrateInSmallSurfs(vector<shared_ptr<RevEngRegion> >& small_sf_reg,
				   vector<RevEngRegion*>& nosf_reg,
				   vector<RevEngRegion*>& include_reg)
//===========================================================================
{
  // Sort small regions according to size (could possibly include also accuracy)
  for (size_t ki=0; ki<small_sf_reg.size(); ++ki)
    for (size_t kj=ki+1; kj<small_sf_reg.size(); ++kj)
      if  (small_sf_reg[kj]->numPoints() > small_sf_reg[ki]->numPoints())
	std::swap(small_sf_reg[ki], small_sf_reg[kj]);

  // Check accuracy of not assigned regions with respect to the small surfaces
  double angtol = 5.0*anglim_;
  for (size_t ki=0; ki<nosf_reg.size(); ++ki)
    {
      BoundingBox bb1 = nosf_reg[ki]->getBbox();
      vector<double> maxdist, avdist;
      vector<int> num_in, num2_in;
      vector<vector<double> > parvals;
      vector<vector<pair<double,double> > > distang;
      vector<int> sfflag;
      vector<size_t> cand_ix;
      for (size_t kj=0; kj<small_sf_reg.size(); ++kj)
	{
	  BoundingBox bb2 = small_sf_reg[kj]->getBbox();
	  if (!bb1.overlaps(bb2, approx_tol_))
	      continue;
	  shared_ptr<ParamSurface> surf = small_sf_reg[kj]->getSurface(0)->surface();
	  bool cyllike = (surf->instanceType() == Class_Cylinder ||
			  surf->instanceType() == Class_Cone);
	  double maxd, avd;
	  int inside, inside2;
	  vector<double> param;
	  vector<pair<double,double> > d_a;
	  vector<RevEngPoint*> in, out;
	  RevEngUtils::distToSurf(nosf_reg[ki]->pointsBegin(),
				  nosf_reg[ki]->pointsEnd(), surf, approx_tol_,
				  maxd, avd, inside, inside2, in, out,
				  param, d_a, angtol);
	  int sf_flag = small_sf_reg[kj]->defineSfFlag(0, approx_tol_, inside,
						       inside2, avd, cyllike);
	  if (sf_flag < ACCURACY_POOR)
	    {
	      // Check "closeness"
	      double dom1[4];
	      small_sf_reg[kj]->getDomain(dom1);
	      RectDomain sf_dom(Vector2D(dom1[0],dom1[2]), Vector2D(dom1[1],dom1[3]));
	      Vector2D c1(param[0], param[1]);
	      Vector2D c2(param[0], param[1]);
	      for (size_t kr=2; kr<param.size(); kr+=2)
		{
		  c1[0] = std::min(c1[0], param[kr]);
		  c2[0] = std::max(c2[0], param[kr]);
		  c1[1] = std::min(c1[1], param[kr+1]);
		  c2[1] = std::max(c2[1], param[kr+1]);
		}
	      RectDomain pnt_dom(c1, c2);
	      double udel = 0.25*(dom1[1]-dom1[0]);
	      double vdel = 0.25*(dom1[3]-dom1[2]);
	      BoundingBox bb1 = small_sf_reg[kj]->boundingBox();
	      double ddel = 0.25*bb1.low().dist(bb1.high());
	      BoundingBox bb2 = nosf_reg[ki]->boundingBox();
	      if (sf_dom.overlap(pnt_dom, udel, vdel) &&
		  bb1.overlaps(bb2, ddel))
		{
		  maxdist.push_back(maxd);
		  avdist.push_back(avd);
		  num_in.push_back(inside);
		  num2_in.push_back(inside2);
		  parvals.push_back(param);
		  distang.push_back(d_a);
		  sfflag.push_back(sf_flag);
		  cand_ix.push_back(kj);
		}
	    }
	}
      if (cand_ix.size() > 0)
	{
	  // Sort candidates
	  double eps = 0.01*approx_tol_;
	  vector<int> perm(cand_ix.size());
	  for (size_t kj=0; kj<perm.size(); ++kj)
	    perm[kj] = kj;
	  for (size_t kj=0; kj<perm.size(); ++kj)
	    for (size_t kr=kj+1; kr<perm.size(); ++kr)
	      if (sfflag[perm[kr]] < sfflag[perm[kj]] ||
		  (sfflag[perm[kr]] == sfflag[perm[kj]] &&
		   avdist[perm[kr]] < avdist[perm[kj]]-eps))
		std::swap(perm[kj], perm[kr]);

	  vector<RevEngRegion*> dummy;
	  small_sf_reg[cand_ix[perm[0]]]->includeAdjacentRegion(nosf_reg[ki],
								maxdist[perm[0]],
								avdist[perm[0]],
								num_in[perm[0]],
								num2_in[perm[0]],
								parvals[perm[0]],
								distang[perm[0]],
								dummy);
	  include_reg.push_back(nosf_reg[ki]);
	}
      
      int stop_break1 = 1;
    }
  int stop_break2 = 1;
}

bool checkAccuracy(RevEngRegion* reg, int min_num_pt, double tol, double angtol,
		   shared_ptr<ElementarySurface> surf,
		   vector<vector<RevEngPoint*> >& points, vector<BoundingBox>& bb,
		   size_t& num_points, double& maxdist, double& avdist, int& num_in,
		   int& num2_in, vector<vector<double> >& parvals,
		   vector<vector<pair<double,double> > >& dist_ang,
		   vector<RevEngPoint*>& non_assigned_pts)
{
  num_points = 0;
  maxdist = 0.0;
  avdist = 0.0;
  num_in = 0;
  num2_in = 0;
  bool type_cyl = (surf->instanceType() == Class_Cylinder ||
		   surf->instanceType() == Class_Cone);
  for (size_t kj=0; kj<points.size(); )
    {
      if (points[kj].size() == 0)
	{
	  points.erase(points.begin()+kj);
	  continue;
	}
      double maxd, avd;
      int inside, inside2;
      vector<double> param;
      vector<pair<double,double> > d_a;
      vector<RevEngPoint*> in, out;
      RevEngUtils::distToSurf(points[kj].begin(), points[kj].end(),
			      surf, tol, maxd, avd, inside,
			      inside2, in, out, param, d_a, angtol);
      int sf_flag = reg->defineSfFlag((int)points[kj].size(), 0,
				      tol, inside, inside2, avd, type_cyl);
      if (sf_flag >= ACCURACY_POOR)
	{
	  non_assigned_pts.insert(non_assigned_pts.end(), points[kj].begin(),
				  points[kj].end());
	  points.erase(points.begin()+kj);
	  bb.erase(bb.begin()+kj);
	}
      else
	{
	  maxdist = std::max(maxdist, maxd);
	  avdist += avd;
	  num_in += inside;
	  num2_in += inside2;
	  parvals.push_back(param);
	  dist_ang.push_back(d_a);
	  num_points += points[kj].size();
	  ++kj;
	}
    }
  avdist /= (double)parvals.size();
  if ((int)num_points < min_num_pt)
    {
      for (size_t kj=0; kj<points.size(); ++kj)
	if (points[kj].size() > 0)
	  non_assigned_pts.insert(non_assigned_pts.end(), points[kj].begin(),
				  points[kj].end());
      return false;
    }
  return true;
}

//===========================================================================
void RevEng::extractSmallSurfs(SmallSurface& small_surf,
			       vector<shared_ptr<RevEngRegion> >& small_sf_reg,
			       vector<shared_ptr<HedgeSurface> >& small_sf_hedge,
			       vector<RevEngRegion*>& nosf_reg,
			       vector<RevEngPoint*>& non_assigned_pts)
//===========================================================================
{
  double eps = 1.0e-6;
  double angtol = 5.0*anglim_;
  int min_num_pt = min_point_region_/10;
  int classtype = CLASSIFICATION_UNDEF;
  
  vector<shared_ptr<ElementarySurface> > surfs = small_surf.surfs_;
  vector<vector<vector<RevEngPoint*> > > points(surfs.size());
  vector<vector<BoundingBox> > bb(surfs.size());
  if (surfs.size() == 1)
    {
      points[0] = small_surf.assos_points_;
      bb[0] = small_surf.bbox_;
    }
  else
    {
      // Distribute points according to surfaces
      for (size_t ki=0; ki<surfs.size(); ++ki)
	{
	  points[ki].resize(small_surf.assos_points_.size());
	  bb[ki].resize(small_surf.assos_points_.size());
	  for (size_t kr=0; kr<bb[ki].size(); ++kr)
	    bb[ki][kr] = BoundingBox(3);
	}

      for (size_t ki=0; ki<small_surf.assos_points_.size(); ++ki)
	{
	  for (size_t kj=0; kj<small_surf.assos_points_[ki].size(); ++kj)
	    {
	      Vector3D xyz = small_surf.assos_points_[ki][kj]->getPoint();
	      Point pos(xyz[0], xyz[1], xyz[2]);

	      double mind = std::numeric_limits<double>::max();
	      int min_ix = -1;
	      for (size_t kr=0; kr<surfs.size(); ++kr)
		{
		  double upar, vpar, dist;
		  Point close;
		  surfs[kr]->closestPoint(pos, upar, vpar, close, dist, eps);
		  if (dist < mind)
		    {
		      mind = dist;
		      min_ix = (int)kr;
		    }
		}
	      if (min_ix >= 0)
		{
		  points[min_ix][ki].push_back(small_surf.assos_points_[ki][kj]);
		  bb[min_ix][ki].addUnionWith(pos);
		}
	    }
	}
    }

  for (size_t ki=0; ki<surfs.size(); ++ki)
    {
      size_t npt = 0;
      for (size_t kj=0; kj<points[ki].size(); ++kj)
	npt += points[ki][kj].size();
      if ((int)npt < min_num_pt)
	{
	  for (size_t kj=0; kj<points[ki].size(); ++kj)
	    if (points[ki][kj].size() > 0)
	      non_assigned_pts.insert(non_assigned_pts.end(), points[ki][kj].begin(),
				      points[ki][kj].end());
	  continue;
	}
      
      shared_ptr<ElementarySurface> curr_sf;
      bool type_cyl = (surfs[ki]->instanceType() == Class_Cylinder ||
		       surfs[ki]->instanceType() == Class_Cone);
      if (points[ki].size() > 1)
	{
	  // Recompute free surface parameters
	  vector<RevEngPoint*> all_points;
	  for (size_t kj=0; kj<points[ki].size(); ++kj)
	    all_points.insert(all_points.end(), points[ki][kj].begin(),
			      points[ki][kj].end());
	  double diag = bbox_.low().dist(bbox_.high());
	  curr_sf = RevEngUtils::elemsurfWithAxis(surfs[ki], all_points, mainaxis_, diag);
	}
      else
	curr_sf = surfs[ki];

      // Check accuracy and extract point groups with a large distance.
      size_t num_points = 0;
      double maxdist = 0.0, avdist = 0.0;
      int num_in = 0, num2_in = 0;
      vector<vector<double> > parvals;
      vector<vector<pair<double,double> > > dist_ang;
      bool OK = checkAccuracy(nosf_reg[0], min_num_pt, approx_tol_, angtol,
			      curr_sf, points[ki], bb[ki],
			      num_points, maxdist, avdist, num_in, num2_in,
			      parvals, dist_ang, non_assigned_pts);
		    
      if (!OK)
	continue;

#ifdef DEBUG_SMALL
      std::ofstream of1("cand_sf_points.g2");
      curr_sf->writeStandardHeader(of1);
      curr_sf->write(of1);
      for (size_t kv=0; kv<points[ki].size(); ++kv)
	{
	  if (points[ki][kv].size() > 0)
	    {
	      of1 << "400 1 0 0" << std::endl;
	      of1 << points[ki][kv].size() << std::endl;
	      for (size_t kw=0; kw<points[ki][kv].size(); ++kw)
		of1 << points[ki][kv][kw]->getPoint() << std::endl;
	    }
	}
#endif
      double distlim = 2.0*approx_tol_;
      if (maxdist > distlim)
	{
	  // Remove the most distant points and try again
	  int num2 = 0;
	  for (size_t kj=0; kj<points[ki].size(); ++kj)
	    {
	      bb[ki][kj] = BoundingBox(3);
	      for (size_t kr=0; kr<points[ki][kj].size(); )
		{
		  if (dist_ang[kj][kr].first > distlim)
		    points[ki][kj].erase(points[ki][kj].begin()+kr);
		  else
		    {
		      Vector3D xyz = points[ki][kj][kr]->getPoint();
		      Point pos(xyz[0], xyz[1], xyz[2]);
		      bb[ki][kj].addUnionWith(pos);
		      num2++;
		      ++kr;
		    }
		}
	    }
	  if (num2 < min_num_pt)
	    continue;

	  // Recompute free surface parameters
	  vector<RevEngPoint*> all_points;
	  for (size_t kj=0; kj<points[ki].size(); ++kj)
	    all_points.insert(all_points.end(), points[ki][kj].begin(),
			      points[ki][kj].end());
	  double diag = bbox_.low().dist(bbox_.high());
	  shared_ptr<ElementarySurface> curr_sf2;
	  curr_sf2 = RevEngUtils::elemsurfWithAxis(curr_sf, all_points, mainaxis_, diag);
	  std::swap(curr_sf, curr_sf2);

	  OK = checkAccuracy(nosf_reg[0], min_num_pt, approx_tol_, angtol,
			     curr_sf, points[ki], bb[ki],
			     num_points, maxdist, avdist, num_in, num2_in,
			     parvals, dist_ang, non_assigned_pts);
		    
	  if (!OK)
	    continue;
	  
#ifdef DEBUG_SMALL
	  std::ofstream of3("cand_sf_points2.g2");
	  curr_sf->writeStandardHeader(of3);
	  curr_sf->write(of3);
	  for (size_t kv=0; kv<points[ki].size(); ++kv)
	    {
	      if (points[ki][kv].size() > 0)
		{
		  of3 << "400 1 0 0" << std::endl;
		  of3 << points[ki][kv].size() << std::endl;
		  for (size_t kw=0; kw<points[ki][kv].size(); ++kw)
		    of3 << points[ki][kv][kw]->getPoint() << std::endl;
		}
	    }
#endif
	}
      
      // Compute parameter domains
      vector<RectDomain> dom(parvals.size());
      for (size_t kj=0; kj<parvals.size(); ++kj)
	{
	  Vector2D c1(parvals[kj][0], parvals[kj][1]);
	  Vector2D c2(parvals[kj][0], parvals[kj][1]);
	  for (size_t kr=2; kr<parvals[kj].size(); kr+=2)
	    {
	      c1[0] = std::min(c1[0], parvals[kj][kr]);
	      c2[0] = std::max(c2[0], parvals[kj][kr]);
	      c1[1] = std::min(c1[1], parvals[kj][kr+1]);
	      c2[1] = std::max(c2[1], parvals[kj][kr+1]);
	    }
	  dom[kj] = RectDomain(c1, c2);
	}

      // Sort groups according to size
      RectDomain totdom = dom[0];
      vector<double> urange(dom.size());
      vector<double> vrange(dom.size());
      vector<double> diag(dom.size());
      vector<size_t> perm(points[ki].size());
      for (size_t kj=0; kj<perm.size(); ++kj)
	{
	  perm[kj] = kj;
	  totdom.addUnionWith(dom[kj]);
	  urange[kj] = dom[kj].umax() - dom[kj].umin();
	  vrange[kj] = dom[kj].vmax() - dom[kj].vmin();
	  diag[kj] = bb[ki][kj].low().dist(bb[ki][kj].high());
	}
      for (size_t kj=0; kj<perm.size(); ++kj)
	for (size_t kr=kj+1; kr<perm.size(); ++kr)
	  if (points[ki][perm[kr]].size() > points[ki][perm[kj]].size())
	    std::swap(perm[kr], perm[kj]);

      vector<size_t> first_group;
      size_t ix = 0;
      first_group.push_back(ix);
      size_t kj=0, kr=0, kh=0, kv=0;
      for (kj=0; kj<perm.size(); kj=kr)
	{
	  for (kr=kj+1, kv=kj+2; kr<perm.size(); )
	    {
	      for (kh=ix; kh<kr; ++kh)
		{
		  double udel = 0.125*(urange[perm[kr]]+urange[perm[kh]]);
		  double vdel = 0.125*(vrange[perm[kr]]+vrange[perm[kh]]);
		  double ddel = 0.125*(diag[perm[kr]]+diag[perm[kh]]);
		  if (dom[perm[kh]].overlap(dom[perm[kr]], udel, vdel) &&
		      bb[ki][perm[kh]].overlaps(bb[ki][perm[kr]], ddel))
		    break;
		}
	      if (kh >= kr)
		{
		  // Not close
		  if (kv >= perm.size())
		    {
		      ix = kr;
		      first_group.push_back(ix);
		      break;
		    }
		  if (kr < perm.size()-1)
		    std::swap(perm[kr], perm[kv]);
		  ++kv;
		}
	      else
		{
		  ++kr;
		  kv = kr+1;
		}
	    }
	}
      first_group.push_back(perm.size());
      
#ifdef DEBUG_SMALL
      std::ofstream of2("cand_sf_points4.g2");
      curr_sf->writeStandardHeader(of2);
      curr_sf->write(of2);
      for (size_t kj=1; kj<first_group.size(); ++kj)
	{
	  size_t num = 0;
	  for (size_t kv=first_group[kj-1]; kv<first_group[kj]; ++kv)
	    num += points[ki][perm[kv]].size();
	  of2 << "400 1 0 0" << std::endl;
	  of2 << num << std::endl;
	  for (size_t kv=first_group[kj-1]; kv<first_group[kj]; ++kv)
	    for (size_t kw=0; kw<points[ki][perm[kv]].size(); ++kw)
	      of2 << points[ki][perm[kv]][kw]->getPoint() << std::endl;
	}
#endif
       for (size_t kj=1; kj<first_group.size(); ++kj)
	{
	  size_t num = 0;
	  for (size_t kv=first_group[kj-1]; kv<first_group[kj]; ++kv)
	    num += points[ki][perm[kv]].size();

	  if (num < min_num_pt)
	    {
	      for (size_t kv=first_group[kj-1]; kv<first_group[kj]; ++kv)
		if (points[ki][perm[kv]].size() > 0)
	      non_assigned_pts.insert(non_assigned_pts.end(), 
				      points[ki][perm[kv]].begin(),
				      points[ki][perm[kv]].end());
	    }
	  else
	    {
	      // Create region with associated surface
	      // First remove the associated points from their previous region
	      // Collect points with respect to regions
	      vector<vector<RevEngPoint*> > reg_points;
	      for (size_t kv=first_group[kj-1]; kv<first_group[kj]; ++kv)
		for (size_t kw=0; kw<points[ki][perm[kv]].size(); ++kw)
		  {
		    RevEngPoint *currpt = points[ki][perm[kv]][kw];
		    RevEngRegion *reg = currpt->region();
		    size_t kr;
		    for (kr=0; kr<reg_points.size(); ++kr)
		      if (reg_points[kr][0]->region() == reg)
			break;
		    if (kr < reg_points.size())
		      reg_points[kr].push_back(currpt);
		    else
		      {
			vector<RevEngPoint*> curr_reg_pt;
			curr_reg_pt.push_back(currpt);
			reg_points.push_back(curr_reg_pt);
		      }
		  }

	      for (size_t kr=0; kr<reg_points.size(); ++kr)
		{
		  RevEngRegion *reg = reg_points[kr][0]->region();
		  reg->removePoints(reg_points[kr]);
		}

	      // Collect points and set parameter value
	      vector<RevEngPoint*> surf_pts;
	      for (size_t kv=first_group[kj-1]; kv<first_group[kj]; ++kv)
		for (size_t kw=0; kw<points[ki][perm[kv]].size(); ++kw)
		  {
		    RevEngPoint *currpt = points[ki][perm[kv]][kw];
		    double u1 = parvals[perm[kv]][2*kw];
		    double v1 = parvals[perm[kv]][2*kw+1];
		    currpt->setPar(Vector2D(u1,v1));
		    double dd = dist_ang[perm[kv]][kw].first;
		    double ang = dist_ang[perm[kv]][kw].second;
		    currpt->setSurfaceDist(dd, ang);
		    surf_pts.push_back(currpt);
		  }

	      shared_ptr<RevEngRegion> small_reg(new RevEngRegion(classtype,
								  edge_class_type_,
								  surf_pts));
	      shared_ptr<HedgeSurface> hedge(new HedgeSurface(curr_sf, small_reg.get()));
	      small_reg->setHedge(hedge.get());
	      small_reg->updateInfo(approx_tol_, angtol);
	      int sf_flag = small_reg->defineSfFlag(0, approx_tol_, num_in, num2_in,
						   avdist, type_cyl);
	      small_reg->setSurfaceFlag(sf_flag);
	      small_sf_reg.push_back(small_reg);
	      small_sf_hedge.push_back(hedge);

	      
	      int stop_remove = 1;
	    }
	}
     int stop_break = 1;
    }
}



//===========================================================================
bool RevEng::identifySmallRotational(vector<RevEngPoint*>& points,
				     Point midp, Point loc0, Point axis, Point Cx,
				     double ppar1, double ppar2,
				     vector<shared_ptr<ElementarySurface> >& sfs)
//===========================================================================
{
  double eps = 1.0e-9;
  bool found = false;
  double num_fac = 0.5;
  double anglim2 = 5*anglim_; //10*anglim_;
  int num_pt_lim = min_point_region_/20;
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group_points;
  group_points.push_back(std::make_pair(points.begin(), points.end()));

  if (loc0.dimension() == 0)
    {
      // Define point from cylinder algorithm
      Point low = bbox_.low();
      Point high = bbox_.high();
      shared_ptr<Cylinder> tmp_cyl =
	RevEngUtils::cylinderWithAxis(points, axis, low, high, mainaxis_);
      loc0 = tmp_cyl->location();
    }
  Point loc = loc0 + ((midp - loc0)*axis)*axis;

  vector<Point> rotated;
  RevEngUtils::rotateToPlane(group_points, Cx, axis, loc, rotated);
  Point loc1 = loc + ppar1*axis;
  Point loc2  = loc + ppar2*axis;
#ifdef DEBUG_SMALL
  std::ofstream of4("axis_rotate.g2");
  of4 << "400 1 0 4 255 0 0 255" << std::endl;
  of4 << rotated.size() << std::endl;
  for (size_t kr=0; kr<rotated.size(); ++kr)
    of4 << rotated[kr] << std::endl;

  of4 << "410 1 0 0" << std::endl;
  of4 << "1" << std::endl;
  of4 << loc1 << " " << loc2 << std::endl;
#endif

  // Parameterize rotated points on axis
  shared_ptr<SplineCurve> line_cv(new SplineCurve(loc1, loc2));
  vector<double> pts;
  vector<double> param;
  vector<double> distance;
  double tmin = line_cv->startparam();
  double tmax = line_cv->endparam();
  double tmin2 = tmax;
  double tmax2 = tmin;
  for (size_t ki=0; ki<rotated.size(); ++ki)
    {
      double tpar, dist;
      Point close;
      line_cv->closestPoint(rotated[ki], tmin, tmax, tpar, close, dist);
      pts.insert(pts.end(), rotated[ki].begin(), rotated[ki].end());
      param.push_back(tpar);
      distance.push_back(dist);
      tmin2 = std::min(tmin2, tpar);
      tmax2 = std::max(tmax2, tpar);
    }

  if (tmax2 - tmin2 < eps)
    return false;

  // Statistics
  int num_del = 10;
  double tdel = (tmax2 - tmin2)/(double)(num_del);
  vector<int> num_pts(num_del, 0);
  vector<double> min_dist(num_del, std::numeric_limits<double>::max());
  vector<double> max_dist(num_del, std::numeric_limits<double>::lowest());
  vector<double> avdist(num_del, 0.0);
  for (size_t kj=0; kj<param.size(); ++kj)
    {
      int ka = (int)((param[kj]-tmin2)/tdel);
      if (ka >= (int)num_pts.size())
	ka = (int)num_pts.size() - 1;
      num_pts[ka]++;
      min_dist[ka] = std::min(min_dist[ka], distance[kj]);
      max_dist[ka] = std::max(max_dist[ka], distance[kj]);
      avdist[ka] += distance[kj];
    }

  vector<double> range(num_del);
  vector<double> var(num_del);
  vector<double> bddist(num_del);
  vector<double> rfrac(num_del);
  double avnum = 0.0, avrange = 0.0, avbdd = 0.0, avrfrac = 0.0;
  double fac = 1.0/(double)num_del;
  for (size_t kj=0; kj<avdist.size(); ++kj)
    {
      if (num_pts[kj] > 0)
	avdist[kj] /= (double)num_pts[kj];
      range[kj] = max_dist[kj] - min_dist[kj];
      var[kj] = 0.5*(min_dist[kj]+max_dist[kj]) - avdist[kj];
      bddist[kj] = std::min(avdist[kj]-min_dist[kj], max_dist[kj]-avdist[kj]);
      rfrac[kj] = range[kj]/(double)num_pts[kj];
      avnum += fac*(double)num_pts[kj];
      avrange += fac*range[kj];
      avbdd += fac*bddist[kj];
      avrfrac += fac*rfrac[kj];
    }
  double bdfrac = avbdd/avrange;

  double numfac = 1.5;
  double rfac = 1.5;
  double rffac = 5.0;  // Need another criterion

  // Shorten curve
  double tmin3 = tmin2, tmax3 = tmax2;
  for (int ka=0; ka<num_del; ++ka)
    {
      if (((double)num_pts[ka] > numfac*avnum || rfrac[ka] > rffac*avrfrac) &&
	  range[ka] > rfac*avrange && bddist[ka]/range[ka] > bdfrac)
	tmin3 += tdel;
      else
	break;
    }
  for (int ka=num_del-1; ka>=0; --ka)
    {
      if (((double)num_pts[ka] > numfac*avnum || rfrac[ka] > rffac*avrfrac) &&
	  range[ka] > rfac*avrange && bddist[ka]/range[ka] > bdfrac)
	tmax3 -= tdel;
      else
	break;
    }

  if (tmin3 >= tmax3)
    return found;

  vector<Point> rotated2;
  vector<double> param2;
  if (tmin3 > tmin2 || tmax3 < tmax2)
    {
      for (size_t ki=0; ki<param.size(); ++ki)
	if (param[ki] >= tmin3 && param[ki] <= tmax3)
	  {
	    rotated2.push_back(rotated[ki]);
	    param2.push_back(param[ki]);
	  }
    }
  else
    {
      rotated2 = rotated;
      param2 = param;
    }

  int in = 12;
  int ik = 3;
  shared_ptr<SplineCurve> approx_cv;
  RevEngUtils::curveApprox(rotated2, param2, ik, in, approx_cv);
  
#ifdef DEBUG_SMALL
  approx_cv->writeStandardHeader(of4);
  approx_cv->write(of4);
#endif

  shared_ptr<SplineCurve> approx_line0;
  RevEngUtils::curveApprox(rotated2, param2, 2, 2, approx_line0);
  
#ifdef DEBUG_SMALL
  approx_line0->writeStandardHeader(of4);
  approx_line0->write(of4);
#endif

  // Check accuracy of line
  vector<Point> der0(2);
  approx_line0->point(der0,
		     0.5*(approx_line0->startparam()+approx_line0->endparam()), 1);
  shared_ptr<Line> curr_line0(new Line(der0[0], der0[1]));

  double maxdist0, avdist0;
  int num_in0;
  vector<double> curr_dist0;
  RevEngUtils::distToCurve(rotated, curr_line0, approx_tol_, maxdist0,
			   avdist0, num_in0, curr_dist0);
  if (num_in0 > (int)rotated.size()/2 && avdist0 <= approx_tol_)
    {
      // Define surface
      double ang = der0[1].angle(axis);
      ang = std::min(ang, M_PI-ang);
      Point axis_pt = loc + ((der0[0]-loc)*axis)*axis;
      double rad = der0[0].dist(axis_pt);
      if (ang < anglim2)
	{
	  // Create cylinder
	  shared_ptr<Cylinder> cyl(new Cylinder(rad, axis_pt, axis, Cx));
	  sfs.push_back(cyl);
	}
      else
	{
	  // Create cone
	  // Sign of angle
	  Point pos3;
	  approx_line0->point(pos3, approx_line0->endparam());
	  Point axis_pt2 = loc + ((pos3-loc)*axis)*axis;
	  double rad2 = pos3.dist(axis_pt2);
	  int sgn = ((pos3 - der0[0])*axis < 0.0) ? -1 : 1;
	  if (sgn*rad2 < sgn*rad)
	    ang *= -1;
	  
	  shared_ptr<Cone> cone(new Cone(rad, axis_pt, axis, Cx, ang));
	  sfs.push_back(cone);
	}
      return true;
    }

  // Angle between consecutive control segments
  vector<double> seg_ang(in-2);
  vector<double> coefs(approx_cv->coefs_begin(), approx_cv->coefs_end());
  int dim = 3;
  Point pt1(coefs.begin(), coefs.begin()+dim);
  Point pt2(coefs.begin()+dim, coefs.begin()+2*dim);
  for (int ka=2; ka<in; ++ka)
    {
      Point pt3(coefs.begin()+ka*dim,coefs.begin()+(ka+1)*dim);
      Point vec1 = pt2 - pt1;
      Point vec2 = pt3 - pt2;
      seg_ang[ka-2] = vec1.angle(vec2);
      pt1 = pt2;
      pt2 = pt3;
    }

  // Potential lines
  vector<vector<Point> > line_coefs;
  vector<vector<double> > line_par;
  vector<Point> currc;
  vector<double> currp;
  for (int ka=0; ka<(int)seg_ang.size(); ++ka)
    {
      if (seg_ang[ka] > anglim2)
	{
	  if (currc.size() > 0)
	    {
	      line_coefs.push_back(currc);
	      line_par.push_back(currp);
	      currc.clear();
	      currp.clear();
	    }
	}
      else
	{
	  if (currc.size() == 0)
	    {
	      currc.push_back(Point(coefs.begin()+ka*dim,coefs.begin()+(ka+1)*dim));
	      currc.push_back(Point(coefs.begin()+(ka+1)*dim,coefs.begin()+(ka+2)*dim));
	      currp.push_back(approx_cv->basis().grevilleParameter(ka));
	      currp.push_back(approx_cv->basis().grevilleParameter(ka+1));
	    }
	  currc.push_back(Point(coefs.begin()+(ka+2)*dim,coefs.begin()+(ka+3)*dim));
	  currp.push_back(approx_cv->basis().grevilleParameter(ka+2));

	}
    }

  if (line_par.size() == 0)
    {
      if (currp.size() < 2 || currp[currp.size()-1] - currp[0] < eps)
	return false;
      vector<double> tmp_par(2);
      tmp_par[0] = currp[0];
      tmp_par[1] = currp[currp.size()-1];
      line_par.push_back(tmp_par);
    }
  
  for (size_t kj=0; kj<line_par.size(); ++kj)
    {
      double tmin4 = line_par[kj][0];
      double tmax4 = line_par[kj][line_par[kj].size()-1];
      vector<Point> rotated3;
      vector<double> param3;
      if (tmin4 > tmin3 || tmax4 < tmax3)
	{
	  for (size_t kr=0; kr<param2.size(); ++kr)
	    if (param2[kr] >= tmin4 && param2[kr] <= tmax4)
	      {
		rotated3.push_back(rotated2[kr]);
		param3.push_back(param2[kr]);
	      }
	}
      else
	{
	  rotated3 = rotated2;
	  param3 = param2;
	}

      if ((int)rotated3.size() < num_pt_lim)
	continue;
      
      shared_ptr<SplineCurve> approx_line;
      RevEngUtils::curveApprox(rotated3, param3, 2, 2, approx_line);
      
#ifdef DEBUG_SMALL
      approx_line->writeStandardHeader(of4);
      approx_line->write(of4);
#endif

      // Check accuracy for all points
      vector<Point> der(2);
      approx_line->point(der,
			 0.5*(approx_line->startparam()+approx_line->endparam()), 1);
      shared_ptr<Line> curr_line(new Line(der[0], der[1]));

      double maxdist, avdist;
      int num_in;
      vector<double> curr_dist;
      RevEngUtils::distToCurve(rotated, curr_line, approx_tol_, maxdist,
			       avdist, num_in, curr_dist);

      // double rad = 0.0;
      // for (size_t kr=0; kr<curr_dist.size(); ++kr)
      // 	{
      // 	  if (curr_dist[kr] <= approx_tol_)
      // 	    {
      // 	      Point tmp = loc + ((rotated[kr]-loc)*axis)*axis;
      // 	      double rad0 = rotated[kr].dist(tmp);
      // 	      rad += rad0;
      // 	    }
      // 	}
      // rad /= (double)num_in;
	  
#ifdef DEBUG_SMALL
      std::ofstream of1("in_out_curr.g2");
      vector<Point> in_pt, out_pt;
      for (size_t kr=0; kr<curr_dist.size(); ++kr)
	{
	  if (curr_dist[kr] <= approx_tol_)
	    in_pt.push_back(rotated[kr]);
	  else
	    out_pt.push_back(rotated[kr]);
	}
      if (in_pt.size() > 0)
	{
	  of1 << "400 1 0 4 255 0 0 255" << std::endl;
	  of1 << in_pt.size() << std::endl;
	  for (size_t kr=0; kr<in_pt.size(); ++kr)
	    of1 << in_pt[kr] << std::endl;
	}
      if (out_pt.size() > 0)
	{
	  of1 << "400 1 0 4 0 255 0 255" << std::endl;
	  of1 << out_pt.size() << std::endl;
	  for (size_t kr=0; kr<out_pt.size(); ++kr)
	    of1 << out_pt[kr] << std::endl;
	}
#endif

      // Remove most distant points and try again
      double dlim = 2.0*approx_tol_;
      vector<Point> rotated4;
      vector<double> param4;
      for (size_t kr=0; kr<curr_dist.size(); ++kr)
	{
	  if (curr_dist[kr] <= dlim)
	    {
	      rotated4.push_back(rotated[kr]);
	      param4.push_back(param[kr]);
	    }
	}

      shared_ptr<SplineCurve> approx_line2;
      RevEngUtils::curveApprox(rotated4, param4, 2, 2, approx_line2);
      
#ifdef DEBUG_SMALL
      std::ofstream of2("in_out2_curr.g2");
      approx_line2->writeStandardHeader(of2);
      approx_line2->write(of2);
#endif

      // Check accuracy for all points
      vector<Point> der2(2);
      approx_line2->point(der2,
			 0.5*(approx_line2->startparam()+approx_line2->endparam()), 1);
      shared_ptr<Line> curr_line2(new Line(der2[0], der2[1]));

      double maxdist2, avdist2;
      int num_in2;
      vector<double> curr_dist2;
      RevEngUtils::distToCurve(rotated, curr_line2, approx_tol_, maxdist2,
			       avdist2, num_in2, curr_dist2);
     // double rad2 = 0.0;
     //  for (size_t kr=0; kr<curr_dist.size(); ++kr)
     // 	{
     // 	  if (curr_dist2[kr] <= approx_tol_)
     // 	    {
     // 	      Point tmp = loc + ((rotated[kr]-loc)*axis)*axis;
     // 	      double rad0 = rotated[kr].dist(tmp);
     // 	      rad2 += rad0;
     // 	    }
     // 	}
     //  rad2 /= (double)num_in2;
      
#ifdef DEBUG_SMALL
      vector<Point> in_pt2, out_pt2;
      for (size_t kr=0; kr<curr_dist2.size(); ++kr)
	{
	  if (curr_dist2[kr] <= approx_tol_)
	    in_pt2.push_back(rotated[kr]);
	  else
	    out_pt2.push_back(rotated[kr]);
	}
      if (in_pt.size() > 0)
	{
	  of2 << "400 1 0 4 255 0 0 255" << std::endl;
	  of2 << in_pt2.size() << std::endl;
	  for (size_t kr=0; kr<in_pt2.size(); ++kr)
	    of2 << in_pt2[kr] << std::endl;
	}
      if (out_pt2.size() > 0)
	{
	  of2 << "400 1 0 4 0 255 0 255" << std::endl;
	  of2 << out_pt2.size() << std::endl;
	  for (size_t kr=0; kr<out_pt2.size(); ++kr)
	    of2 << out_pt2[kr] << std::endl;
	}
#endif

      // Define surface
      if ((double)(std::max(num_in0, std::max(num_in, num_in2))) <
	  num_fac*(double)rotated2.size())
	continue;
      if (num_in0 > std::max(num_in, num_in2))
	{
	  //std::swap(rad, rad2);
	  std::swap(der0[0], der2[0]);
	  std::swap(der0[1], der2[1]);
	}
      else if (num_in > num_in2)
	{
	  //std::swap(rad, rad2);
	  std::swap(der[0], der2[0]);
	  std::swap(der[1], der2[1]);
	}
      if (std::max(num_in0, std::max(num_in, num_in2)) > num_pt_lim)
	{
	  double ang = der2[1].angle(axis);
	  ang = std::min(ang, M_PI-ang);
	  Point axis_pt = loc + ((der2[0]-loc)*axis)*axis;
	  double rad = der2[0].dist(axis_pt);
	  if (ang < anglim2)
	    {
	      // Create cylinder
	      shared_ptr<Cylinder> cyl(new Cylinder(rad, axis_pt, axis, Cx));
	      sfs.push_back(cyl);
	    }
	  else
	    {
	      // Create cone
	      // Sign of angle
	      Point pos3;
	      if (num_in2 > num_in) 
		approx_line2->point(pos3, approx_line2->endparam());
	      else
		approx_line->point(pos3, approx_line->endparam());
	      Point axis_pt2 = loc + ((pos3-loc)*axis)*axis;
	      double rad2 = pos3.dist(axis_pt2);
	      int sgn = ((pos3 - der[0])*axis < 0.0) ? -1 : 1;
	      if (sgn*rad2 < sgn*rad)
		ang *= -1;
	      
	      shared_ptr<Cone> cone(new Cone(rad, axis_pt, axis, Cx, ang));
	      sfs.push_back(cone);
	    }
	  found = true;
	}
    }

  return found;
}

//===========================================================================
bool RevEng::identifySmallPlanar(vector<RevEngPoint*>& points,
				 Point loc, Point axis, Point Cx,
				 double ppar1, double ppar2, double delta,
				 vector<shared_ptr<ElementarySurface> >& sfs)
//===========================================================================
{
  // Sort groups with respect to axis
  vector<Point> axis_pt(points.size());
  vector<double> axis_par(points.size());
  double tmin = std::numeric_limits<double>::max();
  double tmax = std::numeric_limits<double>::lowest();
  double avpar = 0.0;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      Vector3D xyz = points[ki]->getPoint();
      Point pos(xyz[0], xyz[1], xyz[2]);
      double tpar = (pos - loc)*axis;
      tmin = std::min(tmin, tpar);
      tmax = std::max(tmax, tpar);
      axis_par[ki] = tpar;
      axis_pt[ki] = loc + axis_par[ki]*axis;
    }

  size_t ki;
  double tlim = 1.5*approx_tol_;
  int nlim = 100;  // Currently
  double angtol = 5.0*anglim_;
  bool found = false;
  int num = 0;
  Point mid(0.0, 0.0, 0.0);
  for (ki=0; ki<axis_par.size(); ++ki)
    {
      if (axis_par[ki] <= ppar1+delta || axis_par[ki] >= ppar2-delta)
	continue;

      Point vec = axis_pt[ki] - loc;
      Point pos2 = loc + (vec*axis)*axis;
      mid += pos2;
      ++num;
    }

  if (num == 0)
    return false;
  mid /= (double)num;
      
  shared_ptr<Plane> plane(new Plane(mid, axis, Cx));
#ifdef DEBUG_SMALL
  std::ofstream of("small_plane.g2");
  plane->writeStandardHeader(of);
  plane->write(of);
#endif

  // Check accuracy
  double maxd, avd;
  int n_in, n2_in;
  vector<double> parvals;
  vector<pair<double,double> > distang;
  vector<RevEngPoint*> in, out;
  RevEngUtils::distToSurf(points, plane, approx_tol_, angtol,
			  maxd, avd, n_in, n2_in, distang);

  if (avd <= approx_tol_ && n2_in > points.size()/2)
    {
      sfs.push_back(plane);
      found = true;
    }
  return found;
}
 
//===========================================================================
bool RevEng::identifySmallPlanar(vector<RevEngRegion*>& groups,
				 Point loc, Point axis, Point Cx,
				 double ppar1, double ppar2, double delta,
				 vector<shared_ptr<ElementarySurface> >& sfs)
//===========================================================================
{
  // Sort groups with respect to axis
  vector<Point> axis_pt(groups.size());
  vector<double> axis_par(groups.size());
  vector<pair<double,double> > axis_range(groups.size());
  for (size_t ki=0; ki<groups.size(); ++ki)
    {
      int num = groups[ki]->numPoints();
      double tmin = std::numeric_limits<double>::max();
      double tmax = std::numeric_limits<double>::lowest();
      double avpar = 0.0;
      double fac = 1.0/(double)num;
      for (int ka=0; ka<num; ++ka)
	{
	  Vector3D xyz = groups[ki]->getPoint(ka)->getPoint();
	  Point pos(xyz[0], xyz[1], xyz[2]);
	  double tpar = (pos - loc)*axis;
	  tmin = std::min(tmin, tpar);
	  tmax = std::max(tmax, tpar);
	  avpar += fac*tpar;
	}
      axis_par[ki] = avpar;
      axis_range[ki] = std::make_pair(tmin,tmax);
      axis_pt[ki] = loc + axis_par[ki]*axis;
    }

  for (size_t ki=0; ki<groups.size(); ++ki)
    for (size_t kj=ki+1; kj<groups.size(); ++kj)
      {
	if (axis_par[kj] < axis_par[ki])
	  {
	    std::swap(axis_par[ki], axis_par[kj]);
	    std::swap(axis_pt[ki], axis_pt[kj]);
	    std::swap(axis_range[ki], axis_range[kj]);
	    std::swap(groups[ki], groups[kj]);
	  }
      }

  size_t ki, kj;
  double tlim = 1.5*approx_tol_;
  int nlim = 100;  // Currently
  double angtol = 5.0*anglim_;
  bool found = false;
  for (ki=0; ki<axis_par.size(); ki=kj)
    {
      for (kj=ki+1; kj<axis_par.size() && axis_par[kj]-axis_par[ki]<tlim; ++kj);
      int num = 0;
      for (size_t kr=ki; kr<kj; ++kr)
	num += groups[kr]->numPoints();
      if (num < nlim)
	continue;
      if (axis_range[ki].second <= ppar1+delta)
	continue;
      if (axis_range[ki].first >= ppar2-delta)
	continue;

      Point mid(0.0, 0.0, 0.0);
      int nmb = 0;
      for (size_t kr=ki; kr<kj; ++kr)
	{
	  int num2 = groups[kr]->numPoints();
	  for (int ka=0; ka<num2; ++ka)
	    {
	      Vector3D pos0 = groups[kr]->getPoint(ka)->getPoint();
	      Point pos(pos0[0], pos0[1], pos0[2]);
	      double tpar = (pos - loc)*axis;
	      if (tpar <= ppar1+delta || tpar >= ppar2-delta)
		continue;
	      Point vec = pos - loc;
	      Point pos2 = loc + (vec*axis)*axis;
	      mid += pos2;
	      ++nmb;
	    }
	}

      if (nmb == 0)
	continue;
      mid /= (double)nmb;
      
      shared_ptr<Plane> plane(new Plane(mid, axis, Cx));
#ifdef DEBUG_SMALL
      std::ofstream of("small_plane.g2");
      plane->writeStandardHeader(of);
      plane->write(of);
#endif

      // Check accuracy
      double maxdist = 0.0, avdist = 0.0;
      int num_in = 0, num2_in = 0;
      double fac = 1.0/(double)(kj-ki);
      for (size_t kr=ki; kr<kj; ++kr)
	{
	  double maxd, avd;
	  int n_in, n2_in;
	  vector<double> parvals;
	  vector<pair<double,double> > distang;
	  vector<RevEngPoint*> in, out;
	  RevEngUtils::distToSurf(groups[kr]->pointsBegin(), groups[kr]->pointsEnd(),
				  plane, approx_tol_, maxd, avd, n_in,
				  n2_in, in, out, parvals, distang, angtol);
	  maxdist = std::max(maxdist, maxd);
	  avdist += fac*avd;
	  num_in += n_in;
	  num2_in += n2_in;
	}
      int surf_flag = groups[ki]->defineSfFlag(num, 0, approx_tol_, num_in,
					       num2_in, avdist, false);
      if (surf_flag < NOT_SET)
	{
	  sfs.push_back(plane);
	  found = true;
	}
    }
  return found;
}
 

//===========================================================================
void RevEng::adaptToMainAxis()
//===========================================================================
{
  doAdaptToAxis();

#ifdef DEBUG
  std::cout << "Post merge similar. Number of regions: " << regions_.size() << std::endl;

  checkConsistence("Regions11_1");

  if (regions_.size() > 0)
    {
      std::cout << "Regions 11_1" << std::endl;
      std::ofstream of("regions11_1.g2");
      std::ofstream ofm("mid_regions11_1.g2");
      std::ofstream ofs("small_regions11_1.g2");
      writeRegionStage(of, ofm, ofs);
      std::ofstream of0("regions11_1_helix.g2");
      for (int ki=0; ki<(int)regions_.size(); ++ki)
	{
	  if (regions_[ki]->getSurfaceFlag() == PROBABLE_HELIX)
	    {
	      regions_[ki]->writeRegionInfo(of0);
	      if (regions_[ki]->hasSurface())
		regions_[ki]->writeSurface(of0);
	    }
	}
      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf11_1.g2");
	  writeRegionWithSurf(of);
	}
    }
#endif
	  
#ifdef DEBUG
  std::ofstream ofu1("unresolved2.g2");
  for (size_t kr=0; kr<regions_.size(); ++kr)
    {
      if (regions_[kr]->hasSurface())
	continue;
      if (regions_[kr]->hasAssociatedBlend())
	continue;
      regions_[kr]->writeRegionPoints(ofu1);
    }
#endif

  // Merge adjacent edges if possible
  double eps = 1.0e-9;
  for (int ka=0; ka<(int)edges_.size(); ++ka)
    {
      RevEngRegion* adj1[2];
      edges_[ka]->getAdjacent(adj1[0], adj1[1]);
      for (int kb=ka+1; kb<(int)edges_.size(); ++kb)
	{
	  RevEngRegion* adj2[2];
	  edges_[kb]->getAdjacent(adj2[0], adj2[1]);

	  if ((adj1[0] == adj2[0] || adj1[0] == adj2[1]) &&
	      (adj1[1] == adj2[0] || adj1[1] == adj2[1]))
	    {
	      double par1, par2;
	      bool adjacent = edges_[ka]->isAdjacent(edges_[kb].get(),
						     approx_tol_, par1, par2);
	      if (adjacent)
		{
#ifdef DEBUG_BLEND
		  std::ofstream of("adj_blend.g2");
		  adj1[0]->writeRegionPoints(of);
		  adj1[1]->writeRegionPoints(of);
		  vector<shared_ptr<ParamCurve> > cvs1 = edges_[ka]->getSpaceCurves();
		  for (size_t kh=0; kh<cvs1.size(); ++kh)
		    {
		      cvs1[kh]->writeStandardHeader(of);
		      cvs1[kh]->write(of);
		    }
		  vector<shared_ptr<ParamCurve> > cvs2 = edges_[kb]->getSpaceCurves();
		  for (size_t kh=0; kh<cvs2.size(); ++kh)
		    {
		      cvs2[kh]->writeStandardHeader(of);
		      cvs2[kh]->write(of);
		    }
	      
#endif
		  // Check for seam of adjacent surface
		  shared_ptr<ParamSurface> surf1 = adj1[0]->getSurface(0)->surface();
		  shared_ptr<ParamSurface> surf2 = adj1[1]->getSurface(0)->surface();
		  shared_ptr<ElementarySurface> elem1 =
		    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf1);
		  shared_ptr<ElementarySurface> elem2 =
		    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf2);
		  bool close1=false, close2=false, close3=false, close4=false;
		  if (elem1.get())
		    elem1->isClosed(close1, close2);
		  if (elem2.get())
		    elem2->isClosed(close3, close4);
		  Point pos = edges_[ka]->point(par1);
		  if (close1 || close2)
		    {
		      RectDomain dom = elem1->containingDomain();
		      double upar, vpar, dist;
		      Point close_pt;
		      elem1->closestBoundaryPoint(pos, upar, vpar, close_pt, dist, eps);
		      if (dist < approx_tol_)
			{
			  if (close1 && (fabs(upar-dom.umin()) < eps || fabs(dom.umax()-upar) < eps))
			    adjacent = false;
			  if (close2 && (fabs(vpar-dom.vmin()) < eps || fabs(dom.vmax()-vpar) < eps))
			    adjacent = false;
			}
		    }
		  
		  if (close3 || close4)
		    {
		      RectDomain dom = elem2->containingDomain();
		      double upar, vpar, dist;
		      Point close_pt;
		      elem2->closestBoundaryPoint(pos, upar, vpar, close_pt, dist, eps);
		      if (dist < approx_tol_)
			{
			  if (close3 && (fabs(upar-dom.umin()) < eps || fabs(dom.umax()-upar) < eps))
			    adjacent = false;
			  if (close4 && (fabs(vpar-dom.vmin()) < eps || fabs(dom.vmax()-vpar) < eps))
			    adjacent = false;
			}
		    }
		}

	      if (adjacent)
		{
		  double t1 = edges_[ka]->startparam();
		  double t2 = edges_[ka]->endparam();
		  if (fabs(par1-t1) < fabs(par2-t2))
		    std::swap(edges_[ka], edges_[kb]);
		  bool done = edges_[ka]->append(edges_[kb].get(), approx_tol_);
		  if (done)
		    {
		      edges_.erase(edges_.begin()+kb);
		      kb--;
		    }
		  int stop_break0 = 1;
		}
	    }
	}
    }

   int stop_break = 1;
  
}

//===========================================================================
void RevEng::doAdaptToAxis()
//===========================================================================
{
  vector<SurfaceProperties> sfprop;
  collectAxis(sfprop);

  // Sort surfaces into groups with roughly the same axis
  double epsang = 0.05;
  double eps = 1.0e-4;
  vector<DirectionCone> axis_cone;
  vector<vector<size_t> > group_ixs;
  vector<int> num_pts;
  for (size_t ki=0; ki<sfprop.size(); ++ki)
    {
      if (sfprop[ki].type_ == Class_Sphere)
	continue;   // Sphere centre is relevant, axis is derived
      size_t kr;
      for (kr=0; kr<axis_cone.size(); ++kr)
	{
	  Point centre = axis_cone[kr].centre();
	  if (sfprop[ki].dir_*centre < 0.0)
	    sfprop[ki].dir_ *= -1;
	  Point axis = sfprop[ki].dir_;
	  double angle = axis_cone[kr].unionAngle(axis);
	  if (angle <= epsang)
	    {
	      axis_cone[kr].addUnionWith(axis);
	      group_ixs[kr].push_back(ki);
	      num_pts[kr] += sfprop[ki].num_points_;
	      break;
	    }
	}
      if (kr == axis_cone.size())
	{
	  DirectionCone cone(sfprop[ki].dir_);
	  axis_cone.push_back(cone);
	  vector<size_t> ixs;
	  ixs.push_back(ki);
	  group_ixs.push_back(ixs);
	  num_pts.push_back(sfprop[ki].num_points_);
	}
    }

  // Sort according to number of points associated to the surface
  for (size_t ki=0; ki<num_pts.size(); ++ki)
    for (size_t kj=ki+1; kj<num_pts.size(); ++kj)
      {
	if (num_pts[kj] > num_pts[ki])
	  {
	    std::swap(num_pts[ki], num_pts[kj]);
	    std::swap(axis_cone[ki], axis_cone[kj]);
	    std::swap(group_ixs[ki], group_ixs[kj]);
	  }
      }

  // For each group, ecompute axis direction based on reliability of surface 
  // type and number of points associated to the surface
  double planefac = 1.0;
  double cylfac = 0.75;
  double conefac = 0.6;
  double torfac = 0.4;
  vector<Point> axes(group_ixs.size());
  for (size_t ki=0; ki<group_ixs.size(); ++ki)
    {
      Point curr(0.0, 0.0, 0.0);
      for (size_t kj=0; kj<group_ixs[ki].size(); ++kj)
	{
	  Point curr2 = sfprop[group_ixs[ki][kj]].dir_;
	  double fac = planefac;
	  if (sfprop[group_ixs[ki][kj]].type_ == Class_Cylinder)
	    fac = cylfac;
	  else if (sfprop[group_ixs[ki][kj]].type_ == Class_Cone)
	    fac = conefac;
	  else if (sfprop[group_ixs[ki][kj]].type_ == Class_Torus)
	    fac = torfac;
	  int num = sfprop[group_ixs[ki][kj]].num_points_;
	  curr += (double)num*fac*curr2;
	}
      if (curr.length() < eps)
	axes[ki] = axis_cone[ki].centre();
      else
	{
	  curr.normalize();
	  axes[ki] = curr;
	}
    }

  adjustWithMainAxis(axes, num_pts); 

  size_t num_axis = model_axis_.size();
  vector<RevEngEdge*> nopar_edgs;
  for (size_t ki=0; ki<num_pts.size(); ++ki)
    {
      int ix=-1;
      double min_ang = M_PI;
      for (int ka=0; ka<3; ++ka)
	{
	  double ang = axes[ki].angle(mainaxis_[ka]);
	  ang = std::min(ang, M_PI-ang);
	  if (ang < min_ang)
	    {
	      min_ang = ang;
	      ix = ka;
	    }
	}

      Point Cx = mainaxis_[(ix+2)%3].cross(axes[ki]);
      Cx.normalize();

      AxisInfo curr_axis(axes[ki]);
      model_axis_.push_back(curr_axis);
      
     // Identify surfaces with the "same" axis including location
      vector<vector<int> > plane_adapt;
      vector<int> num_adapt1;
      for (size_t kj=0; kj<group_ixs[ki].size(); ++kj)
	{
	  if (sfprop[group_ixs[ki][kj]].type_  != Class_Plane)
	    continue;

	  shared_ptr<ParamSurface> surf = surfaces_[sfprop[group_ixs[ki][kj]].sfix_]->surface();
	  shared_ptr<Plane> plane = dynamic_pointer_cast<Plane,ParamSurface>(surf);
	  Point loc = plane->location();

	  size_t kr, kh;
	  for (kr=0; kr<plane_adapt.size(); ++kr)
	    {
	      for (kh=0; kh<plane_adapt[kr].size(); ++kh)
		{
		  shared_ptr<ParamSurface> surf2 = surfaces_[plane_adapt[kr][kh]]->surface();
		  shared_ptr<Plane> plane2 = dynamic_pointer_cast<Plane,ParamSurface>(surf2);
		  Point loc2 = plane2->location();
		  double dd = fabs((loc-loc2)*axes[ki]);
		  if (dd <= approx_tol_)
		    break;
		}
	      if (kh < plane_adapt[kr].size())
		break;
	    }
	  if (kr == plane_adapt.size())
	    {
	      vector<int> curr_adapt;
	      curr_adapt.push_back(sfprop[group_ixs[ki][kj]].sfix_);
	      plane_adapt.push_back(curr_adapt);
	      num_adapt1.push_back(sfprop[group_ixs[ki][kj]].num_points_);
	    }
	  else
	    {
	      plane_adapt[kr].push_back(sfprop[group_ixs[ki][kj]].sfix_);
	      num_adapt1[kr] += sfprop[group_ixs[ki][kj]].num_points_;
	    }
	}

      for (size_t kj=0; kj<plane_adapt.size(); ++kj)
	{
	  Point location = planarFit(plane_adapt[kj], axes[ki]);
	  model_axis_[num_axis+ki].addPlaneLocation(location, num_adapt1[kj]);
	}
      
      // The other surfaces
      vector<vector<int> > rotational_adapt;
       vector<int> num_adapt2;
      for (size_t kj=0; kj<group_ixs[ki].size(); ++kj)
	{
	  if (sfprop[group_ixs[ki][kj]].type_ == Class_Plane)
	    continue;

	  Point loc = sfprop[group_ixs[ki][kj]].loc_;
	  size_t kr, kh;
	  for (kr=0; kr<rotational_adapt.size(); ++kr)
	    {
	      for (kh=0; kh<rotational_adapt[kr].size(); ++kh)
		{
		  shared_ptr<ParamSurface> surf2 =
		    surfaces_[rotational_adapt[kr][kh]]->surface();
		  shared_ptr<ElementarySurface> elem =
		    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf2);
		  Point loc2 = elem->location();
		  Point tmp = loc + ((loc2-loc)*axes[ki])*axes[ki];
		  double dd = loc2.dist(tmp);
		  if (dd <= approx_tol_)
		    break;
		}
	      if (kh < rotational_adapt[kr].size())
		break;
	    }
	  if (kr == rotational_adapt.size())
	    {
	      vector<int> curr_adapt;
	      curr_adapt.push_back(sfprop[group_ixs[ki][kj]].sfix_);
	      rotational_adapt.push_back(curr_adapt);
	      num_adapt2.push_back(sfprop[group_ixs[ki][kj]].num_points_);
	    }
	  else
	    {
	      rotational_adapt[kr].push_back(sfprop[group_ixs[ki][kj]].sfix_);
	      num_adapt2[kr] += sfprop[group_ixs[ki][kj]].num_points_;
	    }
	}

      for (size_t kj=0; kj<rotational_adapt.size(); ++kj)
	{
	  Point location = rotationalFit(rotational_adapt[kj], axes[ki], Cx,
					 nopar_edgs);
	  if (location.dimension() > 0)
	    model_axis_[num_axis+ki].addRotationalLocation(location, num_adapt2[kj]);
	}
    }

  // Sort axes
  vector<int> num_axis_pt(model_axis_.size());
  for (size_t ki=0; ki<model_axis_.size(); ++ki)
    {
      int num_pt = 0;
      for (size_t kj=0; kj<model_axis_[ki].plane_loc_.size(); ++kj)
	num_pt += model_axis_[ki].plane_loc_[kj].second;
      for (size_t kj=0; kj<model_axis_[ki].rotational_loc_.size(); ++kj)
	num_pt += model_axis_[ki].rotational_loc_[kj].second;
      num_axis_pt[ki] = num_pt;
    }
  
   for (size_t ki=0; ki<model_axis_.size(); ++ki)
    {
      for (size_t kj=ki+1; kj<model_axis_.size(); ++kj)
	if (num_axis_pt[ki] < num_axis_pt[kj])
	  {
	    std::swap(model_axis_[ki], model_axis_[kj]);
	    std::swap(num_axis_pt[ki], num_axis_pt[kj]);
	  }
    }
 

  if (nopar_edgs.size() > 0)
    {
      // Must split edge at seam
      vector<shared_ptr<RevEngEdge> > added_edgs;
      vector<shared_ptr<RevEngRegion> > added_regs;
      vector<shared_ptr<HedgeSurface> > added_sfs;
      for (size_t ki=0; ki<nopar_edgs.size(); ++ki)
	{
	  if (nopar_edgs[ki]->isClosed(approx_tol_))
	    continue;  // Will be handled later
	  nopar_edgs[ki]->splitAtSeam(approx_tol_, added_edgs, added_regs,
				      added_sfs);
	}
      if (added_edgs.size() > 0)
	edges_.insert(edges_.end(), added_edgs.begin(), added_edgs.end());
      if (added_regs.size() > 0)
	regions_.insert(regions_.end(), added_regs.begin(), added_regs.end());
      if (added_sfs.size() > 0)
	surfaces_.insert(surfaces_.end(), added_sfs.begin(), added_sfs.end());
    }
  
  int stop_break = 1;
}

//===========================================================================
void RevEng::adjustWithMainAxis(vector<Point>& axes, vector<int>& num_pts)
//===========================================================================
{
  double angtol1 = 5.0*anglim_;
  double angtol2 = 5.0*angtol1;
  double pihalf = 0.5*M_PI;
  double eps = 1.0e-6;

  // Change sign of axis if in conflict with existing main axis
  for (size_t ki=0; ki<axes.size(); ++ki)
    {
      for (int ka=0; ka<3; ++ka)
	{
	  double ang = axes[ki].angle(mainaxis_[ka]);
	  ang = std::min(ang, M_PI-ang);
	  double scpr = axes[ki]*mainaxis_[ka];
	  if (scpr < 0.0 && ang<angtol2)
	    axes[ki] *= -1.0;
	}
    }

  vector<Point> axes2(axes.size());
  for (size_t ki=0; ki<axes.size(); ++ki)
    axes2[ki] = axes[ki];
  
  for (int ix1=(int)axes2.size()-1; ix1>=0; --ix1)
    {
      Vector3D axis1(axes2[ix1][0], axes2[ix1][1], axes2[ix1][2]);
      for (int ix2=ix1-1; ix2>=0; --ix2)
	{
	  Vector3D axis2(axes2[ix2][0], axes2[ix2][1], axes2[ix2][2]);
	  double ang = axes2[ix1].angle(axes2[ix2]);
	  double ang2 = fabs(pihalf-ang);
	  if (ang2 > eps && ang2 < angtol1)
	    {
	      // Adjust axes2 to achieve orthogonality
	      Point vec = axes2[ix1].cross(axes2[ix2]);
	      int sgn = (ang < pihalf) ? 1 : -1;
	      double fac = 1.0/(double)(num_pts[ix1] + num_pts[ix2]);
	      double ang3 = (double)num_pts[ix2]*fac*ang2;
	      double ang4 = (double)num_pts[ix1]*fac*ang2;
	      Matrix3D mat1, mat2;
	      mat1.setToRotation(-sgn*ang3, vec[0], vec[1], vec[2]);
	      mat2.setToRotation(sgn*ang4, vec[0], vec[1], vec[2]);
	      Vector3D axis3 = mat1*axis1;
	      Vector3D axis4 = mat2*axis2;
	      axis3.normalize_checked();
	      axis4.normalize_checked();
	      axes2[ix1] = Point(axis3[0], axis3[1], axis3[2]);
	      axes2[ix2] = Point(axis4[0], axis4[1], axis4[2]);
	      int stop_break = 1;
	    }
	}
    }

  // Ensure orthogonality
  int ix1 = 0;
  int ix2 = -1, ix3 = -1;
  for (int ka=1; ka<(int)axes2.size(); ++ka)
    {
      double ang = axes2[ix1].angle(axes2[ka]);
      if (fabs(pihalf-ang) < angtol1)
	{
	  if (ix2 < 0)
	    ix2 = ka;
	  else if (ix3 < 0)
	    {
	      ix3 = ka;
	      break;
	    }
	}
    }

  if (ix2 < 0)
    {
      ix2 = (int)axes2.size();
      ix3 = (int)axes2.size()+1;
      axes2.resize(axes2.size()+2);
    }
  else if (ix3 < 0)
    {
      ix3 = (int)axes2.size();
      axes2.resize(axes2.size()+1);
    }

  if (ix3 > (int)axes.size())
    {
      int min_ix = 0;
      double min_ang = fabs(axes2[ix1].angle(mainaxis_[0]));
      for (int kb=1; kb<3; ++kb)
	{
	  double ang = fabs(axes2[ix1].angle(mainaxis_[kb]));
	  if (ang < min_ang)
	    {
	      min_ix = kb;
	      min_ang = ang;
	    }
	}

      int min_ix2 = -1;
      if (ix2 < (int)axes.size())
	{
	  min_ix2 = (min_ang == 0) ? 1 : 0;
	  double min_ang2 = fabs(axes2[ix2].angle(mainaxis_[min_ix2]));
	  for (int kb=0; kb<3; ++kb)
	    {
	      if (kb == min_ix || kb == min_ix2)
		continue;
	      double ang = fabs(axes2[ix2].angle(mainaxis_[kb]));
	      if (ang < min_ang)
		{
		  min_ix2 = kb;
		  min_ang2 = ang;
		}
	    }
	}
      else
	{
	  for (int kb=0; kb<3; ++kb)
	    if (kb != min_ix)
	      {
		axes2[ix2] = mainaxis_[kb];
		min_ix2 = kb;
		break;
	      }
	}

      for (int kb=0; kb<3; ++kb)
	if (kb != min_ix && kb != min_ix2)
	  {
	    axes2[ix3] = mainaxis_[kb];
	    break;
	  }
    }

  axes2[ix3] = axes2[ix1].cross(axes2[ix2]);
  axes2[ix2] = axes2[ix3].cross(axes2[ix1]);

  for (size_t ki=0; ki<axes.size(); ++ki)
    axes[ki] = axes2[ki];

  mainaxis_[0] = axes2[ix1];
  mainaxis_[1] = axes2[ix2];
  mainaxis_[2] = axes2[ix3];
}

//===========================================================================
Point RevEng::planarFit(vector<int>& sf_ix, Point axis)
//===========================================================================
{
  double angtol = 5.0*anglim_;
  
  // Collect data points
  vector<RevEngPoint*> points;
  Point loc(0.0, 0.0, 0.0);
  int num = 0;
  for (size_t ki=0; ki<sf_ix.size(); ++ki)
    {
      HedgeSurface* surf = surfaces_[sf_ix[ki]].get();
      shared_ptr<Plane> plane =
	dynamic_pointer_cast<Plane,ParamSurface>(surf->surface());
      if (!plane.get())
	continue;
      loc += plane->location();
      num++;
      vector<RevEngRegion*> reg = surf->getRegions();
      for (size_t kj=0; kj<reg.size(); ++kj)
	{
	  points.insert(points.end(), reg[kj]->pointsBegin(),
			reg[kj]->pointsEnd());
	}
    }
  if (num > 0)
    loc /= (double)num;
  
  shared_ptr<Plane> plane2 = RevEngUtils::planeWithAxis(points, axis,
							loc, mainaxis_);

  Point centre = plane2->location();
  Point Cx = plane2->direction2();

  for (size_t ki=0; ki<sf_ix.size(); ++ki)
    {
      HedgeSurface* surf = surfaces_[sf_ix[ki]].get();
      shared_ptr<Plane> plane =
	dynamic_pointer_cast<Plane,ParamSurface>(surf->surface());
      if (!plane.get())
	continue;
      Point loc2 = plane->location();
      loc2 -= ((loc2-centre)*axis)*axis;
      shared_ptr<Plane> plane3(new Plane(loc2, axis, Cx));
      if (plane->direction()*axis < 0.0)
	plane3->swapParameterDirection();

      vector<RevEngRegion*> reg = surf->getRegions();
      double maxdist = 0.0, avdist = 0.0;
      int num_in = 0, num2_in = 0, num_pt = 0;
      double maxdist_init = 0.0, avdist_init = 0.0;
      int num_in_init = 0, num2_in_init = 0;
      vector<vector<double> > parvals(reg.size());
      vector<vector<pair<double,double> > > dist_ang(reg.size());
      vector<int> init_flag(reg.size());;
      for (size_t kj=0; kj<reg.size(); ++kj)
	{
	  double maxd0, avd0;
	  int num_in0, num2_in0;
	  vector<RevEngPoint*> in0, out0;
	  vector<pair<double,double> > distang0;
	  vector<double> parvals0;
	  RevEngUtils::distToSurf(reg[kj]->pointsBegin(), reg[kj]->pointsEnd(),
				  plane3, approx_tol_, maxd0, avd0, num_in0, num2_in0,
				  in0, out0, parvals0, distang0, angtol);
	  int num0 = reg[kj]->numPoints();
	  maxdist = std::max(maxdist, maxd0);
	  avdist = (double)num_pt*avdist + (double)num0*avd0;
	  num_in += num_in0;
	  num2_in += num2_in0;
	  parvals[kj] = parvals0;
	  dist_ang[kj] = distang0;

	  double maxd1, avd1;
	  int num_in1, num2_in1;
	  reg[kj]->getAccuracy(maxd1, avd1, num_in1, num2_in1);
	  maxdist_init = std::max(maxdist_init, maxd1);
	  avdist_init = (double)num_pt*avdist_init + (double)num0*avd1;
	  num_in_init += num_in1;
	  num2_in_init += num2_in1;
	  num_pt += num0;
	  init_flag[kj] = reg[kj]->getSurfaceFlag();
	}
      avdist /= (double)num_pt;
      avdist_init /= (double)num_pt;
      int surf_flag = reg[0]->defineSfFlag(num_pt, 0, approx_tol_, num_in,
					   num2_in, avdist, false);
      int surf_flag_init = reg[0]->defineSfFlag(num_pt, 0, approx_tol_, num_in_init,
						num2_in_init, avdist_init, false);

      double av_fac = 1.1;
      double in_fac = 0.9;
      if (surf_flag == 0 || surf_flag <= surf_flag_init ||
	  (avdist < av_fac*avdist_init &&
	   (double)num2_in > in_fac*(double)num2_in_init &&
	   (double)num_in > in_fac*(double)num_in_init))
	{
	  // Replace surface
	  for (size_t kj=0; kj<reg.size(); ++kj)
	    {
	      vector<RevEngEdge*> dummy;
	      reg[kj]->updateSurfaceAndInfo(plane3, approx_tol_, angtol,
					    parvals[kj], dist_ang[kj], dummy);
	    }
	}
      
      int stop_break = 1;
    }
  return centre;
}

//===========================================================================
Point RevEng::rotationalFit(vector<int>& sf_ix, Point axis, Point Cx,
			    vector<RevEngEdge*>& nopar_edgs)
//===========================================================================
{
  double angtol = 5.0*anglim_;
  Point dummy_loc;
  
  // Compute point on axis
  Point loc(0.0, 0.0, 0.0);
  int num = 0;
  for (size_t ki=0; ki<sf_ix.size(); ++ki)
    {
      HedgeSurface* surf = surfaces_[sf_ix[ki]].get();
      shared_ptr<ElementarySurface> elem =
	dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf->surface());
      if (!elem.get())
	continue;
      loc += elem->location();
      num++;
    }
  if (num > 0)
    loc /= (double)num;

  // Add spheres with centre on axis
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      shared_ptr<ParamSurface> sf = surfaces_[ki]->surface();
      if (sf->instanceType() == Class_Sphere)
	{
	  shared_ptr<Sphere> sphere = dynamic_pointer_cast<Sphere,ParamSurface>(sf);
	  Point centre = sphere->location();
	  Point tmp = loc + ((centre-loc)*axis)*axis;
	  if (centre.dist(tmp) < approx_tol_)
	    sf_ix.push_back((int)ki);
	}
    }

  double udomain[2];
  udomain[0] = std::numeric_limits<double>::max();
  udomain[1] = std::numeric_limits<double>::lowest();
  int num_pts = 0;
  for (size_t ki=0; ki<sf_ix.size(); ++ki)
    {
      HedgeSurface* surf = surfaces_[sf_ix[ki]].get();
      shared_ptr<ElementarySurface> elem =
	dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf->surface());
      if (!elem.get())
	continue;
      Point loc2 = elem->location();
      Point loc2_2 = loc + ((loc2 - loc)*axis)*axis;

      vector<RevEngPoint*> points;
      vector<RevEngRegion*> reg = surf->getRegions();
      for (size_t kj=0; kj<reg.size(); ++kj)
	{
	  points.insert(points.end(), reg[kj]->pointsBegin(),
			reg[kj]->pointsEnd());
	}

      shared_ptr<ElementarySurface> elem2;
      if (elem->instanceType() == Class_Cylinder)
	elem2 = RevEngUtils::cylinderWithAxis(points, axis, Cx, loc2_2);
      else if (elem->instanceType() == Class_Cone)
	{
	  double len = bbox_.low().dist(bbox_.high());
	  elem2 = RevEngUtils::coneWithAxis(points, axis, Cx, loc2_2, len);
	  if (!elem2.get())
	    elem2 = elem;
	}
      else if (elem->instanceType() == Class_Torus)
	elem2 = RevEngUtils::torusWithAxis(points, axis, Cx, loc2_2);
      else if (elem->instanceType() == Class_Sphere)
	elem2 = RevEngUtils::sphereWithAxis(points, axis, Cx, loc2_2);

      double maxdist = 0.0, avdist = 0.0;
      int num_in = 0, num2_in = 0, num_pt = 0;
      double maxdist_init = 0.0, avdist_init = 0.0;
      int num_in_init = 0, num2_in_init = 0;
      vector<vector<double> > parvals(reg.size());
      vector<vector<pair<double,double> > > dist_ang(reg.size());
      vector<int> init_flag(reg.size());;
      for (size_t kj=0; kj<reg.size(); ++kj)
	{
	  double maxd0, avd0;
	  int num_in0, num2_in0;
	  vector<RevEngPoint*> in0, out0;
	  vector<pair<double,double> > distang0;
	  vector<double> parvals0;
	  RevEngUtils::distToSurf(reg[kj]->pointsBegin(), reg[kj]->pointsEnd(),
				  elem2, approx_tol_, maxd0, avd0, num_in0, num2_in0,
				  in0, out0, parvals0, distang0, angtol);
	  int num0 = reg[kj]->numPoints();
	  maxdist = std::max(maxdist, maxd0);
	  avdist = (double)num_pt*avdist + (double)num0*avd0;
	  num_in += num_in0;
	  num2_in += num2_in0;
	  parvals[kj] = parvals0;
	  dist_ang[kj] = distang0;

	  double maxd1, avd1;
	  int num_in1, num2_in1;
	  reg[kj]->getAccuracy(maxd1, avd1, num_in1, num2_in1);
	  maxdist_init = std::max(maxdist_init, maxd1);
	  avdist_init = (double)num_pt*avdist_init + (double)num0*avd1;
	  num_in_init += num_in1;
	  num2_in_init += num2_in1;
	  num_pt += num0;
	  init_flag[kj] = reg[kj]->getSurfaceFlag();
	}
      avdist /= (double)num_pt;
      avdist_init /= (double)num_pt;

      bool cyl_like = (elem->instanceType() == Class_Cylinder ||
		       elem->instanceType() == Class_Cone);
      int surf_flag = reg[0]->defineSfFlag(num_pt, 0, approx_tol_, num_in,
					   num2_in, avdist, cyl_like);
      int surf_flag_init = reg[0]->defineSfFlag(num_pt, 0, approx_tol_, num_in_init,
						num2_in_init, avdist_init, cyl_like);

      double av_fac = 1.1;
      double in_fac = 0.9;
      if (surf_flag == 0 || surf_flag <= surf_flag_init ||
	  (avdist < av_fac*avdist_init &&
	   (double)num2_in > in_fac*(double)num2_in_init &&
	   (double)num_in > in_fac*(double)num_in_init))
	{
	  // Replace surface
	  for (size_t kj=0; kj<reg.size(); ++kj)
	    {
	      reg[kj]->updateSurfaceAndInfo(elem2, approx_tol_, angtol,
					    parvals[kj], dist_ang[kj],
					    nopar_edgs);
	    }
	}
      for (size_t kj=0; kj<reg.size(); ++kj)
	{
	  num_pts += reg[kj]->numPoints();
	  double dom[4];
	  reg[kj]->getDomain(dom);
	  udomain[0] = std::min(udomain[0], dom[0]);
	  udomain[1] = std::max(udomain[1], dom[1]);
	}      
      int stop_break = 1;
    }

  double seclim = 0.125*M_PI;
  if (num_pts < min_point_region_ && udomain[1]-udomain[0] < seclim)
    return dummy_loc;
  return loc;
}


//===========================================================================
void RevEng::collectAxis(vector<SurfaceProperties>& sfprop)
//===========================================================================
{
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      int code;
      ClassType type = surfaces_[ki]->instanceType(code);
      ClassType type2 = Class_Unknown;
      shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
      shared_ptr<ElementarySurface> elemsf =
	dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf);

      RevEngRegion *reg0 = surfaces_[ki]->getRegion(0);
      if (!reg0)
	continue;  // Surface should have been remove prior to this!!!
      
      if (reg0->hasBlendEdge())
	continue;  // Derived information
      
      if (reg0->hasAssociatedBlend())
	continue;  // Outdated information
      
      int nreg = surfaces_[ki]->numRegions();
      int num_pts = 0;
      shared_ptr<ParamSurface> primary;
;
      int num_pt_primary = 0;
      bool prefer_base = false;  // Maybe check accuracy between surface
      // and primary surface later. Can also use history information about
      // surface updates
      bool type_cyl = (type == Class_Cylinder || type == Class_Cone);
      int surfflag = 0;
      for (int ka=0; ka<nreg; ++ka)
	{
	  RevEngRegion *reg = surfaces_[ki]->getRegion(ka);
	  primary = reg->getBase();
	  int num = reg->numPoints();
	  num_pts += num;
	  double maxd, avd;
	  int num_in, num2_in;
	  reg->getAccuracy(maxd, avd, num_in, num2_in);
	  int surfflag0 = reg->defineSfFlag(num, 0, approx_tol_, num2_in,
					   num2_in, avd, type_cyl);
	  surfflag = std::max(surfflag, surfflag0);
	  vector<shared_ptr<ElementarySurface> > sfs(primary.get() ? 2: 1);
	  sfs[0] = elemsf;
	  if (primary.get())
	    sfs[1] = dynamic_pointer_cast<ElementarySurface,ParamSurface>(primary);
	  double maxd2, avd2;
	  int num_in2, num2_in2, surfflag2;
	  int ix = reg->checkSurfaceAccuracy(sfs, approx_tol_, 5.0*anglim_, maxd2,
					     avd2, num_in2, num2_in2, surfflag2);
	  
	  if (reg->hasBaseSf() && num > num_pt_primary)
	    {
	      double maxdp, avdp;
	      int num_inp, num2_inp;
	      reg->getBaseDist(maxdp, avdp, num_inp, num2_inp);
	      if (num_inp > num/2 && avdp < approx_tol_)
		{
		  num_pt_primary = num;
		}
	    }
	}
      if (primary.get())
	type2 = primary->instanceType();
      
      if ((!elemsf.get()) || prefer_base)
	{
	  if (primary.get())
	    elemsf =
	      dynamic_pointer_cast<ElementarySurface,ParamSurface>(primary);
	}
      if (!elemsf.get())
	continue;
      
      Point loc, dir;
      double rad1, rad2;
      loc = elemsf->location();
      dir = elemsf->direction();
      rad1 = elemsf->radius(0.0, 0.0);   // Not good enough for cones
      rad2 = elemsf->radius2(0.0, 0.0);   // Not good enough for cones
      SurfaceProperties currprop((int)ki, type, num_pts, surfflag, dir, loc, 
				 type2, rad1, rad2);
      sfprop.push_back(currprop);
    }
}


//===========================================================================
void RevEng::trimSurfaces()
//===========================================================================
{
  double angtol = 5.0*anglim_;
#ifdef DEBUG_TRIM
  std::ofstream of0("trimmedsfs.g2");
#endif  
  // Missing edges between surfaces
  size_t num_edg = edges_.size();
  recognizeEdges(true);

  vector<RevEngRegion*> out_regs;
  vector<HedgeSurface*> out_sfs;
  for (size_t kj=0; kj<edges_.size(); ++kj)
    {
      int type = edges_[kj]->getType();
      if (type == NOT_BLEND)
	edges_[kj]->setTrimCurves(approx_tol_, angtol, out_regs, out_sfs);
    }
  if (out_regs.size() > 0 || out_sfs.size() > 0)
    {
      int ix = -1;
      updateRegionsAndSurfaces(ix, out_regs, out_sfs);
    }
#ifdef DEBUG_TRIM
      vector<shared_ptr<ftEdge> > trim_edgs1;
      for (size_t kr=0; kr<regions_.size(); ++kr)
	{
	  // if (regions_[kr]->numPoints() == 0)
	  //   std::cout << "Finished set blend boundaries, empty region, ki=" << kr << ", region: " << regions_[kr].get() << std::endl;
	  int num = regions_[kr]->numTrimEdges();
	  if (num > 0)
	    {
	      vector<shared_ptr<ftEdge> > curr_edgs = regions_[kr]->getTrimEdges();
	      trim_edgs1.insert(trim_edgs1.end(), curr_edgs.begin(), curr_edgs.end());
	    }
	}
      if (trim_edgs1.size() > 0)
	{
	  std::ofstream oft1("trim_edges_1.g2");
	  for (size_t kr=0; kr<trim_edgs1.size(); ++kr)
	    {
	      shared_ptr<ParamCurve> tmp = trim_edgs1[kr]->geomCurve();
	      shared_ptr<CurveOnSurface> tmp2 =
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp);
	      shared_ptr<ParamCurve> tmp3 = tmp2->spaceCurve();
	      tmp3->writeStandardHeader(oft1);
	      tmp3->write(oft1);
	    }
	}
#endif
      //#if 0
   for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      if (!regions_[ki]->hasSurface())
	continue;
#ifdef DEBUG_TRIM
      std::ofstream of("trimreg.g2");
      regions_[ki]->writeRegionPoints(of);
#endif

      if (regions_[ki]->numTrimEdges() == 0)
	continue;
      if (regions_[ki]->hasBlendEdge())
	continue;
      regions_[ki]->extendBoundaries(mean_edge_len_, min_point_region_,
				     approx_tol_, angtol, mainaxis_);
    }
   //#endif      
#ifdef DEBUG_TRIM
      vector<shared_ptr<ftEdge> > trim_edgs2;
      for (size_t kr=0; kr<regions_.size(); ++kr)
	{
	  // if (regions_[kr]->numPoints() == 0)
	  //   std::cout << "Finished set blend boundaries, empty region, ki=" << kr << ", region: " << regions_[kr].get() << std::endl;
	  int num = regions_[kr]->numTrimEdges();
	  if (num > 0)
	    {
	      vector<shared_ptr<ftEdge> > curr_edgs = regions_[kr]->getTrimEdges();
	      trim_edgs2.insert(trim_edgs2.end(), curr_edgs.begin(), curr_edgs.end());
	    }
	}
      if (trim_edgs2.size() > 0)
	{
	  std::ofstream oft2("trim_edges_2.g2");
	  for (size_t kr=0; kr<trim_edgs2.size(); ++kr)
	    {
	      shared_ptr<ParamCurve> tmp = trim_edgs2[kr]->geomCurve();
	      shared_ptr<CurveOnSurface> tmp2 =
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp);
	      shared_ptr<ParamCurve> tmp3 = tmp2->spaceCurve();
	      tmp3->writeStandardHeader(oft2);
	      tmp3->write(oft2);
	    }
	}
#endif

  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      bool trimmed;
      if (!regions_[ki]->hasSurface())
	continue;
#ifdef DEBUG_TRIM
      std::ofstream of2("trimreg2.g2");
      regions_[ki]->writeRegionPoints(of2);
#endif
      if (regions_[ki]->numTrimEdges() == 0)
	{
	  trimmed = true;
	  regions_[ki]->getSurface(0)->trimWithPoints(approx_tol_);
	}
      else
	{
	  try {
	    trimmed = regions_[ki]->trimSurface(approx_tol_);
	  }
	  catch (...)
	    {
	      trimmed = regions_[ki]->getSurface(0)->trimWithPoints(approx_tol_);
#ifdef DEBUG_TRIM
	      std::cout << ki << " trim with points " << std::endl;
#endif
	    }
	  if (!trimmed)
	    {
	      trimmed = regions_[ki]->getSurface(0)->trimWithPoints(approx_tol_);
#ifdef DEBUG_TRIM
	      std::cout << ki << " trim with points " << std::endl;
#endif
	    }
	}
      
#ifdef DEBUG_TRIM
      if (trimmed)
	{
	  shared_ptr<ParamSurface> tsurf = regions_[ki]->getSurface(0)->surface();
	  shared_ptr<BoundedSurface> bdsurf = dynamic_pointer_cast<BoundedSurface,ParamSurface>(tsurf);
	  if (bdsurf.get())
	    {
	      int valid_state;
	      bool valid = bdsurf->isValid(valid_state);
	      std::cout << "BoundedSurf " << ki << " is valid? " << valid << " " << valid_state << std::endl;
	    }
	  tsurf->writeStandardHeader(of2);
	  tsurf->write(of2);
	}
#endif
	
    }
#ifdef DEBUG
   std::cout << "Finished trim surf, regions: " << regions_.size() << ", surfaces: " << surfaces_.size() << std::endl;
   if (regions_.size() > 0)
    {
      std::cout << "Regions13" << std::endl;
      std::ofstream of("regions13.g2");
      std::ofstream ofm("mid_regions13.g2");
      std::ofstream ofs("small_regions13.g2");
      writeRegionStage(of, ofm, ofs);
      std::ofstream of0("regions13_helix.g2");
      for (int ki=0; ki<(int)regions_.size(); ++ki)
	{
	  if (regions_[ki]->getSurfaceFlag() == PROBABLE_HELIX)
	    {
	      regions_[ki]->writeRegionInfo(of0);
	      if (regions_[ki]->hasSurface())
		regions_[ki]->writeSurface(of0);
	    }
	}
      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf13.g2");
	  writeRegionWithSurf(of);
	}
     }
   
   if (edges_.size() > 0)
     {
       std::ofstream ofe("edges13.g2");
       writeEdgeStage(ofe);
     }
#endif
  
#if 0
#ifdef DEBUG
  std::ofstream of1("surfbd.g2");
#endif
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      // Restrict unbounded surfaces
      surfaces_[ki]->ensureSurfaceBounded();
#ifdef DEBUG
      surfaces_[ki]->surface()->writeStandardHeader(of1);
      surfaces_[ki]->surface()->write(of1);
#endif
    }

#ifdef DEBUG
  std::ofstream of3("trimsurfs.g2");
#endif
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      surfaces_[ki]->trimWithPoints(approx_tol_);

#ifdef DEBUG
      surfaces_[ki]->surface()->writeStandardHeader(of3);
      surfaces_[ki]->surface()->write(of3);
#endif
      int stop1 = 1;
    }
#endif
  int stop_break = 1;
}

 //===========================================================================
shared_ptr<SurfaceModel> RevEng::createModel()
//===========================================================================
{
  // Ensure bounded surfaces
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
      if (!surf->isBounded())
	{
	  RevEngRegion *reg = surfaces_[ki]->getRegion(0);
	  double dom[4];
	  reg->getDomain(dom);
	  shared_ptr<Plane> plane = dynamic_pointer_cast<Plane,ParamSurface>(surf);
	  shared_ptr<Cylinder> cyl = dynamic_pointer_cast<Cylinder,ParamSurface>(surf);
	  shared_ptr<Cone> cone = dynamic_pointer_cast<Cone,ParamSurface>(surf);
	  if (plane.get())
	    plane->setParameterBounds(dom[0], dom[2], dom[1], dom[3]);
	  else if (cyl.get())
	    cyl->setParamBoundsV(dom[2], dom[3]);
	  else if (cone.get())
	    cone->setParamBoundsV(dom[2], dom[3]);
 	}
    }

  double gap_tol = std::max(10.0*int_tol_, 0.01*approx_tol_);
  vector<shared_ptr<ftSurface> > tmpsfs(surfaces_.begin(), surfaces_.end());
  sfmodel_ = shared_ptr<SurfaceModel>(new SurfaceModel(approx_tol_, gap_tol,
						       10*gap_tol, anglim_, 10*anglim_,
						       tmpsfs));
#ifdef DEBUG_MODEL
  int num_bd = sfmodel_->nmbBoundaries();
  if (num_bd > 0)
    {
      std::ofstream of1("model_boundaries.g2");
      for (int ka=0; ka<num_bd; ++ka)
	{
	  ftCurve bd = sfmodel_->getBoundary(ka);
	  int num_seg = bd.numSegments();
	  for (int kb=0; kb<num_seg; ++kb)
	    {
	      shared_ptr<ParamCurve> cv = bd.segment(kb).spaceCurve();
	      cv->writeStandardHeader(of1);
	      cv->write(of1);
	    }
	}
    }
  
#endif
  return sfmodel_;
}

 //===========================================================================
void RevEng::initParameters()
//===========================================================================
{
  // Set default parameters
  min_next_ = 10;  // Minimum number of neighbouring points
  max_next_ = 500; //std::min(80, tri_sf_->size()/200); //500;
  rfac_ = 6.0; //3.0;  // Factor for radius in which to search for neighbouring points
  cfac_ = 8.0;  // Edge points from curvature is given by
  // cfac_ times the average length of triangulation edges in a vertex
  norm_plane_lim_= 0.005; // Limit for when the cone angle corresponding
  // to triangle normals indicate an edge
  zero_H_ = 0.005; //0.001; //0.0001;  // When mean curvature is considered zero
  zero_K_ = 0.005; //0.001; //0.0001;  // When Gauss curvature is considered zero
  norm_ang_lim_ = 0.1*M_PI; // Limit for when the cone angle corresponding
    // to triangle normals indicate an edge
  min_point_region_ = 200; //50; //10;  // Should be updated with regard to the total
  // number of points
  approx_tol_ = 0.001;  // Very preliminary
  int_tol_ = 1.0e-6;
  anglim_ = 0.01;
  max_nmb_outlier_ = 3;
  
  model_character_  = 2;
  prefer_elementary_ = 0;
  mainaxis_[0] = Point(1.0, 0.0, 0.0);
  mainaxis_[1] = Point(0.0, 1.0, 0.0);
  mainaxis_[2] = Point(0.0, 0.0, 1.0);
}

 //===========================================================================
void RevEng::updateParameters()
//===========================================================================
{
  if (model_character_ == SMOOTH)
    {
      rfac_ = 4.0;
    }
  else if (model_character_ == MEDIUM_ROUGH)
    {
      zero_H_ = 0.007;
      zero_K_ = 0.007;
    }
  else
    {
      rfac_ = 6.0;
      zero_H_ = 0.01;
      zero_K_ = 0.01;
      anglim_ = 0.02;
    }
}

 //===========================================================================
int RevEng::setSmallRegionNumber()
//===========================================================================
{
  vector<int> nmb_pt_reg(regions_.size());
  for (size_t ki=0; ki<regions_.size(); ++ki)
    nmb_pt_reg[ki] = regions_[ki]->numPoints();

  std::sort(nmb_pt_reg.begin(), nmb_pt_reg.end());
  int tot_num = tri_sf_->size();
  int num_reg = (int)regions_.size();
  int idel = tot_num/num_reg;
  int min_num = std::min(10, tot_num/20);
  int ixmax = (int)(0.99*num_reg);
  int Q4 = 3*num_reg/4;
  int max_num = std::max(min_num, nmb_pt_reg[ixmax]);
  int Q4_num = nmb_pt_reg[Q4];
  max_num = std::min(max_num, 10*idel);
  int ixdel = std::max(num_reg/100, 2);
  int prev = nmb_pt_reg[0], prev0 = 0;
  int ix;
  int fac = 2;
  int min_jump = std::min(idel, Q4_num); //2;
  for (ix=ixdel; ix<num_reg; ix+=ixdel)
    {
      int diff = nmb_pt_reg[ix] - prev;
      if (diff > fac*(std::max(min_jump, prev-prev0)))
	break;
      if (diff > 0)
	prev0 = prev;
      prev = nmb_pt_reg[ix];
    }
  ix = std::min(ix, ixmax);
      
  int num = std::max(min_num, std::min(nmb_pt_reg[ix], max_num/2));
  return num;
}


 //===========================================================================
void RevEng::checkConsistence(std::string text) const
//===========================================================================
{
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      vector<RevEngRegion*> adjacent;
      regions_[ki]->getAdjacentRegions(adjacent);
      for (size_t kj=0; kj<adjacent.size(); ++kj)
	{
	  size_t kr;
	  for (kr=0; kr<regions_.size(); ++kr)
	    if (adjacent[kj] == regions_[kr].get())
	      break;
	  if (kr == regions_.size())
	    std::cout << text << ", Obsolete region pointer, ki=" << ki << ", kj=" << kj << ". Region: " << adjacent[kj] << std::endl;
	}
    }
  for (int ki=0; ki<(int)surfaces_.size(); ++ki)
    {
      int numreg = surfaces_[ki]->numRegions();
      for (int ka=0; ka<numreg; ++ka)
	{
	  RevEngRegion *reg = surfaces_[ki]->getRegion(ka);
	  size_t kr;
	  for (kr=0; kr<regions_.size(); ++kr)
	    if (reg == regions_[kr].get())
	      break;
	  if (kr == regions_.size())
	    std::cout << text << ", surface 1. Obsolete region pointer, ki=" << ki << ", ka=" << ka << ". Region: " << reg << std::endl;
	  vector<RevEngRegion*> adjacent;
	  reg->getAdjacentRegions(adjacent);
	  for (size_t kj=0; kj<adjacent.size(); ++kj)
	    {
	      size_t kr;
	      for (kr=0; kr<regions_.size(); ++kr)
		if (adjacent[kj] == regions_[kr].get())
		  break;
	      if (kr == regions_.size())
		std::cout << text << ", surface. Obsolete region pointer, ki=" << ki << ", kj=" << kj << ". Region: " << adjacent[kj] << std::endl;
	    }
	}
    }
}

 //===========================================================================
void RevEng::storeClassified(ostream& os) const
//===========================================================================
{
  storeParams(os);
  int nmbpts = tri_sf_->size();
  os << nmbpts << std::endl;
  for (int ki=0; ki<nmbpts; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      pt->store(os);
    }
}

 //===========================================================================
void RevEng::readClassified(istream& is)
//===========================================================================
{
  readParams(is);
  int nmbpts;
  is >> nmbpts;
  tri_sf_ = shared_ptr<ftPointSet>(new ftPointSet());
  vector<vector<int> > next_ix(nmbpts);
  for (int ki=0; ki<nmbpts; ++ki)
    {
      shared_ptr<RevEngPoint> vertex(new RevEngPoint());
      vertex->read(is, next_ix[ki]);
      tri_sf_->addEntry(vertex);
    }

  // Add next information
  for (int ki=0; ki<nmbpts; ++ki)
    {
      ftSamplePoint* pt1 = (*tri_sf_)[ki];
      for (size_t kr=0; kr<next_ix[ki].size(); ++kr)
	{
	  int ix = next_ix[ki][kr];
	  ftSamplePoint* pt2 = (*tri_sf_)[ix];
	  pt1->addNeighbour(pt2);
	}
    }

  max_next_ = std::min(80, tri_sf_->size()/200);
  max_next_ = std::max(2*min_next_, max_next_);
  setBoundingBox();
}

 //===========================================================================
void RevEng::storeGrownRegions(ostream& os)
//===========================================================================
{
  storeClassified(os);
  os << surfaces_.size() << std::endl;
  for (int ka=0; ka<(int)surfaces_.size(); ++ka)
    {
      surfaces_[ka]->setId(ka);
      surfaces_[ka]->store(os);
    }
  os << regions_.size() << std::endl;
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setId((int)ki);
      regions_[ki]->store(os);
    }

  os << single_points_.size() << std::endl;
  for (size_t ki=0; ki<single_points_.size(); ++ki)
    os << single_points_[ki]->getIndex() << " ";
  os << std::endl;
  
  os << edges_.size() << std::endl;
  for (int ka=0; ka<(int)edges_.size(); ++ka)
    {
      edges_[ka]->setId(ka);
      edges_[ka]->store(os);
    }

  os << model_axis_.size() << std::endl;
  for (size_t ki=0; ki<model_axis_.size(); ++ki)
    {
      os << model_axis_[ki].axis_ << " " << model_axis_[ki].plane_loc_.size();
      os <<  " " << model_axis_[ki].rotational_loc_.size() << std::endl;
      for (size_t kj=0; kj<model_axis_[ki].plane_loc_.size(); ++kj)
	{
	  os << model_axis_[ki].plane_loc_[kj].first << " ";
	  os << model_axis_[ki].plane_loc_[kj].second << " ";
	}
      os << std::endl;
      for (size_t kj=0; kj<model_axis_[ki].rotational_loc_.size(); ++kj)
	{
	  os << model_axis_[ki].rotational_loc_[kj].first << " ";
	  os << model_axis_[ki].rotational_loc_[kj].second << " ";
	}
      os << std::endl;
    }
      
}

 //===========================================================================
void RevEng::readGrownRegions(istream& is)
//===========================================================================
{
  readClassified(is);
  int nmb = tri_sf_->size();
  vector<ftSamplePoint*> tmp_pts(nmb);
  for (int ka=0; ka<nmb; ++ka)
    tmp_pts[ka] = (*tri_sf_)[ka];
  std::set<ftSamplePoint*> tmp_set(tmp_pts.begin(), tmp_pts.end());
  std::cout << "Read grown, size1 = " << tmp_pts.size() << ", size2 = " << tmp_set.size() << std::endl;
  curvatureFilter();
  int num_sfs;
  is >> num_sfs;
  if (num_sfs > 0)
    surfaces_.resize(num_sfs);
  for (int ki=0; ki<num_sfs; ++ki)
    {
      surfaces_[ki] = shared_ptr<HedgeSurface>(new HedgeSurface());
      surfaces_[ki]->read(is);
    }
  
  int num_regions;
  is >> num_regions;
  regions_.resize(num_regions);
  for (int ki=0; ki<num_regions; ++ki)
    {
      vector<int> sf_id;
      regions_[ki] = shared_ptr<RevEngRegion>(new RevEngRegion(edge_class_type_));
      regions_[ki]->read(is, tri_sf_, sf_id);
      for (size_t kj=0; kj<sf_id.size(); ++kj)
	{
	  for (size_t kr=0; kr<surfaces_.size(); ++kr)
	    {
	      if (sf_id[kj] == surfaces_[kr]->getId())
		{
		  regions_[ki]->addHedge(surfaces_[kr].get());
		  surfaces_[kr]->addRegion(regions_[ki].get());
		  break;
		}
	    }
	}
    }

  for (int ki=0; ki<num_regions; ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

  int num_single;
  is >> num_single;
  single_points_.resize(num_single);
  for (int ki=0; ki<num_single; ++ki)
    {
      int ix;
      is >> ix;
      RevEngPoint* pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ix]);
      single_points_[ki] = pt;
    }
  
  int num_edgs;
  is >> num_edgs;
  if (num_edgs > 0)
    edges_.resize(num_edgs);
  for (int ki=0; ki<num_edgs; ++ki)
    {
      edges_[ki] = shared_ptr<RevEngEdge>(new RevEngEdge());
      int id1, id2, id3;
      vector<int> blend_id;
      edges_[ki]->read(is, id1, id2, id3, blend_id);

      for (size_t kr=0; kr<regions_.size(); ++kr)
	{
	  int id = regions_[kr]->getId();
	  if (id == id1)
	    {
	      edges_[ki]->setReg1(regions_[kr].get());
	      regions_[kr]->addRevEdge(edges_[ki].get());
	      break;
	    }
	}
      
      for (size_t kr=0; kr<regions_.size(); ++kr)
	{
	  int id = regions_[kr]->getId();
	  if (id == id2)
	    {
	      edges_[ki]->setReg2(regions_[kr].get());
	      regions_[kr]->addRevEdge(edges_[ki].get());
	      break;
	    }
	}

      if (id3 >= 0)
	{
	  for (size_t kr=0; kr<regions_.size(); ++kr)
	    {
	      int id = regions_[kr]->getId();
	      if (id == id3)
		{
		  edges_[ki]->setBlendRegSurf(regions_[kr].get());
		  regions_[kr]->setBlendEdge(edges_[ki].get());
		  break;
		}
	    }
	}

      for (size_t kh=0; kh<blend_id.size(); ++kh)
	{
	  for (size_t kr=0; kr<regions_.size(); ++kr)
	    {
	      int id = regions_[kr]->getId();
	      if (blend_id[kh] == id)
		{
		  edges_[ki]->addBlendRegion(regions_[kr].get());
		  regions_[kr]->setAssociatedBlend(edges_[ki].get());
		  break;
		}
	    }
	}
    }

  int num_axis;
  is >> num_axis;
  for (int ki=0; ki<num_axis; ++ki)
    {
      Point axis(3);
      int num_planar, num_rotational;
      is >> axis[0] >> axis[1] >> axis[2]  >> num_planar >> num_rotational;
      model_axis_.push_back(AxisInfo(axis));
      for (int kj=0; kj<num_planar; ++kj)
	{
	  Point loc(3);
	  int num;
	  is >> loc[0] >> loc[1] >> loc[2] >> num;
	  model_axis_[ki].addPlaneLocation(loc, num);
	}
      for (int kj=0; kj<num_rotational; ++kj)
	{
	  Point loc(3);
	  int num;
	  is >> loc[0] >> loc[1] >> loc[2] >> num;
	  model_axis_[ki].addRotationalLocation(loc, num);
	}
    }
  int stop_break = 1;
}
  

 //===========================================================================
void RevEng::storeParams(ostream& os) const
//===========================================================================
{
  os <<  model_character_ << " " << mean_edge_len_ << " " << min_next_;
  os << " " << norm_ang_lim_;
  os << " " << norm_plane_lim_ << " " << zero_H_ << " " << zero_K_;
  os << " " << min_point_region_ << " " << approx_tol_ ;
  os << " " << anglim_ << " " << max_nmb_outlier_ << " ";
  os << edge_class_type_ << " " << classification_type_ << std::endl;
  os << mainaxis_[0] << " " << mainaxis_[1] << " " << mainaxis_[2] << std::endl;
}

 //===========================================================================
void RevEng::readParams(istream& is)
//===========================================================================
{
  is >>  model_character_ >> mean_edge_len_ >> min_next_;
  is >> norm_ang_lim_ >> norm_plane_lim_ >> zero_H_ >> zero_K_;
  is >> min_point_region_ >> approx_tol_ >> anglim_ >> max_nmb_outlier_;
  is >> edge_class_type_ >> classification_type_;
  mainaxis_[0].resize(3);
  mainaxis_[1].resize(3);
  mainaxis_[2].resize(3);
  is >> mainaxis_[0] >> mainaxis_[1] >> mainaxis_[2];
}

//===========================================================================
void RevEng::writeRegionWithSurf(ostream& of) const
//===========================================================================
{
  for (size_t kr=0; kr<surfaces_.size(); ++kr)
    {
      int numreg = surfaces_[kr]->numRegions();
      if (numreg > 0)
	{
	  RevEngRegion *reg = surfaces_[kr]->getRegion(0);
	  reg->writeRegionPoints(of);
	  reg->writeSurface(of);
	}
      else
	{
	  double diag = bbox_.low().dist(bbox_.high());
	  surfaces_[kr]->limitSurf(diag);
	  surfaces_[kr]->surface()->writeStandardHeader(of);
	  surfaces_[kr]->surface()->write(of);
	}
    }
}
  
 //===========================================================================
void RevEng::writeRegionStage(ostream& of, ostream& ofm, ostream& ofs) const
//===========================================================================
{
  std::cout << "Num regions: " << regions_.size() << ", num surfaces: " << surfaces_.size() << " num edges: " << edges_.size() << std::endl;
  
  vector<Vector3D> small;
  int nmb_one = 0;
  int low = min_point_region_/10;
  for (size_t kr=0; kr<regions_.size(); ++kr)
    {
      // BoundingBox bbox = regions_[kr]->boundingBox();
      // if (bbox.low().dist(bbox.high()) < 0.1)
      //   std::cout << "Small bounding box" << std::endl;
      // std::set<RevEngPoint*> tmpset(regions_[kr]->pointsBegin(), regions_[kr]->pointsEnd());
      // if (tmpset.size() != regions_[kr]->numPoints())
      // 	std::cout << "Point number mismatch. " << kr << " " << tmpset.size() << " " << regions_[kr]->numPoints() << std::endl;
      int nmb = regions_[kr]->numPoints();
      if (nmb < low && regions_[kr]->hasSurface() == false)
	{
	  for (int ki=0; ki<nmb; ++ki)
	    small.push_back(regions_[kr]->getPoint(ki)->getPoint());
	}
      else if (nmb < min_point_region_ && regions_[kr]->hasSurface() == false)
	{
	  ofm << "400 1 0 0" << std::endl;
	  ofm << nmb << std::endl;
	  for (int ki=0; ki<nmb; ++ki)
	    {
	      ofm << regions_[kr]->getPoint(ki)->getPoint() << std::endl;
	    }
	  if (regions_[kr]->hasSurface())
	    regions_[kr]->writeSurface(ofm);
	}
      else
	{
	  if (nmb > 0)
	    {
	      of << "400 1 0 0" << std::endl;
	      of << nmb << std::endl;
	      for (int ki=0; ki<nmb; ++ki)
		{
		  of << regions_[kr]->getPoint(ki)->getPoint() << std::endl;
		}
	    }
	  if (regions_[kr]->hasSurface())
	    regions_[kr]->writeSurface(of);
	}
      if (nmb == 1)
	nmb_one++;
    }
  std::cout << "Number of regions with one point: " << nmb_one << std::endl;
  ofs << "400 1 0 4 0 0 0 255" << std::endl;
  ofs << small.size() << std::endl;
  for (size_t kr=0; kr<small.size(); ++kr)
    ofs << small[kr] << std::endl;
}

//===========================================================================
void RevEng::writeEdgeStage(ostream& of) const
//===========================================================================
{
  for (size_t kr=0; kr<edges_.size(); ++kr)
    {
      vector<shared_ptr<ParamCurve> > cvs = edges_[kr]->getSpaceCurves();
      for (size_t kh=0; kh<cvs.size(); ++kh)
	{
	  cvs[kh]->writeStandardHeader(of);
	  cvs[kh]->write(of);
	}
      int num_blend = edges_[kr]->numBlendRegs();
      for (int ka=0; ka<num_blend; ++ka)
	edges_[kr]->getBlendReg(ka)->writeRegionPoints(of);
    }
}
