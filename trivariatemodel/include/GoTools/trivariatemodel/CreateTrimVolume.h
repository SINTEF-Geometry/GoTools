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

#ifndef _CREATETRIMVOLUME_H
#define _CREATETRIMVOLUME_H

#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/utils/DirectionCone.h"

namespace Go
{
  class ftSurface;
  class ftVolume;
  class ParamVolume;
  class ParamSurface;

  enum sf_type
  {
    UNKNOWN = -1,
    FREEFORM = 0,
    PLANAR = 1,
    ROTATIONAL = 2
  };

  /// Create a trimmed volume from a boundary represented solid respecting 
  /// characteristica of the input shape
  /// Currently, only special configurations will give a result different from the
  /// trivial one (a spline volume defined from the axis parallel bounding box of
  /// the model is trimmed by the model itself). Rotational models where the rotational
  /// surfaces are cylinders and cones should give a more boundary fitted result.

  class CreateTrimVolume
  {
  public:
    /// Constructor
    /// model : The trimming shell / outer shell of the brep solid
    /// material: Eventual material specification associated with the brep solid
    ///           that should be maintained in the trivariate model
    CreateTrimVolume(shared_ptr<SurfaceModel> model, int material=-1);

    /// Destructor
    ~CreateTrimVolume();

    /// Add information about voids
    void addVoids(std::vector<shared_ptr<SurfaceModel> >& voids)
    {
      voids_.insert(voids_.end(), voids.begin(), voids.end());
    }

    /// Define a rotational trimmed volume if possible and fetch result
    shared_ptr<ftVolume> fetchRotationalTrimVol(bool create_degen = true,
						bool refine_sharp = false);
 
   /// Define the trimmed volume and fetch result
    shared_ptr<ftVolume> fetchOneTrimVol(bool refine_sharp = false);

  private:
    shared_ptr<SurfaceModel> model_;
    std::vector<shared_ptr<SurfaceModel> > voids_;
    int material_;
    BoundingBox bigbox_;

    vector<vector<shared_ptr<ftSurface> > > face_grp_;
    vector<shared_ptr<ParamSurface> > under_sf_;
    vector<BoundingBox> bbox_;
    vector<DirectionCone> cone_;
    vector<double> sfsize_;
    vector<sf_type> sf_type_; // Unknown, freeform, planar, rotational
    vector<Point> sf_pt_;    
    vector<Point> sf_axis_;  // Only set for rotational surfaces
    vector<Point> sf_centre_;  // Only set for rotational surfaces

    bool identifyBoundaryFaces(std::vector<std::pair<shared_ptr<ftSurface>, shared_ptr<ParamSurface> > >& side_sfs);

    void extractMaxSet(std::vector<shared_ptr<ftSurface> >& bd_faces,
		       std::vector<shared_ptr<ftSurface> >& trim_faces);

    shared_ptr<ftVolume> 
      createTrimVolume(shared_ptr<ParamVolume> vol, 
		       std::vector<std::pair<shared_ptr<ftSurface>, shared_ptr<ParamSurface> > >& side_sfs);
    
    void identifyInnerTrim(std::vector<shared_ptr<ftSurface> >& bd_faces,
			   std::vector<shared_ptr<ftSurface> >& trim_faces);

    void computeGroupInfo(double tol);

    void findSideSfs(double tol, double angtol,
		     std::vector<std::pair<shared_ptr<ftSurface>, shared_ptr<ParamSurface> > >& side_sfs);

    int
      checkCandPair(Point vec, shared_ptr<ParamSurface> sf1, int bd_type1, 
		    BoundingBox& box1, shared_ptr<ParamSurface> sf2,
		    int bd_type2, BoundingBox& box2, 
		    double tol);

    void
      oneSideSf(int bd_type, std::vector<int>& face_group_ix, 
		Point bd_vec, Point dir,
		std::vector<shared_ptr<ftSurface> >& ref_faces,
		double tol, double angtol,
		std::pair<shared_ptr<ftSurface>, shared_ptr<ParamSurface> >& side_sf);
    void extendSurfaces(std::vector<std::pair<shared_ptr<ftSurface>, shared_ptr<ParamSurface> > >& side_sfs);

    void trimSideSurfaces(std::vector<std::pair<shared_ptr<ftSurface>, shared_ptr<ParamSurface> > >& side_sfs,
			  std::vector<bool>& test_inner);
    void refineInSharpEdges(shared_ptr<ParamVolume>& vol);

    bool checkIsoPar(shared_ptr<ParamSurface> surf,
		     shared_ptr<ParamVolume> vol,
		     int pardir, double parval, double tol);

    void repairShell(int degree);

    void analyzePrio(int* prio, int nmb_pri, 
		     Point coord[], Point coord_pos[]);

    bool 
      identifyRotationalAxis(Point& centre, Point& axis, 
			     Point& vec, double& angle,
			     std::vector<shared_ptr<ftSurface> >& rotational_faces, 
			     std::vector<shared_ptr<ftSurface> >& other_faces);

    void
      defineRotationalSurfaces(Point centre, Point axis,
			       Point vec, double angle,
			       std::vector<shared_ptr<ftSurface> >& rot_faces,
			       double& rad1, double& rad2,
			       std::vector<std::pair<shared_ptr<ftSurface>, 
			       shared_ptr<ParamSurface> > >& side_sfs);
    void
      defineEndSurfaces(Point centre, Point axis, double rad,
			std::vector<std::pair<shared_ptr<ftSurface>, 
			shared_ptr<ParamSurface> > >& side_sfs);

    void 
      defineRotationalEndSurfaces(Point centre, Point axis, 
				  Point vec, double angle,
				  double rad1, double rad2,
				  std::vector<std::pair<shared_ptr<ftSurface>, 
				  shared_ptr<ParamSurface> > >& side_sfs);
  void
  orientSurfaces(Point centre, Point axis, Point vec, double angle,
		 std::vector<shared_ptr<SplineSurface> >& sfs);

  void limitUnderlyingSurfaces();

  bool 
    updateSideSfs(shared_ptr<SurfaceModel>& shell,
		  std::vector<std::pair<shared_ptr<ftSurface>, shared_ptr<ParamSurface> > >& side_sfs);

  bool checkRotSfExtent(const Point& centre, const Point& axis,
			const Point& vec, double angle, 
			double radius1, double radius2);

  shared_ptr<ParamVolume> 
  defineDegenRot(Point& pos, Point& axis,
		 double angle, 
		 shared_ptr<ParamSurface> outer_sf);
  };

} // namespace Go


#endif // _CREATETRIMVOLUME_H
