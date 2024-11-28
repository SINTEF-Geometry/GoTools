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

#ifndef _REVENGEDGE_H
#define _REVENGEDGE_H

#include "GoTools/compositemodel/RevEngRegion.h"
#include "GoTools/geometry/CurveOnSurface.h"

namespace Go
{
  /// Flag indicating type of associated blend surface. Only partially in use:
  /// BLEND_NOT_SET = blend surface is expected, NOT_BLEND = sharp edge between
  /// surfaces
  enum
  {
   BLEND_NOT_SET, NOT_BLEND, STRAIGHT_BLEND, CIRCULAR_BLEND, OTHER_BLEND
  };

  /** RevEngEdge - Represents an edge between two surfaces associated RevEngRegions.
      Information about intersection curve, adjacent regions, possible blend
      surfaces/regions and other regions associated to this blend
   *
   */

  class RevEngEdge
  {
  public:
    /// Constructor
    RevEngEdge();

    /// Constructor giving the regions whom associated surfaces is intersected
    /// to produce this edge
    RevEngEdge(RevEngRegion* reg1, RevEngRegion* reg2);

    /// Constructor giving information to produce an associated blend surface
    RevEngEdge(int type, RevEngRegion* reg1, 
	       std::vector<shared_ptr<CurveOnSurface> > cvs1,
	       bool out1, RevEngRegion* reg2,
	       std::vector<shared_ptr<CurveOnSurface> > cvs2,
	       bool out2, double radius, double width);

    /// Destructor
    ~RevEngEdge();

    /// Set entity identity. Used in storing of stage (see RevEng)
    void setId(int Id)
    {
      Id_ = Id;
    }

    /// Enquire entity identity
    int getId()
    {
      return Id_;
    }

    /// Enquire type of blend associated to edge
    int getType()
    {
      return blend_type_;
    }

    /// Change first adjacent region and ensure that related information is updated
    void setReg1(RevEngRegion *reg);
    
    /// Change second adjacent region and ensure that related information is updated
    void setReg2(RevEngRegion *reg);

    /// Increase the pool of regions associated to a future blend surface
    void addBlendRegion(RevEngRegion* reg)
    {
      blend_regs_.push_back(reg);
    }

    /// Increase the pool of regions associated to a future blend surface
   void addBlendRegions(std::vector<RevEngRegion*>& regs)
    {
      blend_regs_.insert(blend_regs_.end(), regs.begin(), regs.end());
    }

    /// Remove all regions associated to a blend surface
    void clearBlendRegions()
    {
      blend_regs_.clear();
    }

    /// Remove one region associated to a blend surface
    void removeBlendReg(RevEngRegion* reg)
    {
      auto it = std::find(blend_regs_.begin(), blend_regs_.end(), reg);
      if (it != blend_regs_.end())
	blend_regs_.erase(it);
    }

    /// Fetch geometry space curves of intersection curves between adjacent region surfaces
    /// Expects one curve
    std::vector<shared_ptr<ParamCurve> > getSpaceCurves();

    /// Number of regions associated to a future blend surface
    int numBlendRegs()
    {
      return (int)blend_regs_.size();
    }

    /// Fetch specified region associated to a future blend surface
    RevEngRegion* getBlendReg(int ix)
    {
      if (ix < 0 || ix >= (int)blend_regs_.size())
	return 0;
      else
	return blend_regs_[ix];
    }

    /// Fetch all regions associated to a future blend surface
    void getAllBlendRegs(std::vector<RevEngRegion*>& blend_regs)
    {
      if (blend_regs_.size() > 0)
	blend_regs.insert(blend_regs.end(), blend_regs_.begin(),
			  blend_regs_.end());
    }
    
   /// Fetch adjacent region/surfaces used to create edge
    void getAdjacent(RevEngRegion*& reg1, RevEngRegion*& reg2)
    {
      reg1 = adjacent1_;
      reg2 = adjacent2_;
    }

    /// Fetch information on how the blend surface should be placed related to the
    /// adjacent region surfaces
    void getOuterInfo(bool& out1, bool& out2)
    {
      out1 = outer1_;
      out2 = outer2_;
    }


    /// Fetch intersection curves between adjacent region surfaces
    /// \param first The curve (parameter curve) to select
    /// Expects one curve
    void getCurve(std::vector<shared_ptr<CurveOnSurface> >& cvs, bool first=true)
    {
      if (first && cvs1_.size() > 0)
	cvs.insert(cvs.end(), cvs1_.begin(), cvs1_.end());
      else if ((!first) && cvs2_.size() > 0)
	cvs.insert(cvs.end(), cvs2_.begin(), cvs2_.end());
    }

    /// Get endpoints of intersection curve in geometry space
    void getCrvEndPoints(Point& pos1, Point& pos2);

    /// Fetch the distance from the intesection curve where points related to the blend
    /// surface is to be fetched. This distance is typically larger than the final size
    /// of the blend surface to make sure to get sufficient information to compute the blend.
    double getDistance()
    {
      return distance_;
    }

    /// Fetch estimate of blend radius
    double getRadius()
    {
      return radius_;
    }

    /// Set estimate of blend radius
    void setRadius(double radius)
    {
      radius_ = radius;
    }

    /// The number of times an associated region surface is changed. Indicates a need for
    /// updating the edge
    int getSurfChangeCount()
    {
      return surfchangecount_;
    }

    /// An adjacent region surface has been changed
    void increaseSurfChangeCount()
    {
      surfchangecount_++;
    }

    /// The edge is updated. Reset change count
    void resetSurfChangeCount()
    {
      surfchangecount_ = 0;
    }

    /// The number of times an associated region is increased with more points. Indicates
    /// a possibility for a longer edge
    int getExtendCount()
    {
      return extendcount_;
    }
    
    /// An adjacent region has been extended with more points
     void increaseExtendCount()
    {
      extendcount_++;
    }

    /// The edge is updated. Reset edge count
   void resetExtendCount()
    {
      extendcount_ = 0;
    }

    /// Set associated blend surface 
    void setBlendRegSurf(RevEngRegion* blend)
    {
      defined_blend_ = blend;
    }

    /// Fetch associated blend surface 
    RevEngRegion* getBlendRegSurf()
    {
      return defined_blend_;
    }

    /// Closest point between the intersection curves and the point pos
    void closestPoint(const Point& pos, double& par, Point& close,
		      double& dist);

    /// Closest point between specified intersection curve and the point pos
    void closestPoint(const Point& pos, double& par, Point& close,
		      double& dist, int& ix);

    /// Update intersection curve information (CurveOnSurface) when the adjacent
    /// region reg has been associated another surface new_surf
    void replaceSurf(RevEngRegion* reg, shared_ptr<ParamSurface>& new_surf,
		     double tol);

    /// Check consistency between geometry and parameter space curves and fix
    /// if necessary
    void fixMismatchCurves(double tol);

    /// Remove intersection curves
    void eraseCurves()
    {
      cvs1_.clear();
      cvs2_.clear();
    }

    /// Remove information about one adjacent regions
    void eraseAdjacent(RevEngRegion *adj)
    {
      if (adjacent1_ == adj)
	adjacent1_ = 0;
      else if (adjacent2_ == adj)
	adjacent2_ = 0;
    }

    /// Replace intersection curves
    void replaceCurves(std::vector<shared_ptr<CurveOnSurface> > int_cvs1,
		       std::vector<shared_ptr<CurveOnSurface> > int_cvs2)
    {
      cvs1_.clear();
      cvs2_.clear();
      cvs1_.insert(cvs1_.end(), int_cvs1.begin(), int_cvs1.end());
      cvs2_.insert(cvs2_.end(), int_cvs2.begin(), int_cvs2.end());
    }

    /// Check if the specified edge endpoint lies at the seam of a closed
    /// adjacent region surface
    int closedSfAtEnd(double tol, double& par, Point& pos, bool at_start);

    /// Check if the edge is closed with respect to the tolerance tol
    bool isClosed(double tol);

    /// Check if this edge and the edge other is adjacent and in that case
    /// return the parameter values of the joint
    bool isAdjacent(RevEngEdge* other, double tol, double& par1, double& par2);

    double startparam()
    {
      return cvs1_[0]->startparam();
    }

    double endparam()
    {
      return cvs1_[cvs1_.size()-1]->endparam();
    }

    /// Evaluate edge
    Point point(double par);

    /// Append the edge other to this one. Check if the conditions allow append, turn
    /// curves if nesessary and update associated information
    bool append(RevEngEdge *other, double tol);

    /// Check if the parameter curve is missing for any intersection curve
    int missingParCrv();

    /// Split edge at the seam of an adjacent region surface and update associated
    /// information accordingly.
    /// \param added_edgs New edges due to split. As the adjacent region surfaces may have
    /// more than one seam more than one new edge may be created
    /// \param added_regs If a region associated to a future blend surface is divided
    /// between two edges
    /// \param added_sfs Not used
    void splitAtSeam(double tol, std::vector<shared_ptr<RevEngEdge> >& added_edgs,
		     std::vector<shared_ptr<RevEngRegion> >& added_regs,
		     std::vector<shared_ptr<HedgeSurface> >& added_sfs);

    /// Update the edge if and adjacent region is changed and extendCurve is
    /// not applied
    bool updateCurve(double int_tol, double tol, double len);

    /// Recompute intersection curve if the adjacent regions are extended with
    /// more points. Keep parts of the associated information
    /// \param int_tol Tolerance used in intersection
    /// \param tol Approximation tolerance
    /// \param anglim Used in estimation of edge extent and to compute angular tolerance
    /// \param len Indicates size of environment
    /// \param lenlim Minimum length of intersection curve
    /// \param blendlim Upper limit of with of blend surface
    /// \param added_regions If regions are split to be associated a future blend surface
    /// \param extract_groups If regions are split to maintain connectivity
    /// \param out_sfs Surface to remove since the associated region has been too small
    bool
    extendCurve(double int_tol, double tol, double anglim, 
		double len, double lenlim, double blendlim,
		std::vector<shared_ptr<RevEngRegion> >& added_regions,
		std::vector<std::vector<RevEngPoint*> >& extract_groups,
		std::vector<HedgeSurface*>& out_sfs);

    /// Define trimming curves from edges that will not lead to a blend surface.
    /// Associate regions connected to this blend surface to the adjacent regions as
    /// appropriate. Freed regions and eventual associated surfaces are reported in
    /// out_regs and out_sfs
    void setTrimCurves(double tol, double angtol,
		       std::vector<RevEngRegion*>& out_regs,
		       std::vector<HedgeSurface*>& out_sfs);

    /// Check if the edge other is contained in this (coincidence)
    bool contains(RevEngEdge *other, double tol);

    /// Include the edge other in this (presupposes coincidence) and update parameters
    /// accordingly
    bool integrate(RevEngEdge *other);

    /// Store edge to file
    void store(std::ostream& os);

    /// Read edge from file and collect information necessary to recreate the data structure
    void read(std::istream& is, int& reg_id1, int& reg_id2,
	      int& reg_id3, std::vector<int>& blend_id);

  private:
    /// Unique id for edge. Used in storing and reading data structure to and from file
    int Id_;

    /// First regions where the region surface is intersected to produce this edge
    RevEngRegion* adjacent1_;

    /// Second regions where the region surface is intersected to produce this edge
    RevEngRegion* adjacent2_;

    /// Blend surface associated to the intersection curve between region surfaces
    /// represented in this edge
    RevEngRegion* defined_blend_;

    /// Intersection curve with parameter curve in region surface one
    std::vector<shared_ptr<CurveOnSurface> > cvs1_;
    
    /// Intersection curve with parameter curve in region surface two
    std::vector<shared_ptr<CurveOnSurface> > cvs2_;

    /// Regions associated to a future blend surface corresponding to this edge
    std::vector<RevEngRegion*> blend_regs_;

    /// Type of blend: currently BLEND_NOT_SET or NOT_BLEND
    int blend_type_;
    
    /// The distance from the interection curve in which to search for points associated
    /// to the blend surface
    double distance_;

    /// Estimate of blend radius
    double radius_;

    /// How a blend surface will be placed relative to the outside of region surface one
    bool outer1_;

    /// How a blend surface will be placed relative to the outside of region surface two
    bool outer2_;

    /// Whether any of the adjacent region surfaces is changed
    int surfchangecount_;

    /// Whether any of the adjacent regions is extended with more points
    int extendcount_;

    shared_ptr<RevEngEdge> doSplit(size_t ix, int side, double par, double tol,
				   std::vector<shared_ptr<RevEngRegion> >& added_regs,
				   std::vector<shared_ptr<HedgeSurface> >& added_sfs);
    
    void updateParCurve(RevEngRegion* adj, double int_tol);

  };
}

#endif
