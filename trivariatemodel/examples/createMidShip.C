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

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/geometry/SplineInterpolator.h"
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/creators/LoftSurfaceCreator.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/SweepVolumeCreator.h"
#include "GoTools/trivariate/LoftVolumeCreator.h"
#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/trivariatemodel/VolumeModel.h"
#include "sislP.h"
#include <fstream>

using namespace Go;
using std::vector;

//===========================================================================
//
/// Description: 
///
///              Build a volume model representing a simplified mid ship,
///              hull with stiffeneres in both directions and deck.
///              The model is represented as a block structured volume
///              model.
///
///              First curves describing the inner and outer boundary
///              of a cross section of the hull is constructed. The curves
///              are C1, cubic non-rational B-splines consisting of
///              linear segments and circle segment approximations.	
///              The curves are lofted to produce a section surface of the
///              hull.
///              The hull section is sweeped along the longitudinal axis
///              of the ship to produce the hull volume.
///              The deck if constructed as a loft between sections of the
///              hull sides. Prior to this construction the hull is divided
///              to obey the conditions for block structuring and to split
///              the hull at the lower corners.
///              The stiffeners in the length direction is created by loft
///              after splitting the hull bottom and the deck volume.
///              The stiffeners across the ship is constructed by linear
///              loft between surfaces constructed as Coons patches to 
///              interpolate the already existing hull, deck and stiffeners.
///              Finally, the volumes are collected into a volume model
///              and corresponding coefficients between adjacent volumes
///              are ensured.
///  
/// 
///
/// Input/Output: 
///
///               The model is described by a set of hard coded parameters.
///               The size of these parameters are rather arbitrary.
///               The final volumes are written to the file midship.g2.
///               Some volumes will be right handed and some will be left
///               handed.
///               The current stage of the model is stored to files several
///               times during the construction.
/// 
/// Note:       
///
///             The combination of requiring C1 continuity of section curves 
///             and corresponding parameters at joints between circular and  
///             linear pieces on the inner and outer side of the hull, leads
///             slightly skewed parameterization on some linear parts of the
///             construction.
///             This model is build by manually splitting at block boundaries.
///             Often this process will be necessary to get the wanted
///             model, but it is a procedure that requires the application
///             to keep good track of parameterizations and parameter 
///             directions.
///
//   
//===========================================================================

int main( int argc, char* argv[] )
{
  double eps = 1.0e-4;  // Tolerance used in circle approximation
  // Ship parameters
  double width = 50.0;  // The with of the ship between circular parts
  double height = 30.0; // The height of the ship above circular parts
  double length = 80.0; // The length of the midship
  double radius = 4.0;  // Radius of the circular arc of the outer hull side
  double hull_thick = 0.2;  // Thickness of hull
  double angle = 0.5*M_PI;  // Angle curresponing to the circular arcs
  int nmb_stiff1 = 3;  // Number of stiffeners across the ship
  int nmb_stiff2 = 4;  // Number of stiffeners along the ship
  double height_frac = 0.9;  // Height of stiffeners compared to ship height
  double stiff_thick = 0.1;   // Thickness of stiffeners
  
  double repar = 1.2*(radius-0.5*hull_thick)*angle; // Factor used to
  // reparametrize the circular arcs to get better input to the construction
  // of linear curve segments. The arcs at the outer and inner hull side has
  // curve length radius*angle and (radis-hull_thick)*angle respectively. 
  // The parameter domain of the circular arcs must be the same to get
  // corresponding parameters at joints for the inner and outer curve.
  // The constant factor is included to get a better correspondance of the
  // parameterization of the linear pieces for the inner and outer curve
  double repar_fac = 0.7;  // Used in Hermite interpolation. Explained later
  double len1 = -1.0, len2;  // Length of derivatives at endpoints of the
                             // circular arcs

  // File containing cross section curves at the inner and outer side of the
  // hull
  std::string output_cvs("data/midship_section_curves.g2");

  // The hull volume before deck and stiffneres are added
  std::string output_hull("data/midship_hull.g2");

  // The deck volume prior to splitting due to stiffeners
  std::string output_deck("data/midship_deck.g2");

  // The hull, the deck and the stiffeners along the ship. Block structuring
  // of the current model is performed
  std::string output_stiff1("data/midship_stiff1.g2");

  // The hull, the deck and the stiffeners in both. Block structuring
  // is performed. All blocks are present, but the coefficient correspondance
  // is not performed
  std::string output_stiff2("data/midship_stiff2.g2");

  // The last output file will contain the final volume model for the midship
  std::string output("data/midship.g2");

  // Prepare output files
  std::ofstream of1(output_cvs.c_str());
  std::ofstream of2(output_hull.c_str());
  std::ofstream of3(output_deck.c_str());
  std::ofstream of4(output_stiff1.c_str());
  std::ofstream of5(output_stiff2.c_str());
  std::ofstream outfile(output.c_str());

  // Centre of constuction
  Point midpt(0.0, 0.0, 0.0);
  Point axis(0.0, 1.0, 0.0);  
  Point axis2(1.0, 0.0, 0.0);  
  Point axis3(0.0, 0.0, 1.0);

  // The cross section curves at the outer and inner side of the hull are
  // constructed. The outer curve first
  std::cout << "Making cross section curves" << std::endl;
  vector<shared_ptr<SplineCurve> > hull_cvs(2);
  for (int kr=0; kr<2; ++kr)
    {
      // Make rounded pieces
      // Define centre and startpoints for the radial arcs of the current
      // curve
      Point startpt1 = midpt - 0.5*width*axis2;
      Point centre1 = startpt1 + radius*axis3;
      Point startpt2 = midpt + 0.5*width*axis2;
      Point centre2 = startpt2 + radius*axis3;

      // Use sisl construction of circular arcs to allow for non-rational
      // B-spline curves. Later these curves will be input for Coons patch
      // constructions which do not handle rational curves
      SISLCurve *qc1 = NULL;
      SISLCurve *qc2 = NULL;
      int status = 0;
      s1303(startpt1.begin(), eps, angle, centre1.begin(), axis.begin(),
	    3, &qc1, &status);
      if (status < 0)
	{
	  std::cout << "Error in curve approximation" << std::endl;
	  exit(-1);
	}

      s1303(startpt2.begin(), eps, -angle, centre2.begin(), axis.begin(),
	    3, &qc2, &status);
      if (status < 0)
	{
	  std::cout << "Error in curve approximation" << std::endl;
	  exit(-1);
	}

      // Convert to GoTools format
      shared_ptr<SplineCurve> circ1(SISLCurve2Go(qc1));
      shared_ptr<SplineCurve> circ2(SISLCurve2Go(qc2));
      freeCurve(qc1);  // Note that the sisl entities are alloced by malloc and
      // must consequently be freed by free and not by delete. Use sisl
      // functionality for freeing
      freeCurve(qc2);

      // The circular arcs are now C1 cubic splines on the parameter
      // interval [0,2] and with one inner double knot. Other sizes
      // of the circular arc and other tolerances may give a different
      // number of inner knots. Reparameterize to the tangent lengths
      // more appropriate for the construction of the linear pieces
      circ1->setParameterInterval(0.0, repar);
      circ2->setParameterInterval(0.0, repar);
      circ2->reverseParameterDirection();  // The curves are to follow
      // each other from tail to head along the hull side

      // Create line segments interpolating the end points of the
      // circular arcs with C1 continuitity
      // Prepare for interpolation
      // We use the GoTools interpolation here, but the sisl interpolation
      // is more flexible and might give a better result
      SplineInterpolator interpol;  // Create an empty SplineInterpolator.

      // Mid segment
      // Evaluate boundary information
      vector<Point> bd1(2);
      vector<Point> bd2(2);
      circ1->point(bd1, circ1->startparam(), 1);
      circ2->point(bd2, circ2->endparam(), 1);

      // Set tangent conditions in endpoints
      interpol.setHermiteConditions(bd2[1], bd1[1]);

      // Collect data points and parameterize
      vector<double> pos;
      vector<double> par(2);
      pos.insert(pos.end(), bd2[0].begin(), bd2[0].end());
      pos.insert(pos.end(), bd1[0].begin(), bd1[0].end());
      par[0] = 0.0;
      par[1] = width;  // Curve length parametrization

      // Create an empty spline curve and create content by interpolation
      shared_ptr<SplineCurve> mid(new SplineCurve());
      mid->interpolate(interpol, 2, 3, &par[0], &pos[0]);

      // Side segments
      vector<Point> bd3(2);
      circ1->point(bd3, circ1->endparam(), 1);
      Point top1 = bd3[0] + height*axis3;

      // Hermite interpolation is used also for the side segments although
      // we do only have tangent conditions in one endpoint. However, we know
      // that directions of the tangents must be the same as we want a
      // linear segment represented as a cubic curve. The free length of the
      // tangent at the top of the side curve is set to get a reasonable
      // parameterization compared to linear segments constructed by linear
      // interpolation
      len2 = bd3[1].length();
      if (len1 > 0.0)
	repar_fac *= (len1/len2); // A partial compensation for the fact that
      // the circular arcs have the same parameterization for the inner and
      // outer curve, but different lengths
      interpol.setHermiteConditions(bd3[1], repar_fac*bd3[1]);
      pos.clear();
      pos.insert(pos.end(), bd3[0].begin(), bd3[0].end());
      pos.insert(pos.end(), top1.begin(), top1.end());
      par[1] = height;
      shared_ptr<SplineCurve> side1(new SplineCurve());
      side1->interpolate(interpol, 2, 3, &par[0], &pos[0]);

      // The second side curve
      vector<Point> bd4(2);
      circ2->point(bd4, circ2->startparam(), 1);
      Point top2 = bd4[0] + height*axis3;

      interpol.setHermiteConditions(repar_fac*bd4[1], bd4[1]);
      pos.clear();
      pos.insert(pos.end(), top2.begin(), top2.end());
      pos.insert(pos.end(), bd4[0].begin(), bd4[0].end());
      shared_ptr<SplineCurve> side2(new SplineCurve());
      side2->interpolate(interpol, 2, 3, &par[0], &pos[0]);

      // Represent the curve segments as one curve. No reparametrization
      // is allowed as that would destroy the correspondance of parameters
      // at joints for the inner and outer curve. The lack of reparametrization
      // is also the cause for constructing the side curves as cubic
      // curve with tangent conditions and not just by linear interpolation
      // of the end points
      double dist;
      side2->appendCurve(circ2.get(), 1, dist, false);
      std::cout << " Append error: " << dist << std::endl;
      side2->appendCurve(mid.get(), 1, dist, false);
      std::cout << " Append error: " << dist << std::endl;
      side2->appendCurve(circ1.get(), 1, dist, false);
      std::cout << " Append error: " << dist << std::endl;
      side2->appendCurve(side1.get(), 1, dist, false);
      std::cout << " Append error: " << dist << std::endl;
  
      side2->writeStandardHeader(of1);
      side2->write(of1);

      hull_cvs[kr] = side2;
      midpt[2] += hull_thick;
      radius -= hull_thick;
      len1 = len2;
    }

  // Create surface interpolating the two hull curves
  std::cout << "Making cross section surface" << std::endl;
  shared_ptr<SplineSurface> hull_section(LoftSurfaceCreator::loftSurface(hull_cvs.begin(), 2));

  // Create mid ship hull
  // The first parameter direction is along the section curves, the second
  // parameter direction is across the hull, and the third along the ship
  std::cout << "Hull volume" << std::endl;
  shared_ptr<SplineCurve> sweep_cv(new SplineCurve(midpt, midpt+length*axis));
  shared_ptr<SplineVolume> hull = 
    shared_ptr<SplineVolume>(SweepVolumeCreator::linearSweptVolume(*hull_section,
								   *sweep_cv,
								   midpt));
  // The parametrization of this volume corresponds to the geometric
  // distances in all parameter directions

  hull->writeStandardHeader(of2);
  hull->write(of2);

  // Make deck
  std::cout << "Deck volume" << std::endl;
  // First split hull volume to create corner-to-corner blocks adjacent to
  // the deck. Do also split the hull in the middle of the circular arcs to
  // prepare for further constructions
  // Compute split parameters
  vector<double> single_knots;  // Knot valuess in the curve direction of
  // the volume
  hull->basis(0).knotsSimple(single_knots);
  vector<double> split_par(6);
  split_par[0] = (1.0 - height_frac)*height;    // Top of deck
  split_par[1] = split_par[0] + stiff_thick;  // Bottom of deck
  split_par[2] = single_knots[2];             // Middle of first arc
  split_par[3] = single_knots[5];             // Middle of second arc
  split_par[5] = single_knots[single_knots.size()-1] - (1.0 - height_frac)*height;
                                       // Top of deck on the other side
  split_par[4] = split_par[5] - stiff_thick;  // Bottom of deck

  // Perform splitting
  vector<shared_ptr<SplineVolume> > sub_hulls = 
    hull->split(split_par, 0);
  
  // Create deck volum. First fetch associated boundary surfaces at each
  // hull side. The relevant surface is associated to the maximum parameter 
  // in the second parameter direction
  vector<shared_ptr<SplineSurface> > deck_bd(2);
  deck_bd[0] = sub_hulls[1]->getBoundarySurface(3);
  deck_bd[1] = sub_hulls[5]->getBoundarySurface(3);

  // The two surfaces are oppositely oriented in the first parameter direction,
  // turn the orientation of the second surface
  deck_bd[1]->reverseParameterDirection(true);

  // Create volume by lofting
  // The first parameter direction is in the direction of the hull curve,
  // the second along the ship, and the third across the ship
  // In the third parameter direction the volume is parameterized on the
  // interval [0,1] as this is default in the loft construction. We could
  // choose to reparameterize the volume in this direciton, but do rather
  // take this into account in later constructions
  shared_ptr<SplineVolume> deck(LoftVolumeCreator::loftVolume(deck_bd.begin(), 2));
  deck->writeStandardHeader(of3);
  deck->write(of3);

  // Create stiffeners in the length direction of the ship
  // Split the deck volume and the middle hull volume
  // Note that the parameterization of these volumes differs. Thus,
  // closest point computations is used to find the correspondance
  // First define the split parameters corresponding to the deck
  std::cout << "Stiffeners along the ship" << std::endl;
  vector<double> split_deck(2*nmb_stiff1);
  double ta = deck->startparam(3); // Start parameter in the lofting direction
  double tb = deck->endparam(3);   // End parameter

  // Fetch curve along the lower deck in the front part and compute the
  // geometrical length of this curve
  shared_ptr<SplineSurface> deck_section = deck->getBoundarySurface(2);
  shared_ptr<SplineCurve> deck_cv(deck_section->constParamCurve(deck_section->startparam_u(),
								false));
  double len = deck_cv->estimatedCurveLength(2); // Linear, two points are enough

  // First compute split parameters according to geometrical length
  double del = (len+stiff_thick)/(double)(nmb_stiff1+1);
  double curr = -0.5*stiff_thick;
  for (int kr=0; kr<nmb_stiff1; ++kr)
    {
      curr += del;
      split_deck[2*kr] = curr - 0.5*stiff_thick;
      split_deck[2*kr+1] = curr + 0.5*stiff_thick;
    }

  // Replace by true parameter values
  for (size_t ki=0; ki<split_deck.size(); ++ki)
    split_deck[ki] = ta + split_deck[ki]*(tb - ta)/len;

  // Fetch the corresponding curve on the ship bottom
  shared_ptr<SplineSurface> bottom_section = sub_hulls[3]->getBoundarySurface(4);
  shared_ptr<SplineCurve> bottom_cv(bottom_section->constParamCurve(bottom_section->endparam_v(), 
								    true));
  // Compute the corresponding parameter values on the bottom curve
  vector<double> split_bottom(split_deck.size());
  ta = bottom_cv->startparam();
  tb = bottom_cv->endparam();
  double tc, tdist;
  Point ptclose;
  for (size_t ki=0; ki<split_deck.size(); ++ki)
    {
      Point deck_pos = deck_cv->ParamCurve::point(split_deck[ki]);
      bottom_cv->closestPoint(deck_pos, ta, tb, tc, ptclose, tdist);
      split_bottom[ki] = tc;
    }

  // Perform splitting
  vector<shared_ptr<SplineVolume> > sub_bottom =
    sub_hulls[3]->split(split_bottom, 0);
  sub_hulls.erase(sub_hulls.begin()+3);
  sub_hulls.insert(sub_hulls.begin()+3, sub_bottom.begin(), sub_bottom.end());

  vector<shared_ptr<SplineVolume> > sub_deck =
    deck->split(split_deck, 2);

  // Create stiffeners
  vector<shared_ptr<SplineVolume> > stiffener1(nmb_stiff1);
  for (int kr=0; kr<nmb_stiff1; ++kr)
    {
      // Fetch boundary surfaces related to the current stiffener
      deck_bd[0] = sub_hulls[4+2*kr]->getBoundarySurface(3);
      // First the sub_hulls vector contains the block above the deck,
      // the block corresponding to the deck, the hull side, and the first
      // part of the bottom. Then the stiffeners start
      // The indexing of boundary surface corresponds to maximum parameter
      // in the second parameter direction
      deck_bd[1]= sub_deck[1+2*kr]->getBoundarySurface(1);
      // The indexing of the boundary surface corresponds to the maximum
      // parameter in the first parameter direction
      deck_bd[1]->swapParameterDirection();  // Since the volumes originating
      // from the hull and from the deck have different ordering of the 
      // parameter directions

      // Create volume
      // The parameter directions will be along the hull curve, in the
      // length direction of the ship, and in the height direction
      stiffener1[kr] = 
	shared_ptr<SplineVolume>(LoftVolumeCreator::loftVolume(deck_bd.begin(), 2));
    }

  // Rearrange volumes to store all vertical sections (including those from
  // the hull) in the same vector. The hull sections will be the first
  // and the last entity
  // Make sure to have consistent parameter directions and orientation
  stiffener1.insert(stiffener1.begin(), sub_hulls.begin()+2, sub_hulls.begin()+3);
  sub_hulls.erase(sub_hulls.begin()+2);  // Only one representation of each volume
  stiffener1[0]->swapParameterDirection(0, 1);
  stiffener1[0]->swapParameterDirection(1, 2);
  stiffener1[0]->reverseParameterDirection(2);
  stiffener1.push_back(sub_hulls[sub_hulls.size()-3]);
  sub_hulls.erase(sub_hulls.end()-3);
  stiffener1[stiffener1.size()-1]->swapParameterDirection(0, 1);
  stiffener1[stiffener1.size()-1]->swapParameterDirection(1, 2);
  stiffener1[stiffener1.size()-1]->reverseParameterDirection(0);

  for (size_t ki=0; ki<sub_hulls.size(); ++ki)
    {
      sub_hulls[ki]->writeStandardHeader(of4);
      sub_hulls[ki]->write(of4);
    }
  for (size_t ki=0; ki<sub_deck.size(); ++ki)
    {
      sub_deck[ki]->writeStandardHeader(of4);
      sub_deck[ki]->write(of4);
    }
  for (size_t ki=0; ki<stiffener1.size(); ++ki)
    {
      stiffener1[ki]->writeStandardHeader(of4);
      stiffener1[ki]->write(of4);
    }
  
  // Prepare for construction of the stiffeners across the ship
  // Split all current volumes. The parameterization in the relevent
  // parameter direction corresponds to the geometric length and all 
  // volumes are linear in this direction. Thus, no matching of
  // parameter values are required
  std::cout << "Stiffeners across the ship" << std::endl;
  vector<double> split_side(2*(nmb_stiff2-1));
  ta = sub_hulls[0]->startparam(3);
  tb = sub_hulls[0]->endparam(3);
  del = (tb - ta - nmb_stiff2*stiff_thick)/(double)(nmb_stiff2-1);
  curr = ta + stiff_thick;
  for (int kr=0; kr<nmb_stiff2-1; ++kr)
    {
      split_side[2*kr] = curr;
      curr += del;
      split_side[2*kr+1] = curr;
      curr += stiff_thick;
    }
  int nmb_split = (int)split_side.size() + 1;

  // Remaining hull volumes
  vector<shared_ptr<SplineVolume> > sub_hulls2;
  for (size_t ki=0; ki<sub_hulls.size(); ++ki)
    {
      vector<shared_ptr<SplineVolume> > curr_sub = 
	sub_hulls[ki]->split(split_side, 2);
      // Split in the third parameter direction as this is along the ship
      sub_hulls2.insert(sub_hulls2.end(), curr_sub.begin(), curr_sub.end());
    }

  // Extended set of stiffeners
  vector<shared_ptr<SplineVolume> > sub_stiffeners;
  for (size_t ki=0; ki<stiffener1.size(); ++ki)
    {
      vector<shared_ptr<SplineVolume> > curr_sub = 
	stiffener1[ki]->split(split_side, 1);
      // Split in the second parameter direction as this is along the ship
      sub_stiffeners.insert(sub_stiffeners.end(), curr_sub.begin(), curr_sub.end());
    }
  
  // Deck
  vector<shared_ptr<SplineVolume> > sub_deck2;
  for (size_t ki=0; ki<sub_deck.size(); ++ki)
    {
      vector<shared_ptr<SplineVolume> > curr_sub = 
	sub_deck[ki]->split(split_side, 1);
      // Split in the second parameter direction as this is along the ship
      sub_deck2.insert(sub_deck2.end(), curr_sub.begin(), curr_sub.end());
    }
   
  // Create the last volumes
  vector<shared_ptr<SplineVolume> > stiffener2;

  int perm[4];  // Translates from ccw numbering of boundary curves to
  // the numbering: umin, umax, vmin, vmax
  perm[0] = 3;
  perm[1] = 1;
  perm[2] = 0;
  perm[3] = 2;
  // Along the ship
  for (int kr=0; kr<nmb_stiff2; ++kr)
    {
      // Across the ship
      for (int kj=0; kj<=nmb_stiff1; ++kj)
  	{
  	  // Since no volume creation method where you can give the
  	  // boundary surfaces in two parameter directions, but not the
  	  // last, exist, we construct the surfaces on both sides of the
  	  // volume and loft in the length direction of the ship
  	  vector<shared_ptr<SplineSurface> > bd_sf(2);

	  // Fetch the surfaces surrounding the new stiffener volumes
	  // 1. stiffener surface
	  shared_ptr<SplineSurface> stiff1 = 
	    sub_stiffeners[kj*nmb_split+2*kr]->getBoundarySurface(1);
	  // Bottom surface
	  shared_ptr<SplineSurface> bottom = 
	    sub_hulls2[2*(1+kj)*nmb_split+2*kr]->getBoundarySurface(3);
	  // 2. stiffener surface
	  shared_ptr<SplineSurface> stiff2 = 
	    sub_stiffeners[(kj+1)*nmb_split+2*kr]->getBoundarySurface(0);
	  // Deck surface
	  shared_ptr<SplineSurface> decksf = 
	    sub_deck2[2*kj*nmb_split+2*kr]->getBoundarySurface(1);

  	  for (int kh=0; kh<2; ++kh)
  	    {
	      // Fetch curves associated with the current boundary surfaces
	      // The curves are expected to have a sequence around the
	      // surface to be constructed and be oriented with the
	      // head of one curve following the tail of the previous curve
	      vector<shared_ptr<ParamCurve> > coons_cvs(4);
	      // umin
	      coons_cvs[0] = 
		shared_ptr<ParamCurve>(stiff1->edgeCurve(perm[kh]));
	      coons_cvs[0]->reverseParameterDirection();
	      // vmin
	      coons_cvs[1] =
		shared_ptr<ParamCurve>(bottom->edgeCurve(perm[2+kh]));
	      // umin
	      coons_cvs[2] =
		shared_ptr<ParamCurve>(stiff2->edgeCurve(perm[kh]));
	      // umin
	      coons_cvs[3] =
		shared_ptr<ParamCurve>(decksf->edgeCurve(perm[kh]));
	      coons_cvs[3]->reverseParameterDirection();

	      // Create surface
	      CurveLoop loop(coons_cvs, eps);
	      bd_sf[kh] = 
		shared_ptr<SplineSurface>(CoonsPatchGen::createCoonsPatch(loop));
  	    }
	  // Create volume
  	  shared_ptr<SplineVolume> vol(LoftVolumeCreator::loftVolume(bd_sf.begin(),
  	  							     2));
  	  stiffener2.push_back(vol);
  	}
    }

  // Collect all constructed volumes
  vector<shared_ptr<ParamVolume> > all_blocks;
  all_blocks.insert(all_blocks.end(), sub_hulls2.begin(), sub_hulls2.end());
  all_blocks.insert(all_blocks.end(), sub_stiffeners.begin(), sub_stiffeners.end());
  all_blocks.insert(all_blocks.end(), stiffener2.begin(), stiffener2.end());
  all_blocks.insert(all_blocks.end(), sub_deck2.begin(), sub_deck2.end());

  for (size_t ki=0; ki<all_blocks.size(); ++ki)
    {
      all_blocks[ki]->writeStandardHeader(of5);
      all_blocks[ki]->write(of5);
    }

  // All volume blocks are now created, but they do not have common
  // spline spaces along all common boundaries. Collect the volumes into a
  // volume model and ensure common spline spaces and corresponding
  // coefficients
  // First add topology stuctures to the volumes
  vector<shared_ptr<ftVolume> > blocks;  // Storage for volumes including
                                         // topology information
  // Topology tolerances
  double neighbour = 0.01;
  double gap = 0.0001;
  double bend = 0.5;
  double kink = 0.01;
  double approxtol = 0.01;
  for (size_t ki=0; ki<all_blocks.size(); ++ki)
    {
      shared_ptr<ftVolume> ftvol = 
	shared_ptr<ftVolume>(new ftVolume(all_blocks[ki], gap, kink));
      blocks.push_back(ftvol);
    }
      
  // Create a volume model. The topology build performs coincidence testing
  // between possible adjacent volumes. This may be time consuming
  std::cout << "Volume topology build"  << std::endl;
  shared_ptr<VolumeModel> volmodel = 
    shared_ptr<VolumeModel>(new VolumeModel(blocks, gap, neighbour, 
					    kink, 10.0*kink));


  // Ensure common spline spaces and corresponding coefficients
  std::cout << "Ensure corresponding coefficients" << std::endl;
  volmodel->makeCommonSplineSpaces();
  volmodel->averageCorrespondingCoefs();

  // Fetch the number of volumes
  int nmb_vol = volmodel->nmbEntities();

  // Write the volumes to the output file
  // With the current file format, the topology information is lost
  for (int kr=0; kr<nmb_vol; ++kr)
    {
      shared_ptr<ParamVolume> vol = volmodel->getVolume(kr);
       vol->writeStandardHeader(outfile);
       vol->write(outfile);
     }

}
